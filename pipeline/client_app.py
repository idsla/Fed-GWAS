# client/main_client.py

import numpy as np
import sys
import os
import math
import time
from datetime import datetime
from typing import Dict, Optional, List
from flwr.client import ClientApp
from flwr.common import Context, ConfigRecord, RecordDict
from pipeline.clients.logger_manager import LoggerManager

# Add parent directory to path for server imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pipeline.clients.base_client import BaseGWASClient
from pipeline.server.prg_masking import create_client_masking_helper
from pipeline.clients.local_qc import (
    compute_genotype_counts,
    compute_missingness_counts,
    run_local_lr,
    parse_insignificant_snps,
    exclude_snps,
    exclude_samples_by_missing_rate,
)
from pipeline.clients.iterative_king import handle_iterative_king
from pipeline.clients.iterative_lr import handle_iterative_lr, compute_lr_chunk_size
from pipeline.clients.data_loder import DataLoader
from pipeline.clients.client_qc_aggregator import _compute_exclusion_list
from pipeline.clients.client_to_client import create_client_messenger

from pipeline.clients.flwr_config import get_all_public_keys
from pipeline.clients.c2c_payloads import pack_typed_payload_uint8, unpack_typed_payload_uint8
from pipeline.clients.lr_privacy import lr_token, resolve_tokens_to_snp_ids, make_dummy_tokens
from pipeline.clients.seed_sync import (
    get_persist_base_dir,
    key_paths,
    load_persisted_seed,
    persist_seed,
    load_persisted_keys_into_masking_helper,
    compute_global_seed_deterministic,
    compute_global_seed_from_encrypted_messages,
    resolve_my_hash_id,
    ensure_global_seed,
)
from pipeline.utils.monitoring_config import resolve_monitoring_settings
from pipeline.utils.performance.monitoring_runtime import (
    ClientPerformanceTracker,
    monitored_stage,
)


def _resolve_simulation_client_config_path(config_path, center_id: int):
    """Resolve center config for both CLI flat layout and legacy nested layout.

    `fedgwas-sim setup-experiment` writes `configs/center_1.yaml`, while older
    repository experiments use `configs/center_1/config.yaml`. Flower only
    receives `config_path` through run-config, so the client app normalizes that
    value at startup before constructing `FedLRClient`.
    """
    import pathlib

    base = pathlib.Path(config_path)
    if base.name == "configs":
        candidates = [
            base / f"center_{center_id}.yaml",
            base / f"center_{center_id}" / "config.yaml",
        ]
    else:
        candidates = [
            base / "configs" / f"center_{center_id}.yaml",
            base / "configs" / f"center_{center_id}" / "config.yaml",
        ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


class FedLRClient(BaseGWASClient):
    # Helper utilities are intentionally implemented in `pipeline/clients/*` modules to keep this file small.
    
    def __init__(
        self, 
        partition_id,
        context: Context,
        config_file="config.yaml", 
        partition_by="samples"
    ):
        
        # Store context and client_state for easier access
        self.context = context
        self.client_state = (context.state)
        self.partition_id = partition_id
        self.partition_by = partition_by
        self.client_id = f"client_{partition_id}"
        self._config_file_path = config_file

        # Determine simulation mode from context (needed early to decide output layout)
        simulation_mode = context.run_config.get("simulation", False)

        # Check if this is first initialization or recovery from state
        if not self._is_client_initialized():
        
            # First initialization - use DataLoader
            loader = DataLoader(config_file)
            loader.transform_data()
            self._apply_run_scoped_intermediate(loader, is_first_init=True)

            # Ensure output directories exist (files use unique prefixes like client_1_log.txt to avoid collisions)
            os.makedirs(loader.intermediate_dir, exist_ok=True)
            os.makedirs(loader.log_dir, exist_ok=True)
            
            # Store all configuration to state
            self.client_state = self._store_config_to_state(
                self.client_state, loader, self.client_id, self.partition_by
            )
            self.is_first_init = True

            self.log_dir = loader.log_dir
            self.logger = LoggerManager.get_logger(self.client_id, self.log_dir, is_first_init=True)

        else:
            loader = self._load_config_from_state(self.client_state.config_records)
            self._apply_run_scoped_intermediate(loader, is_first_init=False)
            self.is_first_init = False

            # Ensure output directories exist (files use unique prefixes like client_1_log.txt to avoid collisions)
            os.makedirs(loader.intermediate_dir, exist_ok=True)
            os.makedirs(loader.log_dir, exist_ok=True)

            self.log_dir = loader.log_dir
            self.logger = LoggerManager.get_logger(self.client_id, self.log_dir, is_first_init=False)
                
        self.plink_prefix = loader.plink_prefix
        self.intermediate_dir = loader.intermediate_dir
        
        # Allow forcing phenotype fix in deployment mode for toy/sim experiments
        phenotype_fix_missing_to_case = context.run_config.get("phenotype_fix_missing_to_case", False)
        
        # Load sample_offset from config parameters
        sample_offset = loader.parameters.get("sample_offset") if hasattr(loader, 'parameters') and loader.parameters else None
        
        # Initialize base client with loaded configuration
        super().__init__(
            self.plink_prefix, 
            client_id=self.client_id, 
            partition_by=partition_by,
            simulation_mode=simulation_mode,
            phenotype_fix_missing_to_case=phenotype_fix_missing_to_case,
            sample_offset=sample_offset,
        )

        # If phenotype-fix is enabled, use a derived prefix (do not mutate original input files)
        try:
            fixed_prefix = os.path.join(self.intermediate_dir, f"phenofixed_{self.client_id}")
            self.plink_prefix = self.ensure_phenotype_fixed_prefix(fixed_prefix)
        except Exception as e:
            # Never hard-fail the whole client just due to phenotype-fix; log and continue
            try:
                self.logger.warning(f"[Client {self.client_id}] Phenotype-fix setup failed: {e}")
            except Exception:
                pass

        # Load thresholds from configuration
        thresholds = loader.get_thresholds()
        # Ensure all thresholds are floats (config values may be strings from YAML)
        self.maf_threshold = float(thresholds.get("maf_threshold", 0.01))
        self.miss_threshold = float(thresholds.get("missing_threshold", 0.05))  # Fixed default from 0.1 to 0.05
        self.hwe_threshold = float(thresholds.get("hwe_threshold", 1e-6))
        # Support separate thresholds for local and global LR, with backward compatibility
        base_p_threshold = float(thresholds.get("p_threshold", 5e-3))
        self.local_lr_threshold = float(thresholds.get("local_lr_threshold", base_p_threshold))
        self.global_lr_threshold = float(thresholds.get("global_lr_threshold", base_p_threshold))
        # Keep p_threshold for backward compatibility (defaults to global_lr_threshold)
        self.p_threshold = float(thresholds.get("p_threshold", self.global_lr_threshold))
        self.king_threshold = float(thresholds.get("king_threshold", 0.4))  # Load KING threshold from config
        # Store thresholds dict for later use (convert values to float)
        self.thresholds = {k: float(v) if isinstance(v, (int, float, str)) else v for k, v in thresholds.items()}
        
        # Load other configuration from loader
        self.flower_config = loader.get_flower_config()
        self.parameters = loader.get_parameters()
        self.participation = loader.get_participation()

        config_path = None
        try:
            config_path = context.run_config.get("config_path")
        except Exception:
            pass
        self.monitoring_settings = resolve_monitoring_settings(
            center_config_file=self._config_file_path,
            config_path=str(config_path) if config_path else None,
        )
        experiment_name = self.monitoring_settings.get("experiment_name", "fedgwas")
        self._perf_tracker = None
        if (
            self.monitoring_settings.get("enable_performance_monitoring")
            or self.monitoring_settings.get("enable_network_monitoring")
        ):
            self._perf_tracker = ClientPerformanceTracker(
                self.client_id,
                self.log_dir,
                self.monitoring_settings,
                experiment_name,
            )
            self.logger.info(
                f"[Client {self.client_id}] Performance monitoring enabled "
                f"(output: {self.log_dir})"
            )
        
        # For final LR significance, if desired
        self.lr_final = {}
        
        # For accumulating partial LR p-values
        self.lr_pvals = {}
        
        # PRG-MASKING helper (will be initialized when we know num_clients)
        self.masking_helper = None
        self.num_clients = 3  # Default fallback; should be sourced from server config.
        
        # Initialize king_accumulator for iterative KING analysis
        self.king_accumulator = {}

    def _store_config_to_state(
        self, 
        client_state: RecordDict, 
        loader: DataLoader, 
        client_id: str, 
        partition_by: str
    ):
        """
        Store all client configuration to state for future rounds.
        
        Args:
            loader: DataLoader instance with configuration
            client_id: Client identifier
            partition_by: Partitioning strategy
        """ 
        # Store in state
        client_state.config_records["client_params"] = ConfigRecord({
            "client_id": client_id,
            "partition_by": partition_by,
            'initialized': True
        })
        
        # Convert DataLoader to config records
        config_dict = loader.to_config_records()
        client_state.config_records.update(config_dict)
        
        return client_state

    def _apply_run_scoped_intermediate(self, loader: DataLoader, is_first_init: bool) -> None:
        """
        Optionally scope intermediate_dir to a per-run subdirectory to avoid map collisions.
        Uses run_id/run_tag from run_config or parameters; otherwise auto-generates on first init.
        """
        run_id = None
        try:
            run_id = self.context.run_config.get("run_id") or self.context.run_config.get("run_tag")
        except Exception:
            run_id = None
        if run_id is None and getattr(loader, "parameters", None):
            run_id = loader.parameters.get("run_id") or loader.parameters.get("run_tag")
        if not run_id and is_first_init:
            run_id = datetime.now().strftime("%Y%m%d_%H%M%S")

        if not run_id:
            return

        base_dir = loader.intermediate_dir
        base_name = os.path.basename(base_dir.rstrip(os.sep))
        if base_name.startswith("run_"):
            base_dir = os.path.dirname(base_dir.rstrip(os.sep))
        run_dir = os.path.join(base_dir, f"run_{run_id}")
        if loader.intermediate_dir != run_dir:
            loader.intermediate_dir = run_dir
        os.makedirs(loader.intermediate_dir, exist_ok=True)

        # Persist updated intermediate_dir when reusing state
        if not is_first_init and getattr(self, "client_state", None) is not None:
            try:
                paths_rec = self.client_state.config_records.get("client_paths")
                paths_dict = dict(paths_rec) if paths_rec else {}
                if paths_dict.get("intermediate_dir") != loader.intermediate_dir:
                    paths_dict["intermediate_dir"] = loader.intermediate_dir
                    self.client_state.config_records["client_paths"] = ConfigRecord(paths_dict)
            except Exception as e:
                try:
                    self.logger.warning(f"[Client {self.client_id}] Failed to persist run-scoped intermediate_dir: {e}")
                except Exception:
                    pass

    def _persist_plink_prefix(self, new_prefix: str):
        """
        Persist the current plink prefix into client_state so it survives across rounds/processes.
        """
        # Ensure absolute path for consistency
        if not os.path.isabs(new_prefix):
            new_prefix = os.path.abspath(new_prefix)
        
        self.plink_prefix = new_prefix
        if not self.client_state:
            self.logger.warning(f"[Client {self.client_id}] No client_state available, cannot persist plink_prefix")
            return
        try:
            paths_rec = self.client_state.config_records.get("client_paths")
            paths_dict = dict(paths_rec) if paths_rec else {}
            paths_dict["plink_prefix"] = new_prefix
            self.client_state.config_records["client_paths"] = ConfigRecord(paths_dict)
            self.logger.info(f"[Client {self.client_id}] Persisted plink_prefix to state: {new_prefix}")
        except Exception as e:
            self.logger.warning(f"[Client {self.client_id}] Failed to persist plink prefix to state: {e}")
    
    def _restore_plink_prefix(self):
        """
        Restore the plink prefix from client_state if it exists (e.g., after QC filtering).
        """
        current_prefix = self.plink_prefix
        if not self.client_state:
            self.logger.warning(f"[Client {self.client_id}] No client_state available, cannot restore plink_prefix")
            return
        try:
            paths_rec = self.client_state.config_records.get("client_paths")
            if paths_rec:
                paths_dict = dict(paths_rec)
                persisted_prefix = paths_dict.get("plink_prefix")
                if persisted_prefix:
                    # Ensure absolute path
                    if not os.path.isabs(persisted_prefix):
                        persisted_prefix = os.path.abspath(persisted_prefix)
                    
                    # Check if the persisted prefix exists
                    if os.path.exists(persisted_prefix + ".bim"):
                        self.plink_prefix = persisted_prefix
                        self.logger.info(f"[Client {self.client_id}] Restored plink_prefix from state: {current_prefix} -> {self.plink_prefix}")
                    else:
                        self.logger.warning(f"[Client {self.client_id}] Persisted prefix {persisted_prefix} does not exist, keeping current prefix {current_prefix}")
                else:
                    self.logger.warning(f"[Client {self.client_id}] No plink_prefix found in state, keeping current prefix {current_prefix}")
            else:
                self.logger.warning(f"[Client {self.client_id}] No client_paths in state, keeping current prefix {current_prefix}")
        except Exception as e:
            self.logger.warning(f"[Client {self.client_id}] Failed to restore plink prefix from state: {e}, keeping current prefix {current_prefix}")
    
    def _load_config_from_state(self, client_config_records: dict):
        """
        Load client configuration from state.
        
        Returns: 
            dict: Dictionary containing all client configuration
        """        
        # Recreate DataLoader from stored config
        loader = DataLoader.from_config_records(client_config_records)
        
        
        return loader

    def _is_client_initialized(self) -> bool:
        """
        Check if client has been initialized before.
        
        Returns:
            bool: True if client was previously initialized
        """
        if not self.client_state:
            return False
            
        return (
            "client_params" in self.client_state.config_records and
            'initialized' in self.client_state.config_records["client_params"] and
            self.client_state.config_records["client_params"]["initialized"]
        )

    def fit(self, parameters, config):

        stage = config.get("stage", "unknown")
        self.logger.info(f"[Client {self.client_id}] Stage: {stage}")

        try:
            # Extract stage from config
            stage = config.get("stage", "unknown")
            # Debug-level logging for verbose information
            self.logger.debug(f"[Client {self.client_id}] FIT CALLED - stage: {stage}")
            self.logger.debug(f"[Client {self.client_id}] Config: {config}")
            self.logger.debug(f"[Client {self.client_id}] Parameters length: {len(parameters)}")
            self.logger.debug(f"[Client {self.client_id}] Config keys: {list(config.keys()) if hasattr(config, 'keys') else 'no keys'}")

            # ------------------------------------------------------------
            # Derive num_clients early (before masking_helper is created)
            # ------------------------------------------------------------
            cfg_n = config.get("num_clients")
            if isinstance(cfg_n, int) and cfg_n > 0:
                self.num_clients = cfg_n
            else:
                keys = get_all_public_keys(config)
                if keys:
                    self.num_clients = len(keys)
            
            # Exit process if this client opts out of the current stage
            if not self.participation.get(stage, True):
                self._maybe_cluster_reply_skew_sleep()
                return [], 0, {"skipped": True, "stage": stage}

            with monitored_stage(self._perf_tracker, stage):
                result = self._fit_stage_body(parameters, config, stage)
            self._maybe_cluster_reply_skew_sleep()
            return result

        except Exception as e:
            self.logger.exception(
                f"[Client {self.client_id}] fit failed at stage={config.get('stage', 'unknown')}: {e}"
            )
            self._maybe_cluster_reply_skew_sleep()
            # Always return a valid NumPyClient 3-tuple so Flower SuperNode can reply.
            return [], 0, {
                "fit_error": str(e),
                "stage": str(config.get("stage", "unknown")),
                "failed": True,
            }

    def _maybe_cluster_reply_skew_sleep(self) -> None:
        """Avoid Flower 1.19 SuperNode 'Invalid out message' when reply.created_at <= instruction.created_at."""
        raw = os.environ.get("FEDGWAS_FLR_REPLY_SKEW_SEC", "")
        if not raw:
            return
        try:
            delay = float(raw)
        except ValueError:
            return
        if delay > 0:
            time.sleep(delay)

    def _fit_stage_body(self, parameters, config, stage):
            # Initialize masking helper if needed
            if self.masking_helper is None:
                # Convert client_id to integer for PRG masking
                try:
                    # Try to extract partition_id from client_id (e.g., "client_1" -> 1)
                    if self.client_id.startswith("client_"):
                        client_id_int = int(self.client_id.split("_")[1])
                    else:
                        client_id_int = int(self.client_id)
                except (ValueError, TypeError):
                    # If conversion fails, use hash of string as integer
                    client_id_int = abs(hash(self.client_id)) % (10**9)
                self.masking_helper = create_client_masking_helper(client_id_int, self.num_clients)
                # If actor was restarted, reload persisted ECC keys/seed so encryption/decryption still works
                base_dir = get_persist_base_dir(getattr(self, "log_dir", None), getattr(self, "intermediate_dir", None))
                pem = load_persisted_keys_into_masking_helper(self.masking_helper, base_dir)
                if pem:
                    self.my_public_key_pem = pem
                s = load_persisted_seed(base_dir)
                if s is not None:
                    self.global_seed = s
            else:
                base_dir = get_persist_base_dir(getattr(self, "log_dir", None), getattr(self, "intermediate_dir", None))
                if getattr(self.masking_helper, "private_key", None) is None:
                    pem = load_persisted_keys_into_masking_helper(self.masking_helper, base_dir)
                    if pem:
                        self.my_public_key_pem = pem
                if getattr(self, "global_seed", None) is None:
                    s = load_persisted_seed(base_dir)
                    if s is not None:
                        self.global_seed = s

            ################################################################################
            # Stage 1: Key Exchange - Generate and send DH public key for secure aggregation
            ################################################################################
            if stage == "key_exchange":
                self.logger.info(f"[Client {self.client_id}] Stage: key_exchange")
                curve_params = {k: v for k, v in config.items() if k == "curve"}
                public_key_pem = self.masking_helper.generate_ecc_keypair(curve_params)
                # Store our public key PEM so we can later resolve our hash-based ID from all_public_keys
                self.my_public_key_pem = public_key_pem

                # Persist keypair for later stages (or actor restarts)
                try:
                    base_dir = get_persist_base_dir(getattr(self, "log_dir", None), getattr(self, "intermediate_dir", None))
                    priv_path, pub_path = key_paths(base_dir)
                    if self.masking_helper.private_key is not None:
                        with open(priv_path, "w") as f:
                            f.write(self.masking_helper.private_key.export_key(format="PEM"))
                    with open(pub_path, "w") as f:
                        f.write(public_key_pem)
                except Exception:
                    pass
                
                # Convert PEM string to bytes and then to uint8 array for transmission
                public_key_bytes = public_key_pem.encode('utf-8')
                public_key_array = np.frombuffer(public_key_bytes, dtype=np.uint8)
                
                self.logger.info(f"[Client {self.client_id}] Generated DH public key")
                return [public_key_array], 1, {}
            
            ################################################################################
            # Stage 2: Sync - Send encrypted seed shares for global random seed generation
            ################################################################################
            elif stage == "sync":
                self.logger.info(f"[Client {self.client_id}] Stage: sync with encrypted client-to-client communication")
                
                # Get all public keys and compute shared secrets
                all_public_keys = get_all_public_keys(config)
                curve_params = {k: v for k, v in config.items() if k == "curve"}
                
                if not curve_params and "curve" not in config:
                    # Extract curve params from the strategy (they should be in config)
                    self.logger.warning(f"[Client {self.client_id}] Missing curve params in sync stage")
                    # Do not leak plaintext seeds; let sync_response use deterministic fallback
                    return [], 1, {"sync_mode": "no_curve_params"}
                
                # correct sync protocol: encrypted broadcast of local_seed
                # Each client encrypts its *local_seed* for every recipient.
                # Server forwards ciphertexts but cannot decrypt; each client decrypts all seeds and sums.

                if not all_public_keys:
                    self.logger.warning(f"[Client {self.client_id}] Missing all_public_keys in sync stage; cannot encrypt broadcast. Falling back to deterministic seed in sync_response.")
                    return [], 1, {"sync_mode": "no_public_keys"}

                my_hash_id = resolve_my_hash_id(all_public_keys, getattr(self, "my_public_key_pem", None))
                if my_hash_id is None:
                    self.logger.warning(f"[Client {self.client_id}] Could not resolve my hash_id in sync stage; cannot encrypt broadcast. Falling back to deterministic seed in sync_response.")
                    return [], 1, {"sync_mode": "no_hash_id"}

                if not (self.masking_helper and self.masking_helper.private_key):
                    self.logger.warning(f"[Client {self.client_id}] No private key available, cannot encrypt broadcast. Falling back to deterministic seed in sync_response.")
                    return [], 1, {"sync_mode": "no_private_key"}

                messenger = create_client_messenger(
                    my_hash_id,
                    self.num_clients,
                    self.masking_helper.private_key
                )

                public_key_ids = sorted([int(k) for k in all_public_keys.keys()])
                seed_array = np.array([float(self.local_seed)], dtype=np.float64)
                encrypted_arrays = []

                for recipient_hash_id in public_key_ids:
                    recipient_key_pem = all_public_keys.get(str(recipient_hash_id))
                    if not recipient_key_pem:
                        # Never downgrade to plaintext/pickle: abort and use deterministic fallback in sync_response
                        self.logger.warning(f"[Client {self.client_id}] No public key for recipient {recipient_hash_id}; aborting encrypted broadcast and falling back to deterministic seed in sync_response.")
                        return [], 1, {"sync_mode": "incomplete_public_keys"}
                    msg = messenger.encrypt_for_recipient(seed_array, recipient_hash_id, recipient_key_pem)
                    encrypted_arrays.append(np.frombuffer(msg, dtype=np.uint8))

                self.logger.info(f"[Client {self.client_id}] ✓ Encrypted-broadcasted local_seed to {len(encrypted_arrays)} recipients (server cannot decrypt)")
                return encrypted_arrays, 1, {}
            

            ################################################################################
            # Stage 3: Sync Response - Compute global seed from encrypted shares
            ################################################################################
            elif stage == "sync_response":
                self.logger.info(f"[Client {self.client_id}] Stage: sync_response - computing global seed")
                
                # Compute global seed from encrypted shares forwarded by server
                # OR use deterministic method if sync stage was skipped (Flower bug)
                # This completes the sync process: clients now have both public keys and global seed
                # Ensure we can parse public keys from JSON config
                if "all_public_keys" not in config:
                    config = dict(config)
                    config["all_public_keys"] = get_all_public_keys(config)
                if not hasattr(self, 'global_seed') or self.global_seed is None:
                    if len(parameters) > 0:
                        # Server forwarded encrypted shares - decrypt and aggregate locally
                        self.logger.info(f"[Client {self.client_id}] Computing global_seed from encrypted shares")
                        all_public_keys = config.get("all_public_keys", {})
                        self.global_seed = compute_global_seed_from_encrypted_messages(
                            parameters=parameters,
                            all_public_keys=all_public_keys,
                            my_public_key_pem=getattr(self, "my_public_key_pem", None),
                            private_key=getattr(self.masking_helper, "private_key", None),
                            num_clients=self.num_clients,
                            local_seed=self.local_seed,
                            server_round=int(config.get("server_round", 0) or 0),
                        )
                    else:
                        # Fallback to deterministic method (sync stage was skipped)
                        self.logger.info(f"[Client {self.client_id}] Computing global_seed deterministically (sync stage was skipped)")
                        all_public_keys = config.get("all_public_keys", {})
                        self.global_seed = compute_global_seed_deterministic(all_public_keys, int(config.get("server_round", 0) or 0), self.local_seed)
                
                # Verify global seed is set
                if self.global_seed is not None:
                    self.logger.info(f"[Client {self.client_id}] ✓ Sync complete: global_seed={self.global_seed} (ready for SNP ID anonymization)")
                    # Persist for later stages/actor restarts
                    base_dir = get_persist_base_dir(getattr(self, "log_dir", None), getattr(self, "intermediate_dir", None))
                    persist_seed(base_dir, int(self.global_seed))
                else:
                    self.logger.warning(f"[Client {self.client_id}] ⚠ Global seed not computed, using fallback")
                    # Final fallback - compute deterministically
                    all_public_keys = config.get("all_public_keys", {})
                    self.global_seed = compute_global_seed_deterministic(all_public_keys, int(config.get("server_round", 0) or 0), self.local_seed)
                
                return [], 1, {}
            
            ################################################################################
            # Stage 4: Local QC - Filter samples by per-sample missing rate threshold
            ################################################################################
            elif stage == "local_qc":
                # Global seed should already be computed in sync_response stage
                # But if sync_response was skipped, compute it here
                if not hasattr(self, 'global_seed') or self.global_seed is None:
                    # Try to load persisted global_seed first (actor restarts)
                    base_dir = get_persist_base_dir(getattr(self, "log_dir", None), getattr(self, "intermediate_dir", None))
                    s = load_persisted_seed(base_dir)
                    if s is not None:
                        self.global_seed = s
                if not hasattr(self, 'global_seed') or self.global_seed is None:
                    self.logger.warning(f"[Client {self.client_id}] Global seed not set (sync_response may have been skipped), computing deterministically")
                    # Try to get public keys from config if available
                    all_public_keys = config.get("all_public_keys", {})
                    if not all_public_keys:
                        # Try to get from masking_helper if available
                        if hasattr(self, 'masking_helper') and self.masking_helper:
                            # Public keys might be stored elsewhere, use deterministic method
                            pass
                    if not all_public_keys:
                        all_public_keys = get_all_public_keys(config)
                    self.global_seed = compute_global_seed_deterministic(all_public_keys, int(config.get("server_round", 0) or 0), self.local_seed)
                    if self.global_seed is not None:
                        self.logger.info(f"[Client {self.client_id}] ✓ Computed global_seed={self.global_seed} in local_qc (fallback)")
                        base_dir = get_persist_base_dir(getattr(self, "log_dir", None), getattr(self, "intermediate_dir", None))
                        persist_seed(base_dir, int(self.global_seed))
                
                local_mind_threshold = config.get("local_mind_threshold", 0.1)
                self.logger.info(f"[Client {self.client_id}] Stage: local_qc (mind={local_mind_threshold})")
                # Exclude samples with missing rate > local_mind_threshold
                new_prefix = exclude_samples_by_missing_rate(
                    self.plink_prefix, local_mind_threshold, log_dir=self.log_dir
                )
                self.plink_prefix = new_prefix
                self.logger.info(f"[Client {self.client_id}] Local QC done => new prefix {self.plink_prefix}")
                return [], 1, {}

            ################################################################################
            # Stage 4: Global QC - Compute and send masked MAF, HWE, and missingness statistics
            ################################################################################
            elif stage == "global_qc":
                
                self.logger.info(f"[Client {self.client_id}] Stage: global_qc with PRG-MASKING")
                
                # Initialize masking helper if needed (may be None if key_exchange was skipped)
                if self.masking_helper is None:
                    self.logger.info(f"[Client {self.client_id}] Initializing masking helper for global_qc stage")
                    try:
                        if self.client_id.startswith("client_"):
                            client_id_int = int(self.client_id.split("_")[1])
                        else:
                            client_id_int = int(self.client_id)
                    except (ValueError, TypeError):
                        client_id_int = abs(hash(self.client_id)) % (10**9)
                    self.masking_helper = create_client_masking_helper(client_id_int, self.num_clients)
                
                # Get all public keys and compute shared secrets
                all_public_keys = get_all_public_keys(config)
                curve_params = {k: v for k, v in config.items() if k == "curve"}
                
                if all_public_keys and self.masking_helper is not None:
                    self.masking_helper.compute_shared_secrets(all_public_keys, curve_params or config)
                elif all_public_keys and self.masking_helper is None:
                    self.logger.warning(f"[Client {self.client_id}] Cannot compute shared secrets: masking_helper is None")
                
                # Compute local QC arrays
                counts_array = compute_genotype_counts(
                    self.plink_prefix, self.client_id, log_dir=self.log_dir
                )
                missing_array = compute_missingness_counts(
                    self.plink_prefix, self.client_id, log_dir=self.log_dir
                )
                
                # Log array shapes for debugging
                self.logger.info(f"[Client {self.client_id}] QC arrays computed: counts shape={counts_array.shape}, missing shape={missing_array.shape}")
                if counts_array.size == 0 or missing_array.size == 0:
                    self.logger.error(f"[Client {self.client_id}] ERROR: Empty QC arrays! counts.size={counts_array.size}, missing.size={missing_array.size}")
                    self.logger.error(f"[Client {self.client_id}] PLINK prefix was: {self.plink_prefix}")
                    # Check if PLINK files exist
                    for ext in ['.bed', '.bim', '.fam']:
                        plink_file = f"{self.plink_prefix}{ext}"
                        exists = os.path.exists(plink_file)
                        self.logger.error(f"[Client {self.client_id}] PLINK file {plink_file} exists: {exists}")

                # Get thresholds from stored instance variables (loaded from config file via DataLoader)
                maf_thresh = getattr(self, 'maf_threshold', 0.01)
                miss_thresh = getattr(self, 'miss_threshold', 0.05)
                hwe_thresh = getattr(self, 'hwe_threshold', 1e-6)
                # Convert string thresholds to float if needed
                if isinstance(maf_thresh, str):
                    maf_thresh = float(maf_thresh)
                if isinstance(miss_thresh, str):
                    miss_thresh = float(miss_thresh)
                if isinstance(hwe_thresh, str):
                    hwe_thresh = float(hwe_thresh)
                threshold_array = np.array([maf_thresh, miss_thresh, hwe_thresh], dtype=np.float64)
                
                # SECURITY (match sync): encrypted broadcast of QC payloads (typed)
                # Each client encrypts its local QC arrays for every recipient.
                # Server forwards ciphertexts but cannot decrypt; each client decrypts all senders and aggregates locally.

                # Align num_clients from server config if provided
                cfg_n = config.get("num_clients")
                if isinstance(cfg_n, int) and cfg_n > 0:
                    self.num_clients = cfg_n

                if not (all_public_keys and self.masking_helper and self.masking_helper.private_key):
                    self.logger.warning(f"[Client {self.client_id}] Missing keys for encrypted broadcast in global_qc; sending plaintext arrays (server can see)")
                    return [counts_array, missing_array, threshold_array], 1, {}

                my_hash_id = resolve_my_hash_id(all_public_keys, getattr(self, "my_public_key_pem", None))
                if my_hash_id is None:
                    self.logger.warning(f"[Client {self.client_id}] Could not resolve my hash_id in global_qc; sending plaintext arrays (server can see)")
                    return [counts_array, missing_array, threshold_array], 1, {}

                messenger = create_client_messenger(
                    my_hash_id,
                    self.num_clients,
                    self.masking_helper.private_key
                )

                public_key_ids = sorted([int(k) for k in all_public_keys.keys()])
                if len(public_key_ids) != self.num_clients:
                    self.logger.warning(f"[Client {self.client_id}] Mismatch: {len(public_key_ids)} public keys but num_clients={self.num_clients}")

                payload_counts_u8 = pack_typed_payload_uint8("counts", counts_array)
                payload_missing_u8 = pack_typed_payload_uint8("missing", missing_array)
                payload_thresh_u8 = pack_typed_payload_uint8("thresholds", threshold_array)

                encrypted_arrays = []
                for recipient_hash_id in public_key_ids:
                    recipient_key_pem = all_public_keys.get(str(recipient_hash_id))
                    if not recipient_key_pem:
                        self.logger.warning(f"[Client {self.client_id}] No public key for recipient {recipient_hash_id}; skipping encrypted QC payloads")
                        continue
                    enc_counts = messenger.encrypt_for_recipient(payload_counts_u8, recipient_hash_id, recipient_key_pem)
                    enc_missing = messenger.encrypt_for_recipient(payload_missing_u8, recipient_hash_id, recipient_key_pem)
                    enc_thresh = messenger.encrypt_for_recipient(payload_thresh_u8, recipient_hash_id, recipient_key_pem)
                    encrypted_arrays.extend([
                        np.frombuffer(enc_counts, dtype=np.uint8),
                        np.frombuffer(enc_missing, dtype=np.uint8),
                        np.frombuffer(enc_thresh, dtype=np.uint8),
                    ])

                self.logger.info(f"[Client {self.client_id}] ✓ Encrypted-broadcasted QC payloads to {len(public_key_ids)} recipients (server cannot decrypt)")
                return encrypted_arrays, 1, {}

            ################################################################################
            # Stage 5: Global QC Response - Receive masked QC data, aggregate locally, compute exclusion list
            ################################################################################
            elif stage == "global_qc_response":
                # Receive encrypted QC payloads from server relay and aggregate locally
                try:
                    from flwr.common import parameters_to_ndarrays
                    all_public_keys = get_all_public_keys(config)
                    cfg_n = config.get("num_clients")
                    if isinstance(cfg_n, int) and cfg_n > 0:
                        self.num_clients = cfg_n

                    if not (all_public_keys and self.masking_helper and self.masking_helper.private_key):
                        self.logger.warning(f"[Client {self.client_id}] Missing keys in global_qc_response; skipping QC aggregation")
                        return [], 1, {}

                    my_hash_id = resolve_my_hash_id(all_public_keys, getattr(self, "my_public_key_pem", None))
                    if my_hash_id is None:
                        self.logger.warning(f"[Client {self.client_id}] Could not resolve my hash_id in global_qc_response; skipping QC aggregation")
                        return [], 1, {}

                    messenger = create_client_messenger(my_hash_id, self.num_clients, self.masking_helper.private_key)

                    ndarrays = parameters_to_ndarrays(parameters) if hasattr(parameters, "tensors") else parameters
                    if not isinstance(ndarrays, list):
                        ndarrays = list(ndarrays)

                    # Decrypt all messages intended for us; aggregate by kind
                    sums: Dict[str, Optional[np.ndarray]] = {"counts": None, "missing": None, "thresholds": None}
                    thresholds_list: List[np.ndarray] = []

                    from pipeline.clients.client_to_client import unpack_envelope
                    for enc_arr in ndarrays:
                        if not isinstance(enc_arr, np.ndarray):
                            continue
                        enc_msg = enc_arr.tobytes()
                        try:
                            meta = unpack_envelope(enc_msg)
                        except Exception:
                            continue
                        if int(meta.get("recipient_id", -1)) != my_hash_id:
                            continue
                        sender_id = int(meta.get("sender_id"))
                        sender_key_pem = all_public_keys.get(str(sender_id))
                        if not sender_key_pem:
                            continue
                        dec_u8 = messenger.decrypt_from_sender(enc_msg, sender_id, sender_key_pem)
                        if not isinstance(dec_u8, np.ndarray) or dec_u8.dtype != np.uint8:
                            continue
                        kind, arr = unpack_typed_payload_uint8(dec_u8)
                        if kind not in sums:
                            continue
                        if kind == "thresholds":
                            thresholds_list.append(arr.astype(np.float64))
                        if sums[kind] is None:
                            sums[kind] = arr.astype(np.float64)
                        else:
                            sums[kind] = sums[kind] + arr.astype(np.float64)

                    if sums["counts"] is None or sums["missing"] is None or not thresholds_list:
                        self.logger.warning(f"[Client {self.client_id}] No decryptable QC payloads found for aggregation")
                        return [], 1, {}

                    thresh_stack = np.stack(thresholds_list, axis=0)  # (n_senders, 3)
                    # Ensure all threshold values are floats (defensive: config values may be strings)
                    maf_final = float(np.min(thresh_stack[:, 0].astype(np.float64)))
                    miss_final = float(np.max(thresh_stack[:, 1].astype(np.float64)))
                    hwe_final = float(np.min(thresh_stack[:, 2].astype(np.float64)))
                    thresholds_final = np.array([maf_final, miss_final, hwe_final], dtype=np.float64)

                    # Debug QC computation (check both environment variable and config)
                    debug_qc_env = os.environ.get('FEDGWAS_DEBUG_QC', '0') == '1'
                    debug_qc_config = config.get('debug_qc', False)
                    debug_qc = debug_qc_env or debug_qc_config
                    exclusion_set = _compute_exclusion_list(sums["counts"], sums["missing"], thresholds_final, debug=debug_qc)
                    exclusion_indices = sorted(list(exclusion_set))
                    
                    # Always log QC statistics when 0 exclusions (helps diagnose issue)
                    if len(exclusion_indices) == 0:
                        self.logger.warning(f"[Client {self.client_id}] QC WARNING: 0 SNPs excluded! Checking aggregated statistics...")
                        self.logger.info(f"[Client {self.client_id}] QC DEBUG: counts_sum shape={sums['counts'].shape}, missing_sum shape={sums['missing'].shape}")
                        self.logger.info(f"[Client {self.client_id}] QC DEBUG: thresholds={thresholds_final}")
                        # Check all SNPs to find which ones should be excluded
                        should_exclude = []
                        for i in range(sums['counts'].shape[0]):
                            N_AA, N_Aa, N_aa = sums['counts'][i]
                            N_obs, N_miss = sums['missing'][i]
                            total_geno = N_AA + N_Aa + N_aa
                            if total_geno == 0:
                                should_exclude.append((i, "zero_geno"))
                                continue
                            freqA = (2 * N_AA + N_Aa) / (2 * total_geno)
                            maf = min(freqA, 1.0 - freqA)
                            if maf < thresholds_final[0]:
                                should_exclude.append((i, f"maf={maf:.6f}"))
                                continue
                            missing_rate = N_miss / (N_obs + N_miss) if (N_obs + N_miss) > 0 else 0
                            if missing_rate > thresholds_final[1]:
                                should_exclude.append((i, f"missing={missing_rate:.6f}"))
                                continue
                            # Check HWE
                            p = freqA
                            q = 1.0 - p
                            E_AA = p * p * total_geno
                            E_Aa = 2.0 * p * q * total_geno
                            E_aa = q * q * total_geno
                            chi_sq = 0.0
                            for obs, exp in [(N_AA, E_AA), (N_Aa, E_Aa), (N_aa, E_aa)]:
                                if exp > 1e-9:
                                    chi_sq += (obs - exp) ** 2 / exp
                            from scipy.stats import chi2
                            pval_hwe = 1 - chi2.cdf(chi_sq, df=1)
                            if pval_hwe < thresholds_final[2]:
                                should_exclude.append((i, f"hwe={pval_hwe:.2e}"))
                        
                        if should_exclude:
                            self.logger.warning(f"[Client {self.client_id}] QC DEBUG: Found {len(should_exclude)} SNPs that SHOULD be excluded but weren't!")
                            for idx, reason in should_exclude[:10]:  # Show first 10
                                self.logger.warning(f"[Client {self.client_id}] QC DEBUG: SNP {idx} should be excluded: {reason}")
                        else:
                            self.logger.info(f"[Client {self.client_id}] QC DEBUG: All SNPs pass QC thresholds (this is unexpected if baseline excluded SNPs)")
                        
                        # Sample a few SNPs to show their values
                        for i in range(min(5, sums['counts'].shape[0])):
                            N_AA, N_Aa, N_aa = sums['counts'][i]
                            N_obs, N_miss = sums['missing'][i]
                            total_geno = N_AA + N_Aa + N_aa
                            if total_geno > 0:
                                freqA = (2 * N_AA + N_Aa) / (2 * total_geno)
                                maf = min(freqA, 1.0 - freqA)
                                missing_rate = N_miss / (N_obs + N_miss) if (N_obs + N_miss) > 0 else 0
                                self.logger.info(f"[Client {self.client_id}] QC DEBUG: SNP {i} sample: MAF={maf:.6f}, Missing={missing_rate:.6f}, Geno=[{N_AA},{N_Aa},{N_aa}], Obs/Miss=[{N_obs},{N_miss}]")
                    
                    self.logger.info(f"[Client {self.client_id}] Computed global exclusion list locally: {len(exclusion_indices)} SNPs (server cannot decrypt)")

                    # Map indices -> SNP IDs and exclude locally
                    bim_file = self.plink_prefix + ".bim"
                    snp_ids: List[str] = []
                    if os.path.exists(bim_file):
                        with open(bim_file, "r") as f:
                            for line in f:
                                parts = line.strip().split()
                                if len(parts) >= 2:
                                    snp_ids.append(parts[1])
                    excluded_snp_ids = [snp_ids[i] for i in exclusion_indices if 0 <= i < len(snp_ids)]
                    if excluded_snp_ids:
                        # Save exclusion list to file for result analysis
                        exclusion_file = os.path.join(self.log_dir, "global_qc_exclusion.txt")
                        try:
                            with open(exclusion_file, 'w') as f:
                                for snp_id in sorted(excluded_snp_ids):
                                    f.write(f"{snp_id}\n")
                            self.logger.info(f"[Client {self.client_id}] Saved {len(excluded_snp_ids)} excluded SNPs to {exclusion_file}")
                        except Exception as e:
                            self.logger.warning(f"[Client {self.client_id}] Failed to save exclusion list: {e}")
                        
                        new_prefix = exclude_snps(self.plink_prefix, excluded_snp_ids, "global_filtered", log_dir=self.log_dir)
                        self.plink_prefix = new_prefix  # Update the instance variable
                        self._persist_plink_prefix(new_prefix)
                        self.logger.info(f"[Client {self.client_id}] Global QC filter => new prefix {self.plink_prefix} (persisted to state)")
                    else:
                        self.logger.info(f"[Client {self.client_id}] No SNPs excluded from global QC")
                    return [], 1, {}
                except Exception as e:
                    self.logger.warning(f"[Client {self.client_id}] Failed to process encrypted global_qc_response: {e}")
                    import traceback
                    self.logger.warning(traceback.format_exc())
                    return [], 1, {}

            ################################################################################
            # Stage 6: Init Chunks (KING) - Use locally computed global seed for data chunking
            ################################################################################
            elif stage == "init_chunks":
                # Restore plink_prefix from state (may have been updated by QC filtering)
                self._restore_plink_prefix()
                
                # KING should chunk by SNPs to accumulate the same pair across SNP subsets
                self.partition_by = "snps"

                # Use global_seed computed locally in sync stage (server doesn't know it)
                if not hasattr(self, 'global_seed') or self.global_seed is None:
                    self.logger.warning(f"[Client {self.client_id}] Global seed not computed, using local seed as fallback")
                    self.global_seed = int(self.local_seed % (10**9))
                
                # Pass global seed to partition_data for deterministic partitioning
                partition_config = dict(config) if hasattr(config, '__iter__') else {}
                partition_config['global_seed'] = self.global_seed
                # Prefer SNP chunk size for KING; fall back to chunk_size
                if hasattr(self, "snp_chunk_size") and self.snp_chunk_size is not None:
                    partition_config["chunk_size"] = self.snp_chunk_size
                elif "snp_chunk_size" in partition_config:
                    partition_config["chunk_size"] = partition_config.get("snp_chunk_size")
                self.partition_data(partition_config)
                self.current_chunk_idx = 0
                # Persist the chunking round so iterative KING can recreate identical chunks
                try:
                    from flwr.common import ConfigRecord
                    chunk_round = int(config.get("server_round", 0) or 0)
                    paths_rec = self.client_state.config_records.get("client_paths")
                    paths_dict = dict(paths_rec) if paths_rec else {}
                    paths_dict["king_chunk_round"] = chunk_round
                    paths_dict["king_chunk_idx"] = 0
                    self.client_state.config_records["client_paths"] = ConfigRecord(paths_dict)
                except Exception:
                    pass
                self.logger.info(f"[Client {self.client_id}] Created {len(self.chunk_data)} chunks for iterative KING using locally-computed global_seed={self.global_seed}.")
                return [], 1, {}

            ################################################################################
            # Stage 7: Iterative KING - Process KING relationship analysis chunk by chunk
            ################################################################################
            elif stage == "iterative_king":
                # Restore plink_prefix from state (may have been updated by QC filtering)
                self._restore_plink_prefix()
                
                self.logger.info(f"[Client {self.client_id}] === ITERATIVE_KING STAGE CALLED ===")
                self.logger.info(f"[Client {self.client_id}] Chunk data length: {len(self.chunk_data) if hasattr(self, 'chunk_data') else 'no chunk_data'}")
                self.logger.info(f"[Client {self.client_id}] Calling handle_iterative_king function...")
                result = handle_iterative_king(self, parameters, config)
                self.logger.info(f"[Client {self.client_id}] handle_iterative_king returned: {result}")
                return result

            ################################################################################
            # Stage 8: Local LR - Run local logistic regression and identify insignificant SNPs
            ################################################################################
            elif stage == "local_lr":
                # Restore plink_prefix from state (may have been updated by QC or KING filtering)
                self._restore_plink_prefix()
                
                self.logger.info(f"[Client {self.client_id}] Stage: local_lr")
                # Use local_lr_threshold from client config (loaded from config file), not from stage config
                # Ensure local_lr_threshold is always a float (config values may be strings)
                local_lr_threshold_raw = getattr(self, "local_lr_threshold", config.get("local_lr_threshold", 1e-3))
                local_lr_threshold = float(local_lr_threshold_raw) if local_lr_threshold_raw is not None else 1e-3
                self.logger.info(f"[Client {self.client_id}] Using local_lr_threshold={local_lr_threshold} for local LR filtering")
                assoc_file = run_local_lr(self.plink_prefix, out_prefix="local_lr_temp", log_dir=self.log_dir)
                insign_snps = parse_insignificant_snps(assoc_file, p_threshold=local_lr_threshold)
                self.logger.info(
                    f"[Client {self.client_id}] Found {len(insign_snps)} insignificant SNPs "
                    f"(p >= {local_lr_threshold})"
                )

                # Send tokenized identifiers to server so it cannot learn plaintext SNP IDs
                lr_tokenized = bool(int(config.get("lr_tokenized", 1))) if str(config.get("lr_tokenized", 1)).isdigit() else True
                if lr_tokenized:
                    # Ensure global_seed exists (should from sync/sync_response); fallback deterministic if missing
                    base_dir = get_persist_base_dir(getattr(self, "log_dir", None), getattr(self, "intermediate_dir", None))
                    all_public_keys = get_all_public_keys(config)
                    self.global_seed = ensure_global_seed(
                        getattr(self, "global_seed", None),
                        base_dir,
                        all_public_keys,
                        int(config.get("server_round", 0) or 0),
                        self.local_seed,
                    )
                    persist_seed(base_dir, int(self.global_seed))
                    token_salt = str(config.get("lr_token_salt", f"r{int(config.get('server_round', 0) or 0)}"))
                    # Store token salt for later resolution (persist to state so it survives client restarts)
                    self.lr_token_salt = token_salt
                    # Also persist to disk in case actor state is lost between rounds
                    try:
                        with open(os.path.join(base_dir, "lr_token_salt.txt"), "w") as f:
                            f.write(token_salt)
                    except Exception:
                        pass
                    # Persist to state
                    if self.client_state:
                        try:
                            paths_rec = self.client_state.config_records.get("client_paths")
                            if paths_rec:
                                paths_dict = dict(paths_rec)
                            else:
                                paths_dict = {}
                            paths_dict["lr_token_salt"] = token_salt
                            self.client_state.config_records["client_paths"] = ConfigRecord(paths_dict)
                        except Exception as e:
                            self.logger.warning(f"[Client {self.client_id}] Failed to persist lr_token_salt: {e}")
                    self.logger.info(f"[Client {self.client_id}] Stored lr_token_salt={token_salt} for later resolution")
                    tokens = [lr_token(getattr(self, "global_seed", None), self.local_seed, s, token_salt=token_salt) for s in insign_snps]
                    # Optional padding: does NOT prevent leakage to an honest-but-curious server that computes intersection,
                    # but can prevent size leakage to observers/logs if server forwards without inspecting.
                    lr_pad_to = int(config.get("lr_pad_to", 0) or 0)
                    if lr_pad_to > 0 and len(tokens) < lr_pad_to:
                        tokens.extend(make_dummy_tokens(getattr(self, "global_seed", None), self.local_seed, lr_pad_to - len(tokens), token_salt=token_salt))
                    joined = "\n".join(tokens)
                    self.logger.info(f"[Client {self.client_id}] Sending {len(tokens)} tokenized insign. SNPs to server ")
                else:
                    joined = "\n".join(insign_snps)
                    self.logger.warning(f"[Client {self.client_id}] Sending plaintext insign. SNP IDs (lr_tokenized=0)")

                data_bytes = joined.encode("utf-8")
                data_array = np.frombuffer(data_bytes, dtype=np.uint8)
                self.logger.info(f"[Client {self.client_id}] Found {len(insign_snps)} insign. SNPs locally (already filtered)")
                return [data_array], 1, {}

            ################################################################################
            # Stage 9: Local LR Filter Response - Receive and apply intersection of insignificant SNPs
            ################################################################################
            elif stage == "local_lr_filter_response":
                # Restore latest filtered prefix (should be KING-filtered before LR filtering)
                try:
                    self._restore_plink_prefix()
                    if not os.path.isabs(self.plink_prefix):
                        self.plink_prefix = os.path.abspath(self.plink_prefix)
                    self.logger.info(f"[Client {self.client_id}] Local LR filter response using plink_prefix={self.plink_prefix}")
                except Exception as e:
                    self.logger.warning(f"[Client {self.client_id}] Failed to restore plink_prefix before LR filter response: {e}")

                try:
                    from flwr.common import parameters_to_ndarrays
                    ndarrays = parameters_to_ndarrays(parameters) if hasattr(parameters, "tensors") else parameters
                    if not isinstance(ndarrays, list):
                        ndarrays = list(ndarrays)

                    if not ndarrays:
                        self.logger.warning(f"[Client {self.client_id}] Local LR filter response: no token arrays received")
                        return [], 1, {}

                    # Decode all token sets
                    token_sets = []
                    for arr in ndarrays:
                        token_str = arr.tobytes().decode("utf-8").strip()
                        toks = [tok for tok in token_str.split() if tok]
                        token_sets.append(set(toks))
                        self.logger.info(f"[Client {self.client_id}] Received token set size={len(toks)} from server relay")

                    # Compute intersection locally (secure intersection)
                    intersection = token_sets[0].copy()
                    for s in token_sets[1:]:
                        intersection &= s
                    self.logger.info(f"[Client {self.client_id}] Local LR: intersection size={len(intersection)}")

                    lr_tokenized = bool(int(config.get("lr_tokenized", 1))) if str(config.get("lr_tokenized", 1)).isdigit() else True
                    if lr_tokenized:
                        # Use the token salt from when tokens were generated (stored in local_lr stage)
                        # Try to restore from state first, then instance variable, then fallback
                        token_salt = None
                        if self.client_state:
                            try:
                                paths_rec = self.client_state.config_records.get("client_paths")
                                if paths_rec:
                                    paths_dict = dict(paths_rec)
                                    token_salt = paths_dict.get("lr_token_salt")
                            except Exception as e:
                                self.logger.warning(f"[Client {self.client_id}] Failed to restore lr_token_salt from state: {e}")
                        if not token_salt:
                            token_salt = getattr(self, "lr_token_salt", None)
                        if not token_salt:
                            # Try to restore from disk if state is missing (actor restart)
                            try:
                                base_dir = get_persist_base_dir(getattr(self, "log_dir", None), getattr(self, "intermediate_dir", None))
                                salt_path = os.path.join(base_dir, "lr_token_salt.txt")
                                if os.path.exists(salt_path):
                                    with open(salt_path, "r") as f:
                                        token_salt = f.read().strip()
                            except Exception:
                                pass
                        if not token_salt:
                            # Fallback: use previous round's salt (server_round - 1)
                            token_salt = str(config.get("lr_token_salt", f"r{int(config.get('server_round', 0) or 0) - 1}"))
                            self.logger.warning(f"[Client {self.client_id}] Using fallback token_salt={token_salt} (stored value not found)")
                        else:
                            self.logger.info(f"[Client {self.client_id}] Using stored token_salt={token_salt} for resolution (tokens were generated with this salt)")
                        
                        # Debug: log what we're trying to resolve
                        self.logger.info(f"[Client {self.client_id}] Resolving {len(intersection)} tokens with global_seed={getattr(self, 'global_seed', None)}, local_seed={self.local_seed}, token_salt={token_salt}")
                        
                        snps_to_exclude = resolve_tokens_to_snp_ids(
                            self.plink_prefix,
                            list(intersection),
                            getattr(self, "global_seed", None),
                            self.local_seed,
                            token_salt=token_salt,
                        )
                        self.logger.info(f"[Client {self.client_id}] Resolved {len(snps_to_exclude)} SNP IDs for exclusion locally (server never saw SNP IDs)")
                        
                        if len(snps_to_exclude) == 0 and len(intersection) > 0:
                            self.logger.warning(f"[Client {self.client_id}] WARNING: {len(intersection)} tokens in intersection but 0 SNP IDs resolved! This suggests a token salt mismatch or token generation issue.")
                        if snps_to_exclude:
                            new_prefix = exclude_snps(self.plink_prefix, snps_to_exclude, "lr_filtered", log_dir=self.log_dir)
                            self.plink_prefix = new_prefix
                            self._persist_plink_prefix(new_prefix)
                            self.logger.info(f"[Client {self.client_id}] Local LR filter => new prefix {self.plink_prefix} (persisted to state)")
                        else:
                            self.logger.info(f"[Client {self.client_id}] No SNPs to exclude from local LR (all significant or no intersection)")
                    else:
                        intersection_snps = list(intersection)
                        self.logger.info(f"[Client {self.client_id}] Received {len(intersection_snps)} SNPs from server relay (plaintext)")
                        if intersection_snps:
                            new_prefix = exclude_snps(self.plink_prefix, intersection_snps, "lr_filtered", log_dir=self.log_dir)
                            self.plink_prefix = new_prefix
                except Exception as e:
                    self.logger.warning(f"[Client {self.client_id}] Local LR filter response failed: {e}")
                return [], 1, {}

            ################################################################################
            # Stage 10: Init Chunks (LR) - Partition filtered data for final LR analysis
            ################################################################################
            elif stage == "init_chunks_lr":
                # Restore plink_prefix from state (should be after QC + KING filtering)
                # This ensures we use the dataset that has been filtered by:
                # 1. Global QC (SNPs excluded)
                # 2. Iterative KING (relatives/samples excluded)
                old_prefix = self.plink_prefix
                self._restore_plink_prefix()
                if self.plink_prefix != old_prefix:
                    self.logger.info(f"[Client {self.client_id}] Restored plink_prefix: {old_prefix} -> {self.plink_prefix}")
                else:
                    self.logger.warning(f"[Client {self.client_id}] plink_prefix unchanged after restore: {self.plink_prefix}")
                
                # Verify the prefix exists and log its properties
                if os.path.exists(f"{self.plink_prefix}.bim"):
                    import subprocess
                    try:
                        # Count variants in the filtered dataset
                        result = subprocess.run(
                            ["wc", "-l", f"{self.plink_prefix}.bim"],
                            capture_output=True, text=True, check=True
                        )
                        variant_count = int(result.stdout.split()[0])
                        self.logger.info(f"[Client {self.client_id}] Filtered dataset has {variant_count} variants")
                    except Exception as e:
                        self.logger.warning(f"[Client {self.client_id}] Could not count variants: {e}")
                else:
                    self.logger.error(f"[Client {self.client_id}] Filtered prefix {self.plink_prefix} does not exist!")
                
                # Clear existing chunk_data to force recreation with the correct (filtered) prefix
                # This is important because chunks from KING stage might have been created with the old prefix
                if hasattr(self, 'chunk_data'):
                    self.logger.info(f"[Client {self.client_id}] Clearing {len(self.chunk_data)} existing chunks to recreate with filtered prefix")
                    self.chunk_data = []
                # LR should chunk by samples
                self.partition_by = "samples"
                
                # Pass global seed to partition_data for deterministic partitioning
                partition_config = dict(config) if hasattr(config, '__iter__') else {}
                # Never silently fall back to 0: ensure global_seed is set (persisted or deterministic fallback)
                base_dir = get_persist_base_dir(getattr(self, "log_dir", None), getattr(self, "intermediate_dir", None))
                all_public_keys = get_all_public_keys(partition_config)
                self.global_seed = ensure_global_seed(
                    getattr(self, "global_seed", None),
                    base_dir,
                    all_public_keys,
                    int(partition_config.get("server_round", 0) or 0),
                    self.local_seed,
                )
                persist_seed(base_dir, int(self.global_seed))
                partition_config['global_seed'] = int(self.global_seed)

                # Optionally compute chunk size from a target number of LR chunks
                chunk_size, sample_count, target_chunks, used_target = compute_lr_chunk_size(self, partition_config)
                partition_config["chunk_size"] = int(chunk_size)
                if used_target:
                    num_chunks = int(math.ceil(sample_count / float(chunk_size))) if sample_count else "unknown"
                    self.logger.info(
                        f"[Client {self.client_id}] LR target chunks={target_chunks} "
                        f"(samples={sample_count}, chunk_size={chunk_size}, chunks={num_chunks})"
                    )
                
                # Log the prefix that will be used for partitioning
                self.logger.info(f"[Client {self.client_id}] Partitioning with plink_prefix={self.plink_prefix} for iterative LR")
                self.partition_data(partition_config)
                self.current_chunk_idx = 0
                # Persist the chunking round so iterative LR can recreate identical chunks
                try:
                    from flwr.common import ConfigRecord
                    chunk_round = int(partition_config.get("server_round", 0) or 0)
                    paths_rec = self.client_state.config_records.get("client_paths")
                    paths_dict = dict(paths_rec) if paths_rec else {}
                    paths_dict["lr_chunk_round"] = chunk_round
                    paths_dict["lr_chunk_idx"] = 0
                    self.client_state.config_records["client_paths"] = ConfigRecord(paths_dict)
                except Exception:
                    pass
                self.logger.info(f"[Client {self.client_id}] Created {len(self.chunk_data)} chunks for iterative LR.")
                return [], 1, {}

            ################################################################################
            # Stage 11: Iterative LR - Process final logistic regression analysis chunk by chunk
            ################################################################################
            elif stage == "iterative_lr":
                # Ensure we are using the latest filtered prefix (after KING + local LR)
                try:
                    self._restore_plink_prefix()
                    if not os.path.isabs(self.plink_prefix):
                        self.plink_prefix = os.path.abspath(self.plink_prefix)
                    self.logger.info(f"[Client {self.client_id}] Iterative LR using plink_prefix={self.plink_prefix}")
                except Exception as e:
                    self.logger.warning(f"[Client {self.client_id}] Failed to restore plink_prefix before iterative LR: {e}")

                return handle_iterative_lr(self, parameters, config)

            # Server sends final LR results for de-anonymization/labeling
            elif stage == "iterative_lr_response":
                return handle_iterative_lr(self, parameters, config)

            else:
                # default fallback
                return [], 1, {}
    
    def __del__(self):
        if getattr(self, "_perf_tracker", None) is not None:
            try:
                self._perf_tracker.finalize()
            except Exception:
                pass
        # clean up logger when client is destroyed
        LoggerManager.close_logger(self.client_id)

def client_fn(context: Context):
    """Client function that creates and returns a FedLRClient instance."""
    
    # Get partition_id from context - handle different modes
    if hasattr(context, 'partition_id'):
        partition_id = context.partition_id
    elif hasattr(context, 'node_config') and 'partition-id' in context.node_config:
        partition_id = context.node_config['partition-id']
    else:
        # In simulation mode, we might not have partition_id, so use a default
        # partition_id = 0
        # print(f"[Client] Warning: No partition_id found in context, using default: {partition_id}")
        raise ValueError("No partition_id found!")

    # Priority: node_config (from SuperNode) > run_config (from simulation)
    # When using local-deployment with SuperNodes, node_config should have the correct config-file
    # Check node_config first, regardless of simulation mode
    config_file_path = None
    simulation = None  # Will be set if needed
    
    if hasattr(context, "node_config") and isinstance(context.node_config, dict) and "config-file" in context.node_config:
        # Local deployment mode: use config-file from SuperNode node_config
        config_file_path = context.node_config["config-file"]
        # Resolve to absolute path if relative
        import pathlib
        config_path_obj = pathlib.Path(config_file_path)
        if not config_path_obj.is_absolute():
            repo_root = pathlib.Path.cwd()
            for parent in [repo_root, *repo_root.parents]:
                if (parent / "pyproject.toml").is_file():
                    repo_root = parent
                    break
            config_path_obj = repo_root / config_path_obj
        config_file_path = str(config_path_obj.resolve())
        simulation = False  # Using node_config means deployment mode
    else:
        # Fallback to simulation mode logic if node_config not available
        if "simulation" in context.run_config:
            simulation = context.run_config["simulation"]
        else:
            # default to simulation mode
            simulation = True

        if simulation:
            # simulation mode use config_path from run_config
            if "config_path" in context.run_config:
                config_path = context.run_config["config_path"]
                config_file_path = str(
                    _resolve_simulation_client_config_path(
                        config_path,
                        partition_id + 1,
                    )
                )
            else:
                raise ValueError("Need to specify the config path in the run_config")
        else:
            # deployment mode fallback: use default config location
            config_file_path = "configs/config.yaml"
    
    # Verify config file exists
    import pathlib
    config_path_check = pathlib.Path(config_file_path)
    if not config_path_check.exists():
        raise FileNotFoundError(f"Client config file not found: {config_file_path} (cwd: {pathlib.Path.cwd()})")
    
    
    # Create and return the client, passing the context
    return FedLRClient(
        partition_id=partition_id+1,
        context=context,
        config_file=config_file_path,
        partition_by="samples"
    ).to_client()

app = ClientApp(client_fn)
    
    
