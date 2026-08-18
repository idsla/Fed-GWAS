from __future__ import annotations

import json
import logging
import os
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
from flwr.common import (
    FitIns,
    FitRes,
    Parameters,
    Scalar,
    ndarrays_to_parameters,
    parameters_to_ndarrays,
)
from flwr.server.client_manager import ClientManager
from flwr.server.client_proxy import ClientProxy
from flwr.server.strategy import Strategy

from pipeline.server.aggregator_king import run_server_king
from pipeline.server.aggregator_lr import run_server_lr
from pipeline.server.prg_masking import create_prg_masking_aggregator
from pipeline.utils.performance.monitoring_runtime import (
    ServerPerformanceTracker,
    merge_performance_csvs_to_results_root,
    resolve_results_root_from_server_output,
)
from pipeline.utils.run_retention import maybe_apply_retention_on_complete


class FederatedGWASStrategy(Strategy):
    """
    Strict, simulation-friendly strategy which:
    - Keeps stage config Scalar-only (JSON strings, ints)
    - Forwards encrypted client-to-client blobs via Parameters (server cannot decrypt)
    - Does NOT use the old "advance stages on failures" workaround
    """

    def __init__(
        self,
        output_dir: Optional[str] = None,
        chunk_size: int = 1000,
        lr_pad_to: int = 0,
        monitoring_settings: Optional[Dict] = None,
        retention_settings: Optional[Dict] = None,
    ):
        self.monitoring_settings = monitoring_settings or {}
        self.retention_settings = retention_settings or {}
        self._server_tracker: Optional[ServerPerformanceTracker] = None
        self._monitoring_finalized = False
        self._retention_applied = False
        self._server_output_base = output_dir

        # Set up file logging if output_dir is provided
        if output_dir:
            self.output_dir = os.path.join(output_dir, "intermediate")
            log_dir = os.path.join(output_dir, "logs")
            os.makedirs(log_dir, exist_ok=True)
            log_file = os.path.join(log_dir, "server_log.txt")
            
            # Configure parent logger for all server modules (pipeline.server.*)
            server_logger = logging.getLogger("pipeline.server")
            server_logger.setLevel(logging.INFO)
            
            # Remove existing handlers to avoid duplicates
            server_logger.handlers = []
            
            # Add file handler
            file_handler = logging.FileHandler(log_file, mode='w')
            file_formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
            file_handler.setFormatter(file_formatter)
            server_logger.addHandler(file_handler)
            
            # Also keep console handler for immediate feedback
            console_handler = logging.StreamHandler()
            console_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
            console_handler.setFormatter(console_formatter)
            server_logger.addHandler(console_handler)
            
            # Prevent propagation to root logger to avoid duplicate console output
            # But allow child loggers to inherit handlers (they will propagate to this logger)
            server_logger.propagate = False
            
            # Get the strategy-specific logger (child of pipeline.server)
            self.logger = logging.getLogger(__name__)
            # Ensure it propagates to parent (pipeline.server) so it uses the file handler
            self.logger.propagate = True
            self.logger.setLevel(logging.INFO)
            self.logger.info(f"[Server] Logging initialized. Server log file: {log_file}")
            if self.monitoring_settings.get("enable_performance_monitoring"):
                experiment_name = self.monitoring_settings.get("experiment_name", "fedgwas")
                self._server_tracker = ServerPerformanceTracker(
                    log_dir, self.monitoring_settings, experiment_name
                )
                self.logger.info("[Server] Performance monitoring enabled (output: %s)", log_dir)
        else:
            # Fallback: only console logging if no output_dir
            self.logger = logging.getLogger(__name__)
            self.logger.setLevel(logging.INFO)
            if not self.logger.handlers:
                handler = logging.StreamHandler()
                formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
                handler.setFormatter(formatter)
                self.logger.addHandler(handler)
            self.output_dir = "./server_intermediate"
        
        os.makedirs(self.output_dir, exist_ok=True)

        self.current_stage = "key_exchange"
        try:
            self.chunk_size = max(1, int(chunk_size))
        except (TypeError, ValueError):
            self.chunk_size = 1000
        try:
            self.lr_pad_to = max(0, int(lr_pad_to))
        except (TypeError, ValueError):
            self.lr_pad_to = 0

        self.connected_clients: set[str] = set()
        self.num_clients: Optional[int] = None

        self.prg_aggregator = create_prg_masking_aggregator(self.num_clients)

        # transient storage for forwarding ciphertext blobs
        self._pending_forward_params: Optional[Parameters] = None
        
        # Store final LR results to send when all clients are done
        self._final_lr_results: Optional[bytes] = None
        self._final_lr_meta: Optional[np.ndarray] = None
        
        # Track if termination has been logged to avoid repeated messages
        self._termination_logged = False

    def _finalize_monitoring(self) -> None:
        if self._monitoring_finalized or self._server_tracker is None or not self._server_tracker.enabled:
            return
        self._server_tracker.finalize()
        if self._server_output_base:
            results_root = resolve_results_root_from_server_output(self._server_output_base)
            merge_performance_csvs_to_results_root(results_root)
            self.logger.info("[Server] Merged performance CSVs into %s", results_root)
        self._monitoring_finalized = True

    def _apply_retention(self) -> None:
        if self._retention_applied or not self._server_output_base:
            return
        results_root = resolve_results_root_from_server_output(self._server_output_base)
        maybe_apply_retention_on_complete(
            results_root,
            self.retention_settings,
            logger=self.logger,
        )
        self._retention_applied = True

    # -------------------------
    # Strategy API
    # -------------------------
    def initialize_parameters(self, client_manager: ClientManager) -> Optional[Parameters]:
        return Parameters(tensors=[], tensor_type="")

    def configure_fit(
        self, server_round: int, parameters: Parameters, client_manager: ClientManager
    ) -> List[Tuple[ClientProxy, FitIns]]:
        # Graceful termination: if stage is "done", stop processing
        if self.current_stage == "done":
            if not self._termination_logged:
                self.logger.info("[Server] Pipeline complete. Gracefully terminating server.")
                self._termination_logged = True
            return []  # Return empty list to stop further rounds
        
        num_available = client_manager.num_available()
        clients = client_manager.sample(num_clients=num_available, min_num_clients=num_available)

        config = self._get_stage_config(server_round)

        # Send forwarded blobs (if any) as Parameters to the next stage
        params_to_send = self._pending_forward_params if self._pending_forward_params is not None else parameters
        self._pending_forward_params = None

        fit_ins = FitIns(parameters=params_to_send, config=config)
        return [(client, fit_ins) for client in clients]

    def aggregate_fit(
        self,
        server_round: int,
        results: List[Tuple[ClientProxy, FitRes]],
        failures: List[Union[Tuple[ClientProxy, FitRes], BaseException]],
    ) -> Tuple[Optional[Parameters], Dict[str, Scalar]]:
        # Graceful termination: if stage is "done", return early
        if self.current_stage == "done":
            self.logger.info("[Server] Stage is 'done'. No further aggregation needed.")
            self._finalize_monitoring()
            self._apply_retention()
            return None, {}
        
        # track connected clients
        for client_proxy, _ in results:
            self.connected_clients.add(str(getattr(client_proxy, "cid", str(client_proxy))))
        if self.num_clients is None and self.connected_clients:
            self.num_clients = len(self.connected_clients)
            if hasattr(self.prg_aggregator, "num_clients"):
                self.prg_aggregator.num_clients = self.num_clients

        if self._server_tracker is not None and self._server_tracker.enabled:
            self._server_tracker.begin_stage(
                self.current_stage, len(results), len(failures)
            )

        # Key exchange: collect ECC public keys
        if self.current_stage == "key_exchange":
            for client_proxy, fit_res in results:
                if not fit_res.parameters:
                    continue
                nds = parameters_to_ndarrays(fit_res.parameters)
                if not nds:
                    continue
                public_key_str = nds[0].tobytes().decode("utf-8")
                cid = getattr(client_proxy, "cid", str(client_proxy))
                try:
                    cid_int = int(cid)
                except Exception:
                    cid_int = abs(hash(cid)) % (10**18)
                self.prg_aggregator.add_client_public_key(cid_int, public_key_str)

            if self.prg_aggregator.is_key_exchange_complete():
                self.current_stage = "sync"
            return None, {}

        # Sync: forward encrypted seed messages (server cannot decrypt)
        if self.current_stage == "sync":
            if len(results) == 0:
                # retry sync next round; do NOT advance stages on failures
                self.logger.warning(f"[Server] Sync stage: no results (failures={len(failures)}), retrying")
                return None, {"stage": "sync"}

            forwarded: List[np.ndarray] = []
            for _, fit_res in results:
                if not fit_res.parameters:
                    continue
                forwarded.extend(parameters_to_ndarrays(fit_res.parameters))

            self.logger.info(
                f"[Server] Sync stage: received {len(forwarded)} encrypted seed messages; forwarding via Parameters (server cannot decrypt)"
            )
            self._pending_forward_params = ndarrays_to_parameters(forwarded)
            self.current_stage = "sync_response"
            return None, {}

        # Sync response: clients compute global seed locally
        if self.current_stage == "sync_response":
            self.logger.info("[Server] Sync Response stage: clients computing global seed (from encrypted shares or deterministically)")
            self.current_stage = "local_qc"
            return None, {}

        if self.current_stage == "local_qc":
            self.current_stage = "global_qc"
            return None, {}

        # Global QC: forward encrypted QC blobs
        if self.current_stage == "global_qc":
            if len(results) == 0:
                self.logger.warning(f"[Server] Global QC stage: no results (failures={len(failures)}), retrying")
                return None, {"stage": "global_qc"}

            forwarded: List[np.ndarray] = []
            for _, fit_res in results:
                if not fit_res.parameters:
                    continue
                forwarded.extend(parameters_to_ndarrays(fit_res.parameters))

            self.logger.info(
                f"[Server] Global QC stage: received {len(forwarded)} encrypted QC messages (server cannot decrypt); forwarding"
            )
            self._pending_forward_params = ndarrays_to_parameters(forwarded)
            self.current_stage = "global_qc_response"
            return None, {}

        if self.current_stage == "global_qc_response":
            # Note: Server does not have access to exclusion list (clients compute it locally from encrypted data)
            # Exclusion lists are saved by clients in their log directories
            # Result analyzer will aggregate them from client logs
            self.current_stage = "init_chunks"
            return None, {}

        if self.current_stage == "init_chunks":
            self.current_stage = "iterative_king"
            return Parameters(tensors=[], tensor_type=""), {}

        if self.current_stage == "iterative_king":
            def _all_clients_done(metric_key: str) -> bool:
                if not results:
                    return False
                done_flags = []
                for _, fit_res in results:
                    done_flags.append(bool(getattr(fit_res, "metrics", {}).get(metric_key, False)))
                return all(done_flags)

            any_params = False
            forwarded_maps: List[np.ndarray] = []
            for _, fit_res in results:
                if not fit_res.parameters:
                    continue
                nds = parameters_to_ndarrays(fit_res.parameters)
                if nds and len(nds) > 0:
                    any_params = True
                if nds and len(nds) > 1:
                    forwarded_maps.extend(nds[1:])

            all_done = _all_clients_done("king_done")

            if any_params:
                try:
                    # Pass full results for deterministic sorting by (client_id, chunk_index)
                    king_results = run_server_king(self, results=results, output_dir=self.output_dir)
                    if king_results is not None and len(king_results) > 0:
                        self.logger.info(f"[Server] Sending KING results back to clients ({len(king_results)} bytes)")
                        payloads = [king_results] + forwarded_maps
                        return ndarrays_to_parameters(payloads), {}
                except Exception as e:
                    self.logger.error(f"[Server] KING analysis failed: {e}")
                # Stay in iterative_king to continue draining chunks
                return None, {}

            if all_done:
                self.logger.info("[Server] All clients reported KING chunks done; moving to local_lr stage")
                self.current_stage = "local_lr"
                return None, {}

            self.logger.info("[Server] No KING chunks received; retrying iterative_king stage")
            return None, {}

        if self.current_stage == "local_lr":
            # Router behavior: just forward all token sets to clients; clients perform secure intersection locally.
            token_arrays: List[np.ndarray] = []
            for client_proxy, fit_res in results:
                if not fit_res.parameters:
                    self.logger.warning(f"[Server] Local LR: Client {client_proxy.cid} sent no parameters")
                    continue
                nds = parameters_to_ndarrays(fit_res.parameters)
                if not nds:
                    self.logger.warning(f"[Server] Local LR: Client {client_proxy.cid} sent empty ndarrays")
                    continue
                token_arrays.extend(nds)
                # Log size of first ndarray (token list)
                try:
                    token_str = nds[0].tobytes().decode("utf-8").strip()
                    toks = [tok for tok in token_str.split() if tok]
                    self.logger.info(f"[Server] Local LR: Client {client_proxy.cid} sent {len(toks)} tokenized insignificant SNPs")
                except Exception as e:
                    self.logger.warning(f"[Server] Local LR: Failed to decode tokens from client {client_proxy.cid}: {e}")

            if not token_arrays:
                self.logger.warning("[Server] Local LR: No valid token sets received from clients")
                self.current_stage = "init_chunks_lr"
                return None, {}

            # Forward all token arrays back; clients will compute intersection themselves
            self.current_stage = "local_lr_filter_response"
            self.logger.info(f"[Server] Local LR: Forwarding {len(token_arrays)} token arrays to clients for local intersection")
            return ndarrays_to_parameters(token_arrays), {}

        if self.current_stage == "local_lr_filter_response":
            self.current_stage = "init_chunks_lr"
            return None, {}

        if self.current_stage == "init_chunks_lr":
            self.current_stage = "iterative_lr"
            return Parameters(tensors=[], tensor_type=""), {}

        if self.current_stage == "iterative_lr":
            all_params = [fit_res.parameters for _, fit_res in results if fit_res.parameters]
            def _all_clients_done(metric_key: str) -> bool:
                if not results:
                    return False
                done_flags = []
                for _, fit_res in results:
                    done_flags.append(bool(getattr(fit_res, "metrics", {}).get(metric_key, False)))
                return all(done_flags)

            any_params = False
            for _, fit_res in results:
                if not fit_res.parameters:
                    continue
                nds = parameters_to_ndarrays(fit_res.parameters)
                if nds and len(nds) > 0:
                    any_params = True
                    break

            all_done = _all_clients_done("lr_done")

            if any_params:
                try:
                    lr_result = run_server_lr(self, all_params, output_dir=self.output_dir)
                    self.logger.info("[Server] LR analysis completed successfully")
                    if lr_result is not None and len(lr_result) > 0:
                        participants = []
                        for client_proxy, fit_res in results:
                            if not fit_res.parameters:
                                continue
                            try:
                                nds = parameters_to_ndarrays(fit_res.parameters)
                                if nds and len(nds) > 0 and len(nds[0]) > 0:
                                    participants.append(client_proxy.cid)
                            except Exception:
                                continue
                        meta = {
                            "round": int(server_round),
                            "participants": participants,
                            "num_clients": len(participants),
                        }
                        meta_arr = np.frombuffer(json.dumps(meta).encode("utf-8"), dtype=np.uint8)
                        
                        # Store results for final send when all clients are done
                        self._final_lr_results = lr_result
                        self._final_lr_meta = meta_arr
                        
                        self.logger.info("[Server] Sending LR p-values back to clients for de-anonymization")
                        return ndarrays_to_parameters([lr_result, meta_arr]), {}
                    self.logger.warning("[Server] No LR results to send back to clients")
                except Exception as e:
                    self.logger.error(f"[Server] LR analysis failed: {e}")
                return None, {}

            if all_done:
                # Send final LR results if we have them stored
                if self._final_lr_results is not None and self._final_lr_meta is not None:
                    self.logger.info("[Server] All clients reported LR chunks done; sending final LR results")
                    self.logger.info("[Server] Sending final LR p-values back to clients for de-anonymization")
                    result_params = ndarrays_to_parameters([self._final_lr_results, self._final_lr_meta])
                    # Clear stored results after sending
                    self._final_lr_results = None
                    self._final_lr_meta = None
                    # Move to response stage to send results, then done in next round
                    self.current_stage = "iterative_lr_response"
                    return result_params, {}
                else:
                    self.logger.info("[Server] All clients reported LR chunks done; moving to done stage")
                    self.logger.info("[Server] Pipeline execution complete. All stages finished successfully.")
                    self.current_stage = "done"
                    self._finalize_monitoring()
                    self._apply_retention()
                    return None, {}

            self.logger.info("[Server] No LR chunks received; retrying iterative_lr stage")
            return None, {}

        if self.current_stage == "iterative_lr_response":
            # After sending final LR results, move to done stage
            self.logger.info("[Server] Final LR results sent to clients; moving to done stage")
            self.logger.info("[Server] Pipeline execution complete. All stages finished successfully.")
            self.current_stage = "done"
            self._finalize_monitoring()
            self._apply_retention()
            return None, {}

        return None, {}

    def configure_evaluate(
        self, server_round: int, parameters: Parameters, client_manager: ClientManager
    ):
        return []

    def evaluate(self, server_round: int, parameters: Parameters) -> Optional[Tuple[float, Dict[str, Scalar]]]:
        # No centralized evaluation for this pipeline
        return None

    def aggregate_evaluate(
        self,
        server_round: int,
        results,
        failures,
    ) -> Tuple[Optional[float], Dict[str, Scalar]]:
        # No evaluation aggregation
        return None, {}

    # -------------------------
    # Stage config (Scalar-only)
    # -------------------------
    def _get_stage_config(self, server_round: int) -> Dict[str, Scalar]:
        config: Dict[str, Scalar] = {"stage": self.current_stage}

        # always include round number for client-side salts / deterministic fallbacks
        config["server_round"] = int(server_round)

        if self.current_stage == "key_exchange":
            config.update(self.prg_aggregator.get_curve_params())
            return config

        all_public_keys = self.prg_aggregator.get_all_public_keys()
        if all_public_keys:
            config["all_public_keys_json"] = json.dumps(all_public_keys)
            config["num_clients"] = int(len(all_public_keys))
        else:
            config["all_public_keys_json"] = "{}"
            config["num_clients"] = int(self.num_clients or 0)

        if self.current_stage in ("sync", "global_qc"):
            config.update(self.prg_aggregator.get_curve_params())

        if self.current_stage in ("init_chunks", "init_chunks_lr"):
            config["chunk_size"] = int(self.chunk_size)

        if self.current_stage in ("local_lr", "local_lr_filter_response"):
            config["lr_tokenized"] = 1
            config["lr_token_salt"] = f"r{server_round}"
            config["lr_pad_to"] = int(self.lr_pad_to)

        if self.current_stage == "iterative_lr_response":
            # This stage sends final LR results to clients
            pass

        return config
