# client/main_client.py

import flwr as fl
import logging
import numpy as np
import sys
import os

# Add parent directory to path for server imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pipeline.src.clients.base_client import BaseGWASClient
from pipeline.src.server.prg_masking import create_client_masking_helper
from pipeline.src.clients.local_qc import (
    compute_genotype_counts,
    compute_missingness_counts,
    run_local_lr,
    parse_insignificant_snps,
    exclude_snps,
    exclude_samples_by_missing_rate,
)
from pipeline.src.clients.iterative_king import handle_iterative_king
from pipeline.src.clients.iterative_lr import handle_iterative_lr
from pipeline.src.clients.data_loder import DataLoader
import os
import uuid

class FedLRClient(BaseGWASClient):
    def __init__(self, config_file="config.yaml", partition_by="samples"):
        # Use DataLoader to load configuration and transform data if necessary
        loader = DataLoader(config_file)
        
        # Clear intermediate and log directories at the start of each run
        import shutil
        if os.path.exists(loader.intermediate_dir):
            shutil.rmtree(loader.intermediate_dir)
        os.makedirs(loader.intermediate_dir, exist_ok=True)
        if os.path.exists(loader.log_dir):
            shutil.rmtree(loader.log_dir)
        os.makedirs(loader.log_dir, exist_ok=True)
        # transform_data() returns the PLINK dataset prefix (e.g., "data/client_data")
        plink_prefix = loader.transform_data()
        client_id = f"client_{uuid.uuid4().hex[:6]}"
        super().__init__(plink_prefix, client_id=client_id, partition_by=partition_by, log_dir=loader.log_dir)
        self.intermediate_dir = loader.intermediate_dir
        self.log_dir = loader.log_dir

       # super().__init__(plink_prefix, client_id="client_1", partition_by=partition_by)
        # Overwrite thresholds from config loaded via DataLoader
        thresholds = loader.get_thresholds()
        self.maf_threshold = thresholds.get("maf_threshold", 0.01)
        self.miss_threshold = thresholds.get("missing_threshold", 0.1)
        self.hwe_threshold = thresholds.get("hwe_threshold", 1e-6)
        self.p_threshold = thresholds.get("p_threshold", 5e-3)
        # Flower configuration, e.g., server address and num_rounds, from DataLoader
        self.flower_config = loader.get_flower_config()
        # Additional parameters from config (e.g., chunk sizes)
        self.parameters = loader.get_parameters()
        # Participation flags per stage from config
        self.participation = loader.participation
        # For final LR significance, if desired
        self.lr_final = {}
        # For accumulating partial LR p-values
        self.lr_pvals = {}
        
        # PRG-MASKING helper (will be initialized when we know num_clients)
        self.masking_helper = None
        self.num_clients = 3  # Default, should be configured

    def fit(self, parameters, config):
        
        stage = config.get("stage", "key_exchange")
        
        # Exit process if this client opts out of the current stage
        if not self.participation.get(stage, True):
            self.logger.info(f"[Client {self.client_id}] Exiting: not participating in stage '{stage}'")
            sys.exit(0)

        # Initialize masking helper if needed
        if self.masking_helper is None:
            client_id = self.client_id
            self.masking_helper = create_client_masking_helper(client_id, self.num_clients)

        ################################################################################
        # Stage 1: Key Exchange - Generate and send DH public key for secure aggregation
        ################################################################################
        if stage == "key_exchange":
            self.logger.info(f"[Client {self.client_id}] Stage: key_exchange")
            curve_params = {k: v for k, v in config.items() if k == "curve"}
            public_key_pem = self.masking_helper.generate_ecc_keypair(curve_params)
            
            # Convert PEM string to bytes and then to uint8 array for transmission
            public_key_bytes = public_key_pem.encode('utf-8')
            public_key_array = np.frombuffer(public_key_bytes, dtype=np.uint8)
            
            self.logger.info(f"[Client {self.client_id}] Generated DH public key")
            return [public_key_array], 1, {"message": "public_key_sent"}
        
        ################################################################################
        # Stage 2: Sync - Send masked local seed for global random seed generation
        ################################################################################
        elif stage == "sync":
            
            self.logger.info(f"[Client {self.client_id}] Stage: sync with PRG-MASKING")
            
            # Get all public keys and compute shared secrets
            all_public_keys = config.get("all_public_keys", {})
            curve_params = {k: v for k, v in config.items() if k == "curve"}
            
            if not curve_params and "curve" not in config:
                # Extract curve params from the strategy (they should be in config)
                self.logger.warning(f"[Client {self.client_id}] Missing curve params in sync stage")
                return [np.array([self.local_seed], dtype=np.float64)], 1, {}
            
            self.masking_helper.compute_shared_secrets(all_public_keys, curve_params or config)
            
            # Apply PRG masking to local seed
            masked_seed = self.masking_helper.mask_data(self.local_seed)
            seed_array = np.array([masked_seed], dtype=np.float64)
            
            self.logger.info(f"[Client {self.client_id}] Original seed: {self.local_seed}, Masked: {masked_seed}")
            return [seed_array], 1, {}

        ################################################################################
        # Stage 3: Local QC - Filter samples by per-sample missing rate threshold
        ################################################################################
        if stage == "local_qc":
            
            local_mind_threshold = config.get("local_mind_threshold", 0.1)
            self.logger.info(f"[Client {self.client_id}] Stage: local_qc (mind={local_mind_threshold})")
            # Exclude samples with missing rate > local_mind_threshold
            new_prefix = exclude_samples_by_missing_rate(self.plink_prefix, local_mind_threshold, log_dir=self.log_dir)
            self.plink_prefix = new_prefix
            self.logger.info(f"[Client {self.client_id}] Local QC done => new prefix {self.plink_prefix}")
            return [], 1, {}

        ################################################################################
        # Stage 4: Global QC - Compute and send masked MAF, HWE, and missingness statistics
        ################################################################################
        elif stage == "global_qc":
            
            self.logger.info(f"[Client {self.client_id}] Stage: global_qc with PRG-MASKING")
            
            # Get all public keys and compute shared secrets
            all_public_keys = config.get("all_public_keys", {})
            curve_params = {k: v for k, v in config.items() if k == "curve"}
            
            if all_public_keys:
                self.masking_helper.compute_shared_secrets(all_public_keys, curve_params or config)
            
            # Compute local QC arrays
            counts_array = compute_genotype_counts(self.plink_prefix, self.client_id, log_dir=self.log_dir)
            missing_array = compute_missingness_counts(self.plink_prefix, self.client_id, log_dir=self.log_dir)
            maf_thresh = config.get("maf_threshold", 0.01)
            miss_thresh = config.get("missing_threshold", 0.1)
            hwe_thresh = config.get("hwe_threshold", 1e-6)
            threshold_array = np.array([maf_thresh, miss_thresh, hwe_thresh], dtype=np.float64)
            
            # Apply PRG masking to arrays
            if all_public_keys:
                masked_counts = self.masking_helper.mask_data(counts_array)
                masked_missing = self.masking_helper.mask_data(missing_array)
                masked_thresholds = self.masking_helper.mask_data(threshold_array)
                
                self.logger.info(f"[Client {self.client_id}] Applied PRG masking to QC arrays")
                return [masked_counts, masked_missing, masked_thresholds], 1, {}
            else:
                # Fallback to unmasked if no keys available
                self.logger.warning(f"[Client {self.client_id}] No public keys available, sending unmasked data")
                return [counts_array, missing_array, threshold_array], 1, {}

        ################################################################################
        # Stage 5: Global QC Response - Receive and apply global SNP exclusion list
        ################################################################################
        elif stage == "global_qc_response":
            
            if len(parameters) > 0:
                excluded_data = parameters[0].tobytes().decode("utf-8").split()
                
                # Suppose these are SNP IDs to drop; exclude them from local dataset
                new_prefix = exclude_snps(self.plink_prefix, excluded_data, "global_filtered", log_dir=self.log_dir)
                self.plink_prefix = new_prefix
                self.logger.info(f"[Client {self.client_id}] Global QC filter => new prefix {self.plink_prefix}")
            
            return [], 1, {}

        ################################################################################
        # Stage 6: Init Chunks (KING) - Receive global seed and partition data for KING analysis
        ################################################################################
        elif stage == "init_chunks":
            
            if len(parameters) > 0:
                self.global_seed = int(parameters[0][0])
            
            self.partition_data(config)
            self.current_chunk_idx = 0
            self.logger.info(f"[Client {self.client_id}] Created {len(self.chunk_files)} chunks for iterative KING.")
            return [], 1, {}

        ################################################################################
        # Stage 7: Iterative KING - Process KING relationship analysis chunk by chunk
        ################################################################################
        elif stage == "iterative_king":
            return handle_iterative_king(self, parameters, config)

        ################################################################################
        # Stage 8: Local LR - Run local logistic regression and identify insignificant SNPs
        ################################################################################
        elif stage == "local_lr":
            self.logger.info(f"[Client {self.client_id}] Stage: local_lr")
            p_threshold = config.get("p_threshold", 1e-3)
            assoc_file = run_local_lr(self.plink_prefix, out_prefix="local_lr_temp", log_dir=self.log_dir)
            insign_snps = parse_insignificant_snps(assoc_file, p_threshold=p_threshold)
            # Send these 'insignificant' SNPs to the server
            joined = "\n".join(insign_snps)
            data_bytes = joined.encode("utf-8")
            data_array = np.frombuffer(data_bytes, dtype=np.uint8)
            self.logger.info(f"[Client {self.client_id}] Found {len(insign_snps)} insign. SNPs locally")
            return [data_array], 1, {}

        ################################################################################
        # Stage 9: Local LR Filter Response - Receive and apply intersection of insignificant SNPs
        ################################################################################
        elif stage == "local_lr_filter_response":
            if len(parameters) > 0:
                intersection_str = parameters[0].tobytes().decode("utf-8")
                intersection_snps = intersection_str.split()
                self.logger.info(f"[Client {self.client_id}] Received intersection of {len(intersection_snps)} SNPs from server")
                if intersection_snps:
                    new_prefix = exclude_snps(self.plink_prefix, intersection_snps, "lr_filtered", log_dir=self.log_dir)
                    self.plink_prefix = new_prefix
            return [], 1, {}

        ################################################################################
        # Stage 10: Init Chunks (LR) - Partition filtered data for final LR analysis
        ################################################################################
        elif stage == "init_chunks_lr":
            self.partition_data(config)
            self.current_chunk_idx = 0
            self.logger.info(f"[Client {self.client_id}] Created {len(self.chunk_files)} chunks for iterative LR.")
            return [], 1, {}

        ################################################################################
        # Stage 11: Iterative LR - Process final logistic regression analysis chunk by chunk
        ################################################################################
        elif stage == "iterative_lr":
            return handle_iterative_lr(self, parameters, config)

        else:
            # default fallback
            return [], 1, {}

def main():
    import sys
    config_file = sys.argv[1] if len(sys.argv) > 1 else "config.yaml"
    client = FedLRClient(config_file=config_file, partition_by="samples")
    fl.client.start_numpy_client(server_address="127.0.0.1:8080", client=client)

if __name__ == "__main__":
    main()