# client/main_client.py

import flwr as fl
import logging
import numpy as np
import sys
import os

# Add parent directory to path for server imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from base_client import BaseGWASClient
from server.prg_masking import create_client_masking_helper
from local_qc import (
    compute_genotype_counts,
    compute_missingness_counts,
    run_local_lr,
    parse_insignificant_snps,
    exclude_snps,
    exclude_samples_by_missing_rate,
)
from iterative_king import handle_iterative_king
from iterative_lr import handle_iterative_lr
from data_loder import DataLoader

class FedLRClient(BaseGWASClient):
    def __init__(self, config_file="config.yaml", partition_by="samples"):
        # Use DataLoader to load configuration and transform data if necessary
        loader = DataLoader(config_file)
        # transform_data() returns the PLINK dataset prefix (e.g., "data/client_data")
        plink_prefix = loader.transform_data()
        super().__init__(plink_prefix, client_id="client_1", partition_by=partition_by)
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
        stage = config.get("stage", "sync")
        # Exit process if this client opts out of the current stage
        if not self.participation.get(stage, True):
            logging.info(f"[Client {self.client_id}] Exiting: not participating in stage '{stage}'")
            sys.exit(0)

        # Initialize masking helper if needed
        if self.masking_helper is None:
            client_id = int(self.client_id.split('_')[1]) if '_' in self.client_id else 0
            self.masking_helper = create_client_masking_helper(client_id, self.num_clients)

        # 0) key_exchange: DH key exchange for PRG-MASKING
        if stage == "key_exchange":
            logging.info(f"[Client {self.client_id}] Stage: key_exchange")
            dh_params = {k: v for k, v in config.items() if k.startswith("dh_")}
            public_key_str = self.masking_helper.generate_dh_keypair(dh_params)
            
            # Convert string to bytes and then to uint8 array for transmission
            public_key_bytes = public_key_str.encode('utf-8')
            public_key_array = np.frombuffer(public_key_bytes, dtype=np.uint8)
            
            logging.info(f"[Client {self.client_id}] Generated DH public key")
            return [public_key_array], 1, {"message": "public_key_sent"}

        # 1) local_qc: Filter samples by per-sample missing rate
        if stage == "local_qc":
            local_mind_threshold = config.get("local_mind_threshold", 0.1)
            logging.info(f"[Client {self.client_id}] Stage: local_qc (mind={local_mind_threshold})")
            # Exclude samples with missing rate > local_mind_threshold
            new_prefix = exclude_samples_by_missing_rate(self.plink_prefix, local_mind_threshold)
            self.plink_prefix = new_prefix
            logging.info(f"[Client {self.client_id}] Local QC done => new prefix {self.plink_prefix}")
            return [], 1, {}

        # 2) global_qc: compute MAF, HWE, and per-SNP missingness with PRG-MASKING
        elif stage == "global_qc":
            logging.info(f"[Client {self.client_id}] Stage: global_qc with PRG-MASKING")
            
            # Get all public keys and compute shared secrets
            all_public_keys = config.get("all_public_keys", {})
            dh_params = {k: v for k, v in config.items() if k.startswith("dh_")}
            
            if all_public_keys:
                self.masking_helper.compute_shared_secrets(all_public_keys, dh_params or config)
            
            # Compute local QC arrays
            counts_array = compute_genotype_counts(self.plink_prefix, self.client_id)
            missing_array = compute_missingness_counts(self.plink_prefix, self.client_id)
            maf_thresh = config.get("maf_threshold", 0.01)
            miss_thresh = config.get("missing_threshold", 0.1)
            hwe_thresh = config.get("hwe_threshold", 1e-6)
            threshold_array = np.array([maf_thresh, miss_thresh, hwe_thresh], dtype=np.float64)
            
            # Apply PRG masking to arrays
            if all_public_keys:
                masked_counts = self.masking_helper.mask_data(counts_array)
                masked_missing = self.masking_helper.mask_data(missing_array)
                masked_thresholds = self.masking_helper.mask_data(threshold_array)
                
                logging.info(f"[Client {self.client_id}] Applied PRG masking to QC arrays")
                return [masked_counts, masked_missing, masked_thresholds], 1, {}
            else:
                # Fallback to unmasked if no keys available
                logging.warning(f"[Client {self.client_id}] No public keys available, sending unmasked data")
                return [counts_array, missing_array, threshold_array], 1, {}

        # 3) global_qc_response: server returns SNPs or samples to exclude globally
        elif stage == "global_qc_response":
            if len(parameters) > 0:
                excluded_data = parameters[0].tobytes().decode("utf-8").split()
                # Suppose these are SNP IDs to drop; exclude them from local dataset
                new_prefix = exclude_snps(self.plink_prefix, excluded_data, "global_filtered")
                self.plink_prefix = new_prefix
                logging.info(f"[Client {self.client_id}] Global QC filter => new prefix {self.plink_prefix}")
            return [], 1, {}

        # 4) sync: combine seeds → self.global_seed using PRG-MASKING
        elif stage == "sync":
            logging.info(f"[Client {self.client_id}] Stage: sync with PRG-MASKING")
            # Get all public keys and compute shared secrets
            all_public_keys = config.get("all_public_keys", {})
            dh_params = {k: v for k, v in config.items() if k.startswith("dh_")}
            
            if not dh_params and "dh_p" not in config:
                # Extract DH params from the strategy (they should be in config)
                logging.warning(f"[Client {self.client_id}] Missing DH params in sync stage")
                return [np.array([self.local_seed], dtype=np.float64)], 1, {}
            
            self.masking_helper.compute_shared_secrets(all_public_keys, dh_params or config)
            
            # Apply PRG masking to local seed
            masked_seed = self.masking_helper.mask_data(self.local_seed)
            seed_array = np.array([masked_seed], dtype=np.float64)
            
            logging.info(f"[Client {self.client_id}] Original seed: {self.local_seed}, Masked: {masked_seed}")
            return [seed_array], 1, {}

        # 5) init_chunks: partition data for iterative KING
        elif stage == "init_chunks":
            if len(parameters) > 0:
                self.global_seed = int(parameters[0][0])
            self.partition_data(config)
            self.current_chunk_idx = 0
            logging.info(f"[Client {self.client_id}] Created {len(self.chunk_files)} chunks for iterative KING.")
            return [], 1, {}

        # 6) iterative_king: filter out highly related samples
        elif stage == "iterative_king":
            return handle_iterative_king(self, parameters, config)

        # 7) local_lr: run local LR to detect insignificant SNPs
        elif stage == "local_lr":
            logging.info(f"[Client {self.client_id}] Stage: local_lr")
            p_threshold = config.get("p_threshold", 1e-3)
            assoc_file = run_local_lr(self.plink_prefix, out_prefix="local_lr_temp")
            insign_snps = parse_insignificant_snps(assoc_file, p_threshold=p_threshold)
            # Send these ‘insignificant’ SNPs to the server
            joined = "\n".join(insign_snps)
            data_bytes = joined.encode("utf-8")
            data_array = np.frombuffer(data_bytes, dtype=np.uint8)
            logging.info(f"[Client {self.client_id}] Found {len(insign_snps)} insign. SNPs locally")
            return [data_array], 1, {}

        # 8) local_lr_filter_response: server intersection => exclude SNPs locally
        elif stage == "local_lr_filter_response":
            if len(parameters) > 0:
                intersection_str = parameters[0].tobytes().decode("utf-8")
                intersection_snps = intersection_str.split()
                logging.info(f"[Client {self.client_id}] Received intersection of {len(intersection_snps)} SNPs from server")
                if intersection_snps:
                    new_prefix = exclude_snps(self.plink_prefix, intersection_snps, "lr_filtered")
                    self.plink_prefix = new_prefix
            return [], 1, {}

        # 9) init_chunks_lr: partition data for iterative LR
        elif stage == "init_chunks_lr":
            self.partition_data(config)
            self.current_chunk_idx = 0
            logging.info(f"[Client {self.client_id}] Created {len(self.chunk_files)} chunks for iterative LR.")
            return [], 1, {}

        # 10) iterative_lr: iterative logistic regression stage
        elif stage == "iterative_lr":
            return handle_iterative_lr(self, parameters, config)

        else:
            # default fallback
            return [], 1, {}

def main():
    import flwr as fl
    client = FedLRClient(config_file="config.yaml", partition_by="samples")
    fl.client.start_numpy_client(server_address="127.0.0.1:8080", client=client)

if __name__ == "__main__":
    main()