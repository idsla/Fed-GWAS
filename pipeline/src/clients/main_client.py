# client/main_client.py

import flwr as fl
import logging
import numpy as np

from base_client import BaseGWASClient
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
from data_loader import DataLoader

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
        # For final LR significance, if desired
        self.lr_final = {}
        # For accumulating partial LR p-values
        self.lr_pvals = {}

    def fit(self, parameters, config):
        stage = config.get("stage", "sync")

        # 1) local_qc: Filter samples by per-sample missing rate
        if stage == "local_qc":
            local_mind_threshold = config.get("local_mind_threshold", 0.1)
            logging.info(f"[Client {self.client_id}] Stage: local_qc (mind={local_mind_threshold})")
            # Exclude samples with missing rate > local_mind_threshold
            new_prefix = exclude_samples_by_missing_rate(self.plink_prefix, local_mind_threshold)
            self.plink_prefix = new_prefix
            logging.info(f"[Client {self.client_id}] Local QC done => new prefix {self.plink_prefix}")
            return [], 1, {}

        # 2) global_qc: compute MAF, HWE, and per-SNP missingness
        elif stage == "global_qc":
            logging.info(f"[Client {self.client_id}] Stage: global_qc")
            # Example local computations for the server to aggregate
            counts_array = compute_genotype_counts(self.plink_prefix, self.client_id)
            missing_array = compute_missingness_counts(self.plink_prefix, self.client_id)
            # Return them to the server for global aggregation
            return [counts_array, missing_array], 1, {}

        # 3) global_qc_response: server returns SNPs or samples to exclude globally
        elif stage == "global_qc_response":
            if len(parameters) > 0:
                excluded_data = parameters[0].tobytes().decode("utf-8").split()
                # Suppose these are SNP IDs to drop; exclude them from local dataset
                new_prefix = exclude_snps(self.plink_prefix, excluded_data, "global_filtered")
                self.plink_prefix = new_prefix
                logging.info(f"[Client {self.client_id}] Global QC filter => new prefix {self.plink_prefix}")
            return [], 1, {}

        # 4) sync: combine seeds → self.global_seed
        elif stage == "sync":
            seed_np = np.array([self.local_seed], dtype=np.int64)
            return [seed_np], 1, {}

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