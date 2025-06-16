# server/strategy.py

import numpy as np
import flwr as fl

from .aggregator_qc import aggregate_global_qc
from .aggregator_king import run_server_king
from .aggregator_lr import run_server_lr, merge_insign_snp_sets
from .prg_masking import create_prg_masking_aggregator
from flwr.common import parameters_to_ndarrays

# Define which stage must precede each stage for client participation
PREREQ_STAGE = {
    "key_exchange": None,  # New: DH key exchange stage
    "global_qc": None,
    "global_qc_response": "global_qc",
    "init_chunks": "global_qc_response",
    "iterative_king": "init_chunks",
    "local_lr": "iterative_king",
    "local_lr_filter_response": "local_lr",
    "init_chunks_lr": "local_lr_filter_response",
    "iterative_lr": "init_chunks_lr",
    # sync is the first stage after key exchange
    "sync": "key_exchange",
}

def secure_sum(seeds, num_clients=None):
    """
    Secure sum using custom secure aggregation.
    """
    if not seeds:
        return 0
    aggregator = create_secure_aggregator(num_clients or len(seeds))
    return aggregator.secure_sum_seeds(seeds)

class FederatedGWASStrategy(fl.server.strategy.FedAvg):
    def __init__(self, num_clients: int = 3):
        super().__init__()
        self.global_seed = 0
        self.current_stage = "key_exchange"  # Start with key exchange
        self.chunk_size = 1000
        self.num_clients = num_clients

        self.global_exclusion = []
        self.lr_data = {}
        self.participants_per_stage = {}
        
        # PRG-MASKING aggregator
        self.prg_aggregator = create_prg_masking_aggregator(num_clients)

    def current_stage_config(self):
        return {"stage": self.current_stage}

    def on_fit_config_fn(self, rnd: int):
        config = {"stage": self.current_stage}
        
        if self.current_stage == "key_exchange":
            # Provide DH parameters for key exchange
            config.update(self.prg_aggregator.get_dh_params())
        elif self.current_stage == "sync":
            # Provide all public keys for secure aggregation
            if self.prg_aggregator.is_key_exchange_complete():
                config["all_public_keys"] = self.prg_aggregator.get_all_public_keys()
            else:
                config = {"stage": "key_exchange"}  # Fallback to key exchange
        elif self.current_stage == "global_qc":
            # Provide public keys for QC data masking
            config["all_public_keys"] = self.prg_aggregator.get_all_public_keys()
        elif self.current_stage == "init_chunks":
            config["chunk_size"] = self.chunk_size
        elif self.current_stage == "init_chunks_lr":
            config["chunk_size"] = self.chunk_size
            
        return config

    def aggregate_fit(self, rnd: int, results, failures):
        
        # Record participants for this stage
        self.participants_per_stage[self.current_stage] = set(cid for cid, _ in results)

        # Enforce prerequisite participation: only keep clients who were in the prereq stage
        prereq = PREREQ_STAGE.get(self.current_stage)
        if prereq is not None:
            allowed = self.participants_per_stage.get(prereq, set())
            # Filter results to only those client IDs
            filtered_results = []
            for cid, fit_res in results:
                if cid in allowed:
                    filtered_results.append((cid, fit_res))
                else:
                    self.logger.info(f"Skipping client {cid} in stage '{self.current_stage}' because they did not participate in '{prereq}'")
            results = filtered_results

        if self.current_stage == "key_exchange":
            # Collect public keys from clients
            print(f"[Server] Collecting public keys from {len(results)} clients...")
            for cid, fit_res in results:
                if fit_res.parameters:
                    ndarrays = parameters_to_ndarrays(fit_res.parameters)
                    # Public key is sent as uint8 array, convert back to string
                    public_key_bytes = ndarrays[0].tobytes()
                    public_key_str = public_key_bytes.decode('utf-8')
                    self.prg_aggregator.add_client_public_key(int(cid), public_key_str)
                    print(f"  > Received public key from client {cid}")
            
            # Check if we can proceed to sync stage
            if self.prg_aggregator.is_key_exchange_complete():
                self.current_stage = "sync"
                print("[Server] Key exchange complete, proceeding to sync stage")
            
            return [], {}

        elif self.current_stage == "sync":
            # Collect masked local seeds from clients
            masked_seeds = []
            for _, fit_res in results:
                if fit_res.parameters:
                    ndarrays = parameters_to_ndarrays(fit_res.parameters)
                    masked_seed = ndarrays[0][0]  # This is already masked by client
                    masked_seeds.append(masked_seed)
            
            # Simple summation - masks cancel out automatically
            self.global_seed = int(sum(masked_seeds) % (10**9))
            print(f"[Server] PRG-MASKING seed aggregation: {len(masked_seeds)} clients -> global_seed={self.global_seed}")
            self.current_stage = "global_qc"
            return [], {}

        elif self.current_stage == "global_qc":
            # Collect masked QC data from clients
            masked_data = []
            for _, fit_res in results:
                ndarrays = parameters_to_ndarrays(fit_res.parameters)
                if len(ndarrays) == 3:
                    # These arrays are already masked by clients using PRG-MASKING
                    masked_c_arr = ndarrays[0]
                    masked_m_arr = ndarrays[1]  
                    masked_t_arr = ndarrays[2]
                    masked_data.append([masked_c_arr, masked_m_arr, masked_t_arr])
            
            # Simple summation - masks cancel out automatically in secure_sum_arrays
            excl_set = aggregate_global_qc(self, masked_data, config={
                "maf_threshold": None,
                "missing_threshold": None,
                "hwe_threshold": None
            })
            self.global_exclusion = sorted(list(excl_set))
            self.current_stage = "global_qc_response"
            print(f"[Server] PRG-MASKING QC aggregation: {len(masked_data)} clients -> {len(excl_set)} SNPs excluded")
            return [], {}

        elif self.current_stage == "global_qc_response":
            if self.global_exclusion:
                excl_str = "\n".join(str(i) for i in self.global_exclusion)
                arr = np.frombuffer(excl_str.encode("utf-8"), dtype=np.uint8)
                self.current_stage = "init_chunks"
                return [arr], {}
            else:
                self.current_stage = "init_chunks"
                return [], {}

        elif self.current_stage == "init_chunks":
            global_seed_np = np.array([self.global_seed], dtype=np.int64)
            self.current_stage = "iterative_king"
            return [global_seed_np], {}

        elif self.current_stage == "iterative_king":
            for _, fit_res in results:
                if fit_res.parameters:
                    run_server_king(self, fit_res.parameters)
            self.current_stage = "local_lr"
            return [], {}

        elif self.current_stage == "local_lr":
            all_sets = []
            for _, fit_res in results:
                if fit_res.parameters:
                    data_bytes = fit_res.parameters[0].numpy.tobytes()
                    snp_set = set(data_bytes.decode("utf-8").split())
                    all_sets.append(snp_set)
            intersection = merge_insign_snp_sets(all_sets)
            if intersection:
                arr = np.frombuffer("\n".join(sorted(intersection)).encode("utf-8"), dtype=np.uint8)
                self.current_stage = "local_lr_filter_response"
                return [arr], {}
            else:
                self.current_stage = "init_chunks_lr"
                return [], {}

        elif self.current_stage == "local_lr_filter_response":
            self.current_stage = "init_chunks_lr"
            return [], {}

        elif self.current_stage == "init_chunks_lr":
            global_seed_np = np.array([self.global_seed], dtype=np.int64)
            self.current_stage = "iterative_lr"
            return [global_seed_np], {}

        elif self.current_stage == "iterative_lr":
            for _, fit_res in results:
                if fit_res.parameters:
                    run_server_lr(self, fit_res.parameters)
            self.current_stage = "done"
            return [], {}

        return super().aggregate_fit(rnd, results, failures)