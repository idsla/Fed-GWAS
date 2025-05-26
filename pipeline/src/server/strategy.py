# server/strategy.py

import numpy as np
import flwr as fl

from .aggregator_qc import aggregate_global_qc
from .aggregator_king import run_server_king
from .aggregator_lr import run_server_lr, merge_insign_snp_sets
from flwr.common import parameters_to_ndarrays

# Define which stage must precede each stage for client participation
PREREQ_STAGE = {
    "global_qc": None,
    "global_qc_response": "global_qc",
    "init_chunks": "global_qc_response",
    "iterative_king": "init_chunks",
    "local_lr": "iterative_king",
    "local_lr_filter_response": "local_lr",
    "init_chunks_lr": "local_lr_filter_response",
    "iterative_lr": "init_chunks_lr",
    # sync is the first stage
    "sync": None,
}

def secure_sum(seeds):
    """
    Placeholder secure sum: simply sum the seeds.
    """
    return int(sum(seeds) % (10**9))

class FederatedGWASStrategy(fl.server.strategy.FedAvg):
    def __init__(self):
        super().__init__()
        self.global_seed = 0
        self.current_stage = "sync"
        self.chunk_size = 1000

        self.global_exclusion = []
        self.lr_data = {}
        self.participants_per_stage = {}

    def current_stage_config(self):
        return {"stage": self.current_stage}

    def on_fit_config_fn(self, rnd: int):
        if self.current_stage == "sync":
            return {"stage": "sync"}
        elif self.current_stage == "global_qc":
            return {"stage": "global_qc"}
        elif self.current_stage == "global_qc_response":
            return {"stage": "global_qc_response"}
        elif self.current_stage == "init_chunks":
            return {"stage": "init_chunks", "chunk_size": self.chunk_size}
        elif self.current_stage == "iterative_king":
            return {"stage": "iterative_king"}
        elif self.current_stage == "local_lr":
            return {"stage": "local_lr"}
        elif self.current_stage == "local_lr_filter_response":
            return {"stage": "local_lr_filter_response"}
        elif self.current_stage == "init_chunks_lr":
            return {"stage": "init_chunks_lr", "chunk_size": self.chunk_size}
        elif self.current_stage == "iterative_lr":
            return {"stage": "iterative_lr"}
        else:
            return {}

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

        if self.current_stage == "sync":
            local_seeds = []
            for _, fit_res in results:
                if fit_res.parameters:
                    ndarrays = parameters_to_ndarrays(fit_res.parameters)
                    seed_val = ndarrays[0][0]
                    #seed_val = fit_res.parameters[0].numpy[0]
                    local_seeds.append(seed_val)
            self.global_seed = secure_sum(local_seeds)
            self.current_stage = "global_qc"
            return [], {}

        elif self.current_stage == "global_qc":
            partial_data = []
            for _, fit_res in results:
                ndarrays = parameters_to_ndarrays(fit_res.parameters)

                if len(ndarrays) == 3:
                    c_arr = ndarrays[0]
                    m_arr = ndarrays[1]
                    t_arr = ndarrays[2]
                    partial_data.append([c_arr, m_arr, t_arr])
                """ if len(fit_res.parameters) == 3:
                    c_arr = fit_res.parameters[0]
                    m_arr = fit_res.parameters[1]
                    t_arr = fit_res.parameters[2] 
                    partial_data.append([c_arr, m_arr, t_arr]) """
            excl_set = aggregate_global_qc(self, partial_data, config={
                "maf_threshold": None,   # These thresholds will be unified from client data.
                "missing_threshold": None,
                "hwe_threshold": None
            })
            self.global_exclusion = sorted(list(excl_set))
            self.current_stage = "global_qc_response"
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