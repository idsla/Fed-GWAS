import flwr as fl
import pandas as pd
from flwr.server.strategy import FedAvg
from flwr.common import NDArrays, FitRes
from typing import List, Tuple, Dict
import os

class CustomFedAvg(FedAvg):
    def aggregate_fit(
        self,
        server_round: int,
        results: List[Tuple[NDArrays, FitRes]],
        failures: List[BaseException],
    ) -> Tuple[NDArrays, Dict[str, float]]:
        """
        Aggregate model weights and metrics from clients, and combine imiss and lmiss files.

        Parameters:
        - server_round (int): The current round of federated learning.
        - results (List[Tuple[NDArrays, FitRes]]): List of client results, each containing weights and metrics.
        - failures (List[BaseException]): List of exceptions raised by clients.

        Returns:
        - Tuple[NDArrays, Dict[str, float]]: Aggregated weights and metrics.
        """
        # Filter out clients with zero examples
        results = [res for res in results if res[1].num_examples > 0]

        if len(results) == 0:
            # If no clients have valid results, return empty
            return [], {}

        # Calculate total number of examples (should now be > 0)
        num_examples_total = sum([res[1].num_examples for res in results])

        if num_examples_total == 0:
            # If after filtering, the total number of examples is still zero, return empty
            return [], {}

        # Proceed with the standard aggregation if there are valid results
        aggregated_weights, aggregated_metrics = super().aggregate_fit(server_round, results, failures)
        
        # Aggregate the imiss and lmiss files from each client
        imiss_files = [res[1].metrics['imiss_file'] for res in results if 'imiss_file' in res[1].metrics]
        lmiss_files = [res[1].metrics['lmiss_file'] for res in results if 'lmiss_file' in res[1].metrics]

        try:
            # Aggregate imiss files
            imiss_dfs = [pd.read_csv(f) for f in imiss_files]
            combined_imiss_df = pd.concat(imiss_dfs, ignore_index=True)
            combined_imiss_df.to_csv("aggregated_imiss.csv", index=False)
            
            # Aggregate lmiss files
            lmiss_dfs = [pd.read_csv(f) for f in lmiss_files]
            combined_lmiss_df = pd.concat(lmiss_dfs, ignore_index=True)
            combined_lmiss_df.to_csv("aggregated_lmiss.csv", index=False)
        except Exception as e:
            print(f"Error during file aggregation: {e}")
            aggregated_metrics['file_aggregation_error'] = str(e)

        return aggregated_weights, aggregated_metrics

# Use this strategy in the server
strategy = CustomFedAvg()
server_config = fl.server.ServerConfig(num_rounds=10)

# Start the Flower server with the custom strategy
if __name__ == "__main__":
    fl.server.start_server(strategy=strategy, config=server_config)