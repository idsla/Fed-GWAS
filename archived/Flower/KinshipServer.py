import flwr as fl
import numpy as np
import sys

class CustomFedAvg(fl.server.strategy.FedAvg):
    
    def __init__(self, num_rounds: int, *args, **kwargs):
        """
        Custom strategy that extends Flower's FedAvg to aggregate kinship matrices
        and track network usage.

        Parameters
        ----------
        num_rounds : int
            Number of federated rounds to run.
        *args : tuple
            Additional arguments to be passed to the FedAvg constructor.
        **kwargs : dict
            Additional keyword arguments to be passed to the FedAvg constructor.
        """
        super().__init__(*args, **kwargs)
        self.global_kinship_matrix: np.ndarray = None
        self.num_rounds: int = num_rounds
        self.total_data_sent: int = 0
        self.total_data_received: int = 0

    def aggregate_fit(self, rnd: int, results, failures):
        """
        Aggregate the kinship matrices received from the clients and track network usage.

        Parameters
        ----------
        rnd : int
            The current round number.
        results : List[Tuple[fl.server.client_proxy.ClientProxy, fl.common.FitRes]]
            List of tuples containing client proxies and fit results.
        failures : List[BaseException]
            List of failures that occurred during the round.

        Returns
        -------
        Tuple[fl.common.Parameters, Dict[str, float]]
            The aggregated global kinship matrix and any additional metrics.
        """
        kinship_matrices = []
        for res in results:
            client_proxy, fit_res = res
            ndarrays = fl.common.parameters_to_ndarrays(fit_res.parameters)
            kinship_matrix = ndarrays[0]
            kinship_matrices.append(kinship_matrix)

            # Track data received
            self.total_data_received += sys.getsizeof(kinship_matrix)

        # Handle case where no valid kinship matrices are received
        if len(kinship_matrices) == 0:
            print("No valid kinship matrices received from clients.")
            return fl.common.ndarrays_to_parameters([self.global_kinship_matrix]), {}

        # Initialize or update the global kinship matrix
        if self.global_kinship_matrix is None:
            self.global_kinship_matrix = np.mean(kinship_matrices, axis=0)
        else:
            self.global_kinship_matrix = (self.global_kinship_matrix + np.mean(kinship_matrices, axis=0)) / 2

        # Track data sent
        global_kinship_params = fl.common.ndarrays_to_parameters([self.global_kinship_matrix])
        self.total_data_sent += sys.getsizeof(self.global_kinship_matrix)

        # Log the network usage
        print(f"Round {rnd} - Data sent: {self.total_data_sent / 1024} KB, Data received: {self.total_data_received / 1024} KB")

        # Return the global kinship matrix as a Parameters object
        return global_kinship_params, {}

def start_server(num_rounds: int) -> None:
    """
    Start the Flower server with the custom FedAvg strategy that tracks network usage.

    Parameters
    ----------
    num_rounds : int
        Number of federated rounds to run.
    """
    strategy = CustomFedAvg(num_rounds=num_rounds)

    # Start the Flower server
    fl.server.start_server(
        server_address="0.0.0.0:8080",  # Server address and port
        strategy=strategy,
        config=fl.server.ServerConfig(num_rounds=num_rounds)
    )

if __name__ == "__main__":
    num_rounds: int = 10  # Define the number of rounds
    start_server(num_rounds=num_rounds)