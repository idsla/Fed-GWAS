# server/main_server.py

import flwr as fl
from .strategy import FederatedGWASStrategy
from flwr.server import ServerConfig

def main():
    num_clients = 3  # Configure expected number of clients
    strategy = FederatedGWASStrategy(num_clients=num_clients)
    strategy.on_fit_config_fn = lambda rnd: strategy.on_fit_config_fn(rnd)

    fl.server.start_server(
        server_address="127.0.0.1:8080",
        strategy=strategy,
        config=ServerConfig(
            num_rounds=50,
            min_fit_clients=1,
            min_available_clients=1,
            min_eval_clients=0,
        ),
    )

if __name__ == "__main__":
    main()
