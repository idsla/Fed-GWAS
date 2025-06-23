# server/main_server.py

import flwr as fl 
from strategy import FederatedGWASStrategy
from flwr.server import ServerConfig

def main():
    strategy = FederatedGWASStrategy(
        min_fit_clients=1,
        min_available_clients=2,
    )
    strategy.on_fit_config_fn = lambda rnd: strategy.current_stage_config()

    fl.server.start_server(
        server_address="127.0.0.1:8080",
        strategy=strategy,
        config=ServerConfig(num_rounds=50),
    )

if __name__ == "__main__":
    main()
