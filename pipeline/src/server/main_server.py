# server/main_server.py

import flwr as fl
from .strategy import FederatedGWASStrategy

def main():
    strategy = FederatedGWASStrategy()
    fl.server.start_server(
        server_address="127.0.0.1:8080",
        strategy=strategy,
        config={"num_rounds": 50},  # or however many max rounds you allow
    )

if __name__ == "__main__":
    main()