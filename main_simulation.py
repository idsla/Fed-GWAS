from pipeline.server_app import app as server
from pipeline.client_app import app as client
from flwr.simulation import run_simulation
import argparse


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--num-clients", type=int, default=2)
    
    args = parser.parse_args()
    NUM_CLIENTS = args.num_clients
    backend_config = {}

    run_simulation(
        server_app=server,
        client_app=client,
        num_supernodes=NUM_CLIENTS,
        backend_config=backend_config,
        verbose_logging=True
    )