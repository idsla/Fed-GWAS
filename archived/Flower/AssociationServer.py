import flwr as fl
import numpy as np
import statsmodels.api as sm
import json

NUM_ITERATIONS_PHASE2 = 5  # Number of iterations to share data

class GWASAggregator:
    def __init__(self):
        self.X_data = []
        self.y_data = []
        self.snp_ids = []  # Store SNP identifiers

    def add_data(self, X_part, y_part, snp_ids_part):
        # Add data part received from a client
        self.X_data.append(X_part)
        self.y_data.append(y_part)

        # Store SNP IDs only once, assuming all clients send the same IDs
        if not self.snp_ids:
            self.snp_ids = snp_ids_part

    def get_aggregated_data(self):
        # Aggregate all collected data parts
        X_combined = np.vstack(self.X_data)
        y_combined = np.concatenate(self.y_data)
        return X_combined, y_combined, self.snp_ids

    def clear(self):
        # Clear data after processing
        self.X_data = []
        self.y_data = []
        self.snp_ids = []

class GWASStrategy(fl.server.strategy.FedAvg):
    def __init__(self, aggregator: GWASAggregator):
        super().__init__()
        self.aggregator = aggregator

    def aggregate_fit(self, rnd, results, failures):
        # Collect shared genotype and phenotype data from clients
        for _, fit_res in results:
            shared_data_json = fit_res.metrics.get("shared_data")
            if shared_data_json:
                try:
                    # Parse the JSON string back into a dictionary
                    shared_data = json.loads(shared_data_json)
                    X_part = np.array(shared_data["X_sample"])
                    y_part = np.array(shared_data["y_sample"])
                    snp_ids = shared_data["snp_ids"]
                    self.aggregator.add_data(X_part, y_part, snp_ids)
                except json.JSONDecodeError as e:
                    print(f"Failed to decode shared data: {e}")
                except KeyError as e:
                    print(f"Missing key in shared data: {e}")

        # Perform regression on the final round (last iteration)
        if rnd == NUM_ITERATIONS_PHASE2:
            X_combined, y_combined, snp_ids = self.aggregator.get_aggregated_data()

            # Add intercept to the combined data
            X_combined = sm.add_constant(X_combined)

            # Perform logistic regression on aggregated data
            model = sm.Logit(y_combined, X_combined)
            try:
                result = model.fit(disp=False)
                p_values = result.pvalues.tolist()  # Convert to list for JSON transport

                # Print SNPs and their corresponding p-values
                print(f"Global SNPs and p-values for round {rnd}:")
                for snp, p_value in zip(snp_ids, p_values[1:]):  # Exclude intercept p-value
                    print(f"SNP: {snp}, P-Value: {p_value}")

                # Return p-values to clients
                return None, {"global_p_values": dict(zip(snp_ids, p_values[1:]))}
            except Exception as e:
                print(f"Logistic regression failed to converge: {e}")

            # Clear aggregator data after regression
            self.aggregator.clear()

        return None, {}

# Initialize the data aggregator
aggregator = GWASAggregator()

# Start Flower server with custom strategy
strategy = GWASStrategy(aggregator=aggregator)

fl.server.start_server(
    strategy=strategy,
    server_address="localhost:8080",
    config=fl.server.ServerConfig(num_rounds=NUM_ITERATIONS_PHASE2),
)
