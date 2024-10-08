import flwr as fl
import numpy as np
import statsmodels.api as sm

NUM_ITERATIONS_PHASE2 = 5  # Number of iterations to share data
class GWASAggregator:
    def __init__(self):
        self.X_data = []
        self.y_data = []

    def add_data(self, X_part, y_part):
        # Add data part received from a client
        self.X_data.append(X_part)
        self.y_data.append(y_part)

    def get_aggregated_data(self):
        # Aggregate all collected data parts
        X_combined = np.vstack(self.X_data)
        y_combined = np.concatenate(self.y_data)
        return X_combined, y_combined

    def clear(self):
        # Clear data after processing
        self.X_data = []
        self.y_data = []

class GWASStrategy(fl.server.strategy.FedAvg):
    def __init__(self, aggregator: GWASAggregator):
        super().__init__()
        self.aggregator = aggregator

    def aggregate_fit(self, rnd, results, failures):
        # Collect shared genotype and phenotype data from clients
        for res in results:
            shared_data = res.metrics.get("shared_data")
            if shared_data:
                X_part = np.array(shared_data[0])
                y_part = np.array(shared_data[1])
                self.aggregator.add_data(X_part, y_part)

        # Perform regression on the final round (last iteration)
        if rnd == NUM_ITERATIONS_PHASE2:
            X_combined, y_combined = self.aggregator.get_aggregated_data()
            
            # Add intercept to the combined data
            X_combined = sm.add_constant(X_combined)

            # Perform logistic regression on aggregated data
            model = sm.Logit(y_combined, X_combined)
            try:
                result = model.fit(disp=False)
                p_values = result.pvalues.tolist()  # Convert to list for JSON transport
                print(f"Global p-values for round {rnd}:\n{p_values}")

                # Return p-values to clients
                return None, {"global_p_values": p_values}
            except Exception as e:
                print(f"Logistic regression failed to converge: {e}")

            # Clear aggregator data after regression
            self.aggregator.clear()

        return None, {}

# Initialize the data aggregator
aggregator = GWASAggregator()

# Start Flower server with custom strategy
strategy = GWASStrategy(aggregator=aggregator)
fl.server.start_server(strategy=strategy, server_address="localhost:8080")
