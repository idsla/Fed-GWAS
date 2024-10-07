import flwr as fl
import numpy as np
import statsmodels.api as sm
from fedgwas.io.reader import IODataHandler as io

# Parameters
P_VALUE_THRESHOLD_PHASE1 = 0.3  # Relaxed threshold for Phase 1 filtering
P_VALUE_THRESHOLD_PHASE2 = 0.05  # Stricter threshold for Phase 2 sharing
NUM_ITERATIONS_PHASE2 = 5  # Number of iterations to share data
SAMPLE_SIZE = 50  # Number of samples to send to the server in each iteration

class GWASClient(fl.client.NumPyClient):
    def __init__(self, bed, fam):
        self.bed_path = bed
        self.fam_path = fam
        self.genotype_data = io._load_genotype_data(self)
        self.phenotype_data = io._load_phenotype_data(self)

    def _handle_missing_values(self, X):
        # Replace NaN and Inf with the column mean
        col_mean = np.nanmean(X, axis=0)
        inds = np.where(np.isnan(X))
        X[inds] = np.take(col_mean, inds[1])

        # Check for any Inf values and replace them with column mean as well
        inf_indices = np.isinf(X)
        X[inf_indices] = np.take(col_mean, np.where(inf_indices)[1])

        return X

    def _perform_local_regression(self):
        # Handle missing or inf values in genotype data
        X_data = self._handle_missing_values(self.genotype_data)
        
        # Add intercept to genotype data for logistic regression
        X_data = sm.add_constant(X_data)
        y_data = self.phenotype_data

        # Perform logistic regression
        model = sm.Logit(y_data, X_data)
        result = model.fit(disp=False)
        p_values = result.pvalues  # Extract p-values for each SNP
        
        return result, p_values

    def _filter_snps_phase1(self, p_values):
        # Identify SNPs that are above the relaxed p-value threshold for Phase 1
        significant_indices = np.where(p_values < P_VALUE_THRESHOLD_PHASE1)[0]
        return significant_indices

    def _sample_data_for_sharing(self, indices):
        # Randomly sample a subset of data for sharing in Phase 2
        num_samples = min(SAMPLE_SIZE, len(self.genotype_data))
        sample_indices = np.random.choice(len(self.genotype_data), num_samples, replace=False)

        X_sample = self.genotype_data[sample_indices, :][:, indices]
        y_sample = self.phenotype_data[sample_indices]
        return X_sample, y_sample

    def fit(self, parameters, config):
        # Phase 1: Perform local regression and filter significant SNPs
        result, p_values = self._perform_local_regression()
        significant_indices = self._filter_snps_phase1(p_values)
        
        # Phase 2: Iteratively share subsets of data with the server
        X_sample, y_sample = self._sample_data_for_sharing(significant_indices)
        
        # Return empty parameters and metrics containing filtered and sampled data
        metrics = {"shared_data": (X_sample.tolist(), y_sample.tolist())}
        return [], len(y_sample), metrics

    def evaluate(self, parameters, config):
        # Return dummy loss and metrics for compatibility
        return 0.0, len(self.phenotype_data), {}

bed_file = '/Users/sonamrathod/Documents/Rutgers/Project/Git_Repo/cc.qc3.bed'
fam_file = '/Users/sonamrathod/Documents/Rutgers/Project/Git_Repo/cc.qc3.fam'

# Create the client using the PLINK data
client = GWASClient(bed_file, fam_file)

# Start Flower client
fl.client.start_numpy_client(server_address="127.0.0.1:8080", client=client)
