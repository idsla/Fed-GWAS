import flwr as fl
import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed

class MissingRateClient(fl.client.NumPyClient):
    def __init__(self, fam_path: str, bim_path: str, bed_path: str):
        """
        Initialize the MissingRateClient with file paths to .fam, .bim, and .bed files.

        Parameters:
        - fam_path (str): Path to the .fam file.
        - bim_path (str): Path to the .bim file.
        - bed_path (str): Path to the .bed file.
        """
        self.fam_path = fam_path
        self.bim_path = bim_path
        self.bed_path = bed_path
    
    def get_parameters(self, config: dict) -> list:
        """
        Get the model parameters.

        Parameters:
        - config (dict): Configuration dictionary.

        Returns:
        - list: An empty list as no parameters are used in this client.
        """
        return []

    def fit(self, parameters: list, config: dict) -> tuple:
        """
        Perform the training of the model.

        Parameters:
        - parameters (list): List of model parameters.
        - config (dict): Configuration dictionary.

        Returns:
        - tuple: Empty list, 0, and an empty dictionary as no training is performed.
        """
        return [], 0, {}

    def evaluate(self, parameters: list, config: dict) -> tuple:
        """
        Evaluate the missing rate per SNP and per individual, generating PLINK-like .imiss and .lmiss files.

        Parameters:
        - parameters (list): List of model parameters.
        - config (dict): Configuration dictionary.

        Returns:
        - tuple: Main metric (average SNP missing rate), number of examples (individuals), and a dictionary of additional metrics.
        """
        # Read genotype data using SnpReader
        snp_reader = Bed(self.bed_path)
        genotype_data = snp_reader.read().val

        n_individuals, n_snps = genotype_data.shape

        # Calculate missing rate per SNP and per individual
        missing_rate_per_snp = np.sum(np.isnan(genotype_data), axis=0) / n_individuals
        missing_rate_per_individual = np.sum(np.isnan(genotype_data), axis=1) / n_snps

        # Generate miss.imiss like file
        fam_df = pd.read_csv(self.fam_path, delim_whitespace=True, header=None)
        fam_df.columns = ['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENOTYPE']
        fam_df['F_MISS'] = missing_rate_per_individual
        imiss_file = f"client_imiss_{self.fam_path.split('/')[-1]}.csv"
        fam_df[['FID', 'IID', 'F_MISS']].to_csv(imiss_file, index=False)

        # Generate miss.lmiss like file
        bim_df = pd.read_csv(self.bim_path, delim_whitespace=True, header=None)
        bim_df.columns = ['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2']
        bim_df['F_MISS'] = missing_rate_per_snp
        lmiss_file = f"client_lmiss_{self.bim_path.split('/')[-1]}.csv"
        bim_df[['CHR', 'SNP', 'F_MISS']].to_csv(lmiss_file, index=False)
        
        # Return the main metric, the number of examples, and the dictionary of additional metrics
        main_metric = np.mean(missing_rate_per_snp)
        return main_metric, genotype_data.shape[0], {
            "snp_missing_rate": np.mean(missing_rate_per_snp),
            "individual_missing_rate": np.mean(missing_rate_per_individual),
            "imiss_file": imiss_file,
            "lmiss_file": lmiss_file,
        }

# Test or separate script for running the client
if __name__ == "__main__":
    client = MissingRateClient(
        fam_path="/Users/sonamrathod/Documents/Rutgers/Project/Git_Repo/cc.qc3.fam",
        bim_path="/Users/sonamrathod/Documents/Rutgers/Project/Git_Repo/cc.qc3.bim",
        bed_path="/Users/sonamrathod/Documents/Rutgers/Project/Git_Repo/cc.qc3.bed"
    )
    fl.client.start_numpy_client(server_address="localhost:8080", client=client)