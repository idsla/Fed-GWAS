import flwr as fl
import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed
import argparse
import logging
from fedgwas.io.reader import IODataHandler
from fedgwas.parameters import QUALITY_CONTROL
from  fedgwas.quality_control.qc import QualityControl 

# Configure logging
logging.basicConfig(level=logging.INFO)

class QCClient(fl.client.NumPyClient):
    def __init__(self, client_id, bed_path, bim_path, fam_path, selected_methods)    :
        self.client_id = client_id
        self.bed_path = bed_path
        self.bim_path = bim_path
        self.fam_path = fam_path
        self.selected_methods = selected_methods
        self.bed_data = IODataHandler._load_genotype_data(self) 
        self.genotype_data, self.bed = IODataHandler.read_bed(self)

    def calculate_missing_rate(self):
        n_individuals, n_snps = self.bed_data.shape

        # Calculate missing rates per SNP
        missing_per_snp = np.sum(np.isnan(self.bed_data), axis=0) / n_individuals
        snp_data = pd.DataFrame({
            'SNP': self.bed.sid,
            #'SNP': [f"rs{i}" for i in range(n_snps)],  # Replace with actual SNP IDs if available
            'N_MISS': np.sum(np.isnan(self.bed_data), axis=0),
            'N_GENO': n_individuals,
            'F_MISS': missing_per_snp
        })

        # Calculate missing rates per individual
        missing_per_ind = np.sum(np.isnan(self.bed_data), axis=1) / n_snps
        ind_data = pd.DataFrame({
            #'FID': [f"FAM{j}" for j in range(n_individuals)],  
            #'IID': [f"IND{j}" for j in range(n_individuals)],  
            'FID': self.bed.iid[:, 0],  # Ensure bed.iid is correctly assigned
            'IID': self.bed.iid[:, 1],
            'N_MISS': np.sum(np.isnan(self.bed_data), axis=1),
            'N_GENO': n_snps,
            'F_MISS': missing_per_ind
        })

        # Save reports locally
        snp_data.to_csv(f"missingness.lmiss", sep='\t', index=False)
        ind_data.to_csv(f"missingness.imiss", sep='\t', index=False)
        

        return snp_data, ind_data

    def calculate_maf(self):
        """
        Calculate Minor Allele Frequency (MAF) for each SNP in a BED file dataset.

        Returns:
            str: JSON serialized MAF data.
        """
        # Extract BIM file data (metadata) and genotype data
        bim_df = pd.DataFrame(IODataHandler.read_bim(self))
        genotype_df = pd.DataFrame(self.bed_data)

        # Initialize lists to store results
        snps, mafs, nchrobs, filtered_bim_indices = [], [], [], []

        # Iterate through genotype data to calculate MAF
        for snp_idx, genotype in genotype_df.iterrows():
            # Total chromosome count
            valid_genotypes = genotype[~genotype.isna()]  # Exclude missing genotypes
            no_of_chzms = len(valid_genotypes) * 2  # Diploid organisms

            # Skip SNPs with no observations
            if no_of_chzms == 0:
                continue

            # Allele counts
            count_a1 = 2 * (valid_genotypes == 0).sum() + (valid_genotypes == 1).sum()
            count_a2 = 2 * (valid_genotypes == 2).sum() + (valid_genotypes == 1).sum()

            # Allele frequencies
            freq_a1 = count_a1 / no_of_chzms
            freq_a2 = count_a2 / no_of_chzms

            # Minor Allele Frequency (MAF)
            maf_coeff = min(freq_a1, freq_a2)

            # Append results for this SNP
            mafs.append(maf_coeff)
            nchrobs.append(no_of_chzms)
            snps.append(bim_df.iloc[snp_idx]["SNP"])  # Use SNP ID from BIM file
            filtered_bim_indices.append(snp_idx)  # Keep track of processed BIM indices

        # Filter BIM file to include only processed SNPs
        filtered_bim_df = bim_df.iloc[filtered_bim_indices].reset_index(drop=True)

        # Create DataFrame with results
        df = pd.DataFrame({
            'CHR': filtered_bim_df['Chromosome'],
            'SNP': snps,
            'POS': filtered_bim_df['BP'],
            'A1': filtered_bim_df['Allele1'],
            'A2': filtered_bim_df['Allele2'],
            'MAF': mafs,
            'NCHROBS': nchrobs
        })
        return df


    def filter_missingness_samples(self):
        n_individuals, n_snps = self.bed_data.shape
        threshold = QUALITY_CONTROL['filter_missingness_samples']['threshold']

        # Calculate missing rates per individual
        n_miss = np.sum(np.isnan(self.bed_data), axis=1)  # Total missing genotypes per individual
        n_geno = n_snps - n_miss  # Total observed genotypes per individual
        f_miss = n_miss / n_snps  # Fraction missing per individual

        # Create DataFrame for results
        mind_data = pd.DataFrame({
           # "FID": [f"FAM{j}" for j in range(n_individuals)],  # Replace with actual family IDs if available
            #"IID": [f"IND{j}" for j in range(n_individuals)],  # Replace with actual individual IDs if available
            'FID': self.bed.iid[:, 0],  # Ensure bed.iid is correctly assigned
            'IID': self.bed.iid[:, 1],
            "N_MISS": n_miss,
            "N_GENO": n_geno,
            "F_MISS": f_miss
        })

        filtered_data = mind_data[mind_data["F_MISS"] < threshold].reset_index(drop=True)
        # Serialize results as JSON
        mind_json = filtered_data.to_json(orient="split")
        logging.info(f"Client {self.client_id}: Calculated individual-level missingness metrics.")
        return mind_json
    
    def filter_maf_variants(self, maf_df):
        """
        Filter SNPs based on MAF thresholds.
        """
        filtered_df = maf_df.copy()
        maf_min =0.01
        maf_max =0.5
        mac_min = None
        mac_max = None
        if maf_min is not None:
            filtered_df = filtered_df[filtered_df['MAF'] >= maf_min]
        if maf_max is not None:
            filtered_df = filtered_df[filtered_df['MAF'] <= maf_max]
        if mac_min is not None or mac_max is not None:
            filtered_df['MAC'] = filtered_df.apply(
                lambda row: min(row['NCHROBS'] * row['MAF'], (1 - row['MAF']) * row['NCHROBS']),
                axis=1
            )
            if mac_min is not None:
                filtered_df = filtered_df[filtered_df['MAC'] >= mac_min]
            if mac_max is not None:
                filtered_df = filtered_df[filtered_df['MAC'] <= mac_max]

        if 'MAC' in filtered_df.columns:
            filtered_df = filtered_df.drop(columns=['MAC'])
        return filtered_df
    
    def calculate_counts(self, genotype_data):
        # Compute local counts for each SNP
        counts = {}
        genotype_data_df = pd.DataFrame(genotype_data)
        for snp in range(genotype_data_df.shape[1]):
            snp_geno = genotype_data_df.iloc[:, snp].values  # NumPy array
            obs_hom1 = int(np.sum(snp_geno == 0))
            obs_hets = int(np.sum(snp_geno == 1))
            obs_hom2 = int(np.sum(snp_geno == 2))

            counts[snp] = {
                "obs_hom1": obs_hom1,
                "obs_hets": obs_hets,
                "obs_hom2": obs_hom2,
            }

        # Serialize counts to JSON
        counts_df = pd.DataFrame(counts)
        return counts_df
    
    def geno(self) -> str:
        """
        Calculate missing genotype rates for each SNP in the dataset (horizontal partitioning).

        Returns:
            str: JSON serialized DataFrame with SNP-level missingness metrics.
        """
        n_individuals, n_snps = self.bed_data.shape
        threshold = QUALITY_CONTROL['geno']['threshold']  
        # Calculate missing rates per SNP
        n_miss = np.sum(np.isnan(self.bed_data), axis=0)  # Total missing genotypes per SNP
        n_geno = n_individuals - n_miss  # Total observed genotypes per SNP
        f_miss = n_miss / n_individuals  # Fraction missing per SNP

        # Create DataFrame for results
        geno_data = pd.DataFrame({
            'SNP': self.bed.sid,
            #"SNP": [f"rs{i}" for i in range(n_snps)],  # Replace with actual SNP IDs if available
            "N_MISS": n_miss,
            "N_GENO": n_geno,
            "F_MISS": f_miss
        })

        filtered_data = geno_data[geno_data["F_MISS"] < threshold].reset_index(drop=True)
        # Serialize results as JSON
        geno_json = filtered_data.to_json(orient="split")
        logging.info(f"Client {self.client_id}: Calculated genotype missingness metrics.")
        return geno_json

    
    def fit(self, parameters, config):
        logging.info(f"Client {self.client_id} is performing fit.")
        results = {}
        num_examples = self.bed_data.shape[0]

        if 'snp_missingness' in self.selected_methods:
            snp_data, ind_data = self.calculate_missing_rate()
            results['snp_missingness'] = snp_data.to_json(orient="split")

        if 'ind_missingness' in self.selected_methods:
            snp_data, ind_data = self.calculate_missing_rate()
            results['ind_missingness'] = ind_data.to_json(orient="split")
       
        if 'mind' in self.selected_methods:
            filtered_data = self.filter_missingness_samples()
            results['mind'] = filtered_data
        
        if 'geno' in self.selected_methods:
            filtered_geno = self.geno()
            results['geno'] = filtered_geno
        
        if 'hwe' in self.selected_methods:
            # Use filtered data if available
            genotype_data = self.bed_data
            if 'filtered_data' in results:
                genotype_data = pd.read_json(results['filtered_data'], orient="split").values
            counts_df = self.calculate_counts(genotype_data)
            # Store under 'hwe' key
            results['hwe'] = counts_df.to_json()

        if 'calculate_maf' in self.selected_methods:
            df = self.calculate_maf()
            maf_json = df.to_json(orient="split")
            results['calculate_maf'] = maf_json

        if 'maf_filter' in self.selected_methods:
            maf_df = self.calculate_maf()
            maf_df = self.filter_maf_variants(maf_df)
            results['maf_filter'] = maf_df.to_json(orient="split")    

        return [], num_examples, results

    def evaluate(self, parameters, config):
        # No evaluation in this setup
        return 0.0, len(self.bed_data), {}

def start_client(client_id, bed_file_path, bim_file_path, fam_file_path, selected_methods):
    client = QCClient(client_id, bed_file_path, bim_file_path, fam_file_path, selected_methods)
    fl.client.start_numpy_client(server_address="localhost:8080", client=client)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="QC Pipeline Client")
    parser.add_argument('--client_id', type=int, required=True, help='Client ID')
    parser.add_argument('--bed_file', type=str, required=True, help='Path to BED file')
    parser.add_argument('--bim_file', type=str, required=True, help='Path to BIM file')
    parser.add_argument('--fam_file', type=str, required=True, help='Path to FAM file')
    parser.add_argument('--methods', nargs='+', required=True, help='QC methods to run')
    args = parser.parse_args()

    start_client(
        client_id = args.client_id,
        bed_file_path= args.bed_file,
        bim_file_path= args.bim_file,
        fam_file_path=args.fam_file,
        selected_methods = args.methods

        #client_id=1,
        #bed_file_path="/Users/sonamrathod/Documents/Rutgers/Project/Git_Repo/cc.qc3.bed",
        #bim_file_path="/Users/sonamrathod/Documents/Rutgers/Project/Git_Repo/cc.qc3.bim",
        #fam_file_path="/Users/sonamrathod/Documents/Rutgers/Project/Git_Repo/cc.qc3.fam",
        #selected_methods = ["hwe", "snp_missingness", "ind_missingness", "maf_filter", "mind", "geno", "calculate_maf"],  # Methods to execute
        #output_prefix="Client_Sonam"
    )
