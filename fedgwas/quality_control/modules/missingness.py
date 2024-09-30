import pandas as pd
import numpy as np
import logging

from .base_qc_module import BaseQCModule
from fedgwas.parameters import QUALITY_CONTROL

class MissingDiagnosis(BaseQCModule):

    """
    Perform missing rate check and filter out SNPs with missing rate greater than a threshold.
    """

    def __init__(self, data):
        self.data = data

        super().__init__()

    def calculate_missing_rate(self, genotype_data: pd.DataFrame, output_prefix: str):
        self.log_report("Calculate missing Rate")
        try:
            n_individuals, n_snps = genotype_data.shape

            # Calculate missing rates per SNP
            missing_per_snp = np.sum(np.isnan(genotype_data), axis=0) / n_individuals
            snp_data = pd.DataFrame({
                'SNP': self.bed_snp.sid,  # Ensure bed.sid is correctly assigned
                'N_MISS': np.sum(np.isnan(genotype_data), axis=0),
                'N_GENO': n_individuals,
                'F_MISS': missing_per_snp
            })

            # Calculate missing rates per individual
            missing_per_ind = np.sum(np.isnan(genotype_data), axis=1) / n_snps
            ind_data = pd.DataFrame({
                'FID': self.bed_snp.iid[:, 0],  # Ensure bed.iid is correctly assigned
                'IID': self.bed_snp.iid[:, 1],
                'N_MISS': np.sum(np.isnan(genotype_data), axis=1),
                'N_GENO': n_snps,
                'F_MISS': missing_per_ind
            })

            # Save to files
            snp_data.to_csv(f"{output_prefix}.lmiss", sep='\t', index=False)
            ind_data.to_csv(f"{output_prefix}.imiss", sep='\t', index=False)

            return snp_data, ind_data
        except Exception as e:
            print(f"An error occurred while calculating missing rates: {e}")
            return pd.DataFrame(), pd.DataFrame()

    def filter_missingness_samples(self, genotype_data: pd.DataFrame, fam: pd.DataFrame, output_prefix) -> pd.DataFrame:
        """
        Filter individuals based on missing genotype data and save the filtered FAM file.

        Parameters:
        - genotype_data (pd.DataFrame): Genotype data matrix.
        - fam (pd.DataFrame): FAM file data.
        - output_prefix (str): Prefix for the output file name.

        Returns:
        - pd.DataFrame: Filtered FAM data.
        """
        
        self.log_report(f"Filter individuals based on missing genotype data and save the filtered FAM file")
        
        n_individuals, n_snps = genotype_data.shape
        threshold = QUALITY_CONTROL['filter_missingness_samples']['threshold']
        
        # Calculate the proportion of missing genotypes for each individual
        missing_per_ind = np.sum(np.isnan(genotype_data), axis=1) / n_snps
        
        # Identify individuals to keep
        keep_individuals = missing_per_ind <= threshold
        
        # Filter genotype matrix and individual IDs
        # filtered_genotypes = genotype_data[keep_individuals, :]
        filtered_fam = fam[keep_individuals].reset_index(drop=True)
        individuals_removed = n_individuals - filtered_fam.shape[0]
        
        log_message = (
            f"{n_snps} variants loaded from file.\n"
            f"{individuals_removed} people removed due to missing genotype data (--mind).\n"
            f"{n_snps} variants and {keep_individuals.sum()} people pass filters and QC."
        )
        print(log_message)
        logging.info(log_message)
        
        # Save the filtered data to new PLINK files
        # Write FAM file
        fam.to_csv(f"{output_prefix}.fam", sep=" ", header=False, index=False)
        print(f" filter missingness sample {fam.head()}")
        return filtered_fam
    
    def geno(self, genotypes, bim, threshold, output_prefix):
        """
        Filter SNPs based on missing genotype data and save the filtered BIM file.

        Parameters:
        - genotypes (pd.DataFrame): Genotype data matrix.
        - bim (pd.DataFrame): BIM file data.
        - output_prefix (str): Prefix for the output file name.

        Returns:
        - pd.DataFrame: Filtered BIM data.
        """
        n_individuals, n_snps = genotypes.shape
        threshold = QUALITY_CONTROL['geno']['threshold']       

        # Calculate the proportion of missing genotypes for each SNP using the filtered genotype matrix
        snp_missing_proportions = np.sum(np.isnan(genotypes), axis=0) / n_individuals
        # Identify SNPs to keep
        keep_snps = snp_missing_proportions <= threshold
        
        # Filter genotype matrix and SNPs
        filtered_bim = bim[keep_snps].reset_index(drop=True)

        snp_removed = n_snps - filtered_bim.shape[0]
        
        log_message = (
            f"{n_snps} variants loaded from file.\n"
            f"{snp_removed} SNPs removed due to missing genotype data (--geno).\n"
            f"{filtered_bim.shape[0]} variants and {n_individuals} people pass filters and QC."
        )
        print(log_message)
        logging.info(log_message)
        
        # Save the filtered data to new PLINK files
        # Write BIM file
        bim.to_csv(f"{output_prefix}.bim", sep="\t", header=False, index=False)
        return filtered_bim