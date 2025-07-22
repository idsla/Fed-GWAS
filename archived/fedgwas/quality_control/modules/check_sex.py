import pandas as pd
import numpy as np

from .base_qc_module import BaseQCModule

class CheckSex(BaseQCModule):

    """
    Check sex
    """

    def __init__(self):

        super().__init__()

    
    def perform_diagnosis(self, data):
        """
        check sex
        """
        return data
    
    def perform_filtering(self, data):
        """
        filter
        """
        return data
    
    def check_sex(self):
        """
        Check the sex of individuals based on their homozygosity rate on the X chromosome.

        Parameters:
        - genotype_data: The DataFrame containing genotype data (unused in this function, but could be useful in other contexts).
        - bed: An instance of PyPlink or equivalent for accessing PLINK files.

        Returns:
        - A DataFrame containing the inferred sex and other information for each individual.
        """
        self.log_report("Running Check Sex")
         # Load FAM and BIM data using the bed object
        fam_df = pd.DataFrame(self.bed.get_fam(), columns=['fid', 'iid', 'father', 'mother', 'gender', 'status'])
        bim_df = pd.DataFrame(self.bed.get_bim()).reset_index()
        # Filter for X chromosome SNPs (assuming '23' is the X chromosome)
        x_chromosome_snps = bim_df[bim_df['chrom'] == 23]['snp'].dropna()
        # Process each individual
        results = []
        for index, row in fam_df.iterrows():
            observed_homoz = []
            for snp in x_chromosome_snps:
                # Get the genotype for the current SNP and individual
                genotype = self.bed.get_geno_marker(snp)[index]
                observed_homoz.append(genotype)
            #print(f"observed_homoz -->{observed_homoz}")
            observed_homozygosity = self.calculate_snpsex(np.array(observed_homoz))
            snpsexBn = 1 if observed_homozygosity > 0.8 else 2
            # Determine status (OK/PROBLEM) based on whether SNPSEX matches the reported sex
            status = 'OK' if snpsexBn == row['gender'] and snpsexBn != 0 else 'PROBLEM'
            # Calculate inbreeding coefficient (F)
            F = (observed_homozygosity - 0.5) / 0.5 if observed_homozygosity > 0.5 else 0  # Simple inbreeding coefficient calculation
            results.append({
                'FID': row['fid'],
                'IID': row['iid'],
                'PEDSEX': row['gender'],
                'SNPSEX': snpsexBn,
                'STATUS': status,
                'F': F
            })
        # Convert results to DataFrame
        results_df = pd.DataFrame(results)
        print(results_df.head())
        return results_df

    def calculate_snpsex(self,genotypes: np.ndarray):
        """
        Calculate the homozygosity rate based on the provided genotypes.

        Parameters:
        - genotypes: An array of genotypes for an individual.

        Returns:
        - The homozygosity rate as a float.
        """
        self.log_report("Running snp Sex")
        # Count homozygous genotypes (0/0 and 2/2)
        homozygous = sum(genotypes == 0) + sum(genotypes == 2)
        # Calculate the total number of valid genotypes (ignoring missing data)
        total_valid = len(genotypes[genotypes >= 0])
        # Compute the homozygosity rate
        homo_rate = homozygous / total_valid if total_valid > 0 else 0
        return homo_rate