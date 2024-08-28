import pandas as pd
import numpy as np
from fedgwas.quality_control.qc_utils import QCUtils
import logging

# Configure logging
logging.basicConfig(filename='report_QC.log', filemode='w', level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')
class QualityControl:
    def __init__(self):
        pass

    def check_sex(self, genotype_data: pd.DataFrame, phenotype_data: pd.DataFrame):
        raise NotImplementedError

    def calculate_missing_rate(self, genotype_data: pd.DataFrame, bed, output_prefix: str):
   
        try:
            n_individuals, n_snps = genotype_data.shape

            # Calculate missing rates per SNP
            missing_per_snp = np.sum(np.isnan(genotype_data), axis=0) / n_individuals
            snp_data = pd.DataFrame({
                'SNP': bed.sid,  # Ensure bed.sid is correctly assigned
                'N_MISS': np.sum(np.isnan(genotype_data), axis=0),
                'N_GENO': n_individuals,
                'F_MISS': missing_per_snp
            })

            # Calculate missing rates per individual
            missing_per_ind = np.sum(np.isnan(genotype_data), axis=1) / n_snps
            ind_data = pd.DataFrame({
                'FID': bed.iid[:, 0],  # Ensure bed.iid is correctly assigned
                'IID': bed.iid[:, 1],
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

    def calculate_maf(self, genotype_data: pd.DataFrame):
        raise NotImplementedError

    def hardy_weinberg_test(self, genotype_data: pd.DataFrame, threshold: float):
        """Perform Hardy-Weinberg Equilibrium test and filter SNPs based on the threshold."""
        try:
            # Initialize lists to store results
            hwe_results = []
            snp_indices_to_keep = []

            for snp_idx in range(genotype_data.shape[1]):
                snp_geno = genotype_data.iloc[:, snp_idx].values  # Convert to NumPy array
                obs_hom1 = np.sum(snp_geno == 0)
                obs_hets = np.sum(snp_geno == 1)
                obs_hom2 = np.sum(snp_geno == 2)

                # Calculate p-value using the provided HWE function
                p_value = QCUtils.snphwe(obs_hets, obs_hom1, obs_hom2)
                hwe_results.append((snp_idx, p_value))

                if p_value >= threshold:
                    snp_indices_to_keep.append(snp_idx)

            # Create filtered genotype data
            filtered_geno = genotype_data.iloc[:, snp_indices_to_keep]

            return filtered_geno, snp_indices_to_keep
        except Exception as e:
            print(f"An error occurred while performing Hardy-Weinberg test: {e}")
            return genotype_data, []    # Initialize lists to store results
        
    def filter_missingness_samples(self, genotype_data: pd.DataFrame, fam, threshold: float, output_prefix) -> pd.DataFrame:
        n_individuals, n_snps = genotype_data.shape
    
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
        return filtered_fam
    
    def geno(self, genotypes, bim, threshold, output_prefix):
       
        n_individuals, n_snps = genotypes.shape
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

    def filter_maf_variants(self, genotype_data: pd.DataFrame, maf_min=None, maf_max=None, mac_min=None, mac_max=None) -> pd.DataFrame:
        """
        Filter SNPs based on minor allele frequency (MAF) and minor allele count (MAC) thresholds.

        Parameters:
        - df (DataFrame): DataFrame containing SNPs, MAF, allele counts, and total chromosome counts.
        - maf_min (float): Minimum MAF threshold.
        - maf_max (float): Maximum MAF threshold.
        - mac_min (int): Minimum MAC threshold.
        - mac_max (int): Maximum MAC threshold.

        Returns:
        - DataFrame: Filtered DataFrame based on specified criteria.
        """
         
        if maf_min is not None:
            df = df[df['MAF'] >= maf_min]
        if maf_max is not None:
            df = df[df['MAF'] <= maf_max]
        if mac_min is not None:
            # Calculating minor allele counts
            df['MAC'] = df.apply(lambda row: min(row['NCHROBS'] * row['MAF'], (1 - row['MAF']) * row['NCHROBS']), axis=1)
            df = df[df['MAC'] >= mac_min]
        if mac_max is not None:
            # Ensure MAC is already calculated
            if 'MAC' not in df.columns:
                df['MAC'] = df.apply(lambda row: min(row['NCHROBS'] * row['MAF'], (1 - row['MAF']) * row['NCHROBS']), axis=1)
            df = df[df['MAC'] <= mac_max]

        return df

    # def filter_data():
    #     raise NotImplementedError

    def generate_report(output_path):
        raise NotImplementedError
