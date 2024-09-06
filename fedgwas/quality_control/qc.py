import pandas as pd
import numpy as np
from fedgwas.quality_control.qc_utils import QCUtils
import logging
from fedgwas.parameters import QUALITY_CONTROL
# Configure logging
logging.basicConfig(filename='report_QC.log', filemode='w', level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

class QualityControl:
    def __init__(self):
        pass

    def check_sex(self, genotype_data: pd.DataFrame, bed):
        """
        Check the sex of individuals based on their homozygosity rate on the X chromosome.

        Parameters:
        - genotype_data: The DataFrame containing genotype data (unused in this function, but could be useful in other contexts).
        - bed: An instance of PyPlink or equivalent for accessing PLINK files.

        Returns:
        - A DataFrame containing the inferred sex and other information for each individual.
        """
         # Load FAM and BIM data using the bed object
        fam_df = pd.DataFrame(bed.get_fam(), columns=['fid', 'iid', 'father', 'mother', 'gender', 'status'])
        bim_df = pd.DataFrame(bed.get_bim()).reset_index()
        # Filter for X chromosome SNPs (assuming '23' is the X chromosome)
        x_chromosome_snps = bim_df[bim_df['chrom'] == 23]['snp'].dropna()
        # Process each individual
        results = []
        for index, row in fam_df.iterrows():
            observed_homoz = []
            for snp in x_chromosome_snps:
                # Get the genotype for the current SNP and individual
                genotype = bed.get_geno_marker(snp)[index]
                observed_homoz.append(genotype)
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
    
    def calculate_snpsex(genotypes: np.ndarray):
        """
        Calculate the homozygosity rate based on the provided genotypes.

        Parameters:
        - genotypes: An array of genotypes for an individual.

        Returns:
        - The homozygosity rate as a float.
        """
        # Count homozygous genotypes (0/0 and 2/2)
        homozygous = sum(genotypes == 0) + sum(genotypes == 2)
        # Calculate the total number of valid genotypes (ignoring missing data)
        total_valid = len(genotypes[genotypes >= 0])
        # Compute the homozygosity rate
        homo_rate = homozygous / total_valid if total_valid > 0 else 0
        return homo_rate
        

    def calculate_missing_rate(self, genotype_data: pd.DataFrame, bed, output_prefix: str):
        """
        Calculate missing rates per SNP and per individual.

        Parameters:
        genotype_data (pd.DataFrame): DataFrame containing genotype data where rows represent individuals and columns represent SNPs.
        bed: An object containing 'iid' (individual IDs) and 'sid' (SNP IDs).
        output_prefix (str): Prefix for the output files.

        Returns:
        Two DataFrames, one for SNP missing rates and one for individual missing rates.
        """
   
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
    
    def calculate_maf(self, bed_file: pd.DataFrame):
            '''
            This functions aims to calculate maf taking a DataFrame of bed file as input and returns a dataFrame containing maf column 
            minor allele and major allele and some other parameters
            '''
            bim_df = pd.DataFrame(bed_file.get_bim()).reset_index()
            genotype_df =self.count_genotype(bed_file)
            snps, mafs,nchrobs =[], [], []

            for snp, genotype in genotype_df.iterrows():
                no_of_chzms = 2*(genotype.get(0,0) + genotype.get(1,0) + genotype.get(2,0))
                #print(no_of_chzms,"genotype 0 ", 2*genotype.get(0,0), "genotype 1 ", genotype.get(1,0), "genotype 2 ", genotype.get(2,0))
                #break    
                #allele count
                count_a1 = 2*genotype.get(0,0) +genotype.get(1,0)
                count_a2 = 2*genotype.get(2,0) + genotype.get(1,0)
                #calculate the frequency
                freq_a1= count_a1/no_of_chzms
                freq_a2 = count_a2/no_of_chzms
                #Calculating MAF by taking minium of both frequency
                maf_coeff= min(freq_a1,freq_a2)
                mafs.append(maf_coeff)
                nchrobs.append(no_of_chzms)
                snps.append(snp)
            
            df = pd.DataFrame({
                'CHR': bim_df['chrom'],
                'SNP': snps,
                'POS': bim_df['pos'],
                'A1': bim_df['a1'],
                'A2': bim_df['a2'],
                'MAF': mafs,
                'NCHROBS': nchrobs
                })
            df.to_csv('snp_maf_data.csv', index=False)

            print(df.head())
            return df
    def count_genotype(self, bedfile):
        '''
        Generate count of genotype
        Parameters: bedfile (object)
        Return genotype_df (DataFrame) contains counts of genotype accros snp_ids
        '''
        genotype_counts={}
        for snp_id, genotypes in bedfile.iter_geno():
            unique_value, counts=np.unique(genotypes, return_counts= True)
            counts_= dict(zip(unique_value,counts))
            genotype_counts[snp_id]=counts_
        genotype_df = pd.DataFrame.from_dict(genotype_counts,orient='index')
        genotype_df = genotype_df.fillna(0)
        return genotype_df
    
    def hardy_weinberg_test(self, genotype_data: pd.DataFrame, threshold: float):
        """Perform Hardy-Weinberg Equilibrium test and filter SNPs based on the threshold."""
        try:
            # Initialize lists to store results
            hwe_results = []
            snp_indices_to_keep = []
            threshold = QUALITY_CONTROL['hwe']['threshold']
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
        n_individuals, n_snps = genotype_data.shape
        threshold = QUALITY_CONTROL['missingness']['threshold']
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
        threshold = QUALITY_CONTROL['missingness']['threshold']       

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
       '''
        Generate report based on the input parameters
        Input parameters
        output_path (string): report file name
        start_time (string): Start time of the execution
        end_time (string): Execution end time
        total_variants (int): Total genetic variants
        total_individuals (int): Total number of individuals
        filtered_snps (int): Total filtered snps
        non_missing_genotypes (int): Total genotypes which are not null
        '''
        # with open(output_path, 'w') as report_file:
        #     report_file.write(f"Genetic Analysis Report - Generated on {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        #     report_file.write(f"Analysis Duration: {end_time - start_time}\n")
        #     report_file.write(f"Total Variants Processed: {total_variants}\n")
        #     report_file.write(f"Total Individuals Processed: {total_individuals}\n")
        #     report_file.write(f"Total Non-Missing Genotypes: {total_individuals * total_variants}\n")
        #     report_file.write(f"Genotyping Rate: {non_missing_genotypes / (total_individuals * total_variants):.6f}\n")
        #     report_file.write(f"Variants Remaining after Filtering: {len(filtered_snps)}\n")
        #     report_file.write(f"Variants Removed due to Minor Allele Thresholds: {total_variants - len(filtered_snps)}\n")       
