import pandas as pd
import numpy as np
from . import qc_utils as utils
class QualityControl:
    def __init__(self):
        pass
    
    def check_sex(self, genotype_data: pd.DataFrame, phenotype_data: pd.DataFrame):
        raise NotImplementedError

    def calculate_missing_rate(self, genotype_data: pd.DataFrame):
        raise NotImplementedError


    
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
    

    def hardy_weinberg_test(self, genotype_data: pd.DataFrame):
        raise NotImplementedError

    def filter_missingness_samples(self, genotype_data: pd.DataFrame, threshold: float) -> pd.DataFrame:
        raise NotImplementedError

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

    def filter_data():
        raise NotImplementedError

    def generate_report(output_path:None, start_time: None, end_time: None, total_variants: None, total_individuals:None, filtered_snps:None, non_missing_genotypes: None ):
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
        #Generate report
        with open(output_path, 'w') as report_file:
            report_file.write(f"Genetic Analysis Report - Generated on {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            report_file.write(f"Analysis Duration: {end_time - start_time}\n")
            report_file.write(f"Total Variants Processed: {total_variants}\n")
            report_file.write(f"Total Individuals Processed: {total_individuals}\n")
            report_file.write(f"Total Non-Missing Genotypes: {total_individuals * total_variants}\n")
            report_file.write(f"Genotyping Rate: {non_missing_genotypes / (total_individuals * total_variants):.6f}\n")
            report_file.write(f"Variants Remaining after Filtering: {len(filtered_snps)}\n")
            report_file.write(f"Variants Removed due to Minor Allele Thresholds: {total_variants - len(filtered_snps)}\n")
            # report has been generated.
        
    
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
