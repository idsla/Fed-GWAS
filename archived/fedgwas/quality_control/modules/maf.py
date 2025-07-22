import pandas as pd
import numpy as np

from .base_qc_module import BaseQCModule

class MAF(BaseQCModule):

    """
    Filter data based on Minor Allele Frequency (MAF)
    """

    def __init__(
        self,
        maf_min=None,  # TODO: Add type hint and default numerical values 
        maf_max=None,  # TODO: Add type hint and default numerical values
        mac_min=None,  # TODO: Add type hint and default numerical values
        mac_max=None   # TODO: Add type hint and default numerical values
    ):
        
        super().__init__()
        
        # statistics
        self.qc_stats = {}
        
        # params
        self.maf_min = maf_min
        self.maf_max = maf_max
        self.mac_min = mac_min
        self.mac_max = mac_max
        self.qc_params = {
            'maf_min': self.maf_min,
            'maf_max': self.maf_max,
            'mac_min': self.mac_min,
            'mac_max': self.mac_max
        }

    def perform_diagnosis(self):
        '''
        This functions aims to calculate maf taking a DataFrame of bed file as input and returns a dataFrame containing maf column 
        minor allele and major allele and some other parameters
        '''
        self.log_report("Calculate maf ")
        bim_df = pd.DataFrame(self.bed.get_bim()).reset_index()

        # perform maf calculation
        genotype_df =self.count_genotype(self.bed)
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
        
        # save caculation results to qc_stats
        self.qc_stats = {
            'CHR': bim_df['chrom'],
            'SNP': snps,
            'POS': bim_df['pos'],
            'A1': bim_df['a1'],
            'A2': bim_df['a2'],
            'MAF': mafs,
            'NCHROBS': nchrobs
        }
        
        return self.qc_stats
    
    def perform_filtering(
        self, genotype_data: pd.DataFrame, maf_min=None, maf_max=None, mac_min=None, mac_max=None
    ) -> pd.DataFrame:
        
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
        if maf_min is None:
            maf_min = self.maf_min
        if maf_max is None:
            maf_max = self.maf_max
        if mac_min is None:
            mac_min = self.mac_min
        if mac_max is None:
            mac_max = self.mac_max

        # Filter based on MAF    
        genotype_data = genotype_data[genotype_data['MAF'] >= maf_min]
        genotype_data = genotype_data[genotype_data['MAF'] <= maf_max]
        
        # Calculating minor allele counts
        genotype_data['MAC'] = genotype_data.apply(
            lambda row: min(row['NCHROBS'] * row['MAF'], (1 - row['MAF']) * row['NCHROBS']), axis=1
        )
        genotype_data = genotype_data[genotype_data['MAC'] >= mac_min]
        
        # Ensure MAC is already calculated
        if 'MAC' not in genotype_data.columns:
            genotype_data['MAC'] = genotype_data.apply(
                lambda row: min(row['NCHROBS'] * row['MAF'], (1 - row['MAF']) * row['NCHROBS']), axis=1
            )
        genotype_data = genotype_data[genotype_data['MAC'] <= mac_max]
        
        return genotype_data
    
    @staticmethod
    def count_genotype(bedfile):
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