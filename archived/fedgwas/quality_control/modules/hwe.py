import logging
import numpy as np
import pandas as pd

from .base_qc_module import BaseQCModule
from fedgwas.quality_control.qc_utils import QCUtils
from fedgwas.parameters import QUALITY_CONTROL

class HWE(BaseQCModule):

    """
    perform Hardy-Weinberg Equilibrium (HWE) test using exact p-value calculation.
    """    
    
    def __init__(self, data):
        self.data = data

        super().__init__()

    def hardy_weinberg_test(self, genotype_data: pd.DataFrame, threshold: float):
        """Perform Hardy-Weinberg Equilibrium test and filter SNPs based on the threshold."""
        self.log_report(f"Perform Hardy-Weinberg Equilibrium test and filter SNPs based on the threshold {threshold}")
        try:
            # Initialize lists to store results
            hwe_results = []
            snp_indices_to_keep = []
            threshold = QUALITY_CONTROL['hardy_weinberg_test']['threshold']
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
            print(f"hardy Weinberg test {filtered_geno.head()}")
            return filtered_geno, snp_indices_to_keep
        except Exception as e:
            print(f"An error occurred while performing Hardy-Weinberg test: {e}")
            return genotype_data, []    # Initialize lists to store results