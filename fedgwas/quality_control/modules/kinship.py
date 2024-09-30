import numpy as np
from pandas_plink import read_plink
import logging
from fedgwas.parameters import KINSHIP
from .base_qc_module import BaseQCModule

# Setup logging
logging.basicConfig(filename='king_analysis.log', level=logging.INFO, format='%(message)s')

class KinshipAnalyzer(BaseQCModule):
    """
    Performs KING kinship analysis on genotype data.
    """

    def __init__(self):
        super().__init__()
        self.firstkin = KINSHIP['final_kinship_estimate']['firstkin']
        self.secondkin = KINSHIP['final_kinship_estimate']['secondkin']
        self.thirdkin = KINSHIP['final_kinship_estimate']['thirdkin']

    def calculate_king_coeff(self, genotype_i: np.ndarray, genotype_j: np.ndarray) -> float:
        """
        Calculates the KING kinship coefficient between two individuals.

        Parameters:
        genotype_i (np.ndarray): Genotype data of the first individual.
        genotype_j (np.ndarray): Genotype data of the second individual.

        Returns:
        float: The calculated KING kinship coefficient.
        """
        n11 = np.sum((genotype_i == 1) & (genotype_j == 1))
        n02 = np.sum((genotype_i == 2) & (genotype_j == 0))
        n20 = np.sum((genotype_i == 0) & (genotype_j == 2))
        n1_s = np.sum(genotype_i == 1)
        s_1 = np.sum(genotype_j == 1)

        if n1_s == 0:
            return 0

        phi_ij = (2 * n11 - 4 * (n02 + n20) - n1_s + s_1) / (4 * n1_s)
        return phi_ij

    
    def incremental_analysis(self, Da: np.ndarray, Db: np.ndarray, snps: np.ndarray, iterations: int):
        """
        Performs incremental analysis to estimate kinship coefficients.

        Parameters:
        Da (np.ndarray): Genotype data for the first group of individuals.
        Db (np.ndarray): Genotype data for the second group of individuals.
        snps (np.ndarray): Indices of SNPs to be used in the analysis.
        iterations (int): Number of iterations for incremental analysis.

        Returns:
        tuple: (combined_kinship, kinship_history)
            combined_kinship (dict): Combined kinship coefficients.
            kinship_history (list): Kinship coefficients at each iteration.
        """
        combined_kinship = {}
        weights = []
        kinship_history = []

        for t in range(1, iterations + 1):
            logging.info(f"Iteration {t}")
            selected_snps = snps[:t]

            kinship_current = {}

            for i in range(Da.shape[0]):
                for j in range(Db.shape[0]):
                    genotype_i = Da[i, selected_snps]
                    genotype_j = Db[j, selected_snps]
                    phi_ij = self.calculate_king_coeff(genotype_i, genotype_j)
                    kinship_current[(i, j)] = phi_ij

            weight = 1 / (iterations + 1 - t)
            weights.append(weight)

            for key in kinship_current:
                if key not in combined_kinship:
                    combined_kinship[key] = 0
                combined_kinship[key] += kinship_current[key] * weights[-1]

            logging.info(f"Kinship estimates after iteration {t}: {combined_kinship}")
            kinship_history.append(kinship_current.copy())

        total_weight = sum(weights)
        for key in combined_kinship:
            combined_kinship[key] /= total_weight

        return combined_kinship, kinship_history

    def final_kinship_estimate(self, combined_kinship):
        """
        Prints and logs the kinship history for each iteration.

        Parameters:
        kinship_history (list[dict]): A list of dictionaries with kinship coefficients for each iteration.
        """
        final_estimates = {}
        for key, value in combined_kinship.items():
            if value > self.firstkin:
                final_estimates[key] = '1st degree'
                logging.info(f"1st degree Final Estimate {key} and {value}")
            elif value > self.secondkin:
                final_estimates[key] = '2nd degree'
            elif value > self.thirdkin:
                final_estimates[key] = '3rd degree'
            else:
                final_estimates[key] = 'unrelated'
        return final_estimates

class KinshipLogger:
    @staticmethod
    def monitor_kinship_history(kinship_history):
        """
        Prints and logs the kinship history for each iteration.
        """
        logging.info("Monitoring Kinship History:")
        for t, kinship_current in enumerate(kinship_history):
            logging.info(f"Iteration {t + 1}:")
            for key, value in kinship_current.items():
                logging.info(f"Pair {key}: Kinship Coefficient {value}")