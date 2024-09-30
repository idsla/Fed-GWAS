import flwr as fl
import numpy as np
from pysnptools.snpreader import Bed
from fedgwas.parameters import KINSHIP
from fedgwas.quality_control.modules.kinship import KinshipAnalyzer
from fedgwas.io.reader import IODataHandler
import csv
from typing import List, Tuple, Dict, Any

class KinshipClient(fl.client.NumPyClient):
    def __init__(self, bed_file_path: str, num_snps_to_share: int):
        """
        Initialize the KinshipClient with the path to the .bed file and the number of SNPs to share.

        Parameters:
        -----------
        bed_file_path : str
            Path to the .bed file.
        num_snps_to_share : int
            Number of SNPs to share with the server in each round.
        """
        self.bed_path = bed_file_path
        self.num_snps_to_share = num_snps_to_share
    
    def calculate_kinship_matrix(self, snp_data: np.ndarray) -> np.ndarray:
        """
        Calculate the kinship matrix based on the SNP data.

        Parameters:
        -----------
        snp_data : np.ndarray
            A 2D numpy array where rows represent individuals and columns represent SNPs.

        Returns:
        --------
        np.ndarray
            A symmetric matrix where each element [i, j] represents the kinship coefficient between individual i and j.
        """
        num_individuals = snp_data.shape[0]
        kinship_matrix = np.zeros((num_individuals, num_individuals))
        
        for i in range(num_individuals):
            for j in range(i + 1, num_individuals):
                kinship_matrix[i, j] = KinshipAnalyzer.calculate_king_coeff(self, snp_data[i], snp_data[j])
                kinship_matrix[j, i] = kinship_matrix[i, j]  # Symmetric matrix
        
        return kinship_matrix

    def select_snps_to_share(self, snp_data: np.ndarray) -> np.ndarray:
        """
        Select a subset of SNPs to share with the server.

        Parameters:
        -----------
        snp_data : np.ndarray
            A 2D numpy array where rows represent individuals and columns represent SNPs.

        Returns:
        --------
        np.ndarray
            A subset of SNPs selected randomly from the input SNP data.
        """
        num_snps = snp_data.shape[1]
        selected_snps = np.random.choice(np.arange(num_snps), self.num_snps_to_share, replace=False)
        return snp_data[:, selected_snps]

    def fit(self, parameters: Dict[str, Any], config: Dict[str, Any]) -> Tuple[List[np.ndarray], int, Dict[str, Any]]:
        """
        Compute kinship matrix and send it to the server.

        Parameters:
        -----------
        parameters : Dict[str, Any]
            Parameters passed from the server.
        config : Dict[str, Any]
            Configuration data from the server.

        Returns:
        --------
        Tuple[List[np.ndarray], int, Dict[str, Any]]:
            A tuple containing:
            - The kinship matrix as a list of NumPy arrays.
            - The number of examples (individuals).
            - A dictionary of additional metrics.
        """
        snp_data = IODataHandler._load_genotype_data(self)
        
        if snp_data.size == 0:
            print("No SNP data available for this client.")
            return [], 0, {}

        snp_subset = self.select_snps_to_share(snp_data)
        kinship_matrix = self.calculate_kinship_matrix(snp_subset)
        
        print(f"Kinship matrix shape: {kinship_matrix.shape}")
        print(f"Kinship matrix: {kinship_matrix}")
        
        num_examples = snp_data.shape[0]
        parameters_to_send = fl.common.ndarrays_to_parameters([kinship_matrix])
        return [kinship_matrix], num_examples, {} 

    def evaluate(self, parameters: List[np.ndarray], config: Dict[str, Any]) -> Tuple[float, int, Dict[str, Any]]:
        """
        Receive the global kinship matrix from the server and classify individuals.

        Parameters:
        -----------
        parameters : List[np.ndarray]
            The global kinship matrix received from the server.
        config : Dict[str, Any]
            Configuration data from the server.

        Returns:
        --------
        Tuple[float, int, Dict[str, Any]]:
            A tuple containing:
            - Evaluation loss (currently 0.0).
            - The number of examples (individuals).
            - A dictionary of additional metrics.
        """
        print(f"Received parameters type: {type(parameters)}")
        print(f"Received parameters content: {parameters}")
        global_kinship_matrix = parameters[0]

        related_pairs, unrelated_pairs = self.classify_kinship(global_kinship_matrix)
        self.save_to_csv(related_pairs, "related_pairs.csv")
        self.save_to_csv(unrelated_pairs, "unrelated_pairs.csv", is_related=False)

        return 0.0, len(global_kinship_matrix), {}

    def classify_kinship(self, global_kinship_matrix: np.ndarray) -> Tuple[List[Tuple[int, int, float]], List[Tuple[int, int, float]]]:
        """
        Classify individuals as related or unrelated based on the kinship matrix.

        Parameters:
        -----------
        global_kinship_matrix : np.ndarray
            The global kinship matrix received from the server.

        Returns:
        --------
        Tuple[List[Tuple[int, int, float]], List[Tuple[int, int, float]]]:
            A tuple containing:
            - A list of related pairs (i, j, kinship_value).
            - A list of unrelated pairs (i, j, kinship_value).
        """
        related_pairs = []
        unrelated_pairs = []

        for i in range(global_kinship_matrix.shape[0]):
            for j in range(i + 1, global_kinship_matrix.shape[1]):
                kinship_value = global_kinship_matrix[i, j]
                if kinship_value > KINSHIP['final_kinship_estimate']['firstkin']:
                    related_pairs.append((i, j, kinship_value))
                elif KINSHIP['final_kinship_estimate']['secondkin'] < kinship_value <= KINSHIP['final_kinship_estimate']['firstkin']:
                    related_pairs.append((i, j, kinship_value))
                elif KINSHIP['final_kinship_estimate']['thirdkin'] < kinship_value <= KINSHIP['final_kinship_estimate']['secondkin']:
                    related_pairs.append((i, j, kinship_value))
                else:
                    unrelated_pairs.append((i, j, kinship_value))

        return related_pairs, unrelated_pairs

    def save_to_csv(self, data: List[Tuple[int, int, float]], filename: str, is_related: bool = True):
        """
        Save related or unrelated pairs to a CSV file.

        Parameters:
        -----------
        data : List[Tuple[int, int, float]]
            The data to be saved (related or unrelated pairs).
        filename : str
            The filename for the output CSV file.
        is_related : bool, optional
            Indicates whether the data contains related pairs (default is True).
        """
        header = ["Individual 1", "Individual 2", "Kinship Coefficient"]
        relationship_type = "Related" if is_related else "Unrelated"
        
        with open(filename, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(header)
            writer.writerows(data)

        print(f"{relationship_type} pairs saved to {filename}")

def start_client(bed_file_path: str, num_snps_to_share: int):
    """
    Start the Flower client with the given .bed file and number of SNPs to share.

    Parameters:
    -----------
    bed_file_path : str
        Path to the .bed file.
    num_snps_to_share : int
        Number of SNPs to share with the server in each round.
    """
    client = KinshipClient(bed_file_path, num_snps_to_share)
    
    fl.client.start_numpy_client(
        server_address="0.0.0.0:8080",  # Replace with your actual server address
        client=client
    )

if __name__ == "__main__":
    bed_file_path = "/Users/sonamrathod/Documents/Rutgers/Project/Git_Repo/cc.qc3.bed"
    num_snps_to_share = 100  # Define how many SNPs to share in each round
    start_client(bed_file_path, num_snps_to_share)