import numpy as np
import pandas as pd
import logging as logger

from pyplink import PyPlink
from pandas_plink import read_plink
import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed

# Configure logging
logger.basicConfig(
    filename='report_IO.log',
    filemode='w',
    level=logger.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

class IODataHandler:
    """
    Class to handle loading, processing, and saving genotype data.
    """
    def __init__(self, bed_path: str, bim_path: str, fam_path: str):
        """
        Initialize the GenotypeDataHandler with the path to the .bed file.

        Parameters:
        -----------
        bed_path : str
            Path to the .bed file.
         bim_path : str
            File path to the .bim file.
        fam_path : str
            File path to the .fam file.
        """
        self.bed_path = bed_path
        self.bim_path = bim_path
        self.fam_path = fam_path

    def read_bed(self):
        try:
            bed = Bed(self.bed_path, count_A1=False)
            genotype_data = bed.read().val
            return genotype_data, bed
        except FileNotFoundError:
            print(f"Error: The file {self.bed_path} does not exist.")
            return None
        except Exception as e:
            print(f"An error occurred while reading the file {self.bed_path}: {e}")
            return None
        
    def _load_phenotype_data(self) -> np.ndarray:
        """
        Loads and preprocesses phenotype data from the .fam file.

        Returns:
        --------
        np.ndarray:
        Processed phenotype data as a NumPy array.
        """
        y = np.loadtxt(self.fam_path, usecols=[5])  # Assuming phenotype is in the 6th column
        return np.where(y == 2, 1, 0)  
    
    def _load_genotype_data(self) -> np.ndarray:
        """
        Loads genotype data using pysnptools.

        Returns:
        --------
        np.ndarray:
            Genotype data as a NumPy array.
        """
        snp_reader = Bed(self.bed_path, count_A1=False)
        return snp_reader.read().val

    def save_filtered_data(self, bim: pd.DataFrame, filtered_geno: pd.DataFrame, snp_indices_to_keep: list, output_prefix: str):
        try:
            # Check if SNP indices are valid
            if not all(0 <= idx < bim.shape[0] for idx in snp_indices_to_keep):
                raise IndexError("One or more SNP indices are out of range.")

            # Filter .bim file based on SNP indices to keep
            filtered_bim = bim.iloc[snp_indices_to_keep, :]

            # Write .bim file
            filtered_bim.to_csv(output_prefix + '.bim', sep='\t', header=False, index=False)

            # Convert DataFrame to NumPy array
            filtered_geno_array = filtered_geno.to_numpy()

            # Write .bed file
            with open(output_prefix + '.bed', 'wb') as bed_file:
                # Magic number
                bed_file.write(bytearray([108, 27, 1]))
                # SNP major mode
                for i in range(filtered_geno_array.shape[1]):
                    snp_data = filtered_geno_array[:, i]

                    # Handle NaN values by replacing them with 3 (missing genotype)
                    snp_data = np.nan_to_num(snp_data, nan=3)

                    # Check if genotype data contains invalid values
                    if not np.all(np.isin(snp_data, [0, 1, 2, 3])):
                        raise ValueError(f"Genotype data contains invalid values at SNP index {i}. Unique values: {np.unique(snp_data)}")
                    
                    # Adjust snp_data to PLINK's 2-bit format (0, 1, 2, or missing)
                    snp_data = snp_data.astype(np.uint8)
                    
                    # Ensure snp_data length is a multiple of 4
                    if len(snp_data) % 4 != 0:
                        padding_length = 4 - (len(snp_data) % 4)
                        snp_data = np.append(snp_data, [3] * padding_length)  # Padding with missing values

                    # PLINK .bed files store genotypes in 2-bit format (four genotypes per byte)
                    packed_snp_data = np.packbits(snp_data.reshape(-1, 4), bitorder='little')
                    bed_file.write(packed_snp_data)
                print(f"HWE is saved")
        except FileNotFoundError as fnfe:
            print(f"Error: {fnfe}")
        except IndexError as ie:
            print(f"Error: {ie}")
        except ValueError as ve:
            print(f"Error: {ve}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

    # Function to read .fam file
    def read_fam(self):
        """Read .fam file."""
        logger.info(f"Reading .fam file from {self.fam_path}")
        return pd.read_csv(self.fam_path, sep='\s+', header=None, 
                        names=['FID', 'IID', 'PID', 'MID', 'Sex', 'Phenotype'])

    # Function to read .bim file
    def read_bim(self):
        """Read .bim file."""
        logger.info(f"Reading .bim file from {self.bim_path}")
        return pd.read_csv(self.bim_path, sep='\s+', header=None, 
                        names=['Chromosome', 'SNP', 'Genetic_distance', 'BP', 'Allele1', 'Allele2'])                 


filename = 'data/begin.cc.bed'

def read_file(filename):
    '''
    Returns bed_file
    '''
    bed_file = PyPlink(filename)
    return bed_file

def pandas_plink(filename):
    """
    Reads PLINK genotype data from a specified file and returns the BIM, FAM, and genotype matrices.

    Parameters:
    filename (str): The path to the PLINK file prefix (without the .bed, .bim, or .fam extensions).

    Returns:
    tuple: A tuple containing three pandas DataFrames - (bim, fam, G).
        bim (DataFrame): SNP information from the BIM file.
        fam (DataFrame): Individual information from the FAM file.
        G (DataFrame): Genotype matrix.
    """
    (bim, fam, G) = read_plink(filename)
    return bim, fam, G

def bim_to_bim_df(bedfile):
    '''
    Returns bim Data Frame from .bed file
    '''
    bim_df = pd.DataFrame(bedfile.get_bim()).reset_index()
    return bim_df

def fam_to_fam_df(bedfile):
    '''
    Returns fam Data Frame from .bed file
    '''
    fam_df = pd.DataFrame(bedfile.get_fam(), columns=['fid', 'iid', 'father', 'mother', 'gender', 'status'])
    return fam_df

def snp_reader_Bed(filepath):
    snp_reader= Bed(filepath)
    return snp_reader 
