import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed
import os
import logging as logger
import argparse
from  fedgwas.quality_control.qc import QualityControl 
from  fedgwas.quality_control.qc_utils import QCUtils 

# Configure logging
logger.basicConfig(filename='report_IO.log', filemode='w', level=logger.DEBUG,
                format='%(asctime)s - %(levelname)s - %(message)s')


def read_bed(bfile):
    try:
        bed = Bed(bfile)
        genotype_data = bed.read().val
        return genotype_data, bed
    except FileNotFoundError:
        print(f"Error: The file {bfile} does not exist.")
        return None
    except Exception as e:
        print(f"An error occurred while reading the file {bfile}: {e}")
        return None

def save_filtered_data(bim, filtered_geno, snp_indices_to_keep, output_prefix):
    try:
        
        # Load .bim and .fam data
        #bim = pd.read_csv(bfile + '.bim', delim_whitespace=True, header=None)
        #fam = pd.read_csv(bim + '.fam', delim_whitespace=True, header=None)

        # Check if SNP indices are valid
        if not all(0 <= idx < bim.shape[0] for idx in snp_indices_to_keep):
            raise IndexError("One or more SNP indices are out of range.")

        # Filter .bim file based on SNP indices to keep
        filtered_bim = bim.iloc[snp_indices_to_keep, :]

        # Write filtered data to new .bed, .bim, and .fam files
        # Write .fam file (it's the same as the input .fam file)
       # fam.to_csv(output_prefix + '.fam', sep='\t', header=False, index=False)

        # Write .bim file
        filtered_bim.to_csv(output_prefix + '.bim', sep='\t', header=False, index=False)

        # Write .bed file
        with open(output_prefix + '.bed', 'wb') as bed_file:
            # Magic number
            bed_file.write(bytearray([108, 27, 1]))

            # SNP major mode
            for i in range(filtered_geno.shape[1]):
                snp_data = filtered_geno[:, i]
                
                # Check if genotype data contains invalid values
                if not np.all(np.isin(snp_data, [-1, 0, 1, 2])):
                    raise ValueError(f"Genotype data contains invalid values at SNP index {i}.")

                # Adjust snp_data to PLINK's 2-bit format (0, 1, 2, or missing)
                snp_data[snp_data == -1] = 3  # Missing genotype coded as 3
                snp_data = snp_data.astype(np.uint8)
                
                # PLINK .bed files store genotypes in 2-bit format (four genotypes per byte)
                packed_snp_data = np.packbits(snp_data.reshape(-1, 4), bitorder='little')
                bed_file.write(packed_snp_data)
    except FileNotFoundError as fnfe:
        print(f"Error: {fnfe}")
    except IndexError as ie:
        print(f"Error: {ie}")
    except ValueError as ve:
        print(f"Error: {ve}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
# Function to read .fam file
def read_fam(file_path):
    """Read .fam file."""
    logger.info(f"Reading .fam file from {file_path}")
    return pd.read_csv(file_path, delim_whitespace=True, header=None, 
                    names=['FID', 'IID', 'PID', 'MID', 'Sex', 'Phenotype'])

# Function to read .bim file
def read_bim(file_path):
    """Read .bim file."""
    logger.info(f"Reading .bim file from {file_path}")
    return pd.read_csv(file_path, delim_whitespace=True, header=None, 
                    names=['Chromosome', 'SNP', 'Genetic_distance', 'BP', 'Allele1', 'Allele2'])

def main():
    # Setup argument parser
    parser = argparse.ArgumentParser(description="Perform QC on GWAS data.")
    parser.add_argument('--fam', required=True, help="Path to the .fam file")
    parser.add_argument('--bim', required=True, help="Path to the .bim file")
    parser.add_argument('--bed', required=True, help="Path to the .bed file")
    parser.add_argument('--threshold', type=float, default=0.05, help="Threshold for HWE p-value")

    args = parser.parse_args()

    # Read the .fam and .bim files
    fam_df = read_fam(args.fam)
    bim_df = read_bim(args.bim)
    threshold = args.threshold
    # Read the .bed file
    genotype_data, bed = read_bed(args.bed)
    geno_df = pd.DataFrame(genotype_data)
    # Print the results
    logger.info(f"Genotype data:\n{genotype_data}")
    qcInstance = QualityControl()
    qcInstance.calculate_missing_rate(geno_df, bed, 'missing_rate')
    qcInstance.filter_missingness_samples(geno_df, fam_df, threshold, 'individual_missing')
    qcInstance.geno(genotype_data, bim_df, threshold, 'snp_missing')
    filtered_geno, snp_indices_to_keep = qcInstance.hardy_weinberg_test(geno_df, threshold)
    save_filtered_data(bim_df, filtered_geno, snp_indices_to_keep, 'hwe_test')

if __name__ == "__main__":
    main()