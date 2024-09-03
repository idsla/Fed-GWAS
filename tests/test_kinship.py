import pandas as pd
from fedgwas.kinship.kinship import KinshipAnalyzer, KinshipLogger
import argparse
from pyplink import PyPlink
from pandas_plink import read_plink
import logging
import datetime
from fedgwas.parameters import KINSHIP
import numpy as np

#Setup logging
logging.basicConfig(filename='king_analysis.log', level=logging.INFO, format='%(message)s')
start_time = datetime.datetime.now()
logging.info(f"Start time: {start_time.strftime('%a %b %d %H:%M:%S %Y')}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform QC on GWAS data.")
    parser.add_argument('--file', required=True, help="Path to the bed file")
    args = parser.parse_args()

    plink_file_path = args.file
    firstkin = KINSHIP['final_kinship_estimate']['firstkin']
    secondkin = KINSHIP['final_kinship_estimate']['secondkin']
    thirdkin = KINSHIP['final_kinship_estimate']['thirdkin']

    if not (firstkin > secondkin) & (firstkin > thirdkin):
        print("Wrong values enter for 1st degree Estimate. 1st degree Estimate cannot be smaller than Second degree and third degree")
    elif not (secondkin > thirdkin):
        print("Wrong values enter for 2nd degree Estimate. Second degree Estimate cannot be smaller than third degree")
    
    bed_file = PyPlink(plink_file_path)

    (bim, fam, G) = read_plink(plink_file_path)
    fraction = 0.5  # Use 50% of the data for testing

    print(G.shape)
    logging.info(f"G.shape: {G.shape}")

    # Subset the data for quicker execution
    num_individuals = int(G.shape[1] * fraction)
    num_snps = int(G.shape[0] * fraction)

    G_subset = G[:num_snps, :num_individuals]

    # Convert the first few rows and columns to a pandas DataFrame for easier viewing
    first_few_lines = pd.DataFrame(G_subset[:5, :5].compute().T)
    first_few_lines.fillna(-1, inplace=True)
    print(first_few_lines)
    logging.info(f"First few lines:\n{first_few_lines}")

    # Split data into two sets for researchers A and B
    data_A = G_subset[:, :num_individuals // 2].compute().T
    data_B = G_subset[:, num_individuals // 2:].compute().T

    snps = np.arange(num_snps)
    iterations = 5  # Adjust the number of iterations as needed

    analyzer = KinshipAnalyzer()
    combined_kinship, kinship_history = analyzer.incremental_analysis(data_A, data_B, snps, iterations)
    final_estimates = analyzer.final_kinship_estimate(combined_kinship)

    print("Final kinship estimates:")
    logging.info("Final kinship estimates:")
    for key, value in final_estimates.items():
        logging.info(f"{key}: {value}")

    # Write final kinship estimates to a CSV file
    final_estimates_df = pd.DataFrame(list(final_estimates.items()), columns=['Pair', 'Degree'])
    final_estimates_df.to_csv('final_kinship_estimates.csv', index=False)

    # Monitor the history of kinship estimates
    logger = KinshipLogger()
    logger.monitor_kinship_history(kinship_history)

    end_time = datetime.datetime.now()
    logging.info(f"End time: {end_time.strftime('%a %b %d %H:%M:%S %Y')}")