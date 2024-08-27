from pyplink import PyPlink
import numpy as np
import pandas as pd
from pandas_plink import read_plink
import logging
import datetime
import argparse

# Setup logging
logging.basicConfig(filename='king_analysis.log', level=logging.INFO, format='%(message)s')
start_time = datetime.datetime.now()
logging.info(f"Start time: {start_time.strftime('%a %b %d %H:%M:%S %Y')}")

def calculate_king_coeff(genotype_i, genotype_j):
    n11 = np.sum((genotype_i == 1) & (genotype_j == 1))
    n02 = np.sum((genotype_i == 2) & (genotype_j == 0))
    n20 = np.sum((genotype_i == 0) & (genotype_j == 2))
    n1_s = np.sum(genotype_i == 1)
    s_1 = np.sum(genotype_j == 1)

    if n1_s == 0:
        return 0

    phi_ij = (2 * n11 - 4 * (n02 + n20) - n1_s + s_1) / (4 * n1_s)
    return phi_ij

def incremental_analysis(Da, Db, snps, iterations):
    combined_kinship = {}
    weights = []
    kinship_history = []
    for t in range(1, iterations + 1):
        print(f"Iteration {t}")
        logging.info(f"Iteration {t}")
        selected_snps = snps[:t]

        kinship_current = {}

        for i in range(Da.shape[0]):
            for j in range(Db.shape[0]):
                genotype_i = Da[i, selected_snps]
                genotype_j = Db[j, selected_snps]
                phi_ij = calculate_king_coeff(genotype_i, genotype_j)
                kinship_current[(i, j)] = phi_ij

        weight = 1 / (iterations + 1 - t)
        weights.append(weight)

        for key in kinship_current:
            if key not in combined_kinship:
                combined_kinship[key] = 0
            combined_kinship[key] += kinship_current[key] * weights[-1]

        print(f"Kinship estimates after iteration {t}: {combined_kinship}")
        logging.info(f"Kinship estimates after iteration {t}: {combined_kinship}")
        kinship_history.append(kinship_current.copy())

    total_weight = sum(weights)
    for key in combined_kinship:
        combined_kinship[key] /= total_weight

    return combined_kinship, kinship_history

def monitor_kinship_history(kinship_history):
    print("Monitoring Kinship History:")
    logging.info("Monitoring Kinship History:")
    for t, kinship_current in enumerate(kinship_history):
        print(f"Iteration {t + 1}:")
        logging.info(f"Iteration {t + 1}:")
        for key, value in kinship_current.items():
            print(f"Pair {key}: Kinship Coefficient {value}")
            logging.info(f"Pair {key}: Kinship Coefficient {value}")
        print()

def final_kinship_estimate(combined_kinship):
    final_estimates = {}
    for key, value in combined_kinship.items():
        if value > firstkin:
            final_estimates[key] = '1st degree'
            logging.info(f"1st degree Final Estimate {key} and {value}")
        elif value > secondkin:
            final_estimates[key] = '2nd degree'
        elif value > thirdkin:
            final_estimates[key] = '3rd degree'
        else:
            final_estimates[key] = 'unrelated'
    return final_estimates

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform QC on GWAS data.")
    parser.add_argument('--file', required=True, help="Path to the bed file")
    parser.add_argument('--firstkin', type=float, default=0.35, help="First Degree Kin threshold value")
    parser.add_argument('--secondkin', type=float, default=0.18, help="Second Degree Kin threshold value")
    parser.add_argument('--thirdkin', type=float, default=0.09, help="Third Degree Kin threshold value")
    args = parser.parse_args()

    plink_file_path = args.file
    firstkin =args.firstkin
    secondkin =args.secondkin
    thirdkin =args.thirdkin
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

    combined_kinship, kinship_history = incremental_analysis(data_A, data_B, snps, iterations)
    final_estimates = final_kinship_estimate(combined_kinship)

    print("Final kinship estimates:")
    logging.info("Final kinship estimates:")
    for key, value in final_estimates.items():
        logging.info(f"{key}: {value}")

    # Write final kinship estimates to a CSV file
    final_estimates_df = pd.DataFrame(list(final_estimates.items()), columns=['Pair', 'Degree'])
    final_estimates_df.to_csv('final_kinship_estimates.csv', index=False)

    # Monitor the history of kinship estimates
    monitor_kinship_history(kinship_history)

    end_time = datetime.datetime.now()
    logging.info(f"End time: {end_time.strftime('%a %b %d %H:%M:%S %Y')}")