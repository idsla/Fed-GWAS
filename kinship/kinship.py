import pandas as pd
import numpy as np
import logging
import datetime
from io import reader 
# Setup logging
logging.basicConfig(filename='king_analysis.log', level=logging.INFO, format='%(message)s')

def calculate_king_coeff(genotype_i, genotype_j):
    n11= np.sum((genotype_i == 1) & (genotype_j ==1))
    n02 = np.sum((genotype_i == 2) & (genotype_j == 0))
    n20 = np.sum((genotype_i == 0) & (genotype_j == 2))

    n1_s= np.sum(genotype_i ==1)
    s_1 = np.sum(genotype_j ==1)

    if n1_s ==0:
        return 0   

    phi_ij = (2 * n11 -4 *(n02 + n20)- n1_s +s_1)/ (4 * n1_s)

    return phi_ij

def incremental_analysis(Da, Db, snps, iterations):
    combined_kinship ={}
    weights=[]
    kinship_history=[]
    for t in range(1,iterations+1):
        print(f"Iterations {t}")
        logging.info(f"Iterations {t}")
        selected_snps= snps[:t]
        
        kinship_current ={}

        for i in range(Da.shape[0]):
            for j in range(Db.shape[0]):
                genotype_i= Da[i, selected_snps]
                genotype_j= Db[j,selected_snps]
                phi_ij= calculate_king_coeff(genotype_i,genotype_j)
                kinship_current[(i,j)] = phi_ij
                #print(f"phi_ij {phi_ij}")
        #logging.info(f" kinship_current {t}: {kinship_current}")
        weights.append(1 / iterations + 1 -t)
        logging.info(f" weights {t}: {weights}")
        
        for key in kinship_current:
            if key not in combined_kinship:
                combined_kinship[key] = 0
            combined_kinship[key] += kinship_current[key] * weights[-1]

        print(f"kinship estimates after iteration {t}: {combined_kinship}")
        #logging.info(f"kinship estimates after iteration {t}: {combined_kinship}")
        #print(f"Kinship estimates after iteration {t}: {combined_kinship}")
        kinship_history.append(dict(kinship_current))

    total_weight = sum(weights)
    logging.info(f"Total weights {t}: {total_weight}")
    for key in combined_kinship:
        combined_kinship[key]= combined_kinship[key]/total_weight
    
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

def final_kinship_estimate(combined_kinship):
    """
    Categorizes kinship coefficients into degrees of relatedness.

    Parameters:
    combined_kinship (dict): A dictionary with pairs of individuals as keys and their kinship coefficients as values.

    Returns:
    dict: A dictionary with pairs of individuals as keys and their categorized degrees of relatedness as values.
    """

    final_estimates = {}
    for key, value in combined_kinship.items():
        if value > 0.35:
            final_estimates[key] = '1st degree'
            logging.info(f"1st degree Final Estimate {key} and {value}")
        elif value > 0.18:
            final_estimates[key] = '2nd degree'
            logging.info(f"1st degree Final Estimate {key} and {value}")
        elif value > 0.09:
            final_estimates[key] = '3rd degree'
        else:
            final_estimates[key] = 'unrelated'
        #print("Final Estimate ",final_estimates[key])
        logging.info(f"Final Estimate {final_estimates[key]}")
    return final_estimates

if __name__ == "__main__":
    plink_file_path = 'data/begin.cc'
    fraction = 0.2  # Use 20% of the data for testing
    bim, fam, G = reader.pandas_plink(plink_file_path)
    # Subset the data for quicker execution
    num_individuals = int(G.shape[1] * fraction)
    num_snps = int(G.shape[0] * fraction)
    G[:3]
    G_subset = G[:num_snps, :num_individuals]
  
    # Convert the first few rows and columns to a pandas DataFrame for easier viewing
    first_few_lines = pd.DataFrame(G_subset[:5, :5].compute().T)
    first_few_lines.fillna(-1, inplace=True)
    #print(first_few_lines)

    # Split data into two sets for researchers A and B
    data_A = G_subset[:, :num_individuals // 2].compute().T
    data_B = G_subset[:, num_individuals // 2:].compute().T
    G_subset
    snps = np.arange(num_snps)
    #snps
    iterations = 5  # Adjust the number of iterations as needed

    combined_kinship, kinship_history = incremental_analysis(data_A, data_B, snps, iterations)
    final_estimates = final_kinship_estimate(combined_kinship)
    monitor_kinship_history(kinship_history)
    logging.info("Final kinship estimates:")

    end_time = datetime.datetime.now()
    logging.info(f"End time: {end_time.strftime('%a %b %d %H:%M:%S %Y')}")