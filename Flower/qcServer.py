import flwr as fl
import numpy as np
import pandas as pd
import argparse
import logging
from io import StringIO
from  fedgwas.quality_control.qc import QCUtils 
from fedgwas.parameters import QUALITY_CONTROL

# Configure logging
logging.basicConfig(level=logging.INFO)

def start_server(selected_methods):
    # Initialize aggregated results
    aggregated_results = {method: [] for method in selected_methods}
    num_rounds = 5  # Set to 5 rounds

    # Define custom server-side strategy
    class QCServerStrategy(fl.server.strategy.FedAvg):
        def __init__(self):
            super().__init__(
                min_fit_clients=2,  # Minimum clients for training
                min_available_clients=2  # Minimum clients needed to proceed
            )

        def aggregate_fit(self, server_round, results, failures):
            # Aggregate data from clients
            nonlocal aggregated_results
            logging.info(f"Server: Received results from {len(results)} clients.")
            for _, fit_res in results:
                client_data = fit_res.metrics
                for method in selected_methods:
                    if method in client_data:
                        aggregated_results[method].append(client_data[method])

            # After final round, process the aggregated results
            if server_round == num_rounds:
                if 'hwe' in selected_methods:
                    perform_hwe_test(aggregated_results['hwe'])
                if 'maf_filter' in selected_methods:
                    aggregate_maf_filter(aggregated_results['maf_filter'])
                if 'geno' in selected_methods:
                    aggregate_geno_filter(aggregated_results['geno'])    
                if 'mind' in selected_methods:
                    aggregate_mind_filter(aggregated_results['mind'])   
                if "snp_missingness" in selected_methods:
                    aggregate_snp_missingness(aggregated_results["snp_missingness"])
                if "ind_missingness" in selected_methods:    
                    aggregate_ind_missingness(aggregated_results["ind_missingness"])
                if "calculate_maf" in selected_methods:    
                    aggregate_maf(aggregated_results["calculate_maf"])


            return None, {}
        
    def aggregate_maf(counts_json_list):
        """
        Aggregate MAF results from clients and calculate overall MAF metrics.

        Parameters:
            counts_json_list (list): List of JSON strings from clients, containing MAF data.

        Saves:
            - aggregated_maf_results.csv: CSV file with aggregated MAF data.
        """
        logging.info("Aggregating MAF results from clients.")

        maf_dfs = []
        for res in counts_json_list:
            try:
                maf_df = pd.read_json(StringIO(res), orient="split")
                maf_dfs.append(maf_df)
            except Exception as e:
                logging.error(f"Failed to parse JSON: {e}")
        
        if not maf_dfs:
            logging.error("No valid MAF dataframes received.")
            return

        try:
            # Concatenate all data
            concatenated_df = pd.concat(maf_dfs, ignore_index=True)

            # Aggregate data: Calculate weighted mean MAF
            aggregated_maf = concatenated_df.groupby("SNP").agg(
                CHR=("CHR", "first"),
                POS=("POS", "first"),
                A1=("A1", "first"),
                A2=("A2", "first"),
                TOTAL_NCHROBS=("NCHROBS", "sum"),
                WEIGHTED_MAF=("MAF", lambda x: (x * concatenated_df.loc[x.index, "NCHROBS"]).sum() / concatenated_df.loc[x.index, "NCHROBS"].sum())
            )
            logging.info("MAF filtering results")
            aggregated_maf.head()
            # Save aggregated results to CSV
            aggregated_maf.to_csv("aggregated_calculated_maf_results.csv")
            logging.info("MAF filtering results saved to 'aggregated_maf_results.csv'.")
        except Exception as e:
            logging.error(f"Error during MAF aggregation: {e}")
    
    def aggregate_snp_missingness(counts_json_list):
        """
        Aggregate SNP-level missingness data from clients.

        Parameters:
            counts_json_list (list): List of JSON strings from clients, containing SNP-level missingness data.

        Saves:
            - aggregated_snp_missingness.csv: Aggregated SNP-level missingness data.
        """
        logging.info("Aggregating SNP-level missingness results from clients.")

        snp_dfs = []
        for res in counts_json_list:
            try:
                snp_df = pd.read_json(StringIO(res), orient="split")
                snp_dfs.append(snp_df)
            except Exception as e:
                logging.error(f"Failed to parse JSON: {e}")

        if not snp_dfs:
            logging.error("No valid SNP-level missingness dataframes received.")
            return

        try:
            # Concatenate and aggregate SNP data
            concatenated_df = pd.concat(snp_dfs, ignore_index=True)
            aggregated_snp = concatenated_df.groupby("SNP").agg(
                TOTAL_N_MISS=("N_MISS", "sum"),  # Total missing genotypes across all clients
                TOTAL_N_GENO=("N_GENO", "sum")  # Total observed genotypes across all clients
            )
            aggregated_snp["F_MISS"] = aggregated_snp["TOTAL_N_MISS"] / aggregated_snp["TOTAL_N_GENO"]

            logging.info("SNP-level missingness results")
            aggregated_snp.head()
            # Save aggregated results to CSV
            aggregated_snp.to_csv("aggregated_snp_missingness.csv")
            logging.info("SNP-level missingness results saved to 'aggregated_snp_missingness.csv'.")
        except Exception as e:
            logging.error(f"Error during SNP-level missingness aggregation: {e}")
    
    def aggregate_ind_missingness(counts_json_list):
        """
        Aggregate individual-level missingness data from clients.

        Parameters:
            counts_json_list (list): List of JSON strings from clients, containing individual-level missingness data.

        Saves:
            - aggregated_ind_missingness.csv: Aggregated individual-level missingness data.
        """
        logging.info("Aggregating individual-level missingness results from clients.")

        ind_dfs = []
        for res in counts_json_list:
            try:
                ind_df = pd.read_json(StringIO(res), orient="split")
                ind_dfs.append(ind_df)
            except Exception as e:
                logging.error(f"Failed to parse JSON: {e}")

        if not ind_dfs:
            logging.error("No valid individual-level missingness dataframes received.")
            return

        try:
            # Concatenate and aggregate individual data
            concatenated_df = pd.concat(ind_dfs, ignore_index=True)
            aggregated_ind = concatenated_df.groupby(["FID", "IID"]).agg(
                TOTAL_N_MISS=("N_MISS", "sum"),  # Total missing genotypes across all clients
                TOTAL_N_GENO=("N_GENO", "sum")  # Total observed genotypes across all clients
            )
            aggregated_ind["F_MISS"] = aggregated_ind["TOTAL_N_MISS"] / aggregated_ind["TOTAL_N_GENO"]

            logging.info("Individual-level missingness")
            aggregated_ind.head()
            # Save aggregated results to CSV
            aggregated_ind.to_csv("aggregated_ind_missingness.csv")
            logging.info("Individual-level missingness results saved to 'aggregated_ind_missingness.csv'.")
        except Exception as e:
            logging.error(f"Error during individual-level missingness aggregation: {e}")

    def aggregate_mind_filter(counts_json_list):
        """
        Aggregate missing genotype rates (MIND) results from clients and calculate overall missingness metrics.

        Parameters:
            counts_json_list (list): List of JSON strings from clients, containing missingness information for individuals.

        Saves:
            - aggregated_mind_results.csv: CSV file with aggregated missingness data for individuals.
        """
        logging.info("Aggregating MIND results from clients.")

        # Load data from JSON strings
        mind_dfs = []
        for res in counts_json_list:
            try:
                mind_df = pd.read_json(StringIO(res), orient="split")
                mind_dfs.append(mind_df)
            except Exception as e:
                logging.error(f"Failed to parse JSON: {e}")
        
        if not mind_dfs:
            logging.error("No valid MIND dataframes received.")
            return

        try:
            # Concatenate all data
            concatenated_df = pd.concat(mind_dfs, ignore_index=True) 
            # Aggregate data: Calculate total N_MISS and N_GENO per individual
            aggregated_mind = concatenated_df.groupby(["FID", "IID"]).agg(
                TOTAL_N_MISS=("N_MISS", "sum"),  # Total missing genotypes across all clients
                TOTAL_N_GENO=("N_GENO", "sum")  # Total observed genotypes across all clients
            )
            aggregated_mind["F_MISS"] = aggregated_mind["TOTAL_N_MISS"] / aggregated_mind["TOTAL_N_GENO"]  # Fraction missing

            logging.info("MIND filtering results")
            aggregated_mind.head()
            # Save aggregated results to CSV
            aggregated_mind.to_csv("aggregated_mind_results.csv")
            logging.info("MIND filtering results saved to 'aggregated_mind_results.csv'.")
        except Exception as e:
            logging.error(f"Error during MIND aggregation: {e}")

    def aggregate_maf_filter(counts_json_list):
        maf_dfs = []
        for c in counts_json_list:
            try:
                # Ensure JSON strings are handled correctly
                maf_df = pd.read_json(StringIO(c), orient="split")
                maf_dfs.append(maf_df)
            except ValueError as e:
                logging.error(f"Failed to parse JSON: {e}")
                continue

        if not maf_dfs:
            logging.error("No valid MAF dataframes received.")
            return

        try:
            # Concatenate all data
            concatenated_df = pd.concat(maf_dfs, ignore_index=True)

            # Aggregate data: Calculate total NCHROBS and weighted mean MAF
            aggregated_maf = concatenated_df.groupby("SNP").agg(
                TOTAL_NCHROBS=("NCHROBS", "sum"),  # Total chromosome observations
                WEIGHTED_MAF=("MAF", lambda x: (x * concatenated_df.loc[x.index, "NCHROBS"]).sum() / concatenated_df.loc[x.index, "NCHROBS"].sum())
            )
            logging.info("MAF filtering results")
            aggregated_maf.head()
            # Save aggregated results to CSV
            aggregated_maf.to_csv("aggregated_maf_results.csv")
            logging.info("MAF filtering results saved to 'aggregated_maf_results.csv'.")
        except Exception as e:
            logging.error(f"Error during MAF aggregation: {e}") 

    def aggregate_geno_filter(counts_json_list):
        """
        Aggregate genotype missingness filtering results from clients and calculate overall missingness metrics.

        Parameters:
            counts_json_list (list): List of JSON strings from clients, containing missingness information.

        Saves:
            - aggregated_geno_results.csv: CSV file with aggregated missingness data.
        """
        logging.info("Aggregating genotype missingness results from clients.")

        geno_dfs = []
        for res in counts_json_list:
            try:
                geno_df = pd.read_json(StringIO(res), orient="split")
                geno_dfs.append(geno_df)
            except Exception as e:
                logging.error(f"Failed to parse JSON: {e}")

        if not geno_dfs:
            logging.error("No valid genotype dataframes received.")
            return

        try:
            # Concatenate all data
            concatenated_df = pd.concat(geno_dfs, ignore_index=True)

            # Aggregate data: Sum `N_MISS` and `N_GENO` for each SNP
            aggregated_geno = concatenated_df.groupby("SNP").agg(
                N_MISS=("N_MISS", "sum"),  # Total missing values across all clients
                N_GENO=("N_GENO", "sum")  # Total genotype counts across all clients
            )
            aggregated_geno["F_MISS"] = aggregated_geno["N_MISS"] / aggregated_geno["N_GENO"]  # Fraction missing

            # Save aggregated results to CSV
            aggregated_geno.to_csv("aggregated_geno_results.csv")
            logging.info("Genotype filtering results saved to 'aggregated_geno_results.csv'.")
            aggregated_geno.head()
        except Exception as e:
            logging.error(f"Error during genotype missingness aggregation: {e}")

    def perform_hwe_test(counts_json_list):
        # Aggregate counts from clients and perform HWE test
        counts_dfs = [pd.read_json(StringIO(c)) for c in counts_json_list]
        # Sum counts across clients
        total_counts = counts_dfs[0]
        for df in counts_dfs[1:]:
            total_counts = total_counts.add(df, fill_value=0)

        # Ensure total_counts is of integer type
        total_counts = total_counts.astype(int)

        # Perform HWE test on aggregated counts
        threshold = QUALITY_CONTROL['hardy_weinberg_test']['threshold']
        snps_passing_hwe = []

        hwe_output = []

        for snp in total_counts.columns:
            obs_counts = total_counts[snp]
            obs_hom1 = obs_counts.get("obs_hom1", 0)
            obs_hets = obs_counts.get("obs_hets", 0)
            obs_hom2 = obs_counts.get("obs_hom2", 0)

            # Total number of individuals
            n_individuals = obs_hom1 + obs_hets + obs_hom2

            if n_individuals == 0:
                continue

            # Calculate allele frequencies
            total_alleles = 2 * n_individuals
            p = (2 * obs_hom1 + obs_hets) / total_alleles  # Frequency of allele 1
            q = 1 - p  # Frequency of allele 2

            # Expected genotype counts under HWE
            exp_hom1 = n_individuals * p ** 2
            exp_hets = n_individuals * 2 * p * q
            exp_hom2 = n_individuals * q ** 2

            # Calculate p-value using the HWE function
            p_value = QCUtils.snphwe(obs_hets, obs_hom1, obs_hom2)

            # Determine if SNP passes the HWE threshold
            passes_threshold = p_value >= threshold
            if passes_threshold:
                snps_passing_hwe.append(snp)

            # Append detailed results for this SNP
            hwe_output.append({
                "SNP": snp,
                "OBS(HOM1)": obs_hom1,
                "OBS(HET)": obs_hets,
                "OBS(HOM2)": obs_hom2,
                "EXP(HOM1)": exp_hom1,
                "EXP(HET)": exp_hets,
                "EXP(HOM2)": exp_hom2,
                "P-VALUE": p_value,
                "PASS_HWE": passes_threshold
            })

        logging.info(f"SNPs passing HWE test on aggregated data: {len(snps_passing_hwe)}")

        # Create a DataFrame from hwe_output
        hwe_results_df = pd.DataFrame(hwe_output)

        # Save the detailed HWE results to a CSV file
        output_file = "hwe_results.csv"
        hwe_results_df.to_csv(output_file, index=False)
        logging.info(f"HWE results saved to {output_file}")
        hwe_results_df.head()
    # Start the server
    strategy = QCServerStrategy()
    fl.server.start_server(
        server_address="localhost:8080",
        strategy=strategy,
        config=fl.server.ServerConfig(num_rounds=num_rounds),
    )

if __name__ == "__main__":
    #selected_methods = ["hwe", "snp_missingness", "ind_missingness", "maf_filter", "mind", "geno", "calculate_maf"]  # Methods to execute
    parser = argparse.ArgumentParser(description="QC Pipeline Server")
    parser.add_argument('--methods', nargs='+', required=True, help='QC methods to run')
    args = parser.parse_args()
    start_server(args.methods)
