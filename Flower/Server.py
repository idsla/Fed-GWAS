import flwr as fl
import numpy as np
import pandas as pd
from flwr.common import NDArrays, FitRes
from typing import List, Tuple, Dict
import json

class CustomFedAvg(fl.server.strategy.FedAvg):
    def __init__(self, *args, hwe_threshold: float = 0.05, **kwargs):
        super().__init__(*args, **kwargs)
        self.hwe_threshold = hwe_threshold
        self.global_geno_df = pd.DataFrame()
        self.global_sex_check_df = pd.DataFrame()
        self.global_hwe_df = pd.DataFrame()
        self.global_missingness_df = pd.DataFrame()
        self.global_missingness_snp_df = pd.DataFrame()
        self.global_missingness_ind_df = pd.DataFrame()
        self.global_maf_df = pd.DataFrame()
        self.global_maf_filter_df = pd.DataFrame()
        self.global_kinship_df = pd.DataFrame()

    def aggregate_fit(
        self,
        server_round: int,
        results: List[Tuple[NDArrays, FitRes]],
        failures: List[BaseException],
    ) -> Tuple[NDArrays, Dict[str, float]]:
        """Aggregate results from clients based on function type (geno, check sex, hwd, or missingness)."""
        
        geno_results = []
        sex_check_results = []
        hwd_results = []
        missingness_results_snp = []
        missingness_results_ind = []
        maf_results = []
        maf_filter_results = []
        kinship_results = []

        for _, fit_res in results:
            metrics = fit_res.metrics
            
            if 'geno_results' in metrics:
                print("Aggregating geno results")
                geno_results.append(metrics)
            elif 'sex_check_results' in metrics:
                print("Aggregating check sex results")
                sex_check_results.extend(json.loads(metrics['sex_check_results']))
            elif 'hwd_results' in metrics:
                print("Aggregating hwd results")
                hwd_results.append(metrics)
            elif 'missing_rate_snp' in metrics:
                print("Aggregating SNP missingness results")
                missingness_results_snp.extend(metrics['missing_rate_snp'])
            elif 'missing_rate_ind' in metrics:
                print("Aggregating individual missingness results")
                missingness_results_ind.extend(metrics['missing_rate_ind'])
            elif 'maf_results' in metrics:
                print("Aggregating MAF results")
                maf_results.extend(json.loads(metrics['maf_results']))
            elif 'filtered_maf_results' in metrics:
                print("Aggregating filtered MAF results")
                maf_filter_results.extend(json.loads(metrics['filtered_maf_results']))
            elif 'kinship_results' in metrics:
                print("Aggregating kinship results")
                kinship_results.append(json.loads(metrics['kinship_results']))

        # Aggregate geno results if present
        if geno_results:
            combined_chr = []
            combined_snp = []
            combined_pos = []
            combined_nchrobs = []
            combined_missingness = []
            for res in geno_results:
                combined_chr.extend(json.loads(res['CHR']))
                combined_snp.extend(json.loads(res['SNP']))
                combined_pos.extend(json.loads(res['POS']))
                combined_nchrobs.extend(json.loads(res['NCHROBS']))
                combined_missingness.extend(json.loads(res['MISSINGNESS']))
            
            combined_geno_df = pd.DataFrame({
                'CHR': combined_chr,
                'SNP': combined_snp,
                'POS': combined_pos,
                'MISSINGNESS': combined_missingness,
                'NCHROBS': combined_nchrobs
            })
            self.global_geno_df = pd.concat([self.global_geno_df, combined_geno_df], ignore_index=True)
            print(f"Aggregated geno results at round {server_round}:")
            print(self.global_geno_df.tail())

        # Aggregate check sex results if present
        if sex_check_results:
            combined_sex_check_df = pd.DataFrame(sex_check_results)
            self.global_sex_check_df = pd.concat([self.global_sex_check_df, combined_sex_check_df], ignore_index=True)
            print(f"Aggregated sex check results at round {server_round}:")
            print(self.global_sex_check_df.tail())

        # Aggregate hwd results if present
        if hwd_results:
            combined_snp = []
            combined_p_value = []
            for res in hwd_results:
                combined_snp.extend(json.loads(res['SNP']))
                combined_p_value.extend(json.loads(res['p_value']))
            
            combined_hwe_df = pd.DataFrame({
                'SNP': combined_snp,
                'p_value': combined_p_value
            })
            
            # Filter based on the HWD threshold
            filtered_hwe_df = combined_hwe_df[combined_hwe_df['p_value'] >= self.hwe_threshold]
            self.global_hwe_df = pd.concat([self.global_hwe_df, filtered_hwe_df], ignore_index=True)
            print(f"Aggregated hwd results at round {server_round}:")
            print(filtered_hwe_df.head())

        # Aggregate missingness results for SNPs if present
        if missingness_results_snp:
            combined_missingness_snp_df = pd.DataFrame(missingness_results_snp)
            self.global_missingness_snp_df = pd.concat([self.global_missingness_snp_df, combined_missingness_snp_df], ignore_index=True)
            print(f"Aggregated SNP missingness results at round {server_round}:")
            print(self.global_missingness_snp_df.tail())

        # Aggregate missingness results for individuals if present
        if missingness_results_ind:
            combined_missingness_ind_df = pd.DataFrame(missingness_results_ind)
            self.global_missingness_ind_df = pd.concat([self.global_missingness_ind_df, combined_missingness_ind_df], ignore_index=True)
            print(f"Aggregated individual missingness results at round {server_round}:")
            print(self.global_missingness_ind_df.tail())

        if maf_results:
            combined_maf_df = pd.DataFrame(maf_results)
            self.global_maf_df = pd.concat([self.global_maf_df, combined_maf_df], ignore_index=True)
            print(f"Aggregated MAF results:")
            print(self.global_maf_df.head())
        
        if maf_filter_results:
            combined_maf_filter_df = pd.DataFrame(maf_filter_results)
            self.global_maf_filter_df = pd.concat([self.global_maf_filter_df, combined_maf_filter_df], ignore_index=True)
            print(f"Aggregated filtered MAF results:")
            print(self.global_maf_filter_df.head())
        if kinship_results:
            aggregated_kinship = {}  
            for kinship in kinship_results:
                for pair, details in kinship.items():
                    if pair not in aggregated_kinship:
                        aggregated_kinship[pair] = details
                    else:
                        # Merge or reconcile results if needed
                        aggregated_kinship[pair]['value'] += details['value']
                        # You can adjust logic based on specific reconciliation needs

            # Convert aggregated kinship results to DataFrame for further processing
            kinship_df = pd.DataFrame.from_dict(aggregated_kinship, orient='index')
            self.global_kinship_df = kinship_df
            print(f"Aggregated kinship results at round {server_round}:")
            print(kinship_df.head())

        # if kinship_results:
        #     combined_kinship_pairs = []
        #     combined_kinship_values = []
        #     combined_kinship_degrees = []

        #     for kinship in kinship_results:
        #         for pair, degree in kinship.items():
        #             combined_kinship_pairs.append(pair)
        #             combined_kinship_values.append(kinship[pair]['value'])
        #             combined_kinship_degrees.append(degree)

        #     combined_kinship_df = pd.DataFrame({
        #         'Pair': combined_kinship_pairs,
        #         'Kinship_Value': combined_kinship_values,
        #         'Degree': combined_kinship_degrees
        #     })
        #     self.global_kinship_df = pd.concat([self.global_kinship_df, combined_kinship_df], ignore_index=True)
        #     print(f"Aggregated kinship results at round {server_round}:")
        #     print(self.global_kinship_df.tail())

        return [], {}
    
    

    def on_conclude(self):
        # Save aggregated geno results if they exist
        if not self.global_geno_df.empty:
            print("Final Aggregated Geno DataFrame:")
            print(self.global_geno_df.head())
            self.global_geno_df.to_csv('pipeline_final_geno_aggregated.csv', index=False)

        # Save aggregated sex check results if they exist
        if not self.global_sex_check_df.empty:
            print("Final Aggregated Sex Check DataFrame:")
            print(self.global_sex_check_df.tail())
            self.global_sex_check_df.to_csv('pipeline_final_sex_check_results.csv', index=False)

        # Save aggregated hwd results if they exist
        if not self.global_hwe_df.empty:
            print("Final Aggregated HWD DataFrame:")
            print(self.global_hwe_df.head())
            self.global_hwe_df.to_csv('pipeline_final_hwd_results.csv', index=False)
        
        if not self.global_missingness_df.empty:
            print("Final Aggregated Missingness DataFrame:")
            print(self.global_missingness_df.head())
            self.global_missingness_df.to_csv('pipeline_final_missingness_results.csv', index=False)

        # Save aggregated SNP missingness results if they exist
        if not self.global_missingness_snp_df.empty:
            print("Final Aggregated SNP Missingness DataFrame:")
            print(self.global_missingness_snp_df.head())
            self.global_missingness_snp_df.to_csv('pipeline_final_missingness_snp_results.csv', index=False)

        # Save aggregated individual missingness results if they exist
        if not self.global_missingness_ind_df.empty:
            print("Final Aggregated Individual Missingness DataFrame:")
            print(self.global_missingness_ind_df.head())
            self.global_missingness_ind_df.to_csv('pipeline_final_missingness_ind_results.csv', index=False)

        if not self.global_maf_df.empty:
            print("Final Aggregated MAF DataFrame:")
            print(self.global_maf_df.head())
            self.global_maf_df.to_csv('pipeline_final_maf_results.csv', index=False)

        if not self.global_maf_filter_df.empty:
            print("Final Aggregated Filtered Filter MAF DataFrame:")
            print(self.global_maf_filter_df.head())
            self.global_maf_filter_df.to_csv('pipeline_final_filtered_maf_results.csv', index=False)
        if not self.global_kinship_df.empty:
            print("Final Aggregated Kinship DataFrame:")
            print(self.global_kinship_df.head())
            self.global_kinship_df.to_csv('pipeline_final_kinship_results.csv', index=False)

if __name__ == "__main__":
    num_rounds = 20
    hwe_threshold = 0.05

    strategy = CustomFedAvg(
        hwe_threshold=hwe_threshold,
        fraction_fit=0.1,
        min_fit_clients=1,
        min_available_clients=1,
    )

    # Configure and start the federated server
    server_config = fl.server.ServerConfig(num_rounds=num_rounds)
    fl.server.start_server(server_address="localhost:8080", strategy=strategy, config=server_config)
    strategy.on_conclude()