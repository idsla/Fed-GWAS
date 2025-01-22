import flwr as fl
import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed
from pyplink import PyPlink
import json
import logging
import datetime
logging.basicConfig(filename='king_analysis.log', level=logging.INFO, format='%(message)s')
start_time = datetime.datetime.now()
# Hardy-Weinberg Equilibrium (HWE) function
def snphwe(obs_hets: int, obs_hom1: int, obs_hom2: int) -> float:
    if obs_hom1 < 0 or obs_hom2 < 0 or obs_hets < 0:
        raise ValueError("Observed counts must be non-negative integers.")
    N = obs_hom1 + obs_hom2 + obs_hets
    rare = min(obs_hom1, obs_hom2) * 2 + obs_hets
    probs = np.zeros(rare + 1)
    mid = rare * (2 * N - rare) // (2 * N)
    if (mid % 2) != (rare % 2):
        mid += 1
    probs[mid] = 1.0
    mysum = 1.0
    curr_hets = mid
    curr_homr = (rare - mid) // 2
    curr_homc = N - curr_hets - curr_homr
    while curr_hets >= 2:
        probs[curr_hets - 2] = probs[curr_hets] * curr_hets * (curr_hets - 1) / (4 * (curr_homr + 1) * (curr_homc + 1))
        mysum += probs[curr_hets - 2]
        curr_hets -= 2
        curr_homr += 1
        curr_homc += 1
    curr_hets = mid
    curr_homr = (rare - mid) // 2
    curr_homc = N - curr_hets - curr_homr
    while curr_hets <= rare - 2:
        probs[curr_hets + 2] = probs[curr_hets] * 4 * curr_homr * curr_homc / ((curr_hets + 2) * (curr_hets + 1))
        mysum += probs[curr_hets + 2]
        curr_hets += 2
        curr_homr -= 1
        curr_homc -= 1
    target = probs[obs_hets]
    p_value = min(1.0, np.sum(probs[probs <= target]) / mysum)
    return p_value


class GenotypeClient(fl.client.NumPyClient):
    def __init__(self, bed_file, num_rounds=5,kinship_thresholds= None, threshold=0.05, modes=None):
        self.bed_file = bed_file
        self.num_rounds = num_rounds
        self.threshold = threshold
        self.kinship_thresholds = kinship_thresholds or {'firstkin': 0.35, 'secondkin': 0.18, 'thirdkin':  0.09}
        self.modes = modes or []
        self.round_counter = {mode: 0 for mode in self.modes}
        self.current_mode_index = 0
        #self.current_mode_round = 0
        # Initialize SNP data if any mode is selected
        if "geno" in self.modes or "hwd" in self.modes or "missingness" in self.modes or "calculate_missing_rate" in self.modes or "calculate_maf" in self.modes or "filter_maf" in self.modes or "kinship":
            snp_reader = Bed(bed_file)
            self.num_snps = snp_reader.sid_count
            self.snps_per_round = self.num_snps // self.num_rounds
            self.snp_start = 0
            self.snp_end = 0
            self.current_round = 0

        if "check_sex" in self.modes or "missingness" in self.modes or "calculate_missing_rate" in self.modes or "calculate_maf" in self.modes or "filter_maf" in self.modes:
            self.bed = PyPlink(bed_file)
            fam_data = self.bed.get_fam()  # Load FAM data to determine total number of individuals
            self.num_individuals = len(fam_data)
            self.individuals_per_round = self.num_individuals // self.num_rounds
            self.current_round = 0
    def calculate_kinship(self, genotype_i, genotype_j):
        """
        Calculate the KING kinship coefficient between two individuals.
        """
        n11 = np.sum((genotype_i == 1) & (genotype_j == 1))
        n02 = np.sum((genotype_i == 2) & (genotype_j == 0))
        n20 = np.sum((genotype_i == 0) & (genotype_j == 2))
        n1_s = np.sum(genotype_i == 1)
        s_1 = np.sum(genotype_j == 1)

        if n1_s == 0:
            return 0.0

        phi_ij = (2 * n11 - 4 * (n02 + n20) - n1_s + s_1) / (4 * n1_s)
        return phi_ij
    def  kinship_analysis(self):
        """
        Perform kinship analysis on the current subset of SNPs.
        """
        
        self.snp_start = self.current_round * self.snps_per_round
        self.snp_end = self.snp_start + self.snps_per_round
        if self.current_round == self.num_rounds - 1:
            self.snp_end = self.num_snps

        # Subset genotype data for the current round
        snp_reader = Bed(self.bed_file, count_A1=False)
        partition = snp_reader[:, self.snp_start:self.snp_end]
        print(f"Partition shape: {partition.shape}")
        print(f"Partition values: {partition.read().val}")
        print(f"Total SNPs in file: {snp_reader.sid_count}")
        print(f"SNP range requested: {self.snp_start} to {self.snp_end}")
        print(f"snp start {self.snp_start} and snp end {self.snp_end}")
        genotype_matrix = partition.read().val
        
        genotype_matrix = np.nan_to_num(genotype_matrix, nan=-1).astype(np.int8)
        print(f"genotype_matrix {genotype_matrix}")
        genotype_subset = genotype_matrix
        print(f"genotype_subset {genotype_subset}")
        num_individuals = genotype_subset.shape[0]

        # Calculate pairwise kinship coefficients
        kinship_results = {}
        for i in range(num_individuals):
            for j in range(i + 1, num_individuals):
                phi_ij = self.calculate_kinship(genotype_subset[i, :], genotype_subset[j, :])
                pair_key = f"{i}_{j}"  # Convert tuple (i, j) to string
                kinship_results[pair_key] = phi_ij

        # Classify kinship degree
        classified_results = {}
        for pair, value in kinship_results.items():
            logging.info(f" Final Estimate {pair} and {value} and threshold {self.kinship_thresholds['firstkin']}")
            if value > self.kinship_thresholds['firstkin']:
                classified_results[pair] = '1st degree'
                logging.info(f"1st degree Final Estimate {pair} and {value}")
            elif value > self.kinship_thresholds['secondkin']:
                classified_results[pair] = '2nd degree'
            elif value > self.kinship_thresholds['thirdkin']:
                classified_results[pair] = '3rd degree'
            else:
                classified_results[pair] = 'unrelated'

        # Increment the round counter
        self.current_round += 1
        return classified_results
    def count_genotype(self, genotype_matrix):
        genotype_counts = {}
        for snp_index in range(genotype_matrix.shape[1]):
            genotypes = genotype_matrix[:, snp_index]
            unique_value, counts = np.unique(genotypes, return_counts=True)
            counts_ = dict(zip(unique_value, counts))
            genotype_counts[snp_index] = counts_
        
        genotype_df = pd.DataFrame.from_dict(genotype_counts, orient='index').fillna(0)
        return genotype_df
    def geno(self):
        print("Running geno function")
        try:
            self.snp_start = self.current_round * self.snps_per_round
            self.snp_end = self.snp_start + self.snps_per_round

            if self.current_round == self.num_rounds - 1:
                self.snp_end = self.num_snps

            snp_reader = Bed(self.bed_file, count_A1=False)
            partition = snp_reader[:, self.snp_start:self.snp_end]
            genotype_matrix = partition.read().val
            genotype_matrix = np.nan_to_num(genotype_matrix, nan=-1).astype(np.int8)

            bed = PyPlink(self.bed_file)
            bim_df = pd.DataFrame(bed.get_bim().reset_index())
            bim_df = bim_df.iloc[self.snp_start:self.snp_end].reset_index(drop=True)

            n_individuals, n_snps = genotype_matrix.shape
            self.snp_missing_proportions = np.sum(np.isnan(genotype_matrix), axis=0) / n_individuals
            keep_snps = self.snp_missing_proportions <= self.threshold
            genotype_matrix_filtered = genotype_matrix[:, keep_snps]
            filtered_bim = bim_df[keep_snps].reset_index(drop=True)

            genotype_df = self.count_genotype(genotype_matrix)
            snps, nchrobs = [], []
            for snp, genotype in genotype_df.iterrows():
                no_of_chzms = 2 * (genotype.get(0, 0) + genotype.get(1, 0) + genotype.get(2, 0))
                nchrobs.append(no_of_chzms)
                snps.append(snp)

            df = pd.DataFrame({
                'CHR': bim_df['chrom'],
                'SNP': snps,
                'POS': bim_df['pos'],
                'NCHROBS': nchrobs,
            })

            metrics = {
                'geno_results': True,  # Flag for server to detect geno results
                'CHR': json.dumps(df['CHR'].tolist()),
                'SNP': json.dumps(df['SNP'].tolist()),
                'POS': json.dumps(df['POS'].tolist()),
                'NCHROBS': json.dumps(df['NCHROBS'].tolist()),
                'MISSINGNESS': json.dumps(self.snp_missing_proportions.tolist()),
                'num_individuals': genotype_matrix.shape[0]
            }
            self.current_round += 1
            return [], len(filtered_bim), metrics

        except Exception as e:
            print(f"Error in geno function: {e}")
            return [], 0, {}
    def hwd(self):
        """Executes the Hardy-Weinberg Disequilibrium test and filters SNPs based on threshold."""
        print("Running hwd function")
        try:
            # Calculate start and end indices for current round's SNPs
            self.snp_start = self.current_round * self.snps_per_round
            self.snp_end = self.snp_start + self.snps_per_round

            # Ensure that the last round processes any remaining SNPs
            if self.current_round == self.num_rounds - 1:
                self.snp_end = self.num_snps

            # Load genotype data for the current partition of SNPs
            snp_reader = Bed(self.bed_file, count_A1=False)
            partition = snp_reader[:, self.snp_start:self.snp_end]
            genotype_matrix = partition.read().val
            genotype_matrix = np.nan_to_num(genotype_matrix, nan=-1).astype(np.int8)

            # Load BIM file to get SNP identifiers
            bed = PyPlink(self.bed_file)
            bim_df = pd.DataFrame(bed.get_bim().reset_index())
            bim_df = bim_df.iloc[self.snp_start:self.snp_end].reset_index(drop=True)

            # Perform HWE test on each SNP
            p_values = []
            snps = []
            for i in range(genotype_matrix.shape[1]):
                snp_data = genotype_matrix[:, i]
                obs_hom1 = np.sum(snp_data == 0)
                obs_hets = np.sum(snp_data == 1)
                obs_hom2 = np.sum(snp_data == 2)

                # Skip SNPs with no observations (all missing)
                if obs_hom1 + obs_hets + obs_hom2 == 0:
                    continue

                # Calculate HWE p-value
                p_value = snphwe(obs_hets, obs_hom1, obs_hom2)
                if p_value >= self.threshold:
                    p_values.append(p_value)
                    snps.append(bim_df.iloc[i]['snp'])  # Use SNP ID from BIM file

            # Increment the round counter for the next call
            self.current_round += 1

            # Package metrics
            metrics = {
                'hwd_results': True,
                'SNP': json.dumps(snps),
                'p_value': json.dumps(p_values),
            }
            return [], len(snps), metrics

        except Exception as e:
            print(f"Error in hwd function: {e}")
            return [], 0, {}


    def calculate_snpsex(self, genotypes):
        """
        Calculate the homozygosity rate based on the provided genotypes.

        Parameters:
        - genotypes: An array of genotypes for an individual.

        Returns:
        - The homozygosity rate as a float.
        """
        # Count homozygous genotypes (0/0 and 2/2)
        homozygous = sum(genotypes == 0) + sum(genotypes == 2)
        # Calculate the total number of valid genotypes (ignoring missing data)
        total_valid = len(genotypes[genotypes >= 0])
        # Compute the homozygosity rate
        homo_rate = homozygous / total_valid if total_valid > 0 else 0
        return homo_rate
    
    def check_sex(self):
        """Executes the check sex functionality."""
        print("Running check sex function")
        # Calculate start and end indices for current round's subset
        start_idx = self.current_round * self.individuals_per_round
        end_idx = start_idx + self.individuals_per_round

        # Ensure that the last round processes any remaining individuals
        if self.current_round == self.num_rounds - 1:
            end_idx = self.num_individuals

        # Load FAM data and filter to the current subset of individuals
        fam_df = pd.DataFrame(self.bed.get_fam(), columns=['fid', 'iid', 'father', 'mother', 'gender', 'status'])
        subset_fam_df = fam_df.iloc[start_idx:end_idx].reset_index(drop=True)

        # Load BIM data to get X chromosome SNPs
        bim_df = pd.DataFrame(self.bed.get_bim()).reset_index()
        x_chromosome_snps = bim_df[bim_df['chrom'] == 23]['snp'].dropna()
        
        results = []
        for index, row in subset_fam_df.iterrows():
            observed_homoz = []
            for snp in x_chromosome_snps:
                genotype = self.bed.get_geno_marker(snp)[index]
                observed_homoz.append(genotype)

            # Calculate homozygosity and infer sex
            observed_homozygosity = self.calculate_snpsex(np.array(observed_homoz))
            snpsexBn = 1 if observed_homozygosity > 0.8 else 2
            status = 'OK' if snpsexBn == row['gender'] and snpsexBn != 0 else 'PROBLEM'
            F = (observed_homozygosity - 0.5) / 0.5 if observed_homozygosity > 0.5 else 0

            results.append({
                'FID': row['fid'],
                'IID': row['iid'],
                'PEDSEX': row['gender'],
                'SNPSEX': snpsexBn,
                'STATUS': status,
                'F': F
            })

        # Increment the round counter for the next call
        self.current_round += 1

        # Package metrics
        metrics = {
            'sex_check_results': json.dumps(results),
            'num_individuals': len(results)
        }
        return [], len(results), metrics

    def filter_missingness_samples(self):
        """Filter individuals based on missing genotype data."""
        print("Running missingness filtering function")
        try:
            # Calculate start and end indices for the current round's individuals
            start_idx = self.current_round * self.individuals_per_round
            end_idx = start_idx + self.individuals_per_round

            # Ensure that the last round processes any remaining individuals
            if self.current_round == self.num_rounds - 1:
                end_idx = self.num_individuals

            # Load genotype data for the current partition of individuals
            snp_reader = Bed(self.bed_file, count_A1=False)
            partition = snp_reader[start_idx:end_idx, :]
            genotype_matrix = partition.read().val

            # Calculate the missingness for each individual
            n_individuals, n_snps = genotype_matrix.shape
            missing_per_ind = np.sum(np.isnan(genotype_matrix), axis=1) / n_snps
            keep_individuals = missing_per_ind <= self.threshold

            # Load and filter FAM data based on missingness
            fam_df = pd.DataFrame(self.bed.get_fam(), columns=['fid', 'iid', 'father', 'mother', 'gender', 'status'])
            subset_fam_df = fam_df.iloc[start_idx:end_idx].reset_index(drop=True)
            filtered_fam = subset_fam_df[keep_individuals].reset_index(drop=True)
            individuals_removed = len(subset_fam_df) - len(filtered_fam)

            # Increment the round counter for the next call
            self.current_round += 1

            # Package metrics
            metrics = {
                'missingness_results': json.dumps(filtered_fam.to_dict(orient='records')),
                'individuals_removed': int(individuals_removed),
                'remaining_individuals': int(keep_individuals.sum())
            }
            return [], len(filtered_fam), metrics

        except Exception as e:
            print(f"Error in filter_missingness_samples function: {e}")
            return [], 0, {}
    def calculate_missing_rate(self):
        """
        Calculate missing rates for both SNPs and individuals with partition logic.
        """
        print("Running calculate_missing_rate function")
        try:
            # Calculate start and end indices for the current round's partition of individuals
            start_idx = self.current_round * self.individuals_per_round
            end_idx = start_idx + self.individuals_per_round

            # Ensure the last round processes any remaining individuals
            if self.current_round == self.num_rounds - 1:
                end_idx = self.num_individuals

            # Load genotype data for the current partition of individuals
            snp_reader = Bed(self.bed_file, count_A1=False)
            partition = snp_reader[start_idx:end_idx, :]
            genotype_matrix = partition.read().val

            # Calculate missing rates for SNPs
            n_individuals, n_snps = genotype_matrix.shape
            missing_per_snp = np.sum(np.isnan(genotype_matrix), axis=0) / n_individuals
            snp_data = pd.DataFrame({
                'SNP': [f"SNP_{i}" for i in range(n_snps)],  # SNP IDs from the Bed partition
                'N_MISS': np.sum(np.isnan(genotype_matrix), axis=0),
                'N_GENO': n_individuals,
                'F_MISS': missing_per_snp
            })

            # Calculate missing rates for individuals
            missing_per_ind = np.sum(np.isnan(genotype_matrix), axis=1) / n_snps
            fam_df = pd.DataFrame(self.bed.get_fam(), columns=['FID', 'IID', 'father', 'mother', 'gender', 'status'])
            subset_fam_df = fam_df.iloc[start_idx:end_idx].reset_index(drop=True)
            ind_data = pd.DataFrame({
                'FID': [f"FID_{i}" for i in range(n_individuals)],
                'IID': [f"IID_{i}" for i in range(n_individuals)],
                'N_MISS': np.sum(np.isnan(genotype_matrix), axis=1),
                'N_GENO': [n_snps] * n_individuals,
                'F_MISS': missing_per_ind
            })

            # Increment the round counter for the next call
            self.current_round += 1

            # Package metrics with Python-native types for serialization
            metrics = {
                'missing_rate_snp': json.dumps(snp_data.to_dict(orient='records')),
                'missing_rate_ind': json.dumps(ind_data.to_dict(orient='records')),
                'num_snps': n_snps,
                'num_individuals': int(n_individuals),
            }
            return [], len(snp_data), metrics

        except Exception as e:
            print(f"Error in calculate_missing_rate function: {e}")
            return [], 0, {}
    def calculate_maf(self):
        """
        Calculate MAF (Minor Allele Frequency) with partition logic for SNPs.
        """
        print("Running calculate_maf function")
        try:
            # Calculate start and end indices for the current round's SNPs
            self.snp_start = self.current_round * self.snps_per_round
            self.snp_end = self.snp_start + self.snps_per_round

            # Ensure the last round processes any remaining SNPs
            if self.current_round == self.num_rounds - 1:
                self.snp_end = self.num_snps

            # Load genotype data for the current partition of SNPs
            snp_reader = Bed(self.bed_file, count_A1=False)
            partition = snp_reader[:, self.snp_start:self.snp_end]
            genotype_matrix = partition.read().val
            genotype_matrix = np.nan_to_num(genotype_matrix, nan=-1).astype(np.int8)

            # Count genotype data
            genotype_df = self.count_genotype(genotype_matrix)

            # Load BIM data for SNP information
            bim_df = pd.DataFrame(self.bed.get_bim().reset_index())
            bim_df = bim_df.iloc[self.snp_start:self.snp_end].reset_index(drop=True)

            snps, mafs, nchrobs = [], [], []
            for snp, genotype in genotype_df.iterrows():
                no_of_chzms = 2 * (genotype.get(0, 0) + genotype.get(1, 0) + genotype.get(2, 0))

                # Allele counts
                count_a1 = 2 * genotype.get(0, 0) + genotype.get(1, 0)
                count_a2 = 2 * genotype.get(2, 0) + genotype.get(1, 0)

                # Calculate allele frequencies
                freq_a1 = count_a1 / no_of_chzms
                freq_a2 = count_a2 / no_of_chzms

                # MAF calculation
                maf_coeff = min(freq_a1, freq_a2)
                mafs.append(maf_coeff)
                nchrobs.append(no_of_chzms)
                snps.append(snp)

            # Create DataFrame for MAF results
            maf_df = pd.DataFrame({
                'CHR': bim_df['chrom'],
                'SNP': snps,
                'POS': bim_df['pos'],
                'A1': bim_df['a1'],
                'A2': bim_df['a2'],
                'MAF': mafs,
                'NCHROBS': nchrobs,
            })

            # Increment the round counter for the next call
            self.current_round += 1

            # Package metrics for the server
            metrics = {
                'maf_results': json.dumps(maf_df.to_dict(orient='records'))
            }
            return [], len(maf_df), metrics

        except Exception as e:
            print(f"Error in calculate_maf function: {e}")
            return [], 0, {}
    def filter_maf_variants(self, maf_min=None, maf_max=0.5, mac_min=None, mac_max=None):
        """
        Filter SNPs based on MAF and MAC thresholds with partition logic.

        Parameters:
        - maf_min (float): Minimum MAF threshold.
        - maf_max (float): Maximum MAF threshold.
        - mac_min (int): Minimum MAC threshold.
        - mac_max (int): Maximum MAC threshold.

        Returns:
        - Tuple: Filtered SNP DataFrame and metrics for server aggregation.
        """
        print("Running filter_maf_variants function")
        try:
            # Calculate start and end indices for the current round's SNPs
            self.snp_start = self.current_round * self.snps_per_round
            self.snp_end = self.snp_start + self.snps_per_round

            # Ensure the last round processes any remaining SNPs
            if self.current_round == self.num_rounds - 1:
                self.snp_end = self.num_snps

            # Load genotype data for the current partition of SNPs
            snp_reader = Bed(self.bed_file, count_A1=False)
            partition = snp_reader[:, self.snp_start:self.snp_end]
            genotype_matrix = partition.read().val
            genotype_matrix = np.nan_to_num(genotype_matrix, nan=-1).astype(np.int8)

            # Count genotype data
            genotype_df = self.count_genotype(genotype_matrix)

            # Load BIM data for SNP information
            bim_df = pd.DataFrame(self.bed.get_bim().reset_index())
            bim_df = bim_df.iloc[self.snp_start:self.snp_end].reset_index(drop=True)

            # Calculate MAF and NCHROBS
            snps, mafs, nchrobs = [], [], []
            for snp, genotype in genotype_df.iterrows():
                no_of_chzms = 2 * (genotype.get(0, 0) + genotype.get(1, 0) + genotype.get(2, 0))
                count_a1 = 2 * genotype.get(0, 0) + genotype.get(1, 0)
                count_a2 = 2 * genotype.get(2, 0) + genotype.get(1, 0)
                freq_a1 = count_a1 / no_of_chzms
                freq_a2 = count_a2 / no_of_chzms
                maf_coeff = min(freq_a1, freq_a2)
                mafs.append(maf_coeff)
                nchrobs.append(no_of_chzms)
                snps.append(snp)

            # Create DataFrame for MAF and NCHROBS
            maf_df = pd.DataFrame({
                'CHR': bim_df['chrom'],
                'SNP': snps,
                'POS': bim_df['pos'],
                'A1': bim_df['a1'],
                'A2': bim_df['a2'],
                'MAF': mafs,
                'NCHROBS': nchrobs,
            })

            # Apply MAF and MAC filtering
            if maf_min is not None:
                maf_df = maf_df[maf_df['MAF'] >= maf_min]
            if maf_max is not None:
                maf_df = maf_df[maf_df['MAF'] <= maf_max]
            if mac_min is not None:
                maf_df['MAC'] = maf_df.apply(lambda row: min(row['NCHROBS'] * row['MAF'], (1 - row['MAF']) * row['NCHROBS']), axis=1)
                maf_df = maf_df[maf_df['MAC'] >= mac_min]
            if mac_max is not None:
                if 'MAC' not in maf_df.columns:
                    maf_df['MAC'] = maf_df.apply(lambda row: min(row['NCHROBS'] * row['MAF'], (1 - row['MAF']) * row['NCHROBS']), axis=1)
                maf_df = maf_df[maf_df['MAC'] <= mac_max]

            # Increment the round counter for the next call
            self.current_round += 1

            # Package metrics for the server
            metrics = {
                'filtered_maf_results': json.dumps(maf_df.to_dict(orient='records')),
                'num_filtered_snps': len(maf_df),
            }
            return [], len(maf_df), metrics

        except Exception as e:
            print(f"Error in filter_maf_variants function: {e}")
            return [], 0, {}

    def fit(self, parameters, config):
        metrics = {}

        # Ensure all modes have been processed
        if self.current_mode_index >= len(self.modes):
            print("All modes completed.")
            return [], 1, metrics  # No more modes to execute

        # Get the current mode and round counter
        current_mode = self.modes[self.current_mode_index]
        current_round = self.round_counter[current_mode]
        print(f"Running round {current_round + 1} for mode {current_mode}")
        print(f"Current round: {current_round + 1}, Total rounds: {self.num_rounds}")

        # Execute the current mode's function
        try:
            if current_mode == "geno":
                _, _, geno_metrics = self.geno()
                metrics.update(geno_metrics)
            elif current_mode == "check_sex":
                _, _, sex_metrics = self.check_sex()
                metrics.update(sex_metrics)
            elif current_mode == "hwd":
                _, _, hwd_metrics = self.hwd()
                metrics.update(hwd_metrics)
            elif current_mode == "missingness":
                _, _, missingness_metrics = self.filter_missingness_samples()
                metrics.update(missingness_metrics)
            elif current_mode == "calculate_missing_rate":
                _, _, missing_rate_metrics = self.calculate_missing_rate()
                metrics.update(missing_rate_metrics)
            elif current_mode == "calculate_maf":
                _, _, maf_metrics = self.calculate_maf()
                metrics.update(maf_metrics)
            elif current_mode == "filter_maf":
                _, _, maf_filter_metrics = self.filter_maf_variants(
                    maf_min=config.get("maf_min", None),
                    maf_max=config.get("maf_max", None),
                    mac_min=config.get("mac_min", None),
                    mac_max=config.get("mac_max", None)
                )
                metrics.update(maf_filter_metrics)
            elif current_mode == "kinship":
                kinship_results = self.kinship_analysis()
                metrics["kinship_results"] = json.dumps(kinship_results)
        except Exception as e:
            print(f"Error in {current_mode} function: {e}")
            metrics[f"{current_mode}_error"] = str(e)

        # Increment the round counter for the current mode
        self.round_counter[current_mode] += 1

        # Check if the current mode has completed all rounds
        if self.round_counter[current_mode] >= self.num_rounds:
            print(f"Completed all rounds for mode: {current_mode}")
            self.current_mode_index += 1  # Move to the next mode

            # Check if all modes are completed
            if self.current_mode_index < len(self.modes):
                next_mode = self.modes[self.current_mode_index]
                print(f"Switched to mode: {next_mode}")
                if next_mode in ["kinship","hwd", "geno", "calculate_maf", "filter_maf"]:
                    self.snp_start = 0
                    self.snp_end = 0
                    self.current_round = 0
                if next_mode in ["check_sex", "missingness", "calculate_missing_rate"]:
                    self.current_round = 0
            # Reset the round counter for the next mode
            #next_mode = self.modes[self.current_mode_index]
            #self.round_counter[next_mode] = 0  # Reset rounds for the new mode
            #print(f"Switched to mode: {next_mode}")

        return [], 1, metrics

    def evaluate(self, parameters, config):
        return 0.0, 1, {}
    
    

if __name__ == "__main__":
    print("Select mode:\n1. Geno\n2. Check Sex\n3. HWD\n4. Missingness Filtering\n5. Calculate Missingness Sample\n6. Calculate MAF\n7. Filter MAF\n8. Kinship")
    choice = input("Enter choice (1, 2, 3, 4, 5, 6, 7, 8 or combination separated by space): ")

    modes = []
    if '1' in choice:
        modes.append("geno")
    if '2' in choice:
        modes.append("check_sex")
    if '3' in choice:
        modes.append("hwd")
    if '4' in choice:
        modes.append("missingness")
    if '5' in choice:
        modes.append("calculate_missing_rate")
    if '6' in choice:
        modes.append("calculate_maf")
    if '7' in choice:
        modes.append("filter_maf")
    if '8' in choice:
        modes.append("kinship")


    client1 = GenotypeClient(bed_file="data/begin.cc", modes=modes)
    fl.client.start_numpy_client(server_address="localhost:8080", client=client1)

