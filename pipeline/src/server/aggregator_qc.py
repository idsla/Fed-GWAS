# server/aggregator_qc.py

import numpy as np

def secure_sum_arrays(arrays):
    """
    Placeholder for a secure sum aggregator.
    Currently performs a plain sum.
    """
    total = np.zeros_like(arrays[0])
    for arr in arrays:
        total += arr
    return total

def aggregate_global_qc(server_strategy, partial_data, config):
    """
    Aggregate partial QC data from clients.
    Each client sends [counts_array, missing_array, thresh_array] where:
      - counts_array: shape (n_snps, 3) for [N_AA, N_Aa, N_aa]
      - missing_array: shape (n_snps, 2) for [N_obs, N_miss]
      - thresh_array: shape (3,) for [maf_threshold, missing_threshold, hwe_threshold]
    Assumes SNP ordering is identical across clients.

    Returns a single set of SNP indices to exclude.
    """
    if not partial_data:
        print("[Server QC] No partial data received.")
        return set()

    all_counts = []
    all_missing = []
    all_thresholds = []
    for item in partial_data:
        counts_arr = item[0]
        missing_arr = item[1]
        thresh_arr = item[2]
        all_counts.append(counts_arr)
        all_missing.append(missing_arr)
        all_thresholds.append(thresh_arr)

    thresh_array = np.stack(all_thresholds)  # shape (n_clients, 3)
    # We use the min MAF threshold, max missing threshold, and min HWE threshold
    maf_final = float(np.min(thresh_array[:, 0]))
    miss_final = float(np.max(thresh_array[:, 1]))
    hwe_final = float(np.min(thresh_array[:, 2]))
    print(f"[Server QC] Unified thresholds => MAF={maf_final}, Missing={miss_final}, HWE={hwe_final}")

    n_snps = all_counts[0].shape[0]
    counts_sum = secure_sum_arrays([arr for arr in all_counts])  # shape (n_snps, 3)
    missing_sum = secure_sum_arrays([arr for arr in all_missing])  # shape (n_snps, 2)

    exclude_indices = set()

    for i in range(n_snps):

        # MAF
        N_AA, N_Aa, N_aa = counts_sum[i]
        N_obs, N_miss = missing_sum[i]

        total_geno = N_AA + N_Aa + N_aa
        if total_geno == 0:
            exclude_indices.add(i)
            continue

        freqA = (2*N_AA + N_Aa)/(2*total_geno)
        maf = min(freqA, 1.0 - freqA)
        if maf < maf_final:
            exclude_indices.add(i)
            continue

        # Missing rate (per-SNP)
        denom = N_obs + N_miss
        if denom == 0:
            exclude_indices.add(i)
            continue
        missing_rate = N_miss / denom
        if missing_rate > miss_final:
            exclude_indices.add(i)
            continue

        # chi-squared based HWE
        p = freqA
        q = 1.0 - p
        E_AA = p*p*total_geno
        E_Aa = 2.0*p*q*total_geno
        E_aa = q*q*total_geno
        chi_sq = 0.0
        for obs, exp in [(N_AA, E_AA), (N_Aa, E_Aa), (N_aa, E_aa)]:
            if exp > 1e-9:
                chi_sq += (obs - exp)**2 / exp
        pval_hwe = 1 - chi2_cdf(chi_sq, df=1)
        if pval_hwe < hwe_final:
            exclude_indices.add(i)
            continue

    return exclude_indices