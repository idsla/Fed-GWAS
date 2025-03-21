# client/local_qc.py
import os
import numpy as np
import logging
from .base_client import run_plink_command

def exclude_samples_by_missing_rate(plink_prefix, mind_threshold=0.1, new_prefix="filtered_data_by_sample"):
    """
    Exclude samples whose missing rate exceeds 'mind_threshold'.
    Produces a new .bed set in 'new_prefix'.

    mind_threshold: float in [0.0, 1.0].
      Example: 0.1 means filter out any sample with >10% missing rate.
    """
    cmd = [
        "plink",
        "--bfile", plink_prefix,
        "--mind", str(mind_threshold),
        "--make-bed",
        "--out", new_prefix
    ]
    run_plink_command(cmd)
    return new_prefix

def compute_genotype_counts(plink_prefix, client_id):
    """
    Run PLINK to compute genotype counts.
    Return Nx3 array of [N_AA, N_Aa, N_aa] for each SNP.
    """
    out_prefix = f"qc_counts_{client_id}"
    cmd = ["plink", "--bfile", plink_prefix, "--freq", "counts", "--out", out_prefix]
    run_plink_command(cmd)

    file_name = out_prefix + ".frq.counts"
    counts = []
    if os.path.exists(file_name):
        with open(file_name, "r") as fin:
            fin.readline()  # header
            for line in fin:
                parts = line.strip().split()
                if len(parts) >= 6:
                    try:
                        N_AA = int(parts[3])
                        N_Aa = int(parts[4])
                        N_aa = int(parts[5])
                        counts.append([N_AA, N_Aa, N_aa])
                    except ValueError:
                        pass
        os.remove(file_name)
    return np.array(counts, dtype=np.int64)

def compute_missingness_counts(plink_prefix, client_id):
    """
    Run PLINK to compute the per-SNP missing rate.
    Return Nx2 array of [N_obs, N_miss].
    """
    out_prefix = f"mr_{client_id}"
    cmd = ["plink", "--bfile", plink_prefix, "--missing", "--out", out_prefix]
    run_plink_command(cmd)

    file_name = out_prefix + ".lmiss"
    mr_counts = []
    if os.path.exists(file_name):
        with open(file_name, "r") as fin:
            fin.readline()  # header
            for line in fin:
                parts = line.strip().split()
                if len(parts) >= 6:
                    try:
                        N_obs = int(parts[3])
                        N_miss = int(parts[4])
                        mr_counts.append([N_obs, N_miss])
                    except ValueError:
                        pass
        os.remove(file_name)
    # remove .log
    log_file = out_prefix + ".log"
    if os.path.exists(log_file):
        os.remove(log_file)

    return np.array(mr_counts, dtype=np.int64)

def run_local_lr(plink_prefix, out_prefix="local_lr"):
    """
    Run PLINK logistic regression.
    Return path to the .assoc.logistic file.
    """
    cmd = [
        "plink",
        "--bfile", plink_prefix,
        "--logistic",
        "--out", out_prefix
    ]
    run_plink_command(cmd)
    return f"{out_prefix}.assoc.logistic"

def parse_insignificant_snps(assoc_file, p_threshold=1e-3):
    """
    Parse .assoc.logistic, return list of SNP IDs with p-value >= p_threshold.
    """
    if not os.path.exists(assoc_file):
        return []
    snps = []
    with open(assoc_file, "r") as fin:
        header = fin.readline().strip().split()
        try:
            i_snp = header.index("SNP")
            i_p = header.index("P")
        except ValueError:
            return []
        for line in fin:
            parts = line.strip().split()
            if len(parts) <= max(i_snp, i_p):
                continue
            try:
                pval = float(parts[i_p])
                if pval >= p_threshold:
                    snps.append(parts[i_snp])
            except ValueError:
                pass
    return snps

def exclude_snps(plink_prefix, snp_list, new_prefix="filtered_data"):
    """
    Exclude the given SNPs from plink_prefix -> produce new_prefix .bed set.
    """
    if not snp_list:
        return plink_prefix
    exclude_file = "temp_exclude_snps.txt"
    with open(exclude_file, "w") as f:
        for snp in snp_list:
            f.write(f"{snp}\n")
    cmd = [
        "plink",
        "--bfile", plink_prefix,
        "--exclude", exclude_file,
        "--make-bed",
        "--out", new_prefix
    ]
    run_plink_command(cmd)
    os.remove(exclude_file)
    return new_prefix