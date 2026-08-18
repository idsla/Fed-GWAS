# client/local_qc.py
import os
import numpy as np
import logging
from pipeline.clients.base_client import run_plink_command

def exclude_samples_by_missing_rate(plink_prefix, mind_threshold=0.1, new_prefix="filtered_data_by_sample", log_dir="logs"):
    """
    Exclude samples whose missing rate exceeds 'mind_threshold'.
    Produces a new .bed set in 'new_prefix'.

    mind_threshold: float in [0.0, 1.0].
      Example: 0.1 means filter out any sample with >10% missing rate.
    """
    new_prefix = os.path.join(log_dir, new_prefix)
    os.makedirs(log_dir, exist_ok=True)
    cmd = [
        "plink",
        "--bfile", plink_prefix,
        "--mind", str(mind_threshold),
        "--make-bed",
        "--allow-no-sex",
        "--out", new_prefix
    ]
    run_plink_command(cmd)
    return new_prefix

def compute_genotype_counts(plink_prefix, client_id, log_dir="logs", plink_binary=None):
    """
    Compute per-SNP genotype counts [N_AA, N_Aa, N_aa] using PLINK outputs.

    We combine:
      - --freq (allele counts via A1_FREQ * NCHROBS)
      - --hardy (observed heterozygote counts)

    Genotype reconstruction:
      N       = NCHROBS / 2
      het     = observed heterozygotes from .hwe (O(HET))
      a1_cnt  = A1_FREQ * NCHROBS  (approximate integer)
      hom1    = (a1_cnt - het) / 2
      hom2    = N - hom1 - het
    """
    os.makedirs(log_dir, exist_ok=True)

    freq_prefix = os.path.join(log_dir, f"qc_freq_{client_id}")
    hwe_prefix = os.path.join(log_dir, f"qc_hwe_{client_id}")

    # Run PLINK --freq to get allele counts (NCHROBS, A1_FREQ)
    cmd_freq = ["plink", "--bfile", plink_prefix, "--freq", "--out", freq_prefix]
    run_plink_command(cmd_freq, plink_binary=plink_binary)

    # Run PLINK --hardy to get observed heterozygotes
    cmd_hwe = ["plink", "--bfile", plink_prefix, "--hardy", "--out", hwe_prefix]
    run_plink_command(cmd_hwe, plink_binary=plink_binary)

    hwe_file = hwe_prefix + ".hwe"

    if not os.path.exists(hwe_file):
        logging.warning(f"[QC] File not found: {hwe_file} for client {client_id}")
        return np.array([], dtype=np.int64).reshape(0, 3)

    # Parse GENO column directly from .hwe file (format: hom1/het/hom2, e.g., "15/57/175")
    counts = []
    with open(hwe_file, "r") as fin:
        header = fin.readline().strip().split()
        try:
            idx_snp = header.index("SNP")
            idx_geno = header.index("GENO")
            idx_test = header.index("TEST")
        except ValueError as e:
            logging.warning(f"[QC] Unexpected header in {hwe_file}: {header}, err={e}")
            return np.array([], dtype=np.int64).reshape(0, 3)

        for line_num, line in enumerate(fin, start=2):
            parts = line.strip().split()
            if len(parts) <= max(idx_snp, idx_geno, idx_test):
                continue
            try:
                test = parts[idx_test]
                # PLINK may output "ALL" or "ALL(NP)" for combined sample stats
                if not test.startswith("ALL"):
                    continue  # use combined sample stats (ALL or ALL(NP))
                
                # Parse GENO column: format is "hom1/het/hom2" (e.g., "15/57/175")
                geno_str = parts[idx_geno]
                geno_counts = geno_str.split('/')
                if len(geno_counts) != 3:
                    logging.warning(f"[QC] Unexpected GENO format at line {line_num}: {geno_str}")
                    continue
                
                hom1 = int(geno_counts[0])  # N_AA (homozygous A1A1)
                het = int(geno_counts[1])   # N_Aa (heterozygous A1A2)
                hom2 = int(geno_counts[2]) # N_aa (homozygous A2A2)
                
                counts.append([hom1, het, hom2])
            except (ValueError, IndexError) as e:
                logging.warning(f"[QC] Failed to parse line {line_num} in {hwe_file}: {e}, parts={parts}")

    if len(counts) == 0:
        logging.warning(f"[QC] No genotype counts reconstructed for client {client_id}")
    else:
        logging.info(f"[QC] Computed {len(counts)} genotype counts for client {client_id}")

    return np.array(counts, dtype=np.int64)

def compute_missingness_counts(plink_prefix, client_id, log_dir="logs", plink_binary=None):
    """
    Run PLINK to compute the per-SNP missing rate.
    Return Nx2 array of [N_obs, N_miss].
    """
    out_prefix = os.path.join(log_dir, f"mr_{client_id}")
    os.makedirs(log_dir, exist_ok=True)
    cmd = ["plink", "--bfile", plink_prefix, "--missing", "--out", out_prefix]
    run_plink_command(cmd, plink_binary=plink_binary)

    file_name = out_prefix + ".lmiss"
    mr_counts = []
    if not os.path.exists(file_name):
        logging.warning(f"[QC] File not found: {file_name} for client {client_id}")
        logging.warning(f"[QC] PLINK prefix was: {plink_prefix}")
        return np.array([], dtype=np.int64).reshape(0, 2)
    
    with open(file_name, "r") as fin:
        header = fin.readline()  # header
        if not header.strip():
            logging.warning(f"[QC] Empty header in {file_name} for client {client_id}")
            return np.array([], dtype=np.int64).reshape(0, 2)
        
        # Parse header to find column indices (PLINK .lmiss format: CHR SNP N_MISS N_GENO F_MISS)
        header_parts = header.strip().split()
        try:
            idx_chr = header_parts.index('CHR') if 'CHR' in header_parts else 0
            idx_snp = header_parts.index('SNP') if 'SNP' in header_parts else 1
            idx_n_miss = header_parts.index('N_MISS') if 'N_MISS' in header_parts else 2
            idx_n_geno = header_parts.index('N_GENO') if 'N_GENO' in header_parts else 3
        except (ValueError, IndexError):
            # Fallback to fixed positions if header parsing fails
            idx_n_miss = 2
            idx_n_geno = 3
        
        for line_num, line in enumerate(fin, start=2):
            parts = line.strip().split()
            if len(parts) >= max(idx_n_miss + 1, idx_n_geno + 1):
                try:
                    # N_MISS is number of missing genotypes, N_GENO is number of non-missing genotypes
                    N_miss = int(parts[idx_n_miss])
                    N_obs = int(parts[idx_n_geno])
                    mr_counts.append([N_obs, N_miss])
                except (ValueError, IndexError) as e:
                    logging.warning(f"[QC] Failed to parse line {line_num} in {file_name}: {e}, parts={parts}")
            elif line.strip():  # Non-empty line that doesn't have enough parts
                logging.warning(f"[QC] Line {line_num} in {file_name} has only {len(parts)} parts, expected >= {max(idx_n_miss + 1, idx_n_geno + 1)}")
    
    if len(mr_counts) == 0:
        logging.warning(f"[QC] No missingness counts parsed from {file_name} for client {client_id}")
        logging.warning(f"[QC] File exists: {os.path.exists(file_name)}, size: {os.path.getsize(file_name) if os.path.exists(file_name) else 0} bytes")
    
    # Don't delete the file - keep it for debugging
    # os.remove(file_name)
    
    # remove .log
    log_file = out_prefix + ".log"
    if os.path.exists(log_file):
        os.remove(log_file)

    result = np.array(mr_counts, dtype=np.int64)
    logging.info(f"[QC] Computed {len(result)} missingness counts for client {client_id}")
    return result

def run_local_lr(plink_prefix, out_prefix="local_lr", log_dir="logs"):
    """
    Run PLINK logistic regression.
    Return path to the .assoc.logistic file.
    """
    out_prefix = os.path.join(log_dir, out_prefix)
    os.makedirs(log_dir, exist_ok=True)
    cmd = [
        "plink",
        "--bfile", plink_prefix,
        "--logistic",
        "--allow-no-sex",  # Allow analysis when sex information is missing/ambiguous
        "--out", out_prefix
    ]
    run_plink_command(cmd)
    return f"{out_prefix}.assoc.logistic"

def parse_insignificant_snps(assoc_file, p_threshold=1e-3):
    """
    Parse .assoc.logistic, return list of SNP IDs with p-value >= p_threshold.
    """
    # Ensure p_threshold is a float (defensive: config values may be strings)
    p_threshold = float(p_threshold) if p_threshold is not None else 1e-3
    
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

def exclude_snps(plink_prefix, snp_list, new_prefix="filtered_data", log_dir="logs"):
    """
    Exclude the given SNPs from plink_prefix -> produce new_prefix .bed set.
    """
    if not snp_list:
        logging.debug(f"[QC] exclude_snps: empty SNP list; keeping prefix={plink_prefix}")
        return plink_prefix

    logging.info(f"[QC] exclude_snps: excluding {len(snp_list)} SNPs from {plink_prefix}")
    
    os.makedirs(log_dir, exist_ok=True)
    exclude_file = os.path.join(log_dir, "temp_exclude_snps.txt")
    
    try:
        with open(exclude_file, "w") as f:
            for snp in snp_list:
                f.write(f"{snp}\n")
        
        # Verify file was created and has content
        if not os.path.exists(exclude_file):
            logging.error(f"[QC] exclude_snps: file not created: {exclude_file}")
            return plink_prefix
            
        with open(exclude_file, "r") as f:
            content = f.read().strip()
            if not content:
                logging.error(f"[QC] exclude_snps: file is empty: {exclude_file}")
                return plink_prefix
            
    except Exception as e:
        logging.error(f"[QC] exclude_snps: failed to write exclude file: {e}")
        return plink_prefix
    
    new_prefix = os.path.join(log_dir, new_prefix)
    cmd = [
        "plink",
        "--bfile", plink_prefix,
        "--exclude", exclude_file,
        "--make-bed",
        "--allow-no-sex",
        "--out", new_prefix
    ]
    
    logging.debug(f"[QC] exclude_snps: running PLINK command with out={new_prefix}")
    
    try:
        run_plink_command(cmd)
        logging.info(f"[QC] exclude_snps: PLINK output prefix={new_prefix}")
    except Exception as e:
        logging.error(f"[QC] exclude_snps: PLINK command failed: {e}")
        # Clean up and return original prefix
        if os.path.exists(exclude_file):
            os.remove(exclude_file)
        return plink_prefix
    
    if os.path.exists(exclude_file):
        os.remove(exclude_file)
        logging.debug(f"[QC] exclude_snps: removed temp file {exclude_file}")
    
    return new_prefix
