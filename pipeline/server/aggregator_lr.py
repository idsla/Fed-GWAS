import os
import uuid
import subprocess
import platform
import shutil
import logging
from pathlib import Path
import numpy as np
import pandas as pd
from flwr.common import parameters_to_ndarrays

logger = logging.getLogger(__name__)

# Cache for PLINK binary path
_PLINK_BINARY = None

def find_plink_binary():
    """Find PLINK binary path"""
    global _PLINK_BINARY
    if _PLINK_BINARY is not None:
        return _PLINK_BINARY
    
    # Find project root - prioritize current working directory
    project_root = None
    cwd = Path.cwd()
    if (cwd / "plink").exists() or (cwd / "experiments").exists() or (cwd / "pyproject.toml").exists():
        project_root = cwd
    
    if project_root is None:
        current_path = Path(__file__).resolve()
        for parent in list(current_path.parents):
            if (parent / "plink").exists() or (parent / "experiments").exists() or (parent / "pyproject.toml").exists():
                project_root = parent
                break
    
    if project_root is None:
        project_root = Path(__file__).parent.parent.parent
    
    # Determine OS-specific PLINK path
    system = platform.system().lower()
    if system == "darwin":
        plink_dir = project_root / "plink" / "plink_mac"
        plink_name = "plink"
    elif system == "linux":
        plink_dir = project_root / "plink" / "plink_linux"
        plink_name = "plink"
    elif system == "windows":
        plink_dir = project_root / "plink" / "plink_win"
        plink_name = "plink.exe"
    else:
        plink_dir = project_root / "plink" / "plink_mac"
        plink_name = "plink"
    
    plink_path = plink_dir / plink_name
    if plink_path.exists() and os.access(plink_path, os.X_OK):
        _PLINK_BINARY = str(plink_path)
        return _PLINK_BINARY
    
    bin_plink = project_root / "bin" / plink_name
    if bin_plink.exists() and os.access(bin_plink, os.X_OK):
        _PLINK_BINARY = str(bin_plink)
        return _PLINK_BINARY
    
    if shutil.which("plink"):
        _PLINK_BINARY = "plink"
        return _PLINK_BINARY
    
    _PLINK_BINARY = str(plink_path)
    return _PLINK_BINARY

def run_server_lr(server_strategy, parameters_list, output_dir=None):
    """
    Server-based LR aggregator with configurable output directory.
    """
    if output_dir is None:
        output_dir = getattr(server_strategy, 'output_dir', '.')
    session_id = uuid.uuid4().hex
    session_path = os.path.join(output_dir, session_id)
    os.makedirs(session_path, exist_ok=True)
    merged_prefix = f"{session_path}/merged_lr"
    bed_files = []

    if not isinstance(parameters_list, list):
        parameters_list = [parameters_list]

    # Process each client's chunk data
    for param in parameters_list:
        nds = parameters_to_ndarrays(param)
        if not nds:
            continue
        chunk_array = nds[0]
            
        if len(chunk_array) == 0:
            continue

        if len(chunk_array) < 3:
            continue
        bed_size = int(chunk_array[0])
        bim_size = int(chunk_array[1])
        fam_size = int(chunk_array[2])
        total_bytes = bed_size + bim_size + fam_size

        # Clients may send uint32 arrays where each element after the 3-size metadata is a single byte (0..255).
        if chunk_array.dtype == np.uint32:
            if len(chunk_array) < 3 + total_bytes:
                continue
            data_u8 = chunk_array[3 : 3 + total_bytes].astype(np.uint8, copy=False)
            chunk_bytes = data_u8.tobytes()
        else:
            data_u8 = chunk_array[3 : 3 + total_bytes].astype(np.uint8, copy=False)
            chunk_bytes = data_u8.tobytes()
        
        # Split the data using the stored sizes
        bed_data = chunk_bytes[:bed_size]
        bim_data = chunk_bytes[bed_size : bed_size + bim_size]
        fam_data = chunk_bytes[bed_size + bim_size : bed_size + bim_size + fam_size]
        
        # Verify we have all the data
        if len(bed_data) != bed_size or len(bim_data) != bim_size or len(fam_data) != fam_size:
            logger.debug(f"[Server LR] Data size mismatch: expected bed={bed_size}, bim={bim_size}, fam={fam_size}, "
                  f"got bed={len(bed_data)}, bim={len(bim_data)}, fam={len(fam_data)}")
            continue
                
        # Create temporary files for this chunk
        chunk_uuid = uuid.uuid4().hex
        bed_file = os.path.join(session_path, f"chunk_{chunk_uuid}")
        
        with open(f"{bed_file}.bed", "wb") as f:
            f.write(bed_data)
        with open(f"{bed_file}.bim", "wb") as f:
            f.write(bim_data)
        with open(f"{bed_file}.fam", "wb") as f:
            f.write(fam_data)
            
        bed_files.append(bed_file)

    if not bed_files:
        logger.warning("[Server LR] No bed files received.")
        return

    first_bed = bed_files[0]
    os.rename(f"{first_bed}.bed", f"{merged_prefix}.bed")
    os.rename(f"{first_bed}.bim", f"{merged_prefix}.bim")
    os.rename(f"{first_bed}.fam", f"{merged_prefix}.fam")
    bed_files = bed_files[1:]

    plink_binary = find_plink_binary()
    for bf in bed_files:
        # First attempt: standard merge
        merge_cmd = [
            plink_binary,
            "--bfile", merged_prefix,
            "--bmerge", bf,
            "--make-bed",
            "--allow-no-sex",
            "--out", f"{merged_prefix}_tmp"
        ]
        
        merge_success = False
        try:
            # Capture output but discard verbose logs (keep stderr for errors)
            # Note: PLINK creates .log file, so we clean it up after
            result = subprocess.run(merge_cmd, check=True, capture_output=True, text=True)
            # Immediately clean up verbose log file if merge succeeds
            log_file = f"{merged_prefix}_tmp.log"
            if os.path.exists(log_file):
                try:
                    os.remove(log_file)
                except OSError:
                    # Ignore file removal errors (file may not exist or be locked)
                    pass
            os.rename(f"{merged_prefix}_tmp.bed", f"{merged_prefix}.bed")
            os.rename(f"{merged_prefix}_tmp.bim", f"{merged_prefix}.bim")
            os.rename(f"{merged_prefix}_tmp.fam", f"{merged_prefix}.fam")
            merge_success = True
        except subprocess.CalledProcessError as e:
            # Check if error is due to strand inconsistency or multiallelic variants
            error_output = e.stderr if e.stderr else (e.stdout if hasattr(e, 'stdout') else "")
            log_file = f"{merged_prefix}_tmp.log"
            missnp_file = f"{merged_prefix}_tmp-merge.missnp"
            
            # Try to handle strand flip issues
            if os.path.exists(missnp_file):
                logger.debug(f"[Server LR] Detected strand inconsistency, attempting --flip")
                try:
                    # First, flip the problematic variants
                    flip_cmd = [
                        plink_binary,
                        "--bfile", bf,
                        "--flip", missnp_file,
                        "--make-bed",
                        "--allow-no-sex",
                        "--out", f"{bf}_flipped"
                    ]
                    subprocess.run(flip_cmd, check=True, capture_output=True, text=True)
                    # Clean up flip log immediately
                    flip_log = f"{bf}_flipped.log"
                    if os.path.exists(flip_log):
                        try:
                            os.remove(flip_log)
                        except OSError:
                            # Ignore file removal errors (file may not exist or be locked)
                            pass
                    
                    # Now try merge again with flipped file
                    merge_cmd_flip = [
                        plink_binary,
                        "--bfile", merged_prefix,
                        "--bmerge", f"{bf}_flipped",
                        "--make-bed",
                        "--allow-no-sex",
                        "--out", f"{merged_prefix}_tmp"
                    ]
                    subprocess.run(merge_cmd_flip, check=True, capture_output=True, text=True)
                    # Immediately clean up merge log
                    merge_log_flip = f"{merged_prefix}_tmp.log"
                    if os.path.exists(merge_log_flip):
                        try:
                            os.remove(merge_log_flip)
                        except OSError:
                            # Ignore file removal errors (file may not exist or be locked)
                            pass
                    os.rename(f"{merged_prefix}_tmp.bed", f"{merged_prefix}.bed")
                    os.rename(f"{merged_prefix}_tmp.bim", f"{merged_prefix}.bim")
                    os.rename(f"{merged_prefix}_tmp.fam", f"{merged_prefix}.fam")
                    merge_success = True
                    logger.debug(f"[Server LR] Merge succeeded after strand flip")
                except subprocess.CalledProcessError:
                    logger.debug(f"[Server LR] Strand flip failed, trying variant ID disambiguation")

            # If strand flip didn't work, try excluding problematic variants on both datasets
            if not merge_success and os.path.exists(missnp_file):
                try:
                    logger.debug(f"[Server LR] Attempting merge after excluding missnp variants")
                    merged_filt = f"{merged_prefix}_filt"
                    bf_filt = f"{bf}_filt"
                    # Exclude missnp from merged dataset
                    subprocess.run(
                        [
                            plink_binary,
                            "--bfile", merged_prefix,
                            "--exclude", missnp_file,
                            "--make-bed",
                            "--allow-no-sex",
                            "--out", merged_filt,
                        ],
                        check=True,
                        capture_output=True,
                        text=True,
                    )
                    # Exclude missnp from incoming chunk
                    subprocess.run(
                        [
                            plink_binary,
                            "--bfile", bf,
                            "--exclude", missnp_file,
                            "--make-bed",
                            "--allow-no-sex",
                            "--out", bf_filt,
                        ],
                        check=True,
                        capture_output=True,
                        text=True,
                    )
                    merge_cmd_excl = [
                        plink_binary,
                        "--bfile", merged_filt,
                        "--bmerge", bf_filt,
                        "--make-bed",
                        "--allow-no-sex",
                        "--out", f"{merged_prefix}_tmp",
                    ]
                    subprocess.run(merge_cmd_excl, check=True, capture_output=True, text=True)
                    os.rename(f"{merged_prefix}_tmp.bed", f"{merged_prefix}.bed")
                    os.rename(f"{merged_prefix}_tmp.bim", f"{merged_prefix}.bim")
                    os.rename(f"{merged_prefix}_tmp.fam", f"{merged_prefix}.fam")
                    merge_success = True
                    logger.debug(f"[Server LR] Merge succeeded after excluding missnp variants")
                except subprocess.CalledProcessError:
                    logger.debug(f"[Server LR] Missnp exclusion merge failed, trying variant ID disambiguation")
                finally:
                    for suffix in [".bed", ".bim", ".fam", ".log"]:
                        for base in [f"{merged_prefix}_filt", f"{bf}_filt"]:
                            path = f"{base}{suffix}"
                            if os.path.exists(path):
                                try:
                                    os.remove(path)
                                except Exception:
                                    pass
            
            # If strand flip didn't work, try making variant IDs unique
            if not merge_success:
                try:
                    logger.debug(f"[Server LR] Attempting to resolve duplicate variant IDs")
                    # Read the BIM file to check for duplicates
                    bim_file = f"{bf}.bim"
                    if os.path.exists(bim_file):
                        # Read BIM file (CHR, SNP, CM, BP, A1, A2)
                        bim_df = pd.read_csv(bim_file, sep='\t', header=None, 
                                            names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
                        
                        # Read merged BIM to check for existing SNP IDs
                        merged_bim = f"{merged_prefix}.bim"
                        if os.path.exists(merged_bim):
                            merged_bim_df = pd.read_csv(merged_bim, sep='\t', header=None,
                                                        names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
                            existing_snps = set(merged_bim_df['SNP'].values)
                            
                            # Make duplicate SNP IDs unique by appending chunk identifier
                            chunk_id = os.path.basename(bf)
                            bim_df['SNP'] = bim_df['SNP'].apply(
                                lambda x: f"{x}_{chunk_id}" if x in existing_snps else x
                            )
                            
                            # Write updated BIM
                            bim_df.to_csv(bim_file, sep='\t', header=False, index=False)
                            
                            # Try merge again with unique IDs
                            merge_cmd_unique = [
                                plink_binary,
                                "--bfile", merged_prefix,
                                "--bmerge", bf,
                                "--make-bed",
                                "--allow-no-sex",
                                "--out", f"{merged_prefix}_tmp"
                            ]
                            subprocess.run(merge_cmd_unique, check=True, capture_output=True)
                            os.rename(f"{merged_prefix}_tmp.bed", f"{merged_prefix}.bed")
                            os.rename(f"{merged_prefix}_tmp.bim", f"{merged_prefix}.bim")
                            os.rename(f"{merged_prefix}_tmp.fam", f"{merged_prefix}.fam")
                            merge_success = True
                            logger.debug(f"[Server LR] Merge succeeded after variant ID disambiguation")
                except Exception as e2:
                    logger.debug(f"[Server LR] Variant ID disambiguation failed: {e2}")
            
            # If all else fails, try excluding problematic variants
            if not merge_success:
                try:
                    logger.debug(f"[Server LR] Attempting merge with --merge-equal-pos flag")
                    merge_cmd_pos = [
                        plink_binary,
                        "--bfile", merged_prefix,
                        "--bmerge", bf,
                        "--merge-equal-pos",
                        "--make-bed",
                        "--allow-no-sex",
                        "--out", f"{merged_prefix}_tmp"
                    ]
                    # Capture output but discard verbose logs (keep stderr for errors)
                    subprocess.run(merge_cmd_pos, check=True, capture_output=True, text=True)
                    # Immediately clean up verbose log file
                    log_file_pos = f"{merged_prefix}_tmp.log"
                    if os.path.exists(log_file_pos):
                        try:
                            os.remove(log_file_pos)
                        except OSError:
                            # Ignore file removal errors (file may not exist or be locked)
                            pass
                    os.rename(f"{merged_prefix}_tmp.bed", f"{merged_prefix}.bed")
                    os.rename(f"{merged_prefix}_tmp.bim", f"{merged_prefix}.bim")
                    os.rename(f"{merged_prefix}_tmp.fam", f"{merged_prefix}.fam")
                    merge_success = True
                    logger.debug(f"[Server LR] Merge succeeded with --merge-equal-pos")
                except subprocess.CalledProcessError:
                    logger.warning(f"[Server LR] All merge attempts failed for chunk {bf}, skipping this chunk")
                    continue
        
        # Clean up temporary files (including verbose PLINK logs)
        # Standard GWAS practice: remove verbose intermediate logs, keep only errors/summaries
        tmp_files_to_remove = [
            f"{merged_prefix}_tmp.log",
            f"{merged_prefix}_tmp-merge.missnp",
            f"{bf}_flipped.bed",
            f"{bf}_flipped.bim",
            f"{bf}_flipped.fam",
            f"{bf}_flipped.log",
            f"{merged_prefix}_tmp-merge.log",  # Additional merge log files
        ]
        for tmp_file in tmp_files_to_remove:
            if os.path.exists(tmp_file):
                try:
                    os.remove(tmp_file)
                except OSError:
                    # Ignore file removal errors (file may not exist or be locked)
                    pass

    # Server runs LR on merged data from clients (clients have already applied QC filtering)
    # No additional filtering needed - clients are responsible for filtering before sending chunks
    lr_cmd = [
        plink_binary,
        "--bfile", merged_prefix,
        "--logistic",
        "--allow-no-sex",  # Allow analysis when sex information is missing/ambiguous
        "--out", f"{session_path}/lr_results"
    ]
    
    try:
        subprocess.run(lr_cmd, check=True, capture_output=True, text=True)
        logger.info("[Server LR] Logistic regression complete.")
    except subprocess.CalledProcessError as e:
        logger.warning(f"[Server LR] PLINK LR failed: {e}")

    # Parse the .assoc.logistic file and return p-values as text.
    assoc_file = f"{session_path}/lr_results.assoc.logistic"
    pvals_text = ""
    if os.path.exists(assoc_file):
        with open(assoc_file, "r") as f:
            header = f.readline()  # skip header
            # The header is expected to contain column names; we assume SNP ID is in the second column and the p-value in the last column
            for line in f:
                parts = line.strip().split()
                if len(parts) < 2:
                    continue
                # Extract SNP ID and p-value
                snp_id = parts[1]
                pval = parts[-1]
                pvals_text += f"{snp_id} {pval}\n"

        # Persist a stable copy of p-values for downstream visualization
        try:
            pvals_path = os.path.join(session_path, "lr_results.pvals.txt")
            with open(pvals_path, "w", encoding="utf-8") as fout:
                fout.write(pvals_text)
        except Exception:
            pass

        # Also write/refresh a stable pointer to the latest session for post-analysis tooling
        try:
            latest_ptr = Path(output_dir) / "latest_lr_session.txt"
            latest_ptr.write_text(session_path, encoding="utf-8")
            # Convenience copy (best-effort)
            latest_assoc = Path(output_dir) / "lr_results_latest.assoc.logistic"
            latest_assoc.write_text(Path(assoc_file).read_text(encoding="utf-8"), encoding="utf-8")
        except Exception:
            pass


    # Convert the concatenated string to a numpy array of uint8
    return np.frombuffer(pvals_text.encode("utf-8"), dtype=np.uint8)

def merge_insign_snp_sets(snp_sets):
    """
    Compute the intersection of sets of insignificant SNPs from clients.
    Returns SNPs that are insignificant in ALL clients.
    """
    if not snp_sets:
        return set()
    if len(snp_sets) == 1:
        return snp_sets[0]
    # Start with first set, then intersect with all others
    merged = snp_sets[0].copy()
    for s in snp_sets[1:]:
        merged = merged.intersection(s)
    return merged
