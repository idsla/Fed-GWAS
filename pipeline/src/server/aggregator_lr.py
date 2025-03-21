import os
import subprocess
import tarfile
import uuid
import numpy as np

def run_server_lr(server_strategy, parameters_list):
    """
    Server-based LR aggregator:
    1. For each parameter (tar file bytes) from clients, save as a tar file.
    2. Unpack each tar file, merge resulting PLINK binary files using --bmerge.
    3. Run PLINK logistic regression on the merged dataset.
    4. Parse the .assoc.logistic file and return p-values as text (e.g., lines of "anonSNPID p_value").
    """
    session_id = uuid.uuid4().hex
    os.makedirs(session_id, exist_ok=True)
    merged_prefix = f"{session_id}/merged_lr"
    bed_files = []

    for param in parameters_list:
        chunk_bytes = param.numpy.tobytes()
        tar_path = os.path.join(session_id, f"{uuid.uuid4().hex}.tar")
        with open(tar_path, "wb") as f:
            f.write(chunk_bytes)
        with tarfile.open(tar_path, "r") as tf:
            tf.extractall(session_id)
        chunk_uuid = uuid.uuid4().hex
        bed_file = os.path.join(session_id, f"chunk_{chunk_uuid}")
        os.rename(os.path.join(session_id, "chunk.bed"), f"{bed_file}.bed")
        os.rename(os.path.join(session_id, "chunk.bim"), f"{bed_file}.bim")
        os.rename(os.path.join(session_id, "chunk.fam"), f"{bed_file}.fam")
        bed_files.append(bed_file)
        os.remove(tar_path)

    if not bed_files:
        print("[Server LR] No bed files received.")
        return

    first_bed = bed_files[0]
    os.rename(f"{first_bed}.bed", f"{merged_prefix}.bed")
    os.rename(f"{first_bed}.bim", f"{merged_prefix}.bim")
    os.rename(f"{first_bed}.fam", f"{merged_prefix}.fam")
    bed_files = bed_files[1:]

    for bf in bed_files:
        merge_cmd = [
            "plink",
            "--bfile", merged_prefix,
            "--bmerge", f"{bf}.bed", f"{bf}.bim", f"{bf}.fam",
            "--make-bed",
            "--out", f"{merged_prefix}_tmp"
        ]
        try:
            subprocess.run(merge_cmd, check=True)
            os.rename(f"{merged_prefix}_tmp.bed", f"{merged_prefix}.bed")
            os.rename(f"{merged_prefix}_tmp.bim", f"{merged_prefix}.bim")
            os.rename(f"{merged_prefix}_tmp.fam", f"{merged_prefix}.fam")
            if os.path.exists(f"{merged_prefix}_tmp.log"):
                os.remove(f"{merged_prefix}_tmp.log")
        except subprocess.CalledProcessError as e:
            print(f"[Server LR] Merge failed: {e}")

    lr_cmd = [
        "plink",
        "--bfile", merged_prefix,
        "--logistic",
        "--out", f"{session_id}/lr_results"
    ]
    try:
        subprocess.run(lr_cmd, check=True)
        print("[Server LR] Logistic regression complete.")
    except subprocess.CalledProcessError as e:
        print(f"[Server LR] PLINK LR failed: {e}")

    # Parse the .assoc.logistic file and return p-values as text.
    assoc_file = f"{session_id}/lr_results.assoc.logistic"
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

    # Convert the concatenated string to a numpy array of uint8
    return np.frombuffer(pvals_text.encode("utf-8"), dtype=np.uint8)

def merge_insign_snp_sets(snp_sets):
    """
    Placeholder: compute the intersection of sets of insignificant SNPs from clients.
    """
    if not snp_sets:
        return set()
    merged = snp_sets[0]
    for s in snp_sets[1:]:
        merged = merged.intersection(s)
    return merged