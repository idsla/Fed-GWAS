# server/aggregator_king.py

import os
import subprocess
import tarfile
import uuid
import numpy as np


def run_server_king(server_strategy, parameters_list, output_dir=None):
    """
    Server-based KING aggregator with configurable output directory.
    """
    if output_dir is None:
        output_dir = getattr(server_strategy, 'output_dir', '.')
    session_id = uuid.uuid4().hex
    session_path = os.path.join(output_dir, session_id)
    os.makedirs(session_path, exist_ok=True)
    merged_prefix = f"{session_path}/merged_king"
    bed_files = []

    # Unpack each client's tar file and collect bed file prefixes
    for param in parameters_list:
        chunk_bytes = param.numpy.tobytes()
        tar_path = os.path.join(session_path, f"{uuid.uuid4().hex}.tar")
        with open(tar_path, "wb") as f:
            f.write(chunk_bytes)
        with tarfile.open(tar_path, "r") as tf:
            tf.extractall(session_path)
        # Assume each tar contains files named exactly: chunk.bed, chunk.bim, chunk.fam
        chunk_uuid = uuid.uuid4().hex
        bed_file = os.path.join(session_path, f"chunk_{chunk_uuid}")
        os.rename(os.path.join(session_path, "chunk.bed"), f"{bed_file}.bed")
        os.rename(os.path.join(session_path, "chunk.bim"), f"{bed_file}.bim")
        os.rename(os.path.join(session_path, "chunk.fam"), f"{bed_file}.fam")
        bed_files.append(bed_file)
        os.remove(tar_path)

    if not bed_files:
        print("[Server KING] No bed files received.")
        return np.array([], dtype=np.uint8)

    # Merge bed files:
    # Use the first bed set as the starting point
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
            print(f"[Server KING] Merge failed: {e}")

    # Compute heterozygosity using PLINK --het on the merged data
    het_cmd = [
        "plink",
        "--bfile", merged_prefix,
        "--het",
        "--out", merged_prefix
    ]
    try:
        subprocess.run(het_cmd, check=True)
        print("[Server KING] Heterozygosity computation complete.")
    except subprocess.CalledProcessError as e:
        print(f"[Server KING] PLINK --het failed: {e}")

    # Parse the .het file to get n1 for each sample
    het_file = f"{merged_prefix}.het"
    n1_dict = {}  # key: sample ID, value: n1 (number of heterozygous SNPs)
    if os.path.exists(het_file):
        with open(het_file, "r") as f:
            header = f.readline()  # skip header
            for line in f:
                parts = line.strip().split()
                if len(parts) < 6:
                    continue
                # Expected columns: FID, IID, O(HOM), E(HOM), N(NM), F
                fid = parts[0]
                iid = parts[1]
                try:
                    n_nm = float(parts[4])
                    o_hom = float(parts[2])
                    n1 = n_nm - o_hom  # number of heterozygous SNPs = N(NM) - O(HOM)
                    n1_dict[fid] = n1  # Assuming FID and IID are identical in our anonymized data
                except ValueError:
                    continue
        # Optionally remove the .het file after parsing
        os.remove(het_file)
    else:
        print("[Server KING] No heterozygosity file found; using placeholder for n1.")

    # Run PLINK --king robust
    king_cmd = [
        "plink",
        "--bfile", merged_prefix,
        "--king", "robust",
        "--out", f"{session_path}/king_results"
    ]
    try:
        subprocess.run(king_cmd, check=True)
        print("[Server KING] KING analysis complete.")
    except subprocess.CalledProcessError as e:
        print(f"[Server KING] PLINK KING failed: {e}")

    # Parse the king_results file (.kin0)
    kin0_file = f"{session_path}/king_results.kin0"
    result_str = ""
    if os.path.exists(kin0_file):
        with open(kin0_file, "r") as f:
            header = f.readline()  # skip header, if any
            for line in f:
                parts = line.strip().split()
                if len(parts) < 5:
                    continue
                # Assume:
                # Column 1: FID1 (anon) ; Column 3: FID2 (anon); Column 5: KING coefficient
                sampleA = parts[0]
                sampleB = parts[2]
                partial_phi = parts[4]
                # Lookup n1 for sampleA using our heterozygosity dictionary
                n1_star = n1_dict.get(sampleA, 0)  # 0 if not found
                result_str += f"{sampleA} {sampleB} {partial_phi} {n1_star}\n"
    else:
        print("[Server KING] No KING result file found.")

    return np.frombuffer(result_str.encode("utf-8"), dtype=np.uint8)