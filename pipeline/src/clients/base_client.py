# client/base_client.py

import flwr as fl
import numpy as np
import random
import subprocess
import os
import logging
import tarfile
import hashlib


logging.basicConfig(
    filename="iteration_log.txt",
    level=logging.INFO,
    format="%(asctime)s %(message)s"
)

def run_plink_command(cmd):
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"PLINK command failed: {e}")


def anonymize_snp_id(old_snp: str, global_seed: int) -> str:
    """
    Deterministically map old_snp -> new_snp using global_seed so that
    all clients produce the SAME mapping for each SNP ID.
    """
    m = hashlib.sha256()
    # Combine the seed and the old SNP ID
    seed_str = f"{global_seed}-{old_snp}"
    m.update(seed_str.encode("utf-8"))
    # Let's take the first 4 bytes from the hash => up to 2^32 range
    short_int = int.from_bytes(m.digest()[:4], "big")
    # Example: newSnp = f"rs{short_int:09d}"
    new_snp = f"rs{short_int}"
    return new_snp


def anonymize_bed_chunk(
        prefix: str,
        new_prefix: str,
        sample_offset: int,
        global_seed: int
):
    """
    1) For each sample in .fam, produce new IDs based on sample_offset.
    2) For each SNP in .bim, produce new IDs using the global_seed-based hash.
    3) Return two dictionaries for chunk-based ID maps:
       sample_map and snp_map, each used for inverting new->old IDs.

    The final .bed/.bim/.fam are written to new_prefix.*,
    with the old->new transformation in place.
    """

    fam_file = prefix + ".fam"
    bim_file = prefix + ".bim"
    bed_file = prefix + ".bed"

    sample_map = {}  # newSampleID -> oldSampleID
    snp_map = {}  # newSnpID -> oldSnpID

    # 1) Anonymize .fam
    with open(fam_file, "r") as fin, open(new_prefix + ".fam", "w") as fout:
        for line in fin:
            parts = line.strip().split()
            if len(parts) < 2:
                continue

            old_fid, old_iid = parts[0], parts[1]

            # Convert old IDs to integers if possible; fallback to hashing
            try:
                fid_int = int(old_fid)
            except ValueError:
                fid_int = abs(hash(old_fid)) % (10 ** 9)
            try:
                iid_int = int(old_iid)
            except ValueError:
                iid_int = abs(hash(old_iid)) % (10 ** 9)

            new_fid = str(sample_offset + fid_int)
            new_iid = str(sample_offset + iid_int)

            # Store map so we can invert later (server sees new, we get old)
            sample_map[new_fid] = old_fid
            sample_map[new_iid] = old_iid

            # Keep rest columns the same
            new_line = [new_fid, new_iid] + parts[2:]
            fout.write("\t".join(new_line) + "\n")

    # 2) Anonymize .bim
    with open(bim_file, "r") as fin, open(new_prefix + ".bim", "w") as fout:
        for line in fin:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            old_snp = parts[1]
            new_snp = anonymize_snp_id(old_snp, global_seed)
            snp_map[new_snp] = old_snp  # new->old
            parts[1] = new_snp
            fout.write("\t".join(parts) + "\n")

    # 3) Rename .bed
    os.rename(bed_file, new_prefix + ".bed")

    return sample_map, snp_map


def create_tar(anon_prefix: str) -> str:
    """
    Archive anonymized .bed/.bim/.fam into a single .tar file for transport.
    """
    tar_file = f"{anon_prefix}.tar"
    with tarfile.open(tar_file, "w") as tf:
        tf.add(f"{anon_prefix}.bed", arcname="chunk.bed")
        tf.add(f"{anon_prefix}.bim", arcname="chunk.bim")
        tf.add(f"{anon_prefix}.fam", arcname="chunk.fam")
    # Remove the anonymized bed/bim/fam
    os.remove(f"{anon_prefix}.bed")
    os.remove(f"{anon_prefix}.bim")
    os.remove(f"{anon_prefix}.fam")
    return tar_file


class BaseGWASClient(fl.client.NumPyClient):
    """
    Base client containing shared attributes & methods.
    Specific stages (QC, iterative LR, etc.) can be implemented in separate modules
    and imported into the final client class.
    """
    def __init__(self, plink_prefix, client_id="client", partition_by="samples"):
        self.plink_prefix = plink_prefix
        self.client_id = client_id
        self.partition_by = partition_by

        # random offset for sample ID anonymization
        self.sample_offset = random.randint(10**12, 10**13)

        self.local_seed = random.randint(0, 10**8)
        self.global_seed = None

        self.chunk_files = []
        self.current_chunk_idx = 0
        self.iteration_results = []

        # dictionaries to store chunk-specific ID maps
        # e.g. chunk_sample_map[i][newID] = oldID
        #      chunk_snp_map[i][newSNP]   = oldSNP
        self.chunk_sample_map = {}
        self.chunk_snp_map = {}

    def get_parameters(self, config):
        return []

    def evaluate(self, parameters, config):
        return 0.0, 1, {}

    def run_make_bed(self, in_prefix, out_prefix, extra_args=None):
        """
        Run PLINK --make-bed with optional extra_args (e.g. --keep, --extract).
        """
        cmd = [
            "plink",
            "--bfile", in_prefix,
            "--make-bed",
            "--out", out_prefix
        ]
        if extra_args:
            cmd[1:1] = extra_args  # Insert after 'plink'
        run_plink_command(cmd)

    def anonymize_and_tar(self, chunk_prefix: str, chunk_index: int) -> str:
        """
        1) Anonymize the BED chunk at chunk_prefix using sample_offset and global_seed.
        2) Capture the sample_map and snp_map for chunk_index.
        3) Tar the anonymized .bed/.bim/.fam and return the path to that .tar file.
        """
        anon_prefix = f"{chunk_prefix}_anon"
        sample_map, snp_map = anonymize_bed_chunk(
            prefix=chunk_prefix,
            new_prefix=anon_prefix,
            sample_offset=self.sample_offset,
            global_seed=self.global_seed or 0
        )

        # store these maps for later ID inversion
        self.chunk_sample_map[chunk_index] = sample_map
        self.chunk_snp_map[chunk_index] = snp_map

        # create the tar file
        tar_file = create_tar(anon_prefix)
        return tar_file

    def partition_data(self, config):
        """
        Partition local data into .bed chunks.
        Then anonymize the .bed chunk (renaming sample and SNP IDs).
        Finally, archive the .bed, .bim, .fam into a single .tar for each chunk.
        """

        partition_by = self.partition_by
        chunk_size = config.get("chunk_size", 100)

        if partition_by == "samples":
            # Partition by samples using the .fam file
            fam_file = self.plink_prefix + ".fam"
            if not os.path.exists(fam_file):
                raise FileNotFoundError(f"{fam_file} not found.")

            with open(fam_file, "r") as f:
                lines = f.readlines()

            sample_ids = []
            for line in lines:
                parts = line.strip().split()
                if len(parts) >= 2:
                    sample_ids.append(parts[1])

            random.seed(self.global_seed)
            random.shuffle(sample_ids)
            chunks = [sample_ids[i: i + chunk_size] for i in range(0, len(sample_ids), chunk_size)]

            chunk_files = []
            for idx, chunk_sids in enumerate(chunks):
                keep_file = f"temp_keep_{self.client_id}_{idx}.txt"
                with open(keep_file, "w") as f:
                    for sid in chunk_sids:
                        f.write(f"{sid}\n")

                chunk_prefix = f"chunk_{self.client_id}_{idx}"
                cmd = [
                    "plink",
                    "--bfile", self.plink_prefix,
                    "--keep", keep_file,
                    "--make-bed",
                    "--out", chunk_prefix
                ]
                run_plink_command(cmd)
                os.remove(keep_file)

                # anonymize and tar using our base_client method
                tar_file = self.anonymize_and_tar(chunk_prefix, idx)
                chunk_files.append(tar_file)

            self.chunk_files = chunk_files

        elif partition_by == "snps":
            # Partition by SNPs using the .bim file
            bim_file = self.plink_prefix + ".bim"
            if not os.path.exists(bim_file):
                raise FileNotFoundError(f"{bim_file} not found.")

            with open(bim_file, "r") as f:
                lines = f.readlines()

            snp_ids = []
            for line in lines:
                parts = line.strip().split()
                if len(parts) >= 2:
                    snp_ids.append(parts[1])

            random.seed(self.global_seed)
            random.shuffle(snp_ids)
            chunks = [snp_ids[i: i + chunk_size] for i in range(0, len(snp_ids), chunk_size)]

            chunk_files = []
            for idx, chunk_snps in enumerate(chunks):
                extract_file = f"temp_extract_{self.client_id}_{idx}.txt"
                with open(extract_file, "w") as f:
                    for snp in chunk_snps:
                        f.write(f"{snp}\n")

                chunk_prefix = f"chunk_{self.client_id}_{idx}"
                cmd = [
                    "plink",
                    "--bfile", self.plink_prefix,
                    "--extract", extract_file,
                    "--make-bed",
                    "--out", chunk_prefix
                ]
                run_plink_command(cmd)
                os.remove(extract_file)

                tar_file = self.anonymize_and_tar(chunk_prefix, idx)
                chunk_files.append(tar_file)

            self.chunk_files = chunk_files

        elif partition_by == "both":
            # Partition both samples and SNPs.
            sample_chunk_size = config.get("sample_chunk_size", chunk_size)
            snp_chunk_size = config.get("snp_chunk_size", chunk_size)

            fam_file = self.plink_prefix + ".fam"
            if not os.path.exists(fam_file):
                raise FileNotFoundError(f"{fam_file} not found.")

            with open(fam_file, "r") as f:
                fam_lines = f.readlines()

            sample_ids = []
            for line in fam_lines:
                parts = line.strip().split()
                if len(parts) >= 2:
                    sample_ids.append(parts[1])
            random.seed(self.global_seed)
            random.shuffle(sample_ids)
            sample_chunks = [
                sample_ids[i: i + sample_chunk_size]
                for i in range(0, len(sample_ids), sample_chunk_size)
            ]

            bim_file = self.plink_prefix + ".bim"
            if not os.path.exists(bim_file):
                raise FileNotFoundError(f"{bim_file} not found.")

            with open(bim_file, "r") as f:
                bim_lines = f.readlines()

            snp_ids = []
            for line in bim_lines:
                parts = line.strip().split()
                if len(parts) >= 2:
                    snp_ids.append(parts[1])
            random.seed(self.global_seed)
            random.shuffle(snp_ids)
            snp_chunks = [
                snp_ids[i: i + snp_chunk_size]
                for i in range(0, len(snp_ids), snp_chunk_size)
            ]

            # Pair them using the minimum number of chunks
            num_pairs = min(len(sample_chunks), len(snp_chunks))
            chunk_files = []
            for idx in range(num_pairs):
                keep_file = f"temp_keep_{self.client_id}_{idx}.txt"
                with open(keep_file, "w") as f:
                    for sid in sample_chunks[idx]:
                        f.write(f"{sid}\n")

                extract_file = f"temp_extract_{self.client_id}_{idx}.txt"
                with open(extract_file, "w") as f:
                    for snp in snp_chunks[idx]:
                        f.write(f"{snp}\n")

                chunk_prefix = f"chunk_{self.client_id}_{idx}"
                cmd = [
                    "plink",
                    "--bfile", self.plink_prefix,
                    "--keep", keep_file,
                    "--extract", extract_file,
                    "--make-bed",
                    "--out", chunk_prefix
                ]
                run_plink_command(cmd)
                os.remove(keep_file)
                os.remove(extract_file)

                tar_file = self.anonymize_and_tar(chunk_prefix, idx)
                chunk_files.append(tar_file)

            self.chunk_files = chunk_files