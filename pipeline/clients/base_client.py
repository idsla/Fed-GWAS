# client/base_client.py

import flwr as fl
import numpy as np
import random
import subprocess
import os
import tarfile
import hashlib
import shutil
import platform
import logging
from pathlib import Path

# Cache for PLINK binary path
_PLINK_BINARY = None


def find_plink_binary():
    """Find PLINK binary path"""
    global _PLINK_BINARY
    if _PLINK_BINARY is not None:
        return _PLINK_BINARY

    # Find project root - prioritize current working directory (experiment runner runs from project root)
    project_root = None

    # First: try current working directory (most reliable when running via experiment runner)
    cwd = Path.cwd()
    if (
        (cwd / "plink").exists()
        or (cwd / "experiments").exists()
        or (cwd / "pyproject.toml").exists()
    ):
        project_root = cwd

    # Second: try to find project root by walking up from file location
    if project_root is None:
        current_path = Path(__file__).resolve()
        for parent in list(current_path.parents):
            if (
                (parent / "plink").exists()
                or (parent / "experiments").exists()
                or (parent / "pyproject.toml").exists()
            ):
                project_root = parent
                break

    # Last fallback: try relative to file location (for development)
    if project_root is None:
        project_root = Path(__file__).parent.parent.parent

    # Determine OS-specific PLINK path
    system = platform.system().lower()
    if system == "darwin":  # macOS
        plink_dir = project_root / "plink" / "plink_mac"
        plink_name = "plink"
    elif system == "linux":
        plink_dir = project_root / "plink" / "plink_linux"
        plink_name = "plink"
    elif system == "windows":
        plink_dir = project_root / "plink" / "plink_win"
        plink_name = "plink.exe"
    else:
        plink_dir = project_root / "plink" / "plink_mac"  # Default to macOS
        plink_name = "plink"

    # Check OS-specific location first
    plink_path = plink_dir / plink_name
    if plink_path.exists() and os.access(plink_path, os.X_OK):
        _PLINK_BINARY = str(plink_path)
        return _PLINK_BINARY

    # Check project bin directory
    bin_plink = project_root / "bin" / plink_name
    if bin_plink.exists() and os.access(bin_plink, os.X_OK):
        _PLINK_BINARY = str(bin_plink)
        return _PLINK_BINARY

    # Try in PATH
    if shutil.which("plink"):
        _PLINK_BINARY = "plink"
        return _PLINK_BINARY
    if shutil.which("plink2"):
        _PLINK_BINARY = "plink2"
        return _PLINK_BINARY

    # Last resort: return the expected path (will fail with clear error)
    _PLINK_BINARY = str(plink_path)
    return _PLINK_BINARY


def run_plink_command(cmd, plink_binary=None):
    """Run PLINK command, automatically replacing 'plink' with a resolved binary path.

    Args:
        cmd: PLINK command arguments. Any literal ``"plink"`` entry is replaced.
        plink_binary: Optional explicit PLINK executable path. When omitted, the
            project-level resolver is used for backwards compatibility.
    """
    # Replace 'plink' in command with actual path
    plink_binary = plink_binary or find_plink_binary()
    cmd = [plink_binary if arg == "plink" else arg for arg in cmd]

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"PLINK command failed: {e}")


def anonymize_snp_id(
    old_snp: str, global_seed: int, chunk_index: int = 0, iteration_id: int = 0
) -> str:
    """
    Deterministically map old_snp -> new_snp using global_seed plus chunk/iteration
    to avoid cross-chunk/iteration linkability.
    """
    m = hashlib.sha256()
    # Combine the seed, chunk_index, iteration_id, and the old SNP ID
    seed_str = f"{global_seed}-{chunk_index}-{iteration_id}-{old_snp}"
    m.update(seed_str.encode("utf-8"))
    # Take the first 4 bytes from the hash => up to 2^32 range
    short_int = int.from_bytes(m.digest()[:4], "big")
    new_snp = f"rs{short_int}"
    return new_snp


def anonymize_bed_chunk(
    prefix: str,
    new_prefix: str,
    sample_offset: int,
    global_seed: int,
    simulation_mode: bool = False,
    anonymize_snps: bool = True,
    chunk_index: int = 0,
    iteration_id: int = 0,
):
    """
    1) For each sample in .fam, produce new IDs based on sample_offset, chunk_index, and iteration_id.
       This ensures the same sample gets different anonymized IDs in different chunks/iterations (better privacy).
    2) For each SNP in .bim, produce new IDs using the global_seed-based hash.
    3) Return two dictionaries for chunk-based ID maps:
       sample_map and snp_map, each used for inverting new->old IDs.

    The final .bed/.bim/.fam are written to new_prefix.*,
    with the old->new transformation in place.

    Security Note: chunk_index and iteration_id are included in anonymization to prevent
    tracking the same sample across different chunks or iterations.
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
            # Use deterministic hashing for non-numeric IDs to avoid process-dependent hashes
            def _stable_id(val: str) -> int:
                try:
                    return int(val)
                except ValueError:
                    import hashlib

                    return int(hashlib.sha256(val.encode()).hexdigest()[:9], 16)

            fid_int = _stable_id(old_fid)
            iid_int = _stable_id(old_iid)

            # Include chunk_index and iteration_id in anonymization to prevent cross-chunk/iteration tracking
            # Use large multipliers to ensure different chunks/iterations produce non-overlapping ID spaces
            chunk_salt = chunk_index * (10**9)  # 1 billion per chunk
            iteration_salt = iteration_id * (10**12)  # 1 trillion per iteration

            # Create a hash that includes all components for better security
            # This ensures deterministic mapping within a chunk/iteration, but different across chunks/iterations
            def _anonymize_sample_id(
                orig_id_int: int, chunk_idx: int, iter_id: int
            ) -> str:
                """Anonymize sample ID with chunk and iteration-specific variation."""
                # Combine all components into a hash for better security
                combined = (
                    f"{sample_offset}_{orig_id_int}_{chunk_idx}_{iter_id}".encode()
                )
                hash_val = int(
                    hashlib.sha256(combined).hexdigest()[:12], 16
                )  # 12 hex digits = 48 bits
                # Use a large base offset to ensure IDs don't collide across clients
                return str(
                    sample_offset + hash_val % (10**12)
                )  # Keep within reasonable range

            new_fid = _anonymize_sample_id(fid_int, chunk_index, iteration_id)
            new_iid = _anonymize_sample_id(iid_int, chunk_index, iteration_id)

            # Store map so we can invert later (server sees new, we get old)
            sample_map[new_fid] = old_fid
            sample_map[new_iid] = old_iid

            # Keep rest columns the same, but optionally adjust phenotype for simulation
            if simulation_mode and len(parts) >= 6:
                # PLINK: 1=control, 2=case, -9/0=missing
                if parts[5] in ("-9", "0"):
                    parts[5] = "2"
            new_line = [new_fid, new_iid] + parts[2:]
            fout.write("\t".join(new_line) + "\n")

    # 2) Anonymize .bim (optionally keep SNP IDs for downstream alignment/comparison)
    with open(bim_file, "r") as fin, open(new_prefix + ".bim", "w") as fout:
        for line in fin:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            old_snp = parts[1]
            if anonymize_snps:
                new_snp = anonymize_snp_id(
                    old_snp, global_seed, chunk_index, iteration_id
                )
            else:
                new_snp = old_snp
            snp_map[new_snp] = old_snp  # new->old
            parts[1] = new_snp
            fout.write("\t".join(parts) + "\n")

    # 3) Rename .bed
    # os.rename(bed_file, new_prefix + ".bed")

    shutil.copy(bed_file, new_prefix + ".bed")

    # Persist sample map to a TSV for downstream de-anonymization (e.g., KING)
    # Include iteration_id in filename so maps can be matched to the correct round.
    try:
        sample_map_path = (
            Path(new_prefix).parent
            / f"{Path(new_prefix).name}_iter{iteration_id}_sample_map.tsv"
        )
        with open(sample_map_path, "w") as f:
            f.write("anon_sample\torig_sample\n")
            for anon_id, orig_id in sample_map.items():
                f.write(f"{anon_id}\t{orig_id}\n")
    except Exception:
        pass

    return sample_map, snp_map


def _write_snp_map_file(map_path: str, snp_map: dict):
    """Persist anonymized->original SNP map to a TSV for downstream de-anonymization."""
    try:
        with open(map_path, "w") as f:
            f.write("anon_snp\torig_snp\n")
            for anon, orig in snp_map.items():
                f.write(f"{anon}\t{orig}\n")
    except Exception as e:
        logging.getLogger().warning(
            f"[BaseClient] Failed to write snp_map to {map_path}: {e}"
        )


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

    Attributes:
        plink_prefix (str): Path to the PLINK dataset prefix (e.g., "data/client_data")
        client_id (str): Unique identifier for the client
        partition_by (str): How to partition the data ("samples", "snps", or "both")
    """

    def __init__(
        self,
        plink_prefix: str,
        client_id: str,
        partition_by="samples",
        simulation_mode: bool = False,
        phenotype_fix_missing_to_case: bool = False,
        sample_offset: int = None,
    ):

        self.plink_prefix = plink_prefix
        self.client_id = client_id
        self.partition_by = partition_by
        self.simulation_mode = simulation_mode
        # When enabled, convert -9/0 (missing) to 2 (case) for toy/simulated phenotypes.
        # This is intentionally opt-in outside simulation mode.
        self.phenotype_fix_missing_to_case = phenotype_fix_missing_to_case

        # Sample offset for ID anonymization (from config, or random if not provided)
        if sample_offset is not None:
            self.sample_offset = sample_offset
        else:
            # Fallback to random if not provided (for backward compatibility)
            self.sample_offset = random.randint(10**12, 10**13)

        self.local_seed = random.randint(0, 10**8)
        self.global_seed = None

        self.chunk_files = []
        self.current_chunk_idx = 0
        self.iteration_results = []
        self.chunk_data = []  # Initialize chunk_data attribute

        # dictionaries to store chunk-specific ID maps
        # e.g. chunk_sample_map[i][newID] = oldID
        #      chunk_snp_map[i][newSNP]   = oldSNP
        self.chunk_sample_map = {}
        self.chunk_snp_map = {}

    def _fix_phenotype_encoding(self, fam_file: str) -> None:
        """
        Fix phenotype encoding in .fam file: convert -9/0 (missing) to 2 (case).
        NOTE: Applied when simulation_mode is on OR when phenotype_fix_missing_to_case is explicitly enabled.
        """
        if not self._should_fix_phenotype():
            return
        if not os.path.exists(fam_file):
            return
        with open(fam_file, "r") as f:
            lines = f.readlines()
        fixed_lines = []
        for line in lines:
            parts = line.strip().split()
            if len(parts) >= 6 and parts[5] in ("-9", "0"):
                parts[5] = "2"
                fixed_lines.append("\t".join(parts) + "\n")
            else:
                fixed_lines.append(line if line.endswith("\n") else (line + "\n"))
        with open(fam_file, "w") as f:
            f.writelines(fixed_lines)

    def _should_fix_phenotype(self) -> bool:
        return bool(
            getattr(self, "simulation_mode", False)
            or getattr(self, "phenotype_fix_missing_to_case", False)
        )

    def ensure_phenotype_fixed_prefix(self, out_prefix: str) -> str:
        """
        If phenotype fixing is enabled, create a derived PLINK prefix with a rewritten .fam
        (missing phenotype -9/0 -> 2) without modifying the original dataset.
        Returns the prefix to use (either self.plink_prefix or out_prefix).
        """
        if not self._should_fix_phenotype():
            return self.plink_prefix

        src_prefix = self.plink_prefix
        src_bed, src_bim, src_fam = (
            src_prefix + ".bed",
            src_prefix + ".bim",
            src_prefix + ".fam",
        )
        dst_bed, dst_bim, dst_fam = (
            out_prefix + ".bed",
            out_prefix + ".bim",
            out_prefix + ".fam",
        )

        if (
            os.path.exists(dst_bed)
            and os.path.exists(dst_bim)
            and os.path.exists(dst_fam)
        ):
            return out_prefix

        os.makedirs(os.path.dirname(os.path.abspath(dst_fam)), exist_ok=True)
        shutil.copyfile(src_bed, dst_bed)
        shutil.copyfile(src_bim, dst_bim)

        # Rewrite .fam with adjusted phenotype encoding
        if os.path.exists(src_fam):
            with open(src_fam, "r") as fin, open(dst_fam, "w") as fout:
                for line in fin:
                    parts = line.strip().split()
                    if len(parts) >= 6 and parts[5] in ("-9", "0"):
                        parts[5] = "2"
                        fout.write("\t".join(parts) + "\n")
                    else:
                        fout.write(line if line.endswith("\n") else (line + "\n"))
        else:
            # If fam missing, just create empty placeholder to avoid crashes
            open(dst_fam, "a").close()

        return out_prefix

    def get_parameters(self, config):
        return []

    def evaluate(self, parameters, config):
        return 0.0, 1, {}

    def run_make_bed(self, in_prefix, out_prefix, extra_args=None):
        """
        Run PLINK --make-bed with optional extra_args (e.g. --keep, --extract).
        """
        cmd = ["plink", "--bfile", in_prefix, "--make-bed", "--allow-no-sex", "--out", out_prefix]
        if extra_args:
            cmd[1:1] = extra_args  # Insert after 'plink'
        run_plink_command(cmd)

    def anonymize_and_tar(
        self, chunk_prefix: str, chunk_index: int, iteration_id: int = 0
    ) -> str:
        """
        1) Anonymize the BED chunk at chunk_prefix using sample_offset, global_seed, chunk_index, and iteration_id.
        2) Capture the sample_map and snp_map for chunk_index.
        3) Tar the anonymized .bed/.bim/.fam and return the path to that .tar file.
        """
        anon_prefix = f"{chunk_prefix}_anon"
        sample_map, snp_map = anonymize_bed_chunk(
            prefix=chunk_prefix,
            new_prefix=anon_prefix,
            sample_offset=self.sample_offset,
            global_seed=self.global_seed or 0,
            simulation_mode=self._should_fix_phenotype(),
            anonymize_snps=True,
            chunk_index=chunk_index,
            iteration_id=iteration_id,
        )

        # store these maps for later ID inversion
        self.chunk_sample_map[chunk_index] = sample_map
        self.chunk_snp_map[chunk_index] = snp_map
        # Persist snp map for downstream de-anonymization of global LR results
        try:
            map_filename = f"{Path(anon_prefix).name}_snp_map.tsv"
            if iteration_id is not None:
                map_filename = f"{Path(anon_prefix).name}_iter{iteration_id}_snp_map.tsv"
            map_path = str(Path(anon_prefix).parent / map_filename)
            _write_snp_map_file(map_path, snp_map)
        except Exception:
            pass

        # create the tar file
        tar_file = create_tar(anon_prefix)
        return tar_file

    def partition_data(self, config):
        """
        Partition local data into .bed chunks.
        Then anonymize the .bed chunk (renaming sample and SNP IDs).
        Finally, return the chunk data as numpy arrays for transmission.
        """

        partition_by = self.partition_by
        chunk_size = config.get("chunk_size", 100)
        intermediate_dir = getattr(self, "intermediate_dir", "intermediate")
        os.makedirs(intermediate_dir, exist_ok=True)

        if partition_by == "samples":
            # Ensure plink_prefix is an absolute path
            if not os.path.isabs(self.plink_prefix):
                self.plink_prefix = os.path.abspath(self.plink_prefix)

            fam_file = self.plink_prefix + ".fam"
            if not os.path.exists(fam_file):
                raise FileNotFoundError(f"{fam_file} not found.")

            with open(fam_file, "r") as f:
                lines = f.readlines()

            fid_map = {}
            sample_ids = []
            for line in lines:
                parts = line.strip().split()
                if len(parts) >= 2:
                    fid, iid = parts[0], parts[1]
                    sample_ids.append(iid)
                    fid_map[iid] = fid

            # Deterministic chunking with iteration variation:
            # - Same global_seed + same sample_ids + same server_round = same chunks
            # - Different server_round = different chunks (even with same seed and samples)
            # Note: If samples are filtered (e.g., after KING), the sample_ids list changes,
            # so chunks will be different even with the same seed (this is expected behavior)
            server_round = config.get("server_round", 0)
            # Combine global_seed with server_round to get different chunks per iteration
            chunk_seed = (self.global_seed + int(server_round)) % (
                2**31
            )  # Keep within int32 range
            random.seed(chunk_seed)
            random.shuffle(sample_ids)
            chunks = [
                sample_ids[i : i + chunk_size]
                for i in range(0, len(sample_ids), chunk_size)
            ]

            chunk_data = []
            for idx, chunk_sids in enumerate(chunks):
                keep_file = os.path.join(
                    intermediate_dir, f"temp_keep_{self.client_id}_{idx}.txt"
                )
                with open(keep_file, "w") as f:
                    f.write("FID IID\n")
                    for sid in chunk_sids:
                        fid = fid_map.get(sid, "0")
                        f.write(f"{fid} {sid}\n")

                chunk_prefix = os.path.join(
                    intermediate_dir, f"chunk_{self.client_id}_{idx}"
                )
                cmd = [
                    "plink",
                    "--bfile",
                    self.plink_prefix,
                    "--keep",
                    keep_file,
                    "--make-bed",
                    "--allow-no-sex",
                    "--out",
                    chunk_prefix,
                ]
                run_plink_command(cmd)
                os.remove(keep_file)

                # Simulation-only phenotype fix (helps toy LR have cases)
                self._fix_phenotype_encoding(f"{chunk_prefix}.fam")

                # Anonymize the chunk before sending to server
                # Include chunk_index and server_round to prevent cross-chunk/iteration tracking
                anon_prefix = f"{chunk_prefix}_anon"
                sample_map, snp_map = anonymize_bed_chunk(
                    prefix=chunk_prefix,
                    new_prefix=anon_prefix,
                sample_offset=self.sample_offset,
                global_seed=self.global_seed or 0,
                simulation_mode=self._should_fix_phenotype(),
                anonymize_snps=True,
                chunk_index=idx,
                iteration_id=server_round,
                )

                # Store ID maps for later inversion (if needed)
                self.chunk_sample_map[idx] = sample_map
                self.chunk_snp_map[idx] = snp_map
                # Persist snp map for downstream de-anonymization of global LR results
                try:
                    map_filename = f"{Path(anon_prefix).name}_snp_map.tsv"
                    if server_round is not None:
                        map_filename = f"{Path(anon_prefix).name}_iter{server_round}_snp_map.tsv"
                    map_path = str(Path(anon_prefix).parent / map_filename)
                    _write_snp_map_file(map_path, snp_map)
                except Exception:
                    pass
                # Debug: log mapping sizes to help KING filtering mapping
                try:
                    map_size = len(sample_map) if hasattr(sample_map, "__len__") else -1
                    logging.getLogger().info(
                        f"[Client {self.client_id}] chunk_sample_map[{idx}] size={map_size}"
                    )
                except Exception:
                    pass

                # Read anonymized chunk files and convert to numpy arrays
                chunk_data.append(self._read_chunk_as_array(anon_prefix, idx))

                # Clean up chunk files (both original and anonymized)
                os.remove(f"{chunk_prefix}.bed")
                os.remove(f"{chunk_prefix}.bim")
                os.remove(f"{chunk_prefix}.fam")
                os.remove(f"{anon_prefix}.bed")
                os.remove(f"{anon_prefix}.bim")
                os.remove(f"{anon_prefix}.fam")

            self.chunk_data = chunk_data
            return chunk_data

        elif partition_by == "snps":
            # Ensure plink_prefix is an absolute path
            if not os.path.isabs(self.plink_prefix):
                self.plink_prefix = os.path.abspath(self.plink_prefix)

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

            server_round = config.get("server_round", 0)
            random.seed(self.global_seed)
            random.shuffle(snp_ids)
            chunks = [
                snp_ids[i : i + chunk_size] for i in range(0, len(snp_ids), chunk_size)
            ]

            chunk_data = []
            for idx, chunk_snps in enumerate(chunks):
                extract_file = os.path.join(
                    intermediate_dir, f"temp_extract_{self.client_id}_{idx}.txt"
                )
                with open(extract_file, "w") as f:
                    for snp in chunk_snps:
                        f.write(f"{snp}\n")

                chunk_prefix = os.path.join(
                    intermediate_dir, f"chunk_{self.client_id}_{idx}"
                )
                cmd = [
                    "plink",
                    "--bfile",
                    self.plink_prefix,
                    "--extract",
                    extract_file,
                    "--make-bed",
                    "--allow-no-sex",
                    "--out",
                    chunk_prefix,
                ]
                run_plink_command(cmd)
                os.remove(extract_file)

                # Simulation-only phenotype fix (helps toy LR have cases)
                self._fix_phenotype_encoding(f"{chunk_prefix}.fam")

                # Anonymize the chunk before sending to server
                # Include chunk_index and server_round to prevent cross-chunk/iteration tracking
                anon_prefix = f"{chunk_prefix}_anon"
                sample_map, snp_map = anonymize_bed_chunk(
                    prefix=chunk_prefix,
                    new_prefix=anon_prefix,
                    sample_offset=self.sample_offset,
                    global_seed=self.global_seed or 0,
                    simulation_mode=self._should_fix_phenotype(),
                    anonymize_snps=True,
                    chunk_index=idx,
                    iteration_id=server_round,
                )

                # Store ID maps for later inversion (if needed)
                self.chunk_sample_map[idx] = sample_map
                self.chunk_snp_map[idx] = snp_map
                # Persist snp map for downstream de-anonymization of global LR results
                try:
                    map_filename = f"{Path(anon_prefix).name}_snp_map.tsv"
                    if server_round is not None:
                        map_filename = f"{Path(anon_prefix).name}_iter{server_round}_snp_map.tsv"
                    map_path = str(Path(anon_prefix).parent / map_filename)
                    _write_snp_map_file(map_path, snp_map)
                except Exception:
                    pass

                # Read anonymized chunk files and convert to numpy arrays
                chunk_data.append(self._read_chunk_as_array(anon_prefix, idx))

                # Clean up chunk files (both original and anonymized)
                os.remove(f"{chunk_prefix}.bed")
                os.remove(f"{chunk_prefix}.bim")
                os.remove(f"{chunk_prefix}.fam")
                os.remove(f"{anon_prefix}.bed")
                os.remove(f"{anon_prefix}.bim")
                os.remove(f"{anon_prefix}.fam")

            self.chunk_data = chunk_data
            return chunk_data

        elif partition_by == "both":
            sample_chunk_size = config.get("sample_chunk_size", chunk_size)
            snp_chunk_size = config.get("snp_chunk_size", chunk_size)

            # Ensure plink_prefix is an absolute path
            if not os.path.isabs(self.plink_prefix):
                self.plink_prefix = os.path.abspath(self.plink_prefix)

            fam_file = self.plink_prefix + ".fam"
            if not os.path.exists(fam_file):
                raise FileNotFoundError(f"{fam_file} not found.")

            # This is a simplified implementation for "both" partitioning
            # In practice, you might want to implement a more sophisticated approach
            chunk_data = []
            chunk_data.append(self._read_chunk_as_array(self.plink_prefix, 0))
            self.chunk_data = chunk_data
            return chunk_data

    def _read_chunk_as_array(self, chunk_prefix, idx):
        """
        Read .bed, .bim, and .fam files as bytes and concatenate them.
        Returns a numpy array with the concatenated data and metadata for splitting.
        """
        bed_file = f"{chunk_prefix}.bed"
        bim_file = f"{chunk_prefix}.bim"
        fam_file = f"{chunk_prefix}.fam"

        # Read all three files
        with open(bed_file, "rb") as f:
            bed_data = f.read()
        with open(bim_file, "rb") as f:
            bim_data = f.read()
        with open(fam_file, "rb") as f:
            fam_data = f.read()

        # Store file sizes for proper reconstruction on server side
        bed_size = len(bed_data)
        bim_size = len(bim_data)
        fam_size = len(fam_data)

        # Concatenate all data
        combined_data = bed_data + bim_data + fam_data

        # Create metadata array with sizes
        metadata = np.array([bed_size, bim_size, fam_size], dtype=np.uint32)

        # Combine metadata and data
        result = np.concatenate(
            [metadata, np.frombuffer(combined_data, dtype=np.uint8)]
        )

        return result
