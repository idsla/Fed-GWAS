# server/aggregator_king.py

import os
import subprocess
import uuid
import platform
import shutil
import logging
from pathlib import Path
import numpy as np
from flwr.common import parameters_to_ndarrays

logger = logging.getLogger(__name__)

# Cache for PLINK binary paths
_PLINK_BINARY = None
_PLINK2_BINARY = None

def find_plink_binary():
    """Find PLINK 1.9 binary path"""
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

def find_plink2_binary():
    """Find PLINK 2.0 binary path"""
    global _PLINK2_BINARY
    if _PLINK2_BINARY is not None:
        return _PLINK2_BINARY
    
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
    
    # Determine OS-specific PLINK2 path
    system = platform.system().lower()
    if system == "darwin":
        plink2_dir = project_root / "plink" / "plink_mac"
        plink2_name = "plink2"
    elif system == "linux":
        plink2_dir = project_root / "plink" / "plink_linux"
        plink2_name = "plink2"
    elif system == "windows":
        plink2_dir = project_root / "plink" / "plink_win"
        plink2_name = "plink2.exe"
    else:
        plink2_dir = project_root / "plink" / "plink_mac"
        plink2_name = "plink2"
    
    plink2_path = plink2_dir / plink2_name
    if plink2_path.exists() and os.access(plink2_path, os.X_OK):
        _PLINK2_BINARY = str(plink2_path)
        return _PLINK2_BINARY
    
    if shutil.which("plink2"):
        _PLINK2_BINARY = "plink2"
        return _PLINK2_BINARY
    
    _PLINK2_BINARY = str(plink2_path)
    return _PLINK2_BINARY


def run_server_king(server_strategy, results=None, parameters_list=None, output_dir=None):
    """
    Run KING analysis on the server side.
    
    Args:
        server_strategy: The server strategy object
        results: List of (ClientProxy, FitRes) tuples (preferred, for deterministic sorting)
        parameters_list: List of parameters from clients (fallback, for backward compatibility)
        output_dir: Directory to store output files
    """
    if output_dir is None:
        output_dir = getattr(server_strategy, 'output_dir', '.')
    session_id = uuid.uuid4().hex
    session_path = os.path.join(output_dir, session_id)
    os.makedirs(session_path, exist_ok=True)
    merged_prefix = f"{session_path}/merged_king"
    
    # Prepare list of (client_id, chunk_index, parameters) tuples for deterministic sorting
    chunk_items = []
    
    if results is not None:
        # Extract client_id and chunk_index from results for deterministic sorting
        for client_proxy, fit_res in results:
            if not fit_res.parameters:
                continue
            
            # Extract client_id
            client_id_str = getattr(client_proxy, "cid", str(client_proxy))
            try:
                client_id = int(client_id_str) if client_id_str.isdigit() else hash(client_id_str) % 10000
            except (ValueError, TypeError, AttributeError):
                # Fallback if client_id_str conversion fails
                client_id = hash(str(client_proxy)) % 10000
            
            # Extract chunk_index from metrics
            chunk_index = fit_res.metrics.get("king_chunk_index", 0) if hasattr(fit_res, "metrics") else 0
            
            chunk_items.append((client_id, chunk_index, fit_res.parameters))
        
        # Sort deterministically by (client_id, chunk_index)
        chunk_items.sort(key=lambda x: (x[0], x[1]))
        logger.debug(f"[Server KING] Sorted {len(chunk_items)} chunks deterministically by (client_id, chunk_index)")
        parameters_list = [item[2] for item in chunk_items]
    elif parameters_list is not None:
        # Fallback: use parameters_list directly (non-deterministic order)
        if not isinstance(parameters_list, list):
            parameters_list = [parameters_list]
    else:
        logger.warning("[Server KING] No parameters or results provided")
        return np.array([], dtype=np.uint8)
    
    bed_files = []
    
    # Process each client's chunk data (now in deterministic order)
    for param in parameters_list:
        logger.debug(f"[Server KING] Processing parameter: type={type(param)}")

        nds = parameters_to_ndarrays(param)
        if not nds:
            logger.debug("[Server KING] No ndarrays extracted from Parameters, skipping")
            continue
        chunk_array = nds[0]
        logger.debug(f"[Server KING] Extracted numpy array from Parameters: type={type(chunk_array)}, size={len(chunk_array)}")

        if len(chunk_array) == 0:
            logger.debug("[Server KING] Empty chunk array, skipping")
            continue

        # Recover bytes:
        # - Some clients send uint32 arrays where each element after the 3-size metadata is a single byte (0..255).
        # - Other clients might send uint8 arrays directly.
        if len(chunk_array) < 3:
            logger.debug("[Server KING] Chunk too short for metadata, skipping")
            continue

        bed_size = int(chunk_array[0])
        bim_size = int(chunk_array[1])
        fam_size = int(chunk_array[2])
        total_bytes = bed_size + bim_size + fam_size
        logger.debug(f"[Server KING] Metadata extracted - bed: {bed_size}, bim: {bim_size}, fam: {fam_size}")

        if chunk_array.dtype == np.uint32:
            logger.debug("[Server KING] chunk_array dtype is uint32; interpreting elements [3:] as byte values")
            if len(chunk_array) < 3 + total_bytes:
                logger.debug(
                    f"[Server KING] Data too short: have {len(chunk_array)-3} bytes, need {total_bytes} bytes"
                )
                continue
            data_u8 = chunk_array[3 : 3 + total_bytes].astype(np.uint8, copy=False)
            chunk_bytes = data_u8.tobytes()
        else:
            # Assume uint8-ish array which already represents raw bytes after metadata
            data_u8 = chunk_array[3 : 3 + total_bytes].astype(np.uint8, copy=False)
            chunk_bytes = data_u8.tobytes()

        # Split the data using the stored sizes
        bed_data = chunk_bytes[:bed_size]
        bim_data = chunk_bytes[bed_size : bed_size + bim_size]
        fam_data = chunk_bytes[bed_size + bim_size : bed_size + bim_size + fam_size]
        
        # Verify we have all the data
        if len(bed_data) != bed_size or len(bim_data) != bim_size or len(fam_data) != fam_size:
            logger.debug(f"[Server KING] Data size mismatch: expected bed={bed_size}, bim={bim_size}, fam={fam_size}, "
                  f"got bed={len(bed_data)}, bim={len(bim_data)}, fam={len(fam_data)}")
            continue
                
        logger.debug(f"[Server KING] Data split successfully, creating temporary files")
        
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
        logger.debug(f"[Server KING] Created temporary files for chunk {chunk_uuid}")

    if not bed_files:
        logger.warning("[Server KING] No bed files received.")
        return np.array([], dtype=np.uint8)

    # Merge bed files:
    # Use the first bed set as the starting point
    first_bed = bed_files[0]
    os.rename(f"{first_bed}.bed", f"{merged_prefix}.bed")
    os.rename(f"{first_bed}.bim", f"{merged_prefix}.bim")
    os.rename(f"{first_bed}.fam", f"{merged_prefix}.fam")
    bed_files = bed_files[1:]

    plink_binary = find_plink_binary()
    for bf in bed_files:
        merge_cmd = [
            plink_binary,
            "--bfile", merged_prefix,
            "--bmerge", bf,
            "--make-bed",
            "--allow-no-sex",
            "--out", f"{merged_prefix}_tmp"
        ]
        try:
            subprocess.run(merge_cmd, check=True, capture_output=True, text=True)
            os.rename(f"{merged_prefix}_tmp.bed", f"{merged_prefix}.bed")
            os.rename(f"{merged_prefix}_tmp.bim", f"{merged_prefix}.bim")
            os.rename(f"{merged_prefix}_tmp.fam", f"{merged_prefix}.fam")
            if os.path.exists(f"{merged_prefix}_tmp.log"):
                os.remove(f"{merged_prefix}_tmp.log")
        except subprocess.CalledProcessError as e:
            logger.warning(f"[Server KING] Merge failed: {e}")

    # Compute heterozygosity using PLINK --het on the merged data
    het_cmd = [
        plink_binary,
        "--bfile", merged_prefix,
        "--het",
        "--out", merged_prefix
    ]
    try:
        subprocess.run(het_cmd, check=True, capture_output=True, text=True)
        logger.info("[Server KING] Heterozygosity computation complete.")
    except subprocess.CalledProcessError as e:
        logger.warning(f"[Server KING] PLINK --het failed: {e}")

    # Parse the .het file to get n1 for each sample (kept for backward compatibility / optional use)
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
        logger.debug("[Server KING] No heterozygosity file found; n1 will be unavailable.")

    # Run PLINK --king robust
    # PLINK 2.0 is required for KING table output
    plink2_binary = find_plink2_binary()
    king_cmd = [plink2_binary, "--bfile", merged_prefix, "--make-king-table", "--out", f"{session_path}/king_results"]
    
    try:
        # Run PLINK2 and log output directly to server logger
        logger.info(f"[Server KING] Running PLINK2 KING command: {' '.join(king_cmd)}")
        
        result = subprocess.run(king_cmd, check=True, capture_output=True, text=True)
        
        # Log STDERR (usually contains important information)
        if result.stderr:
            logger.info("[Server KING] PLINK2 STDERR:")
            for line in result.stderr.strip().split('\n'):
                if line.strip():
                    logger.info(f"[Server KING]   {line}")
        
        # Log STDOUT with better filtering - include progress and summary information
        if result.stdout:
            lines = result.stdout.split('\n')
            # Include lines with: errors, warnings, progress indicators, summaries, completion status
            important_keywords = ['error', 'warning', 'note', 'complete', 'pass', 'variants', 'samples', 
                                 'people', 'pairs', 'writing', 'reading', 'loaded', 'remaining', 
                                 'total', 'done', 'finished', 'summary', 'rate']
            important_lines = []
            for line in lines:
                line_lower = line.lower().strip()
                # Include non-empty lines that contain important keywords or are progress indicators
                if line_lower and (
                    any(keyword in line_lower for keyword in important_keywords) or
                    line_lower.startswith('%') or  # Progress percentages
                    any(char.isdigit() for char in line_lower[:20])  # Lines starting with numbers (statistics)
                ):
                    important_lines.append(line)
            
            if important_lines:
                logger.info("[Server KING] PLINK2 STDOUT (important lines):")
                for line in important_lines:
                    if line.strip():
                        logger.info(f"[Server KING]   {line}")
            else:
                # If no important lines found, log first 50 lines to ensure we capture something
                logger.info("[Server KING] PLINK2 STDOUT (first 50 lines):")
                for line in lines[:50]:
                    if line.strip():
                        logger.info(f"[Server KING]   {line}")
                if len(lines) > 50:
                    logger.info(f"[Server KING]   ... ({len(lines) - 50} more lines omitted)")
        
        logger.info("[Server KING] KING analysis complete.")
    except subprocess.CalledProcessError as e:
        logger.error(f"[Server KING] PLINK2 KING failed with exit code {e.returncode}")
        if e.stderr:
            logger.error("[Server KING] STDERR:")
            stderr_text = e.stderr.decode() if isinstance(e.stderr, bytes) else e.stderr
            for line in stderr_text.strip().split('\n'):
                if line.strip():
                    logger.error(f"[Server KING]   {line}")
        if e.stdout:
            logger.error("[Server KING] STDOUT:")
            stdout_text = e.stdout.decode() if isinstance(e.stdout, bytes) else e.stdout
            for line in stdout_text.strip().split('\n'):
                if line.strip():
                    logger.error(f"[Server KING]   {line}")

    # Parse the king_results file (.kin0)
    kin0_file = f"{session_path}/king_results.kin0"
    result_str = ""
    debug_str = ""
    have_n1 = bool(n1_dict)
    if os.path.exists(kin0_file):
        with open(kin0_file, "r") as f:
            header = f.readline()  # skip header, if any
            for line in f:
                parts = line.strip().split()
                if len(parts) < 8:
                    continue
                # PLINK2 KING table (.kin0) columns:
                # 0: FID1, 1: IID1, 2: FID2, 3: IID2, 4: NSNP, 5: HETHET, 6: IBS0, 7: KINSHIP
                sampleA = parts[0]
                sampleB = parts[2]
                nsnp_str = parts[4]
                hethet_str = parts[5]
                ibs0_str = parts[6]
                kinship_str = parts[7]  # kinship coefficient (may be "NA" or "nan")

                # Skip pairs with NaN/NA kinship values from PLINK2
                if kinship_str.lower() in ["na", "nan", ""]:
                    continue

                try:
                    kinship = float(kinship_str)
                    nsnp = int(nsnp_str)
                    hethet = float(hethet_str)
                    ibs0 = float(ibs0_str)
                    # Check for NaN after conversion
                    if np.isnan(kinship) or np.isnan(hethet) or np.isnan(ibs0):
                        continue
                except (ValueError, OverflowError):
                    # Skip invalid values
                    continue

                # Skip pairs with invalid NSNP or HETHET (would cause division issues)
                if nsnp <= 0 or hethet <= 0:
                    continue

                if have_n1:
                    n1_A = float(n1_dict.get(sampleA, 0.0))
                    n1_B = float(n1_dict.get(sampleB, 0.0))
                    # Wire format: sampleA sampleB ibs0 hethet nsnp kinship n1_A n1_B
                    # Client will accumulate ibs0/hethet/nsnp/n1 across chunks and recompute kinship.
                    result_str += f"{sampleA} {sampleB} {ibs0} {hethet} {nsnp} {kinship} {n1_A} {n1_B}\n"

                    # Additional debug format (8 columns) for offline investigation:
                    # sampleA sampleB ibs0 hethet nsnp kinship n1_A n1_B
                    debug_str += f"{sampleA} {sampleB} {ibs0} {hethet} {nsnp} {kinship} {n1_A} {n1_B}\n"
                else:
                    # Wire format: sampleA sampleB ibs0 hethet nsnp kinship
                    # Client will accumulate ibs0/hethet/nsnp across chunks and recompute kinship.
                    result_str += f"{sampleA} {sampleB} {ibs0} {hethet} {nsnp} {kinship}\n"

                    # Additional debug format (6 columns) for offline investigation:
                    # sampleA sampleB ibs0 hethet nsnp kinship
                    debug_str += f"{sampleA} {sampleB} {ibs0} {hethet} {nsnp} {kinship}\n"
    else:
        logger.warning("[Server KING] No KING result file found.")

    # Persist debug partial results for this KING session (optional offline analysis)
    if debug_str:
        debug_file = f"{session_path}/king_partial_debug.txt"
        try:
            with open(debug_file, "w") as df:
                if have_n1:
                    df.write("# sampleA sampleB IBS0 HETHET NSNP KINSHIP N1_A N1_B\n")
                else:
                    df.write("# sampleA sampleB IBS0 HETHET NSNP KINSHIP\n")
                df.write(debug_str)
        except OSError:
            # Don't fail the KING stage if debug writing fails
            pass

    return np.frombuffer(result_str.encode("utf-8"), dtype=np.uint8)
