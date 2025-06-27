# client/iterative_king.py

import logging
import numpy as np
import os
import subprocess

def handle_iterative_king(client, parameters, config):
    """
    Iterative King Approach:
      - Each iteration:
        1) Send the next chunk (anonymized .tar) to the server.
        2) Receive partial KING results of form:
             sampleID1_ano sampleID2_ano partial_phi n1_star
          Accumulate them in `client.king_accumulator`.
      - Once no more chunks remain:
        1) Compute final phi = sum_phiN / sum_n1.
        2) Filter out samples whose final phi > threshold from local dataset.
    """

    # If we've already processed all chunks, finalize local kinship.
    if client.current_chunk_idx >= len(client.chunk_files):
        logging.info(f"[Client {client.client_id}] No more KING chunks. Performing final local kinship update.")

        # 1) Compute final phi for each pair
        for pair_key, accum in client.king_accumulator.items():
            sum_phiN = accum.get("sum_phiN", 0.0)
            sum_n1   = accum.get("sum_n1", 0.0)
            if sum_n1 > 0:
                accum["phi"] = sum_phiN / sum_n1
            else:
                accum["phi"] = 0.0

        logging.info(f"[Client {client.client_id}] Final local KING update complete. "
                     f"{len(client.king_accumulator)} sample pairs stored.")

        # 2) Filter out samples beyond threshold
        king_threshold = config.get("king_threshold", 0.4)
        _filter_samples_by_king(client, king_threshold)

        return [], 1, {}

    # Otherwise, proceed with the next chunk
    chunk_index = client.current_chunk_idx
    tar_file = client.chunk_files[chunk_index]
    client.current_chunk_idx += 1

    # Read chunk file as bytes
    with open(tar_file, "rb") as f:
        chunk_data = f.read()
    data_array = np.frombuffer(chunk_data, dtype=np.uint8)

    logging.info(f"[Client {client.client_id}] Sending KING chunk "
                 f"{chunk_index+1}/{len(client.chunk_files)}")

    # If the server returned partial KING results for this chunk
    if parameters and len(parameters) > 0:
        try:
            partial_str = parameters[0].tobytes().decode("utf-8").strip()
            lines = partial_str.splitlines()
            logging.info(f"[Client {client.client_id}] Received {len(lines)} partial KING records from server.")
        
            for line in lines:
                # Expect: sampleA_ano sampleB_ano partial_phi n1_star
                parts = line.strip().split()
                if len(parts) < 4:
                    continue
                sampleA_ano, sampleB_ano, phi_str, n1_str = parts
                phi_chunk = float(phi_str)
                n1_star   = float(n1_str)

                # Convert anonymized -> real ID if you keep them in chunk_sample_map.
                # If not, store as is.
                if chunk_index in client.chunk_sample_map:
                    realA = client.chunk_sample_map[chunk_index].get(sampleA_ano, sampleA_ano)
                    realB = client.chunk_sample_map[chunk_index].get(sampleB_ano, sampleB_ano)
                else:
                    realA, realB = sampleA_ano, sampleB_ano

                key = tuple(sorted([realA, realB]))
                if key not in client.king_accumulator:
                    client.king_accumulator[key] = {"sum_phiN": 0.0, "sum_n1": 0.0, "phi": 0.0}

                client.king_accumulator[key]["sum_phiN"] += phi_chunk * n1_star
                client.king_accumulator[key]["sum_n1"]   += n1_star
        except Exception as e:
                print(f"Error in KING records from server parameters: {e}")    

    return [data_array], 1, {}

def _filter_samples_by_king(client, threshold):
    """
    Example approach: remove both samples in any pair whose phi > threshold.
    Adjust logic as needed.
    """
    samples_to_remove = set()
    for (sA, sB), accum in client.king_accumulator.items():
        phi_val = accum.get("phi", 0.0)
        if phi_val > threshold:
            samples_to_remove.add(sA)
            samples_to_remove.add(sB)

    if not samples_to_remove:
        logging.info(f"[Client {client.client_id}] No samples exceed KING threshold {threshold}.")
        return

    logging.info(f"[Client {client.client_id}] Removing {len(samples_to_remove)} samples with phi > {threshold}.")

    temp_remove = "temp_remove_king.txt"
    with open(temp_remove, "w") as f:
        for sid in samples_to_remove:
            f.write(f"{sid}\t{sid}\n")

    filtered_prefix = "king_filtered"
    cmd = [
        "plink",
        "--bfile", client.plink_prefix,
        "--remove", temp_remove,
        "--make-bed",
        "--out", filtered_prefix
    ]
    try:
        subprocess.run(cmd, check=True)
        client.plink_prefix = filtered_prefix
    except subprocess.CalledProcessError as e:
        logging.error(f"[Client {client.client_id}] PLINK filtering failed: {e}")

    os.remove(temp_remove)
    # Optionally reload local samples
    if hasattr(client, 'local_samples'):
        client.local_samples.clear()
        client.load_local_samples()