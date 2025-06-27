# client/iterative_lr.py

import logging
import numpy as np


def handle_iterative_lr(client, parameters, config):
    """
    Iterative LR Approach:
      - On each iteration, the client sends one data_chunk (as a .tar with anonymized bed/bim/fam).
      - The server merges chunks, runs PLINK logistic regression, then returns lines of:
            anonSNP_ID p-value
      - The client accumulates these p-values in a dictionary, converting anonSNP_ID -> realSNP_ID
        using local chunk_snp_map if needed.
      - Once no more chunks remain, the client performs a final local step:
         1. For each realSNP_ID, we have multiple p-values from different chunk passes.
         2. Compare each p-value to a threshold to get a binary significance result.
         3. Take the majority vote of these chunk-based significance flags to produce a final label.
    """

    # If we've sent all chunks, do the final local significance update
    if client.current_chunk_idx >= len(client.chunk_files):
        logging.info(f"[Client {client.client_id}] No more LR chunks. Performing final local significance update.")

        # Example: for each SNP, we have a list of pvals in client.lr_pvals[snpID].
        # We do a majority vote based on threshold.
        threshold = config.get("p_threshold", 5e-8)  # typical GWAS threshold

        final_significance = {}
        for snp_id, pval_list in client.lr_pvals.items():
            # Convert each p-value to a significance label (True/False)
            sig_flags = [(p < threshold) for p in pval_list]
            # Majority vote
            num_sig = sum(sig_flags)
            majority_is_sig = (num_sig > (len(sig_flags) // 2))
            final_significance[snp_id] = majority_is_sig

        logging.info(f"[Client {client.client_id}] Final LR significance updated for {len(final_significance)} SNPs.")
        # You might store final_significance in client.lr_final or some structure for later.
        client.lr_final = final_significance

        return [], 1, {}

    # Otherwise, we still have a chunk to send
    chunk_index = client.current_chunk_idx
    tar_file = client.chunk_files[chunk_index]
    client.current_chunk_idx += 1

    # Read the chunk file as bytes
    with open(tar_file, "rb") as f:
        chunk_data = f.read()
    data_array = np.frombuffer(chunk_data, dtype=np.uint8)

    logging.info(f"[Client {client.client_id}] sending LR chunk {chunk_index + 1}/{len(client.chunk_files)}")

    # If server returned p-values, parse them and store
    if parameters and len(parameters) > 0:
        try:
            # e.g., each line: "anonSNPID p_value"
            lr_results_str = parameters[0].tobytes().decode("utf-8").strip()
            lines = lr_results_str.splitlines()
            logging.info(f"[Client {client.client_id}] Received {len(lines)} partial LR p-values from server.")

            # let us find the local id map for 'chunk_index'
            # chunk_snp_map[chunk_index][anonSNP] -> realSNP
            # For example: client.chunk_snp_map[chunk_index][anon_id] = old_id

            for line in lines:
                parts = line.strip().split()
                if len(parts) < 2:
                    continue
                anon_snp_id, pval_str = parts[0], parts[1]
                pval = float(pval_str)

                # Convert anon_snp_id -> real_snp_id if we have a mapping
                if chunk_index in client.chunk_snp_map:
                    local_map = client.chunk_snp_map[chunk_index]
                    real_snp_id = local_map.get(anon_snp_id, anon_snp_id)
                else:
                    # If no mapping, we store as-is (less ideal if we need real ID)
                    real_snp_id = anon_snp_id

                # Accumulate p-values
                if real_snp_id not in client.lr_pvals:
                    client.lr_pvals[real_snp_id] = []
                client.lr_pvals[real_snp_id].append(pval)

                # Example:
                # client.lr_pvals = {
                #    "rs12345": [1e-7, 2e-5],
                #    "rs99999": [0.3, 0.15],
                #    ...
                # }
        except Exception as e:
                print(f"Error in Iterative LR records from server parameters: {e}")  
    return [data_array], 1, {}