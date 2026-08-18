# client/iterative_lr.py

import os
import math
import json
import numpy as np


def compute_lr_chunk_size(client, config):
    """
    Compute LR chunk size with optional target chunk count.
    Priority:
      1) lr_target_chunks (or lr_target_rounds) -> derived chunk size by sample count
      2) lr_sample_chunk_size
      3) sample_chunk_size
      4) config["sample_chunk_size"]
      5) config["chunk_size"] (server)
      6) default 1000
    Returns: (chunk_size, sample_count, target_chunks, used_target)
    """
    # Base chunk size fallbacks
    chunk_size = None
    if hasattr(client, "lr_sample_chunk_size") and client.lr_sample_chunk_size is not None:
        chunk_size = client.lr_sample_chunk_size
    elif hasattr(client, "sample_chunk_size") and client.sample_chunk_size is not None:
        chunk_size = client.sample_chunk_size
    elif hasattr(client, "parameters") and client.parameters:
        if "lr_sample_chunk_size" in client.parameters:
            chunk_size = client.parameters["lr_sample_chunk_size"]
        elif "sample_chunk_size" in client.parameters:
            chunk_size = client.parameters["sample_chunk_size"]
    if chunk_size is None:
        if hasattr(client, "parameters") and client.parameters and "chunk_size" in client.parameters:
            chunk_size = client.parameters["chunk_size"]
        else:
            chunk_size = config.get("chunk_size", 1000)

    # Optional target chunk count for equalized rounds
    target_chunks = None
    if hasattr(client, "parameters") and client.parameters:
        if "lr_target_chunks" in client.parameters:
            target_chunks = client.parameters.get("lr_target_chunks")
        elif "lr_target_rounds" in client.parameters:
            target_chunks = client.parameters.get("lr_target_rounds")
    if target_chunks is None:
        if "lr_target_chunks" in config:
            target_chunks = config.get("lr_target_chunks")
        elif "lr_target_rounds" in config:
            target_chunks = config.get("lr_target_rounds")

    try:
        target_chunks = int(target_chunks) if target_chunks is not None else None
    except Exception:
        target_chunks = None

    if target_chunks and target_chunks > 0:
        sample_count = _count_samples_from_fam(getattr(client, "plink_prefix", None))
        if sample_count is not None and sample_count > 0:
            derived = int(math.ceil(sample_count / float(target_chunks)))
            derived = max(1, derived)
            return derived, sample_count, target_chunks, True

    return int(chunk_size) if chunk_size is not None else 1000, None, None, False


def _count_samples_from_fam(plink_prefix):
    if not plink_prefix:
        return None
    fam_file = f"{plink_prefix}.fam"
    if not os.path.exists(fam_file):
        return None
    count = 0
    try:
        with open(fam_file, "r") as f:
            for line in f:
                if not line.strip():
                    continue
                parts = line.strip().split()
                if len(parts) >= 2:
                    count += 1
    except Exception:
        return None
    return count


def handle_iterative_lr(client, parameters, config):
    """
    Handle iterative LR analysis stage.
    """
    def _get_state_value(key, default=None):
        try:
            state = getattr(client, "client_state", None)
            if not state:
                return default
            paths_rec = state.config_records.get("client_paths")
            paths_dict = dict(paths_rec) if paths_rec else {}
            return paths_dict.get(key, default)
        except Exception:
            return default

    def _persist_state_value(key, value):
        try:
            state = getattr(client, "client_state", None)
            if not state:
                return
            from flwr.common import ConfigRecord
            paths_rec = state.config_records.get("client_paths")
            paths_dict = dict(paths_rec) if paths_rec else {}
            paths_dict[key] = value
            state.config_records["client_paths"] = ConfigRecord(paths_dict)
        except Exception:
            return

    # If server already sent LR results, process them (even in iterative_lr stage)
    if parameters and len(parameters) > 0 and config.get("stage") in ("iterative_lr", "iterative_lr_response"):
        _handle_lr_response(client, parameters, config)
        if config.get("stage") == "iterative_lr_response":
            _persist_state_value("lr_chunk_idx", len(getattr(client, "chunk_data", [])))
            return [], 1, {
                "lr_done": True,
                "lr_chunks_remaining": 0,
                "lr_total_chunks": len(getattr(client, "chunk_data", [])),
            }

    # Check if we have chunks to process
    if not hasattr(client, 'chunk_data') or len(client.chunk_data) == 0:
        client.logger.info(f"[Client {client.client_id}] No chunk data found, recreating chunks for LR analysis")
        
        # LR chunks by samples
        client.partition_by = "samples"
        # LR partitions by samples (different sample subsets across iterations)
        chunk_size, sample_count, target_chunks, used_target = compute_lr_chunk_size(client, config)
        if used_target:
            num_chunks = int(math.ceil(sample_count / float(chunk_size))) if sample_count else "unknown"
            client.logger.info(
                f"[Client {client.client_id}] Using target LR chunks={target_chunks} "
                f"(samples={sample_count}, chunk_size={chunk_size}, chunks={num_chunks})"
            )
        else:
            client.logger.info(f"[Client {client.client_id}] Using chunk_size: {chunk_size} (for sample partitioning)")
        
        # Partition data and store in chunk_data
        stored_round = _get_state_value("lr_chunk_round", None)
        server_round = int(stored_round) if stored_round is not None else config.get('server_round', 0)
        partition_config = {'chunk_size': chunk_size, 'partition_by': 'samples', 'server_round': server_round}
        if hasattr(client, 'global_seed') and client.global_seed is not None:
            partition_config['global_seed'] = client.global_seed
            client.logger.info(f"[Client {client.client_id}] Using global_seed: {client.global_seed}, server_round: {server_round}")
        else:
            client.logger.warning(f"[Client {client.client_id}] No global_seed available, using default")
            
        chunk_data = client.partition_data(partition_config)
        client.chunk_data = chunk_data
        restored_idx = _get_state_value("lr_chunk_idx", None)
        client.current_chunk_idx = int(restored_idx) if restored_idx is not None else 0
        _persist_state_value("lr_chunk_round", int(server_round))
        
        client.logger.info(f"[Client {client.client_id}] Recreated {len(client.chunk_data)} chunks for LR analysis")
    
    if not hasattr(client, "current_chunk_idx"):
        restored_idx = _get_state_value("lr_chunk_idx", None)
        client.current_chunk_idx = int(restored_idx) if restored_idx is not None else 0

    # Check if we've processed all chunks
    if client.current_chunk_idx >= len(client.chunk_data):
        client.logger.info(f"[Client {client.client_id}] All LR chunks processed, awaiting server LR results")
        _persist_state_value("lr_chunk_idx", len(client.chunk_data))
        return [], 1, {
            "lr_done": True,
            "lr_chunks_remaining": 0,
            "lr_total_chunks": len(client.chunk_data),
        }

    # Otherwise, proceed with the next chunk
    chunk_index = client.current_chunk_idx
    chunk_array = client.chunk_data[chunk_index]
    client.current_chunk_idx += 1
    _persist_state_value("lr_chunk_idx", int(client.current_chunk_idx))

    client.logger.info(f"[Client {client.client_id}] Sending LR chunk "
                 f"{chunk_index+1}/{len(client.chunk_data)}")
    
    debug_chunk_payloads = str(
        config.get("debug_chunk_payloads", os.environ.get("FEDGWAS_DEBUG_CHUNK_PAYLOADS", "0"))
    ).lower() in ("1", "true", "yes")
    if debug_chunk_payloads:
        client.logger.info(f"[Client {client.client_id}] Chunk array type: {type(chunk_array)}")
        client.logger.info(f"[Client {client.client_id}] Chunk array shape: {chunk_array.shape if hasattr(chunk_array, 'shape') else 'no shape'}")
        client.logger.info(f"[Client {client.client_id}] Chunk array size: {len(chunk_array)}")
        if len(chunk_array) > 0:
            client.logger.info(f"[Client {client.client_id}] First few bytes: {chunk_array[:20]}")
            client.logger.info(f"[Client {client.client_id}] Metadata (first 12 bytes): {chunk_array[:12]}")
            if len(chunk_array) >= 12:
                metadata = np.frombuffer(chunk_array[:12], dtype=np.uint32)
                client.logger.info(f"[Client {client.client_id}] File sizes - bed: {metadata[0]}, bim: {metadata[1]}, fam: {metadata[2]}")

    remaining = max(len(client.chunk_data) - client.current_chunk_idx, 0)
    return [chunk_array], 1, {
        "lr_done": remaining == 0,
        "lr_chunks_remaining": remaining,
        "lr_total_chunks": len(client.chunk_data),
    }


def _handle_lr_response(client, parameters, config):
    """
    Process server-returned LR p-values (anonymized), de-anonymize using chunk_snp_map,
    and maintain per-iteration significance labels.
    """
    if not parameters:
        client.logger.warning(f"[Client {client.client_id}] LR response: no parameters received")
        return [], 1, {}

    # Convert Flower Parameters to numpy arrays
    try:
        from flwr.common import parameters_to_ndarrays
        ndarrays = parameters_to_ndarrays(parameters) if hasattr(parameters, "tensors") else parameters
        if not isinstance(ndarrays, list):
            ndarrays = list(ndarrays) if ndarrays is not None else []
    except Exception as e:
        client.logger.error(f"[Client {client.client_id}] Failed to convert parameters: {e}")
        return [], 1, {}

    if not ndarrays or len(ndarrays) == 0:
        client.logger.warning(f"[Client {client.client_id}] LR response: no data in parameters")
        return [], 1, {}

    lr_results_str = ndarrays[0].tobytes().decode("utf-8").strip()
    lines = lr_results_str.splitlines()
    client.logger.info(f"[Client {client.client_id}] Received {len(lines)} LR p-values from server (iterative_lr_response)")

    # Build combined SNP map anon->orig from all chunk maps
    combined_map = {}
    if hasattr(client, "chunk_snp_map"):
        for _, snp_map in client.chunk_snp_map.items():
            combined_map.update(snp_map)
    # Fallback: load SNP maps from disk if in-memory map is empty
    if not combined_map:
        try:
            import glob
            from pathlib import Path
            # Prefer the persisted chunking round if available
            stored_round = None
            try:
                state = getattr(client, "client_state", None)
                if state:
                    paths_rec = state.config_records.get("client_paths")
                    paths_dict = dict(paths_rec) if paths_rec else {}
                    stored_round = paths_dict.get("lr_chunk_round")
            except Exception:
                stored_round = None
            iteration_id = int(stored_round) if stored_round is not None else config.get("server_round", 0)
            interm_dir = getattr(client, 'intermediate_dir', None)
            if interm_dir is None:
                log_dir = getattr(client, 'log_dir', 'logs')
                interm_dir = Path(log_dir).parent / "intermediate"
            else:
                interm_dir = Path(interm_dir)
            interm_dir = interm_dir.resolve()
            map_files = glob.glob(str(interm_dir / f"*iter{iteration_id}_snp_map.tsv"))
            allow_mixed = str(config.get("lr_allow_mixed_maps", "0")).lower() in ("1", "true", "yes")
            if not map_files and allow_mixed:
                map_files = glob.glob(str(interm_dir / "*snp_map.tsv"))
            client.logger.info(f"[Client {client.client_id}] LR: Looking for SNP maps in: {interm_dir}")
            client.logger.info(f"[Client {client.client_id}] LR: Found {len(map_files)} SNP map files: {[Path(mf).name for mf in map_files]}")
            if not map_files:
                client.logger.warning(
                    f"[Client {client.client_id}] LR: No SNP map files for iter={iteration_id}; de-anonymization may be incomplete"
                )
            for mf in map_files:
                with open(mf, "r") as f:
                    _ = f.readline()
                    for line in f:
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            anon_id, orig_id = parts[0], parts[1]
                            if anon_id not in combined_map:
                                combined_map[anon_id] = orig_id
            client.logger.info(f"[Client {client.client_id}] LR: Loaded SNP maps from disk: {len(map_files)} files, {len(combined_map)} entries")
        except Exception as e:
            client.logger.warning(f"[Client {client.client_id}] LR: Failed to load SNP maps from disk: {e}")

    # Ensure lr_pvals and lr_significance dicts exist (load from disk if present)
    if not hasattr(client, "lr_pvals"):
        client.lr_pvals = {}
    if not hasattr(client, "lr_significance"):
        client.lr_significance = {}
    if not hasattr(client, "lr_round_participation"):
        client.lr_round_participation = {}

    # Load persisted participation metadata (if present)
    try:
        part_path = os.path.join(getattr(client, "log_dir", "logs"), "lr_round_participation.json")
        if os.path.exists(part_path):
            with open(part_path, "r") as f:
                client.lr_round_participation = json.load(f)
    except Exception as e:
        client.logger.warning(f"[Client {client.client_id}] LR: Failed to load participation metadata: {e}")

    # Load persisted LR state if in-memory state is empty (actor restarts)
    def _load_lr_state():
        try:
            import json
            log_dir = getattr(client, "log_dir", "logs")
            state_path = os.path.join(log_dir, "lr_pvals_state.json")
            if not os.path.exists(state_path):
                return
            with open(state_path, "r") as f:
                data = json.load(f)
            pvals = data.get("lr_pvals", {})
            sigs = data.get("lr_significance", {})
            if isinstance(pvals, dict):
                client.lr_pvals.update({k: list(v) for k, v in pvals.items()})
            if isinstance(sigs, dict):
                client.lr_significance.update({k: list(v) for k, v in sigs.items()})
            client.logger.info(
                f"[Client {client.client_id}] LR: Loaded persisted state ({len(client.lr_pvals)} SNPs)"
            )
        except Exception as e:
            client.logger.warning(f"[Client {client.client_id}] LR: Failed to load persisted state: {e}")

    def _persist_lr_state():
        try:
            import json
            log_dir = getattr(client, "log_dir", "logs")
            os.makedirs(log_dir, exist_ok=True)
            state_path = os.path.join(log_dir, "lr_pvals_state.json")
            with open(state_path, "w") as f:
                json.dump(
                    {
                        "lr_pvals": client.lr_pvals,
                        "lr_significance": client.lr_significance,
                    },
                    f,
                )
        except Exception as e:
            client.logger.warning(f"[Client {client.client_id}] LR: Failed to persist state: {e}")

    if not client.lr_pvals and not client.lr_significance:
        _load_lr_state()

    def _persist_lr_participation():
        try:
            log_dir = getattr(client, "log_dir", "logs")
            os.makedirs(log_dir, exist_ok=True)
            part_path = os.path.join(log_dir, "lr_round_participation.json")
            with open(part_path, "w") as f:
                json.dump(client.lr_round_participation, f)
        except Exception as e:
            client.logger.warning(f"[Client {client.client_id}] LR: Failed to persist participation: {e}")

    # Optional: per-round participation metadata (JSON) in ndarrays[1]
    if len(ndarrays) > 1:
        try:
            meta_str = ndarrays[1].tobytes().decode("utf-8").strip()
            if meta_str:
                meta = json.loads(meta_str)
                round_id = int(meta.get("round", config.get("server_round", 0) or 0))
                participants = meta.get("participants", [])
                if isinstance(participants, list):
                    client.lr_round_participation[str(round_id)] = participants
                    _persist_lr_participation()
        except Exception as e:
            client.logger.warning(f"[Client {client.client_id}] LR: Failed to parse participation meta: {e}")

    # Use global_lr_threshold for significance classification in global LR
    # Fallback to p_threshold for backward compatibility
    global_lr_threshold = getattr(client, "global_lr_threshold", None)
    if global_lr_threshold is None:
        global_lr_threshold = getattr(client, "p_threshold", 1e-3)
    p_threshold = float(global_lr_threshold) if global_lr_threshold is not None else 1e-3
    iteration = int(config.get("server_round", 0))

    for line in lines:
        parts = line.strip().split()
        if len(parts) < 2:
            continue
        anon_snp_id, pval_str = parts[0], parts[1]
        try:
            pval = float(pval_str)
        except ValueError:
            continue

        real_snp_id = combined_map.get(anon_snp_id, anon_snp_id)

        # Accumulate p-values per SNP
        client.lr_pvals.setdefault(real_snp_id, []).append(pval)
        # Track significance per iteration
        sig = pval < p_threshold
        client.lr_significance.setdefault(real_snp_id, []).append((iteration, sig))

    # Persist updated LR state for cross-iteration aggregation
    _persist_lr_state()
    
    # Log completion of LR response processing
    processed_count = len([line for line in lines if len(line.strip().split()) >= 2])
    client.logger.info(
        f"[Client {client.client_id}] LR response processing complete. Processed {processed_count} p-values, "
        f"total unique SNPs in state: {len(client.lr_pvals)}"
    )

    # Always write final LR results when in iterative_lr_response stage (final results)
    # emit_lr_intermediate controls intermediate results during iterative processing
    is_final_response = config.get("stage") == "iterative_lr_response"
    emit_lr_intermediate = str(
        config.get("emit_lr_intermediate", os.environ.get("FEDGWAS_EMIT_LR_INTERMEDIATE", "0"))
    ).lower() in ("1", "true", "yes")
    
    if is_final_response or emit_lr_intermediate:
        try:
            screening_rule = getattr(client, "lr_screening_rule", config.get("lr_screening_rule", "any_iter"))
            screening_rule = str(screening_rule).strip().lower()
            out_path = os.path.join(client.log_dir, "lr_results_client_deanon.txt")
            with open(out_path, "w") as f:
                f.write("SNP\tp_values\titerations_significant\tp_min\tp_geo_mean\tsig_any\tsig_majority\tscreening_rule\tscreening_sig\n")
                for snp, pvals in client.lr_pvals.items():
                    sig_list = client.lr_significance.get(snp, [])
                    p_min = min(pvals) if pvals else float("nan")
                    if pvals and all(p > 0 for p in pvals):
                        p_geo = math.exp(sum(math.log(p) for p in pvals) / len(pvals))
                    else:
                        p_geo = float("nan")
                    sig_only = [bool(s) for _, s in sig_list] if sig_list else []
                    sig_any = any(sig_only) if sig_only else False
                    sig_majority = False
                    if sig_only:
                        sig_majority = sum(1 for s in sig_only if s) >= ((len(sig_only) + 1) // 2)
                    if screening_rule in ("any_iter", "any", "min"):
                        screening_sig = sig_any or (p_min < p_threshold if p_min == p_min else False)
                    elif screening_rule in ("majority", "maj"):
                        screening_sig = sig_majority
                    elif screening_rule in ("geo", "geo_mean", "geomean"):
                        screening_sig = (p_geo < p_threshold) if p_geo == p_geo else False
                    elif screening_rule in ("mean", "avg", "average"):
                        p_avg = sum(pvals) / len(pvals) if pvals else float("nan")
                        screening_sig = (p_avg < p_threshold) if p_avg == p_avg else False
                    else:
                        screening_sig = sig_any

                    f.write(
                        f"{snp}\t{pvals}\t{sig_list}\t{p_min}\t{p_geo}\t{sig_any}\t{sig_majority}\t{screening_rule}\t{screening_sig}\n"
                    )

            screening_path = os.path.join(client.log_dir, "lr_screening_significant.txt")
            with open(screening_path, "w") as f:
                for snp, pvals in client.lr_pvals.items():
                    sig_list = client.lr_significance.get(snp, [])
                    sig_only = [bool(s) for _, s in sig_list] if sig_list else []
                    sig_any = any(sig_only) if sig_only else False
                    p_min = min(pvals) if pvals else float("nan")
                    if screening_rule in ("any_iter", "any", "min"):
                        screening_sig = sig_any or (p_min < p_threshold if p_min == p_min else False)
                    elif screening_rule in ("majority", "maj"):
                        screening_sig = sum(1 for s in sig_only if s) >= ((len(sig_only) + 1) // 2) if sig_only else False
                    elif screening_rule in ("geo", "geo_mean", "geomean"):
                        p_geo = math.exp(sum(math.log(p) for p in pvals) / len(pvals)) if pvals and all(p > 0 for p in pvals) else float("nan")
                        screening_sig = (p_geo < p_threshold) if p_geo == p_geo else False
                    elif screening_rule in ("mean", "avg", "average"):
                        p_avg = sum(pvals) / len(pvals) if pvals else float("nan")
                        screening_sig = (p_avg < p_threshold) if p_avg == p_avg else False
                    else:
                        screening_sig = sig_any
                    if screening_sig:
                        f.write(f"{snp}\n")
            result_type = "final results" if is_final_response else "intermediate artifacts"
            client.logger.info(
                f"[Client {client.client_id}] LR: wrote {result_type} ({out_path}, {screening_path})"
            )
        except Exception as e:
            client.logger.warning(f"[Client {client.client_id}] LR: Failed to write intermediate artifacts: {e}")

    return [], 1, {}
