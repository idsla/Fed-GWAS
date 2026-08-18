# client/iterative_king.py

import numpy as np
import os
import subprocess
from pathlib import Path
import hashlib
import hmac
import json

def handle_iterative_king(client, parameters, config):
    """
    Handle iterative KING analysis stage.
    """
    client.logger.info(f"[Client {client.client_id}] === ITERATIVE_KING FUNCTION CALLED ===")

    def _stable_id(shared_secret: bytes, orig_id: str) -> str:
        """Derive a stable pseudonym for an original ID using a shared secret (client-only)."""
        if not shared_secret:
            return orig_id
        digest = hmac.new(shared_secret, orig_id.encode("utf-8"), hashlib.sha256).hexdigest()
        # 128-bit hex is enough and keeps IDs compact
        return f"s{digest[:32]}"

    def _ensure_shared_secrets(all_public_keys: dict):
        """Ensure client has shared secrets for peer HMAC derivations."""
        try:
            if not all_public_keys:
                return
            if getattr(client, "masking_helper", None) is None:
                return
            if not getattr(client.masking_helper, "shared_secrets", None):
                # compute shared secrets once per stage if missing
                client.masking_helper.compute_shared_secrets(all_public_keys, config)
        except Exception:
            return
    
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

    def _compute_linfit_coeffs():
        """Compute linear coefficients for kinship formula: y = a + b*(ibs0/hethet).
        Kept as a fallback when per-sample n1 counts are unavailable.
        """
        lf = getattr(client, "king_linfit", None) or {}
        n = float(lf.get("n", 0))
        if n < 2:
            return None
        sum_x = float(lf.get("sum_x", 0.0))
        sum_y = float(lf.get("sum_y", 0.0))
        sum_xx = float(lf.get("sum_xx", 0.0))
        sum_xy = float(lf.get("sum_xy", 0.0))
        denom = n * sum_xx - sum_x * sum_x
        if denom == 0:
            return None
        b = (n * sum_xy - sum_x * sum_y) / denom
        a = (sum_y - b * sum_x) / n
        return a, b
    
    # Load persisted KING accumulator if in-memory state is empty (actor restarts)
    def _load_king_accumulator():
        try:
            import json
            import pickle
            log_dir = getattr(client, "log_dir", "logs")
            state_path = os.path.join(log_dir, "king_accumulator_state.pkl")
            if os.path.exists(state_path):
                with open(state_path, "rb") as f:
                    state_data = pickle.load(f)
                    # Handle both old format (just accumulator) and new format (dict with both)
                    if isinstance(state_data, dict):
                        client.king_accumulator = state_data.get('accumulator', {})
                        client.king_position_to_pair_key = state_data.get('position_to_pair_key', {})
                        client.king_stable_to_orig = state_data.get('stable_to_orig', {})
                        client.king_linfit = state_data.get('linfit', {"n": 0, "sum_x": 0.0, "sum_y": 0.0, "sum_xx": 0.0, "sum_xy": 0.0})
                    else:
                        # Old format: just the accumulator
                        client.king_accumulator = state_data
                        client.king_position_to_pair_key = {}
                        client.king_stable_to_orig = {}
                        client.king_linfit = {"n": 0, "sum_x": 0.0, "sum_y": 0.0, "sum_xx": 0.0, "sum_xy": 0.0}
                client.logger.info(
                    f"[Client {client.client_id}] KING: Loaded persisted accumulator ({len(client.king_accumulator)} pairs), "
                    f"position_to_pair_key ({len(getattr(client, 'king_position_to_pair_key', {}))} position mappings)"
                )
        except Exception as e:
            client.logger.warning(f"[Client {client.client_id}] KING: Failed to load persisted accumulator: {e}")

    def _persist_king_accumulator():
        try:
            import pickle
            log_dir = getattr(client, "log_dir", "logs")
            os.makedirs(log_dir, exist_ok=True)
            state_path = os.path.join(log_dir, "king_accumulator_state.pkl")
            # Persist both accumulator and stable key map
            state_data = {
                'accumulator': client.king_accumulator,
                'position_to_pair_key': getattr(client, 'king_position_to_pair_key', {}),
                'stable_to_orig': getattr(client, 'king_stable_to_orig', {}),
                'linfit': getattr(client, 'king_linfit', {"n": 0, "sum_x": 0.0, "sum_y": 0.0, "sum_xx": 0.0, "sum_xy": 0.0}),
            }
            with open(state_path, "wb") as f:
                pickle.dump(state_data, f)
        except Exception as e:
            client.logger.warning(f"[Client {client.client_id}] KING: Failed to persist accumulator: {e}")

    # Load accumulator from disk if in-memory state is empty
    if not hasattr(client, 'king_accumulator') or len(getattr(client, 'king_accumulator', {})) == 0:
        _load_king_accumulator()
    
    # Initialize accumulator if it still doesn't exist
    if not hasattr(client, 'king_accumulator'):
        client.king_accumulator = {}
    
    # Ensure position map is loaded (it should be loaded by _load_king_accumulator, but initialize if missing)
    if not hasattr(client, 'king_position_to_pair_key'):
        client.king_position_to_pair_key = {}
    # Map stable pseudonyms -> original IDs for local samples (used for filtering and display)
    if not hasattr(client, 'king_stable_to_orig'):
        client.king_stable_to_orig = {}
    if not hasattr(client, 'king_linfit'):
        # Fallback linear fit (only used if per-sample n1 counts are unavailable)
        client.king_linfit = {"n": 0, "sum_x": 0.0, "sum_y": 0.0, "sum_xx": 0.0, "sum_xy": 0.0}
    
    # If we received KING results from server, update accumulator
    if parameters and len(parameters) > 0:
        try:
            from flwr.common import parameters_to_ndarrays
            ndarrays = parameters_to_ndarrays(parameters) if hasattr(parameters, "tensors") else parameters
            if not isinstance(ndarrays, list):
                ndarrays = list(ndarrays)
            
            if ndarrays and len(ndarrays) > 0:
                king_results_bytes = ndarrays[0].tobytes()
                king_results_str = king_results_bytes.decode("utf-8").strip()
                extra_ndarrays = ndarrays[1:] if len(ndarrays) > 1 else []

                # Determine iteration_id for this KING round
                stored_round = _get_state_value("king_chunk_round", None)
                iteration_id = int(stored_round) if stored_round is not None else config.get("server_round", 0)

                # Decrypt per-peer anon->stable maps (if provided)
                peer_anon_index = {}  # anon_id -> (peer_id, stable_id)
                payloads_seen = 0
                payload_entries = 0
                mismatched_iters = 0
                if extra_ndarrays:
                    try:
                        from pipeline.clients.flwr_config import get_all_public_keys
                        from pipeline.clients.seed_sync import resolve_my_hash_id
                        from pipeline.clients.client_to_client import create_client_messenger, unpack_envelope
                        from pipeline.clients.c2c_payloads import unpack_typed_payload_uint8

                        all_public_keys = get_all_public_keys(config)
                        _ensure_shared_secrets(all_public_keys)

                        my_hash_id = resolve_my_hash_id(all_public_keys, getattr(client, "my_public_key_pem", None))
                        if my_hash_id is not None and getattr(client, "masking_helper", None) is not None and getattr(client.masking_helper, "private_key", None) is not None:
                            messenger = create_client_messenger(my_hash_id, client.num_clients, client.masking_helper.private_key)
                            for enc_arr in extra_ndarrays:
                                if not isinstance(enc_arr, np.ndarray):
                                    continue
                                enc_msg = enc_arr.tobytes()
                                try:
                                    meta = unpack_envelope(enc_msg)
                                except Exception:
                                    continue
                                if int(meta.get("recipient_id", -1)) != my_hash_id:
                                    continue
                                sender_id = int(meta.get("sender_id"))
                                sender_key_pem = all_public_keys.get(str(sender_id))
                                if not sender_key_pem:
                                    continue
                                dec_u8 = messenger.decrypt_from_sender(enc_msg, sender_id, sender_key_pem)
                                if not isinstance(dec_u8, np.ndarray) or dec_u8.dtype != np.uint8:
                                    continue
                                kind, arr = unpack_typed_payload_uint8(dec_u8)
                                if kind != "king_map":
                                    continue
                                try:
                                    payload = json.loads(arr.tobytes().decode("utf-8"))
                                except Exception:
                                    continue
                                payload_iter = int(payload.get("iter_id", -1))
                                if payload_iter != int(iteration_id):
                                    mismatched_iters += 1
                                mapping_list = payload.get("map", [])
                                payloads_seen += 1
                                payload_entries += len(mapping_list)
                                for pair in mapping_list:
                                    if not isinstance(pair, (list, tuple)) or len(pair) < 2:
                                        continue
                                    anon_id, stable_id = str(pair[0]), str(pair[1])
                                    peer_anon_index[anon_id] = (sender_id, stable_id)
                    except Exception as e:
                        client.logger.warning(f"[Client {client.client_id}] KING: Failed to decrypt peer maps: {e}")
                if extra_ndarrays:
                    client.logger.info(
                        f"[Client {client.client_id}] KING: peer map payloads={payloads_seen}, "
                        f"entries={payload_entries}, mismatched_iters={mismatched_iters}, "
                        f"peer_anon_index={len(peer_anon_index)}"
                    )
                
                if king_results_str:
                    
                    # Build a combined sample map across all chunks to avoid index mismatches
                    sample_map = {}
                    total_entries = 0
                    if hasattr(client, 'chunk_sample_map') and isinstance(client.chunk_sample_map, dict):
                        for idx_map, chunk_map in client.chunk_sample_map.items():
                            if isinstance(chunk_map, dict):
                                sample_map.update(chunk_map)
                                total_entries += len(chunk_map)
                    # Fallback: load sample maps from disk if in-memory map is empty
                    # NOTE: Only load from this client's intermediate directory (privacy: clients can't access other clients' data)
                    if total_entries == 0:
                        try:
                            import glob
                            from pathlib import Path
                            stored_round = _get_state_value("king_chunk_round", None)
                            iteration_id = int(stored_round) if stored_round is not None else config.get("server_round", 0)
                            # Use intermediate_dir directly if available, otherwise infer from log_dir
                            interm_dir = getattr(client, 'intermediate_dir', None)
                            if interm_dir is None:
                                # Fallback: infer intermediate directory from log_dir
                                log_dir = getattr(client, 'log_dir', 'logs')
                                interm_dir = Path(log_dir).parent / "intermediate"
                            else:
                                interm_dir = Path(interm_dir)
                            # Resolve to absolute path to handle relative paths correctly
                            interm_dir = interm_dir.resolve()
                            # Prefer sample maps matching the current iteration id
                            map_files = glob.glob(str(interm_dir / f"*iter{iteration_id}_sample_map.tsv"))
                            allow_mixed = str(config.get("king_allow_mixed_maps", "0")).lower() in ("1", "true", "yes")
                            if not map_files and allow_mixed:
                                # Fallback to any sample maps if iteration-specific maps not found
                                map_files = glob.glob(str(interm_dir / "*sample_map.tsv"))
                            client.logger.info(f"[Client {client.client_id}] Looking for sample maps in: {interm_dir}")
                            client.logger.info(f"[Client {client.client_id}] Found {len(map_files)} sample map files: {[Path(mf).name for mf in map_files]}")
                            if not map_files:
                                client.logger.warning(
                                    f"[Client {client.client_id}] No sample map files for iter={iteration_id}; "
                                    f"set king_allow_mixed_maps=1 to allow mixed-map fallback"
                                )
                            for mf in map_files:
                                with open(mf, "r") as f:
                                    _ = f.readline()  # header
                                    for line in f:
                                        parts = line.strip().split()
                                        if len(parts) >= 2:
                                            anon_id, orig_id = parts[0], parts[1]
                                            if anon_id not in sample_map:
                                                sample_map[anon_id] = orig_id
                            total_entries = len(sample_map)
                            client.logger.info(f"[Client {client.client_id}] Loaded sample maps from disk: {len(map_files)} files, {total_entries} entries")
                            client.logger.info(f"[Client {client.client_id}] Loaded sample maps from disk: {len(map_files)} files, {total_entries} entries")
                        except Exception as e:
                            client.logger.warning(f"[Client {client.client_id}] Failed to load sample maps from disk: {e}")
                    if total_entries > 0:
                        client.logger.info(
                            f"[Client {client.client_id}] Using combined sample map with total entries={total_entries}"
                        )
                    else:
                        client.logger.warning(f"[Client {client.client_id}] No sample maps found (in-memory or disk), cannot de-anonymize KING results")
                    
                    # Parse KING results: format is "sampleA sampleB ibs0 hethet nsnp kinship\n"
                    # sampleA and sampleB are ANONYMIZED IDs from the server
                    lines = king_results_str.splitlines()
                    client.logger.info(f"[Client {client.client_id}] Received {len(lines)} KING result pairs from server (anonymized)")
                    
                    deanon_count = 0
                    skipped_count = 0
                    # Debug previews
                    try:
                        sample_map_keys_preview = list(sample_map.keys())[:5]
                        client.logger.info(f"[Client {client.client_id}] Sample map key preview: {sample_map_keys_preview}")
                    except Exception:
                        pass
                    
                    for line in lines:
                        parts = line.strip().split()
                        if len(parts) < 6:
                            continue
                        anon_sampleA, anon_sampleB = parts[0], parts[1]
                        ibs0_str, hethet_str, nsnp_str, kinship_str = parts[2], parts[3], parts[4], parts[5]
                        try:
                            ibs0 = float(ibs0_str)
                            hethet = float(hethet_str)
                            nsnp = float(nsnp_str)
                            kinship = float(kinship_str)
                            n1_A = float(parts[6]) if len(parts) >= 8 else 0.0
                            n1_B = float(parts[7]) if len(parts) >= 8 else 0.0
                            
                            # De-anonymize sample IDs using the local sample map
                            # If not found in local map, keep the anonymized ID (it's from another client)
                            orig_sampleA = sample_map.get(anon_sampleA, anon_sampleA)
                            orig_sampleB = sample_map.get(anon_sampleB, anon_sampleB)
                            local_A = anon_sampleA in sample_map
                            local_B = anon_sampleB in sample_map

                            # Same-client pairs: use original IDs (stable and readable)
                            swap_pair = False
                            if local_A and local_B:
                                if str(orig_sampleA) <= str(orig_sampleB):
                                    pair_key = (orig_sampleA, orig_sampleB)
                                    swap_pair = False
                                else:
                                    pair_key = (orig_sampleB, orig_sampleA)
                                    swap_pair = True
                            else:
                                # Cross-client pairs: use stable pseudonyms (client-only) from peer maps
                                peer_A = peer_anon_index.get(anon_sampleA)
                                peer_B = peer_anon_index.get(anon_sampleB)

                                stable_A = None
                                stable_B = None
                                peer_id_A = None
                                peer_id_B = None

                                if peer_A:
                                    peer_id_A, stable_A = peer_A
                                if peer_B:
                                    peer_id_B, stable_B = peer_B

                                # If one sample is local, derive its stable ID using the peer's shared secret
                                if local_A:
                                    if peer_id_B is None:
                                        skipped_count += 1
                                        continue
                                    secret = getattr(getattr(client, "masking_helper", None), "shared_secrets", {}).get(peer_id_B)
                                    if not secret:
                                        skipped_count += 1
                                        continue
                                    stable_A = _stable_id(secret, str(orig_sampleA))
                                    client.king_stable_to_orig[stable_A] = str(orig_sampleA)
                                if local_B:
                                    if peer_id_A is None:
                                        skipped_count += 1
                                        continue
                                    secret = getattr(getattr(client, "masking_helper", None), "shared_secrets", {}).get(peer_id_A)
                                    if not secret:
                                        skipped_count += 1
                                        continue
                                    stable_B = _stable_id(secret, str(orig_sampleB))
                                    client.king_stable_to_orig[stable_B] = str(orig_sampleB)

                                if stable_A is None or stable_B is None:
                                    skipped_count += 1
                                    continue

                                if stable_A <= stable_B:
                                    pair_key = (stable_A, stable_B)
                                    swap_pair = False
                                else:
                                    pair_key = (stable_B, stable_A)
                                    swap_pair = True
                            
                            # Align per-sample n1 counts to the ordered pair_key
                            if swap_pair:
                                n1_A, n1_B = n1_B, n1_A
                                anon_A_aligned, anon_B_aligned = anon_sampleB, anon_sampleA
                            else:
                                anon_A_aligned, anon_B_aligned = anon_sampleA, anon_sampleB
                            
                            if pair_key not in client.king_accumulator:
                                client.king_accumulator[pair_key] = {
                                    "sum_ibs0": 0.0,
                                    "sum_hethet": 0.0,
                                    "sum_nsnp": 0.0,
                                    "sum_n1A": 0.0,
                                    "sum_n1B": 0.0,
                                    "sum_phi_n1": 0.0,
                                    "sum_n1": 0.0,
                                    "sum_phi_hethet": 0.0,
                                    "anon_A": anon_A_aligned,  # Store aligned anonymized IDs for consistent saving
                                    "anon_B": anon_B_aligned
                                }
                            else:
                                # Backward compatibility for accumulators created before n1 fields existed
                                accum_existing = client.king_accumulator[pair_key]
                                if "sum_n1A" not in accum_existing:
                                    accum_existing["sum_n1A"] = 0.0
                                if "sum_n1B" not in accum_existing:
                                    accum_existing["sum_n1B"] = 0.0
                                if "sum_phi_n1" not in accum_existing:
                                    accum_existing["sum_phi_n1"] = 0.0
                                if "sum_n1" not in accum_existing:
                                    accum_existing["sum_n1"] = 0.0
                                if "sum_phi_hethet" not in accum_existing:
                                    accum_existing["sum_phi_hethet"] = 0.0
                            
                            # Accumulate sufficient stats across chunks (additive, NSNP-weighted)
                            # IBS0/HETHET are per-SNP rates in PLINK2 .kin0, so weight by NSNP.
                            client.king_accumulator[pair_key]["sum_ibs0"] += ibs0 * nsnp
                            client.king_accumulator[pair_key]["sum_hethet"] += hethet * nsnp
                            client.king_accumulator[pair_key]["sum_nsnp"] += nsnp
                            client.king_accumulator[pair_key]["sum_n1A"] += n1_A
                            client.king_accumulator[pair_key]["sum_n1B"] += n1_B
                            # Combine per-iteration phi using n1* weights (per paper Eq. 2)
                            # Use n1 of the first sample in the ordered pair as the weight.
                            if n1_A > 0:
                                client.king_accumulator[pair_key]["sum_phi_n1"] += kinship * n1_A
                                client.king_accumulator[pair_key]["sum_n1"] += n1_A
                            # Combine per-iteration phi using HETHET (n11) weights
                            # HETHET is a per-SNP rate; weight by NSNP to get counts.
                            if hethet > 0 and nsnp > 0:
                                client.king_accumulator[pair_key]["sum_phi_hethet"] += kinship * (hethet * nsnp)

                            # Update global linear fit for kinship formula: y = a + b*(ibs0/hethet)
                            # Use chunk-level kinship from PLINK2 as target
                            if hethet > 0:
                                x = ibs0 / hethet
                                lf = client.king_linfit
                                lf["n"] = int(lf.get("n", 0)) + 1
                                lf["sum_x"] = float(lf.get("sum_x", 0.0)) + x
                                lf["sum_y"] = float(lf.get("sum_y", 0.0)) + kinship
                                lf["sum_xx"] = float(lf.get("sum_xx", 0.0)) + x * x
                                lf["sum_xy"] = float(lf.get("sum_xy", 0.0)) + x * kinship
                            
                            # Track how many samples were successfully de-anonymized
                            if orig_sampleA != anon_sampleA:
                                deanon_count += 1
                            if orig_sampleB != anon_sampleB:
                                deanon_count += 1
                        except ValueError:
                            skipped_count += 1
                            continue
                    
                    pairs_processed = len(lines) - skipped_count
                    client.logger.info(f"[Client {client.client_id}] Processed {pairs_processed} KING pairs: de-anonymized {deanon_count} sample IDs from local chunks (kept anonymized IDs for samples from other clients)")
                    
                    # Persist accumulator after updating it
                    _persist_king_accumulator()
                    
                    # After receiving results, check if we're done and can finalize
                    # Note: We need to check the persisted state, not just in-memory state,
                    # because the client might have been restarted (actor restart)
                    restored_idx = _get_state_value("king_chunk_idx", None)
                    persisted_idx = int(restored_idx) if restored_idx is not None else 0
                    
                    # Also check if chunk_data exists (might be lost on actor restart)
                    has_chunk_data = hasattr(client, 'chunk_data') and len(getattr(client, 'chunk_data', [])) > 0
                    if has_chunk_data:
                        total_chunks = len(client.chunk_data)
                        client.logger.info(f"[Client {client.client_id}] After receiving results: persisted_idx={persisted_idx}, total_chunks={total_chunks}, has_chunk_data=True")
                    else:
                        # Chunk data lost, need to recreate - but we can't finalize yet
                        client.logger.info(f"[Client {client.client_id}] After receiving results: persisted_idx={persisted_idx}, chunk_data lost (will recreate), has_chunk_data=False")
                        # Don't finalize - fall through to recreate chunks and continue
                    
                    # Only finalize if we've truly processed all chunks (check persisted state)
                    # persisted_idx is the NEXT chunk to send, so if it equals total_chunks, we're done
                    if has_chunk_data and persisted_idx >= len(client.chunk_data):
                        # All chunks processed, finalize now
                        client.logger.info(f"[Client {client.client_id}] All KING chunks processed after receiving results (persisted_idx={persisted_idx} >= total_chunks={len(client.chunk_data)}), finalizing")
                        
                        # Check if accumulator has data
                        if not hasattr(client, 'king_accumulator') or len(client.king_accumulator) == 0:
                            client.logger.warning(f"[Client {client.client_id}] WARNING: KING accumulator is empty after receiving results!")
                            return [], 1, {}
                        
                        # Compute final phi for each pair using HETHET-weighted combination of per-iteration phi.
                        # If HETHET counts are missing, fall back to the legacy linear fit.
                        coeffs = _compute_linfit_coeffs()
                        if coeffs is None:
                            client.logger.warning(
                                f"[Client {client.client_id}] KING: insufficient data to fit linear formula fallback"
                            )
                        a_b = coeffs if coeffs is not None else (0.0, 0.0)
                        a, b = a_b
                        for pair_key, accum in client.king_accumulator.items():
                            sum_nsnp = float(accum.get("sum_nsnp", 0.0))
                            if sum_nsnp <= 0:
                                accum["phi"] = 0.0
                                continue
                            sum_hethet = float(accum.get("sum_hethet", 0.0))
                            sum_phi_hethet = float(accum.get("sum_phi_hethet", 0.0))
                            if sum_hethet > 0:
                                accum["phi"] = sum_phi_hethet / sum_hethet
                            else:
                                # Fallback to linear fit if HETHET counts are unavailable
                                sum_ibs0 = float(accum.get("sum_ibs0", 0.0))
                                ibs0_rate = sum_ibs0 / sum_nsnp
                                hethet_rate = sum_hethet / sum_nsnp if sum_nsnp > 0 else 0.0
                                if hethet_rate > 0:
                                    x = ibs0_rate / hethet_rate
                                    accum["phi"] = a + b * x
                                else:
                                    accum["phi"] = 0.0
                        
                        client.logger.info(f"[Client {client.client_id}] Final local KING update complete. "
                                     f"{len(client.king_accumulator)} sample pairs stored.")
                        
                        # Optional: Save KING results to file
                        try:
                            king_results_text = _save_king_results(client)
                        except Exception as e:
                            client.logger.warning(f"[Client {client.client_id}] KING results save failed: {e}")
                        
                        # Filter out samples beyond threshold
                        king_threshold = getattr(client, 'thresholds', {}).get("king_threshold", 0.4)
                        client.logger.info(f"[Client {client.client_id}] KING finalization: {len(client.king_accumulator)} pairs, threshold={king_threshold}")
                        try:
                            _filter_samples_by_king(client, king_threshold)
                        except Exception as filter_error:
                            client.logger.error(f"[Client {client.client_id}] Failed to filter samples by KING: {filter_error}")
                        
                        _persist_state_value("king_chunk_idx", len(getattr(client, "chunk_data", [])))
                        return [], 1, {
                            "king_done": True,
                            "king_chunks_remaining": 0,
                            "king_total_chunks": len(getattr(client, "chunk_data", [])),
                        }
                    else:
                        # More chunks to process, continue to send next chunk
                        if has_chunk_data:
                            client.logger.info(f"[Client {client.client_id}] More chunks to process (persisted_idx={persisted_idx}/{len(client.chunk_data)}), will send next chunk in next round")
                        else:
                            client.logger.info(f"[Client {client.client_id}] More chunks to process (persisted_idx={persisted_idx}), will recreate chunks and continue")
                        # Don't return here - fall through to send next chunk
        except Exception as e:
            client.logger.warning(f"[Client {client.client_id}] Failed to parse KING results from server: {e}")
    
    # Check if we have chunks to process
    if not hasattr(client, 'chunk_data') or len(client.chunk_data) == 0:
        client.logger.info(f"[Client {client.client_id}] No chunk data found, recreating chunks for KING analysis")
        
        # KING should chunk by SNPs to accumulate the same pair across SNP subsets
        client.partition_by = "snps"
        # Prefer snp_chunk_size from client config, fallback to server config, then default
        if hasattr(client, 'snp_chunk_size') and client.snp_chunk_size is not None:
            chunk_size = client.snp_chunk_size
        elif hasattr(client, 'parameters') and client.parameters and 'snp_chunk_size' in client.parameters:
            chunk_size = client.parameters['snp_chunk_size']
        else:
            # Fallback to server config or default
            chunk_size = config.get('snp_chunk_size', config.get('chunk_size', 1000))
        stored_round = _get_state_value("king_chunk_round", None)
        server_round = int(stored_round) if stored_round is not None else config.get('server_round', 0)
        client.logger.info(
            f"[Client {client.client_id}] Using chunk_size: {chunk_size}, server_round: {server_round} (for SNP partitioning)"
        )
        
        # Partition data and store in chunk_data
        # Use server_round to get different chunks at different iterations
        partition_config = {
            'chunk_size': chunk_size,
            'partition_by': 'snps',
            'server_round': server_round,
        }
        if hasattr(client, 'global_seed') and client.global_seed is not None:
            partition_config['global_seed'] = client.global_seed
            client.logger.info(f"[Client {client.client_id}] Using global_seed: {client.global_seed}, server_round: {server_round} (SNP chunks, stable sample IDs)")
        else:
            client.logger.warning(f"[Client {client.client_id}] No global_seed available, using default")
            
        chunk_data = client.partition_data(partition_config)
        client.chunk_data = chunk_data
        restored_idx = _get_state_value("king_chunk_idx", None)
        # Restore the index, but ensure it's valid (not beyond chunk_data length)
        if restored_idx is not None:
            restored_idx = int(restored_idx)
            if restored_idx > len(chunk_data):
                client.logger.warning(f"[Client {client.client_id}] Restored king_chunk_idx ({restored_idx}) > total chunks ({len(chunk_data)}), resetting to 0")
                restored_idx = 0
        else:
            restored_idx = 0
        client.current_chunk_idx = restored_idx
        _persist_state_value("king_chunk_round", int(server_round))
        
        client.logger.info(f"[Client {client.client_id}] Recreated {len(client.chunk_data)} chunks for KING analysis, restored_idx={restored_idx}")
    else:
        client.logger.info(f"[Client {client.client_id}] Using existing chunk data: {len(client.chunk_data)} chunks")
        if not hasattr(client, "current_chunk_idx"):
            restored_idx = _get_state_value("king_chunk_idx", None)
            client.current_chunk_idx = int(restored_idx) if restored_idx is not None else 0
    
    # Check if we've processed all chunks
    if client.current_chunk_idx >= len(client.chunk_data):
        client.logger.info(f"[Client {client.client_id}] All KING chunks processed, finalizing results")
        
        # Check if accumulator exists and has data
        if not hasattr(client, 'king_accumulator') or len(client.king_accumulator) == 0:
            client.logger.warning(f"[Client {client.client_id}] WARNING: KING accumulator is empty! No pairs to filter.")
            client.logger.info(f"[Client {client.client_id}] KING filtering skipped (no accumulator data). This may be expected if server doesn't send aggregated results back.")
            _persist_state_value("king_chunk_idx", len(client.chunk_data))
            return [], 1, {
                "king_done": True,
                "king_chunks_remaining": 0,
                "king_total_chunks": len(client.chunk_data),
            }
        
        # 1) Compute final phi for each pair using HETHET-weighted combination of per-iteration phi.
        # If HETHET counts are missing, fall back to the legacy linear fit.
        coeffs = _compute_linfit_coeffs()
        if coeffs is None:
            client.logger.warning(
                f"[Client {client.client_id}] KING: insufficient data to fit linear formula fallback"
            )
        a_b = coeffs if coeffs is not None else (0.0, 0.0)
        a, b = a_b
        for pair_key, accum in client.king_accumulator.items():
            sum_nsnp = float(accum.get("sum_nsnp", 0.0))
            if sum_nsnp <= 0:
                accum["phi"] = 0.0
                continue
            sum_hethet = float(accum.get("sum_hethet", 0.0))
            sum_phi_hethet = float(accum.get("sum_phi_hethet", 0.0))
            if sum_hethet > 0:
                accum["phi"] = sum_phi_hethet / sum_hethet
            else:
                # Fallback to linear fit if HETHET counts are unavailable
                sum_ibs0 = float(accum.get("sum_ibs0", 0.0))
                ibs0_rate = sum_ibs0 / sum_nsnp
                hethet_rate = sum_hethet / sum_nsnp if sum_nsnp > 0 else 0.0
                if hethet_rate > 0:
                    x = ibs0_rate / hethet_rate
                    accum["phi"] = a + b * x
                else:
                    accum["phi"] = 0.0

        client.logger.info(f"[Client {client.client_id}] Final local KING update complete. "
                     f"{len(client.king_accumulator)} sample pairs stored.")

        # Optional: Save KING results to file for post-analysis
        try:
            king_results_text = _save_king_results(client)
        except Exception as e:
            client.logger.warning(f"[Client {client.client_id}] KING results save failed: {e}")

        # 2) Filter out samples beyond threshold
        # Use client's thresholds (from config file) instead of server config
        king_threshold = getattr(client, 'thresholds', {}).get("king_threshold", 0.4)
        client.logger.info(f"[Client {client.client_id}] KING finalization: {len(client.king_accumulator)} pairs, threshold={king_threshold}")
        _filter_samples_by_king(client, king_threshold)

        _persist_state_value("king_chunk_idx", len(client.chunk_data))
        return [], 1, {
            "king_done": True,
            "king_chunks_remaining": 0,
            "king_total_chunks": len(client.chunk_data),
        }

    # Otherwise, proceed with the next chunk
    chunk_index = client.current_chunk_idx
    chunk_array = client.chunk_data[chunk_index]
    is_last_chunk = (chunk_index + 1) >= len(client.chunk_data)
    client.current_chunk_idx += 1
    _persist_state_value("king_chunk_idx", int(client.current_chunk_idx))

    client.logger.info(f"[Client {client.client_id}] Sending KING chunk "
                 f"{chunk_index+1}/{len(client.chunk_data)} {'(LAST CHUNK)' if is_last_chunk else ''}")
    
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
    
    # Build encrypted per-peer anon->stable maps for this chunk (client-only, server cannot decrypt)
    encrypted_maps = []
    try:
        from pipeline.clients.flwr_config import get_all_public_keys
        from pipeline.clients.seed_sync import resolve_my_hash_id
        from pipeline.clients.client_to_client import create_client_messenger
        from pipeline.clients.c2c_payloads import pack_typed_payload_uint8

        all_public_keys = get_all_public_keys(config)
        _ensure_shared_secrets(all_public_keys)
        my_hash_id = resolve_my_hash_id(all_public_keys, getattr(client, "my_public_key_pem", None))
        if all_public_keys and my_hash_id is not None and getattr(client, "masking_helper", None) is not None and getattr(client.masking_helper, "private_key", None) is not None:
            messenger = create_client_messenger(my_hash_id, client.num_clients, client.masking_helper.private_key)
            # Use the per-chunk sample map (anon -> original) to build peer-specific stable IDs
            chunk_sample_map = client.chunk_sample_map.get(chunk_index, {})
            if chunk_sample_map:
                map_entries = len(chunk_sample_map)
                payload_count = 0
                stored_round = _get_state_value("king_chunk_round", None)
                iteration_id = int(stored_round) if stored_round is not None else config.get("server_round", 0)
                for peer_id_str, peer_key_pem in all_public_keys.items():
                    peer_id = int(peer_id_str)
                    if peer_id == my_hash_id:
                        continue
                    secret = getattr(client.masking_helper, "shared_secrets", {}).get(peer_id)
                    if not secret:
                        continue
                    mapping = [[anon_id, _stable_id(secret, str(orig_id))] for anon_id, orig_id in chunk_sample_map.items()]
                    payload = {"iter_id": int(iteration_id), "chunk_index": int(chunk_index), "map": mapping}
                    payload_bytes = json.dumps(payload, separators=(",", ":")).encode("utf-8")
                    payload_u8 = np.frombuffer(payload_bytes, dtype=np.uint8)
                    typed_payload_u8 = pack_typed_payload_uint8("king_map", payload_u8)
                    enc = messenger.encrypt_for_recipient(typed_payload_u8, peer_id, peer_key_pem)
                    encrypted_maps.append(np.frombuffer(enc, dtype=np.uint8))
                    payload_count += 1
                client.logger.info(
                    f"[Client {client.client_id}] KING: built {payload_count} peer map payloads "
                    f"(entries={map_entries}, iter_id={iteration_id}, chunk_index={chunk_index})"
                )
    except Exception as e:
        client.logger.warning(f"[Client {client.client_id}] KING: Failed to build encrypted peer maps: {e}")

    remaining = max(len(client.chunk_data) - client.current_chunk_idx, 0)
    return [chunk_array] + encrypted_maps, 1, {
        "king_done": remaining == 0,
        "king_chunks_remaining": remaining,
        "king_total_chunks": len(client.chunk_data),
        "king_chunk_index": chunk_index,
    }

def _filter_samples_by_king(client, threshold):
    """
    Filter out samples from LOCAL dataset that exceed KING threshold.
    Only filters samples that are in the client's local dataset.
    """
    # Ensure threshold is a float (defensive: config values may be strings)
    threshold = float(threshold) if threshold is not None else 0.4
    
    # Load local sample IDs from .fam file and compute their anonymized versions
    # KING results contain anonymized FID values (from server's merged anonymized chunks)
    fam_file = client.plink_prefix + ".fam"
    if not os.path.exists(fam_file):
        client.logger.error(f"[Client {client.client_id}] Cannot find .fam file: {fam_file}")
        return
    
    # Load original FID/IID pairs from .fam file
    # The accumulator now uses ORIGINAL sample IDs as keys, so we can directly check against local samples
    local_original_samples = set()  # Set of original FIDs in local dataset
    original_fid_to_iid = {}  # orig_fid -> orig_iid mapping for PLINK filtering
    with open(fam_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                orig_fid, orig_iid = parts[0], parts[1]
                local_original_samples.add(orig_fid)
                original_fid_to_iid[orig_fid] = orig_iid
    
    client.logger.info(f"[Client {client.client_id}] Loaded {len(local_original_samples)} local samples from .fam file")
    
    # Now filter: accumulator keys may be original IDs (local pairs), stable pseudonyms (cross-client),
    # or anonymized IDs. We need to de-anonymize them to check against local samples.
    # Load sample maps to de-anonymize accumulator keys and stable pseudonym map for local samples.
    sample_map = {}
    for center_dir in Path(client.log_dir).parent.glob("center_*/intermediate"):
        for sm_file in center_dir.glob("*_sample_map.tsv"):
            with open(sm_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        sample_map[parts[0]] = parts[1]
    stable_to_orig = getattr(client, "king_stable_to_orig", {}) or {}
    
    samples_to_remove = set()  # Set of original FIDs to remove
    debug_pairs_logged = 0
    debug_king_pairs = str(
        getattr(client, "parameters", {}).get(
            "debug_king_pairs", os.environ.get("FEDGWAS_DEBUG_KING_PAIRS", "0")
        )
    ).lower() in ("1", "true", "yes")
    for (anon_sampleA, anon_sampleB), accum in client.king_accumulator.items():
        phi_val = accum.get("phi", 0.0)
        if phi_val > threshold:
            # De-anonymize the accumulator keys (stable pseudonyms -> original, then anon map)
            orig_sampleA = stable_to_orig.get(anon_sampleA, sample_map.get(anon_sampleA, anon_sampleA))
            orig_sampleB = stable_to_orig.get(anon_sampleB, sample_map.get(anon_sampleB, anon_sampleB))
            
            # Only filter if the de-anonymized sample is in the local original samples
            if orig_sampleA in local_original_samples:
                samples_to_remove.add(orig_sampleA)
            if orig_sampleB in local_original_samples:
                samples_to_remove.add(orig_sampleB)
            if debug_king_pairs and debug_pairs_logged < 5:
                is_local_A = orig_sampleA in local_original_samples
                is_local_B = orig_sampleB in local_original_samples
                client.logger.info(
                    f"[Client {client.client_id}] KING pair over threshold: phi={phi_val:.4f} "
                    f"sampleA={orig_sampleA} local={is_local_A}, "
                    f"sampleB={orig_sampleB} local={is_local_B}"
                )
                debug_pairs_logged += 1

    if not samples_to_remove:
        client.logger.info(f"[Client {client.client_id}] No LOCAL samples exceed KING threshold {threshold}.")
        return

    client.logger.info(f"[Client {client.client_id}] Removing {len(samples_to_remove)} LOCAL samples with phi > {threshold} (out of {len(local_original_samples)} local samples).")

    # Write removal file with original FID/IID pairs for PLINK filtering
    log_dir = getattr(client, 'log_dir', 'logs')
    os.makedirs(log_dir, exist_ok=True)
    temp_remove = os.path.join(log_dir, "temp_remove_king.txt")
    with open(temp_remove, "w") as f:
        for orig_fid in samples_to_remove:
            orig_iid = original_fid_to_iid.get(orig_fid, orig_fid)  # Use FID as IID if not found
            # PLINK --remove format: FID IID
            f.write(f"{orig_fid}\t{orig_iid}\n")

    client.logger.info(f"[Client {client.client_id}] Wrote {len(samples_to_remove)} original FID/IID pairs to removal file for PLINK filtering")

    filtered_prefix = os.path.join(log_dir, "king_filtered")
    
    # Find PLINK binary path
    from pipeline.clients.base_client import find_plink_binary
    plink_binary = find_plink_binary()
    
    cmd = [
        plink_binary,
        "--bfile", client.plink_prefix,
        "--remove", temp_remove,
        "--make-bed",
        "--allow-no-sex",
        "--out", filtered_prefix
    ]
    try:
        subprocess.run(cmd, check=True)
        client.plink_prefix = filtered_prefix
        # Persist the filtered prefix to state so subsequent stages use it
        if hasattr(client, '_persist_plink_prefix'):
            client._persist_plink_prefix(filtered_prefix)
            client.logger.info(f"[Client {client.client_id}] KING filter => new prefix {filtered_prefix} (persisted to state)")
        else:
            client.logger.info(f"[Client {client.client_id}] KING filter => new prefix {filtered_prefix}")
    except subprocess.CalledProcessError as e:
        client.logger.error(f"[Client {client.client_id}] PLINK filtering failed: {e}")

    if os.path.exists(temp_remove):
        os.remove(temp_remove)
    # Optionally reload local samples
    if hasattr(client, 'local_samples'):
        client.local_samples.clear()
        client.load_local_samples()


def _save_king_results(client) -> str:
    """Save KING results to a text file and return the text for visualization.
    
    Format: sample1 sample2 ibs0 hethet nsnp phi
    (ibs0/hethet are NSNP-weighted rates)
    
    Uses de-anonymized IDs for local samples (if available), anonymized IDs for other clients' samples.
    This ensures each client's own samples are de-anonymized for readability.
    """
    log_dir = getattr(client, 'log_dir', 'logs')
    os.makedirs(log_dir, exist_ok=True)
    
    king_file = os.path.join(log_dir, f"king_results_{client.client_id}.txt")
    lines = []
    
    # Load local sample FIDs for checking which samples are local
    local_fids = set()
    fam_file = getattr(client, 'plink_prefix', '') + ".fam"
    if os.path.exists(fam_file):
        with open(fam_file, "r") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 1:
                    local_fids.add(parts[0])
    
    # Load sample maps to de-anonymize accumulator keys (which may be anonymized IDs)
    sample_map = {}
    for center_dir in Path(log_dir).parent.glob("center_*/intermediate"):
        for sm_file in center_dir.glob("*_sample_map.tsv"):
            with open(sm_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        sample_map[parts[0]] = parts[1]
    stable_to_orig = getattr(client, "king_stable_to_orig", {}) or {}
    
    for pair_key, accum in client.king_accumulator.items():
        phi = accum.get("phi", 0.0)
        sum_nsnp = float(accum.get("sum_nsnp", 0.0))
        sum_ibs0 = float(accum.get("sum_ibs0", 0.0))
        sum_hethet = float(accum.get("sum_hethet", 0.0))
        ibs0_rate = (sum_ibs0 / sum_nsnp) if sum_nsnp > 0 else 0.0
        hethet_rate = (sum_hethet / sum_nsnp) if sum_nsnp > 0 else 0.0
        
        # The pair_key is already de-anonymized (original IDs) if both samples were de-anonymized,
        # or a mix of original/anonymized if only one was de-anonymized, or both anonymized if neither was.
        # We want to use the pair_key directly since it already contains the best available IDs.
        if isinstance(pair_key, tuple) and len(pair_key) == 2:
            save_A, save_B = pair_key
        else:
            # Fallback: try to get from accumulator or use anonymized IDs
            save_A = accum.get("anon_A", str(pair_key[0]) if isinstance(pair_key, (tuple, list)) else str(pair_key))
            save_B = accum.get("anon_B", str(pair_key[1]) if isinstance(pair_key, (tuple, list)) and len(pair_key) > 1 else "UNKNOWN")
        
        # If the IDs in pair_key are stable pseudonyms for local samples, de-anonymize first
        if save_A in stable_to_orig:
            save_A = stable_to_orig.get(save_A, save_A)
        if save_B in stable_to_orig:
            save_B = stable_to_orig.get(save_B, save_B)

        # If the IDs in pair_key are anonymized numeric IDs, try to de-anonymize using the sample map
        if isinstance(save_A, str) and save_A.isdigit() and len(save_A) > 12:
            # Looks like an anonymized ID, try to de-anonymize
            deanon_A = sample_map.get(save_A, save_A)
            # Only use de-anonymized if it's a local sample
            if deanon_A in local_fids:
                save_A = deanon_A
        
        if isinstance(save_B, str) and save_B.isdigit() and len(save_B) > 12:
            # Looks like an anonymized ID, try to de-anonymize
            deanon_B = sample_map.get(save_B, save_B)
            # Only use de-anonymized if it's a local sample
            if deanon_B in local_fids:
                save_B = deanon_B
        
        lines.append(f"{save_A} {save_B} {ibs0_rate} {hethet_rate} {sum_nsnp} {phi}")
    
    king_text = "\n".join(lines)
    with open(king_file, "w") as f:
        f.write(king_text)
    
    client.logger.info(f"[Client {client.client_id}] Saved {len(lines)} KING results to {king_file}")
    return king_text
