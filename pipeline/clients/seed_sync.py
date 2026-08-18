from __future__ import annotations

import hashlib
import os
from typing import Dict, Optional

import numpy as np
from Crypto.PublicKey import ECC

from pipeline.clients.client_to_client import unpack_envelope, create_client_messenger


def get_persist_base_dir(log_dir: Optional[str], intermediate_dir: Optional[str]) -> str:
    base_dir = log_dir or intermediate_dir or os.getcwd()
    os.makedirs(base_dir, exist_ok=True)
    return base_dir


def key_paths(base_dir: str) -> tuple[str, str]:
    """Return (private_key_path, public_key_path) for persisting ECC keys across actor restarts."""
    os.makedirs(base_dir, exist_ok=True)
    return (
        os.path.join(base_dir, "ecc_private_key.pem"),
        os.path.join(base_dir, "ecc_public_key.pem"),
    )


def seed_path(base_dir: str) -> str:
    """Path for persisting computed global_seed across actor restarts."""
    os.makedirs(base_dir, exist_ok=True)
    return os.path.join(base_dir, "global_seed.txt")


def load_persisted_seed(base_dir: str) -> Optional[int]:
    try:
        p = seed_path(base_dir)
        if not os.path.exists(p):
            return None
        with open(p, "r") as f:
            s = f.read().strip()
        return int(s) if s else None
    except Exception:
        return None


def persist_seed(base_dir: str, seed: int) -> None:
    try:
        p = seed_path(base_dir)
        with open(p, "w") as f:
            f.write(str(int(seed)))
    except Exception:
        return


def load_persisted_keys_into_masking_helper(masking_helper, base_dir: str) -> Optional[str]:
    """
    Load persisted ECC private key into `masking_helper.private_key` if present.
    Returns our public key PEM if it can be derived/loaded.
    """
    try:
        priv_path, pub_path = key_paths(base_dir)
        if not os.path.exists(priv_path):
            return None

        with open(priv_path, "r") as f:
            priv_pem = f.read()
        priv_key = ECC.import_key(priv_pem)

        if masking_helper is not None:
            masking_helper.private_key = priv_key

        try:
            return priv_key.public_key().export_key(format="PEM")
        except Exception:
            if os.path.exists(pub_path):
                with open(pub_path, "r") as f:
                    return f.read()
    except Exception:
        return None

    return None


def resolve_my_hash_id(all_public_keys: Dict[str, str], my_public_key_pem: Optional[str]) -> Optional[int]:
    """
    Resolve this client's hash-based ID by matching our own public key PEM
    against the broadcast all_public_keys mapping (hash_id -> public_key_pem).
    """
    try:
        if not my_public_key_pem or not all_public_keys:
            return None
        my_norm = my_public_key_pem.strip()
        for k, v in all_public_keys.items():
            if isinstance(v, str) and v.strip() == my_norm:
                return int(k)
    except Exception:
        return None
    return None


def compute_global_seed_deterministic(all_public_keys: Dict[str, str], round_num: int, local_seed: int) -> int:
    """
    Fallback: Compute global seed deterministically from public information.
    Used when encrypted seeds are not available.
    """
    try:
        if not all_public_keys:
            return int(local_seed % (10**9))

        client_ids = sorted([int(k) for k in all_public_keys.keys()])
        client_ids_str = "-".join(str(cid) for cid in client_ids)
        seed_input = f"{client_ids_str}-round{int(round_num)}-fedgwas-seed"

        seed_hash = hashlib.sha256(seed_input.encode()).digest()
        return int.from_bytes(seed_hash[:4], "big") % (10**9)
    except Exception:
        return int(local_seed % (10**9))


def compute_global_seed_from_encrypted_messages(
    *,
    parameters,
    all_public_keys: Dict[str, str],
    my_public_key_pem: Optional[str],
    private_key: Optional[ECC.EccKey],
    num_clients: int,
    local_seed: int,
    server_round: int = 0,
) -> int:
    """
    Compute global seed by decrypting encrypted seed messages received via server relay.
    Server forwards ciphertexts but cannot decrypt (no client private keys).
    """
    my_hash_id = resolve_my_hash_id(all_public_keys, my_public_key_pem)
    if my_hash_id is None or private_key is None:
        return compute_global_seed_deterministic(all_public_keys, server_round, local_seed)

    # Extract encrypted messages from parameters (Flower Parameters or list of ndarrays)
    try:
        from flwr.common import parameters_to_ndarrays
        if hasattr(parameters, "tensors"):
            ndarrays = parameters_to_ndarrays(parameters)
        else:
            ndarrays = parameters if isinstance(parameters, list) else [parameters]
    except Exception:
        ndarrays = parameters if isinstance(parameters, list) else [parameters]

    messenger = create_client_messenger(my_hash_id, num_clients, private_key)

    seeds_received: Dict[int, float] = {}
    for encrypted_array in ndarrays:
        if not isinstance(encrypted_array, np.ndarray):
            continue
        encrypted_msg = encrypted_array.tobytes()
        try:
            meta = unpack_envelope(encrypted_msg)
        except Exception:
            continue

        recipient_id = int(meta.get("recipient_id", -1))
        if recipient_id != my_hash_id:
            continue

        sender_id = int(meta.get("sender_id"))
        sender_key_pem = all_public_keys.get(str(sender_id))
        if not sender_key_pem:
            continue

        decrypted = messenger.decrypt_from_sender(encrypted_msg, sender_id, sender_key_pem)
        if isinstance(decrypted, np.ndarray):
            seeds_received[sender_id] = float(decrypted.ravel()[0])
        else:
            seeds_received[sender_id] = float(decrypted)

    if seeds_received:
        return int(float(sum(seeds_received.values())) % (10**9))

    return compute_global_seed_deterministic(all_public_keys, server_round, local_seed)


def ensure_global_seed(
    current_seed: Optional[int],
    base_dir: str,
    all_public_keys: Dict[str, str],
    server_round: int,
    local_seed: int,
) -> int:
    """
    Ensure a non-None global_seed:
    - prefer current_seed if already set
    - else try persisted seed from disk
    - else deterministic fallback from public info
    """
    if current_seed is not None:
        return int(current_seed)

    s = load_persisted_seed(base_dir)
    if s is not None:
        return int(s)

    return int(compute_global_seed_deterministic(all_public_keys, server_round, local_seed))


