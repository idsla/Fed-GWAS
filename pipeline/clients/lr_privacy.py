from __future__ import annotations

import hashlib
import os
from typing import List, Optional


def lr_token(global_seed: Optional[int], local_seed: int, snp_id: str, token_salt: str = "") -> str:
    """
    Tokenize an SNP identifier so the server can't learn plaintext SNP IDs.

    Notes:
    - Keyed by global_seed (server does not know it).
    - token_salt is typically round-scoped to prevent linkability across rounds.
    """
    gs = global_seed
    if gs is None:
        # last resort: still hides SNPs from server, but may break intersection across clients
        gs = int(local_seed % (10**9))
    msg = f"fedgwas-lr|{int(gs)}|{token_salt}|{snp_id}".encode("utf-8")
    return hashlib.sha256(msg).hexdigest()


def make_dummy_tokens(global_seed: Optional[int], local_seed: int, n: int, token_salt: str = "") -> List[str]:
    """Deterministically generate dummy tokens for padding (do not map to any SNP ID)."""
    gs = global_seed
    if gs is None:
        gs = int(local_seed % (10**9))
    out: List[str] = []
    for i in range(int(n)):
        msg = f"fedgwas-lr-dummy|{int(gs)}|{token_salt}|{i}".encode("utf-8")
        out.append(hashlib.sha256(msg).hexdigest())
    return out


def resolve_tokens_to_snp_ids(
    plink_prefix: str, tokens: List[str], global_seed: Optional[int], local_seed: int, token_salt: str = ""
) -> List[str]:
    """Map token set back to local SNP IDs by scanning the current .bim."""
    if not tokens:
        return []

    token_set = set(tokens)
    bim_file = plink_prefix + ".bim"
    if not os.path.exists(bim_file):
        return []

    snp_ids: List[str] = []
    with open(bim_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                snp_ids.append(parts[1])

    to_exclude: List[str] = []
    for snp in snp_ids:
        if lr_token(global_seed, local_seed, snp, token_salt=token_salt) in token_set:
            to_exclude.append(snp)
    return to_exclude


