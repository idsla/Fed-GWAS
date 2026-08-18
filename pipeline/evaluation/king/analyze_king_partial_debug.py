#!/usr/bin/env python3
"""
Analyze KING partial debug outputs vs centralized KING for specific sample pairs.

Usage:
  python -m pipeline.evaluation.king.analyze_king_partial_debug \\
      --results-dir experiments/correctness/tiny_even/results \\
      --baseline-dir experiments/correctness/tiny_even/data/tiny/centralized_baseline \\
      --pair F000300 F000320
"""

import argparse
import pickle
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd


def load_baseline_kinship(
    baseline_dir: Path,
) -> Dict[Tuple[str, str], float]:
    """Load centralized KING kinship (king.kin0) into a dict keyed by sorted (FID1, FID2)."""
    king_file = baseline_dir / "king.kin0"
    if not king_file.exists():
        raise FileNotFoundError(f"Baseline KING file not found: {king_file}")

    baseline: Dict[Tuple[str, str], float] = {}
    with king_file.open("r") as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split()
            if len(parts) < 8:
                continue
            fid1, fid2 = parts[0], parts[2]
            try:
                kinship = float(parts[7])
            except ValueError:
                continue
            key = tuple(sorted((fid1, fid2)))
            baseline[key] = kinship
    return baseline


def load_sample_map(sample_map_file: Path) -> Dict[str, str]:
    """Load a sample map file (anon_sample -> orig_sample) into a dict."""
    sample_map: Dict[str, str] = {}
    if not sample_map_file.exists():
        return sample_map
    with sample_map_file.open("r") as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                anon_id, orig_id = parts[0], parts[1]
                sample_map[anon_id] = orig_id
    return sample_map


def find_sample_maps_for_session(
    results_dir: Path,
    session_id: str,
) -> Dict[str, str]:
    """
    Find and load all sample maps that might correspond to a given server session.
    
    Returns a combined mapping from anonymized -> original IDs.
    """
    combined_map: Dict[str, str] = {}
    
    # Search in all center intermediate directories (prefer latest run_* subdir if present)
    for center_root in results_dir.glob("center_*"):
        interm_root = center_root / "intermediate"
        if not interm_root.exists():
            continue
        run_dirs = [p for p in interm_root.iterdir() if p.is_dir() and p.name.startswith("run_")]
        interm = max(run_dirs, key=lambda p: p.stat().st_mtime) if run_dirs else interm_root
        # Look for sample map files (they might be named with iteration info)
        # Try to find any sample maps that could correspond to this session
        for sample_map_file in interm.glob("*sample_map.tsv"):
            smap = load_sample_map(sample_map_file)
            # Merge into combined map (later entries override if there are conflicts)
            combined_map.update(smap)
    
    return combined_map


def load_partial_debug_for_pair(
    results_dir: Path,
    orig_fid1: str,
    orig_fid2: str,
) -> pd.DataFrame:
    """
    Load all partial debug entries for a given original pair (fid1, fid2) across all KING sessions.
    
    This function:
    1. Loads all sample maps from client intermediate directories
    2. For each server session's debug file, maps anonymized IDs back to original IDs
    3. Filters for rows where the original pair matches the target pair
    
    Each debug file has lines (new format):
      sampleA sampleB IBS0 HETHET NSNP KINSHIP
    (Old format supported for backward compatibility.)
    """
    server_intermediate = results_dir / "server" / "intermediate"
    if not server_intermediate.exists():
        raise FileNotFoundError(f"Server intermediate dir not found: {server_intermediate}")

    # Pre-load all sample maps (they're shared across sessions)
    print("Loading sample maps from client intermediate directories...")
    all_sample_maps: Dict[str, str] = {}
    for center_root in results_dir.glob("center_*"):
        interm_root = center_root / "intermediate"
        if not interm_root.exists():
            continue
        run_dirs = [p for p in interm_root.iterdir() if p.is_dir() and p.name.startswith("run_")]
        interm = max(run_dirs, key=lambda p: p.stat().st_mtime) if run_dirs else interm_root
        for sample_map_file in interm.glob("*sample_map.tsv"):
            smap = load_sample_map(sample_map_file)
            all_sample_maps.update(smap)
    print(f"  Loaded {len(all_sample_maps)} anonymized -> original ID mappings")

    rows: List[Dict[str, object]] = []
    target_pair = set((orig_fid1, orig_fid2))
    
    for session_dir in server_intermediate.iterdir():
        if not session_dir.is_dir():
            continue
        debug_file = session_dir / "king_partial_debug.txt"
        if not debug_file.exists():
            continue
        
        # For this session, try to find session-specific sample maps if they exist
        # (Some implementations might store per-session maps)
        session_sample_map = find_sample_maps_for_session(results_dir, session_dir.name)
        # Combine with global maps
        combined_map = {**all_sample_maps, **session_sample_map}
        
        with debug_file.open("r") as f:
            header = f.readline().strip()
            is_new_format = "IBS0" in header and "HETHET" in header
            for line in f:
                parts = line.strip().split()
                if is_new_format:
                    if len(parts) < 6:
                        continue
                    anonA, anonB = parts[0], parts[1]
                    try:
                        ibs0 = float(parts[2])
                        hethet = float(parts[3])
                        nsnp = float(parts[4])
                        kinship = float(parts[5])
                    except ValueError:
                        continue
                else:
                    # Old format (backward compatibility): sampleA sampleB partial_phi n1_A n1_B NSNP
                    if len(parts) < 5:
                        continue
                    anonA, anonB = parts[0], parts[1]
                    try:
                        kinship = float(parts[2])
                    except ValueError:
                        continue

                # Map anonymized IDs back to original IDs
                origA = combined_map.get(anonA, anonA)  # Fallback to anon if not found
                origB = combined_map.get(anonB, anonB)

                # Check if this row matches our target original pair (unordered)
                if set((origA, origB)) != target_pair:
                    continue

                if is_new_format:
                    rows.append(
                        {
                            "session": session_dir.name,
                            "anon_FID1": anonA,
                            "anon_FID2": anonB,
                            "orig_FID1": origA,
                            "orig_FID2": origB,
                            "ibs0": ibs0,
                            "hethet": hethet,
                            "nsnp": nsnp,
                            "kinship": kinship,
                        }
                    )
                else:
                    # Old format parsing (backward compatibility)
                    try:
                        if len(parts) >= 6:
                            n1_A = float(parts[3])
                            n1_B = float(parts[4])
                            nsnp = int(parts[5])
                            if str(origA) < str(origB):
                                n1_first = n1_A
                            else:
                                n1_first = n1_B
                        else:
                            n1_first = float(parts[3])
                            nsnp = int(parts[4])
                            n1_A = n1_first
                            n1_B = n1_first
                    except ValueError:
                        continue

                    rows.append(
                        {
                            "session": session_dir.name,
                            "anon_FID1": anonA,
                            "anon_FID2": anonB,
                            "orig_FID1": origA,
                            "orig_FID2": origB,
                            "partial_phi": kinship,
                            "n1_first": n1_first,
                            "n1_A": n1_A if len(parts) >= 6 else n1_first,
                            "n1_B": n1_B if len(parts) >= 6 else n1_first,
                            "NSNP": nsnp,
                        }
                    )

    if not rows:
        raise ValueError(
            f"No partial debug entries found for original pair ({orig_fid1}, {orig_fid2}) "
            f"across all sessions"
        )

    df = pd.DataFrame(rows)
    return df


def compute_linfit_from_debug(results_dir: Path) -> Optional[Tuple[float, float]]:
    """Compute linear fit coefficients for kinship from debug files.

    Uses: kinship = a + b * (IBS0 / HETHET).
    Returns (a, b) or None if insufficient data.
    """
    server_intermediate = results_dir / "server" / "intermediate"
    if not server_intermediate.exists():
        return None

    n = 0.0
    sum_x = sum_y = sum_xx = sum_xy = 0.0
    for debug_file in server_intermediate.glob("*/king_partial_debug.txt"):
        with debug_file.open("r") as f:
            header = f.readline().strip()
            if "IBS0" not in header or "HETHET" not in header:
                continue
            for line in f:
                parts = line.strip().split()
                if len(parts) < 6:
                    continue
                try:
                    ibs0 = float(parts[2])
                    hethet = float(parts[3])
                    kinship = float(parts[5])
                except ValueError:
                    continue
                if hethet <= 0:
                    continue
                x = ibs0 / hethet
                n += 1.0
                sum_x += x
                sum_y += kinship
                sum_xx += x * x
                sum_xy += x * kinship

    if n < 2:
        return None
    denom = n * sum_xx - sum_x * sum_x
    if denom == 0:
        return None
    b = (n * sum_xy - sum_x * sum_y) / denom
    a = (sum_y - b * sum_x) / n
    return a, b


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze KING partial debug vs centralized KING.")
    parser.add_argument(
        "--results-dir",
        required=True,
        help="Federated results directory (e.g., experiments/correctness/tiny_even/results)",
    )
    parser.add_argument(
        "--baseline-dir",
        required=True,
        help="Centralized baseline directory (e.g., experiments/correctness/tiny_even/data/tiny/centralized_baseline)",
    )
    parser.add_argument(
        "--pair",
        nargs=2,
        metavar=("FID1", "FID2"),
        required=True,
        help="Sample pair (FID1 FID2) to analyze",
    )
    parser.add_argument(
        "--center-id",
        type=int,
        default=1,
        help="Center ID whose accumulator to use for de-anonymization (default: 1)",
    )

    args = parser.parse_args()
    results_dir = Path(args.results_dir)
    baseline_dir = Path(args.baseline_dir)
    center_id: int = args.center_id
    fid1, fid2 = args.pair
    pair_key = tuple(sorted((fid1, fid2)))

    print("=" * 80)
    print(f"Analyzing KING partials for pair {pair_key}")
    print("=" * 80)

    # Load baseline KING
    baseline = load_baseline_kinship(baseline_dir)
    baseline_phi = baseline.get(pair_key, None)
    if baseline_phi is None:
        print(f"Baseline KING: no entry found for pair {pair_key}")
    else:
        print(f"Baseline KING φ: {baseline_phi:.8f}")

    # Load partial debug entries by mapping anonymized IDs back to original IDs
    # This will find all iterations where this original pair appears, regardless of anonymization
    print(f"\nSearching for original pair {pair_key} across all server sessions...")
    df = load_partial_debug_for_pair(results_dir, fid1, fid2)
    df = df.sort_values(by=["session"])
    print(f"\nFound {len(df)} partial entries across {df['session'].nunique()} unique sessions")

    if "ibs0" in df.columns:
        print("\nPer-iteration partials (session, IBS0, HETHET, NSNP, kinship):")
        for _, row in df.iterrows():
            print(
                f"  {row['session']}: IBS0={row['ibs0']:.3f}, HETHET={row['hethet']:.3f}, "
                f"NSNP={int(row['nsnp'])}, kinship={row['kinship']:.8f}, "
                f"anon_ids=({row['anon_FID1']}, {row['anon_FID2']})"
            )

        coeffs = compute_linfit_from_debug(results_dir)
        if coeffs is None:
            print("\nWARNING: insufficient data to fit linear formula; kinship aggregation will be NaN.")
            a, b = float("nan"), float("nan")
        else:
            a, b = coeffs
            print(f"\nLinear fit (global): kinship = a + b * (IBS0/HETHET)")
            print(f"  a={a:.6f}, b={b:.6f}")

        sum_ibs0 = float((df["ibs0"] * df["nsnp"]).sum())
        sum_hethet = float((df["hethet"] * df["nsnp"]).sum())
        sum_nsnp = float(df["nsnp"].sum())
        kinship_mean = float(df["kinship"].mean())
        kinship_nsnp = float((df["kinship"] * df["nsnp"]).sum() / sum_nsnp) if sum_nsnp > 0 else float("nan")
        if sum_hethet > 0 and not np.isnan(a):
            kinship_fit = a + b * (sum_ibs0 / sum_hethet)
        else:
            kinship_fit = float("nan")

        print("\nAggregated estimates for this pair:")
        print(f"  Sums: IBS0={sum_ibs0:.3f}, HETHET={sum_hethet:.3f}, NSNP={sum_nsnp:.0f}")
        print(f"  Mean(chunk kinship): {kinship_mean:.8f}")
        print(f"  NSNP-weighted mean: {kinship_nsnp:.8f}")
        print(f"  Linfit on aggregated IBS0/HETHET: {kinship_fit:.8f}")

        if baseline_phi is not None:
            print("\nDifferences vs baseline:")
            print(f"  |mean - baseline|   = {abs(kinship_mean - baseline_phi):.8f}")
            print(f"  |nsnp - baseline|   = {abs(kinship_nsnp - baseline_phi):.8f}")
            print(f"  |linfit - baseline| = {abs(kinship_fit - baseline_phi):.8f}")
    else:
        print("\nPer-iteration partials (session, partial_phi, n1_first, NSNP):")
        for _, row in df.iterrows():
            print(
                f"  {row['session']}: φ={row['partial_phi']:.8f}, "
                f"n1_first={row['n1_first']:.1f}, NSNP={int(row['NSNP'])}, "
                f"anon_ids=({row['anon_FID1']}, {row['anon_FID2']})"
            )

        phi = df["partial_phi"].to_numpy(dtype=float)
        n1_first = df["n1_first"].to_numpy(dtype=float)
        nsnp = df["NSNP"].to_numpy(dtype=float)

        sum_n1 = float(np.sum(n1_first))
        sum_nsnp = float(np.sum(nsnp))

        phi_n1 = float(np.sum(phi * n1_first) / sum_n1) if sum_n1 > 0 else float("nan")
        phi_nsnp = float(np.sum(phi * nsnp) / sum_nsnp) if sum_nsnp > 0 else float("nan")

        print("\nAggregated estimates for this pair:")
        print(f"  Federated rule (φ = Σ_k (φ_k * N_Aa^i_k) / Σ_k N_Aa^i_k): {phi_n1:.8f}")
        print(f"  NSNP-weighted (φ_NSNP = Σ φ_i * NSNP_i / Σ NSNP_i): {phi_nsnp:.8f}")

        if baseline_phi is not None:
            print("\nDifferences vs baseline:")
            print(f"  |φ_n1 - baseline|   = {abs(phi_n1 - baseline_phi):.8f}")
            print(f"  |φ_NSNP - baseline| = {abs(phi_nsnp - baseline_phi):.8f}")

    print("\nDone.")


if __name__ == "__main__":
    main()
