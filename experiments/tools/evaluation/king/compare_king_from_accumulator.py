#!/usr/bin/env python3
"""
Compare federated KING accumulator (client-side) to centralized baseline
without using result_analyzer.

This script uses stable_to_orig maps stored in each client's
king_accumulator_state.pkl to de-anonymize pairs, then computes correlation
vs baseline king.kin0.

Usage:
  python experiments/tools/evaluation/king/compare_king_from_accumulator.py \
    --results-dir experiments/real_world/1000genomes/results \
    --baseline-dir experiments/real_world/1000genomes/data/centralized_baseline \
    --center-id 1

Optionally pass --data-dir to locate id_map_center_*.csv when the data layout
differs (e.g., experiments/correctness/tiny_even/data/tiny).
"""

from __future__ import annotations

import argparse
import math
import pickle
from pathlib import Path
from typing import Dict, Tuple


def _load_accumulator(results_dir: Path, center_id: int):
    acc_path = results_dir / f"center_{center_id}" / "logs" / "king_accumulator_state.pkl"
    if not acc_path.exists():
        raise FileNotFoundError(f"Accumulator not found: {acc_path}")
    with acc_path.open("rb") as f:
        state = pickle.load(f)
    if isinstance(state, dict) and "accumulator" in state:
        acc = state["accumulator"]
        linfit = state.get("linfit", {})
        stable_map = state.get("stable_to_orig", {})
    else:
        acc = state
        linfit = {}
        stable_map = {}
    return acc, linfit, stable_map


def _load_all_stable_maps(results_dir: Path) -> Dict[str, str]:
    """Merge stable_to_orig maps from all centers."""
    merged: Dict[str, str] = {}
    for acc_file in results_dir.glob("center_*/logs/king_accumulator_state.pkl"):
        try:
            with acc_file.open("rb") as f:
                state = pickle.load(f)
            if isinstance(state, dict):
                stable = state.get("stable_to_orig", {})
                for k, v in stable.items():
                    merged.setdefault(str(k), str(v))
        except Exception:
            continue
    return merged


def _compute_linfit(linfit: Dict) -> Tuple[float, float] | None:
    try:
        n = float(linfit.get("n", 0))
        if n < 2:
            return None
        sum_x = float(linfit.get("sum_x", 0.0))
        sum_y = float(linfit.get("sum_y", 0.0))
        sum_xx = float(linfit.get("sum_xx", 0.0))
        sum_xy = float(linfit.get("sum_xy", 0.0))
        denom = n * sum_xx - sum_x * sum_x
        if denom == 0:
            return None
        b = (n * sum_xy - sum_x * sum_y) / denom
        a = (sum_y - b * sum_x) / n
        return a, b
    except Exception:
        return None


def _load_baseline(baseline_dir: Path) -> Dict[Tuple[str, str], float]:
    baseline_file = baseline_dir / "king.kin0"
    if not baseline_file.exists():
        raise FileNotFoundError(f"Baseline KING file not found: {baseline_file}")
    baseline: Dict[Tuple[str, str], float] = {}
    with baseline_file.open("r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 8:
                continue
            fid1, fid2 = parts[0], parts[2]
            try:
                kin = float(parts[7])
            except ValueError:
                continue
            key = (fid1, fid2) if fid1 <= fid2 else (fid2, fid1)
            baseline[key] = kin
    return baseline


def _load_id_map_ids(data_dir: Path, center_id: int) -> set[str]:
    csv_path = data_dir / f"center_{center_id}" / f"id_map_center_{center_id}.csv"
    if not csv_path.exists():
        return set()
    ids = set()
    with csv_path.open("r") as f:
        header = f.readline().strip().split(",")
        try:
            fid_idx = header.index("FID")
        except ValueError:
            fid_idx = None
        for line in f:
            parts = line.strip().split(",")
            if fid_idx is not None and fid_idx < len(parts):
                ids.add(parts[fid_idx])
    return ids


def _normalize_pair(a: str, b: str) -> Tuple[str, str]:
    return (a, b) if a <= b else (b, a)


def _resolve_data_dir(results_dir: Path, explicit: Path | None, center_id: int) -> Path:
    if explicit is not None:
        return explicit
    base = results_dir.parent / "data"
    # Direct layout: data/center_*
    if (base / f"center_{center_id}").exists():
        return base
    # Search for id_map under data/**/center_*
    try:
        matches = list(base.rglob(f"id_map_center_{center_id}.csv"))
    except Exception:
        matches = []
    if matches:
        # id_map_center_*.csv lives under .../<data_root>/center_*/id_map_center_*.csv
        return matches[0].parent.parent
    return base


def compare_from_accumulator(
    results_dir: Path,
    baseline_dir: Path,
    center_id: int,
    data_dir: Path | None = None,
    report_path: Path | None = None,
) -> dict:
    print(f"Loading accumulator from center_{center_id}...")
    acc, linfit, stable_map_local = _load_accumulator(results_dir, center_id)
    print(f"  accumulator pairs: {len(acc)}")

    print("Loading stable_to_orig maps from all centers...")
    stable_map = _load_all_stable_maps(results_dir)
    # Ensure local map entries win (if any conflicts)
    for k, v in stable_map_local.items():
        stable_map[str(k)] = str(v)
    print(f"  stable_to_orig entries: {len(stable_map)}")

    print("Loading baseline...")
    baseline = _load_baseline(baseline_dir)
    print(f"  baseline pairs: {len(baseline)}")

    # Load client ID sets for cross/within breakdown
    # Expected layout: experiments/<category>/<dataset>/results
    # data dir is experiments/<category>/<dataset>/data
    data_root = _resolve_data_dir(results_dir, data_dir, center_id)
    ids1 = _load_id_map_ids(data_root, 1)
    ids2 = _load_id_map_ids(data_root, 2)

    coeffs = _compute_linfit(linfit)
    total_pairs = 0
    mapped_pairs = 0
    missing_map_pairs = 0
    missing_baseline_pairs = 0
    cross_pairs = 0
    within_pairs = 0

    # Correlation accumulators
    n = 0
    sum_x = sum_y = sum_xx = sum_yy = sum_xy = 0.0
    sum_abs = 0.0
    max_abs = 0.0

    for (a, b), v in acc.items():
        total_pairs += 1
        sum_nsnp = float(v.get("sum_nsnp", 0.0))
        sum_hethet = float(v.get("sum_hethet", 0.0))
        sum_ibs0 = float(v.get("sum_ibs0", 0.0))
        phi = v.get("phi", None)
        if phi is None:
            sum_phi_hethet = float(v.get("sum_phi_hethet", 0.0))
            sum_hethet = float(v.get("sum_hethet", 0.0))
            if sum_hethet > 0:
                phi = sum_phi_hethet / sum_hethet
        if phi is None and coeffs is not None and sum_hethet > 0:
            a_coef, b_coef = coeffs
            phi = a_coef + b_coef * (sum_ibs0 / sum_hethet)
        if phi is None:
            continue

        a_str = str(a)
        b_str = str(b)
        a_m = stable_map.get(a_str)
        b_m = stable_map.get(b_str)
        # If ID is already original (e.g., HG*/NA*), keep it as-is.
        if a_m is None:
            if a_str.startswith("s"):
                missing_map_pairs += 1
                continue
            a_m = a_str
        if b_m is None:
            if b_str.startswith("s"):
                missing_map_pairs += 1
                continue
            b_m = b_str

        mapped_pairs += 1
        if (a_m in ids1 and b_m in ids2) or (a_m in ids2 and b_m in ids1):
            cross_pairs += 1
        else:
            within_pairs += 1

        key = _normalize_pair(a_m, b_m)
        baseline_phi = baseline.get(key)
        if baseline_phi is None:
            missing_baseline_pairs += 1
            continue

        n += 1
        sum_x += phi
        sum_y += baseline_phi
        sum_xx += phi * phi
        sum_yy += baseline_phi * baseline_phi
        sum_xy += phi * baseline_phi
        diff = phi - baseline_phi
        ad = diff if diff >= 0 else -diff
        sum_abs += ad
        if ad > max_abs:
            max_abs = ad

    def _pearson(n, sx, sy, sxx, syy, sxy):
        denom = (n * sxx - sx * sx) * (n * syy - sy * sy)
        if n < 2 or denom <= 0:
            return float("nan")
        return (n * sxy - sx * sy) / math.sqrt(denom)

    r = _pearson(n, sum_x, sum_y, sum_xx, sum_yy, sum_xy)
    mae = (sum_abs / n) if n else float("nan")

    summary = {
        "total_pairs": total_pairs,
        "mapped_pairs": mapped_pairs,
        "missing_map_pairs": missing_map_pairs,
        "missing_baseline_pairs": missing_baseline_pairs,
        "cross_pairs": cross_pairs,
        "within_pairs": within_pairs,
        "overlapping_pairs": n,
        "pearson_r": r,
        "mae": mae,
        "max_abs_diff": max_abs,
    }

    print("\n=== KING Comparison (Accumulator vs Baseline) ===")
    print(f"Total accumulator pairs: {summary['total_pairs']}")
    print(f"Mapped pairs: {summary['mapped_pairs']}")
    print(f"Missing map pairs: {summary['missing_map_pairs']}")
    print(f"Missing baseline pairs: {summary['missing_baseline_pairs']}")
    print(f"Cross-client pairs (mapped): {summary['cross_pairs']}")
    print(f"Within-client pairs (mapped): {summary['within_pairs']}")
    print(f"Overlapping pairs used: {summary['overlapping_pairs']}")
    print(f"Pearson r: {summary['pearson_r']:.4f}")
    print(f"MAE: {summary['mae']:.6f}")
    print(f"Max abs diff: {summary['max_abs_diff']:.6f}")

    if report_path is not None:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        with report_path.open("w") as f:
            f.write("# KING Evaluation Report\n\n")
            f.write("## Accumulator vs Baseline\n\n")
            f.write(f"- **Total accumulator pairs**: {summary['total_pairs']}\n")
            f.write(f"- **Mapped pairs**: {summary['mapped_pairs']}\n")
            f.write(f"- **Missing map pairs**: {summary['missing_map_pairs']}\n")
            f.write(f"- **Missing baseline pairs**: {summary['missing_baseline_pairs']}\n")
            f.write(f"- **Cross-client pairs (mapped)**: {summary['cross_pairs']}\n")
            f.write(f"- **Within-client pairs (mapped)**: {summary['within_pairs']}\n")
            f.write(f"- **Overlapping pairs used**: {summary['overlapping_pairs']}\n")
            f.write(f"- **Pearson r**: {summary['pearson_r']:.4f}\n")
            f.write(f"- **MAE**: {summary['mae']:.6f}\n")
            f.write(f"- **Max abs diff**: {summary['max_abs_diff']:.6f}\n")

    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--baseline-dir", required=True, type=Path)
    parser.add_argument("--center-id", type=int, default=1)
    parser.add_argument("--data-dir", type=Path, default=None, help="Optional data root containing center_* dirs")
    parser.add_argument("--report", type=Path, default=None, help="Optional report path")
    args = parser.parse_args()
    compare_from_accumulator(
        args.results_dir,
        args.baseline_dir,
        args.center_id,
        data_dir=args.data_dir,
        report_path=args.report,
    )


if __name__ == "__main__":
    main()
