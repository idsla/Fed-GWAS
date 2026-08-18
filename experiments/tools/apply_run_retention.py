#!/usr/bin/env python3
"""Apply product retention tier to a completed federated run results directory."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from pipeline.utils.retention_config import resolve_retention_settings
from pipeline.utils.run_retention import apply_retention_policy


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prune run outputs per retention tier (minimal/standard/research)"
    )
    parser.add_argument("results_dir", type=Path, help="Experiment results root (e.g. .../results_2)")
    parser.add_argument(
        "--config-path",
        type=str,
        default=None,
        help="Flower config_path (e.g. experiments/correctness/tiny_even/configs)",
    )
    parser.add_argument(
        "--experiment-config",
        type=Path,
        default=None,
        help="Direct path to experiment config.yaml",
    )
    parser.add_argument(
        "--tier",
        choices=("minimal", "standard", "research"),
        default=None,
        help="Override tier from config",
    )
    parser.add_argument("--dry-run", action="store_true", help="List actions without deleting")
    args = parser.parse_args()

    settings = resolve_retention_settings(
        config_path=args.config_path,
        experiment_config_file=str(args.experiment_config) if args.experiment_config else None,
        tier_override=args.tier,
    )

    summary = apply_retention_policy(args.results_dir, settings, dry_run=args.dry_run)
    print(json.dumps(summary.to_dict(), indent=2))
    if args.dry_run:
        print(f"\n(dry-run) Would delete {len(summary.deleted_paths)} path(s).")
    else:
        print(f"\nDeleted {len(summary.deleted_paths)} path(s); manifest: {args.results_dir}/retention_manifest.json")


if __name__ == "__main__":
    main()
