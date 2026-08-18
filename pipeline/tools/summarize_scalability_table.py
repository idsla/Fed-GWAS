#!/usr/bin/env python3
"""Summarize Table S1 rows from Matpool-exported stage CSVs under results/<scale>_even/."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def _load_stage_csv(path: Path) -> list[dict[str, float | str]]:
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        row["start_time"] = float(row["start_time"])
        row["end_time"] = float(row["end_time"])
        row["memory_peak_mb"] = float(row["memory_peak_mb"])
    return rows


def summarize_scale(results_root: Path, scale: str) -> dict[str, str | int | float | None]:
    run_dir = results_root / f"{scale}_even"
    server_path = run_dir / "server_stage.csv"
    if not server_path.exists():
        raise FileNotFoundError(server_path)

    server = _load_stage_csv(server_path)
    t0 = min(r["start_time"] for r in server)
    t1 = max(r["end_time"] for r in server)
    wall_min = (t1 - t0) / 60.0
    peak_server = max(r["memory_peak_mb"] for r in server)

    peaks: dict[str, float | None] = {"c1": None, "c2": None}
    for client_id in (1, 2):
        client_path = run_dir / f"client_{client_id}_stage.csv"
        if not client_path.exists():
            continue
        client_rows = _load_stage_csv(client_path)
        window = [
            r
            for r in client_rows
            if r["start_time"] >= t0 and r["end_time"] <= t1 + 1.0
        ]
        if window:
            peaks[f"c{client_id}"] = max(r["memory_peak_mb"] for r in window)

    return {
        "scale": scale.capitalize(),
        "rounds": len(server),
        "wall_min": round(wall_min, 1),
        "peak_server_mb": round(peak_server),
        "peak_c1_mb": None if peaks["c1"] is None else round(peaks["c1"]),
        "peak_c2_mb": None if peaks["c2"] is None else round(peaks["c2"]),
        "results_dir": str(run_dir),
    }


def _fmt_peak(value: float | None) -> str:
    return "---" if value is None else str(int(value))


def print_markdown_row(row: dict[str, str | int | float | None], samples: int, snps: int) -> None:
    print(
        f"| {row['scale']} | {samples:,} | {snps:,} | {row['rounds']} | "
        f"{row['wall_min']} | {row['peak_server_mb']} | "
        f"{_fmt_peak(row['peak_c1_mb'])} | {_fmt_peak(row['peak_c2_mb'])} |"
    )


def print_latex_row(row: dict[str, str | int | float | None], samples: int, snps: int) -> None:
    snp_str = f"{snps:,}".replace(",", "{,}")
    print(
        f"{row['scale']}   & {samples:,}   & {snp_str}  & {row['rounds']}  & "
        f"{row['wall_min']}  & {row['peak_server_mb']} & "
        f"{_fmt_peak(row['peak_c1_mb'])} & {_fmt_peak(row['peak_c2_mb'])} \\\\"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-root",
        type=Path,
        default=Path("results"),
        help="Directory containing <scale>_even/ folders",
    )
    parser.add_argument(
        "--scale",
        choices=("tiny", "small", "medium"),
        action="append",
        default=["tiny", "small"],
    )
    parser.add_argument("--format", choices=("md", "tex"), default="md")
    args = parser.parse_args()

    specs = {
        "tiny": (500, 5_000),
        "small": (2_000, 20_000),
        "medium": (5_000, 50_000),
    }
    for scale in args.scale:
        row = summarize_scale(args.results_root.resolve(), scale)
        samples, snps = specs[scale]
        print(f"# {row['results_dir']}")
        if args.format == "md":
            print_markdown_row(row, samples, snps)
        else:
            print_latex_row(row, samples, snps)


if __name__ == "__main__":
    main()
