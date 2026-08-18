#!/usr/bin/env python3
"""Aggregate wall-clock time and resource metrics after a federated run.

Reads (when present):
  - stage_metrics.csv / client_metrics.csv (PerformanceMonitor output)
  - GNU /usr/bin/time -v text files (*_time.txt)
  - Client/server log directories under results/

Writes:
  - run_summary.json
  - run_summary.md (paper-ready one-liners)
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd


def _parse_gnu_time(path: Path) -> dict[str, Any]:
    text = path.read_text(errors="replace")
    out: dict[str, Any] = {"source": str(path)}
    patterns = {
        "elapsed_seconds": r"Elapsed \(wall clock\) time \(h:m:s or s\):\s*([\d:.]+)",
        "user_seconds": r"User time \(seconds\):\s*([\d.]+)",
        "sys_seconds": r"System time \(seconds\):\s*([\d.]+)",
        "max_rss_kb": r"Maximum resident set size \(kbytes\):\s*(\d+)",
        "avg_rss_kb": r"Average resident set size \(kbytes\):\s*(\d+)",
    }
    for key, pat in patterns.items():
        m = re.search(pat, text)
        if m:
            out[key] = m.group(1)
    if "max_rss_kb" in out:
        out["max_rss_mb"] = round(int(out["max_rss_kb"]) / 1024, 2)
    return out


def _read_stage_client_csvs(stage_csv: Path, client_csv: Path) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    if stage_csv.exists():
        df = pd.read_csv(stage_csv)
        summary["stage_rows"] = len(df)
        if "duration" in df.columns:
            by_stage = df.groupby("stage_name")["duration"].agg(["sum", "max", "mean"])
            summary["stage_duration_sum"] = {
                k: round(float(v), 3) for k, v in by_stage["sum"].items()
            }
        if "memory_peak_mb" in df.columns:
            summary["stage_memory_peak_mb_max"] = round(
                float(df["memory_peak_mb"].max()), 2
            )
    if client_csv.exists():
        df = pd.read_csv(client_csv)
        summary["client_metrics"] = df.to_dict(orient="records")
    return summary


def _load_csv_metrics(results_dir: Path) -> dict[str, Any]:
    """Load merged or per-node PerformanceMonitor CSV outputs."""
    summary = _read_stage_client_csvs(
        results_dir / "stage_metrics.csv",
        results_dir / "client_metrics.csv",
    )
    if summary:
        return summary

    stage_frames = []
    client_frames = []
    for logs_dir in [results_dir / "server" / "logs", *results_dir.glob("center_*/logs")]:
        if not logs_dir.is_dir():
            continue
        if (logs_dir / "stage_metrics.csv").exists():
            stage_frames.append(pd.read_csv(logs_dir / "stage_metrics.csv"))
        if (logs_dir / "client_metrics.csv").exists():
            client_frames.append(pd.read_csv(logs_dir / "client_metrics.csv"))

    if not stage_frames and not client_frames:
        return summary

    merged_dir = results_dir
    if stage_frames:
        pd.concat(stage_frames, ignore_index=True).to_csv(
            merged_dir / "stage_metrics.csv", index=False
        )
    if client_frames:
        pd.concat(client_frames, ignore_index=True).to_csv(
            merged_dir / "client_metrics.csv", index=False
        )
    return _read_stage_client_csvs(
        merged_dir / "stage_metrics.csv",
        merged_dir / "client_metrics.csv",
    )


def _dir_size_mb(path: Path) -> float:
    if not path.exists():
        return 0.0
    total = sum(f.stat().st_size for f in path.rglob("*") if f.is_file())
    return round(total / (1024 * 1024), 2)


def collect(results_dir: Path, time_files: list[Path] | None = None) -> dict[str, Any]:
    results_dir = results_dir.resolve()
    payload: dict[str, Any] = {
        "results_dir": str(results_dir),
        "collected_at": datetime.now().isoformat(),
    }

    payload["csv_metrics"] = _load_csv_metrics(results_dir)

    time_files = time_files or list(results_dir.glob("*_time.txt"))
    if time_files:
        payload["gnu_time"] = [_parse_gnu_time(p) for p in time_files]

    nodes = {}
    for name in ("server", "center_1", "center_2"):
        node_dir = results_dir / name
        if not node_dir.exists():
            continue
        nodes[name] = {
            "intermediate_mb": _dir_size_mb(node_dir / "intermediate"),
            "logs_mb": _dir_size_mb(node_dir / "logs"),
            "log_files": sorted(p.name for p in (node_dir / "logs").glob("*") if p.is_file())
            if (node_dir / "logs").exists()
            else [],
        }
    payload["nodes"] = nodes

    return payload


def to_markdown(summary: dict[str, Any], dataset_label: str) -> str:
    lines = [
        f"## Run metrics — {dataset_label}",
        "",
        f"- Results: `{summary['results_dir']}`",
    ]
    csv = summary.get("csv_metrics", {})
    if "stage_duration_sum" in csv:
        total = sum(csv["stage_duration_sum"].values())
        lines.append(f"- Total staged runtime (sum): **{total:.1f} s**")
    if "stage_memory_peak_mb_max" in csv:
        lines.append(
            f"- Peak memory (stage monitor): **{csv['stage_memory_peak_mb_max']} MB**"
        )
    for entry in summary.get("gnu_time", []):
        if "elapsed_seconds" in entry:
            lines.append(
                f"- Wall time (`{Path(entry['source']).name}`): **{entry['elapsed_seconds']}**"
            )
        if "max_rss_mb" in entry:
            lines.append(f"- Peak RSS: **{entry['max_rss_mb']} MB**")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description="Collect federated run metrics")
    parser.add_argument("results_dir", type=Path, help="Experiment results directory")
    parser.add_argument(
        "--time-file",
        action="append",
        default=[],
        help="GNU time -v output file (repeatable)",
    )
    parser.add_argument("--label", default="experiment", help="Dataset label for markdown")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Write run_summary.json/md here (default: results_dir)",
    )
    args = parser.parse_args()

    out_dir = args.output_dir or args.results_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    summary = collect(args.results_dir, [Path(p) for p in args.time_file])
    json_path = out_dir / "run_summary.json"
    md_path = out_dir / "run_summary.md"

    json_path.write_text(json.dumps(summary, indent=2))
    md_path.write_text(to_markdown(summary, args.label))
    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()
