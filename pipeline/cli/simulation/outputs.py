"""Result collection helpers for simulation CLI."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any

import typer


def collect_metrics(results_dir: Path, time_files: list[Path], label: str, output_dir: Path) -> None:
    """Collect lightweight run metadata into `run_summary.json` and markdown.

    Args:
        results_dir: Project results directory to inspect.
        time_files: Optional GNU time output files supplied by the user.
        label: Human-readable label for the markdown summary.
        output_dir: Destination directory for `run_summary.json` and
            `run_summary.md`.

    Returns:
        None. Writes summary files and prints their paths.
    """
    summary: dict[str, Any] = {
        "label": label,
        "results_dir": str(results_dir.resolve()),
        "collected_at": datetime.now().isoformat(),
        "time_files": [str(path) for path in time_files],
        "nodes": {},
    }
    for node_dir in [results_dir / "server", *sorted(results_dir.glob("center_*"))]:
        if not node_dir.is_dir():
            continue
        summary["nodes"][node_dir.name] = {
            "logs": sorted(path.name for path in (node_dir / "logs").glob("*") if path.is_file())
            if (node_dir / "logs").is_dir()
            else [],
            "intermediate_files": sum(1 for path in (node_dir / "intermediate").rglob("*") if path.is_file())
            if (node_dir / "intermediate").is_dir()
            else 0,
        }

    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "run_summary.json"
    md_path = output_dir / "run_summary.md"
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    md_path.write_text(
        "\n".join(
            [
                f"## Run metrics - {label}",
                "",
                f"- Results: `{summary['results_dir']}`",
                f"- Nodes: {len(summary['nodes'])}",
                f"- Time files: {len(time_files)}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    typer.echo(f"Wrote {json_path}")
    typer.echo(f"Wrote {md_path}")
