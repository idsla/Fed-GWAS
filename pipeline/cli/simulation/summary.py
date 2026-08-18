"""Summary helpers for simulation projects and prepared data."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from pipeline.cli.simulation import paths


def _format_bytes(size: int) -> str:
    """Format byte counts for human-readable CLI reports.

    Args:
        size: Raw byte count.

    Returns:
        Compact string using binary units.
    """
    value = float(size)
    for unit in ("B", "KiB", "MiB", "GiB"):
        if value < 1024 or unit == "GiB":
            return f"{value:.1f} {unit}" if unit != "B" else f"{int(value)} B"
        value /= 1024
    return f"{size} B"


def _plink_triplets(root: Path) -> list[Path]:
    """Return complete PLINK triplet prefixes below a directory.

    Args:
        root: Directory to inspect.

    Returns:
        Prefix paths whose `.bed`, `.bim`, and `.fam` files all exist.
    """
    files = [path for path in root.rglob("*") if path.is_file()] if root.exists() else []
    prefixes = {file.with_suffix("") for file in files if file.suffix in {".bed", ".bim", ".fam"}}
    return sorted(
        prefix
        for prefix in prefixes
        if all(prefix.with_suffix(suffix).exists() for suffix in (".bed", ".bim", ".fam"))
    )


def summarize_data(data_path: Path) -> dict[str, Any]:
    """Collect a lightweight summary of a data directory.

    Args:
        data_path: Directory that contains center subdirectories or PLINK files.

    Returns:
        Mapping with path, center count, file count, PLINK triplet count, and
        total byte size.
    """
    root = data_path.resolve()
    files = [path for path in root.rglob("*") if path.is_file()] if root.exists() else []
    triplets = _plink_triplets(root)
    centers = [path for path in root.iterdir() if path.is_dir() and path.name.startswith("center_")] if root.exists() else []
    center_summaries = []
    for center in sorted(centers):
        center_files = [path for path in center.rglob("*") if path.is_file()]
        center_summaries.append(
            {
                "name": center.name,
                "path": center,
                "files": len(center_files),
                "plink_triplets": len(_plink_triplets(center)),
                "bytes": sum(file.stat().st_size for file in center_files),
            }
        )
    return {
        "path": root,
        "exists": root.exists(),
        "centers": len(centers),
        "files": len(files),
        "plink_triplets": len(triplets),
        "bytes": sum(file.stat().st_size for file in files),
        "human_bytes": _format_bytes(sum(file.stat().st_size for file in files)),
        "center_summaries": center_summaries,
    }


def summarize_experiment(project_path: Path) -> dict[str, Any]:
    """Collect a lightweight summary of a simulation project.

    Args:
        project_path: Simulation project root.

    Returns:
        Mapping with project config, experiment metadata, and key directories.
    """
    root = project_path.resolve()
    project_config = paths.load_yaml(root / "fedgwas.yaml") if (root / "fedgwas.yaml").exists() else {}
    experiment_config = paths.load_yaml(root / "config.yaml") if (root / "config.yaml").exists() else {}
    data_summary = summarize_data(root / str(project_config.get("data_dir", "data")))
    config_dir = root / str(project_config.get("config_dir", "configs"))
    center_configs = sorted(config_dir.glob("center_*.yaml")) if config_dir.exists() else []
    return {
        "path": root,
        "mode": project_config.get("mode"),
        "state": project_config.get("project_state"),
        "preset": project_config.get("preset"),
        "example": project_config.get("example"),
        "experiment_name": experiment_config.get("experiment_name"),
        "experiment_category": experiment_config.get("experiment_category"),
        "scenario": experiment_config.get("scenario"),
        "num_clients": experiment_config.get("num_clients") or project_config.get("num_clients"),
        "config_dir": config_dir,
        "config_files": len(center_configs),
        "data_dir": root / str(project_config.get("data_dir", "data")),
        "data": data_summary,
        "results_dir": root / str(project_config.get("results_dir", "results")),
    }


def render_data_summary(payload: dict[str, Any], console: Console | None = None) -> None:
    """Render a data summary as a Rich terminal report.

    Args:
        payload: Mapping returned by `summarize_data`.
        console: Optional Rich console. When omitted, a default console is used.

    Returns:
        None. Writes formatted report output to the console.
    """
    output = console or Console()
    facts = Table.grid(padding=(0, 2))
    facts.add_column(style="bold")
    facts.add_column()
    facts.add_row("Path", str(payload["path"]))
    facts.add_row("Exists", "yes" if payload["exists"] else "no")
    facts.add_row("Centers", str(payload["centers"]))
    facts.add_row("Files", str(payload["files"]))
    facts.add_row("PLINK triplets", str(payload["plink_triplets"]))
    facts.add_row("Size", payload["human_bytes"])
    output.print(Panel(facts, title="Data Summary", expand=False))

    table = Table(title="Centers")
    table.add_column("Center")
    table.add_column("Files", justify="right")
    table.add_column("PLINK triplets", justify="right")
    table.add_column("Size", justify="right")
    for center in payload["center_summaries"]:
        table.add_row(
            center["name"],
            str(center["files"]),
            str(center["plink_triplets"]),
            _format_bytes(center["bytes"]),
        )
    if not payload["center_summaries"]:
        table.add_row("-", "0", "0", "0 B")
    output.print(table)


def render_experiment_summary(payload: dict[str, Any], console: Console | None = None) -> None:
    """Render an experiment summary as a Rich terminal report.

    Args:
        payload: Mapping returned by `summarize_experiment`.
        console: Optional Rich console. When omitted, a default console is used.

    Returns:
        None. Writes formatted report output to the console.
    """
    output = console or Console()
    facts = Table.grid(padding=(0, 2))
    facts.add_column(style="bold")
    facts.add_column()
    facts.add_row("Path", str(payload["path"]))
    facts.add_row("Mode", str(payload["mode"]))
    facts.add_row("State", str(payload["state"]))
    facts.add_row("Preset", str(payload["preset"]))
    facts.add_row("Example", str(payload["example"]))
    facts.add_row("Experiment", str(payload["experiment_name"]))
    facts.add_row("Category", str(payload["experiment_category"]))
    facts.add_row("Scenario", str(payload["scenario"]))
    facts.add_row("Clients", str(payload["num_clients"]))
    output.print(Panel(facts, title="Experiment Summary", expand=False))

    table = Table(title="Project Readiness")
    table.add_column("Area")
    table.add_column("Detail")
    table.add_column("Status")
    table.add_row("Configs", str(payload["config_dir"]), f"{payload['config_files']} center config(s)")
    table.add_row("Data", str(payload["data_dir"]), f"{payload['data']['plink_triplets']} PLINK triplets")
    table.add_row("Results", str(payload["results_dir"]), "exists" if payload["results_dir"].exists() else "missing")
    output.print(table)


def _status_text(ok: bool):
    """Return a styled Rich status label for a check row.

    Args:
        ok: Whether the check passed.

    Returns:
        Rich renderable containing `PASS` or `FAIL`.
    """
    from rich.text import Text

    return Text("PASS" if ok else "FAIL", style="bold green" if ok else "bold red")


def _next_steps_for_failures(sections: list[dict[str, Any]]) -> list[str]:
    """Build concise recovery hints for failed check rows.

    Args:
        sections: Check section payloads passed to `render_check_report`.

    Returns:
        Ordered command hints for the user.
    """
    failed_titles = {
        section["title"]
        for section in sections
        if any(not check["ok"] for check in section["checks"])
    }
    steps: list[str] = []
    if "Software" in failed_titles:
        steps.append("fedgwas-sim check --software")
    if "Configs" in failed_titles:
        steps.append("fedgwas-sim check --configs")
    if "Data" in failed_titles:
        steps.extend(["fedgwas-sim prepare-data", "fedgwas-sim check --data"])
    if "Outputs" in failed_titles:
        steps.append("fedgwas-sim check --outputs")
    if "Project" in failed_titles:
        steps.append("fedgwas-sim init")
    return list(dict.fromkeys(steps))


def render_check_report(
    sections: list[dict[str, Any]],
    project_dir: Path,
    selected_categories: list[str],
    failures: int,
    console: Console | None = None,
) -> None:
    """Render `fedgwas-sim check` results as a Rich terminal report.

    Args:
        sections: Ordered check sections. Each section contains `title` and a
            `checks` list with `label`, `ok`, and `detail` entries.
        project_dir: Simulation project directory or current working directory
            used as the display anchor.
        selected_categories: Check category keys requested by the CLI.
        failures: Number of failed check rows.
        console: Optional Rich console. When omitted, a default console is used.

    Returns:
        None. Writes formatted report output to the console.
    """
    output = console or Console()
    total = sum(len(section["checks"]) for section in sections)
    passed = total - failures

    facts = Table.grid(padding=(0, 2))
    facts.add_column(style="bold")
    facts.add_column()
    facts.add_row("Project", str(project_dir))
    facts.add_row("Scope", ", ".join(selected_categories))
    facts.add_row("Result", "PASS" if failures == 0 else "FAIL")
    facts.add_row("Checks", f"{passed} passed, {failures} failed")
    output.print(Panel(facts, title="Check Summary", expand=False))

    table = Table()
    table.add_column("Aspect", no_wrap=True)
    table.add_column("Status", no_wrap=True)
    table.add_column("Check")
    table.add_column("Detail")
    for section in sections:
        for check_index, check in enumerate(section["checks"]):
            table.add_row(
                section["title"] if check_index == 0 else "",
                _status_text(check["ok"]),
                str(check["label"]),
                str(check["detail"]),
            )
    output.print(table)

    if failures:
        steps = _next_steps_for_failures(sections)
        message = "One or more checks failed."
        if steps:
            message += "\n\nRun:\n" + "\n".join(f"  {step}" for step in steps)
        output.print(Panel(message, title="Next Steps", expand=False))
