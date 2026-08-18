"""Environment, project, and data validation helpers for simulation CLI."""

from __future__ import annotations

import shutil
import sys
from importlib import metadata
from pathlib import Path
from typing import Any

import typer

from pipeline.cli.simulation import paths


def load_project_config(project_dir: Path | None = None) -> tuple[Path, dict[str, Any]]:
    """Load and validate `fedgwas.yaml` from a simulation project.

    Args:
        project_dir: Optional project root. When omitted, uses `Path.cwd()`.

    Returns:
        Tuple of `(project_root, project_config)`.

    Raises:
        typer.BadParameter: If `fedgwas.yaml` is missing or does not declare
            `mode: simulation`.
    """
    root = (project_dir or Path.cwd()).resolve()
    config_path = root / "fedgwas.yaml"
    if not config_path.exists():
        raise typer.BadParameter(
            f"fedgwas.yaml not found in {root}. Run 'fedgwas-sim setup-experiment' first."
        )
    config = paths.load_yaml(config_path)
    if config.get("mode") != "simulation":
        raise typer.BadParameter(f"{config_path} is not a simulation project.")
    return root, config


def resolve_plink(project_dir: Path, project_config: dict[str, Any]) -> str | None:
    """Resolve PLINK from `fedgwas.yaml`, PATH, or bundled source-tree binaries.

    Args:
        project_dir: Simulation project root.
        project_config: Parsed `fedgwas.yaml` mapping.

    Returns:
        Executable path/name for PLINK or PLINK2, or `None` when unavailable.
    """
    configured = str(project_config.get("plink", "auto"))
    if configured and configured != "auto":
        candidate = Path(configured)
        if not candidate.is_absolute():
            candidate = project_dir / candidate
        return str(candidate) if candidate.exists() else None
    for name in ("plink", "plink2"):
        found = shutil.which(name)
        if found:
            return found
    bundled_root = Path(__file__).resolve().parents[3] / "plink"
    candidates = [
        bundled_root / "plink_win" / "plink.exe",
        bundled_root / "plink_linux" / "plink",
        bundled_root / "plink_mac" / "plink",
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    return None


def flwr_version_ok() -> tuple[bool, str]:
    """Check that the installed Flower version is the pinned 1.19 line.

    Returns:
        Tuple `(ok, detail)`, where `ok` is true only for Flower 1.19.x and
        `detail` is the installed version or a missing-package message.
    """
    try:
        version = metadata.version("flwr")
    except metadata.PackageNotFoundError:
        return False, "Flower is not installed"
    parts: list[int] = []
    for token in version.split(".")[:2]:
        try:
            parts.append(int("".join(ch for ch in token if ch.isdigit())))
        except ValueError:
            parts.append(0)
    ok = parts[:2] == [1, 19]
    return ok, version


def python_version_ok() -> tuple[bool, str]:
    """Return whether the running interpreter satisfies FedGWAS requirements.

    Returns:
        Tuple `(ok, detail)` for Python >=3.11 and the current version string.
    """
    return sys.version_info >= (3, 11), sys.version.split()[0]


def is_writable_dir(path: Path) -> bool:
    """Probe whether a directory can be created and written to.

    Args:
        path: Directory path to test.

    Returns:
        True if a temporary probe file can be written and deleted; false
        otherwise.
    """
    try:
        path.mkdir(parents=True, exist_ok=True)
        probe = path / ".fedgwas-write-test"
        probe.write_text("ok", encoding="utf-8")
        probe.unlink()
        return True
    except OSError:
        return False


def verify_data(project_dir: Path, project_config: dict[str, Any], center: int | None) -> list[Path]:
    """Return missing config or PLINK triplet paths for all requested centers.

    Args:
        project_dir: Simulation project root.
        project_config: Parsed `fedgwas.yaml` mapping.
        center: Optional one-based center id requested by `--center`.

    Returns:
        List of missing paths. An empty list means verification passed.
    """
    missing: list[Path] = []
    for center_id in paths.iter_center_ids(project_config, center):
        config_path = paths.center_config_path(project_dir, project_config, center_id)
        if not config_path.exists():
            missing.append(config_path)
            continue
        center_config = paths.load_yaml(config_path)
        input_data = center_config.get("input_data") or {}
        prefix_value = input_data.get("path")
        if not prefix_value:
            missing.append(config_path)
            continue
        prefix = paths.resolve_prefix(project_dir, str(prefix_value))
        missing.extend(path for path in paths.plink_triplet(prefix) if not path.exists())
    return missing


