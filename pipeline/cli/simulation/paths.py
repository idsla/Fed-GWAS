"""YAML and path helpers for FedGWAS simulation projects."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import typer
import yaml


def write_yaml(path: Path, payload: dict[str, Any]) -> None:
    """Write a YAML mapping while preserving insertion order for readability.

    Args:
        path: Destination YAML file path.
        payload: Mapping to serialize with `yaml.safe_dump`.

    Returns:
        None. The destination file is created or replaced.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")


def display_path(path: Path, project_dir: Path) -> str:
    """Render paths relative to the simulation project when possible.

    Args:
        path: Absolute or relative path to display.
        project_dir: Root directory of the current FedGWAS simulation project.

    Returns:
        POSIX-style project-relative path when `path` is inside `project_dir`;
        otherwise the path as a POSIX-style string.
    """
    try:
        return path.resolve().relative_to(project_dir.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def normalize_config_path(path: str | Path, project_dir: Path) -> str:
    """Normalize a user-supplied path before writing it to YAML.

    Args:
        path: User-provided path, usually a PLINK prefix from `--bfile`.
        project_dir: Root directory of the current FedGWAS simulation project.

    Returns:
        A project-relative POSIX string when possible. Absolute paths outside
        the project remain absolute.
    """
    value = Path(path)
    if value.is_absolute():
        try:
            return value.resolve().relative_to(project_dir.resolve()).as_posix()
        except ValueError:
            return value.as_posix()
    return value.as_posix()


def load_yaml(path: Path) -> dict[str, Any]:
    """Load a YAML file and require the top-level object to be a mapping.

    Args:
        path: YAML file path.

    Returns:
        Parsed top-level mapping. Empty files return an empty dict.

    Raises:
        ValueError: If the YAML top-level value is not a mapping.
        OSError: If the file cannot be opened.
    """
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"YAML file must contain a mapping: {path}")
    return data


def config_dir(project_dir: Path, project_config: dict[str, Any]) -> Path:
    """Return the configured directory that contains center/server YAML files.

    Args:
        project_dir: Simulation project root.
        project_config: Parsed `fedgwas.yaml` mapping.

    Returns:
        Absolute path to the config directory.
    """
    return project_dir / str(project_config.get("config_dir", "configs"))


def data_dir(project_dir: Path, project_config: dict[str, Any]) -> Path:
    """Return the configured project data directory.

    Args:
        project_dir: Simulation project root.
        project_config: Parsed `fedgwas.yaml` mapping.

    Returns:
        Absolute path to the data directory.
    """
    return project_dir / str(project_config.get("data_dir", "data"))


def results_dir(project_dir: Path, project_config: dict[str, Any]) -> Path:
    """Return the configured project results directory.

    Args:
        project_dir: Simulation project root.
        project_config: Parsed `fedgwas.yaml` mapping.

    Returns:
        Absolute path to the results directory.
    """
    return project_dir / str(project_config.get("results_dir", "results"))


def center_config_path(project_dir: Path, project_config: dict[str, Any], center_id: int) -> Path:
    """Resolve a center config path, preferring the new flat CLI layout.

    New projects use `configs/center_1.yaml`. Existing repository experiments
    may still use `configs/center_1/config.yaml`, so the CLI reads both layouts
    and writes the flat path when creating new files.

    Args:
        project_dir: Simulation project root.
        project_config: Parsed `fedgwas.yaml` mapping.
        center_id: One-based center id, for example `1` for center_1.

    Returns:
        Existing config path when found, otherwise the flat path where a new
        center config should be written.
    """
    base = config_dir(project_dir, project_config)
    flat = base / f"center_{center_id}.yaml"
    legacy = base / f"center_{center_id}" / "config.yaml"
    if flat.exists():
        return flat
    if legacy.exists():
        return legacy
    return flat


def server_config_path(project_dir: Path, project_config: dict[str, Any]) -> Path:
    """Resolve the server config path for flat and legacy project layouts.

    Args:
        project_dir: Simulation project root.
        project_config: Parsed `fedgwas.yaml` mapping.

    Returns:
        Existing server config path when found, otherwise the flat path where a
        new server config should be written.
    """
    base = config_dir(project_dir, project_config)
    flat = base / "server.yaml"
    legacy = base / "server" / "config.yaml"
    if flat.exists():
        return flat
    if legacy.exists():
        return legacy
    return flat


def iter_center_ids(project_config: dict[str, Any], center: int | None = None) -> list[int]:
    """Return requested one-based center ids after validating `num_clients`.

    Args:
        project_config: Parsed `fedgwas.yaml` mapping.
        center: Optional one-based center id requested by `--center`.

    Returns:
        List of one-based center ids to process.

    Raises:
        typer.BadParameter: If `num_clients` is invalid or `center` is outside
            the configured client range.
    """
    num_clients = int(project_config.get("num_clients", 0))
    if num_clients < 1:
        raise typer.BadParameter("fedgwas.yaml must set num_clients >= 1.")
    if center is not None:
        if center < 1 or center > num_clients:
            raise typer.BadParameter(f"--center must be between 1 and {num_clients}.")
        return [center]
    return list(range(1, num_clients + 1))


def plink_triplet(prefix: Path) -> list[Path]:
    """Return the three files that make up a PLINK binary dataset prefix.

    Args:
        prefix: PLINK bfile prefix without extension.

    Returns:
        Paths for `<prefix>.bed`, `<prefix>.bim`, and `<prefix>.fam`.
    """
    return [prefix.with_suffix(suffix) for suffix in (".bed", ".bim", ".fam")]


def resolve_prefix(project_dir: Path, prefix: str) -> Path:
    """Resolve a PLINK prefix from project-relative YAML into an absolute path.

    Args:
        project_dir: Simulation project root.
        prefix: Absolute or project-relative PLINK prefix string.

    Returns:
        Absolute `Path` to the PLINK prefix without extension.
    """
    path = Path(prefix)
    if path.is_absolute():
        return path
    return project_dir / path

