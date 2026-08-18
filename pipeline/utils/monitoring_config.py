"""Load monitoring flags from experiment and per-center YAML configs."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

import yaml

DEFAULT_MONITORING: Dict[str, Any] = {
    "enable_performance_monitoring": False,
    "enable_network_monitoring": False,
    "resource_monitoring_interval": 2.0,
    "network_monitoring_interval": 5.0,
    "network_max_samples": 500,
}


def _as_bool(value: Any, default: bool = False) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in ("1", "true", "yes", "on")
    return bool(value)


def _normalize(section: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    if not section:
        return {}
    out: Dict[str, Any] = {}
    if "enable_performance_monitoring" in section:
        out["enable_performance_monitoring"] = _as_bool(
            section["enable_performance_monitoring"]
        )
    if "enable_network_monitoring" in section:
        out["enable_network_monitoring"] = _as_bool(section["enable_network_monitoring"])
    if "resource_monitoring_interval" in section:
        out["resource_monitoring_interval"] = float(section["resource_monitoring_interval"])
    if "network_monitoring_interval" in section:
        out["network_monitoring_interval"] = float(section["network_monitoring_interval"])
    if "network_max_samples" in section:
        out["network_max_samples"] = int(section["network_max_samples"])
    return out


def resolve_experiment_config_path(config_file: Optional[str] = None) -> Optional[Path]:
    """Map a center config path or config_path directory to experiment config.yaml."""
    if not config_file:
        return None
    path = Path(config_file).resolve()
    if path.is_dir():
        candidate = path.parent / "config.yaml"
        return candidate if candidate.exists() else None
    if path.name in ("config.yaml", "config.yml") and path.parent.name == "configs":
        candidate = path.parent.parent / "config.yaml"
        return candidate if candidate.exists() else None
    if path.parent.name == "configs" and path.stem.startswith(("center_", "server")):
        candidate = path.parent.parent / "config.yaml"
        return candidate if candidate.exists() else None
    if path.parent.name.startswith("center_") and path.parent.parent.name == "configs":
        candidate = path.parent.parent.parent / "config.yaml"
        return candidate if candidate.exists() else None
    return None


def load_yaml_section(path: Path, key: str) -> Dict[str, Any]:
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f) or {}
        return data.get(key, {}) if isinstance(data, dict) else {}
    except Exception:
        return {}


def resolve_monitoring_settings(
    *,
    center_config_file: Optional[str] = None,
    config_path: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Merge monitoring settings: defaults < experiment config < center config.
    """
    settings = dict(DEFAULT_MONITORING)

    experiment_path = resolve_experiment_config_path(center_config_file)
    if experiment_path is None and config_path:
        experiment_path = resolve_experiment_config_path(config_path)

    if experiment_path is not None:
        try:
            with open(experiment_path, "r", encoding="utf-8") as f:
                experiment_data = yaml.safe_load(f) or {}
            if isinstance(experiment_data, dict):
                if experiment_data.get("experiment_name"):
                    settings["experiment_name"] = experiment_data["experiment_name"]
                settings.update(_normalize(experiment_data.get("monitoring")))
        except Exception:
            settings.update(_normalize(load_yaml_section(experiment_path, "monitoring")))

    if center_config_file:
        center_path = Path(center_config_file).resolve()
        if center_path.is_file():
            settings.update(_normalize(load_yaml_section(center_path, "monitoring")))

    return settings
