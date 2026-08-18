"""Retention policy for experiment run outputs (product tiers)."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

import yaml

from pipeline.utils.monitoring_config import resolve_experiment_config_path

# Tier presets aligned with product docs: minimal (A), standard (B), research (C)
TIER_PRESETS: Dict[str, Dict[str, bool]] = {
    "minimal": {
        "keep_intermediate": False,
        "keep_network_monitor": False,
        "keep_performance_logs": False,
        "keep_crypto_state": False,
        "keep_king_accumulator": False,
        "keep_text_logs": False,
        "keep_node_performance_csv": False,
        "keep_merged_performance_csv": True,
        "keep_plink_science_outputs": True,
        "keep_evaluation_reports": True,
    },
    "standard": {
        "keep_intermediate": False,
        "keep_network_monitor": False,
        "keep_performance_logs": False,
        "keep_crypto_state": False,
        "keep_king_accumulator": False,
        "keep_text_logs": True,
        "keep_node_performance_csv": True,
        "keep_merged_performance_csv": True,
        "keep_plink_science_outputs": True,
        "keep_evaluation_reports": True,
    },
    "research": {
        "keep_intermediate": True,
        "keep_network_monitor": True,
        "keep_performance_logs": True,
        "keep_crypto_state": True,
        "keep_king_accumulator": True,
        "keep_text_logs": True,
        "keep_node_performance_csv": True,
        "keep_merged_performance_csv": True,
        "keep_plink_science_outputs": True,
        "keep_evaluation_reports": True,
    },
}

DEFAULT_RETENTION: Dict[str, Any] = {
    "tier": "standard",
    "auto_apply_on_complete": True,
}


def _as_bool(value: Any, default: bool = False) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in ("1", "true", "yes", "on")
    return bool(value)


def _normalize_retention(section: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    if not section:
        return {}
    out: Dict[str, Any] = {}
    if "tier" in section and section["tier"]:
        out["tier"] = str(section["tier"]).strip().lower()
    if "auto_apply_on_complete" in section:
        out["auto_apply_on_complete"] = _as_bool(section["auto_apply_on_complete"])
    for key in TIER_PRESETS["standard"]:
        if key in section:
            out[key] = _as_bool(section[key])
    return out


def resolve_retention_settings(
    *,
    config_path: Optional[str] = None,
    experiment_config_file: Optional[str] = None,
    tier_override: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Merge retention: defaults < tier preset < experiment yaml < overrides.
    """
    settings: Dict[str, Any] = dict(DEFAULT_RETENTION)

    exp_path = None
    if experiment_config_file:
        exp_path = Path(experiment_config_file).resolve()
    elif config_path:
        exp_path = resolve_experiment_config_path(config_path)

    yaml_section: Dict[str, Any] = {}
    if exp_path and exp_path.is_file():
        try:
            with open(exp_path, "r", encoding="utf-8") as f:
                data = yaml.safe_load(f) or {}
            if isinstance(data, dict):
                yaml_section = data.get("retention") or {}
        except Exception:
            yaml_section = {}

    settings.update(_normalize_retention(yaml_section))

    tier = (tier_override or settings.get("tier") or "standard").lower()
    if tier not in TIER_PRESETS:
        tier = "standard"
    settings["tier"] = tier

    merged = dict(TIER_PRESETS[tier])
    for key, value in settings.items():
        if key.startswith("keep_"):
            merged[key] = bool(value)
    merged["tier"] = tier
    merged["auto_apply_on_complete"] = settings.get(
        "auto_apply_on_complete", DEFAULT_RETENTION["auto_apply_on_complete"]
    )
    return merged
