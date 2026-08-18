"""Tests for monitoring config resolution from experiment YAML."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline.utils.monitoring_config import resolve_monitoring_settings


def test_tiny_even_experiment_monitoring_enabled():
    repo = Path(__file__).resolve().parents[1]
    center_config = (
        repo / "experiments/correctness/tiny_even/configs/center_1/config.yaml"
    )
    settings = resolve_monitoring_settings(center_config_file=str(center_config))
    assert settings["enable_performance_monitoring"] is True
    assert settings.get("experiment_name") == "tiny_even"
