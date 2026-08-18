"""Tests for run output retention tiers."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline.utils.retention_config import resolve_retention_settings
from pipeline.utils.run_retention import apply_retention_policy


def test_minimal_tier_prunes_intermediate_and_keeps_assoc(tmp_path):
    results = tmp_path / "results"
    logs = results / "center_1" / "logs"
    interm = results / "center_1" / "intermediate"
    logs.mkdir(parents=True)
    interm.mkdir(parents=True)

    (logs / "local_lr_temp.assoc.logistic").write_text("x")
    (logs / "client_1_log.txt").write_text("log")
    (logs / "king_accumulator_state.pkl").write_bytes(b"x" * 100)
    (logs / "performance_monitor.log").write_text("mon")
    (interm / "chunk.bed").write_bytes(b"x")

    settings = resolve_retention_settings(tier_override="minimal")
    summary = apply_retention_policy(results, settings, dry_run=False)

    assert not interm.exists()
    assert (logs / "local_lr_temp.assoc.logistic").exists()
    assert not (logs / "client_1_log.txt").exists()
    assert not (logs / "king_accumulator_state.pkl").exists()
    assert len(summary.deleted_paths) >= 3


def test_standard_keeps_text_logs(tmp_path):
    results = tmp_path / "results"
    logs = results / "center_1" / "logs"
    logs.mkdir(parents=True)
    (logs / "client_1_log.txt").write_text("log")
    (logs / "stage_metrics.csv").write_text("a,b\n")

    settings = resolve_retention_settings(tier_override="standard")
    apply_retention_policy(results, settings, dry_run=False)

    assert (logs / "client_1_log.txt").exists()
    assert (logs / "stage_metrics.csv").exists()
