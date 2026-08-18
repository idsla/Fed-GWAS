"""Tests for server/client performance tracker stage closure."""

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline.utils.performance.monitoring_runtime import ServerPerformanceTracker


def test_server_begin_stage_closes_same_name_rounds():
    with tempfile.TemporaryDirectory() as tmp:
        settings = {"enable_performance_monitoring": True, "experiment_name": "test"}
        tracker = ServerPerformanceTracker(tmp, settings, "test")

        tracker.begin_stage("iterative_king", num_clients=2)
        tracker.begin_stage("iterative_king", num_clients=2)
        tracker.begin_stage("local_lr", num_clients=2)
        tracker.finalize()

        csv_path = Path(tmp) / "stage_metrics.csv"
        assert csv_path.exists()
        text = csv_path.read_text()
        lines = [ln for ln in text.strip().splitlines() if ln.startswith("iterative_king")]
        # Two iterative_king rounds + both should have duration (not empty end_time column)
        assert len(lines) == 2
        for line in lines:
            parts = line.split(",")
            assert parts[3], f"missing end_time in {line}"
            assert parts[4], f"missing duration in {line}"
