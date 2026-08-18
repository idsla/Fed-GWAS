"""Client stage metrics survive incremental append saves (Ray actor restarts)."""

import tempfile
from pathlib import Path

import pandas as pd

from pipeline.utils.performance.performance_monitoring import PerformanceMonitor


def test_save_metrics_append_across_sessions():
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp)

        mon1 = PerformanceMonitor(str(out), "test_client", background_sampling=False)
        mon1.start_stage("init_chunks", "client_1")
        mon1.end_stage("init_chunks", "client_1", success=True)
        mon1.save_metrics(append=True)

        mon2 = PerformanceMonitor(str(out), "test_client", background_sampling=False)
        mon2.start_stage("iterative_king", "client_1")
        mon2.end_stage("iterative_king", "client_1", success=True)
        mon2.save_metrics(append=True)

        df = pd.read_csv(out / "stage_metrics.csv")
        assert len(df) == 2
        assert set(df["stage_name"]) == {"init_chunks", "iterative_king"}

        client_df = pd.read_csv(out / "client_metrics.csv")
        assert len(client_df) == 1
        assert client_df.iloc[0]["client_id"] == "client_1"
        assert client_df.iloc[0]["total_memory_peak_mb"] >= 0
