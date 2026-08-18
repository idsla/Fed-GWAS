"""Runtime hooks wiring PerformanceMonitor into client/server apps."""

from __future__ import annotations

import logging
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Generator, Optional

logger = logging.getLogger(__name__)


def _experiment_name(settings: Dict[str, Any], fallback: str) -> str:
    return str(settings.get("experiment_name") or fallback)


class ClientPerformanceTracker:
    """Per-client performance and optional network monitoring."""

    def __init__(
        self,
        client_id: str,
        output_dir: str,
        settings: Dict[str, Any],
        experiment_name: str,
    ):
        self.client_id = client_id
        self.output_dir = Path(output_dir)
        self.settings = settings
        self._perf = None
        self._network = None
        self._active_stage: Optional[str] = None

        if settings.get("enable_performance_monitoring"):
            from pipeline.utils.performance.performance_monitoring import (
                get_client_monitor,
            )

            name = _experiment_name(settings, experiment_name)
            interval = float(settings.get("resource_monitoring_interval", 2.0))
            self._perf = get_client_monitor(
                client_id,
                str(self.output_dir),
                name,
                sample_interval=interval,
            )

        if settings.get("enable_network_monitoring"):
            from pipeline.utils.performance.network_monitor import NetworkMonitor

            interval = float(settings.get("network_monitoring_interval", 5.0))
            name = _experiment_name(settings, experiment_name)
            self._network = NetworkMonitor(
                str(self.output_dir),
                f"{name}_network",
                max_samples=int(settings.get("network_max_samples", 500)),
            )
            self._network.start_monitoring(interval=interval)

    @property
    def enabled(self) -> bool:
        return self._perf is not None or self._network is not None

    def start_stage(self, stage_name: str) -> None:
        if self._perf is None:
            return
        # Close the previous round even when the logical stage name is unchanged
        if self._active_stage:
            self.end_stage(self._active_stage, success=True)
        self._perf.start_stage(stage_name, self.client_id)
        self._active_stage = stage_name

    def end_stage(
        self,
        stage_name: str,
        *,
        success: bool = True,
        exit_reason: Optional[str] = None,
    ) -> None:
        if self._perf is None:
            return
        self._perf.end_stage(stage_name, self.client_id, success=success, exit_reason=exit_reason)
        if self._active_stage == stage_name:
            self._active_stage = None
        # Persist after each stage so Ray actor restarts do not overwrite history
        self._perf.save_metrics(append=True)

    def finalize(self) -> None:
        if self._active_stage and self._perf is not None:
            self._perf.end_stage(self._active_stage, self.client_id, success=True)
            self._active_stage = None
        if self._perf is not None:
            self._perf.stop_monitoring()
            self._perf.save_metrics(append=True)
        if self._network is not None:
            self._network.stop_monitoring()
            self._network.save_statistics()


class ServerPerformanceTracker:
    """Server-side stage participation and aggregation timing."""

    def __init__(self, output_dir: str, settings: Dict[str, Any], experiment_name: str):
        self.settings = settings
        self._monitor = None
        self._active_stage: Optional[str] = None
        self._stage_started_at: Optional[float] = None

        if not settings.get("enable_performance_monitoring"):
            return

        from pipeline.utils.performance.performance_monitoring import get_server_monitor
        import time

        self._time = time
        name = _experiment_name(settings, experiment_name)
        self._monitor = get_server_monitor(output_dir, name)

    @property
    def enabled(self) -> bool:
        return self._monitor is not None

    def begin_stage(self, stage_name: str, num_clients: int, num_failures: int = 0) -> None:
        if self._monitor is None:
            return
        # Each Flower aggregate_fit round closes the previous interval (incl. iterative_king/lr)
        if self._active_stage:
            self.end_stage(self._active_stage)
        self._monitor.record_stage_participation(stage_name, num_clients, num_failures)
        self._monitor.start_stage(stage_name, "server")
        self._active_stage = stage_name
        self._stage_started_at = self._time.time()

    def end_stage(self, stage_name: str) -> None:
        if self._monitor is None:
            return
        if self._stage_started_at is not None:
            duration = self._time.time() - self._stage_started_at
            self._monitor.record_aggregation_time(stage_name, duration)
        self._monitor.end_stage(stage_name, "server", success=True)
        if self._active_stage == stage_name:
            self._active_stage = None
            self._stage_started_at = None

    def finalize(self) -> None:
        if self._active_stage:
            self.end_stage(self._active_stage)
        if self._monitor is not None:
            self._monitor.stop_monitoring()
            self._monitor.finalize_server_metrics()
            self._monitor.save_metrics()


def merge_performance_csvs_to_results_root(results_root: Path) -> None:
    """
    Concatenate per-node stage_metrics.csv / client_metrics.csv into results_root
    for collect_run_metrics.py and metrics_collector.py.
    """
    try:
        import pandas as pd
    except ImportError:
        logger.warning("pandas not available; skipping performance CSV merge")
        return

    results_root = Path(results_root)
    if not results_root.exists():
        return

    stage_frames = []
    client_frames = []

    search_dirs = [results_root / "server" / "logs", *results_root.glob("center_*/logs")]
    for logs_dir in search_dirs:
        if not logs_dir.is_dir():
            continue
        stage_file = logs_dir / "stage_metrics.csv"
        client_file = logs_dir / "client_metrics.csv"
        if stage_file.exists():
            stage_frames.append(pd.read_csv(stage_file))
        if client_file.exists():
            client_frames.append(pd.read_csv(client_file))

    if stage_frames:
        pd.concat(stage_frames, ignore_index=True).to_csv(
            results_root / "stage_metrics.csv", index=False
        )
    if client_frames:
        pd.concat(client_frames, ignore_index=True).to_csv(
            results_root / "client_metrics.csv", index=False
        )


@contextmanager
def monitored_stage(
    tracker: Optional[ClientPerformanceTracker], stage_name: str
) -> Generator[None, None, None]:
    """Start/end performance monitoring around a single fit() stage."""
    if tracker is None or not tracker.enabled:
        yield
        return
    tracker.start_stage(stage_name)
    try:
        yield
        tracker.end_stage(stage_name, success=True)
    except Exception as exc:
        tracker.end_stage(stage_name, success=False, exit_reason=str(exc))
        raise


def resolve_results_root_from_server_output(server_output_dir: str) -> Path:
    """server output_dir is typically .../results_<tag>/server."""
    path = Path(server_output_dir).resolve()
    if path.name == "server":
        return path.parent
    return path.parent
