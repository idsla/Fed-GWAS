#!/usr/bin/env python3
"""
Performance Monitoring System for Federated GWAS Pipeline

This module provides comprehensive monitoring capabilities for:
- Timing metrics (per-stage, per-client, total runtime)
- Resource usage (memory, CPU, disk I/O)
- Communication metrics (bytes sent/received, network bandwidth)
- Client participation tracking (joins, exits, failures)
- Stage progression and completion rates
"""

import time
import json
import psutil
import logging
import threading
from datetime import datetime
from typing import Dict, List, Optional
from dataclasses import dataclass, asdict
from pathlib import Path
import numpy as np
import pandas as pd

@dataclass
class StageMetrics:
    """Metrics for a single pipeline stage"""
    stage_name: str
    client_id: str
    start_time: float
    end_time: Optional[float] = None
    duration: Optional[float] = None
    memory_peak_mb: Optional[float] = None
    cpu_percent: Optional[float] = None
    bytes_sent: Optional[int] = None
    bytes_received: Optional[int] = None
    exit_reason: Optional[str] = None
    success: bool = True

@dataclass
class ClientMetrics:
    """Aggregated metrics for a single client"""
    client_id: str
    total_runtime: float
    stages_completed: List[str]
    stages_failed: List[str]
    total_memory_peak_mb: float
    total_bytes_sent: int
    total_bytes_received: int
    exit_stage: Optional[str] = None
    exit_reason: Optional[str] = None

@dataclass
class ServerMetrics:
    """Server-side metrics"""
    total_runtime: float
    clients_connected: int
    clients_per_stage: Dict[str, int]
    failures_per_stage: Dict[str, int]
    total_bytes_processed: int
    aggregation_times: Dict[str, float]

class PerformanceMonitor:
    """Main performance monitoring class"""
    
    def __init__(
        self,
        output_dir: str,
        experiment_name: str = "fedgwas_experiment",
        *,
        background_sampling: bool = True,
        sample_interval: float = 1.0,
    ):
        self.output_dir = Path(output_dir)
        self.experiment_name = experiment_name
        # Use the output_dir directly instead of creating a subdirectory
        self.experiment_dir = self.output_dir
        self.experiment_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize logging
        self.logger = self._setup_logger()
        
        # Metrics storage
        self.stage_metrics: List[StageMetrics] = []
        self.client_metrics: Dict[str, ClientMetrics] = {}
        self.server_metrics: Optional[ServerMetrics] = None
        
        # Resource monitoring
        self.process = psutil.Process()
        self.resource_thread = None
        self.monitoring_active = False
        self.background_sampling = background_sampling
        self.sample_interval = max(0.5, float(sample_interval))
        self._rss_peak_mb = 0.0
        
        # Communication tracking (per-stage deltas via baselines)
        self.bytes_sent = 0
        self.bytes_received = 0
        self._comm_baseline_sent = 0
        self._comm_baseline_recv = 0
        
        # Thread safety
        self.lock = threading.Lock()
        # Rows already written to stage_metrics.csv (incremental append saves)
        self._saved_stage_count = 0
        
        self.logger.info(f"Performance monitor initialized for experiment: {experiment_name}")
    
    def _setup_logger(self) -> logging.Logger:
        """Set up logging for the performance monitor"""
        logger = logging.getLogger(f"performance_monitor_{self.experiment_name}")
        logger.setLevel(logging.INFO)
        
        # File handler
        log_file = self.experiment_dir / "performance_monitor.log"
        file_handler = logging.FileHandler(log_file, mode="w")
        file_handler.setFormatter(
            logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        )
        logger.addHandler(file_handler)
        
        return logger

    def _sample_rss_mb(self) -> float:
        return self.process.memory_info().rss / 1024 / 1024
    
    def start_stage(self, stage_name: str, client_id: str) -> str:
        """Start monitoring a pipeline stage"""
        with self.lock:
            start_time = time.time()
            metrics = StageMetrics(
                stage_name=stage_name,
                client_id=client_id,
                start_time=start_time
            )
            
            self._comm_baseline_sent = self.bytes_sent
            self._comm_baseline_recv = self.bytes_received
            self._rss_peak_mb = self._sample_rss_mb()

            # Start resource monitoring if not already active
            if self.background_sampling and not self.monitoring_active:
                self._start_resource_monitoring()
            
            self.stage_metrics.append(metrics)
            stage_id = f"{client_id}_{stage_name}_{start_time}"
            
            self.logger.info(f"Started stage {stage_name} for client {client_id}")
            return stage_id
    
    def end_stage(self, stage_name: str, client_id: str, success: bool = True, 
                  exit_reason: Optional[str] = None) -> None:
        """End monitoring a pipeline stage"""
        with self.lock:
            end_time = time.time()
            
            # Find the most recent stage for this client
            for metrics in reversed(self.stage_metrics):
                if (metrics.stage_name == stage_name and 
                    metrics.client_id == client_id and 
                    metrics.end_time is None):
                    
                    metrics.end_time = end_time
                    metrics.duration = end_time - metrics.start_time
                    metrics.success = success
                    metrics.exit_reason = exit_reason
                    
                    # Resource snapshot (peak during stage if background sampling ran)
                    rss_now = self._sample_rss_mb()
                    metrics.memory_peak_mb = max(self._rss_peak_mb, rss_now)
                    metrics.cpu_percent = self.process.cpu_percent(interval=None)
                    self._rss_peak_mb = 0.0
                    
                    # Per-stage communication delta
                    metrics.bytes_sent = max(0, self.bytes_sent - self._comm_baseline_sent)
                    metrics.bytes_received = max(0, self.bytes_received - self._comm_baseline_recv)
                    
                    self.logger.info(
                        f"Ended stage {stage_name} for client {client_id} "
                        f"(duration: {metrics.duration:.2f}s, success: {success})"
                    )
                    break
    
    def record_communication(self, bytes_sent: int = 0, bytes_received: int = 0) -> None:
        """Record communication metrics"""
        with self.lock:
            self.bytes_sent += bytes_sent
            self.bytes_received += bytes_received
    
    def _start_resource_monitoring(self) -> None:
        """Start background resource monitoring"""
        self.monitoring_active = True
        self.resource_thread = threading.Thread(target=self._monitor_resources, daemon=True)
        self.resource_thread.start()
        self.logger.info("Started resource monitoring")
    
    def _monitor_resources(self) -> None:
        """Background thread for resource monitoring"""
        while self.monitoring_active:
            try:
                # Monitor memory and CPU
                memory_mb = self._sample_rss_mb()
                self._rss_peak_mb = max(self._rss_peak_mb, memory_mb)
                cpu_percent = self.process.cpu_percent()
                
                # Log if resource usage is high
                if memory_mb > 1000:  # > 1GB
                    self.logger.warning(f"High memory usage: {memory_mb:.1f} MB")
                if cpu_percent > 80:  # > 80% CPU
                    self.logger.warning(f"High CPU usage: {cpu_percent:.1f}%")
                
                time.sleep(self.sample_interval)
            except Exception as e:
                self.logger.error(f"Error in resource monitoring: {e}")
                break
    
    def stop_monitoring(self) -> None:
        """Stop all monitoring"""
        self.monitoring_active = False
        if self.resource_thread:
            self.resource_thread.join(timeout=5)
        self.logger.info("Stopped performance monitoring")
    
    def generate_client_summary(self, client_id: str) -> ClientMetrics:
        """Generate summary metrics for a client"""
        client_stages = [m for m in self.stage_metrics if m.client_id == client_id]
        
        if not client_stages:
            return None
        
        total_runtime = sum(m.duration or 0 for m in client_stages)
        stages_completed = [m.stage_name for m in client_stages if m.success]
        stages_failed = [m.stage_name for m in client_stages if not m.success]
        
        total_memory_peak = max(m.memory_peak_mb or 0 for m in client_stages)
        total_bytes_sent = sum(m.bytes_sent or 0 for m in client_stages)
        total_bytes_received = sum(m.bytes_received or 0 for m in client_stages)
        
        # Find exit stage
        exit_stage = None
        exit_reason = None
        for m in reversed(client_stages):
            if not m.success:
                exit_stage = m.stage_name
                exit_reason = m.exit_reason
                break
        
        return ClientMetrics(
            client_id=client_id,
            total_runtime=total_runtime,
            stages_completed=stages_completed,
            stages_failed=stages_failed,
            total_memory_peak_mb=total_memory_peak,
            total_bytes_sent=total_bytes_sent,
            total_bytes_received=total_bytes_received,
            exit_stage=exit_stage,
            exit_reason=exit_reason
        )
    
    def _load_stage_metrics_df(self) -> pd.DataFrame:
        stage_file = self.experiment_dir / "stage_metrics.csv"
        if stage_file.exists():
            try:
                return pd.read_csv(stage_file)
            except Exception as exc:
                self.logger.warning(f"Could not read existing {stage_file}: {exc}")
        return pd.DataFrame()

    def _summaries_from_stage_df(self, stage_df: pd.DataFrame) -> Dict[str, ClientMetrics]:
        summaries: Dict[str, ClientMetrics] = {}
        if stage_df.empty or "client_id" not in stage_df.columns:
            return summaries
        for client_id in stage_df["client_id"].dropna().unique():
            rows = stage_df[stage_df["client_id"] == client_id]
            durations = rows["duration"].fillna(0) if "duration" in rows else pd.Series(dtype=float)
            mem = rows["memory_peak_mb"].fillna(0) if "memory_peak_mb" in rows else pd.Series(dtype=float)
            success_col = rows["success"] if "success" in rows else pd.Series([True] * len(rows))
            completed = rows.loc[success_col.astype(bool), "stage_name"].tolist() if "stage_name" in rows else []
            failed = rows.loc[~success_col.astype(bool), "stage_name"].tolist() if "stage_name" in rows else []
            exit_stage = None
            exit_reason = None
            if failed and "stage_name" in rows:
                fail_rows = rows.loc[~success_col.astype(bool)]
                exit_stage = fail_rows.iloc[-1]["stage_name"]
                if "exit_reason" in fail_rows.columns:
                    exit_reason = fail_rows.iloc[-1].get("exit_reason")
            summaries[str(client_id)] = ClientMetrics(
                client_id=str(client_id),
                total_runtime=float(durations.sum()),
                stages_completed=completed,
                stages_failed=failed,
                total_memory_peak_mb=float(mem.max()) if len(mem) else 0.0,
                total_bytes_sent=int(rows["bytes_sent"].fillna(0).sum()) if "bytes_sent" in rows else 0,
                total_bytes_received=int(rows["bytes_received"].fillna(0).sum())
                if "bytes_received" in rows
                else 0,
                exit_stage=exit_stage,
                exit_reason=exit_reason,
            )
        return summaries

    def save_metrics(self, *, append: bool = False) -> None:
        """Save metrics to files. Use append=True for Ray client actors that may restart."""
        with self.lock:
            pending = self.stage_metrics[self._saved_stage_count :]
            if not pending and not append:
                pending = self.stage_metrics
            if not pending and append:
                stage_df = self._load_stage_metrics_df()
                if not stage_df.empty:
                    self.client_metrics = self._summaries_from_stage_df(stage_df)
                    self._generate_summary_report(stage_df)
                return

            stage_data = [asdict(m) for m in pending]
            new_df = pd.DataFrame(stage_data) if stage_data else pd.DataFrame()
            stage_file = self.experiment_dir / "stage_metrics.csv"

            if append and stage_file.exists() and not new_df.empty:
                existing = self._load_stage_metrics_df()
                stage_df = pd.concat([existing, new_df], ignore_index=True)
            elif append and stage_file.exists() and new_df.empty:
                stage_df = self._load_stage_metrics_df()
            else:
                stage_df = new_df

            if not stage_df.empty:
                stage_df.to_csv(stage_file, index=False)

            self.client_metrics = self._summaries_from_stage_df(stage_df)

            client_data = [asdict(m) for m in self.client_metrics.values()]
            if client_data:
                client_df = pd.DataFrame(client_data)
                client_file = self.experiment_dir / "client_metrics.csv"
                client_df.to_csv(client_file, index=False)

            if self.server_metrics:
                server_file = self.experiment_dir / "server_metrics.json"
                with open(server_file, "w") as f:
                    json.dump(asdict(self.server_metrics), f, indent=2)

            self._generate_summary_report(stage_df if not stage_df.empty else None)

            self.logger.info(
                f"Metrics saved to {self.experiment_dir} "
                f"(append={append}, new_stages={len(pending)})"
            )
            self._saved_stage_count = len(self.stage_metrics)
            if append:
                self.stage_metrics = self.stage_metrics[self._saved_stage_count :]
                self._saved_stage_count = 0
            else:
                self.stage_metrics.clear()
                self._saved_stage_count = 0
            if not append:
                self.client_metrics.clear()
    
    def _generate_summary_report(self, stage_df: Optional[pd.DataFrame] = None) -> None:
        """Generate a comprehensive summary report"""
        report_file = self.experiment_dir / "performance_summary.md"
        if stage_df is None:
            stage_df = self._load_stage_metrics_df()
        use_df = stage_df is not None and not stage_df.empty

        with open(report_file, 'w') as f:
            f.write(f"# Performance Summary Report\n\n")
            f.write(f"**Experiment**: {self.experiment_name}\n")
            f.write(f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # Overall statistics
            total_clients = len(self.client_metrics)
            if use_df:
                total_stages = len(stage_df)
                successful_stages = int(stage_df["success"].fillna(True).astype(bool).sum()) if "success" in stage_df else total_stages
            else:
                total_stages = len(self.stage_metrics)
                successful_stages = sum(1 for m in self.stage_metrics if m.success)
            
            f.write("## Overall Statistics\n\n")
            f.write(f"- **Total Clients**: {total_clients}\n")
            f.write(f"- **Total Stages**: {total_stages}\n")
            f.write(f"- **Successful Stages**: {successful_stages}\n")
            rate = (successful_stages / total_stages * 100) if total_stages else 0.0
            f.write(f"- **Success Rate**: {rate:.1f}%\n\n")
            
            # Per-stage statistics
            f.write("## Per-Stage Statistics\n\n")
            stage_stats = {}
            if use_df and "stage_name" in stage_df.columns:
                for stage_name in stage_df["stage_name"].dropna().unique():
                    rows = stage_df[stage_df["stage_name"] == stage_name]
                    durations = rows["duration"].dropna().tolist() if "duration" in rows else []
                    if durations:
                        ok = rows["success"].fillna(True).astype(bool) if "success" in rows else pd.Series([True] * len(rows))
                        stage_stats[stage_name] = {
                            "count": len(rows),
                            "avg_duration": float(np.mean(durations)),
                            "min_duration": float(np.min(durations)),
                            "max_duration": float(np.max(durations)),
                            "success_rate": float(ok.sum()) / len(rows),
                        }
            else:
                for stage_name in set(m.stage_name for m in self.stage_metrics):
                    stage_metrics = [m for m in self.stage_metrics if m.stage_name == stage_name]
                    durations = [m.duration for m in stage_metrics if m.duration is not None]
                    
                    if durations:
                        stage_stats[stage_name] = {
                            'count': len(stage_metrics),
                            'avg_duration': np.mean(durations),
                            'min_duration': np.min(durations),
                            'max_duration': np.max(durations),
                            'success_rate': sum(1 for m in stage_metrics if m.success) / len(stage_metrics)
                        }
            
            for stage_name, stats in stage_stats.items():
                f.write(f"### {stage_name}\n")
                f.write(f"- **Executions**: {stats['count']}\n")
                f.write(f"- **Avg Duration**: {stats['avg_duration']:.2f}s\n")
                f.write(f"- **Duration Range**: {stats['min_duration']:.2f}s - {stats['max_duration']:.2f}s\n")
                f.write(f"- **Success Rate**: {stats['success_rate']*100:.1f}%\n\n")
            
            # Client statistics
            f.write("## Client Statistics\n\n")
            if self.client_metrics:
                total_runtimes = [m.total_runtime for m in self.client_metrics.values()]
                total_memory = [m.total_memory_peak_mb for m in self.client_metrics.values()]
                total_communication = [m.total_bytes_sent + m.total_bytes_received 
                                     for m in self.client_metrics.values()]
                
                f.write(f"- **Avg Runtime**: {np.mean(total_runtimes):.2f}s\n")
                f.write(f"- **Runtime Range**: {np.min(total_runtimes):.2f}s - {np.max(total_runtimes):.2f}s\n")
                f.write(f"- **Avg Memory Peak**: {np.mean(total_memory):.1f} MB\n")
                f.write(f"- **Total Communication**: {sum(total_communication):,} bytes\n\n")
            
            # Client details
            f.write("## Client Details\n\n")
            for client_id, metrics in self.client_metrics.items():
                f.write(f"### Client {client_id}\n")
                f.write(f"- **Total Runtime**: {metrics.total_runtime:.2f}s\n")
                f.write(f"- **Stages Completed**: {len(metrics.stages_completed)}\n")
                f.write(f"- **Stages Failed**: {len(metrics.stages_failed)}\n")
                f.write(f"- **Memory Peak**: {metrics.total_memory_peak_mb:.1f} MB\n")
                f.write(f"- **Bytes Sent**: {metrics.total_bytes_sent:,}\n")
                f.write(f"- **Bytes Received**: {metrics.total_bytes_received:,}\n")
                if metrics.exit_stage:
                    f.write(f"- **Exit Stage**: {metrics.exit_stage}\n")
                    f.write(f"- **Exit Reason**: {metrics.exit_reason}\n")
                f.write("\n")

class ServerPerformanceMonitor(PerformanceMonitor):
    """Server-specific performance monitor (no background psutil thread)."""
    
    def __init__(self, output_dir: str, experiment_name: str = "fedgwas_server"):
        super().__init__(output_dir, experiment_name, background_sampling=False)
        self.start_time = time.time()
        self.clients_per_stage = {}
        self.failures_per_stage = {}
        self.aggregation_times = {}
    
    def record_stage_participation(self, stage_name: str, num_clients: int, 
                                  num_failures: int = 0) -> None:
        """Record client participation for a stage"""
        self.clients_per_stage[stage_name] = num_clients
        self.failures_per_stage[stage_name] = num_failures
        self.logger.info(f"Stage {stage_name}: {num_clients} clients, {num_failures} failures")
    
    def record_aggregation_time(self, stage_name: str, duration: float) -> None:
        """Record server aggregation time for a stage (sums across rounds)."""
        self.aggregation_times[stage_name] = self.aggregation_times.get(stage_name, 0.0) + duration
        self.logger.info(
            f"Aggregation time for {stage_name}: {duration:.2f}s "
            f"(cumulative {self.aggregation_times[stage_name]:.2f}s)"
        )
    
    def finalize_server_metrics(self) -> None:
        """Finalize server metrics"""
        total_runtime = time.time() - self.start_time
        total_bytes_processed = self.bytes_sent + self.bytes_received
        
        self.server_metrics = ServerMetrics(
            total_runtime=total_runtime,
            clients_connected=len(self.client_metrics),
            clients_per_stage=self.clients_per_stage,
            failures_per_stage=self.failures_per_stage,
            total_bytes_processed=total_bytes_processed,
            aggregation_times=self.aggregation_times
        )

# Global monitor instances
_client_monitors = {}
_server_monitor = None

def get_client_monitor(
    client_id: str,
    output_dir: str,
    experiment_name: str,
    *,
    sample_interval: float = 1.0,
) -> PerformanceMonitor:
    """Get or create a client performance monitor"""
    global _client_monitors
    if client_id not in _client_monitors:
        _client_monitors[client_id] = PerformanceMonitor(
            output_dir,
            f"{experiment_name}_client_{client_id}",
            background_sampling=True,
            sample_interval=sample_interval,
        )
    return _client_monitors[client_id]

def get_server_monitor(output_dir: str, experiment_name: str) -> ServerPerformanceMonitor:
    """Get or create the server performance monitor"""
    global _server_monitor
    if _server_monitor is None:
        _server_monitor = ServerPerformanceMonitor(output_dir, f"{experiment_name}_server")
    return _server_monitor

def save_all_metrics() -> None:
    """Save metrics from all monitors"""
    for monitor in _client_monitors.values():
        monitor.save_metrics()
    if _server_monitor:
        _server_monitor.finalize_server_metrics()
        _server_monitor.save_metrics() 
