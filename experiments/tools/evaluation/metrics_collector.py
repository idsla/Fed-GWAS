#!/usr/bin/env python3
"""
Metrics Collection Framework for Federated GWAS Experiments

This module collects and aggregates metrics from experiment runs including:
- Per-stage timings from logs
- Resource usage (memory, CPU, I/O)
- Network traffic statistics
- PLINK output parsing (QC, KING, LR results)
- Export to structured formats (JSON/CSV)
"""

import os
import re
import json
import csv
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from datetime import datetime
from dataclasses import dataclass, asdict
import pandas as pd

@dataclass
class StageTiming:
    """Timing information for a pipeline stage"""
    stage_name: str
    client_id: str
    start_time: Optional[float] = None
    end_time: Optional[float] = None
    duration: Optional[float] = None
    success: bool = True

@dataclass
class ResourceUsage:
    """Resource usage metrics"""
    client_id: str
    stage_name: str
    memory_peak_mb: Optional[float] = None
    cpu_percent: Optional[float] = None
    disk_io_read_mb: Optional[float] = None
    disk_io_write_mb: Optional[float] = None

@dataclass
class NetworkTraffic:
    """Network traffic metrics"""
    client_id: str
    stage_name: str
    bytes_sent: int = 0
    bytes_received: int = 0
    packets_sent: int = 0
    packets_received: int = 0

@dataclass
class PLINKResults:
    """Parsed PLINK output results"""
    stage: str
    client_id: str
    qc_excluded_snps: Optional[List[str]] = None
    king_estimates: Optional[Dict[str, float]] = None
    lr_pvalues: Optional[Dict[str, float]] = None
    output_file: Optional[str] = None

class MetricsCollector:
    """Main metrics collection class"""
    
    def __init__(self, experiment_dir: str):
        """Initialize metrics collector"""
        self.experiment_dir = Path(experiment_dir)
        self.logger = self._setup_logger()
        
        # Storage for collected metrics
        self.stage_timings: List[StageTiming] = []
        self.resource_usage: List[ResourceUsage] = []
        self.network_traffic: List[NetworkTraffic] = []
        self.plink_results: List[PLINKResults] = []
        
        # Parsed data from existing monitoring files
        self.performance_data: Optional[pd.DataFrame] = None
        self.network_data: Optional[Dict] = None
        
        self.logger.info(f"Metrics collector initialized for: {experiment_dir}")
    
    def _setup_logger(self) -> logging.Logger:
        """Set up logging"""
        logger = logging.getLogger('metrics_collector')
        logger.setLevel(logging.INFO)
        
        handler = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        
        return logger
    
    def collect_all_metrics(self) -> Dict[str, Any]:
        """Collect all available metrics from experiment directory"""
        self.logger.info("Starting comprehensive metrics collection...")
        
        # Collect from existing monitoring files
        self._collect_from_performance_monitor()
        self._collect_from_network_monitor()
        
        # Parse client and server logs
        self._parse_client_logs()
        self._parse_server_logs()
        
        # Parse PLINK outputs
        self._parse_plink_outputs()
        
        # Aggregate metrics
        aggregated = self._aggregate_metrics()
        
        self.logger.info("Metrics collection completed")
        return aggregated
    
    def _collect_from_performance_monitor(self) -> None:
        """Collect metrics from performance monitor CSV files"""
        stage_metrics_file = self.experiment_dir / "stage_metrics.csv"
        client_metrics_file = self.experiment_dir / "client_metrics.csv"
        
        if stage_metrics_file.exists():
            try:
                self.performance_data = pd.read_csv(stage_metrics_file)
                self.logger.info(f"Loaded {len(self.performance_data)} stage metrics from CSV")
            except Exception as e:
                self.logger.warning(f"Could not load stage_metrics.csv: {e}")
        
        if client_metrics_file.exists():
            try:
                client_data = pd.read_csv(client_metrics_file)
                self.logger.info(f"Loaded {len(client_data)} client metrics from CSV")
            except Exception as e:
                self.logger.warning(f"Could not load client_metrics.csv: {e}")
    
    def _collect_from_network_monitor(self) -> None:
        """Collect metrics from network monitor JSON files"""
        network_stats_file = self.experiment_dir / "network_stats.json"
        bandwidth_file = self.experiment_dir / "bandwidth_metrics.json"
        
        if network_stats_file.exists():
            try:
                with open(network_stats_file, 'r') as f:
                    self.network_data = json.load(f)
                self.logger.info(f"Loaded network statistics")
            except Exception as e:
                self.logger.warning(f"Could not load network_stats.json: {e}")
        
        if bandwidth_file.exists():
            try:
                with open(bandwidth_file, 'r') as f:
                    bandwidth_data = json.load(f)
                self.logger.info(f"Loaded bandwidth metrics")
            except Exception as e:
                self.logger.warning(f"Could not load bandwidth_metrics.json: {e}")
    
    def _parse_client_logs(self) -> None:
        """Parse client logs for stage timings"""
        # Look for client log directories
        for center_dir in self.experiment_dir.glob("center_*/logs"):
            client_id = center_dir.parent.name
            log_file = center_dir / "client_*_log.txt"
            
            # Find actual log file
            log_files = list(center_dir.glob("*.txt"))
            if not log_files:
                continue
            
            for log_file in log_files:
                self._parse_log_file(log_file, client_id)
    
    def _parse_log_file(self, log_file: Path, client_id: str) -> None:
        """Parse a single log file for timing information"""
        try:
            with open(log_file, 'r') as f:
                content = f.read()
            
            # Look for stage start/end patterns
            # Pattern: "Stage X started" or "Starting stage X"
            stage_pattern = re.compile(
                r'(?:Stage|stage)\s+(\w+)\s+(?:started|starting|ended|completed)',
                re.IGNORECASE
            )
            
            # Extract timestamps if available
            timestamp_pattern = re.compile(r'(\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2}:\d{2})')
            
            lines = content.split('\n')
            for i, line in enumerate(lines):
                match = stage_pattern.search(line)
                if match:
                    stage_name = match.group(1)
                    # Try to extract timestamp
                    ts_match = timestamp_pattern.search(line)
                    timestamp = None
                    if ts_match:
                        try:
                            timestamp = datetime.strptime(ts_match.group(1), '%Y-%m-%d %H:%M:%S').timestamp()
                        except:
                            pass
                    
                    timing = StageTiming(
                        stage_name=stage_name,
                        client_id=client_id,
                        start_time=timestamp if 'start' in line.lower() else None,
                        end_time=timestamp if 'end' in line.lower() or 'complete' in line.lower() else None
                    )
                    self.stage_timings.append(timing)
        
        except Exception as e:
            self.logger.warning(f"Could not parse log file {log_file}: {e}")
    
    def _parse_server_logs(self) -> None:
        """Parse server logs for aggregation timing"""
        server_log_dir = self.experiment_dir / "server/logs"
        if not server_log_dir.exists():
            return
        
        # Look for server log files
        for log_file in server_log_dir.glob("*.txt"):
            self._parse_log_file(log_file, "server")
    
    def _parse_plink_outputs(self) -> None:
        """Parse PLINK output files for QC, KING, and LR results"""
        # Parse QC results (missing/imiss files)
        self._parse_qc_results()
        
        # Parse KING results (genome files)
        self._parse_king_results()
        
        # Parse LR results (assoc.logistic files)
        self._parse_lr_results()
    
    def _parse_qc_results(self) -> None:
        """Parse QC output files"""
        # Look for .imiss, .lmiss, .hwe files
        for center_dir in self.experiment_dir.glob("center_*/logs"):
            client_id = center_dir.parent.name
            
            # Find QC files
            imiss_files = list(center_dir.glob("*.imiss"))
            for imiss_file in imiss_files:
                try:
                    # Parse PLINK .imiss format
                    df = pd.read_csv(imiss_file, sep=r'\s+', skiprows=1)
                    excluded = df[df['F_MISS'] > 0.05]['IID'].tolist() if 'F_MISS' in df.columns else []
                    
                    result = PLINKResults(
                        stage="global_qc",
                        client_id=client_id,
                        qc_excluded_snps=excluded,
                        output_file=str(imiss_file)
                    )
                    self.plink_results.append(result)
                except Exception as e:
                    self.logger.warning(f"Could not parse QC file {imiss_file}: {e}")
    
    def _parse_king_results(self) -> None:
        """Parse KING kinship estimation results"""
        # Look for .genome files
        server_intermediate = self.experiment_dir / "server/intermediate"
        if server_intermediate.exists():
            genome_files = list(server_intermediate.glob("*.genome"))
            for genome_file in genome_files:
                try:
                    # Parse PLINK .genome format
                    df = pd.read_csv(genome_file, sep=r'\s+')
                    if 'KINSHIP' in df.columns:
                        # Create dictionary of kinship estimates
                        king_dict = {}
                        for _, row in df.iterrows():
                            pair_id = f"{row['IID1']}_{row['IID2']}"
                            king_dict[pair_id] = row['KINSHIP']
                        
                        result = PLINKResults(
                            stage="iterative_king",
                            client_id="server",
                            king_estimates=king_dict,
                            output_file=str(genome_file)
                        )
                        self.plink_results.append(result)
                except Exception as e:
                    self.logger.warning(f"Could not parse KING file {genome_file}: {e}")
    
    def _parse_lr_results(self) -> None:
        """Parse logistic regression results"""
        # Look for .assoc.logistic files
        for center_dir in self.experiment_dir.glob("center_*/logs"):
            client_id = center_dir.parent.name
            lr_files = list(center_dir.glob("*.assoc.logistic"))
            
            for lr_file in lr_files:
                try:
                    # Parse PLINK .assoc.logistic format
                    df = pd.read_csv(lr_file, sep=r'\s+')
                    if 'P' in df.columns and 'SNP' in df.columns:
                        pvalue_dict = {}
                        for _, row in df.iterrows():
                            pvalue_dict[row['SNP']] = row['P']
                        
                        result = PLINKResults(
                            stage="iterative_lr",
                            client_id=client_id,
                            lr_pvalues=pvalue_dict,
                            output_file=str(lr_file)
                        )
                        self.plink_results.append(result)
                except Exception as e:
                    self.logger.warning(f"Could not parse LR file {lr_file}: {e}")
        
        # Also check server intermediate directory
        server_intermediate = self.experiment_dir / "server/intermediate"
        if server_intermediate.exists():
            lr_files = list(server_intermediate.glob("*.assoc.logistic"))
            for lr_file in lr_files:
                try:
                    df = pd.read_csv(lr_file, sep=r'\s+')
                    if 'P' in df.columns and 'SNP' in df.columns:
                        pvalue_dict = {}
                        for _, row in df.iterrows():
                            pvalue_dict[row['SNP']] = row['P']
                        
                        result = PLINKResults(
                            stage="iterative_lr",
                            client_id="server",
                            lr_pvalues=pvalue_dict,
                            output_file=str(lr_file)
                        )
                        self.plink_results.append(result)
                except Exception as e:
                    self.logger.warning(f"Could not parse LR file {lr_file}: {e}")
    
    def _aggregate_metrics(self) -> Dict[str, Any]:
        """Aggregate all collected metrics into a summary"""
        aggregated = {
            'experiment_dir': str(self.experiment_dir),
            'collection_timestamp': datetime.now().isoformat(),
            'stage_timings': len(self.stage_timings),
            'resource_usage': len(self.resource_usage),
            'network_traffic': len(self.network_traffic),
            'plink_results': len(self.plink_results)
        }
        
        # Aggregate stage timings
        if self.performance_data is not None:
            stage_summary = self.performance_data.groupby('stage_name').agg({
                'duration': ['mean', 'min', 'max', 'std'],
                'success': 'sum'
            }).to_dict()
            aggregated['stage_summary'] = stage_summary
        
        # Aggregate network statistics
        if self.network_data:
            total_bytes = sum(stat.get('bytes_sent', 0) + stat.get('bytes_recv', 0) 
                            for stat in self.network_data)
            aggregated['total_network_bytes'] = total_bytes
        
        return aggregated
    
    def export_metrics(self, output_file: Optional[str] = None) -> str:
        """Export all metrics to JSON and CSV files"""
        if output_file is None:
            output_file = str(self.experiment_dir / "collected_metrics.json")
        
        output_path = Path(output_file)
        
        # Export to JSON
        export_data = {
            'stage_timings': [asdict(t) for t in self.stage_timings],
            'resource_usage': [asdict(r) for r in self.resource_usage],
            'network_traffic': [asdict(n) for n in self.network_traffic],
            'plink_results': [
                {
                    'stage': r.stage,
                    'client_id': r.client_id,
                    'output_file': r.output_file,
                    'qc_excluded_count': len(r.qc_excluded_snps) if r.qc_excluded_snps else 0,
                    'king_pairs_count': len(r.king_estimates) if r.king_estimates else 0,
                    'lr_snps_count': len(r.lr_pvalues) if r.lr_pvalues else 0
                }
                for r in self.plink_results
            ],
            'aggregated': self._aggregate_metrics()
        }
        
        with open(output_path, 'w') as f:
            json.dump(export_data, f, indent=2)
        
        # Export stage timings to CSV
        if self.stage_timings:
            csv_file = output_path.with_suffix('.csv')
            with open(csv_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=['stage_name', 'client_id', 
                                                       'start_time', 'end_time', 'duration', 'success'])
                writer.writeheader()
                for timing in self.stage_timings:
                    writer.writerow(asdict(timing))
        
        self.logger.info(f"Metrics exported to {output_path}")
        return str(output_path)


def main():
    """Command-line interface for metrics collection"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Collect metrics from federated GWAS experiments")
    parser.add_argument("experiment_dir", type=str, help="Path to experiment results directory")
    parser.add_argument("--output", type=str, help="Output file path (default: experiment_dir/collected_metrics.json)")
    
    args = parser.parse_args()
    
    collector = MetricsCollector(args.experiment_dir)
    metrics = collector.collect_all_metrics()
    
    output_file = collector.export_metrics(args.output)
    
    print(f"Metrics collection completed. Results saved to: {output_file}")
    print(f"\nSummary:")
    print(f"  - Stage timings: {metrics['stage_timings']}")
    print(f"  - Resource usage records: {metrics['resource_usage']}")
    print(f"  - Network traffic records: {metrics['network_traffic']}")
    print(f"  - PLINK results: {metrics['plink_results']}")


if __name__ == "__main__":
    main()


