#!/usr/bin/env python3
"""
Network Monitoring Utility for Federated GWAS Pipeline

This module provides network monitoring capabilities including:
- Bandwidth usage tracking
- Network interface statistics
- Connection monitoring
- Network performance metrics
"""

import time
import json
import logging
import threading
from datetime import datetime
from typing import Dict, List
from dataclasses import dataclass, asdict
from pathlib import Path
import psutil

@dataclass
class NetworkStats:
    """Network statistics for a time period"""
    timestamp: float
    interface: str
    bytes_sent: int
    bytes_recv: int
    packets_sent: int
    packets_recv: int
    errin: int
    errout: int
    dropin: int
    dropout: int

@dataclass
class BandwidthMetrics:
    """Bandwidth metrics"""
    interface: str
    upload_speed_mbps: float
    download_speed_mbps: float
    total_bytes_sent: int
    total_bytes_recv: int
    connection_count: int

class NetworkMonitor:
    """Network monitoring class"""
    
    def __init__(
        self,
        output_dir: str,
        experiment_name: str = "network_monitor",
        *,
        max_samples: int = 500,
    ):
        self.output_dir = Path(output_dir)
        self.experiment_name = experiment_name
        self.monitor_dir = self.output_dir / experiment_name
        self.monitor_dir.mkdir(parents=True, exist_ok=True)
        self.max_samples = max(50, int(max_samples))
        
        # Initialize logging
        self.logger = self._setup_logger()
        
        # Monitoring state
        self.monitoring_active = False
        self.monitor_thread = None
        
        # Statistics storage (ring buffer cap to limit memory / log churn)
        self.network_stats: List[NetworkStats] = []
        self.baseline_stats = {}
        
        # Thread safety
        self.lock = threading.Lock()
        
        self.logger.info(f"Network monitor initialized for experiment: {experiment_name}")
    
    def _setup_logger(self) -> logging.Logger:
        """Set up logging for the network monitor"""
        logger = logging.getLogger(f"network_monitor_{self.experiment_name}")
        logger.setLevel(logging.INFO)
        
        # File handler
        log_file = self.monitor_dir / "network_monitor.log"
        file_handler = logging.FileHandler(log_file, mode="w")
        file_handler.setFormatter(
            logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        )
        logger.addHandler(file_handler)
        
        return logger
    
    def start_monitoring(self, interval: float = 1.0) -> None:
        """Start network monitoring"""
        if self.monitoring_active:
            self.logger.warning("Network monitoring already active")
            return
        
        self.monitoring_active = True
        self.monitor_thread = threading.Thread(
            target=self._monitor_network, 
            args=(interval,), 
            daemon=True
        )
        self.monitor_thread.start()
        self.logger.info(f"Started network monitoring with {interval}s interval")
    
    def stop_monitoring(self) -> None:
        """Stop network monitoring"""
        self.monitoring_active = False
        if self.monitor_thread:
            self.monitor_thread.join(timeout=5)
        self.logger.info("Stopped network monitoring")
    
    def _monitor_network(self, interval: float) -> None:
        """Background thread for network monitoring"""
        # Get baseline statistics
        self._capture_baseline()
        
        while self.monitoring_active:
            try:
                self._capture_network_stats()
                time.sleep(interval)
            except Exception as e:
                self.logger.error(f"Error in network monitoring: {e}")
                break
    
    def _capture_baseline(self) -> None:
        """Capture baseline network statistics"""
        try:
            net_io = psutil.net_io_counters(pernic=True)
            self.baseline_stats = {
                interface: {
                    'bytes_sent': stats.bytes_sent,
                    'bytes_recv': stats.bytes_recv,
                    'packets_sent': stats.packets_sent,
                    'packets_recv': stats.packets_recv,
                    'errin': stats.errin,
                    'errout': stats.errout,
                    'dropin': stats.dropin,
                    'dropout': stats.dropout
                }
                for interface, stats in net_io.items()
            }
            self.logger.info(f"Captured baseline stats for {len(self.baseline_stats)} interfaces")
        except Exception as e:
            self.logger.error(f"Error capturing baseline stats: {e}")
    
    def _capture_network_stats(self) -> None:
        """Capture current network statistics"""
        try:
            net_io = psutil.net_io_counters(pernic=True)
            timestamp = time.time()
            
            for interface, stats in net_io.items():
                # Skip loopback and virtual interfaces for most monitoring
                if interface.startswith(('lo', 'docker', 'veth', 'br-')):
                    continue
                
                network_stat = NetworkStats(
                    timestamp=timestamp,
                    interface=interface,
                    bytes_sent=stats.bytes_sent,
                    bytes_recv=stats.bytes_recv,
                    packets_sent=stats.packets_sent,
                    packets_recv=stats.packets_recv,
                    errin=stats.errin,
                    errout=stats.errout,
                    dropin=stats.dropin,
                    dropout=stats.dropout
                )
                
                with self.lock:
                    self.network_stats.append(network_stat)
                    overflow = len(self.network_stats) - self.max_samples
                    if overflow > 0:
                        del self.network_stats[:overflow]
                
                # Log significant network activity at debug to avoid huge log files
                if stats.bytes_sent > 1000000 or stats.bytes_recv > 1000000:  # > 1MB
                    self.logger.debug(
                        f"High network activity on {interface}: "
                        f"sent={stats.bytes_sent/1024/1024:.1f}MB, "
                        f"recv={stats.bytes_recv/1024/1024:.1f}MB"
                    )
        
        except Exception as e:
            self.logger.error(f"Error capturing network stats: {e}")
    
    def get_bandwidth_metrics(self, interface: str = None) -> List[BandwidthMetrics]:
        """Calculate bandwidth metrics for specified interface(s)"""
        with self.lock:
            if not self.network_stats:
                return []
            
            # Group stats by interface
            interface_stats = {}
            for stat in self.network_stats:
                if interface and stat.interface != interface:
                    continue
                
                if stat.interface not in interface_stats:
                    interface_stats[stat.interface] = []
                interface_stats[stat.interface].append(stat)
            
            metrics = []
            for iface, stats in interface_stats.items():
                if len(stats) < 2:
                    continue
                
                # Calculate bandwidth over time
                time_span = stats[-1].timestamp - stats[0].timestamp
                if time_span <= 0:
                    continue
                
                bytes_sent_diff = stats[-1].bytes_sent - stats[0].bytes_sent
                bytes_recv_diff = stats[-1].bytes_recv - stats[0].bytes_recv
                
                upload_speed = (bytes_sent_diff * 8) / (time_span * 1000000)  # Mbps
                download_speed = (bytes_recv_diff * 8) / (time_span * 1000000)  # Mbps
                
                # Get connection count
                try:
                    connections = psutil.net_connections()
                    connection_count = len([c for c in connections if c.status == 'ESTABLISHED'])
                except:
                    connection_count = 0
                
                metric = BandwidthMetrics(
                    interface=iface,
                    upload_speed_mbps=upload_speed,
                    download_speed_mbps=download_speed,
                    total_bytes_sent=stats[-1].bytes_sent,
                    total_bytes_recv=stats[-1].bytes_recv,
                    connection_count=connection_count
                )
                metrics.append(metric)
            
            return metrics
    
    def save_network_stats(self) -> None:
        """Save network statistics to files"""
        with self.lock:
            if not self.network_stats:
                self.logger.warning("No network statistics to save")
                return
            
            # Save detailed network stats
            stats_data = [asdict(stat) for stat in self.network_stats]
            stats_file = self.monitor_dir / "network_stats.json"
            with open(stats_file, 'w') as f:
                json.dump(stats_data, f, indent=2)
            
            # Calculate and save bandwidth metrics
            bandwidth_metrics = self.get_bandwidth_metrics()
            bandwidth_data = [asdict(metric) for metric in bandwidth_metrics]
            bandwidth_file = self.monitor_dir / "bandwidth_metrics.json"
            with open(bandwidth_file, 'w') as f:
                json.dump(bandwidth_data, f, indent=2)
            
            # Generate network summary report
            self._generate_network_report()
            
            self.logger.info(f"Network statistics saved to {self.monitor_dir}")
            self.network_stats.clear()
    
    def _generate_network_report(self) -> None:
        """Generate a network performance summary report"""
        report_file = self.monitor_dir / "network_summary.md"
        
        with open(report_file, 'w') as f:
            f.write(f"# Network Performance Summary Report\n\n")
            f.write(f"**Experiment**: {self.experiment_name}\n")
            f.write(f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # Overall statistics
            total_stats = len(self.network_stats)
            interfaces = set(stat.interface for stat in self.network_stats)
            
            f.write("## Overall Statistics\n\n")
            f.write(f"- **Total Measurements**: {total_stats}\n")
            f.write(f"- **Interfaces Monitored**: {len(interfaces)}\n")
            f.write(f"- **Monitoring Duration**: {self._get_monitoring_duration():.1f}s\n\n")
            
            # Per-interface statistics
            f.write("## Per-Interface Statistics\n\n")
            for interface in sorted(interfaces):
                interface_stats = [s for s in self.network_stats if s.interface == interface]
                if not interface_stats:
                    continue
                
                total_sent = sum(s.bytes_sent for s in interface_stats)
                total_recv = sum(s.bytes_recv for s in interface_stats)
                total_errors = sum(s.errin + s.errout for s in interface_stats)
                total_drops = sum(s.dropin + s.dropout for s in interface_stats)
                
                f.write(f"### {interface}\n")
                f.write(f"- **Total Bytes Sent**: {total_sent:,}\n")
                f.write(f"- **Total Bytes Received**: {total_recv:,}\n")
                f.write(f"- **Total Errors**: {total_errors}\n")
                f.write(f"- **Total Drops**: {total_drops}\n")
                
                # Calculate average bandwidth
                if len(interface_stats) > 1:
                    time_span = interface_stats[-1].timestamp - interface_stats[0].timestamp
                    if time_span > 0:
                        avg_upload = (total_sent * 8) / (time_span * 1000000)  # Mbps
                        avg_download = (total_recv * 8) / (time_span * 1000000)  # Mbps
                        f.write(f"- **Average Upload Speed**: {avg_upload:.2f} Mbps\n")
                        f.write(f"- **Average Download Speed**: {avg_download:.2f} Mbps\n")
                
                f.write("\n")
            
            # Current bandwidth metrics
            f.write("## Current Bandwidth Metrics\n\n")
            bandwidth_metrics = self.get_bandwidth_metrics()
            for metric in bandwidth_metrics:
                f.write(f"### {metric.interface}\n")
                f.write(f"- **Upload Speed**: {metric.upload_speed_mbps:.2f} Mbps\n")
                f.write(f"- **Download Speed**: {metric.download_speed_mbps:.2f} Mbps\n")
                f.write(f"- **Active Connections**: {metric.connection_count}\n")
                f.write(f"- **Total Data Sent**: {metric.total_bytes_sent:,} bytes\n")
                f.write(f"- **Total Data Received**: {metric.total_bytes_recv:,} bytes\n\n")
    
    def _get_monitoring_duration(self) -> float:
        """Get the total monitoring duration in seconds"""
        if not self.network_stats:
            return 0.0
        return self.network_stats[-1].timestamp - self.network_stats[0].timestamp

def run_bandwidth_test(duration: int = 30, interface: str = None) -> Dict:
    """Run a quick bandwidth test"""
    monitor = NetworkMonitor("/tmp", "bandwidth_test")
    monitor.start_monitoring(interval=0.5)
    
    print(f"Running bandwidth test for {duration} seconds...")
    time.sleep(duration)
    
    monitor.stop_monitoring()
    metrics = monitor.get_bandwidth_metrics(interface)
    monitor.save_network_stats()
    
    return {
        'metrics': [asdict(m) for m in metrics],
        'duration': duration
    }

if __name__ == "__main__":
    # Example usage
    import argparse
    
    parser = argparse.ArgumentParser(description="Network monitoring utility")
    parser.add_argument("--duration", type=int, default=60, help="Monitoring duration in seconds")
    parser.add_argument("--output", type=str, default="/tmp/network_monitor", help="Output directory")
    parser.add_argument("--experiment", type=str, default="test", help="Experiment name")
    
    args = parser.parse_args()
    
    monitor = NetworkMonitor(args.output, args.experiment)
    monitor.start_monitoring()
    
    print(f"Monitoring network for {args.duration} seconds...")
    time.sleep(args.duration)
    
    monitor.stop_monitoring()
    monitor.save_network_stats()
    
    print("Network monitoring completed. Check the output directory for results.") 
