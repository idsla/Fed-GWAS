#!/usr/bin/env python3
"""
SF-GWAS Comparison Tool

This module compares federated GWAS pipeline performance with SF-GWAS benchmarks:
- Load SF-GWAS published benchmarks (QC ~4.5h, PCA ~44h, LR ~77.8h for UKB-scale)
- Compare runtime per stage
- Generate comparison visualizations
- Highlight advantages/disadvantages
"""

import os
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

# SF-GWAS published benchmarks (from experiments.md)
SF_GWAS_BENCHMARKS = {
    'ukb_scale': {
        'samples': 275812,
        'snps': 98000000,
        'clients': 7,
        'qc_time_hours': 4.5,
        'pca_time_hours': 44.0,  # Note: our pipeline doesn't do PCA, but KING is similar
        'lr_time_hours': 77.8,
        'total_time_hours': 126.3,  # ~5.3 days
        'notes': 'UK Biobank-scale imputed genotypes'
    },
    'small_scale': {
        'samples': 5000,
        'snps': 50000,
        'clients': 3,
        'qc_time_hours': 0.1,  # Estimated
        'pca_time_hours': 0.5,  # Estimated
        'lr_time_hours': 0.5,  # Estimated
        'total_time_hours': 1.1,
        'notes': 'Estimated for small scale'
    }
}

class SFGWASComparator:
    """SF-GWAS comparison class"""
    
    def __init__(self, federated_results_dir: str, scale: str = 'ukb_scale'):
        """Initialize SF-GWAS comparator"""
        self.federated_dir = Path(federated_results_dir)
        self.scale = scale
        self.benchmark = SF_GWAS_BENCHMARKS.get(scale, SF_GWAS_BENCHMARKS['ukb_scale'])
        
        self.logger = self._setup_logger()
        
        # Loaded metrics
        self.federated_metrics: Optional[Dict[str, Any]] = None
        self.comparison_results: Dict[str, Any] = {}
        
        self.logger.info(f"SF-GWAS comparator initialized for scale: {scale}")
    
    def _setup_logger(self) -> logging.Logger:
        """Set up logging"""
        logger = logging.getLogger('sf_gwas_comparison')
        logger.setLevel(logging.INFO)
        
        handler = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        
        return logger
    
    def load_federated_metrics(self) -> None:
        """Load federated experiment metrics"""
        self.logger.info("Loading federated metrics...")
        
        # Try to load from collected_metrics.json
        metrics_file = self.federated_dir / "collected_metrics.json"
        if metrics_file.exists():
            try:
                with open(metrics_file, 'r') as f:
                    data = json.load(f)
                    self.federated_metrics = data.get('aggregated', {})
            except Exception as e:
                self.logger.warning(f"Could not load collected_metrics.json: {e}")
        
        # Try to load from stage_metrics.csv
        stage_metrics_file = self.federated_dir / "stage_metrics.csv"
        if stage_metrics_file.exists():
            try:
                df = pd.read_csv(stage_metrics_file)
                # Calculate per-stage timings
                stage_times = {}
                for stage in df['stage_name'].unique():
                    stage_data = df[df['stage_name'] == stage]
                    stage_times[stage] = {
                        'mean_duration': stage_data['duration'].mean(),
                        'total_duration': stage_data['duration'].sum(),
                        'max_duration': stage_data['duration'].max()
                    }
                
                if not self.federated_metrics:
                    self.federated_metrics = {}
                self.federated_metrics['stage_times'] = stage_times
            except Exception as e:
                self.logger.warning(f"Could not load stage_metrics.csv: {e}")
        
        # Map our stages to SF-GWAS stages
        self._map_stages_to_sf_gwas()
        
        self.logger.info("Federated metrics loaded")
    
    def _map_stages_to_sf_gwas(self) -> None:
        """Map our pipeline stages to SF-GWAS stages"""
        if not self.federated_metrics or 'stage_times' not in self.federated_metrics:
            return
        
        stage_times = self.federated_metrics['stage_times']
        
        # Map our stages to SF-GWAS equivalent stages
        # QC stages -> QC
        qc_stages = ['local_qc', 'global_qc', 'global_qc_response']
        qc_time = sum(stage_times.get(stage, {}).get('total_duration', 0) 
                     for stage in qc_stages) / 3600  # Convert to hours
        
        # KING stages -> Similar to PCA (kinship estimation)
        king_stages = ['init_chunks', 'iterative_king']
        king_time = sum(stage_times.get(stage, {}).get('total_duration', 0) 
                       for stage in king_stages) / 3600
        
        # LR stages -> LR
        lr_stages = ['local_lr', 'local_lr_filter_response', 'init_chunks_lr', 'iterative_lr']
        lr_time = sum(stage_times.get(stage, {}).get('total_duration', 0) 
                     for stage in lr_stages) / 3600
        
        # Total time
        all_stages = list(stage_times.keys())
        total_time = sum(stage_times.get(stage, {}).get('total_duration', 0) 
                        for stage in all_stages) / 3600
        
        self.federated_metrics['mapped_times'] = {
            'qc_time_hours': qc_time,
            'king_time_hours': king_time,  # Maps to PCA in SF-GWAS
            'lr_time_hours': lr_time,
            'total_time_hours': total_time
        }
    
    def compare_performance(self) -> Dict[str, Any]:
        """Compare federated performance with SF-GWAS benchmarks"""
        if not self.federated_metrics or 'mapped_times' not in self.federated_metrics:
            self.logger.warning("Cannot compare: federated metrics not loaded")
            return {}
        
        fed_times = self.federated_metrics['mapped_times']
        
        comparison = {
            'scale': self.scale,
            'benchmark': self.benchmark,
            'federated': fed_times,
            'speedup': {},
            'slowdown': {},
            'relative_performance': {}
        }
        
        # Compare each stage
        stages = ['qc', 'king', 'lr', 'total']
        for stage in stages:
            sf_gwas_key = f'{stage}_time_hours'
            fed_key = f'{stage}_time_hours'
            
            if sf_gwas_key in self.benchmark and fed_key in fed_times:
                sf_gwas_time = self.benchmark[sf_gwas_key]
                fed_time = fed_times[fed_key]
                
                if fed_time > 0:
                    speedup = sf_gwas_time / fed_time
                    comparison['speedup'][stage] = speedup
                    comparison['relative_performance'][stage] = {
                        'sf_gwas_hours': sf_gwas_time,
                        'federated_hours': fed_time,
                        'difference_hours': fed_time - sf_gwas_time,
                        'percent_difference': ((fed_time - sf_gwas_time) / sf_gwas_time) * 100
                    }
        
        self.comparison_results = comparison
        return comparison
    
    def generate_comparison_plot(self, output_file: Optional[str] = None) -> str:
        """Generate comparison bar chart"""
        if not self.comparison_results:
            self.logger.warning("Cannot generate plot: comparison not run")
            return ""
        
        stages = ['QC', 'KING/PCA', 'LR', 'Total']
        sf_gwas_times = [
            self.benchmark.get('qc_time_hours', 0),
            self.benchmark.get('pca_time_hours', 0),  # Use PCA for KING comparison
            self.benchmark.get('lr_time_hours', 0),
            self.benchmark.get('total_time_hours', 0)
        ]
        
        fed_times = [
            self.federated_metrics['mapped_times'].get('qc_time_hours', 0),
            self.federated_metrics['mapped_times'].get('king_time_hours', 0),
            self.federated_metrics['mapped_times'].get('lr_time_hours', 0),
            self.federated_metrics['mapped_times'].get('total_time_hours', 0)
        ]
        
        x = np.arange(len(stages))
        width = 0.35
        
        fig, ax = plt.subplots(figsize=(10, 6))
        bars1 = ax.bar(x - width/2, sf_gwas_times, width, label='SF-GWAS', color='#1f77b4')
        bars2 = ax.bar(x + width/2, fed_times, width, label='Fed-GWAS', color='#ff7f0e')
        
        ax.set_xlabel('Pipeline Stage')
        ax.set_ylabel('Time (hours)')
        ax.set_title(f'Performance Comparison: SF-GWAS vs Fed-GWAS ({self.scale})')
        ax.set_xticks(x)
        ax.set_xticklabels(stages)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add value labels on bars
        for bars in [bars1, bars2]:
            for bar in bars:
                height = bar.get_height()
                if height > 0:
                    ax.text(bar.get_x() + bar.get_width()/2., height,
                           f'{height:.2f}h',
                           ha='center', va='bottom', fontsize=8)
        
        plt.tight_layout()
        
        if output_file is None:
            output_file = str(self.federated_dir / "sf_gwas_comparison.png")
        
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        
        self.logger.info(f"Comparison plot saved to {output_file}")
        return output_file
    
    def generate_speedup_plot(self, output_file: Optional[str] = None) -> str:
        """Generate speedup/slowdown visualization"""
        if not self.comparison_results or 'speedup' not in self.comparison_results:
            self.logger.warning("Cannot generate speedup plot: comparison not run")
            return ""
        
        stages = ['QC', 'KING/PCA', 'LR', 'Total']
        speedups = [
            self.comparison_results['speedup'].get('qc', 1.0),
            self.comparison_results['speedup'].get('king', 1.0),
            self.comparison_results['speedup'].get('lr', 1.0),
            self.comparison_results['speedup'].get('total', 1.0)
        ]
        
        colors = ['green' if s > 1 else 'red' for s in speedups]
        
        fig, ax = plt.subplots(figsize=(10, 6))
        bars = ax.bar(stages, speedups, color=colors, alpha=0.7)
        
        ax.axhline(y=1.0, color='black', linestyle='--', linewidth=1, label='Baseline (1.0x)')
        ax.set_xlabel('Pipeline Stage')
        ax.set_ylabel('Speedup Factor (SF-GWAS / Fed-GWAS)')
        ax.set_title(f'Speedup Comparison: Fed-GWAS vs SF-GWAS ({self.scale})')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add value labels
        for bar, speedup in zip(bars, speedups):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{speedup:.2f}x',
                   ha='center', va='bottom' if speedup > 1 else 'top', fontsize=10)
        
        plt.tight_layout()
        
        if output_file is None:
            output_file = str(self.federated_dir / "speedup_comparison.png")
        
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        
        self.logger.info(f"Speedup plot saved to {output_file}")
        return output_file
    
    def generate_report(self, output_file: Optional[str] = None) -> str:
        """Generate comparison report"""
        if output_file is None:
            output_file = str(self.federated_dir / "sf_gwas_comparison_report.md")
        
        with open(output_file, 'w') as f:
            f.write("# SF-GWAS vs Fed-GWAS Performance Comparison\n\n")
            f.write(f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write(f"**Scale**: {self.scale}\n")
            f.write(f"**Benchmark Dataset**: {self.benchmark['notes']}\n")
            f.write(f"  - Samples: {self.benchmark['samples']:,}\n")
            f.write(f"  - SNPs: {self.benchmark['snps']:,}\n")
            f.write(f"  - Clients: {self.benchmark['clients']}\n\n")
            
            if not self.comparison_results:
                f.write("**Note**: Comparison not yet run. Run compare_performance() first.\n")
                return output_file
            
            f.write("## Performance Comparison\n\n")
            f.write("| Stage | SF-GWAS (hours) | Fed-GWAS (hours) | Difference | Speedup |\n")
            f.write("|-------|----------------|------------------|------------|--------|\n")
            
            stages_map = {
                'qc': 'QC',
                'king': 'KING/PCA',
                'lr': 'LR',
                'total': 'Total'
            }
            
            for stage_key, stage_name in stages_map.items():
                if stage_key in self.comparison_results.get('relative_performance', {}):
                    perf = self.comparison_results['relative_performance'][stage_key]
                    speedup = self.comparison_results['speedup'].get(stage_key, 1.0)
                    
                    f.write(f"| {stage_name} | {perf['sf_gwas_hours']:.2f} | "
                           f"{perf['federated_hours']:.2f} | "
                           f"{perf['difference_hours']:+.2f} | "
                           f"{speedup:.2f}x |\n")
            
            f.write("\n## Analysis\n\n")
            
            # Summary
            total_speedup = self.comparison_results['speedup'].get('total', 1.0)
            if total_speedup > 1.0:
                f.write(f"✓ **Fed-GWAS is {total_speedup:.2f}x faster overall**\n\n")
            elif total_speedup < 1.0:
                f.write(f"✗ **Fed-GWAS is {1/total_speedup:.2f}x slower overall**\n\n")
            else:
                f.write("**Fed-GWAS performance is similar to SF-GWAS**\n\n")
            
            # Stage-by-stage analysis
            f.write("### Stage-by-Stage Analysis\n\n")
            for stage_key, stage_name in stages_map.items():
                if stage_key in self.comparison_results.get('relative_performance', {}):
                    perf = self.comparison_results['relative_performance'][stage_key]
                    speedup = self.comparison_results['speedup'].get(stage_key, 1.0)
                    
                    f.write(f"#### {stage_name}\n")
                    if speedup > 1.0:
                        f.write(f"- Fed-GWAS is **{speedup:.2f}x faster** "
                               f"({perf['sf_gwas_hours']:.2f}h → {perf['federated_hours']:.2f}h)\n")
                    elif speedup < 1.0:
                        f.write(f"- Fed-GWAS is **{1/speedup:.2f}x slower** "
                               f"({perf['sf_gwas_hours']:.2f}h → {perf['federated_hours']:.2f}h)\n")
                    else:
                        f.write(f"- Performance is similar "
                               f"({perf['sf_gwas_hours']:.2f}h vs {perf['federated_hours']:.2f}h)\n")
                    f.write("\n")
            
            # Notes
            f.write("## Notes\n\n")
            f.write("- SF-GWAS benchmarks are from published results for UKB-scale data\n")
            f.write("- KING stage in Fed-GWAS is compared to PCA stage in SF-GWAS (similar functionality)\n")
            f.write("- Times are in hours\n")
            f.write("- Speedup > 1.0 means Fed-GWAS is faster\n")
        
        self.logger.info(f"Comparison report saved to {output_file}")
        return output_file
    
    def run_full_comparison(self) -> Dict[str, Any]:
        """Run complete comparison pipeline"""
        self.logger.info("Running full SF-GWAS comparison...")
        
        self.load_federated_metrics()
        comparison = self.compare_performance()
        
        if comparison:
            self.generate_comparison_plot()
            self.generate_speedup_plot()
            self.generate_report()
        
        self.logger.info("SF-GWAS comparison completed")
        return comparison


def main():
    """Command-line interface for SF-GWAS comparison"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Compare Fed-GWAS with SF-GWAS benchmarks")
    parser.add_argument("federated_results", type=str, help="Path to federated experiment results")
    parser.add_argument("--scale", type=str, default="ukb_scale", 
                       choices=['ukb_scale', 'small_scale'],
                       help="Scale of comparison")
    
    args = parser.parse_args()
    
    comparator = SFGWASComparator(args.federated_results, args.scale)
    comparison = comparator.run_full_comparison()
    
    print("\nSF-GWAS comparison completed!")
    if comparison:
        print(f"\nOverall speedup: {comparison['speedup'].get('total', 1.0):.2f}x")
        print(f"\nDetailed results saved to: {args.federated_results}")


if __name__ == "__main__":
    main()


