#!/usr/bin/env python3
"""
Analyze how varying the federated screening threshold affects coverage and filtering.

This script tests different threshold values to find the optimal balance between:
- Coverage: Preserving baseline significant SNPs
- Filtering: Reducing false positives
"""

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Set publication-quality style
try:
    plt.style.use('seaborn-v0_8-paper')
except OSError:
    try:
        plt.style.use('seaborn-paper')
    except OSError:
        plt.style.use('seaborn')
sns.set_palette("husl")

matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['savefig.bbox'] = 'tight'


def setup_logger() -> logging.Logger:
    logger = logging.getLogger("threshold_analyzer")
    if not logger.handlers:
        logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger


def load_data(fed_file: Path, baseline_file: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load federated and baseline LR results."""
    # Load federated (de-anonymized)
    fed = pd.read_csv(fed_file, sep='\t')
    # Use p_min as the representative p-value
    if 'p_min' in fed.columns:
        fed['P'] = fed['p_min']
    elif 'p_values' in fed.columns:
        # Extract first p-value from list if needed
        fed['P'] = fed['p_values'].apply(lambda x: float(str(x).strip('[]').split(',')[0]) if isinstance(x, str) else x)
    
    # Load baseline
    baseline = pd.read_csv(baseline_file, sep=r'\s+')
    
    return fed, baseline


def compute_coverage_at_thresholds(
    fed: pd.DataFrame,
    baseline: pd.DataFrame,
    baseline_thresholds: List[float],
    fed_thresholds: List[float]
) -> pd.DataFrame:
    """Compute coverage metrics for different threshold combinations."""
    # Merge by SNP ID
    merged = pd.merge(
        fed[['SNP', 'P']].rename(columns={'P': 'P_fed'}),
        baseline[['SNP', 'P']].rename(columns={'P': 'P_baseline'}),
        on='SNP',
        how='inner'
    )
    
    results = []
    
    for base_thresh in baseline_thresholds:
        # Get baseline significant SNPs at this threshold
        baseline_sig = merged[merged['P_baseline'] < base_thresh]
        n_baseline_sig = len(baseline_sig)
        
        if n_baseline_sig == 0:
            continue
        
        for fed_thresh in fed_thresholds:
            # Count how many baseline significant SNPs are preserved
            preserved = baseline_sig[baseline_sig['P_fed'] < fed_thresh]
            n_preserved = len(preserved)
            
            # Count total federated significant
            n_fed_sig = len(merged[merged['P_fed'] < fed_thresh])
            
            # Coverage = preserved / baseline significant
            coverage = n_preserved / n_baseline_sig if n_baseline_sig > 0 else 0.0
            
            # Filtering efficiency = baseline significant / federated significant
            # Higher is better (more focused on true positives)
            filtering_ratio = n_baseline_sig / n_fed_sig if n_fed_sig > 0 else 0.0
            
            results.append({
                'baseline_threshold': base_thresh,
                'federated_threshold': fed_thresh,
                'baseline_significant': n_baseline_sig,
                'preserved': n_preserved,
                'coverage': coverage,
                'federated_significant': n_fed_sig,
                'filtering_ratio': filtering_ratio,
                'false_positives': n_fed_sig - n_preserved,
            })
    
    return pd.DataFrame(results)


def plot_threshold_sensitivity(
    results_df: pd.DataFrame,
    output_path: Path,
    baseline_threshold: float = 5e-8
) -> None:
    """Create plots showing threshold sensitivity."""
    # Filter for the specified baseline threshold
    df = results_df[results_df['baseline_threshold'] == baseline_threshold].copy()
    
    if df.empty:
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: Coverage vs Federated Threshold
    ax1 = axes[0, 0]
    ax1.plot(df['federated_threshold'], df['coverage'], 'o-', linewidth=2, markersize=6)
    ax1.axhline(y=1.0, color='g', linestyle='--', alpha=0.5, label='Perfect coverage (1.0)')
    ax1.set_xlabel('Federated Screening Threshold (p-value)')
    ax1.set_ylabel('Coverage')
    ax1.set_title(f'Coverage vs Federated Threshold\n(Baseline threshold: p < {baseline_threshold})')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_xscale('log')
    
    # Plot 2: Number of SNPs vs Federated Threshold
    ax2 = axes[0, 1]
    ax2.plot(df['federated_threshold'], df['preserved'], 'o-', label='Preserved', linewidth=2, markersize=6)
    ax2.plot(df['federated_threshold'], df['federated_significant'], 's-', label='Total Federated Significant', 
             linewidth=2, markersize=6, alpha=0.7)
    ax2.axhline(y=df['baseline_significant'].iloc[0], color='r', linestyle='--', 
                label=f'Baseline Significant ({int(df["baseline_significant"].iloc[0])})', alpha=0.5)
    ax2.set_xlabel('Federated Screening Threshold (p-value)')
    ax2.set_ylabel('Number of SNPs')
    ax2.set_title('SNP Counts vs Federated Threshold')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_xscale('log')
    
    # Plot 3: Filtering Ratio vs Federated Threshold
    ax3 = axes[1, 0]
    ax3.plot(df['federated_threshold'], df['filtering_ratio'], 'o-', linewidth=2, markersize=6, color='orange')
    ax3.set_xlabel('Federated Screening Threshold (p-value)')
    ax3.set_ylabel('Filtering Ratio (Baseline Sig / Federated Sig)')
    ax3.set_title('Filtering Efficiency vs Federated Threshold\n(Higher = Better filtering)')
    ax3.grid(True, alpha=0.3)
    ax3.set_xscale('log')
    
    # Plot 4: Coverage vs Filtering Trade-off
    ax4 = axes[1, 1]
    scatter = ax4.scatter(df['filtering_ratio'], df['coverage'], 
                         c=df['federated_threshold'], s=100, alpha=0.6, cmap='viridis')
    ax4.set_xlabel('Filtering Ratio')
    ax4.set_ylabel('Coverage')
    ax4.set_title('Coverage vs Filtering Trade-off\n(Color = Federated Threshold)')
    ax4.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax4, label='Federated Threshold')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Analyze threshold sensitivity for federated GWAS screening"
    )
    parser.add_argument("fed_file", type=str, help="Federated LR results file (de-anonymized)")
    parser.add_argument("baseline_file", type=str, help="Baseline LR results file")
    parser.add_argument("--output", type=str, default=None, help="Output directory")
    parser.add_argument("--baseline-thresholds", type=float, nargs='+', 
                       default=[5e-8, 5e-5, 0.0005, 0.005, 0.05],
                       help="Baseline significance thresholds to test")
    parser.add_argument("--fed-thresholds", type=float, nargs='+',
                       default=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
                       help="Federated screening thresholds to test")
    
    args = parser.parse_args()
    logger = setup_logger()
    
    fed_file = Path(args.fed_file)
    baseline_file = Path(args.baseline_file)
    output_dir = Path(args.output) if args.output else fed_file.parent / "threshold_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Loading data...")
    fed, baseline = load_data(fed_file, baseline_file)
    logger.info(f"Federated: {len(fed):,} SNPs, Baseline: {len(baseline):,} SNPs")
    
    logger.info(f"Computing coverage at different thresholds...")
    results_df = compute_coverage_at_thresholds(
        fed, baseline, args.baseline_thresholds, args.fed_thresholds
    )
    
    # Save results table
    results_file = output_dir / "threshold_sensitivity_results.csv"
    results_df.to_csv(results_file, index=False)
    logger.info(f"Results saved to {results_file}")
    
    # Generate plots for genome-wide significance
    plot_file = output_dir / "threshold_sensitivity_5e-8.png"
    plot_threshold_sensitivity(results_df, plot_file, baseline_threshold=5e-8)
    logger.info(f"Plot saved to {plot_file}")
    
    # Generate summary table for manuscript
    summary_df = results_df[results_df['baseline_threshold'] == 5e-8].copy()
    summary_df = summary_df.sort_values('federated_threshold')
    
    summary_table = summary_df[[
        'federated_threshold', 'coverage', 'preserved', 
        'baseline_significant', 'federated_significant', 'filtering_ratio'
    ]].copy()
    summary_table.columns = [
        'Federated Threshold', 'Coverage', 'Preserved', 
        'Baseline Significant', 'Federated Significant', 'Filtering Ratio'
    ]
    summary_table['Coverage'] = summary_table['Coverage'].apply(lambda x: f"{x:.4f}")
    summary_table['Filtering Ratio'] = summary_table['Filtering Ratio'].apply(lambda x: f"{x:.4f}")
    
    summary_file = output_dir / "threshold_sensitivity_summary.md"
    with open(summary_file, 'w') as f:
        f.write("# Threshold Sensitivity Analysis\n\n")
        f.write("## Genome-wide Significance (p < 5×10⁻⁸)\n\n")
        # Create markdown table manually
        f.write("| " + " | ".join(summary_table.columns) + " |\n")
        f.write("| " + " | ".join(["---"] * len(summary_table.columns)) + " |\n")
        for _, row in summary_table.iterrows():
            f.write("| " + " | ".join(str(val) for val in row.values) + " |\n")
        f.write("\n\n")
        f.write("## Key Findings\n\n")
        
        # Find optimal threshold (highest coverage with good filtering)
        optimal = summary_df.loc[summary_df['coverage'].idxmax()]
        f.write(f"- **Maximum Coverage**: {optimal['coverage']:.4f} at threshold p = {optimal['federated_threshold']:.3f}\n")
        f.write(f"  - Preserved: {int(optimal['preserved'])}/{int(optimal['baseline_significant'])} baseline significant SNPs\n")
        f.write(f"  - Filtering ratio: {optimal['filtering_ratio']:.4f}\n\n")
        
        # Find threshold with 100% coverage
        perfect_coverage = summary_df[summary_df['coverage'] >= 0.9999]
        if not perfect_coverage.empty:
            best = perfect_coverage.loc[perfect_coverage['filtering_ratio'].idxmax()]
            f.write(f"- **Perfect Coverage (≥99.99%) with Best Filtering**: Threshold p = {best['federated_threshold']:.3f}\n")
            f.write(f"  - Coverage: {best['coverage']:.4f}\n")
            f.write(f"  - Filtering ratio: {best['filtering_ratio']:.4f}\n")
            f.write(f"  - Total federated significant: {int(best['federated_significant']):,} SNPs\n\n")
    
    logger.info(f"Summary saved to {summary_file}")
    
    # Print recommendations
    logger.info("\n" + "="*60)
    logger.info("RECOMMENDATIONS:")
    logger.info("="*60)
    
    optimal = summary_df.loc[summary_df['coverage'].idxmax()]
    logger.info(f"Maximum coverage ({optimal['coverage']:.4f}) at threshold p = {optimal['federated_threshold']:.3f}")
    
    perfect_coverage = summary_df[summary_df['coverage'] >= 0.9999]
    if not perfect_coverage.empty:
        best = perfect_coverage.loc[perfect_coverage['filtering_ratio'].idxmax()]
        logger.info(f"\nFor 100% coverage with best filtering:")
        logger.info(f"  Threshold: p = {best['federated_threshold']:.3f}")
        logger.info(f"  Coverage: {best['coverage']:.4f}")
        logger.info(f"  Filtering ratio: {best['filtering_ratio']:.4f}")
        logger.info(f"  Total federated significant: {int(best['federated_significant']):,} SNPs")
    
    logger.info("\n" + "="*60)


if __name__ == "__main__":
    main()
