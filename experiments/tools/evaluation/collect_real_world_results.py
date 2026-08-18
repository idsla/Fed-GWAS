#!/usr/bin/env python3
"""
Collect Real-World Experiment Results

Script to collect and analyze results from real-world dataset experiments.
Supports cluster deployment where results are on multiple nodes.
"""

import argparse
import logging
import shutil
import sys
from pathlib import Path
from typing import List, Optional
import json

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def collect_client_results(
    results_dir: Path,
    output_dir: Path,
    client_ids: List[str]
) -> None:
    """
    Collect results from client nodes.
    
    Args:
        results_dir: Base results directory
        output_dir: Output directory for collected results
        client_ids: List of client IDs (e.g., ['center_1', 'center_2'])
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Collecting results from {results_dir}")
    
    for client_id in client_ids:
        client_dir = results_dir / client_id
        if not client_dir.exists():
            logger.warning(f"Client directory not found: {client_dir}")
            continue
        
        # Collect logs
        log_dir = client_dir / "logs"
        if log_dir.exists():
            dest_log_dir = output_dir / client_id / "logs"
            dest_log_dir.mkdir(parents=True, exist_ok=True)
            
            # Copy key result files
            key_files = [
                "lr_results_client_deanon.txt",
                "king_results_client_*.txt",
                "*_log.txt",
                "qc_freq_*.frq",
                "mr_*.imiss",
                "mr_*.lmiss"
            ]
            
            for pattern in key_files:
                for src_file in log_dir.glob(pattern):
                    dest_file = dest_log_dir / src_file.name
                    shutil.copy2(src_file, dest_file)
                    logger.info(f"Copied: {src_file.name}")
        
        # Collect intermediate files (optional, can be large)
        intermediate_dir = client_dir / "intermediate"
        if intermediate_dir.exists():
            logger.info(f"Intermediate files found in {intermediate_dir} (not copied to save space)")


def collect_server_results(
    results_dir: Path,
    output_dir: Path
) -> None:
    """
    Collect results from server node.
    
    Args:
        results_dir: Base results directory
        output_dir: Output directory for collected results
    """
    server_dir = results_dir / "server"
    if not server_dir.exists():
        logger.warning(f"Server directory not found: {server_dir}")
        return
    
    dest_server_dir = output_dir / "server"
    dest_server_dir.mkdir(parents=True, exist_ok=True)
    
    # Collect aggregated results
    intermediate_dir = server_dir / "intermediate"
    if intermediate_dir.exists():
        dest_intermediate = dest_server_dir / "intermediate"
        dest_intermediate.mkdir(parents=True, exist_ok=True)
        
        # Copy latest run results
        run_dirs = sorted(intermediate_dir.glob("*"), key=lambda p: p.stat().st_mtime, reverse=True)
        if run_dirs:
            latest_run = run_dirs[0]
            dest_run = dest_intermediate / latest_run.name
            if latest_run.is_dir():
                shutil.copytree(latest_run, dest_run, dirs_exist_ok=True)
                logger.info(f"Copied server results: {latest_run.name}")
    
    # Collect logs
    log_dir = server_dir / "logs"
    if log_dir.exists():
        dest_log_dir = dest_server_dir / "logs"
        dest_log_dir.mkdir(parents=True, exist_ok=True)
        
        for log_file in log_dir.glob("*.txt"):
            shutil.copy2(log_file, dest_log_dir / log_file.name)
            logger.info(f"Copied server log: {log_file.name}")


def run_evaluators(
    results_dir: Path,
    baseline_dir: Optional[Path] = None
) -> None:
    """
    Run QC and LR evaluators on collected results.
    
    Args:
        results_dir: Results directory
        baseline_dir: Optional baseline directory for comparison
    """
    try:
        qc_path = Path(__file__).parent / "qc"
        lr_path = Path(__file__).parent / "lr"
        sys.path.insert(0, str(qc_path))
        sys.path.insert(0, str(lr_path))

        from qc_evaluator import QCEvaluator  # type: ignore
        from lr_evaluator import LREvaluator  # type: ignore

        if baseline_dir is None:
            logger.error("Baseline directory is required for evaluation")
            return

        logger.info("Running QC evaluator...")
        qc = QCEvaluator(str(results_dir), str(baseline_dir))
        qc.run(report_path=results_dir / "qc_report.md")

        logger.info("Running LR evaluator...")
        lr = LREvaluator(str(results_dir), str(baseline_dir))
        lr.run(report_path=results_dir / "lr_report.md")

        logger.info("QC/LR evaluation complete")
    except ImportError as e:
        logger.error(f"Evaluator import failed: {e}")
    except Exception as e:
        logger.error(f"Evaluation failed: {e}")


def main():
    parser = argparse.ArgumentParser(
        description='Collect and analyze real-world experiment results'
    )
    parser.add_argument(
        'results_dir',
        type=str,
        help='Results directory (e.g., experiments/real_world/1000genomes/results)'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='Output directory for collected results (default: results_dir/collected)'
    )
    parser.add_argument(
        '--clients',
        type=str,
        nargs='+',
        default=['center_1', 'center_2'],
        help='Client IDs to collect (default: center_1 center_2)'
    )
    parser.add_argument(
        '--baseline-dir',
        type=str,
        default=None,
        help='Baseline directory for comparison'
    )
    parser.add_argument(
        '--analyze',
        action='store_true',
        help='Run QC/LR evaluators after collection'
    )
    parser.add_argument(
        '--server-only',
        action='store_true',
        help='Collect only server results'
    )
    parser.add_argument(
        '--clients-only',
        action='store_true',
        help='Collect only client results'
    )
    
    args = parser.parse_args()
    
    results_path = Path(args.results_dir)
    if not results_path.exists():
        logger.error(f"Results directory not found: {results_path}")
        sys.exit(1)
    
    if args.output_dir:
        output_path = Path(args.output_dir)
    else:
        output_path = results_path / "collected"
    
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Collect results
    if not args.server_only:
        collect_client_results(results_path, output_path, args.clients)
    
    if not args.clients_only:
        collect_server_results(results_path, output_path)
    
    logger.info(f"Results collected to: {output_path}")
    
    # Run analyzer if requested
    if args.analyze:
        baseline_path = Path(args.baseline_dir) if args.baseline_dir else None
        run_evaluators(output_path, baseline_path)


if __name__ == '__main__':
    main()
