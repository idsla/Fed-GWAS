#!/usr/bin/env python3
"""
Run QC + LR evaluation (and optional KING) in one shot.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Any, Dict, Optional


def _setup_logger() -> logging.Logger:
    logger = logging.getLogger("evaluate_all")
    if not logger.handlers:
        logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger


def _write_summary(
    report_path: Path,
    qc_results: Optional[Dict[str, Any]],
    lr_results: Optional[Dict[str, Any]],
    king_results: Optional[Dict[str, Any]],
) -> None:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w") as f:
        f.write("# Evaluation Summary\n\n")

        if qc_results:
            f.write("## QC Agreement\n\n")
            f.write(f"- **F1 Score**: {qc_results.get('f1_score', float('nan')):.4f}\n")
            f.write(f"- **Precision**: {qc_results.get('precision', float('nan')):.4f}\n")
            f.write(f"- **Recall**: {qc_results.get('recall', float('nan')):.4f}\n")
            f.write(f"- **Federated Excluded**: {qc_results.get('federated_excluded_count', 0)}\n")
            f.write(f"- **Baseline Excluded**: {qc_results.get('baseline_excluded_count', 0)}\n\n")

        if king_results:
            f.write("## KING (Accumulator vs Baseline)\n\n")
            f.write(f"- **Total accumulator pairs**: {king_results.get('total_pairs', 0)}\n")
            f.write(f"- **Mapped pairs**: {king_results.get('mapped_pairs', 0)}\n")
            f.write(f"- **Missing map pairs**: {king_results.get('missing_map_pairs', 0)}\n")
            f.write(f"- **Missing baseline pairs**: {king_results.get('missing_baseline_pairs', 0)}\n")
            f.write(f"- **Cross-client pairs (mapped)**: {king_results.get('cross_pairs', 0)}\n")
            f.write(f"- **Within-client pairs (mapped)**: {king_results.get('within_pairs', 0)}\n")
            f.write(f"- **Overlapping pairs used**: {king_results.get('overlapping_pairs', 0)}\n")
            pearson = king_results.get("pearson_r")
            if pearson is not None:
                f.write(f"- **Pearson r**: {pearson:.4f}\n")
            mae = king_results.get("mae")
            if mae is not None:
                f.write(f"- **MAE**: {mae:.6f}\n")
            max_abs = king_results.get("max_abs_diff")
            if max_abs is not None:
                f.write(f"- **Max abs diff**: {max_abs:.6f}\n")
            f.write("\n")

        if lr_results:
            f.write("## LR Correlation\n\n")
            local = lr_results.get("local")
            if isinstance(local, dict):
                f.write("### Local LR\n\n")
                f.write(f"- **Pearson r**: {local.get('pearson_r', float('nan')):.4f}\n")
                f.write(f"- **P-value**: {local.get('p_value', float('nan')):.4e}\n")
                f.write(f"- **MSE**: {local.get('mse', float('nan')):.4f}\n")
                f.write(f"- **MAE**: {local.get('mae', float('nan')):.4f}\n")
                f.write(f"- **n_snps**: {local.get('n_snps', 0)}\n")
                f.write(f"- **Top-100 overlap**: {local.get('top_100_overlap', 0.0):.4f}\n\n")

            global_lr = lr_results.get("global")
            if isinstance(global_lr, dict):
                f.write("### Global LR\n\n")
                f.write(f"- **Pearson r**: {global_lr.get('pearson_r', float('nan')):.4f}\n")
                f.write(f"- **P-value**: {global_lr.get('p_value', float('nan')):.4e}\n")
                f.write(f"- **MSE**: {global_lr.get('mse', float('nan')):.4f}\n")
                f.write(f"- **MAE**: {global_lr.get('mae', float('nan')):.4f}\n")
                f.write(f"- **n_snps**: {global_lr.get('n_snps', 0)}\n")
                f.write(f"- **Top-100 overlap**: {global_lr.get('top_100_overlap', 0.0):.4f}\n\n")

            sig = lr_results.get("significance_agreement")
            if isinstance(sig, dict):
                f.write("### LR Significance Agreement\n\n")
                f.write(f"- **Agreement**: {sig.get('agreement', float('nan')):.4f}\n")
                f.write(f"- **TP**: {sig.get('tp', 0)}\n")
                f.write(f"- **TN**: {sig.get('tn', 0)}\n")
                f.write(f"- **FP**: {sig.get('fp', 0)}\n")
                f.write(f"- **FN**: {sig.get('fn', 0)}\n")
                f.write(f"- **n_snps**: {sig.get('n_snps', 0)}\n\n")

            coverage = lr_results.get("coverage")
            if isinstance(coverage, dict) and coverage:
                f.write("### Coverage\n\n")
                for key, stats in coverage.items():
                    threshold = key.replace("threshold_", "")
                    f.write(f"**Threshold {threshold}:**\n")
                    f.write(f"- Coverage: {stats.get('coverage', 0.0):.4f}\n")
                    f.write(f"- Covered: {stats.get('covered', 0)}/{stats.get('total_baseline_sig', 0)} baseline significant SNPs\n")
                    f.write(f"- Total federated positive: {stats.get('total_fed_sig', 0)}\n\n")

            screening = lr_results.get("screening")
            if isinstance(screening, dict) and screening:
                f.write("### Screening Summary\n\n")
                local = screening.get("local")
                if isinstance(local, dict):
                    f.write("**Local Stage:**\n")
                    f.write(f"- Retained: {local.get('retained', 0)}/{local.get('total', 0)} ({local.get('ratio', 0.0):.4f})\n\n")
                global_stage = screening.get("global")
                if isinstance(global_stage, dict):
                    f.write("**Global Stage (Iterative):**\n")
                    f.write(f"- Analyzed: {global_stage.get('iterative', 0)}/{global_stage.get('total', 0)} ({global_stage.get('ratio', 0.0):.4f})\n\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run QC + LR (and optional KING) evaluation")
    parser.add_argument("results_dir", type=str, help="Federated results directory")
    parser.add_argument("--baseline", type=str, required=True, help="Baseline directory")
    parser.add_argument("--report", type=str, default=None, help="Summary report path")
    parser.add_argument("--no-plots", action="store_true", help="Skip LR plot generation")
    parser.add_argument("--king", action="store_true", help="Also run KING evaluation")
    parser.add_argument("--king-center-id", type=int, default=1, help="Center id for KING accumulator")
    parser.add_argument("--king-data-dir", type=str, default=None, help="Optional data root containing center_* dirs")
    parser.add_argument("--qc-only", action="store_true", help="Run QC only")
    parser.add_argument("--lr-only", action="store_true", help="Run LR only")
    parser.add_argument("--king-only", action="store_true", help="Run KING only")

    args = parser.parse_args()

    logger = _setup_logger()
    results_dir = Path(args.results_dir)
    baseline_dir = Path(args.baseline)
    if not results_dir.exists():
        raise FileNotFoundError(f"Results directory not found: {results_dir}")
    if not baseline_dir.exists():
        raise FileNotFoundError(f"Baseline directory not found: {baseline_dir}")

    base_dir = Path(__file__).resolve().parent
    sys.path.insert(0, str(base_dir / "qc"))
    sys.path.insert(0, str(base_dir / "lr"))
    sys.path.insert(0, str(base_dir / "king"))

    from qc_evaluator import QCEvaluator  # type: ignore
    from lr_evaluator import LREvaluator  # type: ignore

    only_flags = [args.qc_only, args.lr_only, args.king_only]
    if sum(1 for flag in only_flags if flag) > 1:
        raise ValueError("Only one of --qc-only, --lr-only, --king-only can be set.")

    run_qc = not any(only_flags) or args.qc_only
    run_lr = not any(only_flags) or args.lr_only
    run_king = args.king or args.king_only
    if args.king_only:
        run_qc = False
        run_lr = False
        run_king = True
    if args.qc_only:
        run_lr = False
        run_king = False
    if args.lr_only:
        run_qc = False
        run_king = False

    qc_results = None
    if run_qc:
        qc = QCEvaluator(str(results_dir), str(baseline_dir))
        qc_results = qc.run(report_path=results_dir / "qc_report.md")
        logger.info("QC evaluation complete")

    lr_results = None
    if run_lr:
        lr = LREvaluator(str(results_dir), str(baseline_dir))
        lr_results = lr.run(report_path=results_dir / "lr_report.md", make_plots=not args.no_plots)
        logger.info("LR evaluation complete")

    king_results = None
    if run_king:
        from compare_king_from_accumulator import compare_from_accumulator  # type: ignore

        king_report = results_dir / "king_report.md"
        king_results = compare_from_accumulator(
            results_dir=results_dir,
            baseline_dir=baseline_dir,
            center_id=args.king_center_id,
            report_path=king_report,
            data_dir=Path(args.king_data_dir) if args.king_data_dir else None,
        )
        logger.info("KING evaluation complete")

    report_path = Path(args.report) if args.report else results_dir / "evaluation_report.md"
    _write_summary(report_path, qc_results, lr_results, king_results)
    logger.info(f"Summary report saved to {report_path}")


if __name__ == "__main__":
    main()
