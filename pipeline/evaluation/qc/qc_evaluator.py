#!/usr/bin/env python3
"""
QC Evaluation Tool for Federated GWAS Experiments

Compares federated QC exclusion lists with centralized PLINK baselines:
- SNP exclusion list agreement (F1/precision/recall)
- Summary report
"""

import argparse
import logging
from pathlib import Path
from typing import Dict, Optional

import pandas as pd


def _setup_logger(name: str) -> logging.Logger:
    logger = logging.getLogger(name)
    if not logger.handlers:
        logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger


class QCEvaluator:
    """Evaluate QC SNP exclusion list agreement."""

    def __init__(self, federated_results_dir: str, baseline_results_dir: Optional[str] = None):
        self.federated_dir = Path(federated_results_dir)
        self.baseline_dir = Path(baseline_results_dir) if baseline_results_dir else None
        self.logger = _setup_logger("qc_evaluator")

        self.federated_qc: Optional[pd.DataFrame] = None
        self.baseline_qc: Optional[pd.DataFrame] = None
        self.analysis_results: Dict[str, Dict] = {}

    def _load_qc_results(self, results_dir: Path, is_baseline: bool) -> None:
        """Load QC exclusion results."""
        excluded_snps = set()

        if is_baseline:
            original_bim = results_dir / "merged.bim"
            filtered_bim = results_dir / "qc.bim"

            if original_bim.exists() and filtered_bim.exists():
                try:
                    original_snps = set()
                    with open(original_bim, 'r') as f:
                        for line in f:
                            parts = line.strip().split()
                            if len(parts) >= 2:
                                original_snps.add(parts[1])

                    filtered_snps = set()
                    with open(filtered_bim, 'r') as f:
                        for line in f:
                            parts = line.strip().split()
                            if len(parts) >= 2:
                                filtered_snps.add(parts[1])

                    excluded_snps = original_snps - filtered_snps
                    self.logger.info(
                        f"Baseline QC: Found {len(original_snps)} original SNPs, "
                        f"{len(filtered_snps)} after filtering, {len(excluded_snps)} excluded"
                    )

                    if excluded_snps:
                        self.baseline_qc = pd.DataFrame({'SNP': list(excluded_snps)})
                        self.logger.info(f"Baseline QC: {len(excluded_snps)} SNPs excluded")
                    else:
                        self.baseline_qc = pd.DataFrame({'SNP': []})
                        self.logger.info("Baseline QC: No SNPs excluded (all passed QC thresholds)")
                except Exception as e:
                    self.logger.warning(f"Could not compare BIM files: {e}")
                    self.baseline_qc = pd.DataFrame({'SNP': []})
            else:
                # Fallback: Try to extract from QC output files (less accurate)
                self.logger.warning("Original/filtered BIM files not found, trying to extract from QC output files")

                maf_threshold = 0.01
                missing_threshold = 0.05
                hwe_threshold = 1e-6

                config_file = results_dir.parent.parent / "config.yaml"
                if config_file.exists():
                    try:
                        import yaml
                        with open(config_file, 'r') as f:
                            config = yaml.safe_load(f)
                            thresholds = config.get('thresholds', {})
                            maf_threshold = thresholds.get('maf_threshold', maf_threshold)
                            missing_threshold = thresholds.get('missing_threshold', missing_threshold)
                            hwe_threshold = thresholds.get('hwe_threshold', hwe_threshold)
                    except Exception as e:
                        self.logger.warning(f"Could not load thresholds from config: {e}")

                lmiss_files = list(results_dir.glob("**/*.lmiss"))
                for lmiss_file in lmiss_files:
                    try:
                        df = pd.read_csv(lmiss_file, sep=r'\s+')
                        df.columns = df.columns.str.strip()
                        if 'F_MISS' in df.columns and 'SNP' in df.columns:
                            excluded = df[df['F_MISS'] > missing_threshold]['SNP'].tolist()
                            excluded_snps.update(excluded)
                    except Exception as e:
                        self.logger.warning(f"Could not parse LMiss file {lmiss_file}: {e}")

                frq_files = list(results_dir.glob("**/*.frq"))
                for frq_file in frq_files:
                    try:
                        df = pd.read_csv(frq_file, sep=r'\s+')
                        df.columns = df.columns.str.strip()
                        if 'MAF' in df.columns and 'SNP' in df.columns:
                            excluded = df[df['MAF'] < maf_threshold]['SNP'].tolist()
                            excluded_snps.update(excluded)
                    except Exception as e:
                        self.logger.warning(f"Could not parse FRQ file {frq_file}: {e}")

                hwe_files = list(results_dir.glob("**/*.hwe"))
                for hwe_file in hwe_files:
                    try:
                        df = pd.read_csv(hwe_file, sep=r'\s+')
                        df.columns = df.columns.str.strip()
                        if 'P' in df.columns and 'SNP' in df.columns:
                            if 'TEST' in df.columns:
                                df_all = df[df['TEST'] == 'ALL']
                            else:
                                df_all = df
                            excluded = df_all[df_all['P'] < hwe_threshold]['SNP'].tolist()
                            excluded_snps.update(excluded)
                    except Exception as e:
                        self.logger.warning(f"Could not parse HWE file {hwe_file}: {e}")

                if excluded_snps:
                    self.baseline_qc = pd.DataFrame({'SNP': list(excluded_snps)})
                else:
                    self.baseline_qc = pd.DataFrame({'SNP': []})
        else:
            exclusion_files = []
            for center_dir in results_dir.glob("center_*/logs"):
                exclusion_file = center_dir / "global_qc_exclusion.txt"
                if exclusion_file.exists():
                    exclusion_files.append(exclusion_file)

            if exclusion_files:
                for exclusion_file in exclusion_files:
                    try:
                        with open(exclusion_file, 'r') as f:
                            client_excluded = set(line.strip() for line in f if line.strip())
                            excluded_snps.update(client_excluded)
                    except Exception as e:
                        self.logger.warning(f"Could not read exclusion file {exclusion_file}: {e}")
                self.logger.info(
                    f"Loaded exclusion lists from {len(exclusion_files)} clients, {len(excluded_snps)} unique SNPs excluded"
                )
            else:
                self.logger.warning("No client exclusion files found, trying to reconstruct from client QC files")
                for center_dir in results_dir.glob("center_*/logs"):
                    lmiss_files = list(center_dir.glob("*.lmiss"))
                    for lmiss_file in lmiss_files:
                        try:
                            df = pd.read_csv(lmiss_file, sep=r'\s+')
                            df.columns = df.columns.str.strip()
                            if 'F_MISS' in df.columns and 'SNP' in df.columns:
                                excluded = df[df['F_MISS'] > 0.05]['SNP'].tolist()
                                excluded_snps.update(excluded)
                        except Exception as e:
                            self.logger.warning(f"Could not parse LMiss file {lmiss_file}: {e}")

            if excluded_snps:
                self.federated_qc = pd.DataFrame({'SNP': list(excluded_snps)})
                self.logger.info(f"Federated QC: {len(excluded_snps)} SNPs excluded")
            else:
                self.federated_qc = pd.DataFrame({'SNP': []})
                self.logger.info("Federated QC: No SNPs excluded (all passed QC thresholds)")

        if not hasattr(self, 'federated_qc') or self.federated_qc is None:
            self.federated_qc = pd.DataFrame({'SNP': []})
            self.logger.info("Federated QC: No exclusion data found, assuming no SNPs excluded")

    def analyze_qc_agreement(self) -> Dict[str, float]:
        if self.federated_qc is None or self.baseline_qc is None:
            self.logger.warning("Cannot analyze QC agreement: missing data")
            return {}

        fed_excluded = set(self.federated_qc['SNP'].tolist())
        baseline_excluded = set(self.baseline_qc['SNP'].tolist())
        self.logger.info(
            f"QC debug: federated_excluded={len(fed_excluded)} baseline_excluded={len(baseline_excluded)}"
        )

        true_positives = len(fed_excluded & baseline_excluded)
        false_positives = len(fed_excluded - baseline_excluded)
        false_negatives = len(baseline_excluded - fed_excluded)

        precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0.0
        recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0.0
        f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0

        agreement = {
            'f1_score': f1_score,
            'precision': precision,
            'recall': recall,
            'true_positives': true_positives,
            'false_positives': false_positives,
            'false_negatives': false_negatives,
            'federated_excluded_count': len(fed_excluded),
            'baseline_excluded_count': len(baseline_excluded)
        }
        self.analysis_results['qc_agreement'] = agreement
        return agreement

    def write_report(self, output_path: Path) -> None:
        agreement = self.analysis_results.get('qc_agreement', {})
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open('w') as f:
            f.write("# QC Evaluation Report\n\n")
            if not agreement:
                f.write("QC agreement could not be computed (missing data).\n")
                return
            f.write("## QC SNP Exclusion Agreement\n\n")
            f.write(f"- **F1 Score**: {agreement['f1_score']:.4f}\n")
            f.write(f"- **Precision**: {agreement['precision']:.4f}\n")
            f.write(f"- **Recall**: {agreement['recall']:.4f}\n")
            f.write(f"- **Federated Excluded**: {agreement['federated_excluded_count']}\n")
            f.write(f"- **Baseline Excluded**: {agreement['baseline_excluded_count']}\n")
            f.write("\n")
            f.write("### Confusion Counts\n\n")
            f.write(f"- **True Positives**: {agreement['true_positives']}\n")
            f.write(f"- **False Positives**: {agreement['false_positives']}\n")
            f.write(f"- **False Negatives**: {agreement['false_negatives']}\n")

    def run(self, report_path: Optional[Path] = None) -> Dict[str, float]:
        if not self.federated_dir.exists():
            raise FileNotFoundError(f"Federated results directory not found: {self.federated_dir}")

        self._load_qc_results(self.federated_dir, is_baseline=False)
        if self.baseline_dir is not None:
            self._load_qc_results(self.baseline_dir, is_baseline=True)

        results = self.analyze_qc_agreement()
        if report_path is not None:
            self.write_report(report_path)
        return results


def main() -> None:
    parser = argparse.ArgumentParser(description='QC evaluation for federated GWAS results')
    parser.add_argument('results_dir', type=str, help='Federated results directory')
    parser.add_argument('--baseline', type=str, required=True, help='Baseline directory')
    parser.add_argument('--report', type=str, default=None, help='Output report path (default: results_dir/qc_report.md)')

    args = parser.parse_args()
    results_dir = Path(args.results_dir)
    report_path = Path(args.report) if args.report else results_dir / 'qc_report.md'

    evaluator = QCEvaluator(str(results_dir), args.baseline)
    evaluator.run(report_path)


if __name__ == '__main__':
    main()
