#!/usr/bin/env python3
"""
LR Evaluation Tool for Federated GWAS Experiments

Compares federated LR results with centralized PLINK baselines:
- Local and global p-value correlation
- Significance agreement and coverage
- Manhattan and correlation plots
"""

import argparse
import csv
import hashlib
import logging
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yaml


def _setup_logger(name: str) -> logging.Logger:
    logger = logging.getLogger(name)
    if not logger.handlers:
        logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger


# Utility to parse average float from list-like string fields such as "[0.4751, 0.5]" -> 0.48755
def _parse_average_float(val: Any) -> Optional[float]:
    try:
        if isinstance(val, float):
            return val
        s = str(val).strip()
        if s.startswith("[") and s.endswith("]"):
            s = s[1:-1]
        if "," in s:
            parts = [p.strip() for p in s.split(",")]
            if parts:
                values = [float(p) for p in parts if p]
                if values:
                    return sum(values) / len(values)
        s = s.strip()
        if s == "":
            return None
        return float(s)
    except Exception:
        return None


def _parse_majority_significant(val: Any) -> Optional[bool]:
    try:
        if isinstance(val, bool):
            return val
        s = str(val).strip()
        if s.startswith("[") and s.endswith("]"):
            s = s[1:-1]
        import re
        matches = re.findall(r'\([^)]+,\s*(True|False)\)', s)
        if matches:
            true_count = sum(1 for m in matches if m.strip() == 'True')
            false_count = sum(1 for m in matches if m.strip() == 'False')
            return true_count > false_count
        return None
    except Exception:
        return None


def anonymize_snp_id(old_snp: str, global_seed: int) -> str:
    m = hashlib.sha256()
    seed_str = f"{global_seed}-{old_snp}"
    m.update(seed_str.encode("utf-8"))
    short_int = int.from_bytes(m.digest()[:4], "big")
    return f"rs{short_int}"


class LREvaluator:
    """Evaluate LR agreement and correlations."""

    def __init__(self, federated_results_dir: str, baseline_results_dir: Optional[str] = None):
        self.federated_dir = Path(federated_results_dir)
        self.baseline_dir = Path(baseline_results_dir) if baseline_results_dir else None
        self.logger = _setup_logger("lr_evaluator")

        self.analysis_results: Dict[str, Dict] = {}
        self.global_seed: Optional[int] = self._load_global_seed()

        self.federated_lr_local: Optional[pd.DataFrame] = None
        self.federated_lr_local_filtered: Optional[pd.DataFrame] = None
        self.federated_lr_global_deanon: Optional[pd.DataFrame] = None
        self.federated_lr_global: Optional[pd.DataFrame] = None
        self.baseline_lr: Optional[pd.DataFrame] = None
        self.baseline_lr_anon: Optional[pd.DataFrame] = None
        self.local_total_snps: Optional[int] = None
        self.local_retained_snps: Optional[int] = None
        self.global_iterative_snps: Optional[int] = None

    def _load_global_seed(self) -> Optional[int]:
        try:
            candidates = [
                self.federated_dir.parent / "config.yaml",
                self.federated_dir.parent.parent / "config.yaml",
            ]
            for cfg_path in candidates:
                if cfg_path.exists():
                    with cfg_path.open("r") as f:
                        cfg = yaml.safe_load(f)
                    seed = cfg.get("global_seed")
                    if seed is not None:
                        self.logger.info(f"Loaded global_seed={seed} from {cfg_path}")
                        return int(seed)
        except Exception as e:
            self.logger.warning(f"Could not load global_seed: {e}")
        return None

    def _load_snp_maps(self, results_dir: Path) -> Dict[str, str]:
        mapping: Dict[str, str] = {}
        for center_dir in results_dir.glob("center_*"):
            interm_root = center_dir / "intermediate"
            if not interm_root.exists():
                continue
            # Prefer latest run-scoped intermediate directory if present
            run_dirs = [p for p in interm_root.iterdir() if p.is_dir() and p.name.startswith("run_")]
            if run_dirs:
                interm = max(run_dirs, key=lambda p: p.stat().st_mtime)
            else:
                interm = interm_root
            for map_file in interm.glob("*snp_map.tsv"):
                try:
                    with open(map_file, "r") as f:
                        header = f.readline()
                        for line in f:
                            parts = line.strip().split()
                            if len(parts) >= 2:
                                anon, orig = parts[0], parts[1]
                                if anon not in mapping:
                                    mapping[anon] = orig
                    self.logger.info(f"Loaded SNP map with {len(mapping)} entries from {map_file}")
                except Exception as e:
                    self.logger.warning(f"Failed to load SNP map {map_file}: {e}")
        return mapping

    def _load_lr_results(self, results_dir: Path, is_baseline: bool) -> None:
        if is_baseline:
            lr_files = list(results_dir.glob("**/*.assoc.logistic"))
            if not lr_files:
                self.logger.warning(f"No LR files found in {results_dir}")
                return
            lr_file = max(lr_files, key=lambda f: f.stat().st_size if f.exists() else 0)
            try:
                df = pd.read_csv(lr_file, sep=r'\s+')
                df.columns = df.columns.str.strip()
                self.baseline_lr = df
                self.logger.info(f"Loaded baseline LR results: {len(df)} SNPs from {lr_file}")
                if self.global_seed is not None and 'SNP' in df.columns:
                    anon_df = df.copy()
                    anon_df['SNP'] = anon_df['SNP'].apply(lambda x: anonymize_snp_id(str(x), self.global_seed))
                    self.baseline_lr_anon = anon_df
                    self.logger.info(f"Generated anonymized baseline LR using global_seed={self.global_seed}")
            except Exception as e:
                self.logger.warning(f"Could not parse LR file {lr_file}: {e}")
            return

        # Federated LR
        client_deanon_files = []
        for center_dir in results_dir.glob("center_*/logs"):
            client_deanon_files.extend(list(center_dir.glob("lr_results_client_deanon.txt")))
        if client_deanon_files:
            deanon_dfs = []
            for fpath in client_deanon_files:
                try:
                    try:
                        csv.field_size_limit(10**7)
                    except (OverflowError, ValueError):
                        csv.field_size_limit(2**20)

                    rows = []
                    with open(fpath, 'r', encoding='utf-8') as f:
                        header = f.readline().strip().split('\t')
                        for line in f:
                            if not line.strip():
                                continue
                            parts = line.strip().split('\t', maxsplit=3)
                            if len(parts) >= 2:
                                snp = parts[0]
                                p_values_str = parts[1] if len(parts) > 1 else ''
                                iterations_sig_str = parts[2] if len(parts) > 2 else ''
                                p_val = _parse_average_float(p_values_str)
                                majority_sig = _parse_majority_significant(iterations_sig_str)
                                if p_val is not None:
                                    row_data = {'SNP': snp, 'P': p_val, 'p_values': p_values_str}
                                    if majority_sig is not None:
                                        row_data['is_significant'] = majority_sig
                                    rows.append(row_data)
                    if rows:
                        df = pd.DataFrame(rows)
                        deanon_dfs.append(df)
                except Exception as e:
                    self.logger.warning(f"Could not parse client deanon LR file {fpath}: {e}")
            if deanon_dfs:
                combined = pd.concat(deanon_dfs, ignore_index=True)
                combined = combined.drop_duplicates(subset=['SNP'], keep='first')
                cols = ['SNP', 'P']
                if 'is_significant' in combined.columns:
                    cols.append('is_significant')
                cols.extend([c for c in combined.columns if c not in ('SNP', 'P', 'is_significant', 'p_values')])
                if 'p_values' in combined.columns:
                    cols.append('p_values')
                combined = combined[cols]
                self.federated_lr_global_deanon = combined
                self.logger.info(
                    f"Loaded client de-anonymized global LR results: {len(combined)} SNPs from {len(client_deanon_files)} files"
                )

        # Local LR from clients
        local_lr_files = []
        retained_snp_sets = []
        for center_dir in results_dir.glob("center_*/logs"):
            local_lr_files.extend(list(center_dir.glob("local_lr*.assoc.logistic")))
            lr_filtered_bim = center_dir / "lr_filtered.bim"
            if lr_filtered_bim.exists():
                try:
                    retained_df = pd.read_csv(lr_filtered_bim, sep=r'\s+', header=None)
                    if retained_df.shape[1] >= 2:
                        retained_snps = set(retained_df.iloc[:, 1].astype(str).tolist())
                        retained_snp_sets.append(retained_snps)
                        self.logger.info(f"Loaded retained SNPs after local LR filter from {lr_filtered_bim}: {len(retained_snps)} SNPs kept")
                except Exception as e:
                    self.logger.warning(f"Could not parse retained SNPs from {lr_filtered_bim}: {e}")

        if local_lr_files:
            local_dfs = []
            for lr_file in local_lr_files:
                try:
                    df = pd.read_csv(lr_file, sep=r'\s+')
                    df.columns = df.columns.str.strip()
                    local_dfs.append(df)
                except Exception as e:
                    self.logger.warning(f"Could not parse local LR file {lr_file}: {e}")

            if local_dfs:
                combined_local = pd.concat(local_dfs, ignore_index=True)
                combined_local = combined_local.drop_duplicates(subset=['SNP'], keep='first')
                self.federated_lr_local = combined_local
                self.local_total_snps = len(combined_local)
                self.logger.info(f"Loaded local LR results: {len(combined_local)} unique SNPs from {len(local_lr_files)} client files")

                if retained_snp_sets:
                    retained_union = set().union(*retained_snp_sets)
                    self.local_retained_snps = len(retained_union)
                    filtered_out_df = combined_local[~combined_local['SNP'].isin(retained_union)]
                    self.federated_lr_local_filtered = filtered_out_df
                    self.logger.info(f"Local LR filtered-out set derived from retained SNPs: {len(filtered_out_df)} SNPs (retained_union={len(retained_union)})")
        else:
            self.logger.warning("No local LR files found in client logs")

    def analyze_lr_correlation(self) -> Dict[str, Any]:
        results: Dict[str, Any] = {}

        baseline_lr_for_merge = self.baseline_lr_anon if self.baseline_lr_anon is not None else self.baseline_lr
        local_lr_threshold = 0.30

        p_threshold = 0.3
        try:
            # Try multiple config paths - check client configs first (most specific)
            config_paths = [
                self.federated_dir.parent.parent / "configs" / "center_1" / "config.yaml",
                self.federated_dir.parent.parent / "configs" / "center_2" / "config.yaml",
                self.federated_dir.parent / "config.yaml",
                self.federated_dir.parent.parent / "config.yaml",
            ]
            for config_file in config_paths:
                if config_file.exists():
                    with open(config_file, 'r') as f:
                        config = yaml.safe_load(f)
                        if 'thresholds' in config and 'p_threshold' in config['thresholds']:
                            p_threshold = float(config['thresholds']['p_threshold'])
                            self.logger.info(f"Using p_threshold={p_threshold} from config: {config_file}")
                            break
        except Exception as e:
            self.logger.warning(f"Could not load p_threshold from config, using default 0.3: {e}")

        significance_threshold = p_threshold

        def compute_significance_agreement(fed_df: pd.DataFrame,
                                           base_df: pd.DataFrame,
                                           fed_label: str,
                                           fed_significance_col: Optional[str] = None) -> Dict[str, Any]:
            if fed_df is None or base_df is None:
                return {}
            if not {'SNP'}.issubset(fed_df.columns) or not {'SNP', 'P'}.issubset(base_df.columns):
                return {}

            fed_cols = ['SNP']
            if fed_significance_col and fed_significance_col in fed_df.columns:
                fed_cols.append(fed_significance_col)
            if 'P' in fed_df.columns:
                fed_cols.append('P')

            merged = pd.merge(
                fed_df[fed_cols],
                base_df[['SNP', 'P']],
                on='SNP',
                how='inner',
                suffixes=('_fed', '_baseline')
            )

            if len(merged) == 0:
                self.logger.warning(f"{fed_label} significance agreement: no overlapping SNPs")
                return {}

            merged['sig_baseline'] = merged['P_baseline'] < significance_threshold

            if fed_significance_col and fed_significance_col in merged.columns:
                merged['sig_fed'] = merged[fed_significance_col]
            elif 'P_fed' in merged.columns:
                merged['sig_fed'] = merged['P_fed'] < significance_threshold
            else:
                self.logger.warning(f"Cannot determine federated significance: missing {fed_significance_col} or P_fed")
                return {}

            tp = int(((merged['sig_fed'] == True) & (merged['sig_baseline'] == True)).sum())
            tn = int(((merged['sig_fed'] == False) & (merged['sig_baseline'] == False)).sum())
            fp = int(((merged['sig_fed'] == True) & (merged['sig_baseline'] == False)).sum())
            fn = int(((merged['sig_fed'] == False) & (merged['sig_baseline'] == True)).sum())

            total = tp + tn + fp + fn
            agreement = (tp + tn) / total if total > 0 else float('nan')

            return {
                'agreement': agreement,
                'tp': tp,
                'tn': tn,
                'fp': fp,
                'fn': fn,
                'n_snps': total,
                'baseline_threshold': significance_threshold,
                'federated_threshold': significance_threshold,
            }

        # Local LR correlation
        if self.federated_lr_local is not None and baseline_lr_for_merge is not None:
            if 'P' in self.federated_lr_local.columns and 'P' in baseline_lr_for_merge.columns:
                if 'SNP' in self.federated_lr_local.columns and 'SNP' in baseline_lr_for_merge.columns:
                    if self.federated_lr_local_filtered is not None:
                        fed_local_filtered = self.federated_lr_local_filtered
                        self.logger.info(f"Local LR comparison (filtered-out from lr_filtered.bim): {len(fed_local_filtered)} SNPs")
                    else:
                        fed_local_filtered = self.federated_lr_local[self.federated_lr_local['P'] >= local_lr_threshold]
                        self.logger.info(f"Local LR comparison (filtered by P >= {local_lr_threshold}): {len(fed_local_filtered)} SNPs")

                    merged = pd.merge(
                        fed_local_filtered[['SNP', 'P']],
                        baseline_lr_for_merge[['SNP', 'P']],
                        on='SNP',
                        how='inner',
                        suffixes=('_fed', '_baseline')
                    )

                    if len(merged) > 0:
                        initial_count = len(merged)
                        merged = merged.dropna(subset=['P_fed', 'P_baseline'])
                        merged = merged[(merged['P_fed'] > 0) & (merged['P_fed'] <= 1) &
                                        (merged['P_baseline'] > 0) & (merged['P_baseline'] <= 1)]
                        valid_count = len(merged)

                        if valid_count == 0:
                            correlation = 0.0
                            p_value = 1.0
                            mse = float('nan')
                            mae = float('nan')
                            top_overlap = 0.0
                        else:
                            merged['neg_log10_fed'] = -np.log10(merged['P_fed'] + 1e-300)
                            merged['neg_log10_baseline'] = -np.log10(merged['P_baseline'] + 1e-300)

                            fed_std = merged['neg_log10_fed'].std()
                            baseline_std = merged['neg_log10_baseline'].std()

                            if fed_std > 0 and baseline_std > 0:
                                correlation, p_value = pearsonr(merged['neg_log10_fed'], merged['neg_log10_baseline'])
                                if np.isnan(correlation):
                                    correlation = 0.0
                                    p_value = 1.0
                            else:
                                correlation = 0.0
                                p_value = 1.0

                            mse = np.mean((merged['P_fed'] - merged['P_baseline']) ** 2)
                            mae = np.mean(np.abs(merged['P_fed'] - merged['P_baseline']))

                            top_k = min(100, valid_count)
                            if top_k > 0:
                                fed_top = set(merged.nsmallest(top_k, 'P_fed')['SNP'].tolist())
                                baseline_top = set(merged.nsmallest(top_k, 'P_baseline')['SNP'].tolist())
                                top_overlap = len(fed_top & baseline_top) / top_k
                            else:
                                top_overlap = 0.0

                        results['local'] = {
                            'pearson_r': correlation,
                            'p_value': p_value,
                            'mse': mse,
                            'mae': mae,
                            'n_snps': len(merged),
                            'top_100_overlap': top_overlap
                        }

        # Global LR correlation
        fed_global_for_merge = self.federated_lr_global_deanon
        if fed_global_for_merge is not None and baseline_lr_for_merge is not None:
            if 'P' in fed_global_for_merge.columns and 'P' in baseline_lr_for_merge.columns:
                if 'SNP' in fed_global_for_merge.columns and 'SNP' in baseline_lr_for_merge.columns:
                    merged = pd.merge(
                        fed_global_for_merge[['SNP', 'P']],
                        baseline_lr_for_merge[['SNP', 'P']],
                        on='SNP',
                        suffixes=('_fed', '_baseline')
                    )

                    if len(merged) > 0:
                        initial_count = len(merged)
                        merged = merged.dropna(subset=['P_fed', 'P_baseline'])
                        merged = merged[(merged['P_fed'] > 0) & (merged['P_fed'] <= 1) &
                                        (merged['P_baseline'] > 0) & (merged['P_baseline'] <= 1)]
                        valid_count = len(merged)

                        if valid_count == 0:
                            correlation = 0.0
                            p_value = 1.0
                            mse = float('nan')
                            mae = float('nan')
                            top_overlap = 0.0
                        else:
                            merged['neg_log10_fed'] = -np.log10(merged['P_fed'] + 1e-300)
                            merged['neg_log10_baseline'] = -np.log10(merged['P_baseline'] + 1e-300)

                            fed_std = merged['neg_log10_fed'].std()
                            baseline_std = merged['neg_log10_baseline'].std()

                            if fed_std > 0 and baseline_std > 0:
                                correlation, p_value = pearsonr(merged['neg_log10_fed'], merged['neg_log10_baseline'])
                                if np.isnan(correlation):
                                    correlation = 0.0
                                    p_value = 1.0
                            else:
                                correlation = 0.0
                                p_value = 1.0

                            mse = np.mean((merged['P_fed'] - merged['P_baseline']) ** 2)
                            mae = np.mean(np.abs(merged['P_fed'] - merged['P_baseline']))

                            top_k = min(100, valid_count)
                            if top_k > 0:
                                fed_top = set(merged.nsmallest(top_k, 'P_fed')['SNP'].tolist())
                                baseline_top = set(merged.nsmallest(top_k, 'P_baseline')['SNP'].tolist())
                                top_overlap = len(fed_top & baseline_top) / top_k
                            else:
                                top_overlap = 0.0

                        results['global'] = {
                            'pearson_r': correlation,
                            'p_value': p_value,
                            'mse': mse,
                            'mae': mae,
                            'n_snps': len(merged),
                            'top_100_overlap': top_overlap
                        }

        # Global LR agreement and coverage
        snp_map = self._load_snp_maps(self.federated_dir)
        fed_agreement_list = []
        iterative_snp_set = set()

        if self.federated_lr_global_deanon is not None and 'is_significant' in self.federated_lr_global_deanon.columns:
            iterative_snps = self.federated_lr_global_deanon[['SNP', 'P', 'is_significant']].copy()

            def normalize_snp_id(snp_id: str) -> str:
                if snp_id.startswith('rs') and snp_id in snp_map:
                    return snp_map[snp_id]
                return snp_id

            iterative_snps['SNP'] = iterative_snps['SNP'].apply(normalize_snp_id)
            iterative_snps = iterative_snps.drop_duplicates(subset=['SNP'], keep='first')
            iterative_snp_set = set(iterative_snps['SNP'].tolist())
            self.global_iterative_snps = len(iterative_snps)
            fed_agreement_list.append(iterative_snps)
        else:
            self.logger.warning("Cannot add iterative_LR SNPs: missing is_significant in client-side results")

        if self.federated_lr_local is not None:
            local_lr_dedup = self.federated_lr_local.drop_duplicates(subset=['SNP'], keep='first')
            local_filtered = local_lr_dedup[local_lr_dedup['P'] >= p_threshold].copy()
            local_filtered = local_filtered[~local_filtered['SNP'].isin(iterative_snp_set)]
            if len(local_filtered) > 0:
                local_filtered['is_significant'] = False
                fed_agreement_list.append(local_filtered[['SNP', 'P', 'is_significant']])

        if fed_agreement_list:
            fed_for_agreement = pd.concat(fed_agreement_list, ignore_index=True)
            fed_for_agreement = fed_for_agreement.drop_duplicates(subset=['SNP'], keep='first')

            if baseline_lr_for_merge is not None:
                sig_agree = compute_significance_agreement(
                    fed_for_agreement,
                    baseline_lr_for_merge,
                    "Global LR Agreement (local_filtered + iterative_LR)",
                    fed_significance_col='is_significant'
                )
                if sig_agree:
                    results['significance_agreement'] = sig_agree

                    coverage_results = {}
                    # Include genome-wide significance and other thresholds
                    for threshold in [5e-8, 5e-5, 0.0005, 0.005, 0.05]:
                        baseline_sig = baseline_lr_for_merge[baseline_lr_for_merge['P'] < threshold]
                        if len(baseline_sig) > 0:
                            # Use both is_significant flag AND p-value threshold for federated SNPs
                            # This captures all SNPs that passed the screening threshold
                            fed_sig_by_flag = fed_for_agreement[fed_for_agreement['is_significant'] == True]
                            # Also include SNPs with p < screening threshold (p_threshold) that may not have is_significant=True
                            fed_sig_by_pvalue = fed_for_agreement[fed_for_agreement['P'] < p_threshold]
                            # Combine both sets
                            fed_sig = pd.concat([fed_sig_by_flag, fed_sig_by_pvalue]).drop_duplicates(subset=['SNP'], keep='first')
                            
                            baseline_sig_snps = set(baseline_sig['SNP'])
                            fed_sig_snps = set(fed_sig['SNP'])
                            covered = len(baseline_sig_snps & fed_sig_snps)
                            coverage = covered / len(baseline_sig_snps) if len(baseline_sig_snps) > 0 else 0.0
                            coverage_results[f'threshold_{threshold}'] = {
                                'coverage': coverage,
                                'covered': covered,
                                'total_baseline_sig': len(baseline_sig_snps),
                                'total_fed_sig': len(fed_sig_snps)
                            }
                    if coverage_results:
                        results['coverage'] = coverage_results

        # Screening summary (local vs global)
        total_snps = self.local_total_snps
        if total_snps is None and baseline_lr_for_merge is not None:
            total_snps = int(baseline_lr_for_merge['SNP'].nunique())
        screening = {}
        if self.local_retained_snps is not None and total_snps:
            screening['local'] = {
                'retained': self.local_retained_snps,
                'total': total_snps,
                'ratio': self.local_retained_snps / total_snps if total_snps > 0 else float('nan'),
            }
        if self.global_iterative_snps is not None and total_snps:
            screening['global'] = {
                'iterative': self.global_iterative_snps,
                'total': total_snps,
                'ratio': self.global_iterative_snps / total_snps if total_snps > 0 else float('nan'),
            }
        if screening:
            results['screening'] = screening

        self.analysis_results['lr_correlation'] = results
        return results

    def generate_manhattan_plot(self, output_file: Optional[str] = None) -> str:
        fed_lr_to_plot = None
        lr_type = "Global"
        if self.federated_lr_global_deanon is not None:
            fed_lr_to_plot = self.federated_lr_global_deanon
            lr_type = "Global (Deanon)"
        elif self.federated_lr_global is not None:
            fed_lr_to_plot = self.federated_lr_global
        elif self.federated_lr_local is not None:
            fed_lr_to_plot = self.federated_lr_local
            lr_type = "Local"

        baseline_lr_for_merge = self.baseline_lr_anon if self.baseline_lr_anon is not None else self.baseline_lr
        if fed_lr_to_plot is None or baseline_lr_for_merge is None:
            self.logger.warning("Cannot generate Manhattan plot: missing data")
            return ""

        fed_df = fed_lr_to_plot.copy()
        baseline_df = baseline_lr_for_merge.copy()

        if 'SNP' not in fed_df.columns:
            if fed_lr_to_plot.index.name == 'SNP':
                fed_df = fed_df.reset_index()
            else:
                self.logger.warning("SNP column not found in federated LR results")
                return ""
        if 'SNP' not in baseline_df.columns:
            if baseline_lr_for_merge.index.name == 'SNP':
                baseline_df = baseline_df.reset_index()
            else:
                self.logger.warning("SNP column not found in baseline LR results")
                return ""

        if ('CHR' not in fed_df.columns or 'BP' not in fed_df.columns) and ('CHR' in baseline_df.columns and 'BP' in baseline_df.columns):
            tmp = pd.merge(
                fed_df[['SNP', 'P']] if 'P' in fed_df.columns else fed_df[['SNP']],
                baseline_df[['SNP', 'CHR', 'BP']],
                on='SNP',
                how='inner'
            )
            if 'P' not in fed_df.columns:
                self.logger.warning("Cannot generate Manhattan plot: federated LR missing P column")
                return ""
            fed_df = tmp

        for col in ['CHR', 'BP', 'P']:
            if col not in fed_df.columns:
                self.logger.warning(f"Cannot generate Manhattan plot: federated missing {col} column")
                return ""
        for col in ['CHR', 'BP', 'P']:
            if col not in baseline_df.columns:
                self.logger.warning(f"Cannot generate Manhattan plot: baseline missing {col} column")
                return ""

        merged = pd.merge(
            fed_df[['SNP', 'CHR', 'BP', 'P']].rename(columns={'P': 'P_fed'}),
            baseline_df[['SNP', 'P']].rename(columns={'P': 'P_baseline'}),
            on='SNP'
        )
        if len(merged) == 0:
            self.logger.warning("No overlapping SNPs for Manhattan plot")
            return ""

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

        for chr_num in merged['CHR'].unique():
            chr_data = merged[merged['CHR'] == chr_num]
            ax1.scatter(chr_data['BP'], -np.log10(chr_data['P_fed'] + 1e-300), alpha=0.6, s=1)
        ax1.set_xlabel('Position')
        ax1.set_ylabel('-log10(P-value)')
        ax1.set_title(f'{lr_type} Federated GWAS Results')
        ax1.grid(True, alpha=0.3)

        for chr_num in merged['CHR'].unique():
            chr_data = merged[merged['CHR'] == chr_num]
            ax2.scatter(chr_data['BP'], -np.log10(chr_data['P_baseline'] + 1e-300), alpha=0.6, s=1)
        ax2.set_xlabel('Position')
        ax2.set_ylabel('-log10(P-value)')
        ax2.set_title('Centralized Baseline Results')
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        if output_file is None:
            output_file = str(self.federated_dir / 'manhattan_plot.png')
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        return output_file

    def generate_correlation_plot(self, output_file: Optional[str] = None) -> str:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        plot_generated = False

        baseline_lr_for_merge = self.baseline_lr_anon if self.baseline_lr_anon is not None else self.baseline_lr

        # Local LR plot
        if self.federated_lr_local is not None and baseline_lr_for_merge is not None:
            fed_local = self.federated_lr_local
            if self.federated_lr_local_filtered is not None:
                fed_local = self.federated_lr_local_filtered
            merged = pd.merge(
                fed_local[['SNP', 'P']],
                baseline_lr_for_merge[['SNP', 'P']],
                on='SNP',
                how='inner',
                suffixes=('_fed', '_baseline')
            )
            if len(merged) > 0:
                merged['neg_log10_fed'] = -np.log10(merged['P_fed'] + 1e-300)
                merged['neg_log10_baseline'] = -np.log10(merged['P_baseline'] + 1e-300)
                axes[0].scatter(merged['neg_log10_baseline'], merged['neg_log10_fed'], alpha=0.5, s=1)
                axes[0].set_xlabel('Baseline -log10(P)')
                axes[0].set_ylabel('Federated -log10(P)')
                axes[0].set_title('Local LR Correlation')
                axes[0].grid(True, alpha=0.3)
                plot_generated = True
            else:
                axes[0].text(0.5, 0.5, 'No overlapping SNPs', ha='center', va='center', transform=axes[0].transAxes)
        else:
            axes[0].text(0.5, 0.5, 'Local LR unavailable', ha='center', va='center', transform=axes[0].transAxes)

        # Global LR plot
        if self.federated_lr_global_deanon is not None and baseline_lr_for_merge is not None:
            merged = pd.merge(
                self.federated_lr_global_deanon[['SNP', 'P']],
                baseline_lr_for_merge[['SNP', 'P']],
                on='SNP',
                how='inner',
                suffixes=('_fed', '_baseline')
            )
            if len(merged) > 0:
                merged['neg_log10_fed'] = -np.log10(merged['P_fed'] + 1e-300)
                merged['neg_log10_baseline'] = -np.log10(merged['P_baseline'] + 1e-300)
                axes[1].scatter(merged['neg_log10_baseline'], merged['neg_log10_fed'], alpha=0.5, s=1)
                axes[1].set_xlabel('Baseline -log10(P)')
                axes[1].set_ylabel('Federated -log10(P)')
                axes[1].set_title('Global LR Correlation')
                axes[1].grid(True, alpha=0.3)
                plot_generated = True
            else:
                axes[1].text(0.5, 0.5, 'No overlapping SNPs', ha='center', va='center', transform=axes[1].transAxes)
        else:
            axes[1].text(0.5, 0.5, 'Global LR unavailable', ha='center', va='center', transform=axes[1].transAxes)

        if output_file is None:
            output_file = str(self.federated_dir / 'lr_correlation_plots.png')
        plt.tight_layout()
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        return output_file

    def write_report(self, output_path: Path) -> None:
        lr = self.analysis_results.get('lr_correlation', {})
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open('w') as f:
            f.write("# LR Evaluation Report\n\n")
            if not lr:
                f.write("LR analysis could not be computed (missing data).\n")
                return

            if 'local' in lr:
                f.write("## Local LR Correlation\n\n")
                f.write(f"- **Pearson r**: {lr['local']['pearson_r']:.4f}\n")
                f.write(f"- **P-value**: {lr['local']['p_value']:.4e}\n")
                f.write(f"- **MSE**: {lr['local']['mse']:.4f}\n")
                f.write(f"- **MAE**: {lr['local']['mae']:.4f}\n")
                f.write(f"- **n_snps**: {lr['local']['n_snps']}\n")
                f.write(f"- **Top-100 overlap**: {lr['local']['top_100_overlap']:.4f}\n\n")

            if 'global' in lr:
                f.write("## Global LR Correlation\n\n")
                f.write(f"- **Pearson r**: {lr['global']['pearson_r']:.4f}\n")
                f.write(f"- **P-value**: {lr['global']['p_value']:.4e}\n")
                f.write(f"- **MSE**: {lr['global']['mse']:.4f}\n")
                f.write(f"- **MAE**: {lr['global']['mae']:.4f}\n")
                f.write(f"- **n_snps**: {lr['global']['n_snps']}\n")
                f.write(f"- **Top-100 overlap**: {lr['global']['top_100_overlap']:.4f}\n\n")

            if 'significance_agreement' in lr:
                sa = lr['significance_agreement']
                f.write("## LR Significance Agreement\n\n")
                f.write(f"- **Agreement**: {sa['agreement']:.4f}\n")
                f.write(f"- **TP**: {sa['tp']}\n")
                f.write(f"- **TN**: {sa['tn']}\n")
                f.write(f"- **FP**: {sa['fp']}\n")
                f.write(f"- **FN**: {sa['fn']}\n")
                f.write(f"- **n_snps**: {sa['n_snps']}\n")
                f.write(f"- **Baseline threshold**: {sa['baseline_threshold']}\n")
                f.write(f"- **Federated threshold**: {sa['federated_threshold']}\n\n")

            if 'coverage' in lr:
                f.write("## Coverage Analysis\n\n")
                for key, stats in lr['coverage'].items():
                    threshold = key.replace('threshold_', '')
                    f.write(f"**Threshold {threshold}:**\n")
                    f.write(f"- Coverage: {stats['coverage']:.4f}\n")
                    f.write(f"- Covered: {stats['covered']}/{stats['total_baseline_sig']} baseline significant SNPs\n")
                    f.write(f"- Total federated positive: {stats['total_fed_sig']}\n\n")

            if 'screening' in lr:
                f.write("## Screening Summary\n\n")
                local = lr['screening'].get('local')
                if local:
                    f.write("**Local Stage:**\n")
                    f.write(f"- Retained: {local['retained']}/{local['total']} ({local['ratio']:.4f})\n\n")
                global_stage = lr['screening'].get('global')
                if global_stage:
                    f.write("**Global Stage (Iterative):**\n")
                    f.write(f"- Analyzed: {global_stage['iterative']}/{global_stage['total']} ({global_stage['ratio']:.4f})\n\n")

    def run(self, report_path: Optional[Path] = None, make_plots: bool = True) -> Dict[str, Any]:
        if not self.federated_dir.exists():
            raise FileNotFoundError(f"Federated results directory not found: {self.federated_dir}")

        self._load_lr_results(self.federated_dir, is_baseline=False)
        if self.baseline_dir is not None:
            self._load_lr_results(self.baseline_dir, is_baseline=True)

        results = self.analyze_lr_correlation()
        if make_plots:
            self.generate_manhattan_plot()
            self.generate_correlation_plot()
        if report_path is not None:
            self.write_report(report_path)
        return results


def main() -> None:
    parser = argparse.ArgumentParser(description='LR evaluation for federated GWAS results')
    parser.add_argument('results_dir', type=str, help='Federated results directory')
    parser.add_argument('--baseline', type=str, required=True, help='Baseline directory')
    parser.add_argument('--report', type=str, default=None, help='Output report path (default: results_dir/lr_report.md)')
    parser.add_argument('--no-plots', action='store_true', help='Skip plot generation')

    args = parser.parse_args()
    results_dir = Path(args.results_dir)
    report_path = Path(args.report) if args.report else results_dir / 'lr_report.md'

    evaluator = LREvaluator(str(results_dir), args.baseline)
    evaluator.run(report_path=report_path, make_plots=not args.no_plots)


if __name__ == '__main__':
    main()
