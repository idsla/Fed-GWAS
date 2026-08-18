"""Apply product retention tiers to a federated run results directory."""

from __future__ import annotations

import json
import shutil
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

# PLINK / GWAS artifacts worth keeping when keep_plink_science_outputs is true
PLINK_SCIENCE_SUFFIXES = (
    ".assoc",
    ".assoc.logistic",
    ".assoc.linear",
    ".kin0",
    ".imiss",
    ".lmiss",
    ".frq",
    ".hwe",
    ".genome",
    ".nosex",
)

PLINK_SCIENCE_PREFIXES = (
    "king_filtered",
    "global_filtered",
    "filtered_data",
    "qc_",
    "king_results",
    "lr_results",
    "local_lr",
    "global_qc",
    "lr_screening",
    "lr_pvals",
    "lr_round",
)

PLINK_SCIENCE_EXACT = {
    "global_qc_exclusion.txt",
    "lr_screening_significant.txt",
}

REPORT_NAMES = {
    "qc_report.md",
    "lr_report.md",
    "king_report.md",
    "evaluation_report.md",
    "experiment_report.md",
    "run_summary.md",
    "run_summary.json",
}

CRYPTO_NAMES = {
    "ecc_private_key.pem",
    "ecc_public_key.pem",
    "global_seed.txt",
}

PERFORMANCE_NODE_FILES = {
    "stage_metrics.csv",
    "client_metrics.csv",
    "performance_summary.md",
}

PERFORMANCE_ROOT_FILES = {
    "stage_metrics.csv",
    "client_metrics.csv",
    "run_summary.md",
    "run_summary.json",
}


@dataclass
class RetentionSummary:
    tier: str
    results_root: str
    deleted_paths: List[str] = field(default_factory=list)
    kept_paths: List[str] = field(default_factory=list)
    skipped_paths: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict:
        return {
            "tier": self.tier,
            "results_root": self.results_root,
            "deleted_count": len(self.deleted_paths),
            "deleted_paths": self.deleted_paths,
            "kept_count": len(self.kept_paths),
            "skipped_count": len(self.skipped_paths),
        }


def _is_plink_science_file(path: Path) -> bool:
    name = path.name
    if name in REPORT_NAMES or name in PLINK_SCIENCE_EXACT:
        return True
    if name.endswith(".json") and (
        name.startswith("lr_") or name.startswith("king_") or "qc" in name
    ):
        return True
    for suffix in PLINK_SCIENCE_SUFFIXES:
        if name.endswith(suffix):
            return True
    for prefix in PLINK_SCIENCE_PREFIXES:
        if name.startswith(prefix):
            return True
    if name.endswith((".bed", ".bim", ".fam")) and "chunk_" not in name:
        # filtered outputs in logs; exclude chunk_* intermediates living under logs
        if "filtered" in name or "king_" in name or "global_" in name:
            return True
    return False


def _delete_path(path: Path, *, dry_run: bool, summary: RetentionSummary) -> None:
    rel = str(path)
    if dry_run:
        summary.deleted_paths.append(rel)
        return
    if path.is_dir():
        shutil.rmtree(path, ignore_errors=True)
    elif path.is_file():
        path.unlink(missing_ok=True)
    summary.deleted_paths.append(rel)


def _should_remove_file(
    path: Path,
    *,
    results_root: Path,
    settings: Dict,
    node_logs: bool,
) -> bool:
    name = path.name

    # Ephemeral PLINK / pipeline scratch (never product deliverables)
    if name.startswith("temp_"):
        return True
    if name.startswith("id_map") or name.startswith("metadata"):
        return True
    if name == "performance_monitor.log":
        return not settings.get("keep_performance_logs")
    if name.endswith(".log"):
        if name == "server_log.txt" or name.endswith("_log.txt"):
            return not settings.get("keep_text_logs")
        return True

    if settings.get("keep_evaluation_reports") and name in REPORT_NAMES:
        return False

    if node_logs and settings.get("keep_text_logs") and name.endswith("_log.txt"):
        return False

    if settings.get("keep_plink_science_outputs") and _is_plink_science_file(path):
        return False

    if (
        not node_logs
        and settings.get("keep_merged_performance_csv")
        and name in PERFORMANCE_ROOT_FILES
    ):
        return False

    if (
        node_logs
        and settings.get("keep_node_performance_csv")
        and name in PERFORMANCE_NODE_FILES
    ):
        return False

    if settings.get("keep_performance_logs") and name == "performance_monitor.log":
        return False

    if settings.get("keep_crypto_state") and name in CRYPTO_NAMES:
        return False

    if settings.get("keep_king_accumulator") and name == "king_accumulator_state.pkl":
        return False

    if settings.get("keep_network_monitor") and path.parent.name.endswith("_network"):
        return False

    return True


def apply_retention_policy(
    results_root: Path,
    settings: Dict,
    *,
    dry_run: bool = False,
) -> RetentionSummary:
    """
    Prune run outputs under results_root (server/, center_*/).

    settings: output of resolve_retention_settings().
    """
    results_root = results_root.resolve()
    summary = RetentionSummary(tier=settings.get("tier", "standard"), results_root=str(results_root))

    if not results_root.exists():
        summary.skipped_paths.append("results_root missing")
        return summary

    # Whole intermediate trees
    if not settings.get("keep_intermediate"):
        for node in sorted(results_root.glob("*/intermediate")):
            if node.is_dir():
                _delete_path(node, dry_run=dry_run, summary=summary)

    # Network monitor subdirs (e.g. logs/tiny_even_network/)
    if not settings.get("keep_network_monitor"):
        for net_dir in results_root.rglob("*_network"):
            if net_dir.is_dir():
                _delete_path(net_dir, dry_run=dry_run, summary=summary)

    # Per-node logs and server logs
    log_dirs = list(results_root.glob("*/logs"))
    for logs in log_dirs:
        if not logs.is_dir():
            continue
        for item in sorted(logs.iterdir()):
            if item.is_dir():
                if not settings.get("keep_network_monitor") and item.name.endswith("_network"):
                    _delete_path(item, dry_run=dry_run, summary=summary)
                continue
            if _should_remove_file(item, results_root=results_root, settings=settings, node_logs=True):
                _delete_path(item, dry_run=dry_run, summary=summary)
            else:
                summary.kept_paths.append(str(item.relative_to(results_root)))

    # Results root files
    for item in results_root.iterdir():
        if item.is_dir():
            continue
        if _should_remove_file(item, results_root=results_root, settings=settings, node_logs=False):
            _delete_path(item, dry_run=dry_run, summary=summary)
        else:
            summary.kept_paths.append(item.name)

    # Write manifest
    manifest = results_root / "retention_manifest.json"
    manifest_data = {
        **summary.to_dict(),
        "settings": {k: v for k, v in settings.items() if k.startswith("keep_") or k == "tier"},
        "applied_at": datetime.now().isoformat(),
        "dry_run": dry_run,
    }
    if not dry_run:
        manifest.write_text(json.dumps(manifest_data, indent=2))

    return summary


def maybe_apply_retention_on_complete(
    results_root: Path,
    settings: Dict,
    *,
    dry_run: bool = False,
    logger=None,
) -> Optional[RetentionSummary]:
    if not settings.get("auto_apply_on_complete"):
        return None
    summary = apply_retention_policy(results_root, settings, dry_run=dry_run)
    if logger:
        logger.info(
            "[Retention] tier=%s deleted=%d kept=%d (dry_run=%s)",
            summary.tier,
            len(summary.deleted_paths),
            len(summary.kept_paths),
            dry_run,
        )
    return summary
