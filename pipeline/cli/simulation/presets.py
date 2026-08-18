"""Preset registry and shared defaults for simulation projects."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


class PresetChoice(str, Enum):
    """Typer-facing preset choices shown in CLI help and validation errors."""

    SYN_TINY = "syn-tiny"
    SYN_SMALL = "syn-small"
    SYN_MEDIUM = "syn-medium"
    GENOMES = "1000genomes"
    HAPMAP = "hapmap"


@dataclass(frozen=True)
class Preset:
    """A supported experiment preset and the defaults used to scaffold it.

    Attributes:
        name: Public preset name accepted by `fedgwas-sim setup-experiment`.
        source: Logical data source label, such as `synthetic` or `1000genomes`.
        scale: Synthetic generator scale name. `None` means the preset is a
            real-data template rather than an immediately generated dataset.
        num_clients: Number of Flower simulation clients/centers to scaffold.
        rounds: Suggested default Flower server rounds for this preset.
        chunk_size: Default server/client chunk size used in generated configs.
        sample_chunk_size: Default per-client sample chunk size.
        snp_chunk_size: Default per-client SNP chunk size.
        lr_target_chunks: Default number of LR chunks targeted by client config.
        synthetic: Whether setup should generate PLINK data immediately.
    description: Human-readable preset description written into config.yaml.
        experiment_name: Optional canonical experiment name written into
            `config.yaml`. Defaults to the preset name with hyphens replaced.
        experiment_category: Experiment category used by reporting and parity
            checks, such as `correctness` or `performance`.
        scenario: Optional scenario label written into `config.yaml`.
        partition_strategy: Data partition strategy written into `config.yaml`.
        server_chunk_size: Optional server chunk size. When omitted, uses
            `chunk_size`.
        analysis: Analysis flags written into `config.yaml`.
        client_thresholds: Thresholds written under `clients.thresholds`.
        center_thresholds: Thresholds written into center configs.
    """

    name: str
    source: str
    scale: Optional[str]
    num_clients: int
    rounds: int
    chunk_size: int
    sample_chunk_size: int
    snp_chunk_size: int
    lr_target_chunks: int
    synthetic: bool
    description: str
    experiment_name: Optional[str] = None
    experiment_category: str = "simulation"
    scenario: Optional[str] = None
    partition_strategy: str = "even"
    server_chunk_size: Optional[int] = None
    analysis: dict[str, bool] = field(
        default_factory=lambda: {
            "generate_baseline": True,
            "compare_results": True,
            "metrics_collection": True,
        }
    )
    client_thresholds: Optional[dict[str, float]] = None
    center_thresholds: Optional[dict[str, float]] = None


# Public preset registry. Keys are CLI-facing names; values carry every default
# needed to create `fedgwas.yaml`, `config.yaml`, Flower config, and center
# configs without consulting repository-only experiment paths.
PRESETS: dict[str, Preset] = {
    "syn-tiny": Preset(
        name="syn-tiny",
        source="synthetic",
        scale="tiny",
        num_clients=2,
        rounds=50,
        chunk_size=200,
        sample_chunk_size=100,
        snp_chunk_size=200,
        lr_target_chunks=2,
        synthetic=True,
        description="Correctness validation: Tiny scale with even partition (2 clients)",
        experiment_name="tiny_even",
        experiment_category="correctness",
        scenario="correctness_tiny",
    ),
    "syn-small": Preset(
        name="syn-small",
        source="synthetic",
        scale="small",
        num_clients=2,
        rounds=80,
        chunk_size=500,
        sample_chunk_size=200,
        snp_chunk_size=500,
        lr_target_chunks=4,
        synthetic=True,
        description="Performance benchmark: Small scale with even partition (2 clients, 2000×20k)",
        experiment_name="small_even",
        experiment_category="performance",
        scenario="performance_small",
        analysis={
            "generate_baseline": False,
            "compare_results": False,
            "metrics_collection": True,
        },
    ),
    "syn-medium": Preset(
        name="syn-medium",
        source="synthetic",
        scale="medium",
        num_clients=2,
        rounds=130,
        chunk_size=1000,
        sample_chunk_size=500,
        snp_chunk_size=500,
        lr_target_chunks=5,
        synthetic=True,
        description="Performance benchmark: Medium scale with even partition (2 clients, 5000×50k)",
        experiment_name="medium_even",
        experiment_category="performance",
        scenario="performance_medium",
        server_chunk_size=500,
        analysis={
            "generate_baseline": False,
            "compare_results": False,
            "metrics_collection": True,
        },
    ),
    "1000genomes": Preset(
        name="1000genomes",
        source="1000genomes",
        scale=None,
        num_clients=2,
        rounds=350,
        chunk_size=500,
        sample_chunk_size=200,
        snp_chunk_size=500,
        lr_target_chunks=4,
        synthetic=False,
        description="Template for a chromosome 22 1000 Genomes simulation.",
    ),
    "hapmap": Preset(
        name="hapmap",
        source="hapmap",
        scale=None,
        num_clients=2,
        rounds=200,
        chunk_size=500,
        sample_chunk_size=200,
        snp_chunk_size=500,
        lr_target_chunks=4,
        synthetic=False,
        description="Template for a HapMap simulation.",
    ),
}

# Default pipeline stage participation for generated center configs. These flags
# mirror the existing tiny/small/medium experiments and keep the full QC, KING,
# and LR workflow enabled unless users edit YAML afterward.
PARTICIPATION_DEFAULTS: dict[str, bool] = {
    "key_exchange": True,
    "sync": True,
    "local_qc": True,
    "global_qc": True,
    "global_qc_response": True,
    "init_chunks": True,
    "iterative_king": True,
    "local_lr": True,
    "local_lr_filter_response": True,
    "init_chunks_lr": True,
    "iterative_lr": True,
}

# Default QC/LR/KING thresholds written into generated center and experiment
# configs. Baseline generation also reads these values to keep comparisons
# aligned with the federated run.
THRESHOLD_DEFAULTS: dict[str, float] = {
    "maf_threshold": 0.01,
    "missing_threshold": 0.05,
    "hwe_threshold": 1e-6,
    "p_threshold": 0.3,
    "local_lr_threshold": 0.3,
    "global_lr_threshold": 0.1,
    "king_threshold": 0.23,
}

# Threshold subset used by experiment-level config files. Per-center configs
# include the LR-specific thresholds, but the committed experiment config files
# only expose the common QC/KING/LR p-value keys under `clients.thresholds`.
EXPERIMENT_THRESHOLD_DEFAULTS: dict[str, float] = {
    "maf_threshold": 0.01,
    "missing_threshold": 0.05,
    "hwe_threshold": 1e-6,
    "p_threshold": 0.3,
    "king_threshold": 0.23,
}
