"""Built-in example experiment templates for the simulation CLI."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any

from pipeline.cli.simulation.presets import PARTICIPATION_DEFAULTS, THRESHOLD_DEFAULTS


class ExampleChoice(str, Enum):
    """Typer-facing choices accepted by `fedgwas-sim init --from-example`."""

    SYN_TINY = "syn-tiny"
    SYN_SMALL = "syn-small"
    SYN_MEDIUM = "syn-medium"
    GENOMES = "1000genomes"


@dataclass(frozen=True)
class ExampleSpec:
    """Normalized example metadata used to create package-local studies.

    Attributes:
        name: CLI-facing example name.
        preset: Optional preset name related to this example.
        experiment_name: Canonical experiment name from the repository example.
        experiment_category: Experiment category from the repository example.
        description: Human-readable scenario description.
        scenario: Scenario label from the repository example.
        scale: Data scale label when available.
        partition_strategy: Data partition strategy.
        num_clients: Number of local simulation clients.
        server_rounds: Number of Flower server rounds.
        server_chunk_size: Chunk size in top-level server config.
        center_chunk_size: Chunk size in center configs.
        sample_chunk_size: Sample chunk size in center configs.
        snp_chunk_size: SNP chunk size in center configs.
        lr_target_chunks: LR target chunks in center configs.
        top_thresholds: Thresholds written under `clients.thresholds`.
        center_thresholds: Thresholds written into center configs.
        analysis: Analysis flags from the repository example.
        data_prefix: PLINK bfile prefix basename for each center.
    """

    name: str
    preset: str | None
    experiment_name: str
    experiment_category: str
    description: str
    scenario: str | None
    scale: str | None
    partition_strategy: str
    num_clients: int
    server_rounds: int
    server_chunk_size: int
    center_chunk_size: int
    sample_chunk_size: int
    snp_chunk_size: int
    lr_target_chunks: int
    top_thresholds: dict[str, float] | None
    center_thresholds: dict[str, float]
    analysis: dict[str, bool]
    data_prefix: str


TOP_THRESHOLDS: dict[str, float] = {
    "hwe_threshold": 1e-6,
    "maf_threshold": 0.01,
    "missing_threshold": 0.05,
    "p_threshold": 0.3,
    "king_threshold": 0.23,
}

REAL_WORLD_THRESHOLDS: dict[str, float] = {
    "maf_threshold": 0.01,
    "missing_threshold": 0.05,
    "hwe_threshold": 1e-6,
    "p_threshold": 0.1,
    "king_threshold": 0.23,
}

EXAMPLES: dict[str, ExampleSpec] = {
    "syn-tiny": ExampleSpec(
        name="syn-tiny",
        preset="syn-tiny",
        experiment_name="tiny_even",
        experiment_category="correctness",
        description="Correctness validation: Tiny scale with even partition (2 clients)",
        scenario="correctness_tiny",
        scale="tiny",
        partition_strategy="even",
        num_clients=2,
        server_rounds=50,
        server_chunk_size=200,
        center_chunk_size=200,
        sample_chunk_size=100,
        snp_chunk_size=200,
        lr_target_chunks=2,
        top_thresholds=dict(TOP_THRESHOLDS),
        center_thresholds=dict(THRESHOLD_DEFAULTS),
        analysis={"generate_baseline": True, "compare_results": True, "metrics_collection": True},
        data_prefix="tiny_center_{center_id}",
    ),
    "syn-small": ExampleSpec(
        name="syn-small",
        preset="syn-small",
        experiment_name="small_even",
        experiment_category="performance",
        description="Performance benchmark: Small scale with even partition (2 clients, 2000×20k)",
        scenario="performance_small",
        scale="small",
        partition_strategy="even",
        num_clients=2,
        server_rounds=80,
        server_chunk_size=500,
        center_chunk_size=500,
        sample_chunk_size=200,
        snp_chunk_size=500,
        lr_target_chunks=4,
        top_thresholds=dict(TOP_THRESHOLDS),
        center_thresholds=dict(THRESHOLD_DEFAULTS),
        analysis={"generate_baseline": False, "compare_results": False, "metrics_collection": True},
        data_prefix="small_center_{center_id}",
    ),
    "syn-medium": ExampleSpec(
        name="syn-medium",
        preset="syn-medium",
        experiment_name="medium_even",
        experiment_category="performance",
        description="Performance benchmark: Medium scale with even partition (2 clients, 5000×50k)",
        scenario="performance_medium",
        scale="medium",
        partition_strategy="even",
        num_clients=2,
        server_rounds=130,
        server_chunk_size=500,
        center_chunk_size=1000,
        sample_chunk_size=500,
        snp_chunk_size=500,
        lr_target_chunks=5,
        top_thresholds=dict(TOP_THRESHOLDS),
        center_thresholds=dict(THRESHOLD_DEFAULTS),
        analysis={"generate_baseline": False, "compare_results": False, "metrics_collection": True},
        data_prefix="medium_center_{center_id}",
    ),
    "1000genomes": ExampleSpec(
        name="1000genomes",
        preset="1000genomes",
        experiment_name="1000genomes",
        experiment_category="real_world",
        description="1000 Genomes Project (chromosome 22) with simulated phenotypes (additive model)",
        scenario=None,
        scale=None,
        partition_strategy="population",
        num_clients=2,
        server_rounds=350,
        server_chunk_size=100,
        center_chunk_size=200,
        sample_chunk_size=200,
        snp_chunk_size=5000,
        lr_target_chunks=6,
        top_thresholds=None,
        center_thresholds=dict(REAL_WORLD_THRESHOLDS),
        analysis={"generate_baseline": True, "compare_results": True, "metrics_collection": True},
        data_prefix="genotypes",
    ),
}


def example_project_config(spec: ExampleSpec) -> dict[str, Any]:
    """Build `fedgwas.yaml` for a built-in example project.

    Args:
        spec: Built-in example specification.

    Returns:
        Project-level simulation config mapping.
    """
    payload: dict[str, Any] = {
        "mode": "simulation",
        "project_state": "configured",
        "example": spec.name,
        "preset": spec.preset,
        "num_clients": spec.num_clients,
        "config_dir": "configs",
        "data_dir": "data",
        "results_dir": "results",
        "logs_dir": "logs",
        "plink": "auto",
    }
    if spec.scale:
        payload["scale"] = spec.scale
    return payload


def example_experiment_config(spec: ExampleSpec) -> dict[str, Any]:
    """Build normalized `config.yaml` for a built-in example.

    Args:
        spec: Built-in example specification.

    Returns:
        Experiment config mapping with project-relative paths.
    """
    clients: dict[str, Any] = {
        "config_files": {center_id - 1: f"configs/center_{center_id}.yaml" for center_id in range(1, spec.num_clients + 1)}
    }
    if spec.top_thresholds is not None:
        clients["participation"] = dict(PARTICIPATION_DEFAULTS)
        clients["thresholds"] = dict(spec.top_thresholds)

    payload: dict[str, Any] = {
        "experiment_name": spec.experiment_name,
        "experiment_category": spec.experiment_category,
        "description": spec.description,
        "clients": clients,
        "data": {
            "data_dir": "data",
            "partition_strategy": spec.partition_strategy,
        },
        "server": {
            "num_server_rounds": spec.server_rounds,
            "chunk_size": spec.server_chunk_size,
            "min_available_clients": 1,
            "min_fit_clients": 1,
        },
        "analysis": dict(spec.analysis),
    }
    if spec.scale:
        payload["data"]["scale"] = spec.scale
    if spec.scenario:
        payload["scenario"] = spec.scenario
    if spec.top_thresholds is not None:
        payload["num_clients"] = spec.num_clients
    return payload


def example_center_config(spec: ExampleSpec, center_id: int) -> dict[str, Any]:
    """Build one normalized center config for a built-in example.

    Args:
        spec: Built-in example specification.
        center_id: One-based center id.

    Returns:
        Center config mapping with project-relative paths.
    """
    return {
        "input_data": {
            "path": f"data/center_{center_id}/{spec.data_prefix.format(center_id=center_id)}",
            "type": "bed",
        },
        "output": {
            "intermediate_dir": f"results/center_{center_id}/intermediate",
            "log_dir": f"results/center_{center_id}/logs",
        },
        "parameters": {
            "sample_offset": 1000000000000,
            "chunk_size": spec.center_chunk_size,
            "sample_chunk_size": spec.sample_chunk_size,
            "snp_chunk_size": spec.snp_chunk_size,
            "lr_target_chunks": spec.lr_target_chunks,
            "run_id": "example",
        },
        "thresholds": dict(spec.center_thresholds),
        "flower": {"server_address": "127.0.0.1:9092", "num_rounds": spec.server_rounds},
        "participation": dict(PARTICIPATION_DEFAULTS),
    }
