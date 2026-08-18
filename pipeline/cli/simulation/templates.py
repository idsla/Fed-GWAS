"""Template builders for generated FedGWAS simulation projects."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any

from pipeline.cli.simulation.presets import (
    EXPERIMENT_THRESHOLD_DEFAULTS,
    PARTICIPATION_DEFAULTS,
    THRESHOLD_DEFAULTS,
    Preset,
)


def project_pyproject(num_clients: int) -> str:
    """Build the minimal Flower `pyproject.toml` needed by `flwr run`.

    Args:
        num_clients: Number of local Flower supernodes/clients to configure.

    Returns:
        TOML string written to the generated project `pyproject.toml`.
    """
    return f"""[project]
name = "fedgwas-simulation-project"
version = "0.1.0"
requires-python = ">=3.11"

[tool.flwr.app]
publisher = "fedgwas"

[tool.flwr.app.components]
serverapp = "pipeline.server_app:app"
clientapp = "pipeline.client_app:app"

[tool.flwr.app.config]
simulation = true
num-server-rounds = 100
config_path = "configs"
phenotype_fix_missing_to_case = false

[tool.flwr.federations]
default = "local-simulation"

[tool.flwr.federations.local-simulation]
options.num-supernodes = {num_clients}
"""


def base_project_config(project_name: str = "fedgwas-simulation-project") -> dict[str, Any]:
    """Build the minimal `fedgwas.yaml` created by `fedgwas-sim init`.

    Args:
        project_name: Human-readable project name stored in the project file.

    Returns:
        Mapping written to `fedgwas.yaml` for an initialized but unconfigured
        simulation project.
    """
    return {
        "mode": "simulation",
        "project_state": "initialized",
        "project_name": project_name,
        "config_dir": "configs",
        "data_dir": "data",
        "results_dir": "results",
        "logs_dir": "logs",
        "plink": "auto",
    }


def base_experiment_config(project_name: str = "fedgwas-simulation-project") -> dict[str, Any]:
    """Build the placeholder `config.yaml` created by `fedgwas-sim init`.

    Args:
        project_name: Human-readable project name stored in the config.

    Returns:
        Minimal experiment config mapping. It is intentionally incomplete until
        `setup-experiment` or `init --from-example` configures a runnable study.
    """
    return {
        "description": "Initialized FedGWAS simulation project.",
        "experiment_name": project_name.replace("-", "_"),
        "experiment_category": "simulation",
        "num_clients": 0,
        "scenario": "initialized",
        "clients": {"config_files": {}, "participation": {}, "thresholds": {}},
        "data": {"data_dir": "data", "partition_strategy": None, "scale": None},
        "server": {"chunk_size": None, "min_available_clients": 1, "min_fit_clients": 1, "num_server_rounds": 0},
        "analysis": {"generate_baseline": False, "compare_results": False, "metrics_collection": False},
    }


def project_config(preset: Preset, seed: int | None) -> dict[str, Any]:
    """Build `fedgwas.yaml`, the CLI-owned project-level configuration.

    Args:
        preset: Preset selected by `setup-experiment`.
        seed: Optional synthetic data seed to persist for reproducibility.

    Returns:
        Mapping written to `fedgwas.yaml`.
    """
    payload: dict[str, Any] = {
        "mode": "simulation",
        "project_state": "configured",
        "preset": preset.name,
        "data_source": preset.source,
        "num_clients": preset.num_clients,
        "config_dir": "configs",
        "data_dir": "data",
        "results_dir": "results",
        "logs_dir": "logs",
        "plink": "auto",
    }
    if preset.scale:
        payload["scale"] = preset.scale
    if seed is not None:
        payload["seed"] = seed
    return payload


def experiment_config(preset: Preset) -> dict[str, Any]:
    """Build `config.yaml`, shared by monitoring, retention, baseline, and evaluators.

    Args:
        preset: Preset selected by `setup-experiment`.

    Returns:
        Mapping written to `config.yaml`.
    """
    config_files = {
        client_id - 1: f"configs/center_{client_id}.yaml"
        for client_id in range(1, preset.num_clients + 1)
    }
    return {
        "description": preset.description,
        "experiment_name": preset.experiment_name or preset.name.replace("-", "_"),
        "experiment_category": preset.experiment_category,
        "num_clients": preset.num_clients,
        "scenario": preset.scenario or preset.name,
        "clients": {
            "config_files": config_files,
            "participation": dict(PARTICIPATION_DEFAULTS),
            "thresholds": dict(preset.client_thresholds or EXPERIMENT_THRESHOLD_DEFAULTS),
        },
        "data": {
            "data_dir": "data",
            "partition_strategy": preset.partition_strategy,
            "scale": preset.scale or preset.source,
        },
        "monitoring": {
            "enable_network_monitoring": False,
            "enable_performance_monitoring": False,
            "network_monitoring_interval": 5.0,
            "resource_monitoring_interval": 2.0,
            "network_max_samples": 500,
        },
        "retention": {"tier": "standard", "auto_apply_on_complete": True},
        "server": {
            "chunk_size": preset.server_chunk_size or preset.chunk_size,
            "min_available_clients": 1,
            "min_fit_clients": 1,
            "num_server_rounds": preset.rounds,
        },
        "analysis": dict(preset.analysis),
    }


def center_config(preset: Preset, center_id: int) -> dict[str, Any]:
    """Build one center's pipeline config with project-relative paths.

    Args:
        preset: Preset selected by `setup-experiment`.
        center_id: One-based center id.

    Returns:
        Mapping written to `configs/center_<id>.yaml`.
    """
    prefix_name = (
        f"{preset.scale}_center_{center_id}"
        if preset.scale
        else f"study_center_{center_id}"
    )
    return {
        "input_data": {
            "path": f"data/center_{center_id}/{prefix_name}",
            "type": "bed",
        },
        "output": {
            "intermediate_dir": f"results/center_{center_id}/intermediate",
            "log_dir": f"results/center_{center_id}/logs",
        },
        "parameters": {
            "sample_offset": 1000000000000,
            "chunk_size": preset.chunk_size,
            "sample_chunk_size": preset.sample_chunk_size,
            "snp_chunk_size": preset.snp_chunk_size,
            "lr_target_chunks": preset.lr_target_chunks,
            "run_id": datetime.now().strftime("%Y%m%d_%H%M%S"),
        },
        "thresholds": dict(preset.center_thresholds or THRESHOLD_DEFAULTS),
        "flower": {"server_address": "127.0.0.1:9092", "num_rounds": preset.rounds},
        "participation": dict(PARTICIPATION_DEFAULTS),
    }


def server_config() -> dict[str, Any]:
    """Build the server output config consumed by `pipeline.server_app`.

    Returns:
        Mapping written to `configs/server.yaml`.
    """
    return {
        "output": {
            "intermediate_dir": "results/server/intermediate",
            "log_dir": "results/server/logs",
        }
    }


def write_real_preset_scripts(project_dir: Path, preset: Preset) -> None:
    """Write placeholder preparation scripts for public-data presets.

    Real genotype datasets can be large and may have access constraints. Setup
    therefore creates a reproducible template, while `prepare-data --download`
    remains the explicit opt-in point for future dataset-specific preparation.

    Args:
        project_dir: Simulation project root.
        preset: Real-data preset selected by `setup-experiment`.

    Returns:
        None. Writes files under `scripts/`.
    """
    scripts_dir = project_dir / "scripts"
    scripts_dir.mkdir(parents=True, exist_ok=True)
    script_path = scripts_dir / "prepare_data.py"
    lines = [
        '"""Prepare public genotype data for this FedGWAS simulation project."""',
        "",
        "from __future__ import annotations",
        "",
        "import argparse",
        "from pathlib import Path",
        "",
        "",
        "def main() -> None:",
        "    parser = argparse.ArgumentParser()",
        "    parser.add_argument('--download', action='store_true')",
        "    parser.add_argument('--seed', type=int, default=42)",
        "    args = parser.parse_args()",
        "    if not args.download:",
        "        print('Template created. Re-run with --download after reviewing dataset access requirements.')",
        "        return",
    ]
    if preset.name == "1000genomes":
        lines.extend(
            [
                "    from pipeline.cli.simulation.real_data import prepare_1000genomes_chr22",
                "",
                "    prepare_1000genomes_chr22(Path.cwd(), seed=args.seed)",
            ]
        )
    else:
        lines.append(f"    raise SystemExit('{preset.name} download/normalization is template-only in this release.')")
    lines.extend(["", "", "if __name__ == '__main__':", "    main()", ""])
    script_path.write_text("\n".join(lines), encoding="utf-8")
    (scripts_dir / "README.md").write_text(
        "\n".join(
            [
                f"# {preset.name} Data Preparation",
                "",
                "This template intentionally does not download public genotype data during setup.",
                "Review access, licensing, filtering, phenotype simulation, and partitioning choices before running preparation.",
                "",
                "Run:",
                "",
                "```bash",
                "python scripts/prepare_data.py --download",
                "```",
                "",
            ]
        ),
        encoding="utf-8",
    )
