"""Typer command handlers for the installed FedGWAS simulation CLI.

The public entry point remains `pipeline.cli.sim:app`. Implementation details
live under `pipeline.cli.simulation` so this module only defines command wiring,
CLI-facing docstrings, and argument/exit behavior.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import Annotated, Optional

import typer

from pipeline.cli.simulation import examples, outputs, paths, setup as setup_mod, summary, validation
from pipeline.cli.simulation.presets import PRESETS, PresetChoice
from pipeline.evaluation import evaluate_all as evaluation_tool
from pipeline.tools import generate_baseline as baseline_tool


#===========================================================================================================
# CLI app wiring
#===========================================================================================================

ROOT_EXAMPLES = """Examples:\n
  mkdir tiny-study\n
  cd tiny-study\n
  fedgwas-sim init\n
  fedgwas-sim setup-experiment syn-tiny\n
  fedgwas-sim check\n
  fedgwas-sim run --rounds 50\n
  fedgwas-sim results collect --label tiny_run\n
"""

app = typer.Typer(
    help="Manage and run local FedGWAS Flower simulations.",
    epilog=ROOT_EXAMPLES,
    no_args_is_help=True,
)
data_app = typer.Typer(help="Configure local center data.", no_args_is_help=True)
baseline_app = typer.Typer(help="Generate centralized comparison baselines.", no_args_is_help=True)
results_app = typer.Typer(help="Collect run outputs and metrics.", no_args_is_help=True)
summarize_app = typer.Typer(help="Summarize configured experiments and data.", no_args_is_help=True)

app.add_typer(data_app, name="data")
app.add_typer(baseline_app, name="baseline")
app.add_typer(results_app, name="results")
app.add_typer(summarize_app, name="summarize")


#===========================================================================================================
# CLI: init
#===========================================================================================================

INIT_HELP = (
    "Initialize a FedGWAS simulation project in the current directory.\n\n"
    "Use --from-example to create one of the built-in experiment layouts."
)

INIT_EXAMPLES = """Examples:\n
  fedgwas-sim init\n
  fedgwas-sim init --from-example syn-tiny\n
  fedgwas-sim init --from-example syn-tiny --seed 42\n
  fedgwas-sim init --from-example syn-small\n
  fedgwas-sim init --from-example syn-medium\n
  fedgwas-sim init --from-example 1000genomes\n
  fedgwas-sim init --from-example 1000genomes --no-prepare-data\n
"""


@app.command(
    "init",
    help=INIT_HELP,
    epilog=INIT_EXAMPLES,
)
def init(
    from_example: Annotated[
        Optional[examples.ExampleChoice],
        typer.Option("--from-example", help="Initialize from a built-in experiment example."),
    ] = None,
    force: Annotated[bool, typer.Option("--force", help="Allow writing into an existing directory.")] = False,
    prepare_data: Annotated[
        bool,
        typer.Option("--prepare-data/--no-prepare-data", help="Prepare data for built-in examples."),
    ] = True,
    seed: Annotated[Optional[int], typer.Option("--seed", help="Synthetic example data random seed.")] = None,
) -> None:
    """Initialize a simulation project directory.

    CLI usage:
        `fedgwas-sim init`
        `fedgwas-sim init --from-example syn-tiny`

    Args:
        from_example: Optional built-in example to materialize into the current
            directory. When omitted, only the base project structure is created.
        force: Whether to allow writing into a non-empty directory that is not
            already marked as a FedGWAS simulation project.
        prepare_data: Whether `--from-example` should prepare the configured
            data immediately after project files are written.
        seed: Optional random seed forwarded to synthetic example generation.

    Returns:
        None. Writes project files and prints the initialized path.

    Raises:
        typer.BadParameter: If the directory is unsafe to initialize or the
            example name is unknown.
    """
    project_dir = Path.cwd().resolve()
    if from_example is None:
        setup_mod.init_project(project_dir, force=force)
        typer.echo(f"Initialized simulation project: {project_dir}")
        return
    setup_mod.create_project_from_example(
        project_dir,
        from_example.value,
        force=force,
        prepare_data=prepare_data,
        seed=seed,
    )
    typer.echo(f"Initialized simulation project from example '{from_example.value}': {project_dir}")


#===========================================================================================================
# CLI: setup-experiment
#===========================================================================================================

SETUP_HELP = (
    "Create a local FedGWAS simulation project from a preset.\n\n"
    "Available presets:\n"
    + "\n".join(f"- {preset.name}: {preset.description}" for preset in PRESETS.values())
)

SETUP_EXAMPLES = """Examples:\n
  fedgwas-sim setup-experiment syn-tiny\n
  fedgwas-sim setup-experiment syn-small --seed 42\n
  fedgwas-sim setup-experiment syn-tiny --out tiny-study\n
  fedgwas-sim setup-experiment hapmap --out hapmap-study\n
"""


@app.command(
    "setup-experiment",
    help=SETUP_HELP,
    epilog=SETUP_EXAMPLES,
)
def setup_experiment(
    preset: Annotated[PresetChoice, typer.Argument(help="Experiment preset name.")],
    out: Annotated[Optional[Path], typer.Option("--out", "-o", help="Project directory to create. Defaults to the current directory.")] = None,
    force: Annotated[bool, typer.Option("--force", help="Overwrite template files in an existing directory.")] = False,
    seed: Annotated[Optional[int], typer.Option("--seed", help="Synthetic data random seed.")] = None,
    download: Annotated[bool, typer.Option("--download/--no-download", help="Download and prepare real-data presets.")] = True,
) -> None:
    """Create a local FedGWAS simulation project from a preset.

    CLI usage:
        `fedgwas-sim setup-experiment syn-tiny --out my-study`
        `fedgwas-sim setup-experiment syn-small --out small-study --seed 42`
        `fedgwas-sim setup-experiment hapmap --out hapmap-study`

    Args:
        preset: The selected experiment preset as a `PresetChoice` enum value.
            Typer validates this argument before the command body runs, so users
            see the supported preset choices in `--help` and validation errors.
        out: Destination project directory. The command resolves this path before
            passing it to the project setup service. When omitted, the current
            directory is configured.
        force: Whether existing template files may be overwritten when `out`
            already exists. Non-empty directories still require explicit
            confirmation through this option.
        seed: Optional random seed persisted into `fedgwas.yaml` and passed to
            synthetic data generation for reproducible generated datasets.
        download: Whether real-data presets should execute their preparation
            script after project files are written.

    Returns:
        None. The command creates the simulation project on disk and prints the
        resolved destination path.

    Raises:
        typer.BadParameter: If the preset is not registered in `PRESETS`.
        OSError: If project files or directories cannot be created.

    The command keeps orchestration in the CLI layer only: it translates Typer
    arguments into the preset name and delegates all project layout, template,
    and optional synthetic-data generation work to `setup_mod.create_project`.
    """
    preset_name = preset.value
    if preset_name not in PRESETS:
        allowed = ", ".join(PRESETS)
        raise typer.BadParameter(f"Unknown preset '{preset_name}'. Choose one of: {allowed}.")
    project_dir = (out or Path.cwd()).resolve()
    setup_mod.create_project(
        project_dir,
        PRESETS[preset_name],
        force=force,
        seed=seed,
        download_real_data=download,
    )
    typer.echo(f"Created simulation project: {project_dir}")


#===========================================================================================================
# CLI: prepare-data
#===========================================================================================================

PREPARE_DATA_HELP = """Prepare or regenerate data for the current simulation project."""

PREPARE_DATA_EXAMPLES = """Examples:\n
  fedgwas-sim prepare-data\n
  fedgwas-sim prepare-data --download\n
"""


@app.command(
    "prepare-data",
    help=PREPARE_DATA_HELP,
    epilog=PREPARE_DATA_EXAMPLES,
)
def prepare_data(
    download: Annotated[bool, typer.Option("--download", help="Allow public dataset download for real-data presets.")] = False,
) -> None:
    """Prepare or regenerate data for the current simulation project.

    CLI usage:
        `fedgwas-sim prepare-data`
        `fedgwas-sim prepare-data --download`

    Args:
        download: Whether real-data preparation scripts are allowed to download
            public datasets. Synthetic presets ignore this flag and generate
            local PLINK data only when required files are missing.

    Returns:
        None. Synthetic presets either report that data is present or regenerate
        missing data. Real-data presets print the template script location unless
        `--download` is supplied.

    Raises:
        typer.BadParameter: If `fedgwas.yaml` references an unknown preset.
        typer.Exit: Exits with the return code from the real-data preparation
            script when `--download` is used.

    The command intentionally gates public dataset download behind `--download`
    so users can inspect the generated preparation script before network and
    storage-heavy work begins.
    """
    project_dir, project_config = validation.load_project_config()
    preset_name = str(project_config.get("preset", ""))
    preset = PRESETS.get(preset_name)
    if not preset:
        raise typer.BadParameter(f"Unknown preset in fedgwas.yaml: {preset_name}")
    if preset.synthetic and preset.scale:
        missing = validation.verify_data(project_dir, project_config, None)
        if missing:
            setup_mod.generate_synthetic_data(project_dir, preset.name, preset.scale, project_config.get("seed"))
            typer.echo("Synthetic data generated.")
        else:
            typer.echo("Synthetic data is already present.")
        return
    script = project_dir / "scripts" / "prepare_data.py"
    if not download:
        typer.echo(f"Preparation template is available at {script}. Re-run with --download after review.")
        return
    result = subprocess.run([sys.executable, str(script), "--download"], cwd=project_dir, check=False)
    raise typer.Exit(result.returncode)


#===========================================================================================================
# CLI: check
#===========================================================================================================

CHECK_HELP = (
    "Check whether the current directory is ready to run.\n\n"
    "Run all checks by default, or select one or more check categories."
)

CHECK_EXAMPLES = """Examples:\n
  fedgwas-sim check\n
  fedgwas-sim check --project\n
  fedgwas-sim check --software\n
  fedgwas-sim check --configs --data\n
  fedgwas-sim check --outputs\n
"""


@app.command(
    "check",
    help=CHECK_HELP,
    epilog=CHECK_EXAMPLES,
)
def check(
    project: Annotated[bool, typer.Option("--project", help="Check project marker and fedgwas.yaml only.")] = False,
    software: Annotated[bool, typer.Option("--software", help="Check Python, FedGWAS, Flower, and PLINK prerequisites.")] = False,
    configs: Annotated[bool, typer.Option("--configs", help="Check server and center config files.")] = False,
    data: Annotated[bool, typer.Option("--data", help="Check configured center PLINK data.")] = False,
    outputs: Annotated[bool, typer.Option("--outputs", help="Check writable results and logs directories.")] = False,
) -> None:
    """Check whether the current directory is ready to run a simulation.

    CLI usage:
        `fedgwas-sim check`
        `fedgwas-sim check --data`
        `fedgwas-sim check --configs --outputs`

    Args:
        project: Whether to check only the project marker/config.
        software: Whether to check only software prerequisites.
        configs: Whether to check only server and center configs.
        data: Whether to check only configured PLINK data.
        outputs: Whether to check only writable output directories.

    Returns:
        None. The command prints one status line per environment, project, data,
        config, and writable-directory check.

    Raises:
        typer.Exit: Exits with code 1 when `fedgwas.yaml` is missing, the
            project is not a simulation project, or any readiness check fails.

    The readiness checks cover Python version, package importability, Flower
    version, PLINK availability, center/server config paths, and writable
    results/logs directories.
    """
    selected = {
        "project": project,
        "software": software,
        "configs": configs,
        "data": data,
        "outputs": outputs,
    }
    failures = _run_checks(selected)
    if failures:
        raise typer.Exit(1)


def _run_checks(selected: dict[str, bool] | None = None) -> int:
    """Run simulation readiness checks and render the Rich check report.

    Args:
        selected: Mapping of check category name to enabled state. When omitted
            or when all values are false, every category is checked.

    Returns:
        Number of failed check rows. Callers decide whether to exit or continue.

    This helper is shared by `fedgwas-sim check` and the pre-flight validation
    inside `fedgwas-sim run`, so the visible report, failure counting, and next
    step guidance stay consistent.
    """
    selected = selected or {
        "project": False,
        "software": False,
        "configs": False,
        "data": False,
        "outputs": False,
    }
    if not any(selected.values()):
        selected = {key: True for key in selected}

    failures = 0
    sections: list[dict] = []
    project_dir: Path | None = None
    project_config: dict | None = None

    def load_config() -> tuple[Path, dict]:
        nonlocal project_dir, project_config
        if project_dir is None or project_config is None:
            project_dir, project_config = validation.load_project_config()
        return project_dir, project_config

    def add_section(title: str, checks: list[tuple[str, bool, str]]) -> None:
        nonlocal failures
        section_checks = []
        for label, ok, detail in checks:
            section_checks.append({"label": label, "ok": ok, "detail": detail})
            if not ok:
                failures += 1
        sections.append({"title": title, "checks": section_checks})

    if selected["project"]:
        try:
            loaded_dir, _loaded_config = load_config()
            add_section("Project", [("fedgwas.yaml", True, paths.display_path(loaded_dir / "fedgwas.yaml", loaded_dir))])
        except Exception as exc:
            add_section("Project", [("fedgwas.yaml", False, str(exc))])

    if selected["software"]:
        current_dir = project_dir or Path.cwd().resolve()
        current_config = project_config or {}
        checks: list[tuple[str, bool, str]] = []
        python_ok, python_detail = validation.python_version_ok()
        checks.append(("Python >=3.11", python_ok, python_detail))
        try:
            import pipeline  # noqa: F401

            checks.append(("FedGWAS importable", True, "pipeline"))
        except Exception as exc:
            checks.append(("FedGWAS importable", False, str(exc)))
        flwr_ok, flwr_detail = validation.flwr_version_ok()
        checks.append(("Flower >=1.19,<1.20", flwr_ok, flwr_detail))
        plink = validation.resolve_plink(current_dir, current_config)
        checks.append(("PLINK available", plink is not None, plink or "not found"))
        add_section("Software", checks)

    if selected["configs"]:
        try:
            loaded_dir, loaded_config = load_config()
            checks = []
            for center_id in paths.iter_center_ids(loaded_config):
                config_path = paths.center_config_path(loaded_dir, loaded_config, center_id)
                checks.append((f"center_{center_id} config", config_path.exists(), paths.display_path(config_path, loaded_dir)))
            server_config = paths.server_config_path(loaded_dir, loaded_config)
            checks.append(("server config", server_config.exists(), paths.display_path(server_config, loaded_dir)))
        except Exception as exc:
            checks = [("configuration", False, str(exc))]
        add_section("Configs", checks)

    if selected["data"]:
        try:
            loaded_dir, loaded_config = load_config()
            missing = validation.verify_data(loaded_dir, loaded_config, None)
            checks = [("configured PLINK triplets", not missing, "all present" if not missing else f"{len(missing)} missing")]
            checks.extend((paths.display_path(path, loaded_dir), False, "missing") for path in missing)
        except Exception as exc:
            checks = [("configured PLINK triplets", False, str(exc))]
        add_section("Data", checks)

    if selected["outputs"]:
        try:
            loaded_dir, loaded_config = load_config()
            checks = []
            for directory_name in ("results_dir", "logs_dir"):
                raw_value = loaded_config.get(directory_name, directory_name.replace("_dir", ""))
                path = loaded_dir / str(raw_value)
                checks.append((f"{directory_name} writable", validation.is_writable_dir(path), paths.display_path(path, loaded_dir)))
        except Exception as exc:
            checks = [("output directories", False, str(exc))]
        add_section("Outputs", checks)

    selected_categories = [name for name, enabled in selected.items() if enabled]
    summary.render_check_report(
        sections,
        project_dir or Path.cwd().resolve(),
        selected_categories,
        failures,
    )
    return failures


#===========================================================================================================
# CLI: reset / clear
#===========================================================================================================

RESET_HELP = (
    "Reset the current FedGWAS simulation project to the base initialized structure.\n\n"
    "This removes only CLI-managed paths and leaves unknown user files untouched."
)

RESET_EXAMPLES = """Examples:\n
  fedgwas-sim reset\n
  fedgwas-sim reset --yes\n
  fedgwas-sim reset --keep-data --yes\n
  fedgwas-sim clear --yes\n
"""


def _reset_project_command(
    *,
    yes: bool,
    keep_data: bool,
    keep_configs: bool,
    keep_results: bool,
) -> None:
    """Run the shared implementation behind `reset` and `clear`.

    Args:
        yes: Whether to skip the interactive confirmation prompt.
        keep_data: Whether to preserve `data/`.
        keep_configs: Whether to preserve `configs/`.
        keep_results: Whether to preserve `results/`.

    Returns:
        None. Removes allowlisted project paths and reinitializes the project.

    Raises:
        typer.Exit: Exits with code 1 when the user declines confirmation.
        typer.BadParameter: If the current directory is not a simulation project.
    """
    project_dir = Path.cwd().resolve()
    if not yes:
        typer.echo("This will reset the current FedGWAS simulation project.")
        typer.echo("CLI-managed configs, data, results, logs, scripts, and generated project files may be removed.")
        confirmed = typer.confirm("Continue?")
        if not confirmed:
            typer.echo("Reset cancelled.")
            raise typer.Exit(1)

    removed = setup_mod.reset_project(
        project_dir,
        keep_data=keep_data,
        keep_configs=keep_configs,
        keep_results=keep_results,
    )
    typer.echo(f"Reset simulation project: {project_dir}")
    if removed:
        typer.echo("Removed:")
        for path in removed:
            typer.echo(f"  {paths.display_path(path, project_dir)}")
    typer.echo("Unknown user files were left untouched.")


@app.command(
    "reset",
    help=RESET_HELP,
    epilog=RESET_EXAMPLES,
)
def reset(
    yes: Annotated[bool, typer.Option("--yes", "-y", help="Skip interactive confirmation.")] = False,
    keep_data: Annotated[bool, typer.Option("--keep-data", help="Preserve the data directory.")] = False,
    keep_configs: Annotated[bool, typer.Option("--keep-configs", help="Preserve the configs directory.")] = False,
    keep_results: Annotated[bool, typer.Option("--keep-results", help="Preserve the results directory.")] = False,
) -> None:
    """Reset the current FedGWAS simulation project.

    CLI usage:
        `fedgwas-sim reset`
        `fedgwas-sim reset --yes`
        `fedgwas-sim reset --keep-data --yes`

    Args:
        yes: Whether to skip the interactive confirmation prompt.
        keep_data: Whether to preserve the `data/` directory.
        keep_configs: Whether to preserve the `configs/` directory.
        keep_results: Whether to preserve the `results/` directory.

    Returns:
        None. The project is restored to the same base state created by
        `fedgwas-sim init`.

    Raises:
        typer.BadParameter: If the current directory is not a simulation project.
        typer.Exit: If the user declines confirmation.
    """
    _reset_project_command(
        yes=yes,
        keep_data=keep_data,
        keep_configs=keep_configs,
        keep_results=keep_results,
    )


@app.command(
    "clear",
    help="Alias for `fedgwas-sim reset`.",
    epilog=RESET_EXAMPLES,
)
def clear(
    yes: Annotated[bool, typer.Option("--yes", "-y", help="Skip interactive confirmation.")] = False,
    keep_data: Annotated[bool, typer.Option("--keep-data", help="Preserve the data directory.")] = False,
    keep_configs: Annotated[bool, typer.Option("--keep-configs", help="Preserve the configs directory.")] = False,
    keep_results: Annotated[bool, typer.Option("--keep-results", help="Preserve the results directory.")] = False,
) -> None:
    """Clear the current FedGWAS simulation project.

    CLI usage:
        `fedgwas-sim clear --yes`

    Args:
        yes: Whether to skip the interactive confirmation prompt.
        keep_data: Whether to preserve the `data/` directory.
        keep_configs: Whether to preserve the `configs/` directory.
        keep_results: Whether to preserve the `results/` directory.

    Returns:
        None. Delegates to the same implementation as `reset`.
    """
    _reset_project_command(
        yes=yes,
        keep_data=keep_data,
        keep_configs=keep_configs,
        keep_results=keep_results,
    )


#===========================================================================================================
# CLI group: data
#===========================================================================================================

DATA_CONFIGURE_HELP = """Point a center config at a PLINK binary prefix."""

DATA_CONFIGURE_EXAMPLES = """Examples:\n
  fedgwas-sim data configure --center 1 --bfile data/center_1/study_center_1\n
"""


@data_app.command(
    "configure",
    help=DATA_CONFIGURE_HELP,
    epilog=DATA_CONFIGURE_EXAMPLES,
)
def data_configure(
    center: Annotated[int, typer.Option("--center", "-c", help="One-based center id.")],
    bfile: Annotated[Path, typer.Option("--bfile", help="PLINK prefix without .bed/.bim/.fam.")],
) -> None:
    """Point a center config at a PLINK binary prefix.

    CLI usage:
        `fedgwas-sim data configure --center 1 --bfile data/center_1/study_center_1`

    Args:
        center: One-based center id whose config should be updated.
        bfile: PLINK binary prefix without the `.bed`, `.bim`, or `.fam`
            extension. Relative paths are normalized relative to the project.

    Returns:
        None. The command updates the selected center YAML file in place and
        prints the project-relative config path.

    Raises:
        typer.BadParameter: If the project config is invalid, `center` is out of
            range, or the selected center config cannot be found.
        OSError: If the config file cannot be read or written.

    The command only changes `input_data.path` and `input_data.type`, preserving
    the rest of the center configuration.
    """
    project_dir, project_config = validation.load_project_config()
    paths.iter_center_ids(project_config, center)
    config_path = paths.center_config_path(project_dir, project_config, center)
    if not config_path.exists():
        raise typer.BadParameter(f"Center config not found: {paths.display_path(config_path, project_dir)}")
    center_config = paths.load_yaml(config_path)
    center_config.setdefault("input_data", {})
    center_config["input_data"]["path"] = paths.normalize_config_path(bfile, project_dir)
    center_config["input_data"]["type"] = "bed"
    paths.write_yaml(config_path, center_config)
    typer.echo(f"Updated {paths.display_path(config_path, project_dir)}")


#===========================================================================================================
# CLI group: summarize
#===========================================================================================================

SUMMARIZE_DATA_HELP = """Summarize a prepared data directory."""

SUMMARIZE_DATA_EXAMPLES = """Examples:\n
  fedgwas-sim summarize data --path data\n
  fedgwas-sim summarize data --path data/center_1\n
"""

SUMMARIZE_EXPERIMENT_HELP = """Summarize a FedGWAS simulation project."""

SUMMARIZE_EXPERIMENT_EXAMPLES = """Examples:\n
  fedgwas-sim summarize experiment --path .\n
  fedgwas-sim summarize experiment --path tiny-study\n
"""


@summarize_app.command(
    "data",
    help=SUMMARIZE_DATA_HELP,
    epilog=SUMMARIZE_DATA_EXAMPLES,
)
def summarize_data(
    path: Annotated[Path, typer.Option("--path", help="Data directory to summarize.")] = Path("data"),
) -> None:
    """Summarize a prepared data directory.

    CLI usage:
        `fedgwas-sim summarize data --path data`

    Args:
        path: Data directory to inspect.

    Returns:
        None. Prints a report-style data inventory.
    """
    payload = summary.summarize_data(path)
    summary.render_data_summary(payload)


@summarize_app.command(
    "experiment",
    help=SUMMARIZE_EXPERIMENT_HELP,
    epilog=SUMMARIZE_EXPERIMENT_EXAMPLES,
)
def summarize_experiment(
    path: Annotated[Path, typer.Option("--path", help="Project directory to summarize.")] = Path("."),
) -> None:
    """Summarize a FedGWAS simulation project.

    CLI usage:
        `fedgwas-sim summarize experiment --path .`

    Args:
        path: Project directory to inspect.

    Returns:
        None. Prints a report-style project summary.
    """
    payload = summary.summarize_experiment(path)
    summary.render_experiment_summary(payload)


#===========================================================================================================
# CLI: run
#===========================================================================================================

RUN_HELP = """Run Flower local simulation after all readiness checks pass."""

RUN_EXAMPLES = """Examples:\n
  fedgwas-sim run --rounds 100\n
  fedgwas-sim run --rounds 100 --no-stream\n
  fedgwas-sim check --data\n
"""


@app.command(
    "run",
    help=RUN_HELP,
    epilog=RUN_EXAMPLES,
)
def run(
    rounds: Annotated[int, typer.Option("--rounds", "-r", help="Number of Flower server rounds.")] = 100,
    stream: Annotated[bool, typer.Option("--stream/--no-stream", help="Stream Flower output.")] = True,
) -> None:
    """Run Flower local simulation after pre-flight checks pass.

    CLI usage:
        `fedgwas-sim run --rounds 100`
        `fedgwas-sim run --rounds 100 --no-stream`

    Args:
        rounds: Number of Flower server rounds to pass as
            `num-server-rounds` in `--run-config`.
        stream: Whether to add Flower's `--stream` flag and stream runtime
            output to the terminal.

    Returns:
        None. The process exits with code 1 if readiness checks fail, or with
        the return code from `flwr run` after the subprocess completes.

    Raises:
        typer.Exit: Exits with code 1 when checks fail. Otherwise exits with
            Flower's return code after the subprocess completes.
        typer.BadParameter: If the current directory is not a valid simulation
            project.

    The command first runs the same project, software, config, data, and output
    readiness checks as `fedgwas-sim check`. The Rich report is printed before
    execution. Flower is launched only when every check passes.
    """
    failures = _run_checks()
    if failures:
        raise typer.Exit(1)

    project_dir, project_config = validation.load_project_config()
    config_dir = str(project_config.get("config_dir", "configs"))
    cmd = ["flwr", "run", ".", "local-simulation"]
    if stream:
        cmd.append("--stream")
    cmd.extend(
        [
            "--run-config",
            f'simulation=true num-server-rounds={rounds} config_path="{config_dir}"',
        ]
    )
    result = subprocess.run(cmd, cwd=project_dir, check=False)
    raise typer.Exit(result.returncode)


#===========================================================================================================
# CLI group: baseline
#===========================================================================================================

BASELINE_GENERATE_HELP = """Generate a centralized PLINK baseline."""

BASELINE_GENERATE_EXAMPLES = """Examples:\n
  fedgwas-sim baseline generate\n
  fedgwas-sim baseline generate --output results/baseline\n
"""


@baseline_app.command(
    "generate",
    help=BASELINE_GENERATE_HELP,
    epilog=BASELINE_GENERATE_EXAMPLES,
)
def baseline_generate(
    output: Annotated[Optional[Path], typer.Option("--output", "-o", help="Baseline output directory.")] = None,
) -> None:
    """Generate a centralized PLINK baseline for comparison.

    CLI usage:
        `fedgwas-sim baseline generate`
        `fedgwas-sim baseline generate --output results/baseline`

    Args:
        output: Optional destination directory for baseline artifacts. When
            omitted, the command writes to `results/baseline` inside the current
            simulation project.

    Returns:
        None. The command prints the generated baseline artifact paths reported
        by the output service.

    Raises:
        typer.BadParameter: If the project config is invalid, PLINK is
            unavailable, or no usable center data is configured.
        subprocess.CalledProcessError: If the underlying PLINK command fails.

    Baseline generation lives in `pipeline.tools.generate_baseline`; the CLI
    only adapts project paths and prints the returned artifact mapping.
    """
    project_dir, project_config = validation.load_project_config()
    if output is None:
        output_dir = paths.results_dir(project_dir, project_config) / "baseline"
    else:
        output_dir = output if output.is_absolute() else project_dir / output
    results = baseline_tool.generate_baseline(project_dir / "config.yaml", output_dir)
    typer.echo("Baseline generation completed.")
    for key, value in results.items():
        typer.echo(f"{key}: {value}")


#===========================================================================================================
# CLI: evaluate
#===========================================================================================================

EVALUATE_HELP = """Run QC, LR, and optional KING evaluation for federated outputs."""

EVALUATE_EXAMPLES = """Examples:\n
  fedgwas-sim evaluate\n
  fedgwas-sim evaluate --baseline results/baseline --king\n
  fedgwas-sim evaluate results --baseline results/baseline --qc-only\n
  fedgwas-sim evaluate --baseline results/baseline --king-only --king-center-id 1\n
"""


@app.command(
    "evaluate",
    help=EVALUATE_HELP,
    epilog=EVALUATE_EXAMPLES,
)
def evaluate(
    results_dir: Annotated[Optional[Path], typer.Argument(help="Federated results directory. Defaults to the project results directory.")] = None,
    baseline: Annotated[Optional[Path], typer.Option("--baseline", help="Baseline results directory.")] = None,
    report: Annotated[Optional[Path], typer.Option("--report", help="Combined evaluation report path.")] = None,
    no_plots: Annotated[bool, typer.Option("--no-plots", help="Skip LR plot generation.")] = False,
    king: Annotated[bool, typer.Option("--king", help="Also run KING evaluation.")] = False,
    king_center_id: Annotated[int, typer.Option("--king-center-id", help="Center id for KING accumulator evaluation.")] = 1,
    king_data_dir: Annotated[Optional[Path], typer.Option("--king-data-dir", help="Optional KING data root containing center_* dirs.")] = None,
    qc_only: Annotated[bool, typer.Option("--qc-only", help="Run QC evaluation only.")] = False,
    lr_only: Annotated[bool, typer.Option("--lr-only", help="Run LR evaluation only.")] = False,
    king_only: Annotated[bool, typer.Option("--king-only", help="Run KING evaluation only.")] = False,
) -> None:
    """Run core evaluation for federated outputs and baseline artifacts.

    CLI usage:
        `fedgwas-sim evaluate`
        `fedgwas-sim evaluate --baseline results/baseline --king`
        `fedgwas-sim evaluate results --baseline results/baseline --qc-only`

    Args:
        results_dir: Optional federated results directory. When omitted, the
            command uses the project `results_dir` from `fedgwas.yaml`.
        baseline: Optional baseline results directory. When omitted, the command
            uses `results/baseline` inside the current simulation project.
        report: Optional combined summary report path. When omitted, the command
            writes `evaluation_report.md` inside the selected results directory.
        no_plots: Whether to skip LR plot generation.
        king: Whether to run KING evaluation in addition to QC and LR.
        king_center_id: Center id for KING accumulator evaluation.
        king_data_dir: Optional data root for KING id-map lookup.
        qc_only: Whether to run QC evaluation only.
        lr_only: Whether to run LR evaluation only.
        king_only: Whether to run KING evaluation only.

    Returns:
        None. The command delegates report generation to
        `pipeline.evaluation.evaluate_all.run_evaluation` and prints generated
        report paths.

    Raises:
        typer.BadParameter: If the current directory is not a simulation project
            or the selected results/baseline paths are invalid.
        OSError: If evaluator reports cannot be written.

    Evaluation logic lives in `pipeline.evaluation.evaluate_all`; the CLI only
    adapts simulation project defaults and user-provided relative paths.
    """
    project_dir, project_config = validation.load_project_config()

    def resolve_project_path(value: Optional[Path], default: Path) -> Path:
        if value is None:
            return default
        return value if value.is_absolute() else project_dir / value

    resolved_results_dir = resolve_project_path(results_dir, paths.results_dir(project_dir, project_config))
    resolved_baseline = resolve_project_path(baseline, resolved_results_dir / "baseline")
    resolved_report = resolve_project_path(report, resolved_results_dir / "evaluation_report.md")
    resolved_king_data_dir = None
    if king_data_dir is not None:
        resolved_king_data_dir = king_data_dir if king_data_dir.is_absolute() else project_dir / king_data_dir

    try:
        result = evaluation_tool.run_evaluation(
            results_dir=resolved_results_dir,
            baseline_dir=resolved_baseline,
            report_path=resolved_report,
            no_plots=no_plots,
            king=king,
            king_center_id=king_center_id,
            king_data_dir=resolved_king_data_dir,
            qc_only=qc_only,
            lr_only=lr_only,
            king_only=king_only,
        )
    except (FileNotFoundError, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from exc

    typer.echo("Evaluation completed.")
    for key in ("report_path", "qc_report_path", "lr_report_path", "king_report_path"):
        value = result.get(key)
        if value:
            typer.echo(f"{key}: {paths.display_path(Path(value), project_dir)}")


#===========================================================================================================
# CLI group: results
#===========================================================================================================

RESULTS_COLLECT_HELP = """Collect metrics and write run summary files."""

RESULTS_COLLECT_EXAMPLES = """Examples:\n
  fedgwas-sim results collect\n
  fedgwas-sim results collect --time-file results/server_app_time.txt --label local_tiny\n
"""


@results_app.command(
    "collect",
    help=RESULTS_COLLECT_HELP,
    epilog=RESULTS_COLLECT_EXAMPLES,
)
def results_collect(
    time_file: Annotated[Optional[list[Path]], typer.Option("--time-file", help="Time output file. Repeatable.")] = None,
    label: Annotated[str, typer.Option("--label", help="Label for run_summary.md.")] = "experiment",
    output_dir: Annotated[Optional[Path], typer.Option("--output-dir", help="Directory for run_summary files.")] = None,
) -> None:
    """Collect metrics and write run_summary.json/md.

    CLI usage:
        `fedgwas-sim results collect`
        `fedgwas-sim results collect --time-file results/server_app_time.txt --label local_tiny`

    Args:
        time_file: Optional repeatable list of GNU time output files to include
            in the collected run summary.
        label: Human-readable run label written into `run_summary.md`.
        output_dir: Optional output directory for `run_summary.json` and
            `run_summary.md`. Defaults to the configured project results
            directory.

    Returns:
        None. The command delegates metric extraction and summary writing to
        `outputs.collect_metrics`.

    Raises:
        typer.BadParameter: If the current directory is not a valid simulation
            project.
        OSError: If summary files cannot be written.

    The command keeps the CLI layer small by resolving project paths here and
    leaving metric parsing and report serialization to the output service.
    """
    project_dir, project_config = validation.load_project_config()
    results_dir = paths.results_dir(project_dir, project_config)
    out_dir = output_dir or results_dir
    outputs.collect_metrics(results_dir, time_file or [], label, out_dir)


#===========================================================================================================
# CLI entrypoint
#===========================================================================================================

def main() -> None:
    """Run the simulation Typer app.

    CLI usage:
        `fedgwas-sim ...`

    Args:
        None.

    Returns:
        None. Typer owns command dispatch and process exit semantics.
    """
    app()
