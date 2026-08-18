from __future__ import annotations

import sys
from pathlib import Path
import re

import pytest
import yaml
from typer.testing import CliRunner

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline.client_app import _resolve_simulation_client_config_path
from pipeline.cli import sim
from pipeline.cli.simulation import setup as setup_mod
from pipeline.evaluation import evaluate_all as evaluation_tool
from pipeline.tools import generate_baseline as baseline_tool
from pipeline.server_app import _resolve_server_config_file


runner = CliRunner()
REPO_ROOT = Path(__file__).resolve().parents[1]


def _help_output(*args: str) -> str:
    result = runner.invoke(sim.app, [*args, "--help"])
    assert result.exit_code == 0, result.output
    return result.output


def _compact(text: str) -> str:
    return re.sub(r"\s+", " ", text)


def _write_triplet(prefix: Path) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".bed", ".bim", ".fam"):
        (prefix.with_suffix(suffix)).write_text("placeholder\n")


def _read_yaml(path: Path) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def _reference_experiment(example_name: str) -> Path:
    return {
        "syn-tiny": REPO_ROOT / "experiments" / "correctness" / "tiny_even",
        "syn-small": REPO_ROOT / "experiments" / "performance" / "small_even",
        "syn-medium": REPO_ROOT / "experiments" / "performance" / "medium_even",
        "1000genomes": REPO_ROOT / "experiments" / "real_world" / "1000genomes",
    }[example_name]


def _config_value(config: dict, dotted_key: str):
    value = config
    for part in dotted_key.split("."):
        if not isinstance(value, dict) or part not in value:
            return None
        value = value[part]
    return value


def _normalize_numbers(value):
    if isinstance(value, dict):
        return {key: _normalize_numbers(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_normalize_numbers(item) for item in value]
    if isinstance(value, str):
        try:
            return float(value)
        except ValueError:
            return value
    return value


def _semantic_config(config: dict) -> dict:
    payload = {
        key: _normalize_numbers(_config_value(config, key))
        for key in (
            "experiment_name",
            "experiment_category",
            "description",
            "scenario",
            "data.scale",
            "data.partition_strategy",
            "server.chunk_size",
            "server.num_server_rounds",
            "clients.participation",
            "clients.thresholds",
        )
        if _config_value(config, key) is not None
    }
    payload["num_clients"] = config.get("num_clients") or len(config.get("clients", {}).get("config_files", {}))
    return payload


def _semantic_center_config(config: dict) -> dict:
    parameters = dict(config.get("parameters", {}))
    parameters.pop("run_id", None)
    return {
        "parameters": _normalize_numbers(parameters),
        "thresholds": _normalize_numbers(config.get("thresholds", {})),
        "participation": config.get("participation", {}),
        "input_type": config.get("input_data", {}).get("type"),
    }


def test_root_help_shows_quickstart_examples():
    output = _help_output()
    compact = _compact(output)

    assert "Examples:" in output
    assert "fedgwas-sim init" in output
    assert "fedgwas-sim setup-experiment syn-tiny" in compact
    assert "fedgwas-sim run --rounds 50" in compact
    assert "fedgwas-sim results collect --label tiny_run" in compact


def test_setup_experiment_help_shows_presets_and_hides_developer_doc_sections():
    output = _help_output("setup-experiment")

    for preset in ("syn-tiny", "syn-small", "syn-medium", "1000genomes", "hapmap"):
        assert preset in output
    assert "1000genomes-chr22" not in output
    assert "Examples:" in output
    assert "fedgwas-sim setup-experiment syn-tiny" in output
    assert "fedgwas-sim setup-experiment syn-tiny --out tiny-study" in output
    assert "fedgwas-sim setup-experiment hapmap --out hapmap-study" in output
    for developer_section in ("Args:", "Returns:", "Raises:"):
        assert developer_section not in output


def test_command_help_shows_examples_without_developer_doc_sections():
    commands = [
        ("data", "configure"),
        ("summarize", "data"),
        ("summarize", "experiment"),
        ("reset",),
        ("run",),
        ("baseline", "generate"),
        ("evaluate",),
        ("results", "collect"),
    ]

    for command in commands:
        output = _help_output(*command)
        assert "Examples:" in output
        for developer_section in ("Args:", "Returns:", "Raises:"):
            assert developer_section not in output


def test_init_help_shows_examples_and_choices():
    output = _help_output("init")

    for example in ("syn-tiny", "syn-small", "syn-medium", "1000genomes"):
        assert example in output
    assert "fedgwas-sim init" in output
    assert "fedgwas-sim init --from-example syn-tiny" in output
    for old_name in ("tiny-even", "small-even", "medium-even"):
        assert old_name not in output
    assert "--prepare-data" in output
    assert "--no-prepare-data" in output


@pytest.mark.parametrize("old_name", ("tiny-even", "small-even", "medium-even", "1000genomes-chr22"))
def test_old_example_and_preset_names_are_removed(tmp_path, monkeypatch, old_name):
    monkeypatch.chdir(tmp_path)

    init_result = runner.invoke(sim.app, ["init", "--from-example", old_name])
    setup_result = runner.invoke(sim.app, ["setup-experiment", old_name])

    assert init_result.exit_code != 0
    assert setup_result.exit_code != 0


def test_init_creates_base_project_structure(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    result = runner.invoke(sim.app, ["init"])

    assert result.exit_code == 0, result.output
    for path in (
        "fedgwas.yaml",
        "pyproject.toml",
        "config.yaml",
        "configs",
        "data",
        "results",
        "logs",
    ):
        assert (tmp_path / path).exists()
    project_config = _read_yaml(tmp_path / "fedgwas.yaml")
    assert project_config["mode"] == "simulation"
    assert project_config["project_state"] == "initialized"


def test_init_from_synthetic_example_prepares_data_by_default(tmp_path, monkeypatch):
    generated = []

    def fake_generate(project_dir: Path, preset_name: str, scale_name: str, seed: int | None) -> None:
        generated.append((project_dir, preset_name, scale_name, seed))
        for center_id in (1, 2):
            _write_triplet(project_dir / "data" / f"center_{center_id}" / f"{scale_name}_center_{center_id}")

    monkeypatch.setattr(setup_mod, "generate_synthetic_data", fake_generate)
    monkeypatch.chdir(tmp_path)

    result = runner.invoke(sim.app, ["init", "--from-example", "syn-tiny", "--seed", "123"])

    assert result.exit_code == 0, result.output
    assert generated == [(tmp_path, "syn-tiny", "tiny", 123)]
    assert (tmp_path / "data" / "center_1" / "tiny_center_1.bed").is_file()


def test_init_from_synthetic_example_can_skip_data_preparation(tmp_path, monkeypatch):
    generated = []

    def fake_generate(project_dir: Path, preset_name: str, scale_name: str, seed: int | None) -> None:
        generated.append((project_dir, preset_name, scale_name, seed))

    monkeypatch.setattr(setup_mod, "generate_synthetic_data", fake_generate)
    monkeypatch.chdir(tmp_path)

    result = runner.invoke(sim.app, ["init", "--from-example", "syn-tiny", "--no-prepare-data"])

    assert result.exit_code == 0, result.output
    assert generated == []
    assert not (tmp_path / "data" / "center_1" / "tiny_center_1.bed").exists()


def test_init_from_real_example_prepares_data_by_default(tmp_path, monkeypatch):
    calls = []

    def fake_prepare(project_dir: Path, preset, download: bool) -> None:
        calls.append((project_dir, preset.name, download))

    monkeypatch.setattr(setup_mod, "prepare_real_data", fake_prepare, raising=False)
    monkeypatch.chdir(tmp_path)

    result = runner.invoke(sim.app, ["init", "--from-example", "1000genomes"])

    assert result.exit_code == 0, result.output
    assert calls == [(tmp_path, "1000genomes", True)]


def test_init_from_real_example_can_skip_data_preparation(tmp_path, monkeypatch):
    calls = []

    def fake_prepare(project_dir: Path, preset, download: bool) -> None:
        calls.append((project_dir, preset.name, download))

    monkeypatch.setattr(setup_mod, "prepare_real_data", fake_prepare, raising=False)
    monkeypatch.chdir(tmp_path)

    result = runner.invoke(sim.app, ["init", "--from-example", "1000genomes", "--no-prepare-data"])

    assert result.exit_code == 0, result.output
    assert calls == []
    assert (tmp_path / "scripts" / "prepare_data.py").is_file()


def test_reset_fails_outside_simulation_project(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "notes.md").write_text("do not touch\n", encoding="utf-8")

    result = runner.invoke(sim.app, ["reset", "--yes"])

    assert result.exit_code != 0
    assert (tmp_path / "notes.md").read_text(encoding="utf-8") == "do not touch\n"


def test_reset_requires_confirmation_by_default(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    init = runner.invoke(sim.app, ["init"])
    assert init.exit_code == 0, init.output
    (tmp_path / "results" / "marker.txt").write_text("run output\n", encoding="utf-8")

    result = runner.invoke(sim.app, ["reset"], input="n\n")

    assert result.exit_code == 1
    assert (tmp_path / "results" / "marker.txt").is_file()


def test_reset_yes_restores_initialized_project_and_preserves_unknown_files(tmp_path, monkeypatch):
    def fake_generate(project_dir: Path, preset_name: str, scale_name: str, seed: int | None) -> None:
        for center_id in (1, 2):
            _write_triplet(project_dir / "data" / f"center_{center_id}" / f"{scale_name}_center_{center_id}")

    monkeypatch.setattr(setup_mod, "generate_synthetic_data", fake_generate)
    monkeypatch.chdir(tmp_path)
    setup = runner.invoke(sim.app, ["setup-experiment", "syn-tiny"])
    assert setup.exit_code == 0, setup.output
    (tmp_path / "notes.md").write_text("keep me\n", encoding="utf-8")

    result = runner.invoke(sim.app, ["reset", "--yes"])

    assert result.exit_code == 0, result.output
    assert (tmp_path / "notes.md").read_text(encoding="utf-8") == "keep me\n"
    assert not (tmp_path / "data" / "center_1" / "tiny_center_1.bed").exists()
    assert not (tmp_path / "configs" / "center_1.yaml").exists()
    assert (tmp_path / "configs").is_dir()
    assert (tmp_path / "data").is_dir()
    assert (tmp_path / "results").is_dir()
    project_config = _read_yaml(tmp_path / "fedgwas.yaml")
    assert project_config["mode"] == "simulation"
    assert project_config["project_state"] == "initialized"
    assert "preset" not in project_config


def test_reset_keep_data_preserves_data_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    init = runner.invoke(sim.app, ["init"])
    assert init.exit_code == 0, init.output
    marker = tmp_path / "data" / "center_1" / "manual.bed"
    marker.parent.mkdir(parents=True)
    marker.write_text("large data\n", encoding="utf-8")

    result = runner.invoke(sim.app, ["reset", "--keep-data", "--yes"])

    assert result.exit_code == 0, result.output
    assert marker.read_text(encoding="utf-8") == "large data\n"
    assert _read_yaml(tmp_path / "fedgwas.yaml")["project_state"] == "initialized"


def test_clear_alias_matches_reset(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    init = runner.invoke(sim.app, ["init"])
    assert init.exit_code == 0, init.output
    (tmp_path / "logs" / "old.log").write_text("old\n", encoding="utf-8")

    result = runner.invoke(sim.app, ["clear", "--yes"])

    assert result.exit_code == 0, result.output
    assert not (tmp_path / "logs" / "old.log").exists()
    assert _read_yaml(tmp_path / "fedgwas.yaml")["project_state"] == "initialized"


@pytest.mark.parametrize("example_name", ("syn-tiny", "syn-small", "syn-medium", "1000genomes"))
def test_init_from_example_matches_reference_config(tmp_path, monkeypatch, example_name):
    monkeypatch.chdir(tmp_path)

    result = runner.invoke(sim.app, ["init", "--from-example", example_name, "--no-prepare-data"])

    assert result.exit_code == 0, result.output
    generated_config = _read_yaml(tmp_path / "config.yaml")
    reference_dir = _reference_experiment(example_name)
    reference_config = _read_yaml(reference_dir / "config.yaml")
    assert _semantic_config(generated_config) == _semantic_config(reference_config)

    generated_center = _read_yaml(tmp_path / "configs" / "center_1.yaml")
    reference_center = _read_yaml(reference_dir / "configs" / "center_1" / "config.yaml")
    assert _semantic_center_config(generated_center) == _semantic_center_config(reference_center)
    assert (tmp_path / "configs" / "server.yaml").is_file()


@pytest.mark.parametrize(
    ("preset_name", "example_name"),
    (("syn-tiny", "syn-tiny"), ("syn-small", "syn-small"), ("syn-medium", "syn-medium")),
)
def test_setup_experiment_presets_match_reference_config(tmp_path, monkeypatch, preset_name, example_name):
    def fake_generate(project_dir: Path, preset_name: str, scale_name: str, seed: int | None) -> None:
        for center_id in (1, 2):
            _write_triplet(project_dir / "data" / f"center_{center_id}" / f"{scale_name}_center_{center_id}")

    monkeypatch.setattr(setup_mod, "generate_synthetic_data", fake_generate)
    monkeypatch.chdir(tmp_path)

    init = runner.invoke(sim.app, ["init"])
    assert init.exit_code == 0, init.output
    result = runner.invoke(sim.app, ["setup-experiment", preset_name])

    assert result.exit_code == 0, result.output
    generated_config = _read_yaml(tmp_path / "config.yaml")
    reference_dir = _reference_experiment(example_name)
    reference_config = _read_yaml(reference_dir / "config.yaml")
    assert _semantic_config(generated_config) == _semantic_config(reference_config)
    generated_center = _read_yaml(tmp_path / "configs" / "center_1.yaml")
    reference_center = _read_yaml(reference_dir / "configs" / "center_1" / "config.yaml")
    assert _semantic_center_config(generated_center) == _semantic_center_config(reference_center)


def test_setup_experiment_configures_current_project(tmp_path, monkeypatch):
    generated = []

    def fake_generate(project_dir: Path, preset_name: str, scale_name: str, seed: int | None) -> None:
        generated.append((project_dir, preset_name, scale_name, seed))
        for center_id in (1, 2):
            _write_triplet(project_dir / "data" / f"center_{center_id}" / f"{scale_name}_center_{center_id}")

    monkeypatch.setattr(setup_mod, "generate_synthetic_data", fake_generate)
    monkeypatch.chdir(tmp_path)

    init = runner.invoke(sim.app, ["init"])
    assert init.exit_code == 0, init.output
    result = runner.invoke(sim.app, ["setup-experiment", "syn-tiny"])

    assert result.exit_code == 0, result.output
    assert generated == [(tmp_path, "syn-tiny", "tiny", None)]
    assert (tmp_path / "configs" / "center_1.yaml").is_file()
    assert (tmp_path / "data" / "center_1" / "tiny_center_1.bed").is_file()


def test_setup_experiment_real_preset_download_flag_can_be_disabled(tmp_path, monkeypatch):
    calls = []

    def fake_prepare(project_dir: Path, preset, download: bool) -> None:
        calls.append((project_dir, preset.name, download))

    monkeypatch.setattr(setup_mod, "prepare_real_data", fake_prepare, raising=False)

    no_download_dir = tmp_path / "hapmap-template"
    no_download = runner.invoke(
        sim.app,
        ["setup-experiment", "hapmap", "--out", str(no_download_dir), "--no-download"],
    )
    assert no_download.exit_code == 0, no_download.output
    assert calls == [(no_download_dir, "hapmap", False)]

    default_download_dir = tmp_path / "hapmap-download"
    default_download = runner.invoke(sim.app, ["setup-experiment", "hapmap", "--out", str(default_download_dir)])
    assert default_download.exit_code == 0, default_download.output
    assert calls[-1] == (default_download_dir, "hapmap", True)


def test_setup_experiment_creates_runnable_project(tmp_path, monkeypatch):
    generated = []

    def fake_generate(project_dir: Path, preset_name: str, scale_name: str, seed: int | None) -> None:
        generated.append((project_dir, preset_name, scale_name, seed))
        for center_id in (1, 2):
            _write_triplet(project_dir / "data" / f"center_{center_id}" / f"{scale_name}_center_{center_id}")

    monkeypatch.setattr(setup_mod, "generate_synthetic_data", fake_generate)

    project_dir = tmp_path / "tiny-study"
    result = runner.invoke(sim.app, ["setup-experiment", "syn-tiny", "--out", str(project_dir)])

    assert result.exit_code == 0, result.output
    assert generated == [(project_dir, "syn-tiny", "tiny", None)]
    assert (project_dir / "fedgwas.yaml").is_file()
    assert (project_dir / "pyproject.toml").is_file()
    assert (project_dir / "configs" / "server.yaml").is_file()
    assert (project_dir / "configs" / "center_1.yaml").is_file()
    assert (project_dir / "configs" / "center_2.yaml").is_file()
    assert (project_dir / "results" / "server" / "logs").is_dir()
    assert (project_dir / "logs").is_dir()

    project_config = yaml.safe_load((project_dir / "fedgwas.yaml").read_text())
    assert project_config["mode"] == "simulation"
    assert project_config["preset"] == "syn-tiny"
    assert project_config["num_clients"] == 2

    center_config = yaml.safe_load((project_dir / "configs" / "center_1.yaml").read_text())
    assert center_config["input_data"]["path"] == "data/center_1/tiny_center_1"
    assert "experiments/" not in (project_dir / "configs" / "center_1.yaml").read_text()


def test_check_help_shows_category_filters():
    output = _help_output("check")

    for flag in ("--project", "--software", "--configs", "--data", "--outputs"):
        assert flag in output


def test_check_renders_rich_report_for_all_categories(tmp_path, monkeypatch):
    def fake_generate(project_dir: Path, preset_name: str, scale_name: str, seed: int | None) -> None:
        for center_id in (1, 2):
            _write_triplet(project_dir / "data" / f"center_{center_id}" / f"{scale_name}_center_{center_id}")

    monkeypatch.setattr(setup_mod, "generate_synthetic_data", fake_generate)
    monkeypatch.setattr(sim.validation, "flwr_version_ok", lambda: (True, "1.19.0"))
    monkeypatch.setattr(sim.validation, "resolve_plink", lambda project_dir, project_config: "plink")
    setup = runner.invoke(sim.app, ["setup-experiment", "syn-tiny", "--out", str(tmp_path)])
    assert setup.exit_code == 0, setup.output
    monkeypatch.chdir(tmp_path)

    result = runner.invoke(sim.app, ["check"])

    assert result.exit_code == 0, result.output
    for label in ("Check Summary", "Project", "Software", "Configs", "Data", "Outputs"):
        assert label in result.output
    assert "Aspect" in result.output
    assert result.output.count("Status") == 1
    assert "PASS" in result.output
    assert "failed" in result.output


def test_check_data_only_ignores_software_failures(tmp_path, monkeypatch):
    def fake_generate(project_dir: Path, preset_name: str, scale_name: str, seed: int | None) -> None:
        for center_id in (1, 2):
            _write_triplet(project_dir / "data" / f"center_{center_id}" / f"{scale_name}_center_{center_id}")

    monkeypatch.setattr(setup_mod, "generate_synthetic_data", fake_generate)
    monkeypatch.setattr(sim.validation, "flwr_version_ok", lambda: (False, "missing"))
    monkeypatch.setattr(sim.validation, "resolve_plink", lambda project_dir, project_config: None)
    setup = runner.invoke(sim.app, ["setup-experiment", "syn-tiny", "--out", str(tmp_path)])
    assert setup.exit_code == 0, setup.output
    monkeypatch.chdir(tmp_path)

    result = runner.invoke(sim.app, ["check", "--data"])

    assert result.exit_code == 0, result.output
    assert "Check Summary" in result.output
    assert "Data" in result.output
    assert "Flower" not in result.output
    assert "PLINK available" not in result.output


def test_check_configs_only_ignores_missing_data(tmp_path, monkeypatch):
    monkeypatch.setattr(sim.validation, "flwr_version_ok", lambda: (False, "missing"))
    monkeypatch.setattr(sim.validation, "resolve_plink", lambda project_dir, project_config: None)
    monkeypatch.chdir(tmp_path)
    init = runner.invoke(sim.app, ["init", "--from-example", "syn-tiny", "--no-prepare-data"])
    assert init.exit_code == 0, init.output

    result = runner.invoke(sim.app, ["check", "--configs"])

    assert result.exit_code == 0, result.output
    assert "Configs" in result.output
    assert "Flower" not in result.output
    assert "tiny_center_1.bed" not in result.output


def test_check_failure_renders_next_steps(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    init = runner.invoke(sim.app, ["init", "--from-example", "syn-tiny", "--no-prepare-data"])
    assert init.exit_code == 0, init.output

    result = runner.invoke(sim.app, ["check", "--data"])

    assert result.exit_code == 1
    assert "FAIL" in result.output
    assert "Next Steps" in result.output
    assert "fedgwas-sim prepare-data" in result.output
    assert "fedgwas-sim check --data" in result.output
    assert "fedgwas-sim data verify" not in result.output


def test_summarize_data_renders_human_readable_report(tmp_path):
    _write_triplet(tmp_path / "data" / "center_1" / "tiny_center_1")
    _write_triplet(tmp_path / "data" / "center_2" / "tiny_center_2")

    result = runner.invoke(sim.app, ["summarize", "data", "--path", str(tmp_path / "data")])

    assert result.exit_code == 0, result.output
    assert "Data Summary" in result.output
    assert "PLINK triplets" in result.output
    assert "center_1" in result.output
    assert "center_2" in result.output


def test_summarize_experiment_renders_human_readable_report(tmp_path, monkeypatch):
    def fake_generate(project_dir: Path, preset_name: str, scale_name: str, seed: int | None) -> None:
        for center_id in (1, 2):
            _write_triplet(project_dir / "data" / f"center_{center_id}" / f"{scale_name}_center_{center_id}")

    monkeypatch.setattr(setup_mod, "generate_synthetic_data", fake_generate)
    setup = runner.invoke(sim.app, ["setup-experiment", "syn-tiny", "--out", str(tmp_path)])
    assert setup.exit_code == 0, setup.output

    result = runner.invoke(sim.app, ["summarize", "experiment", "--path", str(tmp_path)])

    assert result.exit_code == 0, result.output
    assert "Experiment Summary" in result.output
    assert "syn-tiny" in result.output
    assert "tiny_even" in result.output
    assert "Clients" in result.output


def test_check_data_reports_missing_plink_files(tmp_path, monkeypatch):
    def fake_generate(project_dir: Path, preset_name: str, scale_name: str, seed: int | None) -> None:
        _write_triplet(project_dir / "data" / "center_1" / "tiny_center_1")
        _write_triplet(project_dir / "data" / "center_2" / "tiny_center_2")

    monkeypatch.setattr(setup_mod, "generate_synthetic_data", fake_generate)
    project_dir = tmp_path / "tiny-study"
    setup = runner.invoke(sim.app, ["setup-experiment", "syn-tiny", "--out", str(project_dir)])
    assert setup.exit_code == 0, setup.output

    (project_dir / "data" / "center_1" / "tiny_center_1.bim").unlink()
    monkeypatch.chdir(project_dir)

    result = runner.invoke(sim.app, ["check", "--data"])

    assert result.exit_code == 1
    assert "data/center_1/tiny_center_1.bim" in result.output.replace("\\", "/")


def test_data_verify_command_is_removed():
    result = runner.invoke(sim.app, ["data", "verify", "--help"])

    assert result.exit_code != 0
    assert "No such command" in result.output or "No such option" in result.output


def test_data_configure_updates_center_config_and_verify_passes(tmp_path, monkeypatch):
    def fake_generate(project_dir: Path, preset_name: str, scale_name: str, seed: int | None) -> None:
        _write_triplet(project_dir / "data" / "center_1" / "tiny_center_1")
        _write_triplet(project_dir / "data" / "center_2" / "tiny_center_2")

    monkeypatch.setattr(setup_mod, "generate_synthetic_data", fake_generate)
    project_dir = tmp_path / "tiny-study"
    setup = runner.invoke(sim.app, ["setup-experiment", "syn-tiny", "--out", str(project_dir)])
    assert setup.exit_code == 0, setup.output

    custom_prefix = project_dir / "data" / "center_1" / "custom_center_1"
    _write_triplet(custom_prefix)
    monkeypatch.chdir(project_dir)

    configure = runner.invoke(
        sim.app,
        ["data", "configure", "--center", "1", "--bfile", "data/center_1/custom_center_1"],
    )
    assert configure.exit_code == 0, configure.output

    center_config = yaml.safe_load((project_dir / "configs" / "center_1.yaml").read_text())
    assert center_config["input_data"]["path"] == "data/center_1/custom_center_1"

    verify = runner.invoke(sim.app, ["check", "--data"])
    assert verify.exit_code == 0, verify.output


def test_baseline_generate_delegates_to_pipeline_tool_with_default_output(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    init = runner.invoke(sim.app, ["init"])
    assert init.exit_code == 0, init.output
    calls = []

    def fake_generate_baseline(experiment_config, output_dir=None):
        calls.append((Path(experiment_config), Path(output_dir)))
        return {"lr_prefix": str(Path(output_dir) / "lr")}

    monkeypatch.setattr(baseline_tool, "generate_baseline", fake_generate_baseline)

    result = runner.invoke(sim.app, ["baseline", "generate"])

    assert result.exit_code == 0, result.output
    assert calls == [(tmp_path / "config.yaml", tmp_path / "results" / "baseline")]
    assert "lr_prefix" in result.output


def test_baseline_generate_delegates_to_pipeline_tool_with_custom_output(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    init = runner.invoke(sim.app, ["init"])
    assert init.exit_code == 0, init.output
    calls = []
    custom_output = tmp_path / "custom-baseline"

    def fake_generate_baseline(experiment_config, output_dir=None):
        calls.append((Path(experiment_config), Path(output_dir)))
        return {"merged_prefix": str(Path(output_dir) / "merged")}

    monkeypatch.setattr(baseline_tool, "generate_baseline", fake_generate_baseline)

    result = runner.invoke(sim.app, ["baseline", "generate", "--output", str(custom_output)])

    assert result.exit_code == 0, result.output
    assert calls == [(tmp_path / "config.yaml", custom_output)]
    assert "merged_prefix" in result.output


def test_evaluate_help_shows_core_evaluation_options():
    output = _help_output("evaluate")

    for flag in (
        "--baseline",
        "--report",
        "--no-plots",
        "--king",
        "--king-center-id",
        "--king-data-dir",
        "--qc-only",
        "--lr-only",
        "--king-only",
    ):
        assert flag in output


def test_evaluate_delegates_to_core_evaluation_with_project_defaults(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    init = runner.invoke(sim.app, ["init"])
    assert init.exit_code == 0, init.output
    (tmp_path / "results" / "baseline").mkdir(parents=True)
    calls = []

    def fake_run_evaluation(**kwargs):
        calls.append(kwargs)
        return {
            "report_path": str(kwargs["report_path"]),
            "qc_report_path": str(kwargs["results_dir"] / "qc_report.md"),
            "lr_report_path": str(kwargs["results_dir"] / "lr_report.md"),
            "king_report_path": None,
        }

    monkeypatch.setattr(evaluation_tool, "run_evaluation", fake_run_evaluation)

    result = runner.invoke(sim.app, ["evaluate"])

    assert result.exit_code == 0, result.output
    assert calls == [
        {
            "results_dir": tmp_path / "results",
            "baseline_dir": tmp_path / "results" / "baseline",
            "report_path": tmp_path / "results" / "evaluation_report.md",
            "no_plots": False,
            "king": False,
            "king_center_id": 1,
            "king_data_dir": None,
            "qc_only": False,
            "lr_only": False,
            "king_only": False,
        }
    ]
    assert "Evaluation completed." in result.output
    assert "report_path: results/evaluation_report.md" in result.output


def test_evaluate_delegates_to_core_evaluation_with_explicit_options(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    init = runner.invoke(sim.app, ["init"])
    assert init.exit_code == 0, init.output
    calls = []

    def fake_run_evaluation(**kwargs):
        calls.append(kwargs)
        return {
            "report_path": str(kwargs["report_path"]),
            "qc_report_path": None,
            "lr_report_path": None,
            "king_report_path": str(kwargs["results_dir"] / "king_report.md"),
        }

    monkeypatch.setattr(evaluation_tool, "run_evaluation", fake_run_evaluation)

    result = runner.invoke(
        sim.app,
        [
            "evaluate",
            "custom-results",
            "--baseline",
            "baseline-out",
            "--report",
            "reports/evaluation.md",
            "--no-plots",
            "--king",
            "--king-center-id",
            "2",
            "--king-data-dir",
            "data",
            "--king-only",
        ],
    )

    assert result.exit_code == 0, result.output
    assert calls == [
        {
            "results_dir": tmp_path / "custom-results",
            "baseline_dir": tmp_path / "baseline-out",
            "report_path": tmp_path / "reports" / "evaluation.md",
            "no_plots": True,
            "king": True,
            "king_center_id": 2,
            "king_data_dir": tmp_path / "data",
            "qc_only": False,
            "lr_only": False,
            "king_only": True,
        }
    ]
    assert "king_report_path: custom-results/king_report.md" in result.output


def test_evaluate_reports_missing_baseline_as_parameter_error(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    init = runner.invoke(sim.app, ["init"])
    assert init.exit_code == 0, init.output

    result = runner.invoke(sim.app, ["evaluate"])

    assert result.exit_code != 0
    assert "Baseline directory not found" in result.output


def test_run_invokes_flower_with_project_config(tmp_path, monkeypatch):
    project_dir = tmp_path / "tiny-study"
    project_dir.mkdir()
    (project_dir / "fedgwas.yaml").write_text(
        "\n".join(
            [
                "mode: simulation",
                "num_clients: 2",
                "config_dir: configs",
                "data_dir: data",
                "results_dir: results",
                "plink: auto",
            ]
        )
    )
    (project_dir / "configs").mkdir()
    monkeypatch.chdir(project_dir)
    monkeypatch.setattr(sim, "_run_checks", lambda selected=None: 0)

    calls = []

    def fake_run(cmd, cwd=None, check=False):
        calls.append((cmd, cwd, check))

        class Result:
            returncode = 7

        return Result()

    monkeypatch.setattr(sim.subprocess, "run", fake_run)

    result = runner.invoke(sim.app, ["run", "--rounds", "100", "--no-stream"])

    assert result.exit_code == 7
    assert calls == [
        (
            [
                "flwr",
                "run",
                ".",
                "local-simulation",
                "--run-config",
                'simulation=true num-server-rounds=100 config_path="configs"',
            ],
            project_dir,
            False,
        )
    ]


def test_run_exits_before_flower_when_checks_fail(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    init = runner.invoke(sim.app, ["init", "--from-example", "syn-tiny", "--no-prepare-data"])
    assert init.exit_code == 0, init.output

    calls = []

    def fake_run(cmd, cwd=None, check=False):
        calls.append((cmd, cwd, check))

        class Result:
            returncode = 0

        return Result()

    monkeypatch.setattr(sim.subprocess, "run", fake_run)

    result = runner.invoke(sim.app, ["run", "--rounds", "5"])

    assert result.exit_code == 1
    assert calls == []
    assert "Check Summary" in result.output
    assert "FAIL" in result.output
    assert "Next Steps" in result.output


def test_flower_apps_resolve_flat_and_legacy_config_layouts(tmp_path):
    configs_dir = tmp_path / "configs"
    configs_dir.mkdir()
    flat_center = configs_dir / "center_1.yaml"
    flat_server = configs_dir / "server.yaml"
    flat_center.write_text("input_data: {}\n")
    flat_server.write_text("output: {}\n")

    assert _resolve_simulation_client_config_path(configs_dir, 1) == flat_center
    assert _resolve_server_config_file(configs_dir) == flat_server

    flat_center.unlink()
    flat_server.unlink()
    legacy_center = configs_dir / "center_1" / "config.yaml"
    legacy_server = configs_dir / "server" / "config.yaml"
    legacy_center.parent.mkdir()
    legacy_server.parent.mkdir()
    legacy_center.write_text("input_data: {}\n")
    legacy_server.write_text("output: {}\n")

    assert _resolve_simulation_client_config_path(configs_dir, 1) == legacy_center
    assert _resolve_server_config_file(configs_dir) == legacy_server
