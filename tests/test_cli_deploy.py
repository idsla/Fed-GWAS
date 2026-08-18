from __future__ import annotations

import sys
from pathlib import Path

from typer.testing import CliRunner

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline.cli import deploy


runner = CliRunner()


def test_server_start_builds_superlink_command(monkeypatch):
    calls = []

    def fake_run(cmd, cwd=None, check=False):
        calls.append((cmd, cwd, check))

        class Result:
            returncode = 0

        return Result()

    monkeypatch.setattr(deploy.subprocess, "run", fake_run)

    result = runner.invoke(
        deploy.app,
        ["server", "start", "--host", "0.0.0.0", "--fleet-port", "9192", "--exec-port", "9193"],
    )

    assert result.exit_code == 0, result.output
    assert calls == [
        (
            [
                "flower-superlink",
                "--insecure",
                "--fleet-api-address",
                "0.0.0.0:9192",
                "--exec-api-address",
                "0.0.0.0:9193",
            ],
            Path.cwd(),
            False,
        )
    ]


def test_client_start_builds_supernode_command_with_center_defaults(monkeypatch, tmp_path):
    config = tmp_path / "center_1.yaml"
    config.write_text("input_data: {}\n", encoding="utf-8")
    calls = []

    def fake_run(cmd, cwd=None, check=False):
        calls.append((cmd, cwd, check))

        class Result:
            returncode = 0

        return Result()

    monkeypatch.setattr(deploy.subprocess, "run", fake_run)

    result = runner.invoke(
        deploy.app,
        [
            "client",
            "start",
            "--server",
            "192.168.1.88",
            "--center-id",
            "1",
            "--num-centers",
            "2",
            "--config",
            str(config),
        ],
    )

    assert result.exit_code == 0, result.output
    assert calls == [
        (
            [
                "flower-supernode",
                "--insecure",
                "--superlink",
                "192.168.1.88:9092",
                "--clientappio-api-address",
                "0.0.0.0:9094",
                "--node-config",
                f'partition-id=0 num-partitions=2 config-file="{config.resolve()}"',
            ],
            Path.cwd(),
            False,
        )
    ]


def test_server_run_uses_scale_config_path_and_restores_pyproject(monkeypatch, tmp_path):
    pyproject = tmp_path / "pyproject.toml"
    original = """
[tool.flwr.federations.local-deployment]
address = "127.0.0.1:9093"
insecure = true

[tool.flwr.federations.cluster-deployment]
address = "127.0.0.1:9093"
insecure = true
"""
    pyproject.write_text(original, encoding="utf-8")
    config_dir = tmp_path / "experiments" / "correctness" / "tiny_even" / "configs"
    config_dir.mkdir(parents=True)
    monkeypatch.chdir(tmp_path)
    calls = []

    def fake_run(cmd, cwd=None, check=False):
        calls.append((cmd, cwd, check, pyproject.read_text(encoding="utf-8")))

        class Result:
            returncode = 3

        return Result()

    monkeypatch.setattr(deploy.subprocess, "run", fake_run)

    result = runner.invoke(
        deploy.app,
        ["server", "run", "--server", "192.168.1.88", "--rounds", "20", "--scale", "tiny", "--no-stream"],
    )

    assert result.exit_code == 3
    assert pyproject.read_text(encoding="utf-8") == original
    assert calls[0][:3] == (
        [
            "flwr",
            "run",
            ".",
            "local-deployment",
            "--run-config",
            'simulation=false num-server-rounds=20 config_path="experiments/correctness/tiny_even/configs"',
        ],
        tmp_path,
        False,
    )
    assert 'address = "192.168.1.88:9093"' in calls[0][3]


def test_client_check_reports_missing_config(tmp_path):
    result = runner.invoke(
        deploy.app,
        [
            "client",
            "check",
            "--server",
            "192.168.1.88",
            "--center-id",
            "1",
            "--config",
            str(tmp_path / "missing.yaml"),
        ],
    )

    assert result.exit_code == 1
    assert "config file" in result.output
