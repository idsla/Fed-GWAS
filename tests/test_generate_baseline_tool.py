from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline.tools.generate_baseline import BaselineGenerator


def _write_config(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")


def test_baseline_generator_resolves_relative_data_dir_from_config_path(tmp_path, monkeypatch):
    project_dir = tmp_path / "study"
    config_path = project_dir / "config.yaml"
    _write_config(
        config_path,
        {
            "data": {"data_dir": "data"},
            "clients": {"thresholds": {"maf_threshold": 0.02}},
        },
    )
    monkeypatch.setattr(BaselineGenerator, "_find_plink", lambda self: "plink")

    generator = BaselineGenerator(str(config_path), str(project_dir / "results" / "baseline"))

    assert generator.config["data"]["data_dir"] == str(project_dir / "data")


def test_baseline_generator_reads_thresholds_from_clients_thresholds(tmp_path, monkeypatch):
    config_path = tmp_path / "config.yaml"
    _write_config(
        config_path,
        {
            "data": {"data_dir": "data"},
            "clients": {
                "thresholds": {
                    "maf_threshold": 0.03,
                    "missing_threshold": 0.07,
                    "hwe_threshold": 1e-5,
                }
            },
        },
    )
    monkeypatch.setattr(BaselineGenerator, "_find_plink", lambda self: "plink")

    generator = BaselineGenerator(str(config_path), str(tmp_path / "baseline"))

    assert generator._thresholds()["maf_threshold"] == 0.03
    assert generator._thresholds()["missing_threshold"] == 0.07
    assert generator._thresholds()["hwe_threshold"] == 1e-5


def test_baseline_qc_uses_generator_plink_binary(tmp_path, monkeypatch):
    config_path = tmp_path / "config.yaml"
    _write_config(config_path, {"data": {"data_dir": "data"}})
    monkeypatch.setattr(BaselineGenerator, "_find_plink", lambda self: "repo-plink")

    generator = BaselineGenerator(str(config_path), str(tmp_path / "baseline"))
    input_prefix = tmp_path / "merged"
    (tmp_path / "merged.bim").write_text("1 rs1 0 1 A G\n", encoding="utf-8")
    captured: dict[str, str | None] = {}

    def fake_counts(plink_prefix, client_id, log_dir="logs", plink_binary=None):
        captured["counts_plink"] = plink_binary
        return np.array([[10, 2, 8]], dtype=np.int64)

    def fake_missing(plink_prefix, client_id, log_dir="logs", plink_binary=None):
        captured["missing_plink"] = plink_binary
        return np.array([[20, 0]], dtype=np.int64)

    monkeypatch.setattr("pipeline.clients.local_qc.compute_genotype_counts", fake_counts)
    monkeypatch.setattr("pipeline.clients.local_qc.compute_missingness_counts", fake_missing)
    monkeypatch.setattr("pipeline.clients.client_qc_aggregator._compute_exclusion_list", lambda *args: set())
    monkeypatch.setattr(generator, "_run_plink", lambda cmd, error_message: None)
    monkeypatch.setattr("pipeline.tools.generate_baseline.subprocess.run", lambda *args, **kwargs: None)

    generator.run_qc(str(input_prefix))

    assert captured == {"counts_plink": "repo-plink", "missing_plink": "repo-plink"}
