from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline.evaluation import evaluate_all


class FakeQCEvaluator:
    calls: list[tuple[Path, Path, Path | None]] = []

    def __init__(self, federated_results_dir: str, baseline_results_dir: str):
        self.results_dir = Path(federated_results_dir)
        self.baseline_dir = Path(baseline_results_dir)

    def run(self, report_path=None):
        self.__class__.calls.append((self.results_dir, self.baseline_dir, report_path))
        if report_path is not None:
            Path(report_path).write_text("# QC\n", encoding="utf-8")
        return {"f1_score": 1.0, "precision": 1.0, "recall": 1.0}


class FakeLREvaluator:
    calls: list[tuple[Path, Path, Path | None, bool]] = []

    def __init__(self, federated_results_dir: str, baseline_results_dir: str):
        self.results_dir = Path(federated_results_dir)
        self.baseline_dir = Path(baseline_results_dir)

    def run(self, report_path=None, make_plots=True):
        self.__class__.calls.append((self.results_dir, self.baseline_dir, report_path, make_plots))
        if report_path is not None:
            Path(report_path).write_text("# LR\n", encoding="utf-8")
        return {"local": {"pearson_r": 1.0, "p_value": 0.0, "mse": 0.0, "mae": 0.0, "n_snps": 10}}


def _reset_fakes() -> None:
    FakeQCEvaluator.calls = []
    FakeLREvaluator.calls = []


@pytest.fixture(autouse=True)
def patch_evaluators(monkeypatch):
    _reset_fakes()
    king_calls = []

    def fake_king(**kwargs):
        king_calls.append(kwargs)
        report_path = kwargs.get("report_path")
        if report_path is not None:
            Path(report_path).write_text("# KING\n", encoding="utf-8")
        return {"total_pairs": 1, "mapped_pairs": 1}

    monkeypatch.setattr(evaluate_all, "QCEvaluator", FakeQCEvaluator)
    monkeypatch.setattr(evaluate_all, "LREvaluator", FakeLREvaluator)
    monkeypatch.setattr(evaluate_all, "compare_from_accumulator", fake_king)
    yield king_calls


def test_run_evaluation_runs_qc_and_lr_by_default(tmp_path, patch_evaluators):
    results_dir = tmp_path / "results"
    baseline_dir = tmp_path / "baseline"
    results_dir.mkdir()
    baseline_dir.mkdir()

    result = evaluate_all.run_evaluation(results_dir, baseline_dir, report_path=results_dir / "summary.md")

    assert FakeQCEvaluator.calls == [(results_dir, baseline_dir, results_dir / "qc_report.md")]
    assert FakeLREvaluator.calls == [(results_dir, baseline_dir, results_dir / "lr_report.md", True)]
    assert patch_evaluators == []
    assert result["qc"] is not None
    assert result["lr"] is not None
    assert result["king"] is None
    assert (results_dir / "summary.md").is_file()


def test_run_evaluation_runs_king_when_requested(tmp_path, patch_evaluators):
    results_dir = tmp_path / "results"
    baseline_dir = tmp_path / "baseline"
    data_dir = tmp_path / "data"
    results_dir.mkdir()
    baseline_dir.mkdir()

    result = evaluate_all.run_evaluation(
        results_dir,
        baseline_dir,
        king=True,
        king_center_id=2,
        king_data_dir=data_dir,
    )

    assert result["king"] == {"total_pairs": 1, "mapped_pairs": 1}
    assert patch_evaluators[0]["results_dir"] == results_dir
    assert patch_evaluators[0]["baseline_dir"] == baseline_dir
    assert patch_evaluators[0]["center_id"] == 2
    assert patch_evaluators[0]["data_dir"] == data_dir
    assert patch_evaluators[0]["report_path"] == results_dir / "king_report.md"


def test_run_evaluation_honors_only_modes(tmp_path, patch_evaluators):
    results_dir = tmp_path / "results"
    baseline_dir = tmp_path / "baseline"
    results_dir.mkdir()
    baseline_dir.mkdir()

    qc_result = evaluate_all.run_evaluation(results_dir, baseline_dir, qc_only=True)
    assert qc_result["qc"] is not None
    assert qc_result["lr"] is None
    assert qc_result["king"] is None
    assert len(FakeQCEvaluator.calls) == 1
    assert FakeLREvaluator.calls == []

    _reset_fakes()
    lr_result = evaluate_all.run_evaluation(results_dir, baseline_dir, lr_only=True, no_plots=True)
    assert lr_result["qc"] is None
    assert lr_result["lr"] is not None
    assert lr_result["king"] is None
    assert FakeQCEvaluator.calls == []
    assert FakeLREvaluator.calls == [(results_dir, baseline_dir, results_dir / "lr_report.md", False)]

    _reset_fakes()
    king_result = evaluate_all.run_evaluation(results_dir, baseline_dir, king_only=True)
    assert king_result["qc"] is None
    assert king_result["lr"] is None
    assert king_result["king"] is not None


def test_run_evaluation_rejects_multiple_only_modes(tmp_path):
    results_dir = tmp_path / "results"
    baseline_dir = tmp_path / "baseline"
    results_dir.mkdir()
    baseline_dir.mkdir()

    with pytest.raises(ValueError, match="Only one"):
        evaluate_all.run_evaluation(results_dir, baseline_dir, qc_only=True, lr_only=True)
