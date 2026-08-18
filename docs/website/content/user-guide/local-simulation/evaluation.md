---
title: Evaluation and Results
slug: /local-simulation/evaluation
---

# Evaluation and Results

After a local simulation run, use the FedGWAS evaluation tools to generate a
centralized PLINK baseline, compare federated outputs, and collect lightweight
run metadata. The `fedgwas-sim` commands wrap the same evaluation modules used
by repository tools.

## Result Layout

A generated CLI project writes results under the configured `results_dir`,
usually `results/`:

```text
results/
  server/
    logs/
    intermediate/
  center_1/
    logs/
    intermediate/
  center_2/
    logs/
    intermediate/
```

Important artifacts include:

| Path | Purpose |
| --- | --- |
| `results/server/logs/` | Server stage progress and aggregation logs |
| `results/center_*/logs/` | Per-center client logs and PLINK command details |
| `results/server/intermediate/` | Server-side temporary merged data and analysis outputs |
| `results/center_*/intermediate/run_<id>/` | Client-side run-scoped intermediate files |
| `results/evaluation_report.md` | Combined post-run evaluation summary |
| `results/run_summary.json` | Machine-readable run summary from result collection |
| `results/run_summary.md` | Human-readable run summary |

Retention settings can prune large intermediate files after successful
completion while preserving logs and summary artifacts.

## Generate A Centralized Baseline

Generate a centralized baseline after data preparation and before or after the federated run:

```bash
fedgwas-sim baseline generate
```

By default, the baseline is written under `results/baseline`. You can choose a different location:

```bash
fedgwas-sim baseline generate --output data/centralized_baseline
```

Use the same center data and thresholds for baseline generation and federated evaluation. If you edit center configs after generating the baseline, regenerate the baseline.

## Evaluate Federated Outputs

Run QC and LR evaluation with the default baseline path:

```bash
fedgwas-sim evaluate
```

Evaluate against an explicit baseline and include KING comparison:

```bash
fedgwas-sim evaluate results --baseline data/centralized_baseline --king
```

Write the combined report to a chosen path:

```bash
fedgwas-sim evaluate results \
  --baseline data/centralized_baseline \
  --report results/evaluation_report.md \
  --king
```

Evaluation can run selected stages only:

```bash
fedgwas-sim evaluate results --baseline results/baseline --qc-only
fedgwas-sim evaluate results --baseline results/baseline --lr-only --no-plots
fedgwas-sim evaluate results --baseline results/baseline --king-only --king-center-id 1
```

By default, `evaluate` runs QC and LR evaluation. Add `--king` when the run produced KING accumulator artifacts and you want kinship comparison included.

## Reports And Plots

Evaluation may write:

| File | Contents |
| --- | --- |
| `qc_report.md` | QC agreement checks between federated and centralized outputs |
| `lr_report.md` | Logistic regression agreement summary and generated plot references |
| `king_report.md` | KING accumulator comparison when requested |
| `evaluation_report.md` | Combined report with links to stage-specific reports |

Use `--no-plots` when running in a minimal environment or when only tabular agreement metrics are needed.

## Collect Run Metadata

Collect a summary of logs, intermediate file counts, and optional timing files:

```bash
fedgwas-sim results collect --label tiny_run
```

Include GNU `time` output when available:

```bash
fedgwas-sim results collect \
  --time-file results/server_app_time.txt \
  --label tiny_run
```

The command writes:

```text
results/run_summary.json
results/run_summary.md
```

Use these summaries for experiment records, release checks, or lightweight performance notes.

**Success Indicators:**

For a healthy run and evaluation:

- `fedgwas-sim run` exits with status code 0.
- Server and client logs are present.
- Evaluation writes `evaluation_report.md`.
- QC reports do not show unexpected missing outputs.
- LR reports show comparable federated and centralized association outputs for
  the selected scenario.
- KING reports are present when `--king` or `--king-only` is requested.

If evaluation fails because files are missing, confirm the selected results directory, baseline directory, and center output paths match the study that was actually run.

## [Optional] Repository Evaluation Flow

For repository experiments, use the repository tools and experiment paths:

```bash
python experiments/tools/generate_baseline.py \
  experiments/correctness/tiny_even/config.yaml

python experiments/tools/evaluation/evaluate_all.py \
  experiments/correctness/tiny_even/results_2 \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline \
  --king

python experiments/tools/collect_run_metrics.py \
  experiments/correctness/tiny_even/results_2
```

Do not mix these repository paths with a standalone `fedgwas-sim` study unless you intentionally copied all configs, data, and outputs.
