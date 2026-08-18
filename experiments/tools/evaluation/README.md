# Evaluation Tools

This folder contains post-evaluation utilities for federated GWAS experiments, split by stage:

- `qc/`: QC exclusion list agreement
- `king/`: KING evaluation utilities (pair checks, accumulator comparison, debug tools)
- `lr/`: LR agreement, significance coverage, and plots

## Quick Start

All-in-one (QC + LR, optional KING):

```bash
python experiments/tools/evaluation/evaluate_all.py \
  experiments/correctness/tiny_even/results \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline \
  --king --king-center-id 1
```

Single-stage options:

```bash
# QC only
python experiments/tools/evaluation/evaluate_all.py \
  experiments/correctness/tiny_even/results \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline \
  --qc-only

# LR only
python experiments/tools/evaluation/evaluate_all.py \
  experiments/correctness/tiny_even/results \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline \
  --lr-only

# KING only
python experiments/tools/evaluation/evaluate_all.py \
  experiments/correctness/tiny_even/results \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline \
  --king-only --king-center-id 1
```

QC agreement:

```bash
python experiments/tools/evaluation/qc/qc_evaluator.py \
  experiments/correctness/tiny_even/results \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline
```

LR agreement:

```bash
python experiments/tools/evaluation/lr/lr_evaluator.py \
  experiments/correctness/tiny_even/results \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline
```

KING evaluation (accumulator vs baseline):

```bash
python experiments/tools/evaluation/king/compare_king_from_accumulator.py \
  --results-dir experiments/real_world/1000genomes/results \
  --baseline-dir experiments/real_world/1000genomes/data/centralized_baseline \
  --center-id 1
```

## Outputs

- `qc_report.md`: QC agreement summary
- `lr_report.md`: LR correlation + coverage summary
- `evaluation_report.md`: Combined summary (from `evaluate_all.py`)
- `king_report.md`: KING accumulator comparison summary (when using `--king`)
- `manhattan_plot.png`: Manhattan plot from LR evaluator
- `lr_correlation_plots.png`: Local/global LR correlation plots

## Notes

- KING evaluation is intentionally separate from QC/LR evaluation.
- Use the tools in `king/` for cross-client pair checks, accumulator validation, and formula debugging.
