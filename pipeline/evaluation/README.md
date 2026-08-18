# Evaluation Tools

This folder contains post-evaluation utilities for federated GWAS experiments, split by stage:

- `qc/`: QC exclusion list agreement
- `king/`: KING evaluation utilities (pair checks, accumulator comparison, debug tools)
- `lr/`: LR agreement, significance coverage, and plots

## Quick Start

All-in-one (QC + LR, optional KING):

```bash
python -m pipeline.evaluation.evaluate_all \
  experiments/correctness/tiny_even/results \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline \
  --king --king-center-id 1
```

Single-stage options:

```bash
# QC only
python -m pipeline.evaluation.evaluate_all \
  experiments/correctness/tiny_even/results \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline \
  --qc-only

# LR only
python -m pipeline.evaluation.evaluate_all \
  experiments/correctness/tiny_even/results \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline \
  --lr-only

# KING only
python -m pipeline.evaluation.evaluate_all \
  experiments/correctness/tiny_even/results \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline \
  --king-only --king-center-id 1
```

QC agreement:

```bash
python -m pipeline.evaluation.qc.qc_evaluator \
  experiments/correctness/tiny_even/results \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline
```

LR agreement:

```bash
python -m pipeline.evaluation.lr.lr_evaluator \
  experiments/correctness/tiny_even/results \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline
```

KING evaluation (accumulator vs baseline):

```bash
python -m pipeline.evaluation.king.compare_king_from_accumulator \
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
