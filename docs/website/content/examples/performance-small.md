---
title: Small Dataset for Performance Measurement
---

# Performance Small

## Purpose

Measure wall-clock time, stage duration, resource use, and output retention behavior for a two-client small-scale run.

## Data

The small layout is stored under `experiments/performance/small_even/`. Generate synthetic PLINK data locally when needed:

```bash
python pipeline/simulation/simulated_data/generate_synthetic_data.py \
  --scale small \
  --partition-strategy even \
  --seed 42 \
  --output-dir experiments/performance/small_even/data
```

## Command

For cluster benchmarking, use the scripts and guide under `cluster_deployment/`. For local deployment, point SuperNodes at:

```text
experiments/performance/small_even/configs/center_1/config.yaml
experiments/performance/small_even/configs/center_2/config.yaml
```

After the run:

```bash
python experiments/tools/collect_run_metrics.py \
  experiments/performance/small_even/results \
  --label small_even
```

## Expected output

- `stage_metrics.csv` and `client_metrics.csv` under the results root.
- `run_summary.md` and `run_summary.json` from `collect_run_metrics.py`.
- Retention manifest when `retention.auto_apply_on_complete` is enabled.
