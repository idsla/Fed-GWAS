---
title: Tiny Synthetic Dataset
---

# Tiny Synthetic Dataset

## Purpose

Validate that federated QC, relatedness screening, and association screening outputs agree with a centralized PLINK baseline.

## Data

Generate synthetic data if missing:

```bash
python pipeline/simulation/simulated_data/generate_synthetic_data.py \
  --scale tiny \
  --partition-strategy even \
  --seed 42 \
  --output-dir experiments/correctness/tiny_even/data
```

## Command

```bash
python experiments/tools/generate_baseline.py \
  experiments/correctness/tiny_even/config.yaml

flwr run . local-simulation --stream --run-config \
  'simulation=true num-server-rounds=100 config_path="experiments/correctness/tiny_even/configs"'

python experiments/tools/evaluation/evaluate_all.py \
  experiments/correctness/tiny_even/results_2 \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline \
  --king
```

## Expected output

- Server and client logs under `experiments/correctness/tiny_even/results_2/`.
- QC and LR agreement reports from the evaluation tools.
- Successful runs should reach the `done` stage and meet the thresholds documented in `experiments/correctness/tiny_even/README.md`.
