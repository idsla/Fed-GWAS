---
title: Real Data 1000 Genomes
---

# Real Data 1000 Genomes

## Purpose

Run FedGWAS on a real genotype workflow using a prepared 1000 Genomes chr22 subset and generated phenotypes.

## Data

Follow the experiment-specific instructions in `experiments/real_world/1000genomes/README.md`. The raw and generated PLINK data are not committed.

## Command

Use the real-world config path:

```bash
flwr run . local-simulation --stream --run-config \
  'config_path="experiments/real_world/1000genomes/configs"'
```

For deployment or cluster runs, point each SuperNode at the corresponding `center_x/config.yaml` under the same config directory.

## Expected output

- Federated QC, relatedness screening, and association screening outputs under the configured real-world results directory.
- Optional centralized baseline and evaluation reports from `experiments/tools/evaluation/`.
- Manuscript figures and tables under `experiments/real_world/1000genomes/manuscript/` when generated.
