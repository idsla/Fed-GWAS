---
title: Experiments
slug: /experiments
---

# Experiments

`experiments/` stores correctness, performance, and real-world study layouts. Each committed experiment keeps configuration and evaluation tooling in git; PLINK inputs and run outputs are generated locally and ignored.

```text
experiments/
  correctness/
    tiny_even/
      README.md
      config.yaml
      configs/
        server/config.yaml
        center_1/config.yaml
        center_2/config.yaml
      data/
      results/
  performance/
    scales.yaml
    small_even/
    medium_even/
  real_world/
    phenotype_generation/
    1000genomes/
  tools/
    generate_baseline.py
    collect_run_metrics.py
    evaluation/
```

## Scenario catalog

| Scenario | Purpose | Clients | Status |
| --- | --- | --- | --- |
| `correctness/tiny_even` | Default two-center correctness validation | 2 | Configured |
| `performance/small_even` | Small scalability run | 2 | Configured |
| `performance/medium_even` | Medium scalability run | 2 | Configured |
| `real_world/1000genomes` | 1000 Genomes chr22 application experiment | 2 | Configured |

## Generate synthetic data

```bash
python pipeline/simulation/simulated_data/generate_synthetic_data.py \
  --scale tiny \
  --partition-strategy even \
  --seed 42 \
  --output-dir experiments/correctness/tiny_even/data
```

The synthetic data generator lives in the local-only `pipeline/simulation/` tree. It is referenced by the experiment docs, but generated PLINK data is not committed.

## Run a scenario

```bash
flwr run . local-simulation --stream --run-config \
  'simulation=true num-server-rounds=100 config_path="experiments/correctness/tiny_even/configs"'
```

## Result layout

Scenario results are organized under each scenario's configured results directory. The shipped tiny configs currently use `results_2/`:

- `center_x/intermediate/`: client intermediate files.
- `center_x/logs/`: client logs and PLINK outputs.
- `server/logs/`: server logs and stage metrics.
- `server/intermediate/`: temporary merged data, KING outputs, and LR outputs.
- Root-level merged metrics and retention manifests when monitoring and retention are enabled.
