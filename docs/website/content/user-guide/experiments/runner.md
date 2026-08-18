---
title: Experiment Operations
slug: /experiment-operations
---

# Experiment Operations

The committed repository does not include a Python experiment runner module. Experiments are run with Flower CLI commands, the cluster deployment scripts, and the evaluation tools under `experiments/tools/`.

## Local simulation

Use the default tiny correctness config:

```bash
flwr run . local-simulation --stream
```

Override the run config explicitly when needed:

```bash
flwr run . local-simulation --stream --run-config \
  'simulation=true num-server-rounds=100 config_path="experiments/correctness/tiny_even/configs"'
```

## Local deployment

Start one SuperLink, one SuperNode per client, then run:

```bash
flwr run . local-deployment --stream
```

`pyproject.toml` must point `local-deployment` at the SuperLink Exec API:

```toml
[tool.flwr.federations.local-deployment]
address = "127.0.0.1:9093"
insecure = true
options.num-supernodes = 2
```

## Cluster deployment

For three-node runs, use the scripts under `cluster_deployment/scripts/` and the guide at `cluster_deployment/docs/CLUSTER_USER_GUIDE.md`. The cluster scripts wrap setup, data verification, SuperLink/SuperNode startup, run execution, status checks, and log collection.

## Post-run evaluation

Generate a centralized baseline:

```bash
python experiments/tools/generate_baseline.py \
  experiments/correctness/tiny_even/config.yaml
```

Compare federated outputs with the baseline:

```bash
python experiments/tools/evaluation/evaluate_all.py \
  experiments/correctness/tiny_even/results_2 \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline \
  --king
```

Collect monitoring outputs:

```bash
python experiments/tools/collect_run_metrics.py \
  experiments/correctness/tiny_even/results_2
```

## Success indicators

- SuperLink starts the Fleet API.
- SuperNodes register with the federation.
- `flwr run` starts a run successfully.
- Server logs are written under the configured `server/logs/` directory.
- The pipeline progresses through initialization, federated QC, relatedness screening, and association screening.
