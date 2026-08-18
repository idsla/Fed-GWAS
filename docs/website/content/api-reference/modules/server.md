---
title: Server Modules
---

# Server Modules

The server entry point is `pipeline/server_app.py`. `server_fn(context)` reads Flower run configuration and creates `FederatedGWASStrategy`.

## Output directory selection

`pipeline/server_app.py` requires the server output directory to be declared explicitly. It resolves the directory in this order:

1. `server_output_dir` from Flower `run_config`.
2. `configs/server/config.yaml` under the configured `config_path`.

For the default tiny correctness experiment, the server config is:

```text
experiments/correctness/tiny_even/configs/server/config.yaml
```

It declares:

```yaml
output:
  intermediate_dir: "experiments/correctness/tiny_even/results_2/server/intermediate"
  log_dir: "experiments/correctness/tiny_even/results_2/server/logs"
```

## Strategy responsibilities

`FederatedGWASStrategy` extends Flower `Strategy` and implements:

- `initialize_parameters`: returns empty parameters because GWAS does not use a conventional global model.
- `configure_fit`: samples all available clients and sends the current stage configuration.
- `aggregate_fit`: parses client responses, calls aggregators, and advances the stage machine.
- `configure_evaluate`, `aggregate_evaluate`, and `evaluate`: return empty evaluation results because traditional model evaluation is not used.

## Aggregators

| File | Purpose |
| --- | --- |
| `aggregator_king.py` | Reconstructs client chunks, merges PLINK data, runs PLINK `--het`, and runs PLINK2 `--make-king-table`. |
| `aggregator_lr.py` | Reconstructs LR chunks, merges PLINK data, runs `--logistic`, and parses p-values. |
| `prg_masking.py` | Tracks public keys and provides curve parameters for encrypted client-to-client stages. |

## Strict stage behavior

`strategy_strict.py` keeps stage configs Scalar-only and does not advance through encrypted relay stages when no valid client results arrive. It retries `sync`, `global_qc`, KING, and LR drain stages instead of using the older simulation workaround strategy.
