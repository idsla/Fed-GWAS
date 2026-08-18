---
title: API Reference
---

# API Reference

This page summarizes the runtime knobs users most often need when launching FedGWAS through Flower.

## CLI entry point

Most local simulations should start from the packaged CLI:

```bash
fedgwas-sim --help
```

The command creates a local simulation project, configures examples and
presets, runs readiness checks, launches Flower, and evaluates or collects
outputs. See the main [fedgwas-sim CLI](/user-guide/cli-simulation) documentation for
the complete command reference.

## Flower run config

| Key | Type | Default | Purpose |
| --- | --- | --- | --- |
| `simulation` | boolean | `true` | Selects simulation-aware client behavior and output layout. |
| `num-server-rounds` | integer | `300` | Maximum Flower rounds. The strategy stops earlier when it reaches `done`. |
| `config_path` | string | `experiments/correctness/tiny_even/configs` | Directory containing `server/config.yaml` and `center_x/config.yaml` files. |
| `default-partition-by` | string | `samples` | Default partitioning mode for client initialization. |
| `server_output_dir` | string | unset | Optional explicit server output directory override. |
| `chunk_size` / `server_chunk_size` | integer | experiment `server.chunk_size` or `1000` | Chunk size sent during `init_chunks` and `init_chunks_lr`. |
| `lr_pad_to` | integer | experiment `server.lr_pad_to` or `0` | Optional padding target for privacy-preserving tokens in local association filtering. |
| `phenotype_fix_missing_to_case` | boolean | `false` | Creates a derived PLINK prefix with missing phenotypes recoded for toy runs. |

Example:

```bash
flwr run . local-simulation --stream --run-config \
  'simulation=true num-server-rounds=100 config_path="experiments/correctness/tiny_even/configs"'
```

## Center config files

Each `center_x/config.yaml` contains:

| Section | Purpose |
| --- | --- |
| `input_data` | PLINK prefix and input type. |
| `output` | Per-center intermediate and log directories. |
| `parameters` | Chunk sizes, run tags, and anonymization offsets. |
| `thresholds` | Federated QC, relatedness screening, and association screening thresholds. |
| `flower` | Deployment address and nominal round count. |
| `participation` | Per-stage opt-in flags. |

## Server config file

The server config is required at `<config_path>/server/config.yaml` unless `server_output_dir` is provided. It must define `output.log_dir` or `output.intermediate_dir`.

## Main Python entry points

| File | Role |
| --- | --- |
| `pipeline/client_app.py` | Flower `ClientApp`; creates `FedLRClient` and routes stages. |
| `pipeline/server_app.py` | Flower `ServerApp`; resolves config and creates the strategy. |
| `pipeline/cli/sim.py` | Typer application behind the `fedgwas-sim` command. |
| `pipeline/server/strategy_strict.py` | Stage machine and encrypted relay orchestration. |
| `pipeline/clients/base_client.py` | PLINK chunking, anonymization, phenotype helpers. |
| `pipeline/server/aggregator_king.py` | Server-side KING chunk reconstruction and PLINK2 execution. |
| `pipeline/server/aggregator_lr.py` | Server-side LR chunk reconstruction and PLINK execution. |
