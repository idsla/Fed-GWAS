---
title: Simulation Mode
slug: /simulation
---

# Simulation Mode

Simulation mode uses Flower `local-simulation` to run multiple virtual clients
on one machine. It is the fastest way to validate the pipeline and reproduce
the documentation examples after installation.

## Recommended CLI workflow

For installed FedGWAS packages, use `fedgwas-sim` to create a self-contained
study directory, prepare data, validate configuration, run Flower, and collect
results:

```bash
python -m pip install FedGWAS

mkdir my_study
cd my_study

fedgwas-sim init
fedgwas-sim setup-experiment syn-tiny --seed 42
fedgwas-sim check
fedgwas-sim run --rounds 100
fedgwas-sim baseline generate --output data/centralized_baseline
fedgwas-sim evaluate results --baseline data/centralized_baseline --king
fedgwas-sim results collect --label tiny_run
```

See [fedgwas-sim CLI](/user-guide/cli-simulation) for the full command reference,
available presets, project layout, reset behavior, and evaluation options.

## Run the default simulation

Use this repository workflow when developing FedGWAS itself or when running the
checked-in experiment configurations directly:

```bash
flwr run . local-simulation --stream
```

The default settings are defined in `pyproject.toml`:

```toml
[tool.flwr.app.config]
simulation = true
num-server-rounds = 300
config_path = "experiments/correctness/tiny_even/configs"

[tool.flwr.federations.local-simulation]
options.num-supernodes = 2
```

## Scenario layout

Each center needs a local PLINK dataset and a center-specific configuration file:

```text
experiments/correctness/tiny_even/
  configs/
    server/config.yaml
    center_1/config.yaml
    center_2/config.yaml
  data/
    tiny/
      center_1/tiny_center_1.bed
      center_1/tiny_center_1.bim
      center_1/tiny_center_1.fam
      center_2/tiny_center_2.bed
      center_2/tiny_center_2.bim
      center_2/tiny_center_2.fam
  results_2/
```

The committed center configs currently write to `results_2/`; check each center `output.log_dir` and `output.intermediate_dir` if you switch experiments.

## Client initialization

`client_fn` reads the Flower `Context`, maps the partition id to `center_x/config.yaml`, and creates a `FedLRClient`.

On first initialization, the client:

1. Loads YAML settings with `DataLoader`.
2. Binds or converts input data to a PLINK prefix.
3. Stores config records in Flower client state.
4. Creates the configured log and intermediate directories.
5. Initializes the per-client logger.

On later rounds, the client restores the same configuration from Flower state instead of reading the YAML file again.

## Outputs

- Client logs: `output.log_dir` from each center config.
- Client intermediates: `output.intermediate_dir` from each center config.
- Server outputs: the directory declared in `configs/server/config.yaml`, usually under `experiments/correctness/tiny_even/results_2/server/` for the shipped tiny config.
