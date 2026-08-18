---
title: fedgwas-sim CLI
slug: /cli-simulation
---

# fedgwas-sim CLI Usage Guide

This guide will give more usage instructions on `fedgwas-sim`, the command line interface for creating and running local
FedGWAS simulation projects. It is the recommended entry point for users who install FedGWAS from PyPI and want to run examples outside the repository
checkout.

## Quickstart

After installing FedGWAS, confirm that the `fedgwas-sim` console script is available:

```bash
fedgwas-sim --help
```
This command lists the available subcommands and options and details help for each subcommand.

The CLI manages a study directory, writes Flower project files, prepares or verifies data, launches `flwr run . local-simulation`, and collects or evaluates outputs after the run. The recommended workflow starts from a user-owned study directory:

```bash
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

## Subcommands

This section will detail the usage of the each subcommand, but for a quick reference, here is a table of the main subcommands and their purposes:

| Command | Purpose |
| --- | --- |
| `fedgwas-sim init` | Create the base project layout in the current directory. |
| `fedgwas-sim setup-experiment <preset>` | Configure a runnable synthetic or real-world preset. |
| `fedgwas-sim prepare-data` | Generate missing synthetic data or run the generated real-world preparation script. |
| `fedgwas-sim check` | Run readiness checks before a simulation. |
| `fedgwas-sim run` | Launch the local Flower simulation after readiness checks pass. |
| `fedgwas-sim baseline generate` | Generate a centralized baseline for later comparison. |
| `fedgwas-sim evaluate [results_dir]` | Compare federated outputs with centralized baseline outputs. |
| `fedgwas-sim results collect` | Collect run metadata and timing summaries. |
| `fedgwas-sim summarize data` | Inspect prepared data under a data directory. |
| `fedgwas-sim summarize experiment` | Inspect simulation project metadata and readiness signals. |
| `fedgwas-sim data configure` | Point a center config at a user-supplied PLINK prefix. |
| `fedgwas-sim reset` / `clear` | Remove CLI-managed generated files from a simulation project. |

For each subcommnd, run `fedgwas-sim <subcommand> --help` for detailed usage information and options.

## CLI Usage Guide

### Project Initialization

`fedgwas-sim init` creates the base project structure in the current directory. It does not choose an experiment preset or generate data.

```bash
mkdir my_study
cd my_study
fedgwas-sim init
```

Initial structure:

```text
my_study/
  fedgwas.yaml
  pyproject.toml
  config.yaml
  configs/
  data/
  results/
  logs/
```

`fedgwas.yaml` stores CLI project settings such as `mode: simulation`, `config_dir`, `data_dir`, `results_dir`, and PLINK discovery settings.

`pyproject.toml` contains the minimal Flower app wiring:

```toml
[tool.flwr.app.components]
serverapp = "pipeline.server_app:app"
clientapp = "pipeline.client_app:app"

[tool.flwr.federations]
default = "local-simulation"
```

### Built-In Examples

Use `init --from-example` to initialize a project from a packaged example template that mirrors the repository experiments with project-relative paths. By default, examples also prepare their configured data so the project can move directly to `check` and `run`:

```bash
fedgwas-sim init --from-example syn-tiny
fedgwas-sim init --from-example syn-small
fedgwas-sim init --from-example syn-medium
fedgwas-sim init --from-example 1000genomes
```

Use `--no-prepare-data` when working offline, avoiding large downloads, or only inspecting the generated YAML/TOML files:

```bash
fedgwas-sim init --from-example syn-tiny --no-prepare-data
fedgwas-sim init --from-example 1000genomes --no-prepare-data
```

Synthetic examples accept `--seed` for reproducible generated data:

```bash
fedgwas-sim init --from-example syn-tiny --seed 42
```

Supported examples:

| Example | Repository reference | Purpose |
| --- | --- | --- |
| `syn-tiny` | `experiments/correctness/tiny_even` | Small correctness experiment. |
| `syn-small` | `experiments/performance/small_even` | Small performance benchmark. |
| `syn-medium` | `experiments/performance/medium_even` | Medium performance benchmark. |
| `1000genomes` | `experiments/real_world/1000genomes` | 1000 Genomes chr22 application experiment. |

Example initialization writes `config.yaml`, `configs/server.yaml`, and `configs/center_*.yaml` using the same key experiment settings as the repository
example. Paths are normalized to the local project layout. Synthetic examples generate PLINK triplets under `data/center_*`. The `1000genomes` example
downloads chromosome 22 inputs, converts them to PLINK, assigns binary phenotypes, and partitions samples into two center datasets.

### Setup Experiment

`fedgwas-sim setup-experiment <preset>` configures the current project for a named preset. If the current directory is not initialized, the command creates the standard project files as part of setup.

```bash
fedgwas-sim setup-experiment syn-tiny
fedgwas-sim setup-experiment syn-small --seed 42
fedgwas-sim setup-experiment syn-medium
fedgwas-sim setup-experiment 1000genomes
fedgwas-sim setup-experiment hapmap
```

Synthetic presets generate PLINK data automatically. After:

```bash
fedgwas-sim setup-experiment syn-tiny
```

the project has a runnable structure:

```text
my_study/
  fedgwas.yaml
  pyproject.toml
  config.yaml
  configs/
    server.yaml
    center_1.yaml
    center_2.yaml
  data/
    center_1/
      tiny_center_1.bed
      tiny_center_1.bim
      tiny_center_1.fam
    center_2/
      tiny_center_2.bed
      tiny_center_2.bim
      tiny_center_2.fam
  results/
    center_1/
      intermediate/
      logs/
    center_2/
      intermediate/
      logs/
    server/
      intermediate/
      logs/
  logs/
```

Real-world presets default to running their data preparation adapter. Use `--no-download` to generate only configuration and preparation scripts:

```bash
fedgwas-sim setup-experiment hapmap --no-download
fedgwas-sim setup-experiment 1000genomes --no-download
```

### Data Preparation

`fedgwas-sim prepare-data` is data-focused. It can regenerate missing synthetic data for synthetic presets or run the generated preparation script for
real-world presets.

```bash
fedgwas-sim prepare-data
fedgwas-sim prepare-data --download
```

For user-supplied data, center configs should point to PLINK prefixes without file extensions:

```yaml
input_data:
  path: data/center_1/study_center_1
  type: bed
```

### Checks

`fedgwas-sim check` is the unified readiness command. It checks:

- `fedgwas.yaml` exists and declares `mode: simulation`.
- Python, FedGWAS, Flower, and PLINK are available.
- Center and server configs exist.
- Configured center PLINK triplets exist.
- Results and logs directories are writable.

Run all checks:

```bash
fedgwas-sim check
```

Run only selected categories:

```bash
fedgwas-sim check --project
fedgwas-sim check --software
fedgwas-sim check --configs
fedgwas-sim check --data
fedgwas-sim check --outputs
fedgwas-sim check --configs --data
```

Use `check --data` for data verification. The data configuration command remains available for pointing a center config at a user-supplied PLINK prefix:

```bash
fedgwas-sim check --data
fedgwas-sim data configure --center 1 --bfile data/center_1/study_center_1
```

### Summaries

Use `summarize` commands to inspect prepared data or project metadata:

```bash
fedgwas-sim summarize data --path data
fedgwas-sim summarize experiment --path .
```

The output is a Rich terminal report. The data summary reports path existence, center count, file count, PLINK triplet count, human-readable size, and a
per-center table:

```text
Data Summary
Path           my_study/data
Exists         yes
Centers        2
Files          6
PLINK triplets 2
Size           78 B

Centers
Center    Files  PLINK triplets  Size
center_1      3              1   39 B
center_2      3              1   39 B
```

The experiment summary reports preset/example metadata, scenario, client count, key directories, and readiness signals:

```text
Experiment Summary
Path       my_study
Mode       simulation
State      configured
Preset     syn-tiny
Example    None
Experiment tiny_even
Category   correctness
Scenario   correctness_tiny
Clients    2

Project Readiness
Area     Detail              Status
Configs  my_study/configs    2 center config(s)
Data     my_study/data       2 PLINK triplets
Results  my_study/results    exists
```

### Run, Evaluate, And Collect Results

Run the local Flower simulation. The command first runs all readiness checks and only starts Flower when every check passes:

```bash
fedgwas-sim run --rounds 50
fedgwas-sim run --rounds 100 --no-stream
```

Generate a centralized comparison baseline when needed:

```bash
fedgwas-sim baseline generate
fedgwas-sim baseline generate --output results/baseline
```

Evaluate federated outputs against a centralized baseline. This command delegates to the same core evaluator as `python -m pipeline.evaluation.evaluate_all`:

```bash
fedgwas-sim evaluate
fedgwas-sim evaluate --baseline results/baseline --king
fedgwas-sim evaluate results --baseline results/baseline --qc-only
fedgwas-sim evaluate results --baseline results/baseline --lr-only --no-plots
fedgwas-sim evaluate results --baseline results/baseline --king-only --king-center-id 1
```

By default, `evaluate` uses the configured project `results_dir`, reads the baseline from `results/baseline`, writes `evaluation_report.md`, and runs QC +
LR evaluation. Add `--king` to include KING accumulator comparison, or use one of `--qc-only`, `--lr-only`, or `--king-only` for a single stage.

Collect run metadata and timing files:

```bash
fedgwas-sim results collect
fedgwas-sim results collect --time-file results/server_app_time.txt --label tiny_run
```

### Reset A Project

`fedgwas-sim reset` restores the current simulation project to the same base state created by `fedgwas-sim init`.

```bash
fedgwas-sim reset
fedgwas-sim reset --yes
fedgwas-sim clear --yes
```

The command only runs inside a directory with `fedgwas.yaml` and `mode: simulation`. By default it asks for confirmation; use `--yes` in scripts
or tests.

Reset removes only CLI-managed paths:

```text
config.yaml
pyproject.toml
configs/
data/
results/
logs/
scripts/
```

Unknown user files such as `README.md`, notes, notebooks, and custom scripts are left untouched. Use keep flags for large or hand-managed artifacts:

```bash
fedgwas-sim reset --keep-data --yes
fedgwas-sim reset --keep-configs --yes
fedgwas-sim reset --keep-results --yes
```
