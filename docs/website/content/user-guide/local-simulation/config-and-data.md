---
title: Project and Data Setup
slug: /local-simulation/config-and-data
---

# Project, Data and Configurations Setup

The next step after installation is to understand and create the configuration and data files for your local simulated federated GWAS screening run.

## Initialize study directory

In order to start use FedGWAS for local simulation, you need to firstly create a project or study directory, this directory in your local machine will be the place where you will run the local federated GWAS simulation or experiments. For example, you can create a directory called `my_study` and run the following command to initialize the project use `fedgwas-sim` CLI:

```bash
mkdir my_study
cd my_study
fedgwas-sim init
```

This command will initialize the project layout for you for local federated GWAS simulation as follows:

```text
my_study/
  configs/
  data/
  results/
  logs/
  fedgwas.yaml
  pyproject.toml
  config.yaml
```

- `configs/` is the directory where you will place the YAML configuration files for the server and each center.
- `data/` is the directory where you will place the PLINK binary files for each center.
- `results/` is the directory where the simulation results will be stored, with separate subdirectories for the server and each center.
- `logs/` is the directory where logs will be stored, with separate subdirectories for the server and each center.
- `fedgwas.yaml` is the project metadata file that identifies this directory as a FedGWAS simulation project and contains project-level defaults.
- `pyproject.toml` is the Flower project file that configures the simulation app and its connection to the installed FedGWAS package.
- `config.yaml` is the experiment-level configuration file that describes the overall study, including experiment name, category, number of clients, scenario, client config file paths, QC thresholds, and server parameters.

## Configuration Files

In order to run a local federated GWAS simulation, you need to create the configuration files for the project, experiment, the server and each center. The `fedgwas-sim setup-experiment` command generates these files for you with preset values based on the selected scenario and dataset size. You can also create and edit these files manually. 

For example, create a runnable synthetic project with:

```bash
mkdir my_study
cd my_study
fedgwas-sim init
fedgwas-sim setup-experiment syn-tiny --seed 42
```

The generated layout is:

```text
my_study/
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
    server/
      intermediate/
      logs/
    center_1/
      intermediate/
      logs/
    center_2/
      intermediate/
      logs/
  logs/
  fedgwas.yaml
  pyproject.toml
  config.yaml
```

`server.yaml` and `center_*.yaml` are the server and client configuration files that declare input data, output locations, parameters, thresholds, and participation settings. The `data/center_*/tiny_center_*` files are the synthetic PLINK binary datasets for each center.  `results/server/`, `results/center_1/`, and `results/center_2/` are the output directories for the server and each center, with separate subdirectories for intermediate results and logs.


### Project-level Configurations

`fedgwas.yaml` is owned by the simulation CLI. It identifies the directory as a FedGWAS simulation project and records project-level defaults:

```yaml
mode: simulation
project_state: configured
preset: syn-tiny
num_clients: 2
config_dir: configs
data_dir: data
results_dir: results
logs_dir: logs
plink: auto
```

`plink: auto` means FedGWAS searches for `plink`, then `plink2`, then bundled repository PLINK binaries when available. Set this field to a specific executable path only when automatic discovery is not appropriate.

### Flower Project File

`pyproject.toml` wires the generated study directory to the installed FedGWAS Flower app:

```toml
[tool.flwr.app.components]
serverapp = "pipeline.server_app:app"
clientapp = "pipeline.client_app:app"

[tool.flwr.app.config]
simulation = true
num-server-rounds = 100
config_path = "configs"

[tool.flwr.federations.local-simulation]
options.num-supernodes = 2
```

`fedgwas-sim run --rounds N` launches:

```bash
flwr run . local-simulation --stream --run-config \
  'simulation=true num-server-rounds=N config_path="configs"'
```

### Experiment-Level Config

`config.yaml` describes the overall study. It is used by baseline generation, evaluation, monitoring, retention, and some server defaults:

```yaml
experiment_name: tiny_even
experiment_category: correctness
num_clients: 2
scenario: correctness_tiny

clients:
  config_files:
    0: configs/center_1.yaml
    1: configs/center_2.yaml
  thresholds:
    maf_threshold: 0.01
    missing_threshold: 0.05
    hwe_threshold: 1e-6
    p_threshold: 0.3
    king_threshold: 0.23

server:
  chunk_size: 200
  min_available_clients: 1
  min_fit_clients: 1
  num_server_rounds: 50
```

The generated values are safe defaults for the selected preset. Tune them only when you understand the effect on QC thresholds, chunking, runtime, and evaluation comparability.

### Client and Server Configs

#### Server Config

`configs/server.yaml` declares server output locations:

```yaml
output:
  intermediate_dir: results/server/intermediate
  log_dir: results/server/logs
```

The server app requires this output section. During startup it derives the server result directory from either `output.log_dir` or
`output.intermediate_dir`.

#### Clients Center Configs

Each center has its own config file in `configs/center_*.yaml` that declares input data, output locations, parameters, thresholds, and participation settings. An example of a generated center config looks like:

```yaml
input_data:
  path: data/center_1/tiny_center_1
  type: bed

output:
  intermediate_dir: results/center_1/intermediate
  log_dir: results/center_1/logs

parameters:
  sample_offset: 1000000000000
  chunk_size: 200
  sample_chunk_size: 100
  snp_chunk_size: 200
  lr_target_chunks: 2
  run_id: "20260615_120000"

thresholds:
  maf_threshold: 0.01
  missing_threshold: 0.05
  hwe_threshold: 1e-6
  p_threshold: 0.3
  local_lr_threshold: 0.3
  global_lr_threshold: 0.1
  king_threshold: 0.23

participation:
  key_exchange: true
  sync: true
  local_qc: true
  global_qc: true
  global_qc_response: true
  init_chunks: true
  iterative_king: true
  local_lr: true
  local_lr_filter_response: true
  init_chunks_lr: true
  iterative_lr: true
```

`participation` controls whether a client joins each pipeline stage. For normal end-to-end runs, leave all stages enabled.

**PLINK Prefix Convention:**

For binary PLINK input, `input_data.path` is a prefix, not a file name with an extension. This is correct:

```yaml
input_data:
  path: data/center_1/study_center_1
  type: bed
```

The corresponding files must exist:

```text
data/center_1/study_center_1.bed
data/center_1/study_center_1.bim
data/center_1/study_center_1.fam
```

Do not write `.bed`, `.bim`, or `.fam` in `input_data.path`.

## Data

This section describes the expected data format and location for local simulation, and how to point a generated project at your own PLINK data.

### Data Format And Location

FedGWAS accepts genotype data in PLINK 1.9 binary format (`.bed`, `.bim`, `.fam`). VCF input is supported via conversion to PLINK-compatible representations where configured. You can put the PLINK binary files for each client in the `data/` directory. For example, you can have `data/center_1/` and `data/center_2/` directories, each containing the PLINK binary files for that client. 

```text
my_study/
  data/
    center_1/
      study_center_1.bed
      study_center_1.bim
      study_center_1.fam
    center_2/
      study_center_2.bed
      study_center_2.bim
      study_center_2.fam
```

Inside `data/center_*/`, the PLINK binary files for that center should be named with a common prefix (e.g., `study_center_1`) and the appropriate extensions (`.bed`, `.bim`, `.fam`). The `input_data.path` field in each center's config file should point to the prefix of the PLINK files for that center. 

### Synthetic Data Simulation

We provide some functions for simulating PLINK binary datasets for testing and demonstration purposes. These functions can be used to generate synthetic data with specified numbers of samples, SNPs, cases, controls, and random seeds. The generated PLINK files can then be placed in the appropriate `data/center_*/` directories and pointed to by the center config files.

To point a generated project at your own PLINK data:

```bash
fedgwas-sim data configure --center 1 --bfile data/center_1/study_center_1
fedgwas-sim data configure --center 2 --bfile data/center_2/study_center_2
fedgwas-sim check --data
```

### Real-world Data Instructions

We have include some instructions and scripts for preparing some of the real-world datasets used in our paper for local simulation runs. See the `experiments/real-world/` directory in the repository for details.

### Data Summary

We also provide a command to summarize the PLINK data for each center, which can be useful for verifying the data and understanding its characteristics before running the simulation:

```bash
fedgwas-sim summarize data
```

This will output a summary of data as follows:
```
╭────────────────── Data Summary ─────────────────╮
│ Path            XX\my_study\data                │
│ Exists          yes                             │
│ Centers         2                               │
│ Files           20                              │
│ PLINK triplets  3                               │
│ Size            1.6 MiB                         │
╰─────────────────────────────────────────────────╯
                  Centers                     
┏━━━━━━━━━━┳━━━━━━━┳━━━━━━━━━━━━━━━━┳━━━━━━━━━━━┓
┃ Center   ┃ Files ┃ PLINK triplets ┃      Size ┃
┡━━━━━━━━━━╇━━━━━━━╇━━━━━━━━━━━━━━━━╇━━━━━━━━━━━┩
│ center_1 │     7 │              1 │ 449.2 KiB │
│ center_2 │     7 │              1 │ 449.0 KiB │
└──────────┴───────┴────────────────┴───────────┘
```

## Use Setup From Example Experiments

We provide commands to set up the configuration and data for some example experiments that you can run locally. These commands generate the config files and synthetic PLINK data for you based on preset values for different scenarios and dataset sizes. You can use these generated files as templates for setting up your own experiments with custom data and configurations. 

Use the `fedgwas-sim init --from-example <preset>` command to initialize a new project with the example experiment files. We currently support the following presets:
- `syn-tiny`: A tiny synthetic dataset with 100 samples, 100 SNPs, and balanced cases and controls. Good for quick testing and debugging.
- `syn-small`: A small synthetic dataset with 1000 samples, 1000 SNPs, and balanced cases and controls. Good for more extensive testing and debugging.
- `syn-medium`: A medium synthetic dataset with 10,000 samples, 10,000 SNPs, and balanced cases and controls. Good for testing scalability and performance.
- `1000genomes`: A subset of the 1000 Genomes dataset with real-world characteristics. Good for testing with realistic data.

Example usage:

```bash
mkdir my_study
cd my_study
fedgwas-sim init --from-example syn-tiny --seed 42
```

Another way is to firstly setup the project directory layout and use the `fedgwas-sim setup-experiment` command to generate the config and data files for the example experiments:

```bash
mkdir my_study
cd my_study
fedgwas-sim init
fedgwas-sim setup-experiment syn-tiny --seed 42
```

## Reset and Clear Project

We provide a command for you to help you reset a study and keep only the initialized project structure:

```bash
fedgwas-sim reset --yes
```

Preserve expensive data or results when needed:

```bash
fedgwas-sim reset --keep-data --yes
fedgwas-sim reset --keep-results --yes
```

Unknown user files in the project directory are left untouched. You can use it to clean up the generated config and data files and start fresh without having to manually delete files or create a new directory.

## [Optional] Experiments In Repository

Repository examples under `experiments/` use an older nested config layout:

```text
experiments/correctness/tiny_even/
  config.yaml
  configs/
    server/config.yaml
    center_1/config.yaml
    center_2/config.yaml
```

The server and client apps support both this nested layout and the flat `fedgwas-sim` layout. Do not copy paths between the two styles without updating the surrounding config files.

You can refer to these examples for how to set up your own project with custom data and configurations, and run local simulations with the installed package or directly with the repository scripts.
