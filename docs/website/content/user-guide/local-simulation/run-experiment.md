---
title: Run Local FedGWAS Simulation
slug: /local-simulation/run-experiment
---

# Run Local FedGWAS Simulation

After installation and setup local project (configurations and data correctly), you can start to run a local federated GWAS simulation, it will use Flower to execute federated workflow on a single machine with multiple processes simulating the server and centers based on your project configuration and data.


## Check the Run Readiness

Before starting a local simulation run, we recommend you to run the following command to check the run readiness:

```bash
fedgwas-sim check
```

It will check the project configuration, data, and software dependencies to make sure everything is ready for a local simulation run. If any check fails, fix the issue before starting the run. The example output of a successful check is shown below:

```
╭─────────────────────── Check Summary ───────────────────────╮
│ Project  xx\my_study                                        │
│ Scope    project, software, configs, data, outputs          │
│ Result   PASS                                               │
│ Checks   11 passed, 0 failed                                │
╰─────────────────────────────────────────────────────────────╯
┏━━━━━━━━━━┳━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Aspect   ┃ Status ┃ Check                     ┃ Detail                       ┃
┡━━━━━━━━━━╇━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ Project  │ PASS   │ fedgwas.yaml              │ fedgwas.yaml                 │
│ Software │ PASS   │ Python >=3.11             │ 3.11.13                      │
│          │ PASS   │ FedGWAS importable        │ pipeline                     │
│          │ PASS   │ Flower >=1.19,<1.20       │ 1.19.0                       │
│          │ PASS   │ PLINK available           │ xx\plink\plink_win\plink.exe │
│ Configs  │ PASS   │ center_1 config           │ configs/center_1.yaml        │
│          │ PASS   │ center_2 config           │ configs/center_2.yaml        │
│          │ PASS   │ server config             │ configs/server.yaml          │
│ Data     │ PASS   │ configured PLINK triplets │ all present                  │
│ Outputs  │ PASS   │ results_dir writable      │ results                      │
│          │ PASS   │ logs_dir writable         │ logs                         │
└──────────┴────────┴───────────────────────────┴──────────────────────────────┘
```

You can also run individual checks while debugging:

```bash
fedgwas-sim check --software
fedgwas-sim check --configs --data
fedgwas-sim check --outputs
```

We also recomend you to summarize the data and experiment before running the simulation, you can run the following command to get a summary of the data and experiment:

```bash
fedgwas-sim summarize experiment
```

```
╭────────────────────── Experiment Summary ──────────────────────╮
│ Path        XX\my_study                                        │ 
│ Mode        simulation                                         │
│ State       configured                                         │
│ Preset      syn-tiny                                           │
│ Example     tiny-even                                          │
│ Experiment  tiny_even                                          │
│ Category    correctness                                        │
│ Scenario    correctness_tiny                                   │
│ Clients     2                                                  │
╰────────────────────────────────────────────────────────────────╯
```

## Run the Local Federated Simulation

Run local federated GWAS simulation is very simple, you can start Flower local simulation with just one command:

```bash
fedgwas-sim run --rounds 50
```

By default the command streams Flower output. To run without streaming:

```bash
fedgwas-sim run --rounds 50 --no-stream
```

Internally, the CLI launches Flower with the generated project config:

```bash
flwr run . local-simulation --stream --run-config \
  'simulation=true num-server-rounds=50 config_path="configs"'
```

## What Success Looks Like

A successful run:

- Passes all pre-flight checks.
- Starts a Flower local simulation with two simulated clients.
- Writes server logs under `results/server/logs`.
- Writes client logs under `results/center_*/logs`.
- Progresses through privacy-preserving initialization, federated quality control, KING-based relatedness screening, association screening, and completion.

Inspect logs:

```bash
ls results/server/logs
ls results/center_1/logs
ls results/center_2/logs
```

The exact file names can vary by stage, run id, and monitoring settings.

## [Optional] Run Federated Experiments from Repository Fallback Workflow

Use this path only from a repository checkout, usually for development or
reproducing checked-in experiments:

```bash
cd Fed-GWAS
python -m pip install -e .
```

Generate tiny synthetic data for the repository experiment:

```bash
python pipeline/simulation/simulated_data/generate_synthetic_data.py \
  --scale tiny \
  --partition-strategy even \
  --seed 42 \
  --output-dir experiments/correctness/tiny_even/data
```

Run the checked-in config:

```bash
flwr run . local-simulation --stream --run-config \
  'simulation=true num-server-rounds=100 config_path="experiments/correctness/tiny_even/configs"'
```

The repository configs commonly write under
`experiments/correctness/tiny_even/results_2/`. Check the `output` section in each center config before looking for results.
