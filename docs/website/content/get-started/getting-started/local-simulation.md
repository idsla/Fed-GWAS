---
title: Local Federated GWAS Simulation
slug: /local-simulation
---

# Local Federated GWAS Simulation

Use local simulation when you want to try FedGWAS on one machine. This mode
starts multiple virtual clients through Flower `local-simulation`, so it is the
fastest path for installation checks, smoke tests, and example runs.

## Install

Create and activate a Python environment, then install FedGWAS:

```bash
python -m pip install FedGWAS
```

Confirm that the simulation CLI is available:

```bash
fedgwas-sim --help
```

For full GWAS runs, install PLINK and make sure the executable is available on
`PATH`:

```bash
plink --help
plink2 --help
```

## Create a study

Start from an empty directory owned by the user:

```bash
mkdir my_study
cd my_study

fedgwas-sim init
fedgwas-sim setup-experiment syn-tiny --seed 42
```

This writes a local Flower project, center configs, generated synthetic PLINK
data, and output directories under `my_study/`.

## Validate and run

Run readiness checks before starting Flower:

```bash
fedgwas-sim check
```

Start the local simulation:

```bash
fedgwas-sim run --rounds 100
```

## Evaluate outputs

Generate a centralized comparison baseline, evaluate the federated outputs, and
collect run metadata:

```bash
fedgwas-sim baseline generate --output data/centralized_baseline
fedgwas-sim evaluate results --baseline data/centralized_baseline --king
fedgwas-sim results collect --label tiny_run
```

## Next steps

- Use [fedgwas-sim CLI](/user-guide/cli-simulation) for the complete command reference.
- Use [Simulation Mode](/user-guide/simulation) for repository-level Flower settings and
  runtime behavior.
- Use [Examples](/examples/overview) for correctness, performance, and
  1000 Genomes scenarios.
