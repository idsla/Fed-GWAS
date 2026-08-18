---
title: Local Simulation Experiments
slug: /local-simulation/experiments
---

# Local Simulation Experiments

FedGWAS provides two ways to create local simulation studies:

- Presets for installed-package users through `fedgwas-sim setup-experiment`.
- Built-in examples that mirror repository experiment layouts through
  `fedgwas-sim init --from-example`.

Both approaches create a self-contained study directory with project-relative
paths.

## Choose A Workflow

| Workflow | Best for | Command |
| --- | --- | --- |
| Preset | New local simulation project from packaged defaults | `fedgwas-sim setup-experiment <preset>` |
| Built-in example | Reproducing a documented repository scenario with local paths | `fedgwas-sim init --from-example <example>` |
| Repository experiment | Developing or reproducing checked-in `experiments/...` folders | `flwr run . local-simulation` from the repository root |

Use the preset workflow unless you specifically need a repository scenario.

## Presets

Run a preset from an initialized or empty study directory:

```bash
mkdir my_study
cd my_study
fedgwas-sim setup-experiment syn-tiny --seed 42
```

Supported presets:

| Preset | Data source | Clients | Main use case |
| --- | --- | --- | --- |
| `syn-tiny` | Synthetic `tiny` scale | 2 | Correctness checks, smoke tests, quick demonstrations |
| `syn-small` | Synthetic `small` scale | 2 | Small performance and runtime checks |
| `syn-medium` | Synthetic `medium` scale | 2 | Larger performance experiments on a capable workstation |
| `1000genomes` | Public-data template | 2 | 1000 Genomes chromosome 22 simulation workflow |
| `hapmap` | Public-data template | 2 | HapMap-style real-data template |

Synthetic presets generate PLINK `.bed/.bim/.fam` triplets immediately. Real
data presets write configuration and preparation scripts, and may require
dataset-specific access, filtering, phenotype, and storage decisions.

Use `--no-download` for real-data presets when you want to inspect the generated
scripts first:

```bash
fedgwas-sim setup-experiment 1000genomes --no-download
fedgwas-sim setup-experiment hapmap --no-download
```

## Built-In Examples

Built-in examples mirror the repository examples while normalizing paths for a
standalone study directory:

```bash
mkdir tiny_example
cd tiny_example
fedgwas-sim init --from-example syn-tiny --seed 42
```

Supported examples:

| Example | Repository reference | Purpose |
| --- | --- | --- |
| `syn-tiny` | `experiments/correctness/tiny_even` | Two-center correctness validation |
| `syn-small` | `experiments/performance/small_even` | Small scalability run |
| `syn-medium` | `experiments/performance/medium_even` | Medium scalability run |
| `1000genomes` | `experiments/real_world/1000genomes` | 1000 Genomes chromosome 22 application experiment |

Skip data preparation when working offline or when preparing a large public
dataset manually:

```bash
fedgwas-sim init --from-example 1000genomes --no-prepare-data
```

## Correctness Experiments

Use correctness experiments to verify that the complete federated workflow can
run end to end and that outputs agree with centralized baselines within expected
tolerances:

```bash
fedgwas-sim setup-experiment syn-tiny --seed 42
fedgwas-sim check
fedgwas-sim run --rounds 50
fedgwas-sim baseline generate
fedgwas-sim evaluate --king
```

The tiny scale is the best first target because it is fast and exercises federated QC,
relatedness screening, association screening, evaluation, and result collection.

## Performance Experiments

Use `syn-small` and `syn-medium` to evaluate runtime and resource behavior:

```bash
fedgwas-sim setup-experiment syn-small --seed 42
fedgwas-sim check
fedgwas-sim run --rounds 80
fedgwas-sim results collect --label small_run
```

Performance presets are configured primarily for runtime and monitoring. They
may disable baseline comparison by default because centralized comparison can be
more expensive than the smoke-test workflow.

## Real-World Templates

Real-world templates are starting points. Before running them, confirm:

- Dataset access and license requirements.
- Storage and runtime budget.
- Phenotype generation or phenotype file policy.
- Center partitioning strategy.
- QC thresholds and covariates appropriate for the study.

For 1000 Genomes chromosome 22:

```bash
fedgwas-sim setup-experiment 1000genomes --no-download
python scripts/prepare_data.py --download
fedgwas-sim check
```

## Repository Experiments

Repository experiments remain useful for development, regression testing, and
paper/release reproduction. From the repository root:

```bash
python pipeline/simulation/simulated_data/generate_synthetic_data.py \
  --scale tiny \
  --partition-strategy even \
  --seed 42 \
  --output-dir experiments/correctness/tiny_even/data

flwr run . local-simulation --stream --run-config \
  'simulation=true num-server-rounds=100 config_path="experiments/correctness/tiny_even/configs"'
```

The repository layout uses paths under `experiments/...`; the installed CLI
layout uses paths under a user-owned study directory. Keep these workflows
separate when writing commands or interpreting outputs.
