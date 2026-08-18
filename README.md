<p align="center">
  <img src="https://raw.githubusercontent.com/idsla/Fed-GWAS/main/docs/website/images/logo-readme.png" alt="FedGWAS logo" width="640">
</p>

# FedGWAS: a lightweight federated pipeline for privacy-preserving GWAS screening

[![PyPI](https://img.shields.io/pypi/v/FedGWAS.svg)](https://pypi.org/project/FedGWAS/)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-2f8f83)](https://idsla.github.io/Fed-GWAS/)
[![Deploy documentation](https://github.com/idsla/Fed-GWAS/actions/workflows/deploy-docs.yml/badge.svg)](https://github.com/idsla/Fed-GWAS/actions/workflows/deploy-docs.yml)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Overview

FedGWAS is a lightweight federated pipeline for privacy-preserving GWAS screening across distributed genomic datasets. It coordinates federated quality control, KING-based relatedness screening, and logistic-regression association screening under a client–server architecture, while keeping genotype-level computations local to participating clients. FedGWAS uses Flower for federated coordination, PLINK for genetics operations, and encryption, shuffling, anonymization, and lightweight secret-sharing so that the server can relay selected protocol messages without decrypting them.

Use FedGWAS when multiple participating institutions need a shared GWAS screening lifecycle and cannot pool individual-level genotype data.

<p align="center">
  <img src="https://raw.githubusercontent.com/idsla/Fed-GWAS/main/docs/website/images/current_overview.png" alt="FedGWAS workflow illustration" width="100%">
</p>

To get started and learn how to use FedGWAS, use the following resources:

- Documentation site: [Documentation](https://idsla.github.io/Fed-GWAS/)
- Examples gallery: [Examples](https://idsla.github.io/Fed-GWAS/examples/overview)
- API reference: [API Reference](https://idsla.github.io/Fed-GWAS/api-reference/api)
- Technical details: [Technical Details](https://idsla.github.io/Fed-GWAS/user-guide/design/workflow)

## Prerequisites

FedGWAS requires Python 3.11 or later, Flower, and PLINK. Input genotypes use [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/) binary format (`plink` on `PATH` or configured per client). KING-based relatedness screening also needs [PLINK 2](https://www.cog-genomics.org/plink/2.0/) available as `plink2`. VCF input is supported via conversion to PLINK-compatible representations where configured.

```bash
plink --version
plink2 --version
```

For repository-based runs, you can also set the PLINK path in each client config if your environment does not expose `plink` globally.

## FedGWAS Local Simulation Guide

Local simulation mode runs the full FedGWAS screening workflow on one machine by launching multiple simulated clients and a federated server through Flower. Use it to validate an installation, create simulation experiments from preset settings and generated data, prototype client (center) configs, and compare federated outputs against a centralized baseline without setting up a real federated deployment.

You can start local simulation in either of two ways:

1. Recommended: install from PyPI and use `fedgwas-sim` command line interface (CLI)
2. Repository/local workflow: clone this repository and run the old scripts directly

Both workflows require:

- Python 3.11 or later
- PLINK 1.9 (`plink`) and PLINK 2 (`plink2`) available on `PATH` or configured locally
- Flower installed through the package or local environment

### Recommended: PyPI CLI Workflow

Install the package:

```bash
python -m pip install FedGWAS
```

Verify that the simulation CLI is available:

```bash
fedgwas-sim --help
```

Create a standalone study directory and run the tiny two-client simulation:

```bash
mkdir my_study
cd my_study

# initialize study project directory
fedgwas-sim init
# setup data and configurations
fedgwas-sim setup-experiment syn-tiny --seed 42
# validation and run simulation
fedgwas-sim check
fedgwas-sim run --rounds 100
# evaluation and results collection
fedgwas-sim baseline generate --output data/centralized_baseline
fedgwas-sim evaluate results --baseline data/centralized_baseline --king
fedgwas-sim results collect --label tiny_run
```

The usage of the CLI can be found in the [documentation site](https://idsla.github.io/Fed-GWAS/user-guide/cli-simulation).

### Repository/Local Script Workflow

Clone the repository if you want the old direct script workflow, bundled experiment files, cluster deployment scripts, documentation source, or developer tooling:

```bash
git clone https://github.com/idsla/Fed-GWAS.git
cd Fed-GWAS
python -m pip install -e .
```

With `uv`, you can install the local environment with:

```bash
git clone https://github.com/idsla/Fed-GWAS.git
cd Fed-GWAS
uv sync --python 3.11
```

Generate synthetic data:

```bash
python pipeline/simulation/simulated_data/generate_synthetic_data.py \
  --scale tiny \
  --partition-strategy even \
  --seed 42 \
  --output-dir experiments/correctness/tiny_even/data
```

Generate the centralized baseline:

```bash
python experiments/tools/generate_baseline.py \
  experiments/correctness/tiny_even/config.yaml
```

Run the federated simulation:

```bash
flwr run . local-simulation --stream
```

Or run with explicit release-smoke settings:

```bash
flwr run . local-simulation --stream --run-config \
  'simulation=true num-server-rounds=100 config_path="experiments/correctness/tiny_even/configs"'
```

Evaluate the run:

```bash
python experiments/tools/evaluation/evaluate_all.py \
  experiments/correctness/tiny_even/results_2 \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline \
  --king
```

If you changed the active config output paths, pass the results directory from those config files instead.

### Example Simulation Experiments

We have a few preset experiments with generated data and configs in the repository for testing and demonstration. You can find the details of these example experiments in the documentation and the experiment directories:

- Tiny experiment details: [experiments/correctness/tiny_even/README.md](experiments/correctness/tiny_even/README.md)

- Small experiment details: [experiments/performance/small_even/README.md](experiments/performance/small_even/README.md)

- Real-world data experiment details: [experiments/real_world/1000genomes/README.md](experiments/real_world/1000genomes/README.md)

## FedGWAS Cluster Deployment Guide

Use federated deployment when the server and clients should run as separate Flower processes, usually on separate machines. One coordinating **server** runs SuperLink and submits the app. Each participating **client** runs a SuperNode with its own center config and local PLINK data.

FedGWAS provides two equivalent deployment modes. Both wrap `flower-superlink`, `flower-supernode`, and `flwr run . local-deployment`. Use Flower 1.19.x on every node. Pick one mode for a given run; do not mix them.

### CLI

```bash
fedgwas-deploy server start --host 0.0.0.0 --daemon --log-file /tmp/superlink.log

fedgwas-deploy client start \
  --server <server-ip> \
  --center-id <k> \
  --config configs/center_k.yaml \
  --daemon

fedgwas-deploy server run --server <server-ip> --rounds 20 --scale tiny
```

Stop on each node separately: `fedgwas-deploy server stop` on the server, `fedgwas-deploy client stop` on each client.

### Script deployment

Clone this repository on each node and run the wrappers from the repository root:

```bash
nohup cluster_deployment/scripts/cluster-start-server.sh --server-ip 0.0.0.0 \
  > /tmp/superlink.log 2>&1 &

cluster_deployment/scripts/cluster-start-client.sh \
  --server-ip <server-ip> \
  --client-id <k> \
  --config configs/center_k.yaml

cluster_deployment/scripts/cluster-run-app.sh \
  --server-ip <server-ip> \
  --rounds 20 \
  --scale tiny
```

Stop on each node separately: `cluster_deployment/scripts/cluster-stop-all.sh`.

Guides:

- Get Started: [Federated Deployment](https://idsla.github.io/Fed-GWAS/get-started/federated-deployment)
- [CLI Deployment](https://idsla.github.io/Fed-GWAS/user-guide/federated-deployment/cli-deployment)
- [Script Deployment](https://idsla.github.io/Fed-GWAS/user-guide/federated-deployment/script-deployment)
- Worked example: [examples/04three-node-deployment](examples/04three-node-deployment/README.md)

## Federated Protocol Summary

FedGWAS implements a four-stage GWAS screening lifecycle. These conceptual stages are realized as coordinated federated rounds (including initialization, chunking, and aggregation):

1. Privacy-preserving initialization: clients exchange encrypted seed shares via the server, which acts only as a message relay, and derive a shared global seed
2. Federated quality control: local sample-level missingness filtering, then encrypted variant-level summaries for harmonized MAF, missingness, and Hardy–Weinberg filters
3. KING-based relatedness screening: chunked, anonymized kinship estimation
4. Federated association screening: privacy-preserving tokens for local logistic-regression filtering, then federated logistic regression

The current implementation uses encryption, shuffling, anonymization, and lightweight secret-sharing. The server relays encrypted client-to-client payloads for selected stages and does not decrypt those payloads. Logging, metrics, and production- versus research-oriented output retention are configured separately and are not screening stages. See [CURRENT_VERSION.md](docs/CURRENT_VERSION.md) for the current privacy model, stage contracts, and limitations.

Population stratification handling, continuous traits, and association models other than case–control logistic regression are not part of the current screening protocol.

## Troubleshooting (Common Issues)

- `plink` or `plink2` not found: install PLINK 1.9 and PLINK 2, and make sure both are on `PATH` or set in the client config.
- Flower uses the wrong config: pass `--run-config 'config_path="..."'`.
- Empty or missing results: generate the tiny synthetic data and baseline before running.
- TestPyPI or PyPI install fails for a new release: check that the version in `pyproject.toml` has been published and that dependency resolution can reach normal PyPI.

## License

FedGWAS is distributed under the MIT License. See [LICENSE](LICENSE).

## Contributors and Creator

Software development is hosted at the [Rutgers Institute in Data Science, Learning, and Application](https://sites.rutgers.edu/idsla/). The accompanying paper authors are Xinyue Wang (Renmin University of China), Sitao Min, and Jaideep Vaidya (Rutgers University).

**Contributors**:

- Dr. Xinyue Wang <a href="mailto:xinyue.wang@ruc.edu.cn" aria-label="Email Dr. Xinyue Wang">&#9993;</a>
- Dr. Sitao Min <a href="mailto:sitaomin1994@gmail.com" aria-label="Email Dr. Sitao Min">&#9993;</a>
- Dr. Jaideep Vaidya <a href="mailto:jsvaidya@business.rutgers.edu" aria-label="Email Dr. Jaideep Vaidya">&#9993;</a>
