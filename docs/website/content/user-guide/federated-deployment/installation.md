---
title: Installation and Overview
slug: /federated-deployment/installation
---

# Installation and Overview

## Overview

Federated deployment runs FedGWAS as separate Flower processes across machines
or containers. One coordinating **server** starts SuperLink and submits the app
with `flwr run`. Each participating **client** (institution / data holder) starts
a SuperNode with its own center config and local PLINK data. Genotype-level work
stays on the client that holds the data.

Use this mode when simulation on one machine is not enough: the server and
clients are real processes, often on different hosts.

### Two deployment modes

FedGWAS provides two equivalent ways to start that runtime. Both wrap
`flower-superlink`, `flower-supernode`, and `flwr run . local-deployment`. Pick
one mode for a given run; do not mix them.

| Mode | Entry point | When to use |
| --- | --- | --- |
| [CLI](cli-deployment.md) | `fedgwas-deploy` | Package commands on each node |
| [Script deployment](script-deployment.md) | `cluster_deployment/scripts/` | Shell wrappers from a repository checkout |

Shipped examples often use two clients. The layout itself is one server plus
**N** clients: add one SuperNode, one center config, and one local dataset per
participating institution.

### Roles

```text
                    ┌─────────────────────────────┐
                    │  Server                     │
                    │  SuperLink  :9092 / :9093   │
                    │  flwr run (Exec API)        │
                    └──────────┬──────────┬───────┘
                               │          │
              SuperNode        │          │        SuperNode
                               ▼          ▼
                    ┌──────────────┐  ┌──────────────┐
                    │ Client 1     │  │ Client N     │
                    │ ClientAppIo  │  │ ClientAppIo  │
                    │ local PLINK  │  │ local PLINK  │
                    └──────────────┘  └──────────────┘
```

| Role | Runs | Holds |
| --- | --- | --- |
| Server | SuperLink and `flwr run` | App code, server config, server logs, aggregation outputs |
| Client *k* | One SuperNode | That client's center config and that client's PLINK files |

The server relays protocol messages and aggregates selected outputs. It does not
need raw client genotypes to run the workflow. Each client should keep only its
own individual-level data.

All nodes must use the same FedGWAS version and Flower 1.19.x. After
installation, start the federation with [CLI Deployment](cli-deployment.md) or
[Script Deployment](script-deployment.md).

## Installation

Complete these checks on every node before starting SuperLink or SuperNodes.

### Requirements

| Dependency | Where | Requirement |
| --- | --- | --- |
| Python | Server and every client | 3.11 or later |
| FedGWAS | Server and every client | Same package version or same repository checkout |
| Flower | Server and every client | 1.19.x |
| PLINK 1.9 (`plink`) | Every client; optional on server | On `PATH` or set in the center config |
| PLINK 2 (`plink2`) | Every client that runs relatedness screening | On `PATH` or set in the center config |
| Shell utilities | Optional, for script deployment | `bash`, `time`, `ss` or `netstat` |

Install from PyPI:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install FedGWAS
python -m pip install "flwr[simulation]>=1.19.0,<1.20"
```

Or from a repository checkout (needed for script deployment):

```bash
git clone https://github.com/idsla/Fed-GWAS.git
cd Fed-GWAS
python -m pip install -e .
python -m pip install "flwr[simulation]>=1.19.0,<1.20"
```

### Data

Each client stores only its own PLINK 1.9 binary files and center config:

```text
client-k/
  configs/center_k.yaml
  data/center_k/study_center_k.{bed,bim,fam}
```

The server needs the FedGWAS app and a server config. Copy client result
folders to the server later only if you want a centralized baseline or a
merged evaluation tree.

### Config

In deployment mode, each SuperNode gets an explicit center config through Flower
node config. `partition-id` is zero-based; `num-partitions` is the total number
of clients:

```bash
--node-config 'partition-id=<k-1> num-partitions=<N> config-file="configs/center_k.yaml"'
```

The client app uses this `config-file` instead of simulation `config_path`
lookup, so each institution can keep local data and output paths.

Match the client count in the Flower federation as well, for example
`options.num-supernodes = N` under the deployment federation in
`pyproject.toml`. The CLI `--num-centers` flag and the script `--client-id`
values must describe the same N clients.

### Ports

The default deployment uses insecure local-network Flower ports:

| Port | Direction | Purpose |
| --- | --- | --- |
| `9092` | Every client → server | SuperLink Fleet API; SuperNodes connect here |
| `9093` | Operator or server local → server | SuperLink Exec API; `flwr run` submits the app here |
| `9093 + k` | Server → client *k* | ClientAppIo API for that SuperNode (`9094`, `9095`, …) |

Open `9092` from every client to the server, `9093` on the server for app
submission, and one ClientAppIo port from the server to each client. In cloud
deployments, update security groups, network ACLs, and host firewalls.

Connectivity checks:

```bash
# From each client
ping <server-ip>
bash -c 'echo >/dev/tcp/<server-ip>/9092' && echo OK

# From the server, once per client (port 9093+k)
bash -c 'echo >/dev/tcp/<client-k-ip>/<client-k-port>' && echo OK
```

### Clock synchronization

Flower 1.19 message authentication is sensitive to clock skew. Keep server and
client UTC clocks closely aligned before starting SuperNodes.

On Linux nodes:

```bash
timedatectl set-ntp true 2>/dev/null || true
date -u
```

From a repository checkout, check skew from each client:

```bash
bash cluster_deployment/scripts/check-cluster-clock-skew.sh <server-ip>
```

Fix clock skew before debugging `Invalid timestamp` or `UNAUTHENTICATED`.

### Checks

On the server and every client:

```bash
python --version
python -c "import pipeline; import flwr; print('FedGWAS and Flower import OK', flwr.__version__)"
flwr --version
flower-superlink --help
flower-supernode --help
```

On every client:

```bash
plink --version
plink2 --version
```

If `flwr --version` is not 1.19.x:

```bash
python -m pip install --upgrade "flwr[simulation]>=1.19.0,<1.20"
```

With the CLI, you can also run:

```bash
fedgwas-deploy server check
fedgwas-deploy client check --server <server-ip> --center-id <k> --config configs/center_k.yaml --network
```

## Next Steps

After the checks above pass, start the federation with one of the two modes:

1. [CLI Deployment](cli-deployment.md) — `fedgwas-deploy` on the server and on
   each client.
2. [Script Deployment](script-deployment.md) — `cluster_deployment/scripts/`
   from a repository checkout.

A worked three-node example is in
[Three-Node Deployment](/examples/three-node-deployment).
