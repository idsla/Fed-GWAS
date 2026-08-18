---
title: Federated Deployment
slug: /federated-deployment
---

# Federated Deployment

Use federated deployment when the server and clients should run as separate
Flower processes, usually on separate machines or containers. This is closer to
a real multi-site collaboration than [local simulation](/get-started/local-simulation).

One coordinating **server** runs SuperLink and submits the app. Each
participating **client** (institution / data holder) runs a SuperNode with its
own center config and local PLINK data. Genotype-level work stays on the client
that holds the data.

FedGWAS provides two equivalent ways to start that runtime. Both wrap
`flower-superlink`, `flower-supernode`, and `flwr run . local-deployment`. Pick
one mode for a given run; do not mix them.

| Mode | Entry point | When to use |
| --- | --- | --- |
| CLI | `fedgwas-deploy` | Package commands on each node |
| Script deployment | `cluster_deployment/scripts/` | Shell wrappers from a repository checkout |

Full role, port, and check details are in
[Installation and Overview](/user-guide/federated-deployment/installation).

## Prepare each node

Install the same FedGWAS version and Flower 1.19.x on the server and every
client:

```bash
python -m pip install FedGWAS
python -m pip install "flwr[simulation]>=1.19.0,<1.20"
```

Script deployment also needs a repository checkout on each node:

```bash
git clone https://github.com/idsla/Fed-GWAS.git
cd Fed-GWAS
python -m pip install -e .
```

Each client needs PLINK and only its own genotype files:

```bash
plink --version
plink2 --version
```

```text
client-k/
  configs/center_k.yaml
  data/center_k/study_center_k.{bed,bim,fam}
```

Confirm the chosen entry point:

```bash
fedgwas-deploy --help
ls cluster_deployment/scripts/
```

## CLI

On the server:

```bash
fedgwas-deploy server start --host 0.0.0.0 --daemon --log-file /tmp/superlink.log
```

On each client, after SuperLink is up (`<k>` is `1`, `2`, …):

```bash
fedgwas-deploy client start \
  --server <server-ip> \
  --center-id <k> \
  --config configs/center_k.yaml \
  --daemon \
  --log-file /tmp/supernode_k.log
```

On the server, after every SuperNode has registered:

```bash
fedgwas-deploy server run --server <server-ip> --rounds 20 --scale tiny
```

Stop on each node separately:

```bash
fedgwas-deploy server stop    # server node
fedgwas-deploy client stop    # each client node
```

See [CLI Deployment](/user-guide/federated-deployment/cli-deployment)
for checks, `--num-centers`, ports, and troubleshooting.

## Script deployment

Run these from the repository root. On the server:

```bash
nohup cluster_deployment/scripts/cluster-start-server.sh --server-ip 0.0.0.0 \
  > /tmp/superlink.log 2>&1 &
```

On each client:

```bash
cluster_deployment/scripts/cluster-start-client.sh \
  --server-ip <server-ip> \
  --client-id <k> \
  --config configs/center_k.yaml
```

On the server, after every SuperNode has registered:

```bash
cluster_deployment/scripts/cluster-run-app.sh \
  --server-ip <server-ip> \
  --rounds 20 \
  --scale tiny
```

Stop on each node separately:

```bash
cluster_deployment/scripts/cluster-stop-all.sh
```

See [Script Deployment](/user-guide/federated-deployment/script-deployment)
for data checks, clock skew, and troubleshooting.

## Next steps

- [Installation and Overview](/user-guide/federated-deployment/installation):
  roles, requirements, ports, clock sync, and checks.
- [CLI Deployment](/user-guide/federated-deployment/cli-deployment)
- [Script Deployment](/user-guide/federated-deployment/script-deployment)
- [Three-Node Deployment](/examples/three-node-deployment): a worked server plus
  two-client example of both modes.
- [Configuration](/user-guide/configuration): center YAML fields and Flower settings.
- [Pipeline Workflow](/user-guide/design/workflow): screening stages executed
  during the deployment.
