---
title: Three-Node Federated Deployment
---

# Three-Node Deployment

This walkthrough runs the tiny correctness experiment on three machines.

| Role | Example host | Runs |
| --- | --- | --- |
| Server | `192.168.1.88` | SuperLink and `flwr run` |
| Center 1 | `192.168.1.89` | SuperNode for center 1 |
| Center 2 | `192.168.1.90` | SuperNode for center 2 |

Use the same FedGWAS code version on all three nodes.

FedGWAS has two equivalent deployment modes. Both wrap Flower SuperLink,
SuperNode, and `flwr run . local-deployment`. Pick one mode for the whole run.

| Mode | Entry point |
| --- | --- |
| [CLI](/user-guide/federated-deployment/cli-deployment) | `fedgwas-deploy` |
| [Script deployment](/user-guide/federated-deployment/script-deployment) | `cluster_deployment/scripts/` |

The repository example README is
[examples/04three-node-deployment](https://github.com/idsla/Fed-GWAS/tree/main/examples/04three-node-deployment).
The script deployment guide is
[cluster_deployment/docs/CLUSTER_USER_GUIDE.md](https://github.com/idsla/Fed-GWAS/blob/main/cluster_deployment/docs/CLUSTER_USER_GUIDE.md).

## 1. Prepare Each Node

Run on the server and both clients:

```bash
export PATH="/root/miniconda3/envs/fedgwas/bin:$PATH"
cd ~/Fed-GWAS

python -m pip install -e .
python -m pip install "flwr[simulation]>=1.19.0,<1.20"
```

Node setup script (optional; installs PLINK and the `fedgwas` env):

```bash
bash cluster_deployment/scripts/setup-cluster-node.sh
```

On each client, confirm PLINK is available:

```bash
plink --version
```

## 2. Script deployment

Run these commands from the repository root.

Start SuperLink on `192.168.1.88`:

```bash
nohup cluster_deployment/scripts/cluster-start-server.sh --server-ip 0.0.0.0 \
  > /tmp/superlink.log 2>&1 &
```

Start center 1 on `192.168.1.89`:

```bash
cluster_deployment/scripts/cluster-start-client.sh \
  --server-ip 192.168.1.88 \
  --client-id 1 \
  --config experiments/correctness/tiny_even/configs/center_1/config.yaml
```

Start center 2 on `192.168.1.90`:

```bash
cluster_deployment/scripts/cluster-start-client.sh \
  --server-ip 192.168.1.88 \
  --client-id 2 \
  --config experiments/correctness/tiny_even/configs/center_2/config.yaml
```

Submit the run from the server after both SuperNodes are registered:

```bash
cluster_deployment/scripts/cluster-run-app.sh \
  --server-ip 192.168.1.88 \
  --rounds 20 \
  --scale tiny
```

Stop processes on each node:

```bash
cluster_deployment/scripts/cluster-stop-all.sh
cluster_deployment/scripts/cluster-status.sh
```

## 3. CLI (`fedgwas-deploy`)

This mode is equivalent to script deployment. Do not mix the two modes in the same run.

Start SuperLink on `192.168.1.88`:

```bash
fedgwas-deploy server start \
  --host 0.0.0.0 \
  --daemon \
  --log-file /tmp/superlink.log
```

Start center 1 on `192.168.1.89`:

```bash
fedgwas-deploy client start \
  --server 192.168.1.88 \
  --center-id 1 \
  --config experiments/correctness/tiny_even/configs/center_1/config.yaml \
  --daemon \
  --log-file /tmp/supernode1.log
```

Start center 2 on `192.168.1.90`:

```bash
fedgwas-deploy client start \
  --server 192.168.1.88 \
  --center-id 2 \
  --config experiments/correctness/tiny_even/configs/center_2/config.yaml \
  --daemon \
  --log-file /tmp/supernode2.log
```

Submit the run from the server:

```bash
fedgwas-deploy server run \
  --server 192.168.1.88 \
  --rounds 20 \
  --scale tiny
```

The command submits:

```bash
flwr run . local-deployment --stream --run-config \
  'simulation=false num-server-rounds=20 config_path="experiments/correctness/tiny_even/configs"'
```

Stop processes on each node:

```bash
# server node
fedgwas-deploy server stop

# client nodes
fedgwas-deploy client stop
```

## 4. Collect Results

The tiny correctness configs write under:

```text
experiments/correctness/tiny_even/results_2/
```

After the run, copy client result folders back to the server if the nodes do not
share storage. Then collect run metadata:

```bash
python experiments/tools/collect_run_metrics.py \
  experiments/correctness/tiny_even/results_2 \
  --label three_node_tiny
```

## Command mapping

| Cluster script | `fedgwas-deploy` |
| --- | --- |
| `cluster-start-server.sh --server-ip 0.0.0.0` | `server start --host 0.0.0.0` |
| `cluster-start-client.sh --server-ip <ip> --client-id 1 --config <path>` | `client start --server <ip> --center-id 1 --config <path>` |
| `cluster-run-app.sh --server-ip <ip> --rounds 20 --scale tiny` | `server run --server <ip> --rounds 20 --scale tiny` |
| `cluster-status.sh` | `server status` or `client status` |
| `cluster-stop-all.sh` | `server stop` or `client stop` |
