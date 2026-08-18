# Three-Node Deployment

This example runs a three-node cluster (1 SuperLink server + 2 SuperNode clients)
on the tiny correctness experiment.

FedGWAS provides two equivalent deployment modes. Both wrap the same Flower
runtime (`flower-superlink`, `flower-supernode`, `flwr run . local-deployment`).

| Mode | Entry point |
| --- | --- |
| CLI | `fedgwas-deploy` |
| Script deployment | `cluster_deployment/scripts/` |

The scripts live at the repository root:

```text
cluster_deployment/scripts/*.sh
cluster_deployment/docs/CLUSTER_USER_GUIDE.md
```

Full topology, clock-skew, data placement, and teardown notes:
[cluster_deployment/docs/CLUSTER_USER_GUIDE.md](../../cluster_deployment/docs/CLUSTER_USER_GUIDE.md).

## Shared setup

On every node, clone the same revision and install into a Python 3.11 environment:

```bash
git clone https://github.com/idsla/Fed-GWAS.git
cd Fed-GWAS
python -m pip install -e .
python -m pip install "flwr[simulation]>=1.19.0,<1.20"
```

Or use the node setup script (installs PLINK and the `fedgwas` env):

```bash
bash cluster_deployment/scripts/setup-cluster-node.sh
```

Place each center's PLINK files only on that client, then edit
`flower.server_address` in the matching center YAML if you still use that field
for bookkeeping. SuperNodes connect through `--server-ip` / `--server`, not
through the YAML address.

Example hosts used below: server `192.168.1.88`, client 1 `192.168.1.89`,
client 2 `192.168.1.90`.

## Script deployment

Run from the repository root after `chmod +x cluster_deployment/scripts/*.sh`
if needed. Prefix with `bash` if the scripts are not executable.

**Server**

```bash
export PATH="/root/miniconda3/envs/fedgwas/bin:$PATH"
cd ~/Fed-GWAS

nohup cluster_deployment/scripts/cluster-start-server.sh --server-ip 0.0.0.0 \
  > /tmp/superlink.log 2>&1 &
```

**Client 1**

```bash
cluster_deployment/scripts/cluster-start-client.sh \
  --server-ip 192.168.1.88 \
  --client-id 1 \
  --config experiments/correctness/tiny_even/configs/center_1/config.yaml
```

**Client 2**

```bash
cluster_deployment/scripts/cluster-start-client.sh \
  --server-ip 192.168.1.88 \
  --client-id 2 \
  --config experiments/correctness/tiny_even/configs/center_2/config.yaml
```

**Server, after both SuperNodes are registered**

```bash
cluster_deployment/scripts/cluster-run-app.sh \
  --server-ip 192.168.1.88 \
  --rounds 20 \
  --scale tiny
```

**Stop on each node**

```bash
cluster_deployment/scripts/cluster-stop-all.sh
cluster_deployment/scripts/cluster-status.sh
```

## CLI (`fedgwas-deploy`)

**Server**

```bash
fedgwas-deploy server start --host 0.0.0.0 --daemon --log-file /tmp/superlink.log
```

**Client 1**

```bash
fedgwas-deploy client start \
  --server 192.168.1.88 \
  --center-id 1 \
  --config experiments/correctness/tiny_even/configs/center_1/config.yaml \
  --daemon \
  --log-file /tmp/supernode1.log
```

**Client 2**

```bash
fedgwas-deploy client start \
  --server 192.168.1.88 \
  --center-id 2 \
  --config experiments/correctness/tiny_even/configs/center_2/config.yaml \
  --daemon \
  --log-file /tmp/supernode2.log
```

**Server, after both SuperNodes are registered**

```bash
fedgwas-deploy server run \
  --server 192.168.1.88 \
  --rounds 20 \
  --scale tiny
```

**Stop on each node**

```bash
# server node
fedgwas-deploy server stop

# client nodes
fedgwas-deploy client stop
```

## Command mapping

| Cluster script | `fedgwas-deploy` |
| --- | --- |
| `cluster-start-server.sh --server-ip 0.0.0.0` | `server start --host 0.0.0.0` |
| `cluster-start-client.sh --server-ip <ip> --client-id 1 --config <path>` | `client start --server <ip> --center-id 1 --config <path>` |
| `cluster-run-app.sh --server-ip <ip> --rounds 20 --scale tiny` | `server run --server <ip> --rounds 20 --scale tiny` |
| `cluster-status.sh` | `server status` or `client status` |
| `cluster-stop-all.sh` | `server stop` or `client stop` |

Do not mix start/stop tools across nodes in the same run. Pick one mode for the
whole federation.

## Results

Tiny-scale outputs:

```text
experiments/correctness/tiny_even/results_2/
```

Copy `center_1/` and `center_2/` back to the server if the nodes do not share
storage, then:

```bash
python experiments/tools/collect_run_metrics.py \
  experiments/correctness/tiny_even/results_2 \
  --label three_node_tiny
```
