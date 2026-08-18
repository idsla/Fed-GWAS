---
title: Script Deployment
slug: /federated-deployment/script-deployment
---

# Script Deployment

Use the `cluster_deployment/scripts/` wrappers from a repository checkout to
start the federated server and clients. This is one of two equivalent
deployment modes. For package commands, see [CLI Deployment](cli-deployment.md).

Complete [Installation and Overview](installation.md) on every node first. Run
the scripts from the **repository root**. Pick one mode for a given run; do not
mix script and CLI commands.

The scripts wrap the same Flower runtime as the CLI:

- server: `flower-superlink` and `flwr run . local-deployment`
- each client: `flower-supernode`

A longer operator guide with topology and FAQ lives in
[`cluster_deployment/docs/CLUSTER_USER_GUIDE.md`](https://github.com/idsla/Fed-GWAS/blob/main/cluster_deployment/docs/CLUSTER_USER_GUIDE.md).

## Scripts

| Script | Purpose |
| --- | --- |
| `setup-cluster-node.sh` | Install PLINK and a Python 3.11 environment on a node |
| `cluster-start-server.sh` | Start SuperLink on the server (`9092` / `9093`) |
| `cluster-start-client.sh` | Start a SuperNode on one client |
| `cluster-run-app.sh` | Submit `flwr run . local-deployment` |
| `cluster-status.sh` | List Flower processes and listening ports on this node |
| `cluster-stop-all.sh` | Stop SuperLink / SuperNode / stray `flwr run` on this node |
| `check-cluster-clock-skew.sh` | Compare this node's clock with the server |
| `cluster-verify-data.sh` | Check that this client's PLINK files exist |
| `cluster-diagnose.sh` | Optional process, data, and reachability summary |

If the scripts are not executable, prefix them with `bash`.

## Server

Start SuperLink:

```bash
cd Fed-GWAS
export PATH="/path/to/fedgwas/env/bin:$PATH"

nohup cluster_deployment/scripts/cluster-start-server.sh --server-ip 0.0.0.0 \
  > /tmp/superlink.log 2>&1 &

sleep 2
tail -20 /tmp/superlink.log
ss -tln | grep -E '9092|9093' || netstat -tuln | grep -E '9092|9093'
```

`--server-ip 0.0.0.0` binds Fleet API `9092` and Exec API `9093` on all
interfaces. Keep SuperLink running while clients connect and while the app
runs.

## Clients

Start one SuperNode per participating client after SuperLink is up.
`--client-id` is one-based (`1`, `2`, …). Default ClientAppIo ports are
`9094` for client 1 and `9095` for client 2; override with `--port`.

Verify data on client *k*:

```bash
bash cluster_deployment/scripts/cluster-verify-data.sh \
  --scale tiny \
  --client-id <k>
```

Start client *k*:

```bash
cluster_deployment/scripts/cluster-start-client.sh \
  --server-ip <server-ip> \
  --client-id <k> \
  --config configs/center_k.yaml
```

Two-client example:

```bash
cluster_deployment/scripts/cluster-start-client.sh \
  --server-ip 192.168.1.88 \
  --client-id 1 \
  --config experiments/correctness/tiny_even/configs/center_1/config.yaml

cluster_deployment/scripts/cluster-start-client.sh \
  --server-ip 192.168.1.88 \
  --client-id 2 \
  --config experiments/correctness/tiny_even/configs/center_2/config.yaml
```

Background start:

```bash
nohup cluster_deployment/scripts/cluster-start-client.sh \
  --server-ip <server-ip> \
  --client-id 1 \
  --config experiments/correctness/tiny_even/configs/center_1/config.yaml \
  --port 9094 > /tmp/supernode1.log 2>&1 &
```

Each client should hold only its own center config and PLINK files. The
SuperNode connects through `--server-ip`, not through `flower.server_address`
in the YAML.

If a SuperNode log stops at `ClientAppIo gRPC server`, that is usually idle
waiting for `flwr run`. Stage progress appears under the center `output.log_dir`.

### Clock skew

Flower 1.19 rejects SuperNode messages when UTC clocks differ too much. From
each client:

```bash
bash cluster_deployment/scripts/check-cluster-clock-skew.sh <server-ip>
```

If the node cannot change system time, the start-client script supports a
process-local offset when `faketime` is installed:

```bash
export FEDGWAS_FAKETIME="-2s"   # use the value printed by the skew check
cluster_deployment/scripts/cluster-start-client.sh \
  --server-ip <server-ip> \
  --client-id 1 \
  --config configs/center_1.yaml
unset FEDGWAS_FAKETIME
```

## Submit the run

After SuperLink is running and every SuperNode has registered, submit the app
from the server:

```bash
cluster_deployment/scripts/cluster-run-app.sh \
  --server-ip <server-ip> \
  --rounds 20 \
  --scale tiny
```

`--scale` accepts `tiny`, `small`, or `medium`. Pass an explicit config
directory instead:

```bash
cluster_deployment/scripts/cluster-run-app.sh \
  --server-ip <server-ip> \
  --rounds 80 \
  --config-path experiments/performance/small_even/configs
```

`cluster-run-app.sh` temporarily patches the deployment federation addresses in
`pyproject.toml`, runs `flwr run . local-deployment`, and restores the file.

## Outputs

Tiny-scale repository configs commonly write under:

```text
experiments/correctness/tiny_even/results_2/
  server/
  center_1/
  center_2/
```

Copy `center_*` folders from clients to the server if the nodes do not share
storage. Then collect run metadata if needed:

```bash
python experiments/tools/collect_run_metrics.py \
  experiments/correctness/tiny_even/results_2 \
  --label three_node_tiny
```

## Stop

The scripts only affect the current node. Run stop and status on the server and
on every client:

```bash
cluster_deployment/scripts/cluster-stop-all.sh
cluster_deployment/scripts/cluster-status.sh
```

`cluster-stop-all.sh` kills `flower-superlink`, `flower-supernode`, and stray
`flwr run` processes on that node.

## Troubleshooting

| Symptom | Likely cause | Fix |
| --- | --- | --- |
| `cluster_deployment/scripts: No such file` | Not in the repository root | `cd` to the FedGWAS checkout that contains `cluster_deployment/` |
| `Unknown arg: --scale` | Stale `cluster-run-app.sh` | Sync the current `cluster_deployment/scripts/` to the server |
| SuperNode exits immediately | SuperLink down, port in use, or clock skew | Check `/tmp/supernode*.log` and `check-cluster-clock-skew.sh` |
| `config not found` | Relative path resolved from the wrong directory | Run from repository root or pass an absolute `--config` |
| Clients never join | Port `9092` closed, or server cannot reach ClientAppIo | Open `9092`–`9093+k` between hosts |
| `plink` not found | PLINK missing on the client | Run `setup-cluster-node.sh` or install PLINK on `PATH` |

Optional health summary:

```bash
bash cluster_deployment/scripts/cluster-diagnose.sh <server-ip>
```

Before a new run: stop old processes on every node, re-check clock skew, start
SuperLink, start all SuperNodes, then `cluster-run-app.sh`.
