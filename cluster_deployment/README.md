# Cluster Deployment

Scripts and documentation for running the federated GWAS pipeline on a **three-node cluster** (1 SuperLink + 2 SuperNodes). This is the script deployment mode. The equivalent CLI mode is `fedgwas-deploy` (see [examples/04three-node-deployment](../examples/04three-node-deployment/README.md)).

Both modes wrap the same Flower runtime. Use these scripts when you want shell wrappers from the repository; use `fedgwas-deploy` when you prefer package commands on each node.

## Documentation

**[CLUSTER_USER_GUIDE.md](docs/CLUSTER_USER_GUIDE.md)** — topology, code sync, node setup, data placement, start/stop, results, FAQ.

## Scripts (`scripts/`)

| Script | Purpose |
|--------|---------|
| `setup-cluster-node.sh` | PLINK + Python 3.11 `fedgwas` env on each node |
| `cluster-config-template.sh` | Example IPs / ports (source and override) |
| `cluster-start-server.sh` | SuperLink on server (9092 / 9093) |
| `cluster-start-client.sh` | SuperNode on a client (9094 / 9095) |
| `cluster-run-app.sh` | `flwr run . local-deployment` with scale / config-path |
| `cluster-stop-all.sh` | Stop SuperLink / SuperNode / stray `flwr run` on **this** node |
| `cluster-status.sh` | List Flower processes and listening ports |
| `check-cluster-clock-skew.sh` | Fail if \|skew\| > 500 ms vs server (Flower 1.19 auth) |
| `cluster-verify-data.sh` | PLINK `.bed/.bim/.fam` smoke test (`--scale`, `--client-id`) |
| `cluster-diagnose.sh` | Optional health summary (processes, data, reachability) |
| `download-flwr-wheels.sh` | Offline `flwr==1.19.0` wheels when PyPI is unreachable (dev machine) |

Run all scripts from the **repository root** (`cd ~/Fed-GWAS`).

## Quick start

1. Sync the same release tag to all three nodes
2. `bash cluster_deployment/scripts/setup-cluster-node.sh` on each node
3. Copy each center’s PLINK data to the matching client. SuperNodes connect through `--server-ip`, not `flower.server_address`.
4. Verify data: `cluster-verify-data.sh --scale tiny --client-id 1` (client 1), etc.
5. Start SuperLink → SuperNodes → `cluster-run-app.sh`

```bash
export PATH="/root/miniconda3/envs/fedgwas/bin:$PATH"
cd ~/Fed-GWAS

# Server
nohup cluster_deployment/scripts/cluster-start-server.sh --server-ip 0.0.0.0 \
  > /tmp/superlink.log 2>&1 &

# Clients (after editing server_address in center configs)
cluster_deployment/scripts/cluster-start-client.sh \
  --server-ip 192.168.1.88 --client-id 1 \
  --config experiments/correctness/tiny_even/configs/center_1/config.yaml

cluster_deployment/scripts/cluster-start-client.sh \
  --server-ip 192.168.1.88 --client-id 2 \
  --config experiments/correctness/tiny_even/configs/center_2/config.yaml

# App (server, after both SuperNodes registered)
cluster_deployment/scripts/cluster-run-app.sh \
  --server-ip 192.168.1.88 --rounds 20 --scale tiny
```

Optional wall-time + peak RSS log (server, after SuperLink + SuperNodes are up):

```bash
TIME_LOG=experiments/correctness/tiny_even/results_2/server_app_matpool_tiny_even_time.txt
/usr/bin/time -v -o "${TIME_LOG}" \
  cluster_deployment/scripts/cluster-run-app.sh \
  --server-ip 192.168.1.88 --rounds 20 --scale tiny

python experiments/tools/collect_run_metrics.py \
  experiments/correctness/tiny_even/results_2 \
  --time-file "${TIME_LOG}" --label matpool_tiny_even
```

## Default node layout

| Role | Example IP | Ports |
|------|------------|-------|
| Server | 192.168.1.88 | 9092 (Fleet), 9093 (Exec) |
| Client 1 | 192.168.1.89 | 9094 |
| Client 2 | 192.168.1.90 | 9095 |

Source `cluster_deployment/scripts/cluster-config-template.sh` and override IPs as needed.

## Results paths

| Scale | Results directory |
|-------|-------------------|
| tiny | `experiments/correctness/tiny_even/results_2/` |
| small | `experiments/performance/small_even/results/` |
| medium | `experiments/performance/medium_even/results/` |

After a run, rsync client `center_*/logs/` to the server and run `experiments/tools/collect_run_metrics.py`. See the user guide §4.

Performance scales: [experiments/performance/scales.yaml](../experiments/performance/scales.yaml); see `small_even/README.md` and `medium_even/README.md`.
