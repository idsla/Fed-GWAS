---
title: CLI Deployment
slug: /federated-deployment/cli-deployment
---

# CLI Deployment

Use `fedgwas-deploy` to start the federated server and clients with package
commands. This is one of two equivalent deployment modes. For shell wrappers
from a repository checkout, see [Script Deployment](script-deployment.md).

Complete [Installation and Overview](installation.md) on every node first. Pick
one mode for a given run; do not mix CLI and script commands.

The CLI wraps:

- server: `flower-superlink` and `flwr run . local-deployment`
- each client: `flower-supernode`

## Server

Start SuperLink so clients can register:

```bash
fedgwas-deploy server start \
  --host 0.0.0.0 \
  --daemon \
  --log-file /tmp/superlink.log
```

Check that SuperLink is available:

```bash
fedgwas-deploy server check --network
fedgwas-deploy server status
```

`--host 0.0.0.0` binds the Fleet API (`9092`) and Exec API (`9093`) on all
interfaces so remote clients can connect. Keep SuperLink running while clients
start and while the app runs.

## Clients

Start one SuperNode per participating client after SuperLink is up.
`--center-id` is one-based. The CLI maps it to Flower's zero-based
`partition-id` and to ClientAppIo port `9093 + k` unless you pass `--port`.

For N clients, pass the same `--num-centers N` on every client, and set Flower
`options.num-supernodes` to N.

Check a client before starting it:

```bash
fedgwas-deploy client check \
  --server <server-ip> \
  --center-id <k> \
  --config configs/center_k.yaml \
  --network
```

Start client *k*:

```bash
fedgwas-deploy client start \
  --server <server-ip> \
  --center-id <k> \
  --num-centers <N> \
  --config configs/center_k.yaml \
  --daemon \
  --log-file /tmp/supernode_k.log
```

Two-client example:

```bash
fedgwas-deploy client start \
  --server 192.168.1.88 \
  --center-id 1 \
  --config experiments/correctness/tiny_even/configs/center_1/config.yaml \
  --daemon \
  --log-file /tmp/supernode1.log

fedgwas-deploy client start \
  --server 192.168.1.88 \
  --center-id 2 \
  --config experiments/correctness/tiny_even/configs/center_2/config.yaml \
  --daemon \
  --log-file /tmp/supernode2.log
```

| Center id | Partition | Default ClientAppIo port |
| --- | --- | --- |
| `1` | `0` | `9094` |
| `2` | `1` | `9095` |
| `k` | `k - 1` | `9093 + k` |

Each client needs its own center YAML and local PLINK files. SuperNodes connect
through `--server`, not through `flower.server_address` in the YAML.

If a SuperNode log stops at `ClientAppIo gRPC server`, that is usually idle
waiting for `flwr run`. Stage progress appears under the center `output.log_dir`.

## Submit the run

After SuperLink is running and every SuperNode has registered, submit the app
from the server:

```bash
fedgwas-deploy server run \
  --server <server-ip> \
  --rounds 20 \
  --scale tiny
```

Use an explicit config directory instead of a scale shortcut:

```bash
fedgwas-deploy server run \
  --server <server-ip> \
  --rounds 80 \
  --config-path experiments/performance/small_even/configs
```

`--server` patches the deployment federation Exec API address in
`pyproject.toml` for the duration of the run, then restores the file.

## Outputs

Server outputs follow the server config `output` paths, often:

```text
results/server/logs/
results/server/intermediate/
```

Each client writes to its center config paths:

```text
results/center_k/logs/
results/center_k/intermediate/
```

After a multi-machine run, copy client result folders to one tree before
centralized evaluation.

## Stop

Flower does not stop SuperLink or SuperNodes when `flwr run` exits. Stop
processes on **this node only**:

```bash
# server node
fedgwas-deploy server stop

# each client node
fedgwas-deploy client stop
```

Confirm with `fedgwas-deploy server status` or `fedgwas-deploy client status`.

## Equivalent Flower commands

`fedgwas-deploy server start` wraps:

```bash
flower-superlink \
  --insecure \
  --fleet-api-address 0.0.0.0:9092 \
  --exec-api-address 0.0.0.0:9093
```

`fedgwas-deploy client start --center-id 1 --num-centers 2` wraps:

```bash
flower-supernode \
  --insecure \
  --superlink <server-ip>:9092 \
  --clientappio-api-address 0.0.0.0:9094 \
  --node-config 'partition-id=0 num-partitions=2 config-file="..."'
```

`fedgwas-deploy server run --scale tiny` wraps:

```bash
flwr run . local-deployment --stream --run-config \
  'simulation=false num-server-rounds=20 config_path="experiments/correctness/tiny_even/configs"'
```

## Troubleshooting

| Symptom | Likely cause | Fix |
| --- | --- | --- |
| `flwr run` cannot connect | Exec API address is wrong | Pass `--server <server-ip>` or set the federation `address` to `<server-ip>:9093` |
| Clients never join | Fleet API not reachable | Open port `9092`; confirm `--server` is the SuperLink host |
| SuperNode exits immediately | SuperLink down, port in use, or clock skew | Check `/tmp/supernode*.log`, ports, and UTC clocks |
| `config not found` | Relative path resolved from the wrong directory | Run from the app directory or pass an absolute `--config` |
| `plink` not found | PLINK missing on the client | Install PLINK 1.9 / PLINK 2 and re-run `client check` |
| No client logs | Output path wrong, or the run has not reached that client | Check `output.log_dir` and that every SuperNode registered |

Before a new run: stop old processes on every node, confirm clocks and code
versions, start SuperLink, start all SuperNodes, then `server run`.
