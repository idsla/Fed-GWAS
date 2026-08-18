# FedGWAS Three-Node Cluster — User Guide

This guide explains how to deploy and run **synthetic-data federated GWAS** on a **three-node Linux cluster** (1 server + 2 clients): install dependencies, start Flower SuperLink/SuperNodes, run the app, and locate **wall time, peak memory, and per-stage outputs**. It applies to any set of reachable hosts (cloud VMs, on-prem servers, or managed platforms such as Matpool).

### Cluster topology

FedGWAS uses Flower's **SuperLink / SuperNode** layout:

```
                    ┌─────────────────────────────┐
                    │  Server (192.168.1.88)        │
                    │  SuperLink :9092 / :9093      │
                    │  flwr run (Exec API)          │
                    └──────────┬──────────┬─────────┘
                               │          │
              SuperNode        │          │        SuperNode
              (partition 0)    │          │        (partition 1)
                               ▼          ▼
                    ┌──────────────┐  ┌──────────────┐
                    │ Client 1     │  │ Client 2     │
                    │ :9094        │  │ :9095        │
                    │ PLINK data   │  │ PLINK data   │
                    └──────────────┘  └──────────────┘
```

Each client holds **only its own center data**; the server orchestrates rounds and aggregates results. All three nodes run the **same code version** and a shared Python 3.11 `fedgwas` environment.

### Network and connectivity

| Requirement | Detail |
|-------------|--------|
| Reachability | All nodes on a **shared internal network** (same VPC/subnet, LAN, or routed peer network). Hostnames or static IPs must resolve consistently from every node. |
| Client → Server | TCP **9092** — SuperLink Fleet API (`flower.server_address` in client configs) |
| Server → Clients | TCP **9094**, **9095** — ClientAppIO gRPC (server must reach each client during rounds) |
| Server (local) | TCP **9093** — Exec API for `flwr run` on the server node |
| Firewall | Open **9092–9095** between nodes (security groups, `iptables`, cloud ACLs). Insecure mode (no TLS) is used by default in the bundled scripts. |
| Clock sync | UTC clocks within **\|skew\| < 500 ms** across all nodes (Flower 1.19 auth). See §2.0. |
| Verification | From each client: `ping 192.168.1.88` and test port 9092. From server: test 9094/9095 on clients (§2.2). Optional: `cluster_deployment/scripts/cluster-diagnose.sh`. |

**Node layout** (matches `cluster-config-template.sh`; replace IPs with your actual assignments):

| Role | Example IP | Ports |
|------|------------|-------|
| Server (SuperLink + App) | 192.168.1.88 | 9092 (Fleet), 9093 (Exec) |
| Client 1 | 192.168.1.89 | 9094 |
| Client 2 | 192.168.1.90 | 9095 |

Throughout this guide, **192.168.1.88** is the server, **192.168.1.89** is client 1, and **192.168.1.90** is client 2.

---

## 1. Cluster setup

### 1.1 Provision three nodes

Create **3 Linux VMs or servers** (recommended: 8 GB RAM / 4 vCPU / 50 GB disk each) on a network where all nodes can reach each other (same VPC/subnet or equivalent). On managed platforms (e.g. Matpool), use the provider console to create and label nodes as server / client1 / client2.

> **Screenshot 1:** Console or inventory showing three nodes labeled server / client1 / client2.

### 1.2 Get the code (every node)

Fed-GWAS is published as an **open-source repository**. On each cluster node, clone the public repo and check out a **release tag** (recommended) or stable branch. Repeat on all three nodes (server, client 1, client 2).

**Method A — `git clone` on each node (recommended):**

```bash
# On server, client 1, and client 2 (replace URL/tag with the public repo and release)
git clone <PUBLIC_REPO_URL> ~/Fed-GWAS
cd ~/Fed-GWAS
git checkout <release-tag>   # e.g. v0.1.0 — see GitHub Releases

ls ~/Fed-GWAS/pyproject.toml ~/Fed-GWAS/pipeline ~/Fed-GWAS/plink/plink_linux
```

Requires **outbound HTTPS/git** from each node. Use the same tag on all three nodes so code versions match.

**Method B — rsync from a dev machine (optional):**

Use when nodes cannot reach the public git remote but your workstation can.

```bash
# On dev machine (replace user@IP with SSH user and internal IP)
DEV=/path/to/Fed-GWAS

rsync -av --exclude='.git' --exclude='__pycache__' --exclude='.venv' \
  --exclude='experiments/*/results' --exclude='experiments/*/results_*' \
  "${DEV}/" user@192.168.1.88:~/Fed-GWAS/
# Repeat for 192.168.1.89 and 192.168.1.80
```

**Method C — tar + scp or Jupyter Lab (optional):**

Use when neither git nor rsync is available (e.g. restricted managed platforms).

Package a release on a machine that has the repo:

```bash
cd /path/to/Fed-GWAS
git fetch origin
rm -f ./fedgwas-cluster.tar.gz
git archive --format=tar.gz \
  --prefix=Fed-GWAS/ \
  --output=./fedgwas-cluster.tar.gz \
  <release-tag>

ls -lh ./fedgwas-cluster.tar.gz
```

Copy `fedgwas-cluster.tar.gz` to each node (`scp`, provider file upload, or Jupyter Lab), then extract:

```bash
cd ~
tar -xzf fedgwas-cluster.tar.gz
ls ~/Fed-GWAS/pyproject.toml ~/Fed-GWAS/pipeline
```

Repeat on all three nodes. Verify on each:

```bash
ls ~/Fed-GWAS/pyproject.toml ~/Fed-GWAS/pipeline ~/Fed-GWAS/plink/plink_linux
```

### 1.3 Install Fed-GWAS (fedgwas) and dependencies (every node)

Install the **same environment on all three nodes** (server, client 1, client 2). Requires **Python 3.11+**, Flower **1.19.x**, and PLINK on `PATH`. Use a dedicated env named **`fedgwas`** (conda or venv) — do not install into the system Python or an unrelated conda base env.

**Pre-install checks (must pass):**

```bash
/root/miniconda3/bin/conda create -n fedgwas python=3.11 -y   # first time
export PATH="/root/miniconda3/envs/fedgwas/bin:$PATH"
cd /Fed-GWAS

which python pip
python --version    # must be 3.11.x (not 3.8/3.9/3.10)
pip --version       # recommend >= 24.0
```

Confirm `python --version` reports **3.11.x** before installing dependencies.

**System packages (once per node):**

```bash
apt-get update
apt-get install -y git wget curl unzip build-essential time
```

**Method A — setup script (recommended):**

```bash
/root/miniconda3/bin/conda create -n fedgwas python=3.11 -y   # skip if exists
cd /Fed-GWAS   # or ~/Fed-GWAS, depending on extract path
bash cluster_deployment/scripts/setup-cluster-node.sh
```

The script installs PLINK → **pip wheels from PyPI** → **`flwr[simulation] 1.19.x`** (do not install 1.26+) → `pip install -e .`. If PyPI is unreachable, set `PIP_INDEX_URL` to a mirror your network allows, or use offline wheels (below).

**Flower version:** the project pins **1.19.x**. `pip install -U flwr` upgrades to 1.30+ and breaks `pyproject.toml` federation config. Verify: `flwr --version` should show `1.19.x`.

**Method B — manual pip wheels (equivalent to the script):**

```bash
export PATH="/root/miniconda3/envs/fedgwas/bin:$PATH"
cd /Fed-GWAS

pip install --upgrade "pip>=24" setuptools wheel

export PIP_PREFER_BINARY=1
pip install --prefer-binary numpy pandas scipy pyarrow

pip install --prefer-binary "flwr[simulation]==1.19.0" tomli tomli-w psutil

pip install --prefer-binary -e .
```

If Miniconda is not under `/root/miniconda3`, point `PATH` at your `.../envs/fedgwas/bin`. If PyPI is blocked, set `PIP_INDEX_URL` to an accessible mirror before running these commands.

**Verification (all three nodes):**

```bash
export PATH="/root/miniconda3/envs/fedgwas/bin:$PATH"
pip show fedgwas
python -c "import pipeline; import flwr; print('fedgwas OK')"
flwr --version    # must be 1.19.x (not 1.26+ / 1.30+)
plink --version
/usr/bin/time -v true 2>&1 | head -1
```

> **Screenshot 2:** Terminal showing successful `pip show fedgwas` and `python -c "import pipeline"`.

**Troubleshooting:**

| Symptom | Fix |
|---------|-----|
| `No matching distribution found for flwr<1.20,>=1.19.0` with `from versions: ... 1.11.1` | **Wrong env / Python ≠ 3.11:** `export PATH=/root/miniconda3/envs/fedgwas/bin:$PATH`, confirm `python --version` is 3.11.x, then `pip install -U "pip>=24"` and reinstall flwr |
| `PS1: unbound variable` | Use latest `setup-cluster-node.sh`, or Method B |
| `can't find Rust compiler` | Ensure `pip>=24` and `--prefer-binary` |
| conda stuck on `Solving environment` | **Do not conda-install Python packages**; use Method B |
| `flwr: command not found` | `export PATH=.../envs/fedgwas/bin:$PATH` first |
| `No matching distribution found` / `from versions: none` | **pip cannot reach PyPI**; see PyPI connectivity below |
| `Invalid SuperLink connection format` / Flower migration | **Flower too new** (1.30+); downgrade to 1.19.x (below) |
| `Invalid timestamp` / SuperNode exits immediately | Clock skew; run `check-cluster-clock-skew.sh`; if NTP is unavailable use **`FEDGWAS_FAKETIME` + faketime** (§2.0) |
| `ntpdate: step-systime: Operation not permitted` | Container or restricted VM cannot change system clock; **restart nodes from provider console** or enable platform NTP; do not force `date -s` in a web terminal |
| SuperNode idle / other client already in round1 | Check `results_2/center_1/logs/`; on server test reachability to client1:9094 (bash/python below) |
| `Invalid out message` (ClientApp running, SuperNode errors) | Same as above; `cluster-start-client.sh` sets `FEDGWAS_FLR_REPLY_SKEW_SEC=2` — run **`pip install -e .`** and restart SuperNode |

**Downgrade Flower (all nodes if `flwr --version` is not 1.19.x):**

```bash
export PATH="/root/miniconda3/envs/fedgwas/bin:$PATH"
pip install --prefer-binary "flwr[simulation]==1.19.0"
rm -rf ~/.flwr
flwr --version

# If pyproject was corrupted, restore backup or re-upload pyproject.toml from the repo
cd /Fed-GWAS
test -f pyproject.toml.bak.cluster && mv -f pyproject.toml.bak.cluster pyproject.toml
```

**PyPI connectivity (`ERROR: from versions: none`):**

This means pip **could not fetch the index**, not that 1.19 is missing. Diagnose:

```bash
export PATH="/root/miniconda3/envs/fedgwas/bin:$PATH"
which pip python
python --version          # must be 3.11.x
curl -I https://pypi.org/simple/flwr/   # expect HTTP/2 200
```

If `curl` fails or times out, retry Method B with `PIP_INDEX_URL` set to a mirror your network can reach, or use offline wheels:

```bash
# Dev machine
cd /path/to/Fed-GWAS
bash cluster_deployment/scripts/download-flwr-wheels.sh
scp -r cluster_deployment/wheels root@<node>:/tmp/flwr-wheels

# On node (fedgwas + Python 3.11)
export PATH="/root/miniconda3/envs/fedgwas/bin:$PATH"
pip install --no-index --find-links=/tmp/flwr-wheels \
  "flwr[simulation]==1.19.0" tomli tomli-w psutil
pip install --prefer-binary -e .
```

### 1.4 Place experiment data (clients only)

```bash
# Client 1 (192.168.1.89)
ls ./experiments/correctness/tiny_even/data/tiny/center_1/tiny_center_1.{bed,bim,fam}

# Client 2 (192.168.1.90)
ls ./experiments/correctness/tiny_even/data/tiny/center_2/tiny_center_2.{bed,bim,fam}
```

Data is in `.gitignore` — copy it to each client separately (`git archive` does not include it).

**Without rsync: tar + scp or Jupyter upload (recommended)**

```bash
# Dev machine (Fed-GWAS repo root)
cd /path/to/Fed-GWAS

tar -czf /tmp/tiny_center_1.tar.gz \
  -C experiments/correctness/tiny_even/data/tiny center_1
tar -czf /tmp/tiny_center_2.tar.gz \
  -C experiments/correctness/tiny_even/data/tiny center_2

scp /tmp/tiny_center_1.tar.gz user@192.168.1.89:~/
scp /tmp/tiny_center_2.tar.gz user@192.168.1.90:~/
```

**Extract on client 1:**

```bash
mkdir -p ~/Fed-GWAS/experiments/correctness/tiny_even/data/tiny
tar -xzf ~/tiny_center_1.tar.gz \
  -C ~/Fed-GWAS/experiments/correctness/tiny_even/data/tiny
ls ~/Fed-GWAS/experiments/correctness/tiny_even/data/tiny/center_1/tiny_center_1.bed
```

**On client 2:** replace `center_1` / `tiny_center_1` with `center_2` / `tiny_center_2`.

You can also upload `tiny_center_*.tar.gz` via Jupyter Lab and `tar -xzf` on the target VM.

**With rsync (optional):**

```bash
rsync -av experiments/correctness/tiny_even/data/tiny/center_1/ \
  user@192.168.1.89:~/Fed-GWAS/experiments/correctness/tiny_even/data/tiny/center_1/
rsync -av experiments/correctness/tiny_even/data/tiny/center_2/ \
  user@192.168.1.90:~/Fed-GWAS/experiments/correctness/tiny_even/data/tiny/center_2/
```

> **Screenshot 3:** `ls` showing `Fed-GWAS` on all nodes and `.bed/.bim/.fam` on clients.

### 1.5 Configure client SuperLink address

Edit `flower.server_address` on both clients (Fleet API = server IP:9092):

```yaml
# experiments/correctness/tiny_even/configs/center_1/config.yaml
flower:
  server_address: "192.168.1.88:9092"
```

Apply the same change in `center_2/config.yaml`. On the server, `pyproject.toml` lists Exec API `192.168.1.88:9093` under `cluster-deployment`; `cluster-run-app.sh` patches this temporarily at run time.

> **Screenshot 4:** Editor showing `server_address` and the `cluster-deployment` section in `pyproject.toml`.

---

## 2. Start the federation (recommended order)

Run commands on the matching node. Paths assume `~/Fed-GWAS` or `/Fed-GWAS`.

### 2.0 Clock sync (required before SuperNode)

Flower 1.19 requires **client and server UTC clocks to agree** (typically **|skew| < 500 ms**). Large drift causes:

`UNAUTHENTICATED` / `Invalid timestamp`

**On all three VMs (server + both clients):**

```bash
timedatectl set-ntp true 2>/dev/null || true
systemctl restart systemd-timesyncd 2>/dev/null || true
apt-get install -y ntpdate 2>/dev/null || true
ntpdate -u ntp.aliyun.com || ntpdate -u pool.ntp.org
date -u
```

#### `ntpdate` → `step-systime: Operation not permitted`

On some **managed VMs, containers, or restricted shells**, even root cannot change the system clock (missing `CAP_SYS_TIME`). Fix time at the **host or platform layer**:

1. **Restart all three nodes** from your cloud/provider console (often realigns with the host clock)
2. If nodes have **direct SSH** (not a web terminal inside a container), retry `ntpdate` there
3. Confirm with your provider that **NTP / host clock sync** is enabled for the instances
4. After restart, on each client:

```bash
bash cluster_deployment/scripts/check-cluster-clock-skew.sh 192.168.1.88
```

You need **|skew| < 500 ms** before starting SuperNode. If the client is ~5 s **ahead** of the server, the **server clock is slow** — fix the server side, not only the client.

**Compare client vs server:**

```bash
cd /Fed-GWAS
bash cluster_deployment/scripts/check-cluster-clock-skew.sh 192.168.1.88
```

Expect **`OK: clocks aligned`** and exit code 0. **AHEAD** means the client is faster than the server — sync NTP on the **server** too.

Without SSH to the server, run `date -u +%s.%N` on both machines at the same time; the difference should be **< 0.5 s** (absolute).

#### Still cannot fix system time: `faketime`

When `check-cluster-clock-skew.sh` fails, it prints a suggested `FEDGWAS_FAKETIME`. This adjusts **only the SuperNode process**, not system time:

```bash
cd /Fed-GWAS
bash cluster_deployment/scripts/check-cluster-clock-skew.sh 192.168.1.88
# Example: Skew +5798 ms (client ahead) → script suggests:
#   export FEDGWAS_FAKETIME="-5.798s"

apt-get install -y faketime
export FEDGWAS_FAKETIME="-6s"   # use the value from the check script; may differ per client

cluster_deployment/scripts/cluster-start-client.sh \
  --server-ip 192.168.1.88 --client-id 2 \
  --config experiments/correctness/tiny_even/configs/center_2/config.yaml \
  --port 9095
```

Run the check separately on client 1 and client 2 if skew differs. Unset `FEDGWAS_FAKETIME` after the experiment.

### 2.1 Server: start SuperLink

On the **server**, start SuperLink with **`cluster-start-server.sh`**:

```bash
conda activate fedgwas
cd ~/Fed-GWAS
chmod +x cluster_deployment/scripts/*.sh
export PATH="/root/miniconda3/envs/fedgwas/bin:$PATH"

# Background (recommended)
nohup cluster_deployment/scripts/cluster-start-server.sh --server-ip 0.0.0.0 \
  > /tmp/superlink.log 2>&1 &

sleep 2
tail -20 /tmp/superlink.log
ss -tln | grep -E '9092|9093' || netstat -tuln | grep -E '9092|9093'
```

For debugging, run in the foreground (blocks the terminal; Ctrl+C to stop):

```bash
cluster_deployment/scripts/cluster-start-server.sh --server-ip 0.0.0.0
```

> **Screenshot 5:** Server `tail /tmp/superlink.log` showing SuperLink listening on 9092/9093.

### 2.2 Client 1 / Client 2: start SuperNode

Start one SuperNode per client after §2.0 passes and §2.1 is running.

**Client 1 (192.168.1.89):**

```bash
conda activate fedgwas
cd ~/Fed-GWAS
export PATH="/root/miniconda3/envs/fedgwas/bin:$PATH"

# Background (recommended)
nohup cluster_deployment/scripts/cluster-start-client.sh \
  --server-ip 192.168.1.88 \
  --client-id 1 \
  --config experiments/correctness/tiny_even/configs/center_1/config.yaml \
  --port 9094 > /tmp/supernode1.log 2>&1 &

sleep 2
tail -30 /tmp/supernode1.log
```

**Client 2 (192.168.1.90):**

```bash
conda activate fedgwas
cd ~/Fed-GWAS
export PATH="/root/miniconda3/envs/fedgwas/bin:$PATH"

# Background (recommended)
nohup cluster_deployment/scripts/cluster-start-client.sh \
  --server-ip 192.168.1.88 \
  --client-id 2 \
  --config experiments/correctness/tiny_even/configs/center_2/config.yaml \
  --port 9095 > /tmp/supernode2.log 2>&1 &

sleep 2
tail -30 /tmp/supernode2.log
cluster_deployment/scripts/cluster-status.sh
```

For first-time debugging, run **one client in the foreground** (omit `nohup` and `&`) until you see Flower output; then switch to background as above.

If `[1]+ Done` with an almost empty log, SuperNode **exited immediately**. Common causes: SuperLink not running, clock skew, port in use, or `flower-supernode` not on `PATH`.

**SuperNode stops at `ClientAppIo gRPC server` with no new lines:** usually **normal idle** (round logs go to `results_2/center_*/logs/`). If the server is in round 1 but `center_1/logs/` is empty, test connectivity from the **server**:

```bash
bash -c 'echo >/dev/tcp/192.168.1.89/9094' && echo OK || echo FAIL

python3 -c "import socket; socket.create_connection(('192.168.1.89', 9094), timeout=3).close(); print('OK')"
```

When using `faketime`, use the latest `cluster-start-client.sh` (sets `LD_PRELOAD` for ClientApp child processes).

> **Screenshot 6:** SuperNode logs on both clients; on server `ss -tlnp | grep 909` shows connections.

### 2.3 Server: run the app

**Prerequisites:** SuperLink running (§2.1), both SuperNodes connected (§2.2).

```bash
conda activate fedgwas
cd ~/Fed-GWAS
export PATH="/root/miniconda3/envs/fedgwas/bin:$PATH"

cluster_deployment/scripts/cluster-run-app.sh \
  --server-ip 192.168.1.88 \
  --rounds 200 \
  --scale tiny
```

`--scale` accepts `tiny` (default), `small`, or `medium` (see §5). You can also pass an explicit config directory:

```bash
cluster_deployment/scripts/cluster-run-app.sh \
  --server-ip 192.168.1.88 \
  --rounds 200 \
  --config-path experiments/performance/small_even/configs
```

For **small** or **medium**, update each client's `--config` to the matching `experiments/performance/*/configs/center_*/config.yaml` before §2.2.

**Optional — record wall time and peak RSS:**

```bash
/usr/bin/time -v -o experiments/correctness/tiny_even/results_2/server_app_matpool_tiny_even_time.txt \
  cluster_deployment/scripts/cluster-run-app.sh \
  --server-ip 192.168.1.88 \
  --rounds 80 \
  --scale tiny
```

Ensure `output.*_dir` in each config matches the results root (`results_2/` for tiny).

Confirm **`flwr --version` is 1.19.x** (§1.3). Flower 1.19 uses `[tool.flwr.federations.local-deployment]` in `pyproject.toml`. Do **not** install 1.26+ (forces migration to `~/.flwr/config.toml` and breaks this setup).

> **Screenshot 7:** `flwr run` stream showing stage progress (QC → KING → LR).

---

## 3. Results and metrics (tiny_even)

For **tiny_even**, outputs live under:

```
experiments/correctness/tiny_even/results_2/
├── center_1/          # client 1 (on 192.168.1.89)
│   ├── logs/          # per-round ClientApp logs
│   └── intermediate/  # pruned when retention.tier: standard finishes
├── center_2/          # client 2 (on 192.168.1.90)
│   ├── logs/
│   └── intermediate/
└── server/            # server-side aggregation outputs
```

After a successful run, copy `center_1/` and `center_2/` from the clients to the server so all artifacts sit under one tree.

**Wall time and peak memory (optional):**

If you wrapped §2.3 with `/usr/bin/time -v -o ...`, read:

```
experiments/correctness/tiny_even/results_2/server_app_matpool_tiny_even_time.txt
```

Key lines:

- `Elapsed (wall clock time)` — total server app duration
- `Maximum resident set size (kbytes)` — peak RSS

**Per-stage metrics (when monitoring is enabled in config):**

With `monitoring` enabled in `experiments/correctness/tiny_even/config.yaml`, the server merges these to the results root at the `done` stage:

- `stage_metrics.csv` — per-stage timing on the server
- `client_metrics.csv` — per-client timing/memory summaries

Round-by-round detail remains in `center_*/logs/`. With `retention.tier: standard`, large `intermediate/` files are removed automatically; logs and metric CSVs are kept.

> **Screenshot 8:** `ls experiments/correctness/tiny_even/results_2/` and sample lines from `server_app_*_time.txt` or `stage_metrics.csv`.

---

## 4. Standard teardown

When `flwr run` finishes, **Flower does not automatically stop SuperLink or SuperNodes**. That is normal for multi-run cluster workflows, but you should **explicitly stop processes on each node** to free ports (9092–9095) and memory.

### What stops by itself

| Component | After a successful run |
|-----------|-------------------------|
| **`flwr run` (app)** | Exits when the pipeline reaches `done` |
| **SuperLink** | **Keeps running** if you started it manually (§2.1) — stop it on the server when done |
| **SuperNodes (both clients)** | **Never** stopped automatically — clean up on each client node |

### Standard post-run checklist (all three nodes)

Run **once on the server, once on client 1, once on client 2** after §2.3 or §3 completes (success or failure):

```bash
cd /Fed-GWAS
export PATH="/root/miniconda3/envs/fedgwas/bin:$PATH"

cluster_deployment/scripts/cluster-stop-all.sh
cluster_deployment/scripts/cluster-status.sh
```

`cluster-stop-all.sh` kills `flower-superlink`, `flower-supernode`, and any stray `flwr run` on **that** node only. There is no remote shutdown from the server — **each node needs its own stop**.

Expected `cluster-status.sh` output after cleanup:

```
SuperLink:   Not running
SuperNodes:  Not running
App Runner:  Not running
```

### Optional: pull results before teardown

If you need logs or metrics for §3, sync client `logs/` to the server **before** stopping SuperNodes (while processes are still consistent), or copy files with `scp` / your provider's file transfer. Teardown does not delete `results_2/`; it only stops Flower processes.

### Restarting a new run

1. Run `cluster-stop-all.sh` on **all three** nodes (even if one node looks idle).
2. Re-check clock skew (§2.0) if nodes were left running for a long time.
3. Start SuperLink → both SuperNodes → app again (§2.1–§2.3).

Full restart workflow: §4 and §2.1–§2.3 above.

### Manual alternative

```bash
pkill -f "flower-superlink|flower-supernode|flwr run"
```

Prefer `cluster-stop-all.sh` for clearer status messages.

---

## 5. Appendix: scalability (multiple scales)

| Scale | Samples | SNPs | Config directory |
|-------|---------|------|------------------|
| Tiny | 500 | 5,000 | `experiments/correctness/tiny_even` (default) |
| Small | 2,000 | 20,000 | `experiments/performance/small_even` |
| Medium | 5,000 | 50,000 | `experiments/performance/medium_even` |

For each scale, repeat §2 with the matching client configs and `--scale` (or `--config-path`). Small/medium write under `experiments/performance/<scale>_even/results/` instead of `results_2/`.

Parameter definitions: `experiments/performance/scales.yaml`. See `small_even/README.md` and `medium_even/README.md` for data generation and cluster runs.



---

## 6. FAQ

| Symptom | Fix |
|---------|-----|
| Client cannot reach SuperLink | Open firewall/security group for 9092–9095; verify `server_address` is `192.168.1.88:9092` |
| `plink: not found` | Run `which plink` on all nodes; setup script installs PLINK |
| Empty results directory | Check `output.intermediate_dir` / `log_dir` in `config.yaml` match cwd-relative paths |
| No `stage_metrics.csv` | Confirm `monitoring` is enabled in config; otherwise use `*_time.txt` from `/usr/bin/time -v` |
| `Unknown arg: --scale` on server | Re-sync `cluster_deployment/scripts/cluster-run-app.sh` from the dev machine |
