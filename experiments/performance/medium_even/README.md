# Performance Experiment: Medium Scale, Even Partition

## Configuration

- **Scale**: Medium (5,000 samples, 50,000 SNPs)
- **Clients**: 2
- **Partition Strategy**: Even
- **Expected Runtime**: ~240 minutes on Matpool 3-node cluster

## 1. Generate Data (local machine)

```bash
python pipeline/simulation/simulated_data/generate_synthetic_data.py \
  --scale medium \
  --partition-strategy even \
  --seed 42 \
  --output-dir experiments/performance/medium_even/data
```

## 2. Deploy PLINK files to clients (required)

Data is in `.gitignore` — **each client node** must have its own copy:

- Client 1: `data/medium/center_1/medium_center_1.{bed,bim,fam}` (~30 MB)
- Client 2: `data/medium/center_2/medium_center_2.{bed,bim,fam}`

```bash
# On dev machine (repo root)
tar -czf /tmp/medium_center_1.tar.gz \
  -C experiments/performance/medium_even/data/medium center_1
tar -czf /tmp/medium_center_2.tar.gz \
  -C experiments/performance/medium_even/data/medium center_2

scp /tmp/medium_center_1.tar.gz root@<client1>:~/
scp /tmp/medium_center_2.tar.gz root@<client2>:~/
```

**Client 1:**

```bash
mkdir -p /Fed-GWAS/experiments/performance/medium_even/data/medium
tar -xzf ~/medium_center_1.tar.gz \
  -C /Fed-GWAS/experiments/performance/medium_even/data/medium
ls /Fed-GWAS/experiments/performance/medium_even/data/medium/center_1/medium_center_1.bed
plink --bfile /Fed-GWAS/experiments/performance/medium_even/data/medium/center_1/medium_center_1 --freq --out /tmp/plink_check
```

**Client 2:** same with `center_2` / `medium_center_2`.

## 3. Cluster run

**SuperLink (server):**

```bash
nohup cluster_deployment/scripts/cluster-start-server.sh --server-ip 0.0.0.0 \
  > /tmp/superlink.log 2>&1 &
```

**SuperNodes (clients):**

```bash
cluster_deployment/scripts/cluster-start-client.sh \
  --server-ip <SERVER_IP> --client-id 1 \
  --config experiments/performance/medium_even/configs/center_1/config.yaml

cluster_deployment/scripts/cluster-start-client.sh \
  --server-ip <SERVER_IP> --client-id 2 \
  --config experiments/performance/medium_even/configs/center_2/config.yaml
```

**App (server):**

```bash
cluster_deployment/scripts/cluster-run-app.sh \
  --server-ip <SERVER_IP> \
  --rounds 20 \
  --scale medium
```

Verify data before run (on **each client**, pass its id):

```bash
# Client 1
bash cluster_deployment/scripts/cluster-verify-data.sh --scale medium --client-id 1
# Client 2
bash cluster_deployment/scripts/cluster-verify-data.sh --scale medium --client-id 2
```
