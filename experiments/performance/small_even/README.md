# Performance Experiment: Small Scale, Even Partition

## Configuration

- **Scale**: Small (2,000 samples, 20,000 SNPs)
- **Clients**: 2
- **Partition Strategy**: Even
- **Expected Runtime**: ~60 minutes on Matpool 3-node cluster

## 1. Generate Data (local machine or server)

```bash
python pipeline/simulation/simulated_data/generate_synthetic_data.py \
  --scale small \
  --partition-strategy even \
  --seed 42 \
  --output-dir experiments/performance/small_even/data
```

## 2. Deploy PLINK files to clients

Copy per-center `.bed/.bim/.fam` to each client node:

- Client 1: `data/small/center_1/small_center_1.*`
- Client 2: `data/small/center_2/small_center_2.*`

## 3. Cluster run

**Server** (after SuperLink + both SuperNodes are up):

```bash
# §2.1
nohup cluster_deployment/scripts/cluster-start-server.sh --server-ip 0.0.0.0 \
  > /tmp/superlink.log 2>&1 &

# §2.3
cluster_deployment/scripts/cluster-run-app.sh \
  --server-ip 192.168.1.88 \
  --scale small \
  --rounds 20
```

**Clients** (update `flower.server_address` to server IP first):

```bash
# Client 1
cluster_deployment/scripts/cluster-start-client.sh \
  --server-ip 192.168.1.88 --client-id 1 \
  --config experiments/performance/small_even/configs/center_1/config.yaml

# Client 2
cluster_deployment/scripts/cluster-start-client.sh \
  --server-ip 192.168.1.88 --client-id 2 \
  --config experiments/performance/small_even/configs/center_2/config.yaml
```

See `cluster_deployment/docs/CLUSTER_USER_GUIDE.md` §5 for full workflow.
