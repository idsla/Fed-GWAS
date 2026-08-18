# Release Checklist

This document defines what to verify before tagging a release. Research artifacts (manuscript figures, historical `results/`, real-world raw data) stay in the repository but are **not** part of the installable `pipeline` package (`pip install -e .` / wheel only ships `pipeline/`).

## What ships

| Component | Location | Notes |
|-----------|----------|--------|
| Python package | `pipeline/` | Installed via `pip install -e .` / hatch wheel (`packages = ["pipeline"]`) |
| Default Flower config | `pyproject.toml` | `config_path = experiments/correctness/tiny_even/configs` |
| Correctness experiment | `experiments/correctness/tiny_even/` | Recommended smoke test |
| Cluster deployment | `cluster_deployment/` | Scripts + [CLUSTER_USER_GUIDE.md](cluster_deployment/docs/CLUSTER_USER_GUIDE.md) |
| Performance scales | `experiments/performance/` | `scales.yaml`, `small_even/`, `medium_even/` configs |
| Evaluation tools | `experiments/tools/evaluation/` | QC, LR, KING comparison vs baseline |

## Release verification path

### 1. Environment

```bash
uv sync --python 3.11
# PLINK 1.9+ on PATH
```

### 2. Data and baseline (tiny_even)

```bash
python pipeline/simulation/simulated_data/generate_synthetic_data.py \
  --scale tiny --partition-strategy even --seed 42 \
  --output-dir experiments/correctness/tiny_even/data

python experiments/tools/generate_baseline.py \
  experiments/correctness/tiny_even/config.yaml
```

### 3. Unit test

```bash
pytest tests/test_king_federated_unit.py -q
```

### 4. End-to-end simulation

```bash
flwr run . local-simulation --stream --run-config \
  'simulation=true num-server-rounds=100 config_path="experiments/correctness/tiny_even/configs"'
```

### 5. Evaluation (optional, after run completes)

```bash
python experiments/tools/evaluation/evaluate_all.py \
  experiments/correctness/tiny_even/results \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline \
  --king
```

## Paths by use case

| Tier | Purpose | Config path |
|------|---------|-------------|
| **Release / CI** | Fast correctness check | `experiments/correctness/tiny_even/configs` |
| **Local deployment** | Multi-process Flower | Same configs + `flwr run . local-deployment` |
| **Research** | 1000G / manuscript | `experiments/real_world/1000genomes/configs` (override via `--run-config`) |
| **Cluster (3-node)** | Matpool / on-prem SuperLink | `cluster_deployment/` + `flwr run . local-deployment` |

### Cluster smoke (optional, after tiny simulation passes)

On each node: `bash cluster_deployment/scripts/setup-cluster-node.sh`, then follow [cluster_deployment/docs/CLUSTER_USER_GUIDE.md](cluster_deployment/docs/CLUSTER_USER_GUIDE.md). Quick data check on client 1:

```bash
bash cluster_deployment/scripts/cluster-verify-data.sh --scale tiny --client-id 1
```

## Run output retention (product tiers)

After a run completes, the server applies `retention` from the experiment `config.yaml` (default **standard**):

| Tier | Keeps | Prunes |
|------|--------|--------|
| **minimal** | PLINK science outputs, evaluation reports, merged `stage_metrics.csv` / `run_summary.*` | `intermediate/`, logs text, per-node metrics, network monitor, crypto seeds, KING `.pkl` |
| **standard** (default) | minimal + `{client}_log.txt`, per-node performance CSV/summary | `intermediate/`, network monitor, crypto, accumulator |
| **research** | Full tree (debug / paper reproduction) | Nothing |

Manual run:

```bash
python experiments/tools/apply_run_retention.py experiments/correctness/tiny_even/results_2 \
  --config-path experiments/correctness/tiny_even/configs
```

Use `--dry-run` to preview deletions. Writes `retention_manifest.json` under the results root.

## PyPI publishing

Python package releases are published by `.github/workflows/publish-pypi.yml` when a GitHub Release is published. The workflow uses PyPI Trusted Publishing through GitHub Actions OIDC, so no PyPI username, password, API token, or GitHub secret is required in the workflow.

Configure the trusted publisher once in PyPI with these values:

| Field | Value |
|-------|-------|
| PyPI project | `FedGWAS` |
| GitHub owner | `idsla` |
| GitHub repository | `Fed-GWAS` |
| Workflow name | `publish-pypi.yml` |
| Environment name | `pypi` |

Recommended release flow:

```bash
git tag v0.0.1
git push origin v0.0.1
```

Then create and publish a GitHub Release for the tag. The workflow builds `dist/`, checks it with `twine check`, and publishes the distributions to PyPI.

## Not in release scope

- `experiments/**/results/` (gitignored run outputs)
- `experiments/real_world/1000genomes/manuscript/` (figures, tables)
- `pipeline/simulation/`, `pipeline/visualization_archived/`, `pipeline/experiments.md` (local reference)
- Legacy server modules: `strategy.py`, `strategy_new.py`, `aggregator_qc.py`
- Internal dev notes (`plan.md`, `ISSUES_REPORT*.md` if present)

## Pre-tag sanity checks

```bash
# No dashboard remnants
rg -i dashboard --glob '!RELEASE.md' --glob '!*.plan.md' .

# Lockfile matches pyproject (no dash deps)
rg '^name = "dash' uv.lock  # should be empty
```

## Documentation

- Quick start: [README.md](README.md)
- Pipeline behavior and privacy: [CURRENT_VERSION.md](CURRENT_VERSION.md)
- Cluster (3-node): [cluster_deployment/docs/CLUSTER_USER_GUIDE.md](cluster_deployment/docs/CLUSTER_USER_GUIDE.md)
- Tiny experiment details: [experiments/correctness/tiny_even/README.md](experiments/correctness/tiny_even/README.md)
- Performance scales: `experiments/performance/scales.yaml`, `small_even/`, `medium_even/`
