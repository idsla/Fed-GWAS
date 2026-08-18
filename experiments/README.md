# Experiments Directory

Federated GWAS experiment configs, real-world study layouts, and post-run tooling. All paths are relative to the repository root.

## Layout

```
experiments/
├── README.md
├── correctness/
│   └── tiny_even/              # Default smoke test + cluster tiny benchmark
├── performance/
│   ├── scales.yaml             # tiny / small / medium matrix (2 clients each)
│   ├── small_even/             # 2k samples × 20k SNPs (Matpool appendix)
│   └── medium_even/            # 5k samples × 50k SNPs (Matpool appendix)
├── real_world/
│   ├── phenotype_generation/   # Simulated phenotypes on real genotypes
│   └── 1000genomes/            # 1kGP chr22 application experiment
└── tools/
    ├── collect_run_metrics.py  # Wall time / RSS / stage CSV summary
    ├── apply_run_retention.py  # Manual retention tier application
    ├── generate_baseline.py    # Centralized PLINK baseline
    ├── summarize_scalability_table.py
    └── evaluation/             # QC, KING, LR comparison vs baseline
```

Synthetic PLINK data is **not** committed (see `.gitignore`). Generate locally with `pipeline/simulation/simulated_data/generate_synthetic_data.py` (local-only path; not in the git tree).

Run outputs (`results/`, `results_2/`, `data/**/*.bed`) are gitignored.

---

## Experiment matrix

| Study | Directory | Samples (total) | SNPs | Clients | Role |
|-------|-----------|-----------------|------|---------|------|
| Correctness | `correctness/tiny_even` | 500 | 5,000 | 2 | CI / simulation / cluster smoke |
| Performance small | `performance/small_even` | 2,000 | 20,000 | 2 | Scalability |
| Performance medium | `performance/medium_even` | 5,000 | 50,000 | 2 | Scalability |
| Real-world | `real_world/1000genomes` | 1kGP subset | chr22 | 2 | Application |

Details and expected runtime: `performance/scales.yaml` and each experiment’s `README.md`.

---

## 1. Correctness (`tiny_even`)

**Goal:** Federated QC, KING, and LR agree with a centralized PLINK baseline.

**Docs:** [correctness/tiny_even/README.md](correctness/tiny_even/README.md)

```bash
# Data (if missing)
python pipeline/simulation/simulated_data/generate_synthetic_data.py \
  --scale tiny --partition-strategy even --seed 42 \
  --output-dir experiments/correctness/tiny_even/data

# Baseline
python experiments/tools/generate_baseline.py \
  experiments/correctness/tiny_even/config.yaml

# Simulation (default pyproject config)
flwr run . local-simulation --stream --run-config \
  'simulation=true num-server-rounds=100 config_path="experiments/correctness/tiny_even/configs"'

# Evaluate
python experiments/tools/evaluation/evaluate_all.py \
  experiments/correctness/tiny_even/results_2 \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline \
  --king
```

**Results path:** shipped Flower configs under `configs/` write to `experiments/correctness/tiny_even/results_2/` (cluster default). Local simulation or older layouts may use `results/` instead — check each center `log_dir` in `configs/*/config.yaml`.

**Success criteria:** QC F1 > 0.95; LR p-value correlation > 0.99 (see tiny_even README).

---

## 2. Performance & scalability 

**Goal:** Wall-clock time and peak memory on a **three-node** federated deployment (Matpool or equivalent), not single-machine simulation.

**Docs:**

- [performance/small_even/README.md](performance/small_even/README.md)
- [performance/medium_even/README.md](performance/medium_even/README.md)
- Cluster workflow: [../cluster_deployment/docs/CLUSTER_USER_GUIDE.md](../cluster_deployment/docs/CLUSTER_USER_GUIDE.md)

```bash
# Example: small scale data
python pipeline/simulation/simulated_data/generate_synthetic_data.py \
  --scale small --partition-strategy even --seed 42 \
  --output-dir experiments/performance/small_even/data

# After cluster run — rsync client logs to server, then:
python experiments/tools/collect_run_metrics.py \
  experiments/performance/small_even/results \
  --label matpool_small_even

# Optional: markdown/LaTeX row from exported stage CSVs
python experiments/tools/summarize_scalability_table.py \
  --results-root results --scale small --scale medium
```

**Results paths:**

| Scale | Directory |
|-------|-----------|
| tiny | `correctness/tiny_even/results_2/` |
| small | `performance/small_even/results/` |
| medium | `performance/medium_even/results/` |

Experiment `config.yaml` files enable `monitoring` and `retention: standard` by default. Sync the same release tag to all cluster nodes before benchmarking.

---

## 3. Real-world (`1000genomes`)

**Goal:** Application on real genotypes (1000 Genomes chr22) with generated phenotypes.

**Docs:** [real_world/README.md](real_world/README.md), [real_world/1000genomes/README.md](real_world/1000genomes/README.md)

Manuscript figures and tables live under `real_world/1000genomes/manuscript/`. Not required for the default release verification path in [RELEASE.md](../RELEASE.md).

---

## Shared configuration pattern

Each shipped experiment has:

```
<experiment>/
├── config.yaml           # experiment metadata, monitoring, retention, server rounds
├── configs/
│   ├── server/config.yaml
│   ├── center_1/config.yaml
│   └── center_2/config.yaml
├── data/                 # gitignored PLINK inputs
└── results/ or results_2/  # gitignored run outputs
```

Flower `config_path` points at `configs/` (directory), not the top-level `config.yaml`. Monitoring flags are read from the experiment `config.yaml` via `config_path`.

Many experiments also declare an `analysis:` block (for example in `correctness/tiny_even/config.yaml`):

```yaml
analysis:
  generate_baseline: true
  compare_results: true
  metrics_collection: true
```

These three flags are not read by the Flower pipeline (`client_app.py` / `server_app.py`). They document which post-run steps you intend to run manually. The table below maps each flag to the tool or runtime feature to use.

| `analysis` flag | What it means | How to run it |
|-----------------|---------------|---------------|
| `generate_baseline` | Centralized PLINK reference (QC, KING, LR) | `python experiments/tools/generate_baseline.py <experiment>/config.yaml` |
| `compare_results` | Federated vs baseline agreement | `python experiments/tools/evaluation/evaluate_all.py <results_dir> --baseline <baseline_dir>` (see [tools/evaluation/README.md](tools/evaluation/README.md)) |
| `metrics_collection` | Wall time / RSS / stage timing summary | **During run:** `monitoring` in `config.yaml` → `stage_metrics.csv` / `client_metrics.csv`. **After run:** `python experiments/tools/collect_run_metrics.py <results_dir>` |

Do not confuse `analysis` with **`retention`** (`retention.tier` in the same `config.yaml`). Retention runs automatically at server `done` when `retention.auto_apply_on_complete` is true; it prunes logs and intermediates but can keep evaluation reports (`evaluation_report.md`, `qc_report.md`, `lr_report.md`) per tier. It does not generate baselines or run evaluators.

---

## Tools

| Tool | Purpose |
|------|---------|
| [tools/evaluation/evaluate_all.py](tools/evaluation/evaluate_all.py) | QC + LR (+ optional KING) vs baseline |
| [tools/evaluation/metrics_collector.py](tools/evaluation/metrics_collector.py) | Aggregate logs and PLINK outputs |
| [tools/collect_run_metrics.py](tools/collect_run_metrics.py) | `run_summary.json` / `.md` from stage CSVs and GNU `time` |
| [tools/apply_run_retention.py](tools/apply_run_retention.py) | Apply retention tier manually (`--dry-run`) |
| [tools/generate_baseline.py](tools/generate_baseline.py) | Centralized reference run |
| [tools/summarize_scalability_table.py](tools/summarize_scalability_table.py) | Table rows from `results/<scale>_even/` stage exports |

See [tools/evaluation/README.md](tools/evaluation/README.md) for stage-specific evaluators.

---

## Results layout (after a run)

```
results/
├── server/
│   ├── logs/           # server_log.txt, stage_metrics.csv (when monitoring on)
│   └── intermediate/
├── center_1/
│   ├── logs/
│   └── intermediate/
├── center_2/
│   ├── logs/
│   └── intermediate/
├── stage_metrics.csv   # merged at server `done` (optional)
├── run_summary.md      # from collect_run_metrics.py
└── retention_manifest.json
```

Re-running in the same `log_dir` **appends** client `stage_metrics.csv` rows; archive or clear metrics before a fresh benchmark (see [cluster_deployment/README.md](../cluster_deployment/README.md)).