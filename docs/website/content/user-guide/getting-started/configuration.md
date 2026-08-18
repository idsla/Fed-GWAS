---
title: Configuration
slug: /configuration
---

# Configuration

Each participating client has its own `config.yaml` (files use the `center_*` naming convention). The template is stored at `configs/config_template.yaml`. The default run reads from `experiments/correctness/tiny_even/configs`, which contains `server/config.yaml`, `center_1/config.yaml`, and `center_2/config.yaml`.

## Center configuration

```yaml
input_data:
  path: "experiments/correctness/tiny_even/data/tiny/center_1/tiny_center_1"
  type: "bed"

output:
  intermediate_dir: "experiments/correctness/tiny_even/results_2/center_1/intermediate"
  log_dir: "experiments/correctness/tiny_even/results_2/center_1/logs"

parameters:
  sample_offset: 1000000000000
  chunk_size: 100
  sample_chunk_size: 100
  snp_chunk_size: 100
  lr_target_chunks: 2
  run_id: "20260221_2000"

thresholds:
  maf_threshold: 0.01
  missing_threshold: 0.01
  hwe_threshold: 1e-6
  p_threshold: 0.3
  local_lr_threshold: 0.3
  global_lr_threshold: 0.1
  king_threshold: 0.23

flower:
  server_address: "127.0.0.1:8080"
  num_rounds: 10

participation:
  key_exchange: true
  sync: true
  local_qc: true
  global_qc: true
  global_qc_response: true
  init_chunks: true
  iterative_king: true
  local_lr: true
  local_lr_filter_response: true
  init_chunks_lr: true
  iterative_lr: true
```

## Field reference

| Field | Purpose |
| --- | --- |
| `input_data.path` | PLINK 1.9 dataset prefix. Do not include `.bed`, `.bim`, or `.fam`. For VCF input, use the VCF file path; conversion to PLINK-compatible representations is applied where configured. |
| `input_data.type` | Input type. The loader currently handles `bed` (PLINK 1.9 binary) and `vcf` (via conversion). |
| `output.intermediate_dir` | Per-center directory for chunks and temporary files. |
| `output.log_dir` | Per-center directory for logs and PLINK outputs. |
| `parameters.chunk_size` | Default chunk size. KING commonly partitions by samples; final LR commonly partitions by SNPs. |
| `parameters.lr_target_chunks` | Optional target number of LR chunks; used by `compute_lr_chunk_size`. |
| `parameters.run_id` | Optional run tag used to scope intermediate files under `run_<id>`. |
| `thresholds.maf_threshold` | Minimum minor allele frequency used by federated QC. |
| `thresholds.missing_threshold` | Maximum SNP missingness allowed by federated QC. |
| `thresholds.hwe_threshold` | Hardy–Weinberg equilibrium p-value threshold. |
| `thresholds.p_threshold` | Backward-compatible association-screening threshold. |
| `thresholds.local_lr_threshold` | Local logistic-regression threshold for selecting insignificant SNPs to filter before federated association screening. |
| `thresholds.global_lr_threshold` | Threshold for classifying final association-screening p-values on clients. |
| `thresholds.king_threshold` | Kinship threshold used by client-side relatedness screening. |
| `participation` | Per-stage opt-in flags. Use `true` for all stages when running the full pipeline. |

## Server configuration

`pipeline/server_app.py` requires a server config file under the configured experiment path:

```text
<config_path>/server/config.yaml
```

The server config declares the server output directories:

```yaml
output:
  intermediate_dir: "experiments/correctness/tiny_even/results_2/server/intermediate"
  log_dir: "experiments/correctness/tiny_even/results_2/server/logs"
```

## Flower application configuration

Flower entry points and default run settings live in `pyproject.toml`:

```toml
[tool.flwr.app.components]
serverapp = "pipeline.server_app:app"
clientapp = "pipeline.client_app:app"

[tool.flwr.app.config]
simulation = true
num-server-rounds = 300
default-partition-by = "samples"
config_path = "experiments/correctness/tiny_even/configs"
```

In simulation mode, `pipeline/client_app.py` derives the center config path from the partition id:

```text
<config_path>/center_<partition_id + 1>/config.yaml
```
