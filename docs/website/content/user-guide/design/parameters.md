---
title: Parameters
slug: /design/parameters
---

# Parameters

FedGWAS is configured through YAML files at three levels: experiment (`config.yaml`), server (`configs/server.yaml`), and each client (`configs/center_*.yaml`). This page describes the parameters that control chunking, federated QC, relatedness screening, association screening, and post-run output retention.

Client YAML files keep the `center_*` filename convention; each file configures one participating client.

For project layout and setup commands, see [Project and Data Setup](/user-guide/local-simulation/config-and-data).

## Configuration layers

| Layer | Typical file | Controls |
| --- | --- | --- |
| Experiment | `config.yaml` | Study metadata, client config paths, shared thresholds, server chunk / round defaults, retention |
| Server | `configs/server.yaml` | Server output directories used for intermediate files and logs |
| Center | `configs/center_*.yaml` | Local data path, per-center parameters, and thresholds |

Repository examples under `experiments/` may use a nested `configs/center_*/config.yaml` layout. The apps accept both layouts; do not mix path styles without updating surrounding files.

## Center config sections

### `input_data`

| Field | Type | Meaning |
| --- | --- | --- |
| `path` | string | PLINK 1.9 binary prefix (no `.bed` / `.bim` / `.fam` suffix). For VCF, the VCF path; conversion to PLINK format is applied where configured. |
| `type` | string | Input type: `bed` (PLINK 1.9 binary) or `vcf` (converted to PLINK-compatible representations) |

Example:

```yaml
input_data:
  path: data/center_1/study_center_1
  type: bed
```

### `output`

| Field | Type | Meaning |
| --- | --- | --- |
| `intermediate_dir` | string | Chunk files, maps, and other intermediate artifacts |
| `log_dir` | string | PLINK outputs, logs, and stage result text files |

### `parameters`

These fields control anonymization and chunking behavior.

| Field | Typical default | Meaning |
| --- | --- | --- |
| `sample_offset` | `1000000000000` | Offset used when anonymizing numeric sample FID/IID values |
| `chunk_size` | preset-dependent | Default partition size when a more specific chunk size is not set |
| `sample_chunk_size` | preset-dependent | Preferred sample partition size for KING chunks |
| `snp_chunk_size` | preset-dependent | Preferred SNP partition size for association chunks |
| `lr_target_chunks` | optional | Target number of LR chunks; used when deriving LR chunk size |
| `lr_sample_chunk_size` | optional | Explicit sample chunk size override for LR partitioning |
| `run_id` | optional string | Label written into run artifacts for easier result tracking |

Chunk sizing priority for iterative LR is roughly:

1. `lr_sample_chunk_size`
2. `sample_chunk_size`
3. server / config `chunk_size`
4. derivation from `lr_target_chunks` when configured

Larger chunks reduce round count but increase per-round payload size and server merge cost.

### `thresholds`

| Field | Typical default | Used by | Meaning |
| --- | --- | --- | --- |
| `maf_threshold` | `0.01` | Federated QC | Drop SNPs whose minor allele frequency is below this value |
| `missing_threshold` | `0.01`–`0.05` | Federated QC | Sample missingness (`--mind`) and SNP missingness filters |
| `hwe_threshold` | `1e-6` | Federated QC | Drop SNPs whose Hardy–Weinberg equilibrium p-value is below this value |
| `p_threshold` | preset-dependent | Association screening (local filter fallback) | Fallback insignificant-SNP cutoff when `local_lr_threshold` is unset |
| `local_lr_threshold` | preset-dependent | Association screening (local filter) | Local logistic-regression p-value cutoff for privacy-preserving SNP tokens |
| `global_lr_threshold` | preset-dependent | Association screening (federated test) | Cutoff used when classifying final federated association p-values |
| `king_threshold` | `0.23` | Relatedness screening | Relatedness cutoff for optional sample filtering after kinship estimation |

Thresholds may also appear under `config.yaml` → `clients.thresholds`. Keep experiment-level and center-level values aligned unless you intentionally override a center.

## Experiment and server parameters

### Experiment `config.yaml`

| Field | Meaning |
| --- | --- |
| `experiment_name` | Study label used in outputs and evaluation |
| `experiment_category` | Grouping label such as `correctness` or `performance` |
| `num_clients` | Expected number of centers; must match encrypted relay recipient counts |
| `scenario` | Scenario name used by evaluation / bookkeeping |
| `clients.config_files` | Map from client index to center YAML path |
| `clients.thresholds` | Shared QC / LR / KING thresholds |
| `server.chunk_size` | Default chunk size advertised by the strategy |
| `server.min_available_clients` | Flower minimum available clients |
| `server.min_fit_clients` | Flower minimum fit clients |
| `server.num_server_rounds` | Upper bound on Flower server rounds |
| `retention` | Post-run output retention policy (see below) |

`num_clients` is especially important for key exchange and encrypted client-to-client relay: each client encrypts messages for the expected peer set.

### `retention`

Experiment `config.yaml` can set a product retention tier that controls which run artifacts are kept after a successful federated run. This is the implementation of the paper’s **production- versus research-oriented** output modes: production-style runs minimize saved intermediate information; research-style runs retain richer diagnostic outputs. Retention is applied when the stage machine reaches `done` (if auto-apply is enabled). It is not a GWAS screening stage.

```yaml
retention:
  tier: standard          # minimal | standard | research
  auto_apply_on_complete: true
```

| Field | Default | Meaning |
| --- | --- | --- |
| `tier` | `standard` | Preset that selects which artifact classes to keep or prune |
| `auto_apply_on_complete` | `true` | When true, prune outputs automatically at successful completion |

Tier presets:

| Tier | Typical use | Intermediate / debug artifacts | Logs and node metrics | Science outputs / evaluation reports |
| --- | --- | --- | --- | --- |
| `minimal` | Production-style / lean deployment | Pruned | Mostly pruned | Kept |
| `standard` | Default production-style deployment | Pruned | Text logs and node metric CSVs kept | Kept |
| `research` | Research-oriented debugging / validation | Kept (including crypto state and KING accumulators) | Kept | Kept |

Individual `keep_*` flags from the preset can be overridden in the same `retention` block when needed. Manual preview or re-application is available via `experiments/tools/apply_run_retention.py --dry-run`. See [Outputs](/api-reference/outputs) for the resulting layout and `retention_manifest.json`.

### Server `configs/server.yaml`

```yaml
output:
  intermediate_dir: results/server/intermediate
  log_dir: results/server/logs
```

The server derives its working result directory from `output.log_dir` or `output.intermediate_dir`. KING and LR temporary merge directories are created under that server output tree.

## Parameter effects by component

| Component | Most important parameters |
| --- | --- |
| Federated quality control | `missing_threshold`, `maf_threshold`, `hwe_threshold` |
| Relatedness screening | `sample_chunk_size` / `chunk_size`, `sample_offset`, `king_threshold` |
| Association screening | `snp_chunk_size`, `lr_target_chunks`, `local_lr_threshold`, `global_lr_threshold`, `p_threshold` |
| Privacy / identity | `sample_offset`, `num_clients` |
| Output retention | `retention.tier`, `retention.auto_apply_on_complete` |

## Practical guidance

- Start from a generated preset (`fedgwas-sim setup-experiment ...`) and change one group of thresholds at a time.
- Keep client thresholds identical across clients for federated QC and local association filtering so exclusion logic stays consistent.
- Increase `num_server_rounds` when using smaller chunks or more centers; iterative KING/LR need enough rounds to drain all chunks.
- Use smaller chunk sizes for debugging correctness; use larger chunks when measuring performance on bigger datasets.
- Use `retention.tier: research` while validating a new experiment (research-oriented outputs); switch to `standard` or `minimal` for production-style runs that should discard large intermediate trees.

Related pages:

- [GWAS Components](/user-guide/design/gwas-components)
- [Workflow](/user-guide/design/workflow)
- [Privacy and Masking](/user-guide/design/privacy-masking)
- [Outputs](/api-reference/outputs)
