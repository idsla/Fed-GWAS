# Current Version Documentation (Feb 2026, updated 2026-05-23)

Release **v0.0.1** — initial public release.

If you are looking for the quickstart, see `README.md`. For release verification, see [RELEASE.md](RELEASE.md). For three-node deployment, see [cluster_deployment/docs/CLUSTER_USER_GUIDE.md](cluster_deployment/docs/CLUSTER_USER_GUIDE.md).

---

## Manual PyPI publish

The PyPI distribution name is `FedGWAS` (`pip install FedGWAS`). PyPI normalizes package names, so generated artifact filenames use lowercase `fedgwas-...`.

Manual publishing should only happen after the release verification path in [RELEASE.md](RELEASE.md) passes.

```bash
# Optional: verify no stale build artifacts remain.
rm -rf dist/

# Build source and wheel distributions.
uv build

# Validate package metadata and README rendering.
uvx twine check dist/*

# Optional first publish target for a dry run.
uvx twine upload --repository testpypi dist/*

# Install-check from TestPyPI while resolving dependencies from PyPI.
python -m pip install \
  --index-url https://test.pypi.org/simple/ \
  --extra-index-url https://pypi.org/simple \
  FedGWAS==0.0.1

# Publish to production PyPI.
uvx twine upload dist/*
```

For manual upload, use a PyPI API token with username `__token__`. PyPI versions are immutable: after `0.0.1` is uploaded, any fix requires a new version such as `0.0.2`.

The preferred release path is still GitHub Actions Trusted Publishing, configured in `.github/workflows/publish-pypi.yml`, because it does not require storing a PyPI token in GitHub.

---

## v0.2.0 highlights (2026-05-23)

### Cluster deployment (`cluster_deployment/`)
- **User guide**: `cluster_deployment/docs/CLUSTER_USER_GUIDE.md` (setup, clock skew, faketime, teardown, metrics rsync).
- **Scripts shipped** (11): setup, start/stop/status, run-app, clock-skew check, data verify, diagnose, offline Flower wheels, config template.
- **Removed from bundle**: Flower site-packages debug patches (`patch-flwr-*`) — not needed for production runs.
- **`cluster-verify-data.sh`**: `--client-id 1|2|all` so each client checks only its PLINK triplet.
- Wall/RSS logging: wrap `cluster-run-app.sh` with `/usr/bin/time -v -o <file>` (see `cluster_deployment/README.md`).

### Performance experiments (appendix)
- **`experiments/performance/`**: `small_even`, `medium_even` configs + per-scale READMEs; `scales.yaml`; `summarize_scalability_table.py`. Manuscript drafts (`appendix_*.md`, `manuscript_sections.md`, top-level `README.md`) are **local-only** (`.gitignore`).
- **`docker/`**: local Docker workflow; **not** in the git tree (use `cluster_deployment/` on clusters).
- Enable **`monitoring`** + **`retention: standard`** in experiment `config.yaml`; sync latest code to all nodes before cluster runs (metrics CSVs require v0.2.0 pipeline on every node).

### Monitoring & retention (unchanged behavior, documented)
- Per-node `stage_metrics.csv` / `client_metrics.csv` under each `logs/`; merged at server `done`.
- `collect_run_metrics.py` + optional `/usr/bin/time -v` file for paper-ready `run_summary.md`.

---
## Pipeline flow (current version)

The server strategy is **stage-based** and is implemented in `pipeline/server/strategy_strict.py`.

### Stage 1: `key_exchange`
- **Clients**: generate ECC keypair, send public key.
- **Server**: collects public keys.

**Input/Output Formats**:
- **Client → Server**: `Parameters` containing ECC public key as `numpy.ndarray` (uint8, P-256 curve, 65 bytes: 1-byte prefix + 64-byte key).
- **Server → Client**: Config with curve parameters (curve name, generator point, etc.) via `ConfigRecord`.
- **Output**: Server stores all public keys; clients store their own keypair locally.

### Stage 2: `sync` (encrypted seed relay)
- **Clients**: encrypt-broadcast their `local_seed` to all other clients using ECDH-derived AES keys.
- **Server**: forwards encrypted blobs via `Parameters` (**server cannot decrypt**).

**Input/Output Formats**:
- **Client → Server**: `Parameters` containing encrypted binary envelope (AES-256-CBC ciphertext) as `numpy.ndarray` (uint8).
  - Envelope format: `[recipient_hash_id(32B)][iv(16B)][ciphertext][hmac(32B)]`.
  - Each client sends one encrypted message per recipient (N-1 messages for N clients).
- **Server → Client**: `Parameters` containing all encrypted envelopes (server cannot decrypt).
- **Output**: Clients receive encrypted seed messages; server only sees ciphertexts.

### Stage 3: `sync_response`
- **Clients**: decrypt incoming encrypted seed messages addressed to them; compute `global_seed`.
- **Fallback**: if decrypting is impossible, compute a deterministic seed from public information (no security).

**Input/Output Formats**:
- **Input**: Encrypted seed messages from Stage 2.
- **Processing**: Clients decrypt using ECDH-derived AES keys, extract `local_seed` (integer) from each message.
- **Output**: `global_seed` computed as XOR of all decrypted `local_seed` values (or deterministic fallback).
- **Storage**: `global_seed` persisted locally for use in subsequent stages.

### Stage 4: `local_qc`
- **Clients**: run sample missingness filtering locally (PLINK `--mind`). (Note: currently only implement this, can be extended straightforwardly.)

**Input/Output Formats**:
- **Input**: Local PLINK files (`.bed/.bim/.fam`) from client's input data.
- **Processing**: PLINK command `--mind <threshold>` filters samples by missingness rate.
- **Output**: Filtered PLINK files (same format) written to intermediate directory.
- **Network**: No network communication (local-only stage).

### Stage 5: `global_qc` (encrypted QC relay)
- **Clients**: compute QC arrays (genotype counts + missingness + thresholds), pack typed payloads, encrypt-broadcast them to all clients.
- **Server**: forwards encrypted blobs (**server cannot decrypt**).

**Input/Output Formats**:
- **Client → Server**: `Parameters` containing encrypted typed payloads as `numpy.ndarray` (uint8).
  - Payload format: `"FGQP" + kind_len(u8) + kind_ascii + encoded_ndarray_bytes`.
  - Three payload types: `"counts"` (genotype counts per SNP), `"missing"` (missingness rates), `"thresholds"` (QC thresholds).
  - Each payload encrypted with AES-256-CBC using ECDH-derived keys (same envelope format as Stage 2).
- **Server → Client**: `Parameters` containing all encrypted QC payloads (server cannot decrypt).
- **Output**: Clients receive encrypted QC shares from all clients.

### Stage 6: `global_qc_response`
- **Clients**: decrypt QC blobs, aggregate them locally, compute SNP exclusion list, apply SNP filtering.

**Input/Output Formats**:
- **Input**: Encrypted QC payloads from Stage 5.
- **Processing**: 
  - Decrypt and unpack typed payloads (`counts`, `missing`, `thresholds`).
  - Aggregate arrays across all clients (sum counts, aggregate missingness).
  - Compute SNP exclusion list: exclude SNPs failing MAF, missingness, or HWE thresholds.
- **Output**: SNP exclusion list (list of SNP IDs) applied to local dataset via PLINK `--exclude`.
- **Storage**: Filtered PLINK files written to intermediate directory; exclusion list saved for reference.

### Stage 7: `init_chunks`
- **Clients**: chunk the dataset deterministically based on `global_seed`.
- **Anonymization**: clients anonymize sample IDs (FID/IID) and SNP IDs in each chunk:
  - Sample IDs: anonymized using `sample_offset`, `chunk_index`, and `iteration_id` to prevent cross-chunk/iteration tracking.
  - SNP IDs: anonymized using `global_seed`, `chunk_index`, and `iteration_id` via deterministic hashing.
  - ID maps (anonymized → original) are stored locally for later de-anonymization of server results.

**Input/Output Formats**:
- **Input**: Filtered PLINK files from Stage 6, `global_seed` from Stage 3, and chunk controls from stage/client config.
- **Processing**:
  - Deterministic chunking for KING (current implementation uses SNP-based chunking).
  - Chunk size defaults to strategy value (`chunk_size=1000`) and can be overridden by client-side parameters (for example `snp_chunk_size`).
  - For each chunk: extract subset via PLINK `--keep`, then anonymize IDs.
  - Anonymization: create new `.bed/.bim/.fam` files with anonymized IDs.
- **Output**: 
  - Anonymized chunk files (`.bed/.bim/.fam`) stored in intermediate directory.
  - ID maps: `*_sample_map.tsv` (anon_sample → orig_sample) and `*_snp_map.tsv` (anon_snp → orig_snp).
- **Storage**: Chunks and maps persisted locally; ready for transmission in Stage 8.

### Stage 8: `iterative_king`
- **Clients**: send one chunk (as packed array) per round.
- **Clients (privacy add-on)**: encrypt-broadcast a per-round **anon→stable pseudonym** map to each peer (server cannot decrypt).
- **Server**: reconstructs `.bed/.bim/.fam`, merges, runs:
  - PLINK 1.9 `--het`
  - PLINK 2.0 `--make-king-table`

**Input/Output Formats**:
- **Client → Server**: `Parameters` containing packed binary chunk as the **first** `numpy.ndarray` (uint32 or uint8).
  - Format: `[bed_size(uint32)][bim_size(uint32)][fam_size(uint32)][bed_bytes][bim_bytes][fam_bytes]`.
  - Chunk contains anonymized `.bed/.bim/.fam` files serialized as raw bytes.
- **Client → Server (encrypted peer maps)**: additional `numpy.ndarray` payloads (uint8) after the chunk.
  - Each payload is an encrypted typed message (`kind="king_map"`) carrying:
    - `iter_id`, `chunk_index`
    - `map`: list of `[anon_id, stable_pseudonym]` for this chunk
  - **Server forwards these blobs but cannot decrypt them.**
- **Server Processing**:
  - Reconstructs `.bed/.bim/.fam` files from packed arrays.
  - Merges chunks from all clients using PLINK `--bmerge`.
  - Runs PLINK 1.9 `--het` (outputs `.het` file with heterozygosity stats).
  - Runs PLINK 2.0 `--make-king-table` (outputs `.kin0` file with kinship coefficients).
- **Output**: Server produces KING results with anonymized sample ID pairs.

### Stage 8b: KING response payload (within `iterative_king`)
- **Server**: sends KING kinship results back to clients **and forwards encrypted peer maps** during `iterative_king` rounds.
- **Clients**: receive and process KING results, decrypt peer maps, update accumulators, and finalize once all chunks are done.

**Input/Output Formats**:
- **Server → Client**: `Parameters` containing KING results as the **first** `numpy.ndarray` (uint8).
  - Format: UTF-8 encoded text, one line per sample pair: `"sampleA sampleB IBS0 HETHET NSNP KINSHIP\n"`.
  - Sample IDs are anonymized (from Stage 7).
- **Server → Client (forwarded peer maps)**: additional encrypted `numpy.ndarray` payloads (uint8), one per peer.
  - Decrypted by clients only; **server cannot link samples across iterations**.
- **Client Processing**:
  - Parse KING results, extract sample pairs and kinship coefficients.
  - Decrypt peer maps and resolve **anon IDs → stable pseudonyms** for cross-client pairs.
  - Use original IDs for same-client pairs; use stable pseudonyms for cross-client pairs.
  - Accumulate sufficient stats across chunks (NSNP-weighted rates): `sum_ibs0 += IBS0 * NSNP`, `sum_hethet += HETHET * NSNP`, `sum_nsnp += NSNP`.
  - Compute final kinship per pair using HETHET-weighted aggregation across chunks (linear-fit fallback remains for missing HETHET cases).
- **Output**: 
  - KING accumulator updated with de-anonymized local samples and anonymized remote samples.
  - After all chunks: filter related samples (e.g., kinship > 0.2), save results to `.tsv` files.

### Stage 9: `local_lr`
- **Clients**: run local LR; produce list of "insignificant SNPs".
- **Privacy**: clients send **tokens** (not plaintext SNP IDs) to server.

**Input/Output Formats**:
- **Input**: Local PLINK files (post-QC, post-KING filtering), phenotype data from `.fam`.
- **Processing**: Run PLINK 1.9 `--logistic` locally, parse results to identify insignificant SNPs (p-value > threshold).
- **Client → Server**: `Parameters` containing tokenized SNP identifiers as `numpy.ndarray` (uint8).
  - Format: UTF-8 encoded text, space-separated tokens: `"token1 token2 token3 ..."`.
  - Tokens generated via `hash(SNP_ID + salt)` to prevent server from learning SNP identities.
  - Optional padding to hide set sizes.
- **Output**: Server receives tokens only; cannot determine which SNPs are insignificant.

### Stage 10: `local_lr_filter_response`
- **Server**: forwards token sets to clients (router behavior).
- **Clients**: resolve tokens back to SNP IDs locally and filter.

**Input/Output Formats**:
- **Server Processing**: 
  - Receives token sets from all clients.
  - **Router behavior**: forwards all token sets to clients (server does NOT compute intersection).
  - Server only sees individual token sets, not the intersection result.
- **Server → Client**: `Parameters` containing all token sets as `numpy.ndarray` (uint8) list.
  - Format: Each array is UTF-8 encoded text, space-separated tokens (one set per client).
- **Client Processing**:
  - Receives all token sets from server.
  - **Computes intersection locally**: tokens present in all client sets (secure intersection).
  - Resolve intersection tokens back to SNP IDs using local token→SNP map.
  - Apply filtering: exclude SNPs whose tokens are in the intersection.
- **Output**: Filtered PLINK files (SNPs with significant p-values in at least one client retained).
- **Privacy**: Server never sees the intersection result; clients compute it locally for enhanced privacy.

### Stage 11: `init_chunks_lr`
- **Clients**: re-chunk post-filter dataset.

**Input/Output Formats**:
- **Input**: Filtered PLINK files from Stage 10, `global_seed` (same as Stage 7), and LR chunk controls.
- **Processing**: Deterministic LR chunking (sample-based in current implementation) plus anonymization of the post-filter dataset.
- **Output**: 
  - Anonymized LR chunks (`.bed/.bim/.fam`) with updated SNP sets.
  - ID maps: `*_sample_map.tsv` and `*_snp_map.tsv` (may differ from Stage 7 due to SNP filtering).

### Stage 12: `iterative_lr`
- **Clients**: send LR chunks.
- **Server**: reconstructs, merges, runs PLINK 1.9 `--logistic`.

**Input/Output Formats**:
- **Client → Server**: `Parameters` containing packed binary chunk (same format as Stage 8).
  - Format: `[bed_size(uint32)][bim_size(uint32)][fam_size(uint32)][bed_bytes][bim_bytes][fam_bytes]`.
  - Chunks contain anonymized data with filtered SNP sets.
- **Server Processing**:
  - Reconstructs `.bed/.bim/.fam` files from packed arrays.
  - Merges chunks from all clients.
  - Runs PLINK 1.9 `--logistic` (outputs `.assoc.logistic` file with p-values per SNP).
- **Output**: Server produces LR association results with anonymized SNP IDs.

### Stage 13: LR response payload (within `iterative_lr`)
- **Server**: sends LR p-values back to clients during `iterative_lr` rounds.
- **Clients**: receive LR p-values and de-anonymize results locally.

**Input/Output Formats**:
- **Server → Client**: `Parameters` containing LR results as `numpy.ndarray` (uint8).
  - Format: UTF-8 encoded text, one line per SNP: `"SNP_ID p-value\n"` (PLINK `.assoc.logistic` format).
  - SNP IDs are anonymized (from Stage 11).
- **Client Processing**:
  - Parse LR results, extract SNP IDs and p-values.
  - De-anonymize SNP IDs using local SNP maps (from Stage 11).
  - Store results: `{original_SNP_ID: p_value}` dictionary.
- **Output**: 
  - LR state is persisted for iterative accumulation.
  - Optional debug/screening artifacts can be emitted to:
    - `lr_results_client_deanon.txt`
    - `lr_screening_significant.txt`
    when `emit_lr_intermediate=1` (or `FEDGWAS_EMIT_LR_INTERMEDIATE=1`).

### End: `done`
- Strategy transitions to `done`; server does not configure further work.

---

## Repo structure (what matters for the current pipeline)

### Top-level
- **`README.md`**: Environment setup + how to run in simulation/deployment.
- **`pyproject.toml`**: Flower app entrypoints and default run config.
- **`experiments/`**: Experiment configs, datasets, and evaluation tooling.
- **`configs/`**: Config templates and example deployment configs.
- **`pipeline/`**: The Flower app implementation (server + client).
- **`plink/`**: PLINK 1.9 and PLINK 2.0 binaries (project-local).

### `pipeline/`
- **Entry points**
  - `pipeline/server_app.py`: Flower `ServerApp` entry.
  - `pipeline/client_app.py`: Flower `ClientApp` entry.
- **Server**
  - `pipeline/server/strategy_strict.py`: **Current strategy used** (stage orchestration + encrypted relay forwarding).
  - `pipeline/server/aggregator_king.py`: Rebuild/merge chunks + run KING (PLINK2).
  - `pipeline/server/aggregator_lr.py`: Rebuild/merge chunks + run LR (PLINK1.9).
  - `pipeline/server/prg_masking.py`: ECC key exchange helper (public keys).
  - Legacy strategies (`strategy.py`, `strategy_new.py`, `aggregator_qc.py`) are kept locally only (`.gitignore`).
- **Client**
  - `pipeline/clients/base_client.py`: Shared client base, chunk packing, phenotype-fix switch.
  - `pipeline/clients/data_loder.py`: Loads per-client YAML config (input/output/thresholds).
  - `pipeline/clients/client_to_client.py`: **ECDH + AES** encrypted message envelopes (server relay cannot decrypt).
  - `pipeline/clients/c2c_payloads.py`: Typed payload definitions for client-to-client messages.
  - `pipeline/clients/seed_sync.py`: Compute global seed from encrypted relay (or deterministic fallback).
  - `pipeline/clients/local_qc.py`: Local QC + local LR helper commands.
  - `pipeline/clients/client_qc_aggregator.py`: Compute SNP exclusion list client-side after decrypting QC shares.
  - `pipeline/clients/lr_privacy.py`: Tokenize "insignificant SNPs" for privacy-preserving LR filtering.
  - `pipeline/clients/iterative_king.py`, `pipeline/clients/iterative_lr.py`: Chunk loop logic.
  - `pipeline/clients/flwr_config.py`: Flower configuration helpers.
  - `pipeline/clients/logger_manager.py`: Centralized logging management.
- **Monitoring & retention** (`pipeline/utils/`)
  - `monitoring_config.py`, `retention_config.py`, `run_retention.py`
  - `performance/`: runtime monitors and CSV merge (`monitoring_runtime.py`, etc.)
- **Local-only (not in git)**
  - `pipeline/simulation/`: synthetic data generators and local experiment tooling.
  - `pipeline/visualization_archived/`: archived plots/post-analysis (reference).
- **Utilities**
  - `pipeline/utils/client_data_loader.py`: Client data loading utilities.

---

## Documentation site

The Docusaurus site is isolated under `docs/website/`. Site Markdown, images, sidebars, and Docusaurus configuration now live together under that directory.

Install frontend dependencies:

```bash
cd docs/website
npm install
```

Run the documentation site locally:

```bash
cd docs/website
npm run start
```

Build the static documentation site:

```bash
cd docs/website
npm run build
```

Preview the production build:

```bash
cd docs/website
npm run serve
```

With the current GitHub Pages configuration, the local served URL is `http://127.0.0.1:3000/Fed-GWAS/`.

---

## Run modes

### Simulation mode (Flower simulation engine)

Default federation: `local-simulation` (configured in `pyproject.toml`).

Command:

```bash
flwr run . local-simulation --stream
```

Key points:
- Default `config_path` is `experiments/correctness/tiny_even/configs` (from `pyproject.toml`).
- Override `config_path` at runtime for other experiments (for example `experiments/real_world/1000genomes/configs`).
- Clients are typically recreated by the simulator, so we persist seeds/keys opportunistically.

### Local deployment mode (SuperLink + SuperNodes)

Requires three processes:
- 1× `flower-superlink` (Fleet API on 9092; Exec API on 9093 by default)
- 2× `flower-supernode` (one per client partition)
- then `flwr run . local-deployment`

Recommended (ports must match your setup):

```bash
flower-superlink --insecure --fleet-api-address 127.0.0.1:9092
```

```bash
flower-supernode --insecure --superlink 127.0.0.1:9092 --clientappio-api-address 127.0.0.1:9094 \
  --node-config 'partition-id=0 num-partitions=2 config-file="experiments/correctness/tiny_even/configs/center_1/config.yaml"'
```

```bash
flower-supernode --insecure --superlink 127.0.0.1:9092 --clientappio-api-address 127.0.0.1:9095 \
  --node-config 'partition-id=1 num-partitions=2 config-file="experiments/correctness/tiny_even/configs/center_2/config.yaml"'
```

Finally:

```bash
flwr run . local-deployment --stream --run-config "simulation=false num-server-rounds=15"
```

Important deployment notes:
- The app connects to the SuperLink **Exec API** (configured as `address` in `[tool.flwr.federations.local-deployment]` in `pyproject.toml`; typically `127.0.0.1:9093`).
- The `--node-config ... config-file=...` value is now honored by the client (see the "Changes" section).
- PLINK binaries are resolved from the repository (`plink/plink_*`) first, then PATH as fallback.

---

## Experiments

The `experiments/` directory contains experiment scenarios and analysis tools for testing the federated GWAS pipeline. For detailed information about experiment categories, configuration, and usage, see `experiments/README.md`.

### Experiment Structure

```
experiments/
├── README.md            # Detailed experiment documentation
├── correctness/
│   └── tiny_even/       # Tiny scale, even partition (baseline validation)
└── tools/               # Experiment tooling
    ├── evaluation/              # Post-evaluation tooling
    │   ├── metrics_collector.py # Metrics collection framework
    │   ├── collect_real_world_results.py # Result collection for cluster
    │   ├── qc/                  # QC evaluation utilities
    │   ├── king/                # KING evaluation utilities
    │   └── lr/                  # LR evaluation utilities
    ├── generate_baseline.py     # Centralized baseline generator
    └── setup_real_world_experiment.py  # Real-world experiment setup
```

### Correctness Experiment: tiny_even

**Purpose**: Verify federated results match centralized PLINK

**Configuration**:
- 2 clients, even split
- Baseline validation (5-10 min runtime)

**Metrics**:
- QC SNP exclusion list agreement (F1 score)
- KING kinship coefficient correlation (Pearson r)
- LR p-value correlation (Pearson r, Manhattan plots)
- Top significant SNPs overlap


### Running Experiments

#### Using Flower CLI

Experiments are run via Flower CLI and the config sets under `experiments/`.

#### Experiment Configuration

Each experiment has a main `config.yaml` file that specifies:
- **Experiment metadata**: category, purpose, expected runtime
- **Client configurations**: paths to per-client config files
- **Data settings**: data directory, partition strategy, scale
- **Server settings**: rounds, chunk sizes, thresholds
- **Monitoring settings**: enable/disable metrics collection
- **Analysis settings**: baseline generation, result comparison flags

Example config structure:
```yaml
experiment_name: correctness_tiny_even
experiment_category: correctness
description: 'Correctness validation: Tiny scale with even partition'

clients:
  config_files:
    0: experiments/correctness/tiny_even/configs/center_1/config.yaml
    1: experiments/correctness/tiny_even/configs/center_2/config.yaml

data:
  data_dir: experiments/correctness/tiny_even/data/tiny
  partition_strategy: even
  scale: tiny

server:
  num_server_rounds: 15
  chunk_size: 100
  min_available_clients: 1
  min_fit_clients: 1

monitoring:
  enable_network_monitoring: true
  enable_performance_monitoring: true
```

When `enable_performance_monitoring` is true (read from experiment `config.yaml` via `config_path`), the pipeline writes per-node CSVs under each `logs/` directory (`stage_metrics.csv`, `client_metrics.csv`, `performance_summary.md`) and merges them into the experiment results root when the server reaches the `done` stage. Implementation lives in `pipeline/utils/performance/` (`monitoring_runtime.py`, `performance_monitoring.py`, `network_monitor.py`).

Memory-related monitoring defaults: `resource_monitoring_interval` 2s on clients; `network_monitoring_interval` 5s with at most 500 retained samples per client. Disable `enable_network_monitoring` if only CPU/RAM stage CSVs are needed.

#### Analysis Tools

**Metrics Collector** (`experiments/tools/evaluation/metrics_collector.py`):
- Collects and aggregates metrics from experiment runs
- Per-stage timings from `stage_metrics.csv` when present (PerformanceMonitor), else log parsing
- Resource usage (memory, CPU, I/O)
- Network traffic statistics
- PLINK output parsing (QC, KING, LR results)
- Export to structured format (JSON/CSV)

**QC Evaluator** (`experiments/tools/evaluation/qc/qc_evaluator.py`):
- Compares federated QC exclusion lists with centralized baselines
- QC SNP exclusion list agreement (F1 score)
- Summary reporting

**LR Evaluator** (`experiments/tools/evaluation/lr/lr_evaluator.py`):
- Compares federated LR results with centralized baselines
- Local/global p-value correlation
- Significance agreement and coverage
- Generates Manhattan and correlation plots

**Baseline Generator** (`experiments/tools/generate_baseline.py`):
- Creates centralized PLINK baselines for comparison
- Merges federated datasets from all centers
- Runs centralized PLINK QC, KING, and LR
- Stores results for comparison

#### Results Organization

Each experiment's `results/` directory contains:
- **Per-client results**: `center_1/`, `center_2/`, etc.
  - `intermediate/` - GWAS intermediate files (chunks, QC results)
  - `logs/` - Client-specific logs and PLINK outputs
  - `performance/` - Client performance metrics
- **Server results**: `server/`
  - `logs/` - Server logs and federation management
  - `performance/` - Server performance metrics
  - `intermediate/` - Aggregated results (KING, LR)
- **Analysis results**:
  - `collected_metrics.json` - Aggregated metrics
  - `experiment_report.md` - Result comparison report
  - `manhattan_plot.png` - Manhattan plot visualization
  - `correlation_plots.png` - Correlation visualizations

---

## Output layout (current behavior)

### Simulation
- Uses whatever `output.intermediate_dir` / `output.log_dir` are set in the scenario config YAML.

### Deployment
To prevent collisions between clients that share a common base output directory, deployment mode now forces:

- `intermediate_dir = <configured_intermediate_dir>/<client_id>`
- `log_dir = <configured_log_dir>/<client_id>`

This is idempotent (will not keep nesting on repeated state reload).

## Privacy/security model (current implementation)

### What is protected from the server (honest relay assumption)
- **Global seed**: computed client-side from decrypted seed messages; server only forwards ciphertexts.
- **Global QC arrays**: server only forwards ciphertexts; aggregation is done client-side.
- **KING cross-iteration linkability**: server only sees per-round anonymized IDs; stable IDs are exchanged client-to-client under encryption.
- **LR “insignificant SNPs”**: server sees only **tokens**, not SNP IDs.

### Cryptography (current)
- **Key exchange**: ECC keys (P-256).
- **Transport encryption**: ECDH-derived AES-256-CBC for client-to-client messages.
- **Envelope format**: custom binary envelope (no pickle) in `pipeline/clients/client_to_client.py`.

### Known limitations (still true)
- If the server can tamper with public keys, it can mount a MITM by substituting keys (no PKI/pinning by default).
- Tokenization hides SNP IDs but may still leak sizes and overlaps unless padded (optional).

---

## Phenotype encoding switch (`-9/0 -> 2`)

PLINK expects: `1=control`, `2=case`, `-9/0=missing`.

For toy/simulated datasets, missing phenotypes can lead to “phenotype is constant” warnings.

Current behavior:
- In simulation mode (`simulation=true`), phenotype fix can be applied.
- In deployment mode, phenotype fix is **opt-in** via run-config:

```bash
--run-config "phenotype_fix_missing_to_case=true"
```

Implementation detail:
- We do **not mutate** the original input `.fam`.
- We create a derived prefix `phenofixed_<client_id>` in the client intermediate directory and point the pipeline at that prefix.

---

## Major changes made during stabilization (summary)

### Deployment correctness
- Client now honors `--node-config ... config-file=...` (deployment mode).
- Deployment needs SuperLink **Exec API address** (`pyproject.toml` federation config).
- Deployment outputs are now per-client subdirectories to avoid collisions.

### Stage sync + encryption
- Implemented encrypted client-to-client relay for:
  - `sync` seed exchange (AES-256-CBC with ECDH key derivation)
  - `global_qc` QC share exchange (typed payloads with encryption)
- Removed `pickle` network payload fallback (security improvement).
- Added deterministic fallback seed computation if encrypted sync cannot be completed.
- Custom binary envelope format: `[recipient_hash_id][iv][ciphertext][hmac]`.

### Server-side chunk reconstruction fixes
- Fixed server aggregators to correctly extract ndarrays from Flower `Parameters` in both simulation and deployment modes.
- Fixed dtype conversion: clients may send uint32 arrays where values represent bytes; server now reconstructs correctly.
- Robust handling of packed binary chunks: `[bed_size][bim_size][fam_size][bed_bytes][bim_bytes][fam_bytes]`.

### PLINK compatibility fixes
- KING now uses PLINK2 (`--make-king-table`) instead of PLINK 1.9 `--genome`.
- Merge fixed to use `--bmerge <prefix>` (not `.bed .bim .fam` triplet).
- Proper handling of heterozygosity computation (`--het`) before KING analysis.

### Privacy enhancements
- **Anonymization**: Sample and SNP IDs anonymized with chunk/iteration-specific salts to prevent cross-chunk tracking.
- **Tokenization**: Local LR sends **tokenized** insignificant SNP identifiers instead of plaintext SNP IDs.
- **De-anonymization**: Clients can only de-anonymize their own samples/SNPs; remote data remains anonymized.

### Response stages
- Added KING/LR response payload handling so server sends results back to clients during iterative rounds for de-anonymization and accumulation.

### Runtime cleanup and output hygiene (2026-02)
- Removed direct `print(...)` noise from active client/server pipeline paths; runtime output is now logger-driven.
- Reduced chunk payload verbosity by default; detailed payload diagnostics are gated behind debug flags:
  - `debug_chunk_payloads=1` (or `FEDGWAS_DEBUG_CHUNK_PAYLOADS=1`)
  - `debug_king_pairs=1` (or `FEDGWAS_DEBUG_KING_PAIRS=1`)
- Removed redundant LR intermediate persistence paths (`lr_pvals.pkl`, `lr_significance.pkl`).
- LR de-anonymized/screening files are now opt-in (`emit_lr_intermediate=1` or `FEDGWAS_EMIT_LR_INTERMEDIATE=1`).

### Tooling and experiments
- Analysis tools: metrics collector, result analyzer, baseline generator.

---

## Troubleshooting (most common)

- **“No such file or directory: plink” in deployment**:
  - Ensure the process runs from the project root so bundled `plink/plink_*` binaries are discoverable, or set `PATH` to include a PLINK binary.
- **Server log shows “No bed files received”**:
  - Indicates chunk reconstruction decoding mismatch; confirm you are using current `pipeline/server/aggregator_king.py` / `aggregator_lr.py`.
- **LR warns “phenotype is constant”**:
  - Enable `phenotype_fix_missing_to_case=true` for toy runs, or provide real phenotypes with cases/controls.
