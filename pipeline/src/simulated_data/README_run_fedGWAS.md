# Federated GWAS Pipeline (Updated)

This repository implements a federated pipeline for Genome-Wide Association Studies (GWAS) using Flower, PLINK, and custom privacy-preserving protocols. The pipeline supports multi-stage, multi-client GWAS with robust output management, reproducibility, and clear logging.

---

## Environment Setup

### 1. Create and Activate a Conda Environment
We recommend [Miniconda](https://docs.conda.io/en/latest/miniconda/) or [Anaconda](https://www.anaconda.com/products/distribution) to manage your Python environment.

```bash
conda create -n fedgwas python=3.11 -y
conda activate fedgwas
```

Install Poetry
```bash
# Install peotry - in git bash
pipx install poetry
pipx ensurepath
```

#### Key Fields
- **input_data.path**: Path to the input genotype data (PLINK prefix or VCF file).
- **output.intermediate_dir**: Directory for intermediate files (per client, auto-cleared each run).
- **output.log_dir**: Directory for logs (per client, auto-cleared each run).
- **parameters**: Chunking and anonymization settings.
- **thresholds**: QC and association thresholds.
- **flower**: Federated server address and rounds.

---

## Data & Configuration

### Data: Using Simulated Data for Each Center
- Simulated genotype data is generated and partitioned for each center (client) in the federated pipeline, see `pipeline/src/simulated_data/README_synthetic_data.md`
- For each center, you will have a directory such as `pipeline/src/simulated_data/simulated_data/tiny/center_1/` containing:
  - The PLINK files for that center (e.g., `tiny_center_1.bed`, `.bim`, `.fam`)
  - A `config.yaml` file specifying the paths and settings for that center
- Example directory structure for two centers:
  ```
  pipeline/src/simulated_data/simulated_data/tiny/center_1/
    ├── tiny_center_1.bed
    ├── tiny_center_1.bim
    ├── tiny_center_1.fam
    ├── config.yaml
  pipeline/src/simulated_data/simulated_data/tiny/center_2/
    ├── tiny_center_2.bed
    ├── tiny_center_2.bim
    ├── tiny_center_2.fam
    ├── config.yaml
  ```
- Each client runs with its own config and data, and outputs are written to per-client `intermediate` and `logs` directories.

### Client Configuration (`config.yaml`)
Use the following format for each center (client):
```yaml
input_data:
  path: "pipeline/src/simulated_data/simulated_data/tiny/center_1"
  type: "bed"

output:
  intermediate_dir: "pipeline/src/simulated_data/simulated_data/tiny/center_1/intermediate" # all the intermediate data files such as chunks
  log_dir: "pipeline/src/simulated_data/simulated_data/tiny/center_1/logs" # all the intermediate results files such as logs, plink output files

parameters:
  sample_offset: 1000000000000
  chunk_size: 100
  sample_chunk_size: 100
  snp_chunk_size: 100

thresholds:
  maf_threshold: 0.01
  missing_threshold: 0.01
  hwe_threshold: 1e-6
  p_threshold: 1e-3

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

#### Key Fields
- **input_data.path**: Path to the PLINK prefix for this center (no extension).
- **output.intermediate_dir**: Directory for intermediate files (per client, auto-cleared each run).
- **output.log_dir**: Directory for logs (per client, auto-cleared each run).
- **parameters**: Chunking and anonymization settings.
- **thresholds**: QC and association thresholds.
- **flower**: Federated server address and rounds.
- **participation**: Controls which stages this client participates in (all `true` for full participation).

---

## Output & Log Management (NEW)
- **All intermediate and log files are now written to the directories specified in `config.yaml` (`intermediate_dir`, `log_dir`).**
- **At the start of each client run, these directories are automatically cleared and recreated** to avoid redundant or stale outputs.
- **All log messages** (including stage progress, warnings, and errors) are written to `log_dir/iteration_log.txt` using a per-client logger.
- **No outputs are written to the current directory**—everything is organized per client.

---

## Running the Pipeline

### 1. Start the Server
```bash
python -m pipeline.main_server
```

### 2. Start Each Client
```bash
# Edit config.yaml as needed for each client
python -m pipeline.main_client pipeline/src/simulated_data/simulated_data/tiny/center_1/config.yaml
```
- Each client will generate a new unique client ID on each run.
- All outputs and logs will be stored in the specified directories, which are cleared at the start of each run.

---

## Federated Protocol & Workflow

### Stage-by-Stage Summary
1. **Key Exchange**: Clients generate DH keypairs, send public keys to server.
2. **Sync**: Clients receive all public keys and curve params, compute shared secrets, mask local seeds.
3. **Global QC**: Clients compute and mask QC arrays, send to server for aggregation.
4. **Global QC Response**: Server sends SNP exclusion list, clients update datasets.
5. **Init Chunks (KING)**: Server sends global seed, clients partition and anonymize data.
6. **Iterative KING**: Clients process kinship chunks, send results iteratively.
7. **Local LR**: Clients run local logistic regression, identify insignificant SNPs.
8. **Local LR Filter Response**: Server sends intersection of insignificant SNPs, clients filter datasets.
9. **Init Chunks (LR)**: Server sends global seed, clients re-partition for final LR.
10. **Iterative LR**: Clients process LR chunks, send results.
11. **Complete**: Workflow finished, results processed.

### Diagrams

#### Sequence Diagram
```mermaid
sequenceDiagram
    participant C1 as Client 1
    participant C2 as Client 2
    participant S as Server
    ...
```

#### Swimlane Diagram
```mermaid
flowchart TD
    ...
```

---

## Troubleshooting & Tips
- **Missing Curve Params in Sync Stage**: Ensure the server sends the required `curve` parameter in the config for the `sync` stage.
- **No Output in iteration_log.txt**: All logging must use `self.logger` (not the global `logging` module) in client code.
- **Redundant Outputs**: Each run clears the intermediate and log directories, so only current outputs are present.
- **File Not Found**: Check that all paths in `config.yaml` are correct and that PLINK is installed and in your `PATH`.
- **Reproducibility**: All configs, seeds, and outputs are managed per run for full reproducibility.

---

## Contributing
- Please open issues or pull requests for bug fixes, improvements, or new features.

---

## Acknowledgments
- Built with [Flower](https://flower.dev/), [PLINK](https://www.cog-genomics.org/plink/1.9/), and open-source Python tools.
