# Federated GWAS Pipeline (Updated)

This repository implements a federated pipeline for Genome-Wide Association Studies (GWAS) using Flower, PLINK, and custom privacy-preserving protocols. The pipeline supports multi-stage, multi-client GWAS with robust output management, reproducibility, and clear logging.

---
## Environment Setup

### Option 1: Using UV to manage the environment (recommended)

Install uv: https://docs.astral.sh/uv/
```bash
pip install uv
# or
curl -LsSf https://astral.sh/uv/install.sh | sh
```

```bash
uv sync --python 3.11.13
```

Or create a virtual environment and then sync dependencies:

```bash
uv venv --python 3.11.13 # or simply uv venv
# Linux/Mac
source .venv/bin/activate
# Windows PowerShell
.venv\Scripts\activate

# Install dependencies
uv sync
# optional: install dev dependencies
uv sync --dev
```

For adding new packages:
```bash
uv add <package_name>
```

### Option 2: Conda Environment
We recommend [Miniconda](https://docs.conda.io/en/latest/miniconda/) or [Anaconda](https://www.anaconda.com/products/distribution) to manage your Python environment.

```bash
conda create -n fedgwas python=3.11 -y
conda activate fedgwas
```

---

## PLINK

- The pipeline requires [PLINK](https://www.cog-genomics.org/plink/1.9/) (version 1.9 or later).
- Download PLINK from: https://www.cog-genomics.org/plink/1.9/ based on your operating system.
- Place the PLINK binary (`plink` or `plink.exe`) in `bin` folder of the project.

---

## Pipeline Running Instruction (Simulation Mode)

### Data & Configuration Organization

Data and configuration are organized in `centers/center_xx/` folder. Detailed structure is as follows:
```
centers/
├── center_1/
    ├── config.yaml
    ├── data/
    │   ├── center_1.bed
    │   ├── center_1.bim
    │   ├── center_1.fam
    ├── logs/
    ├── intermediate/
├── center_2/
    ├── config.yaml
    ├── data/
    │   ├── center_2.bed
    │   ├── center_2.bim
    │   ├── center_2.fam
    ├── logs/
├── ...
```

Configuration file `config.yaml` for each center contains the following fields:
- `input_data.path`: Path to the input genotype data (PLINK prefix or VCF file).
- `output.intermediate_dir`: Directory for intermediate files (per client, auto-cleared each run).
- `output.log_dir`: Directory for logs (per client, auto-cleared each run).
- `parameters`: Chunking and anonymization settings.
- `thresholds`: QC and association thresholds.
- `flower`: Federated server address and rounds.

For data, the simulated genotype data is generated and partitioned for each center (client) in the federated pipeline, see `prd/README_synthetic_data.md`. Each center has its own data directory and config file and the data format is PLINK .bed/.bim/.fam.

Each client runs with its own config and data, and outputs are written to per-client `intermediate` and `logs` directories.

### Client Configuration (`config.yaml`)

We have put the config template in `configs/config_template.yaml`, you can copy it to `centers/center_xx/config.yaml` and modify the paths and parameters.

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

Key Fields in `config.yaml`
- **input_data.path**: Path to the PLINK prefix for this center (no extension).
- **output.intermediate_dir**: Directory for intermediate files (per client, auto-cleared each run).
- **output.log_dir**: Directory for logs (per client, auto-cleared each run).
- **parameters**: Chunking and anonymization settings.
- **thresholds**: QC and association thresholds.
- **flower**: Federated server address and rounds.
- **participation**: Controls which stages this client participates in (all `true` for full participation).

---

### How to Run the Pipeline

```bash
flwr run . local-simulation --stream
```

Running configurations are set in `pyproject.toml` e.g. `simulation` means running in simulation mode.

---

## Pipeline Running Instruction (Deployment Mode)

### Data & Configuration

Put your data in `user_space` folder as follows:

```
user_space/
├── config.yaml
├── data/
│   ├── center_1.bed
│   ├── center_1.bim
│   ├── center_1.fam
├── logs/
├── intermediate/
├── ...
```

Copy the `configs/config_template.yaml` to `user_space/config.yaml` and modify the paths and parameters.

```yaml
input_data:
  path: "user_data/tiny.bed"
  type: "bed"
```

### Running the Pipeline

#### 1. Lauch Flower Federation

##### 1.1 Start the Server Node in one terminal window
```bash
$ flower-superlink --insecure
```

#### 2. Start each Client in a separate terminal window, client config are in `pipeline/simulated_data/simulated_data/tiny/center_xx/config.yaml`, noted that at least two clients are required.

Initialize client 1
```bash
flower-supernode \
     --insecure \
     --superlink 127.0.0.1:9092 \
     --clientappio-api-address 127.0.0.1:9094 \
     --node-config "partition-id=0 num-partitions=2 config-file=configs/center_1/config.yaml"
```

Initialize client 2
```bash
flower-supernode \
     --insecure \
     --superlink 127.0.0.1:9092 \
     --clientappio-api-address 127.0.0.1:9095 \
     --node-config "partition-id=1 num-partitions=2 config-file=configs/center_2/config.yaml"
```

- Each client will generate a new unique client ID on each run.
- All outputs and logs will be stored in the specified directories, which are cleared at the start of each run.

### Run the Flower App on the Federation

We have set configuration in `pyproject.toml` for the Flower Federation, you can run the following command to start the our FedGWAS Pipeline App on the Federation:

```bash
flwr run . local-deployment --stream
```

### Other Notes:

#### Certificate and HTTP SSL Communication

```bash
flower-superlink \
    --ssl-ca-certfile <your-ca-cert-filepath> \
    --ssl-certfile <your-server-cert-filepath> \
    --ssl-keyfile <your-privatekey-filepath>

flower-supernode \
     --superlink 127.0.0.1:9092 \
     --clientappio-api-address 127.0.0.1:9094 \
     --root-certificates <your-ca-cert-filepath> \
     <other-args>

flower-supernode \
     --superlink 127.0.0.1:9092 \
     --clientappio-api-address 127.0.0.1:9095 \
     --root-certificates <your-ca-cert-filepath> \
     <other-args>
```

---

## Output & Log Management (NEW)
- **All intermediate and log files are now written to the directories specified in `config.yaml` (`intermediate_dir`, `log_dir`).**
- **At the start of each client run, these directories are automatically cleared and recreated** to avoid redundant or stale outputs.
- **All log messages** (including stage progress, warnings, and errors) are written to `log_dir/iteration_log.txt` using a per-client logger.
- **No outputs are written to the current directory**—everything is organized per client.

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
