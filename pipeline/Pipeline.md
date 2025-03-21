# Federated GWAS Pipeline

This repository implements a federated pipeline for Genome-Wide Association Studies (GWAS) using a combination of local processing and iterative outsourcing. The project is divided into two main directories: client/ and server/.

## Overview

The pipeline allows multiple clients to perform local quality control (QC), anonymize and partition their genotype data (in PLINK format), and then share partitial data/results with a central server. The server aggregates data/results, runs statistical analyses (e.g., QC, KING for kinship and logistic regression for association testing), and returns results to clients. Clients then use these results to filter out low-quality SNPs, remove highly related samples, and finalize association statistics.

## Repository Structure

```

├── client/
│   ├── base_client.py         # Base client class with methods for partitioning, anonymization, and file archiving.
│   ├── data_loader.py         # DataLoader class to read configuration from YAML and transform input data.
│   ├── iterative_king.py      # Implements iterative KING stage: sending data chunks and accumulating partial kinship results.
│   ├── iterative_lr.py        # Implements iterative Logistic Regression stage.
│   ├── local_qc.py            # Local quality control functions: computing genotype counts, missingness, and local LR.
│   ├── main_client.py         # Main entry point for the client side.
│   └── config.yaml            # YAML configuration file with parameters (input paths, thresholds, chunk sizes, etc.).
│
└── server/
    ├── aggregator_qc.py       # Aggregator for global QC: merging partial arrays, applying thresholds, and excluding SNPs.
    ├── aggregator_king.py     # Aggregator for iterative KING: merging client .tar chunks, running PLINK KING, and post-processing results.
    ├── aggregator_lr.py       # Aggregator for iterative LR: merging client .tar chunks, running PLINK logistic regression, and parsing results.
    ├── strategy.py            # Custom Flower strategy orchestrating all stages: sync, global QC, iterative KING, local LR, iterative LR.
    └── main_server.py         # Main entry point for the server side.
```

## Pipeline Stages

### 1. Local QC (Client)

- **Client Side:**  
  - Each client computes local per-sample missing rates using functions in `local_qc.py`.
  - Clients filter their local datasets with the threshold. 

### 2. Global QC (Client & Server)

- **Client Side:**  
  - Each client computes local genotype counts (for MAF and HWE) and per-SNP missing rates using functions in `local_qc.py`.
  - Clients load thresholds (e.g., MAF, missing rate, HWE thresholds) from the `config.yaml` file via the DataLoader and return these along with the partial QC arrays.
  
- **Server Side:**  
  - The server aggregates the partial QC arrays using a secure-sum placeholder (see `aggregator_qc.py`).
  - It applies the thresholds (received from the clients) to compute global MAF, missing rate, and HWE p-values.
  - SNP indices failing the thresholds are returned as an exclusion list.
  - This exclusion list is broadcast to the clients in the `global_qc_response` stage so they can filter their local datasets.

### 3. Sync

- **Client Side:**  
  - Each client generates a local seed.
  
- **Server Side:**  
  - The server securely aggregates these seeds (using a secure sum placeholder) to compute a global seed.
  - The global seed is broadcast back to the clients for use in data shuffling and partitioning.

### 4. Iterative KING

- **Client Side:**  
  - Clients partition their PLINK data (using methods in `base_client.py`) into disjoint chunks (in .tar files containing .bed/.bim/.fam) based on samples, SNPs, or both.
  - They send these chunks to the server.
  - Clients receive partial KING results (formatted as: `sampleID1_ano sampleID2_ano partial_phi n1_star`) from the server.
  - The clients accumulate these partial results and later finalize local kinship estimates; then, they filter out samples with high kinship.

- **Server Side:**  
  - The server receives the client .tar files, unpacks them, and merges them using PLINK’s `--bmerge`.
  - The server runs PLINK’s `--het` command on the merged data to compute heterozygosity (n1) for each sample.
  - Next, the server runs PLINK `--king robust` on the merged dataset.
  - It then post-processes the KING results by parsing the `.kin0` file, and using the computed n1 for each sample to produce output lines in the format:  
    `sampleID1_ano sampleID2_ano partial_phi n1_star`
  - These results are returned to the clients.

### 5. Iterative LR

- **Client Side:**  
  - Clients run local logistic regression to identify insignificant SNPs and send these SNP IDs to the server.
  - The server performs a secure intersection of these insignificant SNP sets and broadcasts the intersected set.
  - Clients filter out the insignificant SNPs and re-partition their data for iterative LR.
  - Finally, clients send iterative LR data chunks to the server.

- **Server Side:**  
  - The server merges the LR data chunks from clients (using PLINK’s `--bmerge`).
  - It then runs PLINK logistic regression on the merged data.
  - The server parses the resulting `.assoc.logistic` file to extract SNP IDs and p-values in the format:  
    `anonSNPID p_value`
  - These results are returned to the clients.

## Configuration

The `config.yaml` file (in the `client/` directory) defines parameters such as:

- **Input Data:**  
  - `path`: The path to the input genotype data (e.g., PLINK BED or VCF file).  
  - `type`: Input file type (e.g., `"bed"` or `"vcf"`).

- **Output Directories:**  
  - `intermediate_dir`: Directory for storing intermediate files.  
  - `log_dir`: Directory for storing logs.

- **Parameters:**  
  - `chunk_size`: Default size for partitioning data.  
  - `sample_chunk_size` / `snp_chunk_size`: Chunk sizes when partitioning by both.

- **Thresholds:**  
  - `maf_threshold`: Minimum minor allele frequency.  
  - `missing_threshold`: Maximum allowable missing rate.  
  - `hwe_threshold`: Minimum p-value for the HWE test.  
  - `p_threshold`: Threshold for insignificant SNPs in local LR.

- **Flower Configuration:**  
  - `server_address`: Address for the Flower server.  
  - `num_rounds`: Maximum number of rounds.

## Running the Pipeline

### Server

1. Navigate to the `server/` directory.
2. Start the server:
```bash
python main_server.py
```
### Client

1.  Navigate to the client/ directory.
2.  Ensure config.yaml is configured correctly.
3.  Start the client:
 ```bash
   python main_client.py
```
### Data Loading

The DataLoader class in data_loader.py:
    - reads the configuration from config.yaml.
    - Transforms input data if needed (e.g., converts a VCF file to PLINK binary format).
    - Provides configuration parameters (e.g., thresholds, chunk sizes, intermediate directories) to the client.

### Secure Aggregation

Both the global seed computation and QC data merging use placeholder functions (e.g., secure_sum_arrays in aggregator_qc.py and secure_sum in strategy.py). 

**Replace these placeholders with chosen secure aggregation protocol**.

### Anonymization and Partitioning

Clients anonymize sample IDs (using a client-specific random offset) and SNP IDs (using a global seed for consistent hashing).

Data is partitioned into disjoint chunks (by samples, SNPs, or both) and archived as .tar files.

This ensures privacy while enabling the server to merge data reliably.


