# Federated GWAS Pipeline

This repository implements a federated pipeline for Genome-Wide Association Studies (GWAS) using a combination of local processing and iterative outsourcing. The project is divided into two main directories: client/ and server/.

---

## Environment Setup

### 1. Create and Activate a Conda Environment

We recommend using [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution) to manage your Python environment.

```bash
conda create -n fedgwas python=3.11 -y
conda activate fedgwas
```

### 2. Install Poetry and Project Dependencies
- Install Poetry (if not already installed):
  ```bash
  pip install poetry
  ```
- Install all required packages:
  ```bash
  poetry install
  ```

### 3. (Optional) Install PLINK
- Download [PLINK 1.9+](https://www.cog-genomics.org/plink/1.9/).
- Place the `plink` binary in your `PATH` or specify its location in your scripts.

---

## PLINK

- The pipeline requires [PLINK](https://www.cog-genomics.org/plink/1.9/) (version 1.9 or later).
- Download PLINK from: https://www.cog-genomics.org/plink/1.9/
- Place the PLINK binary (`plink` or `plink.exe`) in a directory included in your system `PATH`, or specify its path in your scripts if needed.

---

## Data

- Place your input genotype data (PLINK .bed/.bim/.fam or VCF files) in a directory such as `data/client_data/`.
- Example data can be downloaded from [PLINK test data](https://www.cog-genomics.org/plink/1.9/resources#test).
- Each client should have its own data directory.

---

## Client Configuration (`config.yaml`)

Example (`pipeline/src/clients/config.yaml`):
```yaml
input_data:
  path: "data/client_data"   # Path to PLINK data (without extension) or VCF file etc.
  type: "bed"               # "bed" or "vcf"

output:
  intermediate_dir: "intermediate"
  log_dir: "logs"

parameters:
  sample_offset: 1000000000000   # if not provided by DataLoader, default offset
  chunk_size: 100                # default chunk size (by samples or SNPs)
  sample_chunk_size: 100
  snp_chunk_size: 100

thresholds:
  maf_threshold: 0.01
  missing_threshold: 0.1
  hwe_threshold: 1e-6
  p_threshold: 1e-3

flower:
  server_address: "127.0.0.1:8080"
  num_rounds: 10
```

### Field Descriptions
- **input_data.path**: Path to the input genotype data (PLINK prefix or VCF file).
- **input_data.type**: File type, either `bed` or `vcf`.
- **output.intermediate_dir**: Directory for intermediate files.
- **output.log_dir**: Directory for logs.
- **parameters.sample_offset**: Offset for anonymizing sample IDs.
- **parameters.chunk_size**: Default chunk size for partitioning.
- **parameters.sample_chunk_size**: Chunk size for sample-based partitioning.
- **parameters.snp_chunk_size**: Chunk size for SNP-based partitioning.
- **thresholds.maf_threshold**: Minimum minor allele frequency.
- **thresholds.missing_threshold**: Maximum missing rate.
- **thresholds.hwe_threshold**: Minimum HWE p-value.
- **thresholds.p_threshold**: P-value threshold for local LR.
- **flower.server_address**: Address of the Flower server.
- **flower.num_rounds**: Number of federated rounds.

---

## Master Configuration
- The pipeline currently uses per-client configuration files. If a master configuration is needed, it should specify global parameters and reference each client’s config.

---

## Running Instructions

### 1. Start the Server
```bash
cd pipeline/src/server
python main_server.py
```

### 2. Start Each Client
```bash
cd pipeline/src/clients
# Edit config.yaml as needed for each client
python main_client.py
```

---

## Results
- Intermediate and final results are stored in the directories specified in each client’s `config.yaml` (`intermediate_dir`, `log_dir`).
- Server-side results (aggregated QC, kinship, association) are stored in the server’s working directory.
- Check logs and output directories for detailed results and logs of each stage.

---

# Pipeline Details

## Requirements

### Set up the environment

```bash
pip install -r requirements.txt
```

### Download Plink

```bash
wget https://www.cog-genomics.org/plink/download/plink2.zip
unzip plink2.zip
```

## Run the Pipeline

### CLient Side

```bash
source venv/bin/activate
python main_client.py
```

### Server Side

```bash
source venv/bin/activate
python main_server.py
```

## Workflow Sequence Diagram

```mermaid
sequenceDiagram
    participant C1 as Client 1
    participant C2 as Client 2
    participant S as Server
    
    Note over C1,S: 1. Key Exchange Stage
    C1->>S: DH Public Key (uint8 array)
    C2->>S: DH Public Key (uint8 array)
    Note over S: Store public keys<br/>Check if exchange complete
    
    Note over C1,S: 2. Sync Stage
    S->>C1: All Public Keys + DH Params
    S->>C2: All Public Keys + DH Params
    C1->>S: Masked Local Seed (PRG-MASKING)
    C2->>S: Masked Local Seed (PRG-MASKING)
    Note over S: Aggregate masked seeds<br/>Compute global seed (masks cancel)
    
    Note over C1,S: 3. Global QC Stage
    C1->>S: Masked QC Arrays<br/>[counts, missing, thresholds]
    C2->>S: Masked QC Arrays<br/>[counts, missing, thresholds]
    Note over S: Aggregate QC data<br/>Determine SNP exclusions
    
    Note over C1,S: 4. Global QC Response Stage
    S->>C1: SNP Exclusion List
    S->>C2: SNP Exclusion List
    Note over C1: Apply exclusions<br/>Update local dataset
    Note over C2: Apply exclusions<br/>Update local dataset
    
    Note over C1,S: 5. Init Chunks (KING) Stage
    S->>C1: Global Seed
    S->>C2: Global Seed
    Note over C1: Partition data into chunks<br/>Anonymize SNP/sample IDs
    Note over C2: Partition data into chunks<br/>Anonymize SNP/sample IDs
    
    Note over C1,S: 6. Iterative KING Stage
    loop For each chunk
        C1->>S: KING chunk results
        C2->>S: KING chunk results
        Note over S: Aggregate KING statistics<br/>Process relationship data
    end
    
    Note over C1,S: 7. Local LR Stage
    C1->>S: Insignificant SNP Set
    C2->>S: Insignificant SNP Set
    Note over S: Compute intersection<br/>of insignificant SNPs
    
    Note over C1,S: 8. Local LR Filter Response Stage
    S->>C1: SNP Intersection Set
    S->>C2: SNP Intersection Set
    Note over C1: Apply LR filtering<br/>Update dataset
    Note over C2: Apply LR filtering<br/>Update dataset
    
    Note over C1,S: 9. Init Chunks (LR) Stage
    S->>C1: Global Seed (for LR)
    S->>C2: Global Seed (for LR)
    Note over C1: Re-partition filtered data<br/>Prepare for final LR
    Note over C2: Re-partition filtered data<br/>Prepare for final LR
    
    Note over C1,S: 10. Iterative LR Stage
    loop For each LR chunk
        C1->>S: LR chunk results
        C2->>S: LR chunk results
        Note over S: Aggregate LR statistics<br/>Process association results
    end
    
    Note over C1,S: 11. Workflow Complete
    Note over S: Final results processed<br/>Stage: done
```

## Detailed Workflow Swimlane Diagram

```mermaid
flowchart TD
    subgraph "Client 1 Operations"
        C1_1[Generate DH Keypair<br/>Send Public Key]
        C1_2[Receive All Public Keys<br/>Compute Shared Secrets<br/>Mask Local Seed]
        C1_3[Load PLINK Data<br/>Compute QC Arrays<br/>Apply PRG Masking]
        C1_4[Receive Exclusion List<br/>Filter SNPs<br/>Update Dataset]
        C1_5[Receive Global Seed<br/>Anonymize IDs<br/>Partition into Chunks]
        C1_6[Process KING Chunks<br/>Send Results]
        C1_7[Run Local LR<br/>Identify Insignificant SNPs]
        C1_8[Receive SNP Intersection<br/>Apply LR Filtering]
        C1_9[Re-partition Filtered Data<br/>Prepare for Final LR]
        C1_10[Process LR Chunks<br/>Send Final Results]
        C1_11[Workflow Complete]
    end

    subgraph "Server Operations"
        S_1[Collect Public Keys<br/>Check Exchange Complete]
        S_2[Distribute Public Keys<br/>Aggregate Masked Seeds<br/>Compute Global Seed]
        S_3[Aggregate Masked QC Data<br/>Apply Thresholds<br/>Generate Exclusion List]
        S_4[Send Exclusion List<br/>Wait for Client Updates]
        S_5[Distribute Global Seed<br/>Coordinate Chunking]
        S_6[Aggregate KING Results<br/>Process Relationships]
        S_7[Collect Insignificant SNPs<br/>Compute Intersection]
        S_8[Send SNP Intersection<br/>Coordinate Filtering]
        S_9[Distribute Global Seed<br/>Prepare LR Phase]
        S_10[Aggregate LR Results<br/>Process Associations]
        S_11[Finalize Results<br/>Set Stage: done]
    end

    subgraph "Client 2 Operations"
        C2_1[Generate DH Keypair<br/>Send Public Key]
        C2_2[Receive All Public Keys<br/>Compute Shared Secrets<br/>Mask Local Seed]
        C2_3[Load PLINK Data<br/>Compute QC Arrays<br/>Apply PRG Masking]
        C2_4[Receive Exclusion List<br/>Filter SNPs<br/>Update Dataset]
        C2_5[Receive Global Seed<br/>Anonymize IDs<br/>Partition into Chunks]
        C2_6[Process KING Chunks<br/>Send Results]
        C2_7[Run Local LR<br/>Identify Insignificant SNPs]
        C2_8[Receive SNP Intersection<br/>Apply LR Filtering]
        C2_9[Re-partition Filtered Data<br/>Prepare for Final LR]
        C2_10[Process LR Chunks<br/>Send Final Results]
        C2_11[Workflow Complete]
    end

    %% Stage 1: Key Exchange
    C1_1 --> S_1
    C2_1 --> S_1
    
    %% Stage 2: Sync
    S_1 --> S_2
    S_2 --> C1_2
    S_2 --> C2_2
    C1_2 --> S_2
    C2_2 --> S_2
    
    %% Stage 3: Global QC
    S_2 --> S_3
    S_3 --> C1_3
    S_3 --> C2_3
    C1_3 --> S_3
    C2_3 --> S_3
    
    %% Stage 4: Global QC Response
    S_3 --> S_4
    S_4 --> C1_4
    S_4 --> C2_4
    
    %% Stage 5: Init Chunks (KING)
    S_4 --> S_5
    S_5 --> C1_5
    S_5 --> C2_5
    
    %% Stage 6: Iterative KING
    C1_5 --> C1_6
    C2_5 --> C2_6
    C1_6 --> S_6
    C2_6 --> S_6
    
    %% Stage 7: Local LR
    S_6 --> S_7
    S_7 --> C1_7
    S_7 --> C2_7
    C1_7 --> S_7
    C2_7 --> S_7
    
    %% Stage 8: Local LR Filter Response
    S_7 --> S_8
    S_8 --> C1_8
    S_8 --> C2_8
    
    %% Stage 9: Init Chunks (LR)
    S_8 --> S_9
    S_9 --> C1_9
    S_9 --> C2_9
    
    %% Stage 10: Iterative LR
    C1_9 --> C1_10
    C2_9 --> C2_10
    C1_10 --> S_10
    C2_10 --> S_10
    
    %% Stage 11: Complete
    S_10 --> S_11
    S_11 --> C1_11
    S_11 --> C2_11

    %% Styling
    classDef client1 fill:#e1f5fe
    classDef server fill:#f3e5f5
    classDef client2 fill:#e8f5e8
    
    class C1_1,C1_2,C1_3,C1_4,C1_5,C1_6,C1_7,C1_8,C1_9,C1_10,C1_11 client1
    class S_1,S_2,S_3,S_4,S_5,S_6,S_7,S_8,S_9,S_10,S_11 server
    class C2_1,C2_2,C2_3,C2_4,C2_5,C2_6,C2_7,C2_8,C2_9,C2_10,C2_11 client2
```

### Stage-by-Stage Operations Summary

**Stage 1 - Key Exchange:**
- *Clients*: Generate DH keypairs using server parameters, send public keys as uint8 arrays
- *Server*: Distribute DH parameters, collect and store client public keys, check completion

**Stage 2 - Sync:**
- *Clients*: Receive all public keys, compute shared secrets, apply PRG masking to local seeds
- *Server*: Distribute public keys, aggregate masked seeds (masks cancel out), compute global seed

**Stage 3 - Global QC:**
- *Clients*: Load PLINK data, compute genotype/missing counts, apply PRG masking to QC arrays
- *Server*: Aggregate masked QC data, apply thresholds (MAF, missing rate, HWE), generate exclusion list

**Stage 4 - Global QC Response:**
- *Clients*: Receive SNP exclusion list, filter local datasets, update PLINK files
- *Server*: Send exclusion list to all participating clients

**Stage 5 - Init Chunks (KING):**
- *Clients*: Receive global seed, anonymize SNP/sample IDs, partition data into chunks
- *Server*: Distribute global seed for deterministic shuffling

**Stage 6 - Iterative KING:**
- *Clients*: Process each chunk for kinship analysis, send results iteratively
- *Server*: Aggregate KING statistics, process relationship data

**Stage 7 - Local LR:**
- *Clients*: Run local logistic regression, identify insignificant SNPs (p > threshold)
- *Server*: Collect insignificant SNP sets, compute intersection across clients

**Stage 8 - Local LR Filter Response:**
- *Clients*: Receive SNP intersection, apply filtering to remove insignificant SNPs
- *Server*: Send intersection of insignificant SNPs to all clients

**Stage 9 - Init Chunks (LR):**
- *Clients*: Re-partition filtered datasets, prepare for final LR analysis
- *Server*: Distribute global seed for LR phase chunking

**Stage 10 - Iterative LR:**
- *Clients*: Process LR chunks, send final association results
- *Server*: Aggregate final LR statistics, process association results

**Stage 11 - Complete:**
- *Clients*: Workflow finished, final results processed
- *Server*: Set stage to "done", finalize all results

