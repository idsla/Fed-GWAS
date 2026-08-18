# Real-World Dataset Experiments

This directory contains experiments using real-world genotype datasets (e.g., 1000 Genomes Project, HapMap) with simulated phenotypes for federated GWAS evaluation.

## Overview

For Application Notes submission, we use **simulated phenotypes on real genotypes**. This approach:
- Uses real genotype data (publicly available, no access restrictions)
- Generates phenotypes using realistic models (see `phenotype_generation/`)
- Provides full reproducibility (reviewers can download same data and reproduce)
- Demonstrates method on real genetic architecture

## Directory Structure

```
experiments/real_world/
├── README.md                    # This file
├── phenotype_generation/        # Phenotype generation tools
│   ├── models.py               # Phenotype models (additive, polygenic, etc.)
│   ├── generate_phenotypes.py  # Main generation tool
│   └── README.md              # Phenotype generation documentation
├── 1000genomes/                # 1000 Genomes Project experiment
│   ├── config.yaml            # Main experiment config
│   ├── configs/               # Per-client configs
│   ├── data/                  # Data directory
│   │   ├── README.md          # Download and setup instructions
│   │   └── phenotypes/        # Generated phenotypes
│   └── results/               # Results directory
└── hapmap/                    # HapMap experiment (optional)
    └── [similar structure]
```

## Adding a New Real-World Dataset

### 1. Create Dataset Directory

```bash
mkdir -p experiments/real_world/<dataset_name>/data
mkdir -p experiments/real_world/<dataset_name>/configs
mkdir -p experiments/real_world/<dataset_name>/results
```

### 2. Download Genotype Data

Place PLINK-format files (`.bed`, `.bim`, `.fam`) in the `data/` directory.

**Requirements**:
- PLINK 1.9/2.0 compatible format
- `.bed`, `.bim`, `.fam` files with matching prefix
- Genotypes encoded as 0, 1, 2 (homozygous ref, heterozygous, homozygous alt)
- Missing values encoded as -9 in `.fam` file

### 3. Filter Multiallelic Variants (Standard GWAS Preprocessing)

**Important**: Multiallelic variants (3+ alleles) must be excluded before analysis because PLINK 1.9 does not handle them well during merging, and most GWAS methods assume biallelic variants.

```bash
python experiments/tools/filter_multiallelic.py \
    --plink-prefix experiments/real_world/<dataset_name>/data/genotypes_raw \
    --output-prefix experiments/real_world/<dataset_name>/data/genotypes
```

This step should be performed **after** downloading/converting data but **before** phenotype generation and data partitioning.

### 4. Generate Phenotypes

Use the phenotype generation tool to create case/control phenotypes:

```bash
python -m experiments.real_world.phenotype_generation.generate_phenotypes \
    --plink-prefix experiments/real_world/<dataset_name>/data/genotypes \
    --output-fam experiments/real_world/<dataset_name>/data/genotypes_pheno.fam \
    --model additive \
    --causal-fraction 0.01 \
    --seed 42 \
    --metadata-file experiments/real_world/<dataset_name>/data/phenotype_metadata.json
```

This will:
- Load genotypes from PLINK files
- Generate binary phenotypes using the specified model
- Update the `.fam` file with phenotypes (1=control, 2=case)
- Save model metadata for reproducibility

### 5. Partition Data Across Clients

Partition the dataset by:
- **Population groups** (e.g., EUR, AFR, ASN for 1000 Genomes)
- **Sample indices** (even split, uneven split, etc.)
- **Geographic regions** (if available)

Create separate PLINK files for each client:
```bash
# Example: Partition by population
plink --bfile genotypes_pheno \
    --keep-fam center_1_samples.txt \
    --make-bed \
    --out center_1/genotypes

plink --bfile genotypes_pheno \
    --keep-fam center_2_samples.txt \
    --make-bed \
    --out center_2/genotypes
```

### 6. Create Configuration Files

Create `config.yaml` and per-client configs in `configs/center_X/config.yaml`:

```yaml
# config.yaml
experiment_name: <dataset_name>
experiment_category: real_world
description: 'Real-world dataset: <dataset_name> with simulated phenotypes'

clients:
  config_files:
    0: experiments/real_world/<dataset_name>/configs/center_1/config.yaml
    1: experiments/real_world/<dataset_name>/configs/center_2/config.yaml

data:
  data_dir: experiments/real_world/<dataset_name>/data
  partition_strategy: population  # or "even", "uneven", etc.

server:
  num_server_rounds: 15
  chunk_size: 100
  min_available_clients: 1
  min_fit_clients: 1

analysis:
  generate_baseline: true
  compare_results: true
  metrics_collection: true
```

### 7. Validate Data

Before running experiments, validate:
- PLINK files are readable
- Phenotypes are present in `.fam` file (column 6)
- Phenotype distribution is reasonable (not all cases or all controls)
- Client partitions are non-overlapping
- Sample IDs are unique within each client

Use the setup script:
```bash
python experiments/tools/setup_real_world_experiment.py \
    --dataset <dataset_name> \
    --validate
```

## Phenotype Requirements

GWAS requires both **genotype** and **phenotype** data:

- **Genotype**: `.bed`, `.bim`, `.fam` files (PLINK format)
- **Phenotype**: Case/control status in `.fam` file column 6:
  - `1` = control
  - `2` = case
  - `-9` or `0` = missing (will be excluded)

**Important**: The phenotype generation tool automatically converts binary phenotypes (0/1) to PLINK format (1/2).

## Reproducibility

For Application Notes submission:

1. **Document phenotype generation**:
   - Model type and parameters
   - Random seed value
   - Causal SNP selection method
   - Effect size distributions

2. **Save metadata**:
   - Model parameters (JSON file)
   - Causal SNP indices
   - Effect sizes
   - Case/control counts

3. **Include in methods**:
   - Clearly state: "Phenotypes were simulated using [model] on real [dataset] genotypes"
   - Provide seed values and model parameters
   - Make phenotype generation script available

## Available Datasets

### 1000 Genomes Project

- **Scale**: ~2,500 samples, ~84M variants
- **Partitioning**: By population groups (EUR, AFR, ASN)
- **Phenotypes**: Simulated using population-stratified model
- **See**: `1000genomes/README.md` for setup instructions

### HapMap (Optional)

- **Scale**: ~1,200 samples, ~1M SNPs
- **Partitioning**: By population groups
- **Phenotypes**: Simulated using additive model
- **Purpose**: Quick validation for reviewers

## Running Experiments

After setup, run experiments using the standard pipeline:

```bash
# Generate baseline (centralized PLINK)
python experiments/tools/generate_baseline.py \
    experiments/real_world/<dataset_name>/config.yaml

# Run federated experiment
flwr run . local-simulation \
    --run-config "config_path=experiments/real_world/<dataset_name>/configs simulation=true num-server-rounds=15"

# Analyze results
python experiments/tools/evaluation/qc/qc_evaluator.py \
    experiments/real_world/<dataset_name>/results \
    --baseline experiments/real_world/<dataset_name>/data/baseline

python experiments/tools/evaluation/lr/lr_evaluator.py \
    experiments/real_world/<dataset_name>/results \
    --baseline experiments/real_world/<dataset_name>/data/baseline
```

## Notes

- Real-world datasets may require significant disk space (1000 Genomes: ~100GB+)
- Download and preprocessing can take time
- Consider using subsets for initial testing
- Validate phenotype generation produces reasonable GWAS results
- Document all preprocessing steps for reproducibility
