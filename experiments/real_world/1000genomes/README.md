# 1000 Genomes Project Experiment

This experiment uses the 1000 Genomes Project genotype data with simulated phenotypes for federated GWAS evaluation.

## Dataset Overview

- **Source**: 1000 Genomes Project Phase 3
- **Scale**: ~2,500 samples, ~84M variants (can use subsets)
- **Access**: Publicly available, no restrictions
- **URL**: https://www.internationalgenome.org/data

## Setup Instructions

### 1. Download Genotype Data

The 1000 Genomes Project provides data in multiple formats. For PLINK format:

**Option A: Download pre-processed PLINK files** (if available)
```bash
# Check 1000 Genomes website for PLINK-format downloads
# Place files in: experiments/real_world/1000genomes/data/
```

**Option B: Convert from VCF to PLINK**
```bash
# Download VCF files from 1000 Genomes
# Convert to PLINK format:
plink --vcf ALL.chr*.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf \
    --make-bed \
    --out experiments/real_world/1000genomes/data/genotypes_raw

# Note: This requires significant disk space (~100GB+)
```

**Option C: Use subset for testing**
```bash
# Download chromosome 22 only (smaller, faster):
# Run: bash experiments/real_world/1000genomes/download_subset.sh
# Then convert and filter as needed
```

### 2. Filter Multiallelic Variants (Standard GWAS Preprocessing)

**Important**: Multiallelic variants (3+ alleles) must be excluded before analysis because:
- PLINK 1.9 does not handle them well during merging
- Most GWAS methods assume biallelic variants
- They can cause merge errors and analysis issues

Filter multiallelic variants using the preprocessing tool:

```bash
python experiments/tools/filter_multiallelic.py \
    --plink-prefix experiments/real_world/1000genomes/data/genotypes_raw \
    --output-prefix experiments/real_world/1000genomes/data/genotypes
```

This tool:
- Uses PLINK2 `--max-alleles 2` if available (most reliable)
- Falls back to PLINK 1.9 `--freq` detection if PLINK2 is not available
- Creates filtered PLINK files with only biallelic variants

**Note**: This step should be performed **after** VCF conversion but **before** phenotype generation and data partitioning.

### 3. Get Population Labels

1000 Genomes samples are labeled by population (EUR, AFR, ASN, etc.).

Download population information:
```bash
# Population labels file (one label per sample, same order as .fam)
# Available from 1000 Genomes project website
# Save as: experiments/real_world/1000genomes/data/population_labels.txt
```

### 4. Generate Phenotypes

Use the phenotype generation tool with population-stratified model:

```bash
python -m experiments.real_world.phenotype_generation.generate_phenotypes \
    --plink-prefix experiments/real_world/1000genomes/data/genotypes \
    --output-fam experiments/real_world/1000genomes/data/genotypes_pheno.fam \
    --model population_stratified \
    --causal-fraction 0.01 \
    --population-file experiments/real_world/1000genomes/data/population_labels.txt \
    --seed 42 \
    --metadata-file experiments/real_world/1000genomes/data/phenotype_metadata.json
```

This will:
- Generate case/control phenotypes using population-stratified model
- Update `.fam` file with phenotypes
- Save model metadata for reproducibility

### 5. Partition by Population Groups

Partition samples by population to simulate multi-site federated scenario:

```bash
# Create sample lists for each population group
# Center 1: European (EUR)
grep -E "EUR" experiments/real_world/1000genomes/data/genotypes_pheno.fam | \
    awk '{print $1, $2}' > experiments/real_world/1000genomes/data/center_1_samples.txt

# Center 2: African (AFR)
grep -E "AFR" experiments/real_world/1000genomes/data/genotypes_pheno.fam | \
    awk '{print $1, $2}' > experiments/real_world/1000genomes/data/center_2_samples.txt

# Optional: Center 3: Asian (ASN)
grep -E "ASN|EAS|SAS" experiments/real_world/1000genomes/data/genotypes_pheno.fam | \
    awk '{print $1, $2}' > experiments/real_world/1000genomes/data/center_3_samples.txt

# Create partitioned PLINK files
plink --bfile experiments/real_world/1000genomes/data/genotypes_pheno \
    --keep experiments/real_world/1000genomes/data/center_1_samples.txt \
    --make-bed \
    --out experiments/real_world/1000genomes/data/center_1/genotypes

plink --bfile experiments/real_world/1000genomes/data/genotypes_pheno \
    --keep experiments/real_world/1000genomes/data/center_2_samples.txt \
    --make-bed \
    --out experiments/real_world/1000genomes/data/center_2/genotypes
```

### 6. Create Configuration Files

Create `config.yaml`:
```yaml
experiment_name: 1000genomes
experiment_category: real_world
description: '1000 Genomes Project with simulated phenotypes (population-stratified model)'

clients:
  config_files:
    0: experiments/real_world/1000genomes/configs/center_1/config.yaml
    1: experiments/real_world/1000genomes/configs/center_2/config.yaml

data:
  data_dir: experiments/real_world/1000genomes/data
  partition_strategy: population

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

Create per-client configs in `configs/center_X/config.yaml`:
```yaml
data:
  plink_prefix: "experiments/real_world/1000genomes/data/center_X/genotypes"

output:
  intermediate_dir: "experiments/real_world/1000genomes/results/center_X/intermediate"
  log_dir: "experiments/real_world/1000genomes/results/center_X/logs"

# ... other config parameters
```

### 7. Validate Setup

```bash
# Check PLINK files
plink --bfile experiments/real_world/1000genomes/data/center_1/genotypes --check

# Verify phenotypes
awk '{print $6}' experiments/real_world/1000genomes/data/center_1/genotypes.fam | sort | uniq -c
# Should show counts for 1 (control) and 2 (case)

# Use setup script
python experiments/tools/setup_real_world_experiment.py \
    --dataset 1000genomes \
    --validate
```

## Running the Experiment

### 1. Generate Baseline

```bash
python experiments/tools/generate_baseline.py \
    experiments/real_world/1000genomes/config.yaml
```

### 2. Run Federated Experiment

**Simulation mode**:
```bash
flwr run . local-simulation \
    --run-config "config_path=experiments/real_world/1000genomes/configs simulation=true num-server-rounds=15"
```

**Deployment mode** (cluster):
```bash
# Follow cluster_deployment/docs/CLUSTER_DEPLOYMENT.md
# Update configs with cluster IPs
```

### 3. Analyze Results

```bash
python experiments/tools/evaluation/qc/qc_evaluator.py \
    experiments/real_world/1000genomes/results \
    --baseline experiments/real_world/1000genomes/data/baseline

python experiments/tools/evaluation/lr/lr_evaluator.py \
    experiments/real_world/1000genomes/results \
    --baseline experiments/real_world/1000genomes/data/baseline
```

## Expected Results

- **Correctness**: Federated results should match centralized baseline
  - QC SNP exclusion agreement: F1 > 0.95
  - KING correlation: Pearson r > 0.99
  - LR correlation: Pearson r > 0.99

- **Population Stratification**: Should be handled correctly
  - No spurious associations due to population structure
  - KING kinship detects relatedness across populations

## Notes

- **Disk Space**: Full 1000 Genomes dataset requires ~100GB+ disk space
- **Memory**: Large datasets may require significant RAM
- **Time**: Full dataset processing can take hours/days
- **Subsets**: Consider using chromosome-specific or sample subsets for testing
- **Documentation**: Document all preprocessing steps for Application Notes

## Phenotype Model Details

For Application Notes, document:
- **Model**: Population-stratified logistic model
- **Causal fraction**: 0.01 (1% of SNPs are causal)
- **Population effects**: Included to test stratification handling
- **Seed**: 42 (for reproducibility)
- **Case rate**: ~0.5 (balanced case/control)

See `phenotype_metadata.json` for full model parameters.
