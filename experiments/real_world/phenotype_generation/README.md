# Phenotype Generation Tools

This directory contains tools for generating phenotypes from real genotype data (e.g., 1000 Genomes Project, HapMap) for federated GWAS evaluation.

## Overview

For Application Notes submission, we use **simulated phenotypes on real genotypes**. This approach:
- Uses real genotype data (publicly available, no access restrictions)
- Generates phenotypes using realistic models
- Provides full reproducibility (reviewers can download same data and reproduce)
- Demonstrates method on real genetic architecture

## Models

### 1. Additive Model (`AdditiveModel`)

Simple additive logistic model:
```
P(Y=1) = logistic(β₀ + Σᵢ βᵢGᵢ + ε)
```

**Parameters**:
- `causal_fraction`: Fraction of SNPs that are causal (default: 0.01)
- `effect_size_std`: Standard deviation of effect sizes (default: 0.1)
- `noise_std`: Standard deviation of random noise (default: 1.0)
- `intercept`: Intercept term (default: 0.0)

### 2. Polygenic Model (`PolygenicModel`)

Polygenic model with explicit heritability control:
```
P(Y=1) = logistic(β₀ + G + E)
```

Where G is genetic component and E is environmental component, scaled to achieve target heritability.

**Parameters**:
- `causal_fraction`: Fraction of SNPs that are causal
- `heritability`: Target heritability h² (default: 0.3)
- `intercept`: Intercept term

### 3. Population-Stratified Model (`PopulationStratifiedModel`)

Adds population-specific effects to test population stratification handling:
```
P(Y=1) = logistic(β₀ + Σᵢ βᵢGᵢ + Pop_effect + ε)
```

**Parameters**:
- `causal_fraction`: Fraction of SNPs that are causal
- `effect_size_std`: Standard deviation of SNP effect sizes
- `population_effect_std`: Standard deviation of population effects (default: 0.2)
- `noise_std`: Standard deviation of noise

## Usage

### Command Line

```bash
python -m experiments.real_world.phenotype_generation.generate_phenotypes \
    --plink-prefix data/1000genomes \
    --output-fam data/1000genomes_pheno.fam \
    --model additive \
    --causal-fraction 0.01 \
    --seed 42 \
    --metadata-file data/phenotype_metadata.json
```

### Python API

```python
from experiments.real_world.phenotype_generation.models import AdditiveModel
import numpy as np

# Load genotypes (n_samples, n_snps) with values 0, 1, 2, or -9 (missing)
genotypes = load_genotypes("data/1000genomes")

# Initialize model
model = AdditiveModel(
    causal_fraction=0.01,
    effect_size_std=0.1,
    seed=42
)

# Generate phenotypes
phenotypes, metadata = model.generate(genotypes)

# phenotypes: binary array (0=control, 1=case)
# metadata: dictionary with model parameters and causal SNPs
```

## Output

The tool generates:
1. **Updated .fam file**: Original .fam file with phenotype column updated (1=control, 2=case)
2. **Metadata file** (optional): JSON file with:
   - Model parameters
   - Causal SNP indices
   - Effect sizes
   - Case/control counts
   - Case rate

## Example: 1000 Genomes Project

```bash
# Download 1000 Genomes data (see 1000genomes/README.md)
# Then generate phenotypes:

python -m experiments.real_world.phenotype_generation.generate_phenotypes \
    --plink-prefix experiments/real_world/1000genomes/data/genotypes \
    --output-fam experiments/real_world/1000genomes/data/genotypes_pheno.fam \
    --model population_stratified \
    --causal-fraction 0.01 \
    --population-file experiments/real_world/1000genomes/data/population_labels.txt \
    --seed 42 \
    --metadata-file experiments/real_world/1000genomes/data/phenotype_metadata.json
```

## Notes

- Phenotypes are generated deterministically with a fixed seed
- Model parameters should be documented in Application Notes methods section
- For reproducibility, include seed values and model parameters in experiment configs
- Generated phenotypes should produce reasonable GWAS results (validate with baseline)
