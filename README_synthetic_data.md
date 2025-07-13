# Synthetic GWAS Data Generation

This directory contains scripts to generate synthetic GWAS datasets for federated analysis testing.

## Overview

The synthetic data generator creates PLINK-format datasets at four scales:

| Scale | Samples | SNPs | Centers | Purpose |
|-------|---------|------|---------|---------|
| Tiny | 500 | 5,000 | 2 | Sanity checks and development |
| Small | 5,000 | 50,000 | 3 | Functional testing |
| Medium | 50,000 | 500,000 | 5 | Performance tuning |
| Large | 275,812 | 98,000,000 | 7 | Scalability benchmarking |

## Features

- **Binary phenotypes** generated using logistic model
- **1% causal SNPs** with effect sizes ~N(0, 0.1)
- **2% missing data** randomly distributed
- **MAFs** drawn from Uniform(0.01, 0.5)
- **Horizontal partitioning** by samples across centers
- **Automatic PLINK integration** for center-specific file generation

## Prerequisites

1. **PLINK** installed and available in PATH
2. **Python dependencies**:
   ```bash
   pip install -r requirements_synthetic.txt
   ```

## Usage

### Generate all scales
```bash
python generate_synthetic_data.py
```

### Generate specific scale
```bash
python generate_synthetic_data.py --scale tiny
python generate_synthetic_data.py --scale small
python generate_synthetic_data.py --scale medium
python generate_synthetic_data.py --scale large
```

### Custom output directory
```bash
python generate_synthetic_data.py --output-dir /path/to/output
```

## Output Structure

```
simulated_data/
├── tiny/
│   ├── tiny.bed
│   ├── tiny.bim
│   ├── tiny.fam
│   ├── metadata.json
│   └── center_1/
│       ├── samples.txt
│       ├── tiny_center_1.bed
│       ├── tiny_center_1.bim
│       └── tiny_center_1.fam
│   └── center_2/
│       └── ...
├── small/
│   └── ...
├── medium/
│   └── ...
└── large/
    └── ...
```

## File Descriptions

- **`*.bed/.bim/.fam`**: PLINK-format genotype files
- **`metadata.json`**: Dataset information including causal SNPs and phenotype counts
- **`center_*/`**: Per-center partitioned datasets
- **`samples.txt`**: Sample ID lists for each center

## Notes

- The Large scale may require significant disk space (~5TB for full UKB-scale)
- Consider using a subset for disk-constrained environments
- All centers share identical SNP sets; only sample partitioning differs
- Phenotype distributions are balanced across centers by default

## Validation

After generation, you can validate the files using PLINK:

```bash
plink --bfile simulated_data/tiny/tiny --freq
plink --bfile simulated_data/tiny/tiny --assoc
``` 