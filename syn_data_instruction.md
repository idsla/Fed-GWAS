# Synthetic Data Generation Instructions

This document describes how to generate synthetic GWAS datasets using our custom Python script and PLINK, and provides links to generated data files.

## Generation Method

We use a custom Python script (`generate_synthetic_data.py`) that:
1. Generates synthetic genotypes and binary phenotypes
2. Creates PLINK-format files (.bed/.bim/.fam)
3. Automatically partitions data across centers
4. Calls PLINK to create center-specific datasets

## PLINK Commands Used

The script automatically calls PLINK for each center:

```bash
plink --bfile <input_prefix> --keep <sample_list> --make-bed --out <output_prefix>
```

Where:
- `<input_prefix>`: Path to the full dataset (e.g., `simulated_data/tiny/tiny`)
- `<sample_list>`: Text file with sample IDs for this center
- `<output_prefix>`: Output path for center-specific files

## Download Links

| Size   | #Samples | #SNPs   | Filename | Data Link | Generation Command |
|--------|----------|---------|----------|-----------|-------------------|
| Tiny   | 500      | 5,000   | tiny.zip | [Google Drive/OneDrive](#) | `python generate_synthetic_data.py --scale tiny` |
| Small  | 5,000    | 50,000  | small.zip | [Google Drive/OneDrive](#) | `python generate_synthetic_data.py --scale small` |
| Medium | 50,000   | 500,000 | medium.zip | [Google Drive/OneDrive](#) | `python generate_synthetic_data.py --scale medium` |
| Large  | 275,812  | 98,000,000 | large.zip | [Google Drive/OneDrive](#) | `python generate_synthetic_data.py --scale large` |

## File Structure

Each dataset contains:
- `*.bed/.bim/.fam`: Full PLINK dataset
- `metadata.json`: Dataset information and causal SNPs
- `center_*/`: Per-center partitioned datasets
  - `samples.txt`: Sample ID list
  - `*_center_*.bed/.bim/.fam`: Center-specific PLINK files

> **Note:** Replace the data links after uploading the generated files to cloud storage. 