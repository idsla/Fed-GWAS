# Synthetic GWAS Data Generation

---

**Table of Contents**
- [Overview](#overview)
- [Simulation Logic](#simulation-logic)
- [Partitioning Scenarios](#partitioning-scenarios)
- [Output Structure](#output-structure)
- [File Descriptions](#file-descriptions)
- [Assigned Samples for Benchmarking](#assigned-samples-for-benchmarking)
- [Centralized Ground Truth Results](#centralized-ground-truth-results)
- [Per-Center Sample ID Mapping](#per-center-sample-id-mapping)
- [Code Structure](#code-structure)
- [Usage](#usage)
- [Reproducibility and Random Seed](#reproducibility-and-random-seed)
- [Test Coverage](#test-coverage)
- [Validation](#validation)
- [Interpreting Relatives and Kinship](#interpreting-relatives-and-kinship)
- [Extending or Customizing the Generator](#extending-or-customizing-the-generator)
- [Troubleshooting & Common Issues](#troubleshooting--common-issues)
- [References](#references)

---

## Overview

This directory contains scripts to generate synthetic GWAS datasets for federated analysis testing, benchmarking, and development. Datasets are generated at four scales:

| Scale  | Samples  | SNPs      | Centers | Purpose                  |
|--------|----------|-----------|---------|--------------------------|
| Tiny   | 500      | 5,000     | 2       | Sanity checks, dev       |
| Small  | 5,000    | 50,000    | 3       | Functional testing       |
| Medium | 50,000   | 500,000   | 5       | Performance tuning       |
| Large  | 275,812  | 98,000,000| 7       | Scalability benchmarking |

---

## Simulation Logic

- **Genotypes**: SNPs are simulated with MAFs drawn from Uniform(0.01, 0.5).
- **Relatives**: ~10-15% of samples are relatives, generated using realistic Mendelian inheritance logic (see `create_relatives`).
  - First-degree, second-degree, and third-degree relatives are included.
  - Relatives are tracked in `relatives_info.csv`.
- **Phenotypes**: Binary phenotypes are generated using a logistic model with 1% causal SNPs (effect sizes ~N(0, 0.1)).
- **Missingness**: 2% of genotype data is set to missing at random.
- **Partitioning**: Samples are split across centers using one of three strategies (see below).
- **PLINK Integration**: PLINK is called automatically to create center-specific .bed/.bim/.fam files.

---

## Partitioning Scenarios

1. **Even Sample Split (`even`)**
   - Equal samples per center, balanced phenotype (~50% cases/controls per center)
   - Use: Functional/correctness testing
2. **Uneven Sample Sizes (`uneven`)**
   - Different sample counts per center, but balanced phenotype in each
   - Use: Sample size imbalance effects
3. **Uneven Size + Phenotype Skew (`skewed`)**
   - Both sample size and phenotype ratios differ by center
   - Use: Phenotype distribution imbalance effects

See [syn_data_instruction.md](syn_data_instruction.md) for detailed PLINK command usage and scenario table.

### PLINK Command Usage for Center-Specific Files
The script automatically calls PLINK to create center-specific datasets:

```bash
plink --bfile <input_prefix> --keep <sample_list> --make-bed --out <output_prefix>
```

Where:
- `<input_prefix>`: Path to the full dataset (e.g., `simulated_data/tiny/tiny`)
- `<sample_list>`: Text file with sample IDs for this center (tab-separated FID/IID)
- `<output_prefix>`: Output path for center-specific files

This ensures each center receives only its assigned samples, with identical SNP sets.

---

## Output Structure

```
simulated_data/
├── tiny/
│   ├── tiny.bed
│   ├── tiny.bim
│   ├── tiny.fam
│   ├── metadata.json
│   ├── relatives_info.csv
│   ├── assigned_samples.txt
│   ├── centralized/
│   │   ├── assoc.log
│   │   ├── assoc.assoc
│   │   ├── king.kin0
│   │   ├── hardy.hwe
│   │   ├── missing.imiss
│   │   ├── ...
│   │   ├── assoc_assigned.*
│   │   ├── king_assigned.*
│   │   ├── hardy_assigned.*
│   │   └── missing_assigned.*
│   └── center_1/
│       ├── samples.txt
│       ├── id_map_center_1.csv
│       ├── metadata_center_1.json
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

---

## File Descriptions

- **`*.bed/.bim/.fam`**: PLINK-format genotype files (full dataset)
- **`metadata.json`**: Dataset info (scale, partition, causal SNPs, case/control counts, family summary)
- **`relatives_info.csv`**: Table of all simulated relatives (columns: `SampleID`, `AncestryID`, `Relatedness`)
- **`assigned_samples.txt`**: List of all FID/IID pairs assigned to any center (union, no duplicates)
- **`center_*/`**: Per-center partitioned datasets
  - **`samples.txt`**: Sample ID list for the center
  - **`id_map_center_{i}.csv`**: Mapping from original sample index to FID/IID for traceability
  - **`metadata_center_{i}.json`**: Per-center metadata (sample/case/control counts, relatives, etc.)
  - **`*_center_*.bed/.bim/.fam`**: Center-specific PLINK files

---

## Assigned Samples for Benchmarking

After partitioning, the file `assigned_samples.txt` in each dataset directory contains the union of all FID/IID pairs assigned to any center (tab-separated, no duplicates). This is essential for:
- **Fair benchmarking**: When comparing federated (partitioned) results to centralized results, use only the assigned samples.
- **Reproducibility**: Ensures you can always reconstruct the set of samples used in federated analysis.

**How to use:**

To run centralized analysis on only the assigned samples, use PLINK’s `--keep` option:

```bash
plink --bfile tiny/tiny --keep tiny/assigned_samples.txt --assoc
```

This ensures the centralized results are directly comparable to the merged federated results.

---

## Centralized Ground Truth Results

After generation, you can run with `--ground-truths` to produce a fully QCed, realistic centralized GWAS analysis for benchmarking and validation. These are saved in a `centralized/` subdirectory for each scale.

### QC and Association Pipeline
The following steps are performed automatically:
1. **Missingness and HWE filtering:**
   - Remove samples with >5% missing genotypes (`--mind 0.05`)
   - Remove SNPs with >5% missingness (`--geno 0.05`)
   - Remove SNPs with HWE p-value < 1e-6 (`--hwe 1e-6`)
2. **Kinship/relatedness estimation:**
   - Compute pairwise IBD/kinship using PLINK `--genome` on the QCed data
3. **Relatives removal:**
   - Remove one sample from each pair with PI_HAT > 0.2 (i.e., keep only unrelateds)
4. **Final association analysis:**
   - Run PLINK `--assoc` on the unrelated, QCed dataset

### Output Files
All intermediate and final files are saved in `centralized/`:
- `qc1.*`: Dataset after missingness and HWE filtering
- `qc1.genome`: Pairwise kinship/IBD estimates
- `unrelated.txt`: List of unrelated samples used for final analysis
- `qc2.*`: Final unrelated, QCed dataset
- `assoc_final.*`: Final association results (gold-standard ground truth)

### Usage

```bash
python generate_synthetic_data.py --scale tiny --ground-truths
# or for all scales:
python generate_synthetic_data.py --scale all --ground-truths
```

**Why is this important?**
- Provides a gold-standard reference for benchmarking federated vs. centralized analysis.
- Ensures that the centralized results are realistic and directly comparable to best-practice GWAS pipelines.

---

## Per-Center Sample ID Mapping

Each center directory contains an `id_map_center_{i}.csv` file, which maps the original sample index in the full dataset to the FID/IID used in that center. This mapping is crucial for:
- Traceability between the original and partitioned datasets
- Joining, auditing, or debugging results across centers
- Reproducibility in federated or distributed analysis

**Example: id_map_center_1.csv**

| original_index | FID      | IID      |
|---------------:|----------|----------|
| 0              | F000001  | I000001  |
| 5              | F000006  | I000006  |
| ...            | ...      | ...      |

- `original_index`: 0-based index in the full dataset
- `FID`/`IID`: Family and individual IDs used in PLINK files for this center

If you shuffle, anonymize, or reindex samples, this file provides the necessary mapping for all downstream analysis.

---

## Code Structure

- **generate_synthetic_data.py**: Main generator script (all logic, CLI)
- **generate_test_scenarios.py**: Script to generate all test scenarios (see below)
- **requirements_synthetic.txt**: Python dependencies ([link](requirements_synthetic.txt))
- **run_tests.py**: Test runner for the suite
- **tests/**: Contains `test_synthetic_data.py` (comprehensive pytest suite)

---

## Usage

### Install dependencies
```

### Generate all test scenarios
```bash
python generate_test_scenarios.py
```

### Custom output directory
```bash
python generate_synthetic_data.py --output-dir /path/to/output
```

### Reproducibility and Random Seed

To ensure reproducibility, you can fix the random seed for all data generation steps:

```bash
python generate_synthetic_data.py --scale tiny --seed 42
```

Using the same seed will generate identical datasets every time. The seed used is also saved in `metadata.json` for each dataset for full traceability.

---

## Test Coverage

- MAF generation (range, shape, dtype)
- Family/relatives simulation (fraction, degree, kinship)
- Genotype generation (with/without relatives)
- Missingness (rate, value)
- Binary phenotype (causal SNPs, effect sizes)
- All partitioning strategies
- PLINK file creation
- Full pipeline integration (tiny/small)
- Output file structure and metadata

Run tests:
```bash
python run_tests.py
# or
pytest tests/test_synthetic_data.py -v
```

---

## Validation

After generation, validate with PLINK:
```bash
plink --bfile simulated_data/tiny/tiny --freq
plink --bfile simulated_data/tiny/tiny --assoc
```

---

## Interpreting Relatives and Kinship

- **relatives_info.csv**: Each row is a (SampleID, AncestryID, Relatedness) triple.
  - `Relatedness` is one of: `first-degree`, `second-degree`, `third-degree`.
  - Use KING or similar tools to verify kinship coefficients:
    - First-degree: ~0.25–0.5
    - Second-degree: ~0.08–0.25
    - Third-degree: ~0–0.15
- Relatives are generated using Mendelian inheritance, so coefficients are realistic but stochastic.

---

## Extending or Customizing the Generator

- **Add new scales**: Edit the `scales` dict in `generate_synthetic_data.py`.
- **Change relative fractions**: Adjust `relative_fraction` and `relative_distribution` in the config.
- **Add new partitioning strategies**: Implement a new method in the generator and add a CLI option.
- **Change phenotype model**: Edit `generate_binary_phenotype`.
- **Add new output files**: Modify `generate_scale` and `create_plink_files`.

---

## Troubleshooting & Common Issues

- **PLINK not found**: Ensure PLINK is installed and in your PATH.
- **Disk space**: Large scale requires significant space (up to several TBs).
- **Test failures**: If tests fail due to stochasticity, rerun or adjust tolerances in the test suite.
- **Slow generation**: Use tiny/small scales for development; large/medium for benchmarking only.
- **Missing dependencies**: Install with `pip install -r requirements_synthetic.txt`.
- **Output files missing**: Check logs for errors, especially PLINK failures.

---

## References

- [requirements_synthetic.txt](requirements_synthetic.txt)
- [syn_data_instruction.md](syn_data_instruction.md)
- [PLINK documentation](https://www.cog-genomics.org/plink/)

--- 
