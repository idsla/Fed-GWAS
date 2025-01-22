# Quality Control Module

## Overview

The **Quality Control (QC)** module ensures that the genomic and phenotype data used for genome-wide association studies (GWAS) is clean, reliable, and ready for analysis. This module includes functions for checking missing data, calculating minor allele frequency (MAF), filtering data based on thresholds, and more.

Each function in this module processes data locally and produces specific outputs required for downstream analysis.



## Detailed User Guide

### 1. `calculate_missing_rate`

**Description**:  
This function calculates the missing genotype rate for each SNP and for each individual. It identifies missing values in the genotype data to determine call rates.

**Parameters**:

- `genotype_data (pd.DataFrame)`: The DataFrame containing genotype data where rows represent individuals and columns represent SNPs.
- `bed`: A `SnpReader` object that provides individual IDs (`iid`) and SNP IDs (`sid`).
- `output_prefix (str)`: A prefix for the output files where missing rate information will be saved.

**Returns**:

- `pd.DataFrame`: Two DataFrames, one for SNP missing rates (`.lmiss`) and one for individual missing rates (`.imiss`).

**Usage Example**:

```python
from pysnptools.snpreader import Bed
bed = Bed('path_to_bed_file')
qc = QualityControl()

snp_missing, ind_missing = qc.calculate_missing_rate(genotype_data, bed, "output_prefix")
print(snp_missing.head())
```

---

### 2. `hardy_weinberg_test`

**Description**:  
This function performs the Hardy-Weinberg Equilibrium (HWE) test on each SNP and filters out those that significantly deviate from HWE based on the given threshold.

**Parameters**:

- `genotype_data (pd.DataFrame)`: The DataFrame containing genotype data.
- `threshold (float)`: The p-value threshold for filtering SNPs based on the Hardy-Weinberg test.

**Returns**:

- `pd.DataFrame`: The filtered genotype data and a list of indices of SNPs that pass the test.

**Usage Example**:

```python
qc = QualityControl()

filtered_geno, snp_indices_to_keep = qc.hardy_weinberg_test(genotype_data, 1e-6)
print(filtered_geno.head())
```

---

### 3. `filter_missingness_samples`

**Description**:  
This function filters individuals based on the proportion of missing genotype data. Individuals with a missing genotype rate higher than the specified threshold are filtered out.

**Parameters**:

- `genotype_data (pd.DataFrame)`: The DataFrame containing genotype data.
- `fam (pd.DataFrame)`: The FAM file data.
- `output_prefix (str)`: A prefix for saving the output FAM file.

**Returns**:

- `pd.DataFrame`: A DataFrame with filtered FAM data.

**Usage Example**:

```python
qc = QualityControl()

filtered_fam = qc.filter_missingness_samples(genotype_data, fam, "output_prefix")
print(filtered_fam.head())
```

---

### 4. `geno`

**Description**:  
This function filters SNPs based on their missing genotype data rate. SNPs that exceed the threshold for missing data are removed, and the results are saved in a BIM file.

**Parameters**:

- `genotypes (pd.DataFrame)`: The genotype data matrix.
- `bim (pd.DataFrame)`: The BIM file data.
- `output_prefix (str)`: A prefix for the output file name.

**Returns**:

- `pd.DataFrame`: The filtered BIM data.

**Usage Example**:

```python
qc = QualityControl()

filtered_bim = qc.geno(genotype_data, bim, "output_prefix")
print(filtered_bim.head())
'''