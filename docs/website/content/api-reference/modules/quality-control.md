---
title: Quality Control
---

# Quality Control

This module implements **federated quality control**: local sample-level missingness filtering, then encrypted variant-level summaries so each client can apply harmonized SNP filters. Flower stages are `local_qc`, `global_qc`, and `global_qc_response`.

## Local sample filtering (`local_qc`)

`exclude_samples_by_missing_rate` runs PLINK `--mind` to remove samples whose missing genotype rate exceeds the configured threshold:

```bash
plink --bfile <plink_prefix> --mind <threshold> --make-bed --out <new_prefix>
```

The output prefix is written under the client's configured `output.log_dir`, and the client updates its active `plink_prefix`.

## Variant filtering (`global_qc`)

During `global_qc`, each client computes:

- Genotype counts and HWE values with `compute_genotype_counts`, using `plink --freq` and `plink --hardy`.
- SNP missingness counts with `compute_missingness_counts`, using `plink --missing`.
- Threshold arrays containing MAF, missingness, and HWE thresholds.

Clients pack these arrays as typed payloads and encrypt them for peers. The server forwards the ciphertext arrays during `global_qc`; it does not compute the exclusion list.

During `global_qc_response`, each client decrypts peer payloads and `_compute_exclusion_list` then:

1. Sums all genotype and missingness arrays.
2. Derives unified thresholds from all clients.
3. Computes MAF, missingness rate, and HWE chi-square p-values for each SNP.
4. Returns the SNP IDs that should be excluded locally.

## SNP filtering

Clients apply the returned exclusion list with `exclude_snps`:

```bash
plink --bfile <plink_prefix> --exclude <exclude_file> --make-bed --out <new_prefix>
```

If the exclusion list is empty, the original `plink_prefix` is kept.
