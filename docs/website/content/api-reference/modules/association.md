---
title: Association Screening
---

# Association Screening

This module implements the **federated association screening** stage: local logistic-regression filtering with privacy-preserving tokens, then federated case–control logistic regression.

## Local logistic-regression filtering

During `local_lr`, each client runs:

```bash
plink --bfile <plink_prefix> --logistic --out <log_dir>/local_lr_temp
```

`parse_insignificant_snps` reads the resulting `.assoc.logistic` file and selects SNPs with `P >= local_lr_threshold` when configured, falling back to `p_threshold`. The client sends privacy-preserving (salted) tokens for those SNPs instead of raw SNP IDs.

The server forwards token arrays back to clients. During `local_lr_filter_response`, clients compute the intersection locally, resolve local SNP IDs, and remove that intersection from local data with `exclude_snps`.

## Federated logistic regression

After filtering, the workflow enters `init_chunks_lr` and `iterative_lr`:

1. Clients repartition filtered data, usually by SNP.
2. Clients send one encoded PLINK chunk per stage call.
3. `run_server_lr` reconstructs and merges chunks.
4. The server runs `plink --logistic`.
5. The server parses `.assoc.logistic` and returns SNP p-values.
6. Clients use persisted SNP map files to map anonymous SNP IDs back to local IDs when possible.
7. Clients classify final p-values with `global_lr_threshold` when configured.

## Current scope

The current implementation focuses on logistic regression for binary phenotypes. Continuous phenotypes, richer covariate handling, population stratification, or alternate association tests should be added through the PLINK command construction, configuration schema, and result parsers.
