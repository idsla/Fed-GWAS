---
title: Examples Overview
---

# Examples Overview

Use these examples as entry points for correctness checks, scalability experiments, and real-world genotype screening runs.

- [Tiny correctness](tiny-correctness): Validate federated QC, relatedness screening, and association screening against a centralized baseline using synthetic tiny PLINK data with 500 samples and 5,000 SNPs.
- [Performance small](performance-small): Measure runtime and memory at small scale using synthetic small PLINK data with 2,000 samples and 20,000 SNPs.
- [1000 Genomes](1000genomes): Run the application workflow on prepared 1000 Genomes chr22 genotypes with generated phenotypes.
- [Three-node deployment](three-node-deployment): Run a server plus two clients with either script deployment (`cluster_deployment/`) or the CLI (`fedgwas-deploy`).
