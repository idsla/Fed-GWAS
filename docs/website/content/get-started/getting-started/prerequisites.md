---
title: Prerequisites
slug: /prerequisites
---

FedGWAS requires Python 3.11+ and PLINK.

## Python

FedGWAS requires Python 3.11+. We recommend using [uv](https://uv.io/) or [conda](https://docs.conda.io/en/latest/) to manage the Python environment and dependencies.

You can install Python from the official Python [website](https://www.python.org/downloads/), or using a package manager like `apt` on Debian-based systems, `brew` on macOS, or `chocolatey` on Windows.

## PLINK

The pipeline calls PLINK for binary dataset creation, missingness checks, genotype counts, KING-based relatedness screening, and logistic-regression association screening.

Input genotypes use **PLINK 1.9** binary format. KING kinship estimation uses **PLINK 2** `--make-king-table`. VCF input is supported via conversion to PLINK-compatible representations where configured.

1. Install [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/) so that `plink` is on `PATH`, or use the bundled project-local binaries under `plink/plink_linux`, `plink/plink_mac`, or `plink/plink_win`.
2. For KING-based relatedness screening, keep [PLINK 2](https://www.cog-genomics.org/plink/2.0/) available as `plink2` in the matching bundled directory or on `PATH`.
3. Verify that client (center) configuration values use the PLINK prefix without `.bed`, `.bim`, or `.fam` extensions.

Remember to add the PLINK executables to your system's PATH variable if you choose to install them globally.

Instructions on how to modify system PATH variables:
- windows: https://www.architectryan.com/2018/03/17/add-to-the-path-on-windows-10/
- mac or linux: https://www.digitalocean.com/community/tutorials/how-to-view-and-update-the-linux-path-environment-variable


## Quick Installation

FedGWAS is available on PyPI. Install it in a clean Python 3.11+ environment, and make sure `plink` and `plink2` are discoverable on `PATH`:

```bash
python -m venv .venv
source .venv/bin/activate # or windows powershell .\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install FedGWAS
```

verify the installation:

```bash
python --version
fedgwas-sim --help
```
