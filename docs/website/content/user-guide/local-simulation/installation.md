---
title: Installation
slug: /local-simulation/installation
---

# Installation

Local simulation runs the FedGWAS pipeline on one workstation with Flower local simulation mode. It is the recommended first workflow for checking an installation, generating synthetic center data, validating configuration, and comparing federated outputs against a centralized baseline.

## Requirements

Use a dedicated Python environment. FedGWAS currently targets:

| Dependency | Requirement |
| --- | --- |
| Python | 3.11 or later |
| FedGWAS | Installed package or editable repository install |
| Flower | 1.19.x, installed through the FedGWAS dependency set |
| PLINK | PLINK 1.9 (`plink`) on `PATH` or from the repository bundle; PLINK 2 (`plink2`) for KING-based relatedness screening |

PLINK 1.9 is required for federated quality control and logistic-regression association screening. KING-based relatedness screening uses PLINK 2 `--make-king-table`. The simulation CLI checks for PLINK before starting a run.

## Install PLINK

Install PLINK 1.9 so that `plink` is on `PATH`, and PLINK 2 as `plink2` for relatedness screening.
Confirm the installation with:

```bash
plink --version
plink2 --version
```

If PLINK is not on `PATH`, local simulation projects can still set `plink` in `fedgwas.yaml`, but the most portable setup is to expose `plink` or `plink2` in the active environment.

## Install From PyPI

For most users, install FedGWAS from PyPI in a clean environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install FedGWAS
```

On Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install FedGWAS
```

Confirm that the simulation CLI is available:

```bash
fedgwas-sim --help
```

**Verify The Runtime**

Run these checks before creating a study:

```bash
python --version
python -c "import pipeline; print('FedGWAS import OK')"
python -c "import flwr; print(flwr.__version__)"
fedgwas-sim --help
plink --version
```

The Flower version should be in the 1.19 line. If it is not, reinstall the supported dependency:

```bash
python -m pip install --upgrade "flwr[simulation]>=1.19.0,<1.20"
```

## [Optional] Install from Repository

You can also install FedGWAS directly from the repository, and run simulations with `flwr run . local-simulation` and local scripts in the repo.

```bash
git clone https://github.com/idsla/Fed-GWAS.git
cd Fed-GWAS
python -m pip install -e .
```

With `uv`:

```bash
git clone https://github.com/idsla/Fed-GWAS.git
cd Fed-GWAS
uv sync --python 3.11
```

The repository includes bundled PLINK binaries under `plink/plink_linux`, `plink/plink_mac`, and `plink/plink_win`. You can add the appropriate directory to your `PATH` or set the PLINK path in `fedgwas.yaml` to use these binaries. And repo contains `experiments` folder with example configurations and data for local simulation runs.

**Difference between PyPI and Repo InstallationL**

The repository workflow and the installed-package workflow use the same core pipeline, but their layouts differ:

| Workflow | Primary command | Config layout |
| --- | --- | --- |
| Installed package | `fedgwas-sim ...` | User study directory with `configs/server.yaml` and `configs/center_*.yaml` |
| Repository checkout | `flwr run . local-simulation` | Repository `experiments/.../configs/server/config.yaml` and `center_*/config.yaml` |

Use one layout consistently within a run.
