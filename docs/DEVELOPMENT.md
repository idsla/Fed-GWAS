# Federated GWAS Pipeline

[![PyPI](https://img.shields.io/pypi/v/FedGWAS.svg)](https://pypi.org/project/FedGWAS/)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-2f8f83)](https://idsla.github.io/Fed-GWAS/)
[![Deploy documentation](https://github.com/idsla/Fed-GWAS/actions/workflows/deploy-docs.yml/badge.svg)](https://github.com/idsla/Fed-GWAS/actions/workflows/deploy-docs.yml)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

This repository implements a federated pipeline for Genome-Wide Association Studies (GWAS) using Flower, PLINK, and custom privacy-preserving protocols. The pipeline supports multi-stage, multi-client GWAS with reproducible outputs and structured logging.

For release verification steps, see [RELEASE.md](RELEASE.md). For implementation details and change history, see [CURRENT_VERSION.md](CURRENT_VERSION.md).

---

## Version Updates

Prepare version bumps on `main` before creating the release tag or publishing a GitHub Release. A normal push to `main` does not publish the package and does not deploy the documentation site.

When bumping a version on `main`:

1. Update the package version in `pyproject.toml` (`[project].version`). Use the plain PEP 440 version, for example `0.0.1`; reserve the leading `v` for Git tags such as `v0.0.1`.
2. Update release-facing docs that mention the current version or install examples:
   - `CURRENT_VERSION.md`
   - `RELEASE.md`
   - `README.md` and `docs/` pages, when user-facing commands or versioned examples changed
3. Keep `CURRENT_VERSION.md` aligned with the code shipped by the version bump. Add the date, release headline, highlights, behavior changes, and any migration notes needed by users.
4. Run the release verification path from `RELEASE.md` before tagging the version bump on `main`. At minimum, run the focused tests and package/docs checks that match the changed surface:

```bash
uv sync --python 3.11
pytest tests/test_king_federated_unit.py -q
python -m pip install --upgrade build twine
python -m build
python -m twine check dist/*
cd docs/website
npm ci
npm run build
```

The GitHub Actions release behavior is:

- `.github/workflows/deploy-docs.yml` deploys GitHub Pages only on tag push (`*`) or manual `workflow_dispatch`.
- `.github/workflows/publish-pypi.yml` builds and publishes to PyPI only when a GitHub Release is published.

After the version bump is on `main`, create the tag on that exact release commit and push it:

```bash
git tag v0.0.1
git push origin v0.0.1
```

Then create and publish the GitHub Release for that tag. The tag push deploys the documentation site, and publishing the GitHub Release triggers the PyPI Trusted Publishing workflow. If a published PyPI version needs a fix, bump to a new version on `main`; PyPI versions are immutable.

---

## Environment Setup

### Option 1: UV (recommended)

Install [uv](https://docs.astral.sh/uv/):

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Sync dependencies (Python 3.11+):

```bash
uv sync --python 3.11
```

Optional dev dependencies:

```bash
uv sync --dev
```

### Option 2: Conda

```bash
conda create -n fedgwas python=3.11 -y
conda activate fedgwas
pip install -e .
pip install -U "flwr[simulation]"
```

---

## PLINK

- Requires [PLINK 1.9+](https://www.cog-genomics.org/plink/1.9/).
- Download the binary for your OS and ensure `plink` is on your `PATH`, or set the path in each client `config.yaml` (`plink.path` if configured).
- Toy reference files are under `plink/`; production runs use experiment data under `experiments/`.

---

## Quick Start (Recommended: tiny_even)

The default Flower config in `pyproject.toml` points to `experiments/correctness/tiny_even/configs` (2 clients, tiny synthetic data).

### Repository layout (experiments)

```
experiments/correctness/tiny_even/
├── config.yaml
├── configs/
│   ├── server/config.yaml
│   ├── center_1/config.yaml
│   └── center_2/config.yaml
├── data/tiny/
│   ├── center_1/          # PLINK .bed/.bim/.fam per client
│   ├── center_2/
│   └── centralized_baseline/   # after generate_baseline
└── results_2/             # gitignored; current shipped config output
```

Config templates: [configs/config_template.yaml](configs/config_template.yaml).

### 1. Generate synthetic data (if not present)

```bash
python pipeline/simulation/simulated_data/generate_synthetic_data.py \
  --scale tiny \
  --partition-strategy even \
  --seed 42 \
  --output-dir experiments/correctness/tiny_even/data
```

### 2. Generate centralized baseline

```bash
python experiments/tools/generate_baseline.py \
  experiments/correctness/tiny_even/config.yaml
```

### 3. Run federated pipeline (simulation)

```bash
flwr run . local-simulation --stream
```

Override rounds or config path:

```bash
flwr run . local-simulation --stream --run-config \
  'simulation=true num-server-rounds=100 config_path="experiments/correctness/tiny_even/configs"'
```

Results are written under each client's `logs/` and `intermediate/` directories (paths set in per-center `config.yaml`). The shipped tiny configs currently write under `experiments/correctness/tiny_even/results_2/`; use the paths in the active center and server config files as the source of truth.

### 4. Retention (optional, automatic)

Experiment `config.yaml` may set `retention.tier` (`minimal` | `standard` | `research`). When `auto_apply_on_complete: true`, the server prunes non-essential artifacts after the run. Manual:

```bash
python experiments/tools/apply_run_retention.py \
  experiments/correctness/tiny_even/results \
  --config-path experiments/correctness/tiny_even/configs \
  --dry-run
```

See [RELEASE.md](RELEASE.md) for tier definitions.

### 5. Evaluate against baseline

```bash
python experiments/tools/evaluation/evaluate_all.py \
  experiments/correctness/tiny_even/results_2 \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline \
  --king
```

See [experiments/correctness/tiny_even/README.md](experiments/correctness/tiny_even/README.md) for expected metrics and success criteria. If you changed the output paths in the active configs, pass that results directory instead.

---

## Documentation Site

The documentation website is a Docusaurus project isolated under `docs/website/`. This keeps Node.js dependencies, build output, Docusaurus configuration, and website content grouped under the repository documentation tree.

The site reads content from three places:

- `docs/website/content/get-started/`: quickstart and onboarding pages under the `Get Started` navigation tab.
- `docs/website/content/user-guide/`: user guide pages under the `User Guide` navigation tab.
- `docs/website/content/examples/`: example workflows under the `Examples` navigation tab.
- `docs/website/content/api-reference/`: API and output reference pages under the `API Reference` navigation tab.

Generated folders are local-only and should not be committed:

- `docs/website/node_modules/`
- `docs/website/.docusaurus/`
- `docs/website/build/`

### Prerequisites

Install Node.js and npm before working on the documentation site. The GitHub Pages workflow uses Node.js 20. Local Node.js 18 also works with the pinned Docusaurus 3.8.1 dependencies, although npm may print engine warnings for transitive packages.

Check your local versions:

```bash
node --version
npm --version
```

### Install website dependencies

Run npm commands from the `docs/website/` directory, not from the repository root:

```bash
cd docs/website
npm install
```

Use `npm ci` when you want to reproduce the exact dependency tree from `package-lock.json`, especially before checking a release or debugging CI:

```bash
cd docs/website
npm ci
```

### Start the local development server

Start Docusaurus in development mode:

```bash
cd docs/website
npm run start
```

The local URL uses the configured GitHub Pages `baseUrl`:

```text
http://localhost:3000/Fed-GWAS/
```

If port `3000` is already in use, pass another port:

```bash
cd docs/website
npm run start -- --port 3001
```

Development mode watches files and refreshes the browser when you edit Markdown, sidebars, or React/CSS files.

### Build the production site

Before opening a pull request or publishing documentation, run a production build:

```bash
cd docs/website
npm run build
```

This writes static files to `docs/website/build/`. A successful build also catches broken Docusaurus routes, invalid MDX, and many broken internal links.

Preview the built site locally:

```bash
cd docs/website
npm run serve
```

The served production build is available at:

```text
http://localhost:3000/Fed-GWAS/
```

### Documentation navigation

The navbar is configured in `docs/website/docusaurus.config.js`:

- `Get Started`: onboarding docs from `docs/website/content/get-started/`.
- `User Guide`: operational docs from `docs/website/content/user-guide/`.
- `Examples`: example workflows from `docs/website/content/examples/`.
- `API Reference`: API and output references from `docs/website/content/api-reference/`.
- Version selector: displayed on the right side of the navbar.

Sidebars are split by section:

- `docs/website/sidebarsGetStarted.js`
- `docs/website/sidebarsUserGuide.js`
- `docs/website/sidebarsExamples.js`
- `docs/website/sidebarsReference.js`

Update the matching sidebar file whenever you add, remove, or rename a page.

### Search

The site uses local search through `@easyops-cn/docusaurus-search-local`. The search index is generated during `npm run build` and includes:

- `docs/`
- `docs/website/content/get-started/`
- `docs/website/content/user-guide/`
- `docs/website/content/examples/`
- `docs/website/content/api-reference/`

If search results look stale, remove generated files and rebuild:

```bash
cd docs/website
npm run clear
npm run build
```

### Troubleshooting

- If `npm run build` reports invalid document IDs, check that the page exists and is listed in the correct sidebar file.
- If a Markdown file fails MDX compilation, check fenced code blocks and unescaped JSX-like syntax such as raw `<tag>` text.
- If the local site opens at `/` but assets fail to load, use `/Fed-GWAS/` in the URL because `baseUrl` is configured for GitHub Pages.
- If npm reports audit warnings, treat them separately from build correctness. Do not run `npm audit fix --force` without checking whether it upgrades Docusaurus or React across major versions.

---

## Three-Node Cluster Deployment

For Matpool or any 3-node layout (1 SuperLink + 2 SuperNodes), use the bundled scripts and guide:

- **Guide:** [cluster_deployment/docs/CLUSTER_USER_GUIDE.md](cluster_deployment/docs/CLUSTER_USER_GUIDE.md)
- **Scripts:** [cluster_deployment/README.md](cluster_deployment/README.md)

```bash
bash cluster_deployment/scripts/setup-cluster-node.sh   # each node
bash cluster_deployment/scripts/cluster-verify-data.sh --scale tiny --client-id 1  # each client
cluster_deployment/scripts/cluster-run-app.sh \
  --server-ip <SERVER_IP> --scale tiny --rounds 20
```

Performance scales (small/medium): `experiments/performance/scales.yaml` and per-scale READMEs under `small_even/`, `medium_even/`.

---

## Local Deployment Mode

Requires SuperLink + two SuperNodes + `flwr run`:

```bash
flower-superlink --insecure
```

```bash
flower-supernode --insecure --superlink 127.0.0.1:9092 --clientappio-api-address 127.0.0.1:9094 \
  --node-config 'partition-id=0 num-partitions=2 config-file="experiments/correctness/tiny_even/configs/center_1/config.yaml"'
```

```bash
flower-supernode --insecure --superlink 127.0.0.1:9092 --clientappio-api-address 127.0.0.1:9095 \
  --node-config 'partition-id=1 num-partitions=2 config-file="experiments/correctness/tiny_even/configs/center_2/config.yaml"'
```

```bash
flwr run . local-deployment --stream
```

---

## Advanced: Real-World Experiments

Larger studies (e.g. 1000 Genomes subset) live under `experiments/real_world/1000genomes/`. These require downloading/preparing data, longer runtime, and overriding `config_path`:

```bash
flwr run . local-simulation --stream --run-config \
  'config_path="experiments/real_world/1000genomes/configs"'
```

Manuscript figures and prior run outputs under `experiments/real_world/1000genomes/manuscript/` are research artifacts and are not required for the default release path.

---

## Output and Logs

- Per-client `intermediate_dir` and `log_dir` are defined in each center `config.yaml`.
- Directories are cleared at the start of each client run to avoid stale artifacts.
- Stage progress and errors go to per-client log files under each configured `output.log_dir`.
- Inspect PLINK outputs (`.assoc.logistic`, `.imiss`, `.frq`, KING kinship files) directly under each client's `logs/`.

---

## Federated Protocol (Summary)

1. **Key exchange** — ECC public keys via server relay  
2. **Sync** — Encrypted seed broadcast (server cannot decrypt)  
3. **Local / global QC** — Encrypted QC shares; exclusion list computed client-side  
4. **Iterative KING** — Chunked kinship with cross-client anonymized IDs  
5. **Local LR + filtering** — Tokenized insignificant SNPs  
6. **Iterative LR** — Chunked association on filtered data  

Full stage contracts and privacy model: [CURRENT_VERSION.md](CURRENT_VERSION.md).

---

## Troubleshooting

- **PLINK not found** — Install PLINK 1.9+ and verify `plink` is on `PATH` or configured in `config.yaml`.  
- **Wrong config** — Check `config_path` in `pyproject.toml` or pass `--run-config`.  
- **Empty results** — Ensure data and baseline exist under `experiments/correctness/tiny_even/data/`.  
- **Reproducibility** — Use fixed seeds in data generation and consistent `config_path` across runs.

---

## Contributing

Open issues or pull requests for bug fixes, improvements, or new features.

## Acknowledgments

Built with [Flower](https://flower.dev/), [PLINK](https://www.cog-genomics.org/plink/1.9/), and open-source Python tools.
