#!/usr/bin/env bash
# Fed-GWAS cluster node setup: PLINK, Python env, cluster script deps, smoke checks.
# Run from anywhere; resolves repo root automatically.
# Non-interactive / Jupyter shells often lack PS1; avoid nounset (set -u) breaking conda.
set -eo pipefail
PS1="${PS1-}"

echo "=========================================="
echo "Fed-GWAS Cluster Node Setup"
echo "=========================================="

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${SCRIPT_DIR}"
while [[ "${REPO_ROOT}" != "/" && ! ( -f "${REPO_ROOT}/pyproject.toml" && -d "${REPO_ROOT}/pipeline" ) ]]; do
  REPO_ROOT="$(cd "${REPO_ROOT}/.." && pwd)"
done
if [[ ! -f "${REPO_ROOT}/pyproject.toml" || ! -d "${REPO_ROOT}/pipeline" ]]; then
  echo "Error: could not find FedGWAS repository root from ${SCRIPT_DIR}" >&2
  exit 1
fi
cd "${REPO_ROOT}"
echo "Repository root: ${REPO_ROOT}"

if [[ "${EUID}" -eq 0 ]]; then
  SUDO=""
else
  SUDO="sudo"
fi

# --- [1/5] PLINK ---
if ! command -v plink &>/dev/null; then
  echo "[1/5] Installing PLINK from local binaries..."

  OS_TYPE="$(uname -s)"
  case "${OS_TYPE}" in
    Linux*) PLINK_SUBDIR="plink_linux" ;;
    Darwin*) PLINK_SUBDIR="plink_mac" ;;
    MINGW*|MSYS*|CYGWIN*) PLINK_SUBDIR="plink_win" ;;
    *) echo "Warning: unknown OS ${OS_TYPE}, trying Linux binary..."; PLINK_SUBDIR="plink_linux" ;;
  esac

  PLINK_ROOT="${REPO_ROOT}/plink/${PLINK_SUBDIR}"

  if [[ -f "${PLINK_ROOT}/plink" || -f "${PLINK_ROOT}/plink.exe" ]]; then
    if [[ -f "${PLINK_ROOT}/plink" ]]; then
      ${SUDO} cp "${PLINK_ROOT}/plink" /usr/local/bin/plink
    else
      ${SUDO} cp "${PLINK_ROOT}/plink.exe" /usr/local/bin/plink
    fi
    ${SUDO} chmod +x /usr/local/bin/plink || true

    if [[ -f "${PLINK_ROOT}/plink2" ]]; then
      ${SUDO} cp "${PLINK_ROOT}/plink2" /usr/local/bin/plink2
      ${SUDO} chmod +x /usr/local/bin/plink2 || true
    elif [[ -f "${PLINK_ROOT}/plink2.exe" ]]; then
      ${SUDO} cp "${PLINK_ROOT}/plink2.exe" /usr/local/bin/plink2
      ${SUDO} chmod +x /usr/local/bin/plink2 || true
    fi
    echo "✓ PLINK installed from ${PLINK_ROOT}"
  else
    echo "Warning: no binary in ${PLINK_ROOT}; downloading PLINK 1.9 (Linux x86_64)..."
    if ! command -v wget &>/dev/null || ! command -v unzip &>/dev/null; then
      echo "Error: install wget and unzip, or add plink binaries under plink/plink_linux/" >&2
      exit 1
    fi
    TMP_DIR="$(mktemp -d)"
    wget -q -O "${TMP_DIR}/plink.zip" \
      https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip
    unzip -q "${TMP_DIR}/plink.zip" -d "${TMP_DIR}"
    ${SUDO} mv "${TMP_DIR}/plink" /usr/local/bin/plink
    ${SUDO} chmod +x /usr/local/bin/plink
    rm -rf "${TMP_DIR}"
    echo "✓ PLINK installed from download"
  fi
else
  echo "[1/5] PLINK already installed: $(plink --version 2>&1 | head -1)"
fi

# --- [2/5] GNU time (benchmark script) ---
if [[ -x /usr/bin/time ]]; then
  if /usr/bin/time -v true &>/dev/null; then
    echo "[2/5] GNU time OK (/usr/bin/time -v)"
  else
    echo "[2/5] Warning: /usr/bin/time exists but -v failed; install package 'time' (Ubuntu)"
  fi
else
  echo "[2/5] Warning: /usr/bin/time not found — cluster-run-benchmark.sh needs it"
  echo "       Ubuntu: sudo apt-get install -y time"
fi

# --- [3/5] Python environment ---
PYTHON="python"
PIP="pip"

if command -v conda &>/dev/null; then
  USE_CONDA=true
  echo "[3/5] Using conda (no conda.sh — safe for Jupyter / set -u shells)"
  # Prefer the real conda executable, not a shell function that may source conda.sh
  CONDA_EXE="$(type -P conda 2>/dev/null || true)"
  if [[ -z "${CONDA_EXE}" && -x "${HOME}/miniconda3/bin/conda" ]]; then
    CONDA_EXE="${HOME}/miniconda3/bin/conda"
  elif [[ -z "${CONDA_EXE}" && -x "/root/miniconda3/bin/conda" ]]; then
    CONDA_EXE="/root/miniconda3/bin/conda"
  fi
  if [[ -z "${CONDA_EXE}" ]]; then
    echo "Error: conda not found on PATH" >&2
    exit 1
  fi
  CONDA_BASE="$("${CONDA_EXE}" info --base)"
  FEDGWAS_ENV="${CONDA_BASE}/envs/fedgwas"
  if [[ ! -d "${FEDGWAS_ENV}" ]]; then
    echo "      Creating conda env 'fedgwas' (Python 3.11)..."
    "${CONDA_EXE}" create -n fedgwas python=3.11 -y
  else
    echo "      Conda env 'fedgwas' already exists"
  fi
  export PATH="${FEDGWAS_ENV}/bin:${PATH}"
  PYTHON="${FEDGWAS_ENV}/bin/python"
  PIP="${FEDGWAS_ENV}/bin/pip"
  echo "      Using ${PYTHON}"
else
  USE_CONDA=false
  echo "[3/5] Conda not found; using venv/uv"
fi

if [[ "${USE_CONDA}" != true ]]; then
  if command -v uv &>/dev/null; then
    [[ -d .venv ]] || uv venv --python 3.11
    # shellcheck disable=SC1091
    source .venv/bin/activate
    PYTHON="python"
    PIP="pip"
  else
    [[ -d .venv ]] || python3.11 -m venv .venv
    # shellcheck disable=SC1091
    source .venv/bin/activate
    PYTHON="python"
    PIP="pip"
  fi
fi

# --- [4/5] Python packages (pip wheels only — recommended on Matpool) ---
echo "[4/5] Installing Python dependencies (pip binary wheels)..."
export PIP_PREFER_BINARY=1
"${PIP}" install --upgrade "pip>=24" setuptools wheel
"${PIP}" install --prefer-binary numpy pandas scipy pyarrow

if [[ ! -f pyproject.toml ]]; then
  echo "Error: pyproject.toml not found in ${REPO_ROOT}" >&2
  exit 1
fi

# flwr before editable install (libcst needs pip>=24 + prebuilt wheel)
"${PIP}" install --prefer-binary "flwr[simulation]>=1.19.0,<1.20" tomli tomli-w psutil
"${PIP}" install --prefer-binary -e .

# --- [5/5] Scripts + verify ---
echo "[5/5] Cluster scripts and smoke checks..."
if [[ -d cluster_deployment/scripts ]]; then
  chmod +x cluster_deployment/scripts/cluster-*.sh 2>/dev/null || true
  chmod +x cluster_deployment/scripts/setup-cluster-node.sh 2>/dev/null || true
  echo "✓ cluster_deployment/scripts/*.sh executable"
fi

"${PYTHON}" - <<'PY'
import sys
for mod in ("pipeline", "flwr", "pandas", "scipy", "psutil", "tomli"):
    __import__(mod)
print("✓ fedgwas (pipeline) + deps OK")
PY

echo ""
echo "=========================================="
echo "Setup complete"
echo "=========================================="
echo ""
echo "Activate:"
if [[ "${USE_CONDA}" == true ]]; then
  echo "  export PATH=\"${FEDGWAS_ENV}/bin:\$PATH\""
  echo "  # or: conda activate fedgwas"
else
  echo "  source ${REPO_ROOT}/.venv/bin/activate"
fi
echo ""
echo "Verify:"
echo "  pip show fedgwas"
echo "  python -c \"import pipeline; import flwr\""
echo "  plink --version"
echo "  /usr/bin/time -v true 2>&1 | head -3"
echo ""
echo "Cluster post-setup (see cluster_deployment/docs/CLUSTER_SETUP.md):"
echo "  • Clients: set flower.server_address to <SERVER_IP>:9092 in configs/center_*/config.yaml"
echo "  • Clients: copy tiny_even .bed/.bim/.fam to data/tiny/center_X/"
echo "  • All nodes: align output paths under experiments/correctness/tiny_even/results_2/"
echo "  • Firewall: open ports 9092-9095 between nodes"
echo ""
