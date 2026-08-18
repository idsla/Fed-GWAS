#!/usr/bin/env bash
set -eo pipefail
PS1="${PS1-}"

# Prefer fedgwas conda env on Matpool nodes
for _bin in \
  "${FEDGWAS_BIN:-}" \
  "${CONDA_PREFIX:-}/bin" \
  "/root/miniconda3/envs/fedgwas/bin" \
  "${HOME}/miniconda3/envs/fedgwas/bin"; do
  if [[ -n "${_bin}" && -x "${_bin}/flwr" ]]; then
    export PATH="${_bin}:${PATH}"
    break
  fi
done
unset _bin

# Run the Flower app on the server node (Flower 1.19.x + pyproject.toml federations).
# Usage:
#   ./scripts/cluster-run-app.sh --server-ip <server_ip> [--rounds N] [--scale tiny|small|medium]
#   ./scripts/cluster-run-app.sh --server-ip <ip> --config-path experiments/performance/small_even/configs

SERVER_IP=""
ROUNDS="15"
SCALE="tiny"
CONFIG_PATH=""
PYPROJECT="pyproject.toml"
BACKUP=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --server-ip)
      SERVER_IP="$2"
      shift 2
      ;;
    --rounds)
      ROUNDS="$2"
      shift 2
      ;;
    --scale)
      SCALE="$2"
      shift 2
      ;;
    --config-path)
      CONFIG_PATH="$2"
      shift 2
      ;;
    *)
      echo "Unknown arg: $1" >&2
      echo "Usage: $0 --server-ip <ip> [--rounds N] [--scale tiny|small|medium] [--config-path <dir>]" >&2
      exit 1
      ;;
  esac
done

if [[ -z "${SERVER_IP}" ]]; then
  echo "Missing required arg: --server-ip <ip>" >&2
  exit 1
fi

if [[ -z "${CONFIG_PATH}" ]]; then
  case "${SCALE}" in
    tiny)
      CONFIG_PATH="experiments/correctness/tiny_even/configs"
      ;;
    small)
      CONFIG_PATH="experiments/performance/small_even/configs"
      ;;
    medium)
      CONFIG_PATH="experiments/performance/medium_even/configs"
      ;;
    *)
      echo "Unknown scale: ${SCALE} (use tiny, small, or medium)" >&2
      exit 1
      ;;
  esac
fi

if [[ ! -d "${CONFIG_PATH}" ]]; then
  echo "Config path not found: ${CONFIG_PATH}" >&2
  exit 1
fi

if [[ ! -f "${PYPROJECT}" ]]; then
  echo "pyproject.toml not found in CWD" >&2
  exit 1
fi

if ! command -v flwr &>/dev/null; then
  echo "Error: flwr not found. export PATH=/root/miniconda3/envs/fedgwas/bin:\$PATH" >&2
  exit 127
fi

FLWR_VER="$(flwr --version 2>&1 || true)"
if [[ "${FLWR_VER}" != *1.19.* ]]; then
  echo "Warning: expected Flower 1.19.x; got: ${FLWR_VER}" >&2
  echo "  pip install 'flwr[simulation]>=1.19.0,<1.20'" >&2
fi

BACKUP="${PYPROJECT}.bak.cluster"
cp "${PYPROJECT}" "${BACKUP}"

restore_pyproject() {
  if [[ -n "${BACKUP}" && -f "${BACKUP}" ]]; then
    mv -f "${BACKUP}" "${PYPROJECT}"
  fi
}
trap restore_pyproject EXIT

echo "[server] Setting Exec API address to ${SERVER_IP}:9093 in pyproject.toml"
python3 - <<PYCODE
import pathlib
import re

addr = "${SERVER_IP}:9093"
path = pathlib.Path("${PYPROJECT}")
text = path.read_text()

for section in ("local-deployment", "cluster-deployment"):
    pattern = (
        rf'(\[tool\.flwr\.federations\.{re.escape(section)}\][^\[]*?'
        rf'address = ")[^"]+(")'
    )
    text, n = re.subn(
        pattern,
        lambda m, a=addr: f"{m.group(1)}{a}{m.group(2)}",
        text,
        count=1,
        flags=re.DOTALL,
    )
    if n == 0:
        raise SystemExit(f"Could not patch address in [tool.flwr.federations.{section}]")

path.write_text(text)
PYCODE

echo "[server] Running flwr app (local-deployment) config=${CONFIG_PATH} for ${ROUNDS} rounds"
flwr run . local-deployment --stream --run-config \
  "simulation=false num-server-rounds=${ROUNDS} config_path=\"${CONFIG_PATH}\""

echo "[server] Restoring pyproject.toml"
restore_pyproject
trap - EXIT
