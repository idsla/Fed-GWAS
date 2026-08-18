#!/usr/bin/env bash
# Start Flower SuperNode on a client node (foreground — use nohup ... & for background).
set -eo pipefail
PS1="${PS1-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${SCRIPT_DIR}"
while [[ "${REPO_ROOT}" != "/" && ! ( -f "${REPO_ROOT}/pyproject.toml" && -d "${REPO_ROOT}/pipeline" ) ]]; do
  REPO_ROOT="$(cd "${REPO_ROOT}/.." && pwd)"
done
if [[ ! -f "${REPO_ROOT}/pyproject.toml" || ! -d "${REPO_ROOT}/pipeline" ]]; then
  echo "Error: could not find FedGWAS repository root from ${SCRIPT_DIR}" >&2
  exit 1
fi

# Prefer fedgwas conda env (Matpool/Jupyter often leave user in base/myconda)
SUPERNODE_BIN=""
for _bin in \
  "${FEDGWAS_BIN:-}" \
  "${CONDA_PREFIX:-}/bin" \
  "/root/miniconda3/envs/fedgwas/bin" \
  "${HOME}/miniconda3/envs/fedgwas/bin"; do
  if [[ -n "${_bin}" && -x "${_bin}/flower-supernode" ]]; then
    SUPERNODE_BIN="${_bin}/flower-supernode"
    export PATH="${_bin}:${PATH}"
    break
  fi
done
unset _bin

SERVER_IP=""
CLIENT_ID=""
CONFIG=""
PORT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --server-ip) SERVER_IP="$2"; shift 2 ;;
    --client-id) CLIENT_ID="$2"; shift 2 ;;
    --config) CONFIG="$2"; shift 2 ;;
    --port) PORT="$2"; shift 2 ;;
    *)
      echo "Unknown arg: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "${SERVER_IP}" || -z "${CLIENT_ID}" || -z "${CONFIG}" ]]; then
  echo "Missing required args. Usage: --server-ip <ip> --client-id <1|2> --config <path> [--port <port>]" >&2
  exit 1
fi

if [[ "${CLIENT_ID}" != "1" && "${CLIENT_ID}" != "2" ]]; then
  echo "client-id must be 1 or 2" >&2
  exit 1
fi

PARTITION_ID=$((CLIENT_ID - 1))
if [[ -z "${PORT}" ]]; then
  if [[ "${CLIENT_ID}" == "1" ]]; then PORT=9094; else PORT=9095; fi
fi

# Resolve config relative to repo root if not absolute
if [[ "${CONFIG}" != /* ]]; then
  if [[ -f "${REPO_ROOT}/${CONFIG}" ]]; then
    CONFIG="${REPO_ROOT}/${CONFIG}"
  fi
fi

cd "${REPO_ROOT}"
export FEDGWAS_REPO_ROOT="${REPO_ROOT}"

_on_exit() {
  local code=$?
  if [[ "${code}" -ne 0 ]]; then
    echo "[client${CLIENT_ID}] ERROR: flower-supernode exited with code ${code}" >&2
    echo "[client${CLIENT_ID}] Check: SuperLink on ${SERVER_IP}:9092, clock skew, port ${PORT} free, flwr 1.19.x" >&2
  fi
  echo "[client${CLIENT_ID}] finished at $(date -u '+%Y-%m-%d %H:%M:%S UTC') exit=${code}"
}
trap _on_exit EXIT

echo "[client${CLIENT_ID}] repo=${REPO_ROOT} cwd=$(pwd) pid=$$"
echo "[client${CLIENT_ID}] superlink ${SERVER_IP}:9092 clientappio 0.0.0.0:${PORT} partition-id=${PARTITION_ID}"
echo "[client${CLIENT_ID}] config=${CONFIG}"

if [[ ! -f "${CONFIG}" ]]; then
  echo "[client${CLIENT_ID}] ERROR: config not found: ${CONFIG}" >&2
  echo "[client${CLIENT_ID}] Run from repo root (cd /Fed-GWAS) or pass an absolute config path." >&2
  exit 1
fi

if [[ -z "${SUPERNODE_BIN}" ]]; then
  echo "[client${CLIENT_ID}] ERROR: flower-supernode not found." >&2
  echo "  export PATH=/root/miniconda3/envs/fedgwas/bin:\$PATH" >&2
  echo "  pip install 'flwr[simulation]==1.19.0'" >&2
  exit 127
fi
echo "[client${CLIENT_ID}] flower-supernode=${SUPERNODE_BIN} ($(flower-supernode --version 2>/dev/null || echo unknown))"

if command -v ss >/dev/null 2>&1; then
  if ss -tln | grep -q ":${PORT} "; then
    echo "[client${CLIENT_ID}] WARNING: port ${PORT} already in use — stop old SuperNode first:" >&2
    echo "  cluster_deployment/scripts/cluster-stop-all.sh" >&2
  fi
fi

if ! (echo >/dev/tcp/"${SERVER_IP}"/9092) 2>/dev/null; then
  if ! python3 -c "import socket; socket.create_connection(('${SERVER_IP}', 9092), timeout=3).close()" 2>/dev/null; then
    echo "[client${CLIENT_ID}] WARNING: cannot reach ${SERVER_IP}:9092 — start SuperLink on server first." >&2
  fi
fi

echo "[client${CLIENT_ID}] Tip: Invalid timestamp → bash cluster_deployment/scripts/check-cluster-clock-skew.sh ${SERVER_IP}"
export FEDGWAS_FLR_REPLY_SKEW_SEC="${FEDGWAS_FLR_REPLY_SKEW_SEC:-2}"

SUPERNODE_LAUNCH=("${SUPERNODE_BIN}")
if [[ -n "${FEDGWAS_FAKETIME:-}" ]]; then
  if command -v faketime >/dev/null 2>&1; then
    echo "[client${CLIENT_ID}] using faketime offset ${FEDGWAS_FAKETIME} (Matpool clock workaround)"
    # Propagate to ClientApp child processes spawned by SuperNode
    export FAKETIME="${FEDGWAS_FAKETIME}"
    export FAKETIME_DONT_RESET=1
    for _lib in \
      /usr/lib/x86_64-linux-gnu/faketime/libfaketime.so.1 \
      /usr/lib/faketime/libfaketime.so.1 \
      /lib/x86_64-linux-gnu/faketime/libfaketime.so.1; do
      if [[ -f "${_lib}" ]]; then
        export LD_PRELOAD="${_lib}${LD_PRELOAD:+:${LD_PRELOAD}}"
        break
      fi
    done
    unset _lib
    export FEDGWAS_FLR_REPLY_SKEW_SEC="${FEDGWAS_FLR_REPLY_SKEW_SEC:-3}"
    SUPERNODE_LAUNCH=(faketime -f "${FEDGWAS_FAKETIME}" "${SUPERNODE_BIN}")
  else
    echo "[client${CLIENT_ID}] ERROR: FEDGWAS_FAKETIME=${FEDGWAS_FAKETIME} but faketime not installed" >&2
    echo "  apt-get install -y faketime" >&2
    exit 1
  fi
fi

"${SUPERNODE_LAUNCH[@]}" \
  --insecure \
  --superlink "${SERVER_IP}:9092" \
  --clientappio-api-address "0.0.0.0:${PORT}" \
  --node-config "partition-id=${PARTITION_ID} num-partitions=2 config-file=\"${CONFIG}\""
