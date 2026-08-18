#!/usr/bin/env bash
set -eo pipefail
PS1="${PS1-}"

for _bin in \
  "${FEDGWAS_BIN:-}" \
  "${CONDA_PREFIX:-}/bin" \
  "/root/miniconda3/envs/fedgwas/bin" \
  "${HOME}/miniconda3/envs/fedgwas/bin"; do
  if [[ -n "${_bin}" && -x "${_bin}/flower-superlink" ]]; then
    export PATH="${_bin}:${PATH}"
    break
  fi
done
unset _bin

# Start Flower SuperLink on the server node.
# Usage: ./scripts/cluster-start-server.sh [--server-ip <ip_or_hostname>]
# Default binds to 0.0.0.0 (all interfaces).

SERVER_IP="0.0.0.0"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --server-ip)
      SERVER_IP="$2"
      shift 2
      ;;
    *)
      echo "Unknown arg: $1" >&2
      exit 1
      ;;
  esac
done

echo "[server] Starting flower-superlink on ${SERVER_IP}:9092 (exec API on 9093)"
if ! command -v flower-superlink &>/dev/null; then
  echo "Error: flower-superlink not found. export PATH=/root/miniconda3/envs/fedgwas/bin:\$PATH" >&2
  exit 127
fi
flower-superlink \
  --insecure \
  --fleet-api-address "${SERVER_IP}:9092" \
  --exec-api-address "${SERVER_IP}:9093"
