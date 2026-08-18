#!/usr/bin/env bash
# Quick cluster health check (run on server and each client).
set -eo pipefail
PS1="${PS1-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SCALE="tiny"
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --scale) SCALE="$2"; shift 2 ;;
    *) POSITIONAL+=("$1"); shift ;;
  esac
done

SERVER_IP="${POSITIONAL[0]:-192.168.1.88}"
CLIENT1_IP="${POSITIONAL[1]:-}"
CLIENT2_IP="${POSITIONAL[2]:-}"

echo "=== Fed-GWAS cluster diagnose ==="
echo "Server SuperLink: ${SERVER_IP}:9092 / Exec 9093"
echo ""

echo "--- Local processes ---"
cluster_deployment/scripts/cluster-status.sh 2>/dev/null || true
echo ""

case "${SCALE}" in
  tiny)
    DATA_BASE="experiments/correctness/tiny_even/data/tiny"
    RESULTS_BASE="experiments/correctness/tiny_even/results_2"
    ;;
  small)
    DATA_BASE="experiments/performance/small_even/data/small"
    RESULTS_BASE="experiments/performance/small_even/results"
    ;;
  medium)
    DATA_BASE="experiments/performance/medium_even/data/medium"
    RESULTS_BASE="experiments/performance/medium_even/results"
    ;;
  *)
    DATA_BASE="experiments/correctness/tiny_even/data/tiny"
    RESULTS_BASE="experiments/correctness/tiny_even/results_2"
    ;;
esac

echo "--- Local PLINK data (${SCALE}) ---"
for c in center_1 center_2; do
  p="${DATA_BASE}/${c}"
  if ls "${p}"/*.bed 1>/dev/null 2>&1; then
    echo "  OK ${p}"
  else
    echo "  MISSING ${p}/*.bed — on this client run: bash cluster_deployment/scripts/cluster-verify-data.sh --scale ${SCALE} --client-id ${c#center_}"
  fi
done
echo ""

echo "--- Local client logs ---"
for c in center_1 center_2; do
  d="${RESULTS_BASE}/${c}/logs"
  if [[ -d "${d}" ]]; then
    n="$(find "${d}" -type f 2>/dev/null | wc -l | tr -d ' ')"
    latest="$(ls -t "${d}"/* 2>/dev/null | head -1 || true)"
    echo "  ${c}: ${n} files, latest: ${latest:-none}"
  else
    echo "  ${c}: no ${d} yet"
  fi
done
echo ""

if [[ -n "${CLIENT1_IP}" ]]; then
  echo "--- Reachability from this host ---"
  for spec in "client1:${CLIENT1_IP}:9094" "client2:${CLIENT2_IP}:9095"; do
    label="${spec%%:*}"
    rest="${spec#*:}"
    ip="${rest%%:*}"
    port="${rest##*:}"
    [[ -z "${ip}" ]] && continue
    if (echo >/dev/tcp/"${ip}"/"${port}") 2>/dev/null || \
       python3 -c "import socket; socket.create_connection(('${ip}', ${port}), timeout=2).close()" 2>/dev/null; then
      echo "  OK ${label} ${ip}:${port}"
    else
      echo "  FAIL ${label} ${ip}:${port} (SuperLink on server must reach ClientAppIo)"
    fi
  done
  echo ""
fi

if command -v ssh >/dev/null 2>&1; then
  echo "--- Clock skew vs server ---"
  bash cluster_deployment/scripts/check-cluster-clock-skew.sh "${SERVER_IP}" || true
fi

echo ""
echo "Notes:"
echo "  • SuperNode stuck after 'ClientAppIo gRPC server' is NORMAL idle until a round starts."
echo "  • Round progress appears in results_2/center_*/logs/, not always in SuperNode stdout."
echo "  • If only one client has logs while round1 runs, the other SuperNode is not in the federation."
echo "  • pyproject expects 2 SuperNodes; both must register with SuperLink before flwr run proceeds fully."
