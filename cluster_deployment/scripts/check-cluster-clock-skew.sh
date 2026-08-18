#!/usr/bin/env bash
# Compare local UTC clock to SuperLink server.
# Flower 1.19 rejects gRPC when skew is too large (either direction).
set -euo pipefail

SERVER_IP="${1:-}"
MAX_SKEW_MS="${MAX_CLOCK_SKEW_MS:-500}"

if [[ -z "${SERVER_IP}" ]]; then
  echo "Usage: $0 <server-ip>" >&2
  echo "Example: $0 192.168.1.88" >&2
  exit 1
fi

LOCAL_TS="$(date -u +%s.%N)"
echo "Local UTC:  $(date -u '+%Y-%m-%d %H:%M:%S') (${LOCAL_TS})"

report_skew() {
  python3 - << PY
local = float("${LOCAL_TS}")
remote = float("${1}")
max_ms = float("${MAX_SKEW_MS}")
delta_ms = (local - remote) * 1000.0
abs_ms = abs(delta_ms)
print(f"Skew (local - server): {delta_ms:+.3f} ms  (|skew| must be < {max_ms:.0f} ms)")
ok = True
if abs_ms > max_ms:
    ok = False
    if delta_ms < 0:
        print("FAIL: local clock is BEHIND server.")
    else:
        print("FAIL: local clock is AHEAD of server.")
    print("  Flower 1.19 may reject SuperNode with: Invalid timestamp / UNAUTHENTICATED")
    print("  1) Sync NTP on ALL nodes (restart VM on Matpool if ntpdate is forbidden)")
    offset_sec = -(delta_ms / 1000.0)
    if offset_sec >= 0:
        ft = f"+{offset_sec:.3f}s"
    else:
        ft = f"{offset_sec:.3f}s"
    print("  2) Matpool workaround — compensate only this SuperNode process:")
    print(f'     apt-get install -y faketime')
    print(f'     export FEDGWAS_FAKETIME="{ft}"')
    print(f'     cluster_deployment/scripts/cluster-start-client.sh ...')
elif delta_ms < -200:
    print("WARNING: local slightly behind server (<500ms but risky).")
elif delta_ms > 200:
    print("WARNING: local slightly ahead of server (<500ms but risky).")
else:
    print("OK: clocks aligned within tolerance.")
import sys
sys.exit(0 if ok else 1)
PY
}

if command -v ssh >/dev/null 2>&1; then
  if REMOTE_TS="$(ssh -o ConnectTimeout=3 -o BatchMode=yes "root@${SERVER_IP}" 'date -u +%s.%N' 2>/dev/null)"; then
    echo "Server UTC: $(ssh -o ConnectTimeout=3 "root@${SERVER_IP}" 'date -u "+%Y-%m-%d %H:%M:%S"') (${REMOTE_TS})"
    report_skew "${REMOTE_TS}"
    exit $?
  fi
fi

echo "Could not SSH to root@${SERVER_IP}. On Server and this node, run simultaneously:"
echo "  date -u +%s.%N"
echo "|local - server| must be < ${MAX_SKEW_MS} ms (both behind and ahead break Flower 1.19)."
exit 1
