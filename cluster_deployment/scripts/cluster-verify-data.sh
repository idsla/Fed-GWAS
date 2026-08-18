#!/usr/bin/env bash
# Verify PLINK triplet exists and is readable on this node.
# On client nodes pass --client-id 1 or 2 (only that center's data is required).
# Omit --client-id (or use --client-id all) to check both centers (dev machine / server).
set -eo pipefail
PS1="${PS1-}"

SCALE="tiny"
CLIENT_ID="all"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${SCRIPT_DIR}"
while [[ "${REPO_ROOT}" != "/" && ! ( -f "${REPO_ROOT}/pyproject.toml" && -d "${REPO_ROOT}/pipeline" ) ]]; do
  REPO_ROOT="$(cd "${REPO_ROOT}/.." && pwd)"
done
if [[ ! -f "${REPO_ROOT}/pyproject.toml" || ! -d "${REPO_ROOT}/pipeline" ]]; then
  echo "Error: could not find FedGWAS repository root from ${SCRIPT_DIR}" >&2
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --scale) SCALE="$2"; shift 2 ;;
    --client-id) CLIENT_ID="$2"; shift 2 ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

case "${SCALE}" in
  tiny)
    BASE="experiments/correctness/tiny_even/data/tiny"
    PFX="tiny_center"
    ;;
  small)
    BASE="experiments/performance/small_even/data/small"
    PFX="small_center"
    ;;
  medium)
    BASE="experiments/performance/medium_even/data/medium"
    PFX="medium_center"
    ;;
  *)
    echo "Unknown scale: ${SCALE} (tiny|small|medium)" >&2
    exit 1
    ;;
esac

case "${CLIENT_ID}" in
  all|both)
    CENTERS=(1 2)
    ;;
  1|2)
    CENTERS=("${CLIENT_ID}")
    ;;
  *)
    echo "Unknown --client-id: ${CLIENT_ID} (use 1, 2, or all)" >&2
    exit 1
    ;;
esac

cd "${REPO_ROOT}"
echo "=== PLINK data check (scale=${SCALE}, client-id=${CLIENT_ID}, cwd=${REPO_ROOT}) ==="

fail=0
for c in "${CENTERS[@]}"; do
  prefix="${BASE}/center_${c}/${PFX}_${c}"
  echo ""
  echo "--- center_${c}: ${prefix} ---"
  ok=true
  for ext in bed bim fam; do
    f="${prefix}.${ext}"
    if [[ ! -f "${f}" ]]; then
      echo "  MISSING ${f}"
      ok=false
      fail=1
    else
      echo "  OK ${f} ($(wc -c <"${f}" | tr -d ' ') bytes)"
    fi
  done
  if [[ "${ok}" == true ]]; then
    if command -v plink &>/dev/null; then
      if plink --bfile "${prefix}" --freq --out "/tmp/fedgwas_data_check_${SCALE}_${c}" >/dev/null 2>&1; then
        echo "  OK plink --bfile reads triplet"
      else
        echo "  FAIL plink --bfile (corrupt or wrong format)"
        fail=1
      fi
    else
      echo "  WARN plink not in PATH — skip read test"
    fi
  fi
done

echo ""
if [[ "${fail}" -eq 0 ]]; then
  echo "All checks passed."
else
  echo "Data missing or unreadable. Copy center_* PLINK files to this node before starting SuperNode." >&2
  exit 1
fi
