#!/usr/bin/env bash
# Download flwr 1.19.x wheels on a dev machine (with PyPI access), then scp to Matpool nodes.
set -eo pipefail

OUT="${1:-./cluster_deployment/wheels}"
mkdir -p "${OUT}"

PYTHON="${PYTHON:-python3.11}"
if ! command -v "${PYTHON}" &>/dev/null; then
  PYTHON=python3
fi

echo "Downloading flwr[simulation]==1.19.0 wheels into ${OUT} ..."
"${PYTHON}" -m pip download -d "${OUT}" \
  "flwr[simulation]==1.19.0" tomli tomli-w psutil

echo ""
echo "Done. Upload to each node, then on the node:"
echo "  export PATH=/root/miniconda3/envs/fedgwas/bin:\$PATH"
echo "  pip install --no-index --find-links=/path/to/wheels 'flwr[simulation]==1.19.0' tomli tomli-w psutil"
echo ""
echo "  scp -r ${OUT} user@192.168.1.88:/tmp/flwr-wheels"
