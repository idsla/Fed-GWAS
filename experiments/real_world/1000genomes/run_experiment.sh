#!/bin/bash
# Run federated GWAS experiment for 1000 Genomes Project
# This script runs the experiment with the optimized threshold (0.1)

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/../../.."

echo "=========================================="
echo "Running 1000 Genomes Federated GWAS Experiment"
echo "Threshold: p = 0.1 (optimized)"
echo "=========================================="
echo ""

# Check if conda environment is activated
if [[ -z "$CONDA_DEFAULT_ENV" ]]; then
    echo "Warning: Conda environment not detected. Make sure to activate your environment first."
    echo "Example: conda activate fedgwas"
    echo ""
fi

# Run the experiment
echo "Starting federated experiment..."
flwr run . local-simulation --stream \
    --run-config "config_path=experiments/real_world/1000genomes/configs simulation=true num-server-rounds=350"

echo ""
echo "=========================================="
echo "Experiment completed!"
echo "=========================================="
