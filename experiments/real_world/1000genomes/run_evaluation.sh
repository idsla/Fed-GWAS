#!/bin/bash
# Run evaluation and regenerate manuscript materials after experiment completes
# Usage: ./run_evaluation.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/../../.."

# Use the Python from the conda environment if available
# Check for conda environment Python first
if [ -n "$CONDA_PREFIX" ] && [ -f "$CONDA_PREFIX/bin/python" ]; then
    PYTHON_CMD="$CONDA_PREFIX/bin/python"
elif command -v python &> /dev/null; then
    PYTHON_CMD="python"
elif command -v python3 &> /dev/null; then
    PYTHON_CMD="python3"
else
    echo "Error: Python not found. Please activate your conda environment first." >&2
    exit 1
fi

RESULTS_DIR="experiments/real_world/1000genomes/results"
BASELINE_DIR="experiments/real_world/1000genomes/data/centralized_baseline"
MANUSCRIPT_DIR="experiments/real_world/1000genomes/manuscript"

echo "=========================================="
echo "Running Evaluation and Generating Manuscript Materials"
echo "=========================================="
echo ""
echo "Using Python: $($PYTHON_CMD --version)"
echo "Python path: $(which $PYTHON_CMD)"
echo ""

# Step 1: Run evaluation
echo "Step 1: Running QC, LR, and KING evaluation..."
$PYTHON_CMD experiments/tools/evaluation/evaluate_all.py \
    "$RESULTS_DIR" \
    --baseline "$BASELINE_DIR" \
    --king \
    --king-center-id 2

echo ""
echo "Step 2: Generating manuscript materials..."
$PYTHON_CMD experiments/tools/generate_manuscript_materials.py \
    "$RESULTS_DIR" \
    --baseline "$BASELINE_DIR" \
    --output "$MANUSCRIPT_DIR"

echo ""
echo "=========================================="
echo "Evaluation and manuscript generation completed!"
echo "=========================================="
echo ""
echo "Results are available in:"
echo "  - Evaluation reports: $RESULTS_DIR/*_report.md"
echo "  - Manuscript materials: $MANUSCRIPT_DIR/"
echo "    - Figures: $MANUSCRIPT_DIR/figures/"
echo "    - Tables: $MANUSCRIPT_DIR/tables/"
echo "    - Text: $MANUSCRIPT_DIR/manuscript_sections.md"
