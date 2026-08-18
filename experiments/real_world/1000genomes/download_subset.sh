#!/bin/bash
# Download a small subset of 1000 Genomes Project data
# This script downloads chromosome 22 (smallest chromosome) for testing

set -e

# Get script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${PROJECT_ROOT}"

DATA_DIR="experiments/real_world/1000genomes/data"
# Use HTTPS instead of FTP for better compatibility
BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
SAMPLE_BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

# Create data directory
mkdir -p "${DATA_DIR}"

echo "Downloading 1000 Genomes Project subset (chromosome 22)..."
echo "This is a small subset for testing (~500MB compressed)"
echo ""

# Download chromosome 22 VCF file
CHR22_VCF="${DATA_DIR}/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
CHR22_VCF_IDX="${DATA_DIR}/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi"

# Use curl if wget is not available
if command -v wget &> /dev/null; then
    DOWNLOAD_CMD="wget -c"
    PROGRESS_FLAG="--progress=bar:force"
elif command -v curl &> /dev/null; then
    DOWNLOAD_CMD="curl -L"
    PROGRESS_FLAG="--progress-bar"
else
    echo "Error: Neither wget nor curl found. Please install one of them."
    exit 1
fi

if [ ! -f "${CHR22_VCF}" ]; then
    echo "Downloading chromosome 22 VCF (~500MB, this may take a few minutes)..."
    ${DOWNLOAD_CMD} "${BASE_URL}/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" \
        ${PROGRESS_FLAG} -o "${CHR22_VCF}"
    echo "✓ Downloaded: ${CHR22_VCF}"
else
    echo "✓ Already exists: ${CHR22_VCF}"
fi

# Download index file
if [ ! -f "${CHR22_VCF_IDX}" ]; then
    echo "Downloading index file..."
    ${DOWNLOAD_CMD} "${BASE_URL}/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi" \
        ${PROGRESS_FLAG} -o "${CHR22_VCF_IDX}"
    echo "✓ Downloaded: ${CHR22_VCF_IDX}"
else
    echo "✓ Already exists: ${CHR22_VCF_IDX}"
fi

# Download sample information
SAMPLE_INFO="${DATA_DIR}/integrated_call_samples_v3.20130502.ALL.panel"
if [ ! -f "${SAMPLE_INFO}" ]; then
    echo "Downloading sample information (population labels)..."
    ${DOWNLOAD_CMD} "${SAMPLE_BASE_URL}/integrated_call_samples_v3.20130502.ALL.panel" \
        ${PROGRESS_FLAG} -o "${SAMPLE_INFO}"
    echo "✓ Downloaded: ${SAMPLE_INFO}"
else
    echo "✓ Already exists: ${SAMPLE_INFO}"
fi

echo ""
echo "Download complete!"
echo ""
echo "Files downloaded:"
echo "  - VCF: ${CHR22_VCF}"
echo "  - Index: ${CHR22_VCF_IDX}"
echo "  - Sample info: ${SAMPLE_INFO}"
echo ""
echo "Next steps:"
echo "  1. Convert VCF to PLINK format:"
echo "     plink --vcf ${CHR22_VCF} --make-bed --out ${DATA_DIR}/genotypes_raw"
echo ""
echo "  2. Filter multiallelic variants (standard GWAS preprocessing):"
echo "     python experiments/tools/filter_multiallelic.py \\"
echo "         --plink-prefix ${DATA_DIR}/genotypes_raw \\"
echo "         --output-prefix ${DATA_DIR}/genotypes"
echo ""
echo "  3. Generate phenotypes:"
echo "     python -m experiments.real_world.phenotype_generation.generate_phenotypes \\"
echo "         --plink-prefix ${DATA_DIR}/genotypes \\"
echo "         --output-fam ${DATA_DIR}/genotypes_pheno.fam \\"
echo "         --model additive --causal-fraction 0.01 --seed 42"
echo ""
echo "  4. Partition data by population:"
echo "     (See README.md for detailed partitioning instructions)"
echo ""