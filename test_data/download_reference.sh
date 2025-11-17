#!/bin/bash
# Download reference genome for Sheriff test data
# This downloads the GRCh38 reference from Ensembl

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "Sheriff Test Data - Reference Genome Download"
echo "=============================================="
echo ""

# Check for required tools
command -v wget >/dev/null 2>&1 || { echo "Error: wget is required but not installed."; exit 1; }

# Reference genome options
ENSEMBL_VERSION="110"
GENOME_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_VERSION}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_VERSION}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${ENSEMBL_VERSION}.gtf.gz"

# Check if already downloaded
if [ -f "reference.fa" ] && [ -f "reference.gtf" ]; then
    echo "Reference files already exist. To re-download, remove them first:"
    echo "  rm reference.fa reference.gtf"
    exit 0
fi

echo "Downloading GRCh38 reference genome (Ensembl release ${ENSEMBL_VERSION})..."
echo "This may take 30-60 minutes depending on your connection."
echo ""

# Download and decompress reference FASTA
if [ ! -f "reference.fa" ]; then
    echo "[1/4] Downloading reference FASTA (~900MB compressed, ~3GB uncompressed)..."
    wget -q --show-progress -O reference.fa.gz "$GENOME_URL"

    echo "[2/4] Decompressing reference FASTA..."
    gunzip reference.fa.gz

    echo "[3/4] Indexing reference FASTA..."
    if command -v samtools >/dev/null 2>&1; then
        samtools faidx reference.fa
    else
        echo "Warning: samtools not found. Skipping FASTA indexing."
        echo "Install samtools and run: samtools faidx reference.fa"
    fi
fi

# Download and decompress GTF
if [ ! -f "reference.gtf" ]; then
    echo "[4/4] Downloading GTF annotation (~50MB compressed)..."
    wget -q --show-progress -O reference.gtf.gz "$GTF_URL"
    gunzip reference.gtf.gz
fi

echo ""
echo "Download complete!"
echo "Files created:"
ls -lh reference.fa reference.gtf 2>/dev/null || true
echo ""
echo "You can now run the Sheriff test pipeline."
