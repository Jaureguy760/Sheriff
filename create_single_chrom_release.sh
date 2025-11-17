#!/bin/bash
# Create single-chromosome test package for GitHub release
# Complete test dataset: BAM + reference + GTF + barcodes + edit sites (~23MB)

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "Creating Single-Chromosome Test Package"
echo "========================================"
echo ""

# Files to include
INCLUDE_FILES=(
    "test_data/test_chr19.bam"
    "test_data/test_chr19.bam.bai"
    "test_data/chr19.fa.gz"
    "test_data/chr19.gtf.gz"
    "test_data/barcodes_chr19.txt"
    "test_data/edit_sites_chr19.txt"
    "test_data/blacklist.bed"
    "test_data/blacklist_seqs.txt"
    "test_data/validate_chr19_test.py"
    "test_data/ci_validation.py"
    "test_data/expected_checksums.json"
    "test_data/generate_checksums.py"
    "test_data/README.md"
)

echo "Files to package:"
for f in "${INCLUDE_FILES[@]}"; do
    if [ -f "$f" ]; then
        SIZE=$(du -h "$f" | cut -f1)
        echo "  $f ($SIZE)"
    else
        echo "  WARNING: $f not found!"
    fi
done

# Create tarball
echo ""
echo "Creating test_data_chr19.tar.gz..."
tar -czf test_data_chr19.tar.gz "${INCLUDE_FILES[@]}"

# Show result
COMPRESSED_SIZE=$(du -h test_data_chr19.tar.gz | cut -f1)
echo ""
echo "========================================"
echo "Package created successfully!"
echo "========================================"
echo "  File: test_data_chr19.tar.gz"
echo "  Size: $COMPRESSED_SIZE"
echo ""
echo "Package contents:"
echo "  - test_chr19.bam: 42,194 reads on chromosome 19"
echo "  - chr19.fa.gz: GRCh38 chromosome 19 (58MB uncompressed)"
echo "  - chr19.gtf.gz: Ensembl 110 annotations for chr19"
echo "  - barcodes_chr19.txt: 4,680 cell barcodes"
echo "  - edit_sites_chr19.txt: 4 known edit sites"
echo "  - Validation scripts and checksums"
echo ""
echo "This package is COMPLETE - no additional downloads needed!"
echo ""
echo "To upload to GitHub release:"
echo "  1. Go to https://github.com/BradBalderson/Sheriff/releases"
echo "  2. Create new release or edit existing"
echo "  3. Attach test_data_chr19.tar.gz as release asset"
echo ""
echo "Users can download and test with:"
echo "  wget -O- -q https://github.com/BradBalderson/Sheriff/releases/latest/download/test_data_chr19.tar.gz | tar xvzf -"
echo "  gunzip test_data/chr19.fa.gz test_data/chr19.gtf.gz"
echo "  python test_data/validate_chr19_test.py"
