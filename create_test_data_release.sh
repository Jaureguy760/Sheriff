#!/bin/bash
# Create test_data.tar.gz for GitHub release
# This packages all test data files (excluding large reference genome)

set -e

echo "Creating Sheriff Test Data Release Package"
echo "==========================================="
echo ""

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Verify test_data directory exists
if [ ! -d "test_data" ]; then
    echo "ERROR: test_data directory not found"
    exit 1
fi

# Files to include (everything except reference genome)
INCLUDE_FILES=(
    "test_data/test_200kb.bam"
    "test_data/test_200kb.bam.bai"
    "test_data/barcodes.txt"
    "test_data/edit_sites.txt"
    "test_data/blacklist.bed"
    "test_data/blacklist_seqs.txt"
    "test_data/download_reference.sh"
    "test_data/run_validation.py"
    "test_data/README.md"
    "test_data/expected_checksums.json"
    "test_data/generate_checksums.py"
    "test_data/ci_validation.py"
)

echo "Files to package:"
TOTAL_SIZE=0
for f in "${INCLUDE_FILES[@]}"; do
    if [ -f "$f" ]; then
        SIZE=$(du -h "$f" | cut -f1)
        echo "  $f ($SIZE)"
        TOTAL_SIZE=$((TOTAL_SIZE + $(stat --printf="%s" "$f")))
    else
        echo "  WARNING: $f not found!"
    fi
done

echo ""
echo "Total uncompressed size: $(numfmt --to=iec $TOTAL_SIZE)"

# Create tarball
echo ""
echo "Creating test_data.tar.gz..."
tar -czf test_data.tar.gz "${INCLUDE_FILES[@]}"

# Show result
COMPRESSED_SIZE=$(du -h test_data.tar.gz | cut -f1)
echo ""
echo "Package created successfully!"
echo "  File: test_data.tar.gz"
echo "  Size: $COMPRESSED_SIZE"
echo ""
echo "To upload to GitHub release:"
echo "  1. Go to https://github.com/BradBalderson/Sheriff/releases"
echo "  2. Create new release or edit existing"
echo "  3. Attach test_data.tar.gz as release asset"
echo ""
echo "Users can download with:"
echo "  wget -O- -q https://github.com/BradBalderson/Sheriff/releases/latest/download/test_data.tar.gz | tar xvzf -"
