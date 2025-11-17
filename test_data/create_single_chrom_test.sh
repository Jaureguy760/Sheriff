#!/bin/bash
# Create single-chromosome test dataset for Sheriff
# Uses chr19 (58MB reference, 42k reads) - smallest reference with meaningful data

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Full reference (you need access to this)
FULL_REF="${1:-/iblm/netapp/data4/bbalderson/edit_capture/data/t7_indel_calling/Homo_sapiens.GRCh38.dna.primary_assembly.fa}"
FULL_GTF="${2:-/iblm/netapp/data4/bbalderson/edit_capture/data/t7_indel_calling/Homo_sapiens.GRCh38.110.gtf}"

if [ ! -f "$FULL_REF" ]; then
    echo "ERROR: Full reference not found: $FULL_REF"
    echo "Usage: $0 [full_reference.fa] [full_gtf]"
    exit 1
fi

TARGET_CHROM="hg38_19"  # 58MB reference, 42k reads

echo "Creating single-chromosome test dataset"
echo "========================================"
echo "Target: $TARGET_CHROM"
echo ""

# Step 1: Extract BAM for single chromosome
echo "[1/5] Extracting BAM for $TARGET_CHROM..."
samtools view -b test_200kb.bam "$TARGET_CHROM" > test_single_chrom.bam
samtools index test_single_chrom.bam

NUM_READS=$(samtools view -c test_single_chrom.bam)
BAM_SIZE=$(du -h test_single_chrom.bam | cut -f1)
echo "  Reads: $NUM_READS"
echo "  BAM size: $BAM_SIZE"

# Step 2: Extract chromosome reference
echo "[2/5] Extracting $TARGET_CHROM from reference..."
samtools faidx "$FULL_REF" "$TARGET_CHROM" > reference_single_chrom.fa
samtools faidx reference_single_chrom.fa

REF_SIZE=$(du -h reference_single_chrom.fa | cut -f1)
echo "  Reference size: $REF_SIZE"

# Step 3: Extract GTF for chromosome
echo "[3/5] Extracting GTF annotations..."
awk -v chr="$TARGET_CHROM" '
    /^#/ {print; next}
    $1 == chr {print}
' "$FULL_GTF" > reference_single_chrom.gtf

GTF_LINES=$(grep -v "^#" reference_single_chrom.gtf | wc -l)
GTF_SIZE=$(du -h reference_single_chrom.gtf | cut -f1)
echo "  GTF entries: $GTF_LINES"
echo "  GTF size: $GTF_SIZE"

# Step 4: Filter barcodes to only those with reads on this chromosome
echo "[4/5] Filtering barcodes..."
samtools view test_single_chrom.bam | \
    awk -F'\t' '{
        for(i=12; i<=NF; i++) {
            if($i ~ /^CB:Z:/) {
                split($i, a, ":")
                print a[3]
            }
        }
    }' | sort | uniq > barcodes_single_chrom.txt

NUM_BARCODES=$(wc -l < barcodes_single_chrom.txt)
echo "  Barcodes with reads: $NUM_BARCODES"

# Step 5: Summary
echo ""
echo "[5/5] Compressing reference (for smaller package)..."
gzip -k reference_single_chrom.fa
COMPRESSED_SIZE=$(du -h reference_single_chrom.fa.gz | cut -f1)
echo "  Compressed reference: $COMPRESSED_SIZE"

echo ""
echo "========================================"
echo "Single-chromosome test dataset created!"
echo "========================================"
echo "Files:"
echo "  test_single_chrom.bam:          $BAM_SIZE ($NUM_READS reads)"
echo "  reference_single_chrom.fa:      $REF_SIZE"
echo "  reference_single_chrom.fa.gz:   $COMPRESSED_SIZE"
echo "  reference_single_chrom.gtf:     $GTF_SIZE"
echo "  barcodes_single_chrom.txt:      $NUM_BARCODES barcodes"
echo ""
echo "Total package size (compressed): ~$(( $(stat --printf="%s" test_single_chrom.bam) / 1000000 + $(stat --printf="%s" reference_single_chrom.fa.gz) / 1000000 ))MB"
echo ""
echo "This is a fully self-contained test that can be:"
echo "- Uploaded to GitHub releases"
echo "- Run without external dependencies"
echo "- Used for CI full-pipeline testing"
