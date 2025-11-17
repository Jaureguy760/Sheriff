#!/bin/bash
# Create minimal reference genome for test_200kb.bam
# Extracts only the chromosomes and regions needed (with padding)

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Full reference (you need access to this to run this script)
FULL_REF="${1:-/iblm/netapp/data4/bbalderson/edit_capture/data/t7_indel_calling/Homo_sapiens.GRCh38.dna.primary_assembly.fa}"
FULL_GTF="${2:-/iblm/netapp/data4/bbalderson/edit_capture/data/t7_indel_calling/Homo_sapiens.GRCh38.110.gtf}"

if [ ! -f "$FULL_REF" ]; then
    echo "ERROR: Full reference not found: $FULL_REF"
    echo "Usage: $0 [full_reference.fa] [full_gtf]"
    exit 1
fi

echo "Creating minimal reference for test_200kb.bam"
echo "=============================================="
echo ""

# Step 1: Get chromosomes with reads
echo "[1/5] Identifying chromosomes with reads..."
CHROMOSOMES=$(samtools idxstats test_200kb.bam | awk '$3 > 0 {print $1}' | sort -V | tr '\n' ' ')
echo "  Chromosomes: $CHROMOSOMES"

# Step 2: Create BED file of regions with padding
echo "[2/5] Creating region BED file (10kb padding around reads)..."
samtools view test_200kb.bam | \
    awk '{print $3"\t"$4-10000"\t"$4+length($10)+10000}' | \
    awk '$2 < 0 {$2 = 0} {print}' OFS="\t" | \
    sort -k1,1 -k2,2n | \
    awk -v OFS="\t" '
        NR==1 {chr=$1; start=$2; end=$3; next}
        $1==chr && $2<=end {if($3>end) end=$3; next}
        {print chr, start, end; chr=$1; start=$2; end=$3}
        END {print chr, start, end}
    ' > regions.bed

NUM_REGIONS=$(wc -l < regions.bed)
TOTAL_BASES=$(awk '{sum += $3-$2} END {print sum}' regions.bed)
echo "  Merged regions: $NUM_REGIONS"
echo "  Total bases (with padding): $TOTAL_BASES (~$((TOTAL_BASES/1000000))MB)"

# Step 3: Extract chromosome sequences
# NOTE: This approach extracts full chromosomes then subsets
# For a truly minimal reference, we'd need to extract just the regions
# But Sheriff needs the full chromosome for coordinate mapping

echo "[3/5] Extracting chromosome FASTAs..."
# Option A: Full chromosomes (larger but simpler)
# for chr in $CHROMOSOMES; do
#     echo "  Extracting $chr..."
#     samtools faidx "$FULL_REF" "$chr" >> reference_mini.fa
# done

# Option B: Just the regions (smaller but requires coordinate remapping)
# This won't work directly with Sheriff because coordinates won't match

# We'll go with a middle ground: extract regions with LARGE padding (1MB)
# This keeps coordinate system intact while being much smaller than full chroms

echo "  Extracting regions with 1MB padding..."
rm -f reference_mini.fa

# Group regions by chromosome and extract extended regions
prev_chr=""
for chr in $CHROMOSOMES; do
    echo "  Processing $chr..."
    # Get min/max positions for this chromosome
    MIN_POS=$(awk -v c="$chr" '$1==c {print $2}' regions.bed | sort -n | head -1)
    MAX_POS=$(awk -v c="$chr" '$1==c {print $3}' regions.bed | sort -n | tail -1)

    # Add 1MB padding on each side
    START=$((MIN_POS - 1000000))
    END=$((MAX_POS + 1000000))
    [ $START -lt 1 ] && START=1

    echo "    Region: $chr:$START-$END"

    # Extract this region
    samtools faidx "$FULL_REF" "${chr}:${START}-${END}" >> reference_mini.fa
done

# Step 4: Extract relevant GTF entries
echo "[4/5] Extracting GTF annotations for regions..."
# Get genes that overlap with our BAM regions
awk 'NR==FNR {chroms[$1]=1; next}
     /^#/ {print; next}
     $1 in chroms {print}' \
    <(echo "$CHROMOSOMES" | tr ' ' '\n') "$FULL_GTF" > reference_mini.gtf

NUM_GTF_LINES=$(grep -v "^#" reference_mini.gtf | wc -l)
echo "  GTF entries: $NUM_GTF_LINES"

# Step 5: Index the mini reference
echo "[5/5] Indexing mini reference..."
samtools faidx reference_mini.fa

# Summary
FA_SIZE=$(du -h reference_mini.fa | cut -f1)
GTF_SIZE=$(du -h reference_mini.gtf | cut -f1)
echo ""
echo "=============================================="
echo "Mini reference created!"
echo "=============================================="
echo "  reference_mini.fa: $FA_SIZE"
echo "  reference_mini.gtf: $GTF_SIZE"
echo "  Total: ~$(($(stat --printf="%s" reference_mini.fa) / 1000000))MB + $(($(stat --printf="%s" reference_mini.gtf) / 1000000))MB"
echo ""
echo "IMPORTANT: This mini reference extracts chromosome:start-end regions."
echo "Coordinates in the BAM will NOT match unless you use the full chromosomes."
echo ""
echo "For a working minimal test set, you need one of:"
echo "1. Full chromosomes (large, but coordinates match)"
echo "2. Subset BAM + remap coordinates (complex)"
echo "3. Skip reference-dependent tests for CI"
