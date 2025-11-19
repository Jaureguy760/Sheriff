#!/usr/bin/env python3
"""
Analyze UMI count distribution in Sheriff BAM files

This script extracts per-cell UMI counts from a BAM file to understand:
1. Distribution of UMI counts across cells
2. Percentage of cells that would benefit from BK-tree (>50 UMIs)
3. Expected performance impact of using BK-tree
"""

import pysam
from collections import defaultdict
import numpy as np

def analyze_umi_distribution(bam_file):
    """Analyze UMI count distribution per cell"""

    print(f"Analyzing BAM file: {bam_file}")
    print("=" * 80)

    # Count UMIs per cell
    cell_umis = defaultdict(set)

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for i, read in enumerate(bam):
            if i % 100000 == 0:
                print(f"Processed {i:,} reads, found {len(cell_umis):,} cells...")

            # Extract cell barcode and UMI
            # Sheriff BAM files use "pN" for UMI instead of "UB"
            if not read.has_tag("CB") or not read.has_tag("pN"):
                continue

            cell_barcode = read.get_tag("CB")
            umi = read.get_tag("pN")

            cell_umis[cell_barcode].add(umi)

    print(f"\nTotal reads processed: {i+1:,}")
    print(f"Total cells found: {len(cell_umis):,}")

    # Calculate statistics
    umi_counts = [len(umis) for umis in cell_umis.values()]
    umi_counts_array = np.array(umi_counts)

    print("\n" + "=" * 80)
    print("UMI COUNT DISTRIBUTION")
    print("=" * 80)

    print(f"\nBasic Statistics:")
    print(f"  Min UMIs per cell:     {np.min(umi_counts_array)}")
    print(f"  Max UMIs per cell:     {np.max(umi_counts_array)}")
    print(f"  Mean UMIs per cell:    {np.mean(umi_counts_array):.2f}")
    print(f"  Median UMIs per cell:  {np.median(umi_counts_array):.2f}")
    print(f"  Std Dev:               {np.std(umi_counts_array):.2f}")

    print(f"\nPercentiles:")
    for p in [10, 25, 50, 75, 90, 95, 99]:
        print(f"  {p}th percentile:       {np.percentile(umi_counts_array, p):.0f}")

    # Analyze by size buckets
    print(f"\nUMI Count Buckets:")
    buckets = [
        (0, 10, "Very small (0-10)"),
        (11, 25, "Small (11-25)"),
        (26, 50, "Medium (26-50)"),
        (51, 100, "Large (51-100)"),
        (101, 200, "Very large (101-200)"),
        (201, float('inf'), "Huge (>200)")
    ]

    for min_umis, max_umis, label in buckets:
        count = sum(1 for n in umi_counts if min_umis <= n <= max_umis)
        pct = 100 * count / len(umi_counts)
        print(f"  {label:20s}: {count:6,} cells ({pct:5.2f}%)")

    # BK-tree analysis
    print("\n" + "=" * 80)
    print("BK-TREE PERFORMANCE ANALYSIS")
    print("=" * 80)

    crossover = 50
    cells_above_crossover = sum(1 for n in umi_counts if n > crossover)
    pct_above = 100 * cells_above_crossover / len(umi_counts)

    print(f"\nCrossover point (BK-tree beats brute force): {crossover} UMIs")
    print(f"  Cells with >{crossover} UMIs:   {cells_above_crossover:,} ({pct_above:.2f}%)")
    print(f"  Cells with <={crossover} UMIs:  {len(umi_counts) - cells_above_crossover:,} ({100-pct_above:.2f}%)")

    # Estimate speedup potential
    print(f"\nEstimated Performance Impact:")

    # For cells > crossover, BK-tree is ~1.5-4x faster depending on size
    # For cells <= crossover, no change (adaptive algorithm uses brute force)

    total_complexity_bf = sum(n * n for n in umi_counts)  # O(n²) for all
    total_complexity_bk = sum(n * n if n <= crossover else n * np.log2(n) * 10 for n in umi_counts)
    # Note: multiply by 10 for log(n) because constant factors

    speedup_ratio = total_complexity_bf / total_complexity_bk

    print(f"  Theoretical complexity reduction: {speedup_ratio:.2f}x")
    print(f"  Expected real-world speedup:      {min(1.5, speedup_ratio):.2f}x")
    print(f"    (limited by constant factors and overhead)")

    # Recommendation
    print("\n" + "=" * 80)
    print("RECOMMENDATION")
    print("=" * 80)

    if pct_above < 5:
        print(f"\nOnly {pct_above:.2f}% of cells would benefit from BK-tree.")
        print("RECOMMENDATION: Stick with brute force (simpler, less overhead)")
    elif pct_above < 20:
        print(f"\n{pct_above:.2f}% of cells would benefit from BK-tree.")
        print("RECOMMENDATION: Use adaptive algorithm (default: crossover={crossover})")
        print("  - Brute force for most cells (simple, fast)")
        print("  - BK-tree for large cells (better asymptotic complexity)")
    else:
        print(f"\n{pct_above:.2f}% of cells would benefit from BK-tree.")
        print("RECOMMENDATION: Use adaptive algorithm (consider lowering crossover)")
        print(f"  - Current crossover: {crossover} UMIs")
        print(f"  - Suggested crossover: {int(np.percentile(umi_counts_array, 70))} UMIs (70th percentile)")

    return umi_counts

if __name__ == "__main__":
    bam_file = "/home/user/Sheriff/example_data/barcode_headAligned_anno.sorted.edit_regions_200kb.bam"
    umi_counts = analyze_umi_distribution(bam_file)
