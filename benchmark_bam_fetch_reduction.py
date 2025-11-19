#!/usr/bin/env python3
"""
Benchmark BAM fetch optimization - measures fetch reduction and performance impact.

This benchmark:
1. Simulates the optimization on realistic edit site distributions
2. Counts fetches before/after optimization
3. Estimates performance improvement
4. Analyzes different scenarios (clustered vs scattered edit sites)

Author: Sheriff Optimization Benchmark
Date: 2025-11-19
"""

import sys
import time
import random
from pathlib import Path
from collections import defaultdict, namedtuple

# Define EditSite namedtuple
EditSite = namedtuple("EditSite", ["chrom", "ref_pos"])


def simulate_edit_sites(n_sites, clustering="mixed"):
    """
    Simulate edit site distributions.

    clustering:
    - "clustered": Sites tend to cluster together (high overlap)
    - "scattered": Sites are well-separated (low overlap)
    - "mixed": Realistic mix of both
    """
    sites = []

    if clustering == "clustered":
        # Create clusters of edit sites
        n_clusters = max(1, n_sites // 20)  # ~20 sites per cluster
        for _ in range(n_clusters):
            cluster_chr = random.choice(["chr1", "chr2", "chr3", "chr4", "chr5"])
            cluster_center = random.randint(10000, 1000000)
            cluster_size = random.randint(10, 30)

            for _ in range(cluster_size):
                # Within 1000bp of cluster center
                pos = cluster_center + random.randint(-500, 500)
                sites.append(EditSite(cluster_chr, max(0, pos)))

    elif clustering == "scattered":
        # Well-separated sites
        for _ in range(n_sites):
            chr_name = random.choice(["chr1", "chr2", "chr3", "chr4", "chr5"])
            pos = random.randint(10000, 10000000)
            sites.append(EditSite(chr_name, pos))

    else:  # mixed
        # Realistic mix: some clusters, some scattered
        n_clustered = n_sites // 2
        n_scattered = n_sites - n_clustered

        # Clustered portion
        n_clusters = max(1, n_clustered // 15)
        for _ in range(n_clusters):
            cluster_chr = random.choice(["chr1", "chr2", "chr3", "chr4", "chr5"])
            cluster_center = random.randint(10000, 1000000)
            cluster_size = random.randint(10, 25)

            for _ in range(cluster_size):
                pos = cluster_center + random.randint(-300, 300)
                sites.append(EditSite(cluster_chr, max(0, pos)))

        # Scattered portion
        for _ in range(n_scattered):
            chr_name = random.choice(["chr1", "chr2", "chr3", "chr4", "chr5"])
            pos = random.randint(10000, 10000000)
            sites.append(EditSite(chr_name, pos))

    return sites


def count_naive_fetches(edit_sites):
    """Count fetches with naive approach (one per site)."""
    return len(edit_sites)


def count_batched_fetches(edit_sites, dist=500):
    """Count fetches with batching optimization."""
    # Group by chromosome
    edit_sites_by_chr = defaultdict(list)
    for edit_site in edit_sites:
        edit_sites_by_chr[edit_site.chrom].append(edit_site)

    # Sort by position
    for chr_name in edit_sites_by_chr:
        edit_sites_by_chr[chr_name].sort(key=lambda x: x.ref_pos)

    # Count batches
    batch_threshold = 2 * dist
    total_batches = 0

    for chr_name, chr_edit_sites in edit_sites_by_chr.items():
        i = 0
        while i < len(chr_edit_sites):
            batch_start_idx = i
            batch_end_idx = i
            batch_region_end = chr_edit_sites[i].ref_pos

            # Extend batch
            while (batch_end_idx + 1 < len(chr_edit_sites) and
                   chr_edit_sites[batch_end_idx + 1].ref_pos <= batch_region_end + batch_threshold):
                batch_end_idx += 1
                batch_region_end = chr_edit_sites[batch_end_idx].ref_pos

            total_batches += 1
            i = batch_end_idx + 1

    return total_batches


def benchmark_scenario(n_sites, clustering, dist=500):
    """Benchmark a specific scenario."""
    print(f"\n{'='*80}")
    print(f"Scenario: {n_sites} edit sites, {clustering} distribution")
    print(f"Distance parameter: {dist}bp")
    print('='*80)

    # Generate edit sites
    edit_sites = simulate_edit_sites(n_sites, clustering)

    # Count fetches
    naive_fetches = count_naive_fetches(edit_sites)
    batched_fetches = count_batched_fetches(edit_sites, dist)

    # Calculate improvement
    fetches_saved = naive_fetches - batched_fetches
    reduction_pct = 100 * (1 - batched_fetches / naive_fetches)

    # Estimate time saved (assuming ~1ms per fetch overhead)
    fetch_overhead_ms = 1.0
    time_saved_ms = fetches_saved * fetch_overhead_ms

    print(f"\n📊 Results:")
    print(f"   Naive approach:    {naive_fetches:,} fetches")
    print(f"   Batched approach:  {batched_fetches:,} fetches")
    print(f"   Fetches saved:     {fetches_saved:,} ({reduction_pct:.1f}% reduction)")
    print(f"   Estimated time saved: {time_saved_ms:.1f}ms")

    return {
        'n_sites': n_sites,
        'clustering': clustering,
        'naive_fetches': naive_fetches,
        'batched_fetches': batched_fetches,
        'reduction_pct': reduction_pct,
        'time_saved_ms': time_saved_ms
    }


def main():
    """Run comprehensive benchmark."""
    print("\n" + "="*80)
    print("BAM Fetch Reduction Benchmark")
    print("="*80)
    print("\nThis benchmark estimates the fetch reduction from batching optimization")
    print("based on different edit site distributions.")

    results = []

    # Scenario 1: Small dataset (100 sites, mixed)
    results.append(benchmark_scenario(100, "mixed", dist=500))

    # Scenario 2: Medium dataset (1000 sites, mixed)
    results.append(benchmark_scenario(1000, "mixed", dist=500))

    # Scenario 3: Large dataset (10000 sites, mixed) - realistic
    results.append(benchmark_scenario(10000, "mixed", dist=500))

    # Scenario 4: Highly clustered (10000 sites)
    results.append(benchmark_scenario(10000, "clustered", dist=500))

    # Scenario 5: Highly scattered (10000 sites)
    results.append(benchmark_scenario(10000, "scattered", dist=500))

    # Scenario 6: Effect of distance parameter (10000 sites, mixed)
    print(f"\n{'='*80}")
    print("Distance Parameter Impact (10,000 sites, mixed distribution)")
    print('='*80)

    edit_sites_fixed = simulate_edit_sites(10000, "mixed")
    for dist in [140, 280, 500, 1000]:
        batched = count_batched_fetches(edit_sites_fixed, dist)
        reduction = 100 * (1 - batched / 10000)
        print(f"   dist={dist:4d}bp: {batched:5,} fetches ({reduction:5.1f}% reduction)")

    # Summary
    print(f"\n{'='*80}")
    print("Summary")
    print('='*80)

    print("\nKey Findings:")
    print(f"1. Mixed distribution (realistic): {results[2]['reduction_pct']:.1f}% fetch reduction")
    print(f"2. Clustered sites (best case):   {results[3]['reduction_pct']:.1f}% fetch reduction")
    print(f"3. Scattered sites (worst case):  {results[4]['reduction_pct']:.1f}% fetch reduction")

    avg_reduction = sum(r['reduction_pct'] for r in results[:3]) / 3
    avg_time_saved = sum(r['time_saved_ms'] for r in results[:3]) / 3

    print(f"\nExpected improvement on realistic data:")
    print(f"   Fetch reduction: {avg_reduction:.1f}%")
    print(f"   Time saved:      {avg_time_saved:.1f}ms per 10K edit sites")

    # Estimate total pipeline impact
    # From analysis: edit site processing is ~8% of runtime (lines 1172-1210)
    # Plus non-barcoded edit processing is ~5% of runtime (lines 411-491)
    # Total: ~13% of runtime affected

    fetch_speedup = 1 / (1 - avg_reduction / 100)
    step_speedup_pct = (fetch_speedup - 1) * 100

    print(f"\nEstimated pipeline impact:")
    print(f"   Edit site processing step: {fetch_speedup:.2f}x faster ({step_speedup_pct:.1f}% improvement)")
    print(f"   Total pipeline (13% of runtime affected): ~{step_speedup_pct * 0.13:.2f}% faster")

    print("\n✅ Benchmark complete!")


if __name__ == "__main__":
    random.seed(42)  # Reproducible results
    main()
