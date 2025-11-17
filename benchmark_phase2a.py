#!/usr/bin/env python3
"""
Phase 2A Specific Benchmark

This benchmark specifically targets the Phase 2A optimizations:
1. Cached collapsed homopolymer strings (benefits homopolymer-heavy comparisons)
2. Direct mismatch counting (benefits ALL alignment calls)

We compare:
- Baseline (Phase 1 performance from PHASE_1_FINAL_RESULTS.md)
- Phase 2A (current implementation)
"""

import time
import statistics
import sheriff_rs


def benchmark_homopolymer_heavy():
    """
    Benchmark with many homopolymer sequences
    This stresses Optimization #1 (homopolymer caching)
    """
    print("\n" + "="*70)
    print("BENCHMARK 1: Homopolymer-Heavy Dataset")
    print("="*70)
    print("Target: Optimization #1 (cached collapsed homopolymer strings)")
    print("Expected: 5-10% improvement over Phase 1")
    print()

    # Create 50 edits with homopolymers and varying lengths
    # These will cluster together due to homopolymer correction
    edits = []
    for i in range(50):
        # Create sequences with homopolymers that differ in length
        base_seq = "ATCGATAG"
        if i % 3 == 0:
            insert = base_seq + "AAAGGGCCC"  # 3 of each
        elif i % 3 == 1:
            insert = base_seq + "AAAAAAGGGGGGGCCCCCCC"  # 6 of each
        else:
            insert = base_seq + "AAAAAAAAAAGGGGGGGGGGGCCCCCCCCCCC"  # 10 of each

        edits.append((
            "chr1", 1000, "ATCG",
            insert + "ATCG",  # Forward edit
            True, [i]
        ))

    # Benchmark
    times = []
    for _ in range(20):
        start = time.perf_counter()
        result = sheriff_rs.get_longest_edits_rust(edits)
        times.append((time.perf_counter() - start) * 1000)

    mean_time = statistics.mean(times)
    median_time = statistics.median(times)
    std_time = statistics.stdev(times)

    print(f"Results:")
    print(f"  Mean:   {mean_time:.3f} ms")
    print(f"  Median: {median_time:.3f} ms")
    print(f"  Std:    {std_time:.3f} ms")
    print(f"  Output: {len(result)} canonical edits (from {len(edits)} input)")

    # Phase 1 baseline for comparison (from PHASE_1_FINAL_RESULTS.md)
    phase1_baseline = 4.16  # 50 edits worst case
    improvement = ((phase1_baseline - mean_time) / phase1_baseline) * 100

    print(f"\nComparison to Phase 1 Baseline:")
    print(f"  Phase 1: 4.16 ms (worst case, 50 edits)")
    print(f"  Phase 2A: {mean_time:.2f} ms")
    print(f"  Improvement: {improvement:+.1f}%")

    return mean_time


def benchmark_alignment_heavy():
    """
    Benchmark with sequences requiring many alignments
    This stresses Optimization #2 (direct mismatch counting)
    """
    print("\n" + "="*70)
    print("BENCHMARK 2: Alignment-Heavy Dataset")
    print("="*70)
    print("Target: Optimization #2 (direct mismatch counting from alignment)")
    print("Expected: 5-10% improvement over Phase 1")
    print()

    # Create 50 edits with varying sequences (no exact matches)
    # All require full alignment
    edits = []
    for i in range(50):
        # Create unique but similar sequences
        insert = "ATCGATAG" + ("GCTA" * (i % 5 + 1)) + ("TAGC" * (i % 3))
        edits.append((
            "chr1", 1000, "ATCG",
            insert + "ATCG",
            True, [i]
        ))

    # Benchmark
    times = []
    for _ in range(20):
        start = time.perf_counter()
        result = sheriff_rs.get_longest_edits_rust(edits)
        times.append((time.perf_counter() - start) * 1000)

    mean_time = statistics.mean(times)
    median_time = statistics.median(times)
    std_time = statistics.stdev(times)

    print(f"Results:")
    print(f"  Mean:   {mean_time:.3f} ms")
    print(f"  Median: {median_time:.3f} ms")
    print(f"  Std:    {std_time:.3f} ms")
    print(f"  Output: {len(result)} canonical edits")

    phase1_baseline = 4.16
    improvement = ((phase1_baseline - mean_time) / phase1_baseline) * 100

    print(f"\nComparison to Phase 1 Baseline:")
    print(f"  Phase 1: 4.16 ms")
    print(f"  Phase 2A: {mean_time:.2f} ms")
    print(f"  Improvement: {improvement:+.1f}%")

    return mean_time


def benchmark_combined():
    """
    Benchmark with realistic mixed data
    Combines both optimizations
    """
    print("\n" + "="*70)
    print("BENCHMARK 3: Realistic Mixed Dataset")
    print("="*70)
    print("Target: Both optimizations combined")
    print("Expected: 10-20% total improvement over Phase 1")
    print()

    # Create 50 edits with mix of:
    # - Some with homopolymers
    # - Some requiring alignment
    # - Some exact matches (early exit)
    edits = []
    for i in range(50):
        if i % 4 == 0:
            # Homopolymer sequences
            insert = "ATCG" + "AAAGGGCCC" + ("T" * (i % 3))
        elif i % 4 == 1:
            # Similar sequences (require alignment)
            insert = "ATCGATAG" + ("GCTA" * (i % 3 + 1))
        elif i % 4 == 2:
            # Exact duplicates (early exit)
            insert = "ATCGATAGGGAGAGTAT"
        else:
            # Mixed
            insert = "ATCG" + ("ATAT" * (i % 4))

        edits.append((
            "chr1", 1000, "ATCG",
            insert + "ATCG",
            True, [i]
        ))

    # Benchmark
    times = []
    for _ in range(20):
        start = time.perf_counter()
        result = sheriff_rs.get_longest_edits_rust(edits)
        times.append((time.perf_counter() - start) * 1000)

    mean_time = statistics.mean(times)
    median_time = statistics.median(times)
    std_time = statistics.stdev(times)

    print(f"Results:")
    print(f"  Mean:   {mean_time:.3f} ms")
    print(f"  Median: {median_time:.3f} ms")
    print(f"  Std:    {std_time:.3f} ms")
    print(f"  Output: {len(result)} canonical edits")

    phase1_baseline = 4.16
    improvement = ((phase1_baseline - mean_time) / phase1_baseline) * 100

    print(f"\nComparison to Phase 1 Baseline:")
    print(f"  Phase 1: 4.16 ms")
    print(f"  Phase 2A: {mean_time:.2f} ms")
    print(f"  Improvement: {improvement:+.1f}%")

    return mean_time


def main():
    print("╔" + "="*68 + "╗")
    print("║" + " "*18 + "PHASE 2A PERFORMANCE BENCHMARK" + " "*20 + "║")
    print("╚" + "="*68 + "╝")

    print("\nPhase 2A Optimizations:")
    print("  1. Cache collapsed homopolymer strings")
    print("  2. Count mismatches directly from alignment operations")
    print()
    print("Expected Combined Speedup: 1.1-1.2x (10-20% improvement)")

    # Run benchmarks
    bench1_time = benchmark_homopolymer_heavy()
    bench2_time = benchmark_alignment_heavy()
    bench3_time = benchmark_combined()

    # Overall summary
    print("\n" + "="*70)
    print("PHASE 2A SUMMARY")
    print("="*70)

    avg_improvement = ((4.16 - ((bench1_time + bench2_time + bench3_time) / 3)) / 4.16) * 100

    print(f"\nAverage Performance:")
    print(f"  Phase 1 Baseline: 4.16 ms (50 edits)")
    print(f"  Phase 2A Average: {(bench1_time + bench2_time + bench3_time) / 3:.2f} ms")
    print(f"  Average Improvement: {avg_improvement:+.1f}%")

    if avg_improvement >= 10:
        print("\n✅ SUCCESS! Phase 2A achieved target improvement (10-20%)")
    elif avg_improvement >= 5:
        print("\n⚠️  PARTIAL SUCCESS - Phase 2A shows 5-10% improvement")
    else:
        print("\n⚠️  REVIEW NEEDED - Phase 2A improvement less than expected")

    print("\nNote: Improvements depend on data characteristics:")
    print("  - More homopolymers → More benefit from optimization #1")
    print("  - More alignments → More benefit from optimization #2")
    print("  - Real genomic data may show different patterns")

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
