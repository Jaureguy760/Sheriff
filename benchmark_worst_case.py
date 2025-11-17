#!/usr/bin/env python3
"""
Worst-case benchmark for edit clustering

This tests the scenario where ALL edits require full alignment:
- Same chromosome
- Same position
- Similar sequences (requiring Smith-Waterman alignment)

This represents the actual bottleneck in edit clustering.
"""

import time
import statistics
import sheriff_rs


def create_worst_case_edits(n_edits=50):
    """
    Create worst-case test edits that require full alignment

    All edits:
    - Same chromosome (chr1)
    - Same position (1000)
    - Same orientation (forward)
    - Similar sequences (require alignment to distinguish)

    This forces the algorithm to perform O(n²) expensive alignments.
    """
    edits = []

    # Base T7 donor template sequence (realistic for Sheriff)
    base_sequence = "GGAGAGTATAGAATGGAGC"

    for i in range(n_edits):
        chrom = "chr1"  # All same chromosome
        ref_pos = 1000  # All same position
        ref_seq = "ATCG"
        forward = True  # All same orientation

        # Create variations of the T7 sequence
        # - Different lengths (15-25bp)
        # - Slight sequence variations (1-3 mismatches)
        # - Some homopolymer differences

        insert_len = 15 + (i % 10)  # Vary length
        alt_seq_insert = base_sequence[:insert_len]

        # Add variation to force alignment
        if i % 5 == 0:
            # Substitute one base
            alt_seq_insert = alt_seq_insert[:10] + "T" + alt_seq_insert[11:]
        elif i % 5 == 1:
            # Add homopolymer
            alt_seq_insert = alt_seq_insert[:10] + "AAA" + alt_seq_insert[10:]
        elif i % 5 == 2:
            # Different homopolymer
            alt_seq_insert = alt_seq_insert[:10] + "GGG" + alt_seq_insert[10:]
        # else: keep similar

        alt_seq = ref_seq + alt_seq_insert
        kmer_matches = [i] if i % 3 == 0 else [i, i+1]

        edits.append((chrom, ref_pos, ref_seq, alt_seq, forward, kmer_matches))

    return edits


def create_best_case_edits(n_edits=50):
    """
    Create best-case test edits (all different, early exit)

    All edits on different chromosomes - hits early exit path.
    """
    edits = []

    for i in range(n_edits):
        chrom = f"chr{i}"  # All different chromosomes
        ref_pos = 1000
        ref_seq = "ATCG"
        alt_seq = ref_seq + "GGAGAGTAT"
        forward = True
        kmer_matches = [i]

        edits.append((chrom, ref_pos, ref_seq, alt_seq, forward, kmer_matches))

    return edits


def benchmark(edits, name, n_runs=10):
    """Run benchmark and return results"""
    print(f"\n⏱️  Benchmarking: {name}")
    print(f"  Input: {len(edits)} edits")

    # Warmup
    result = sheriff_rs.get_longest_edits_rust(edits)

    # Timed runs
    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        result = sheriff_rs.get_longest_edits_rust(edits)
        end = time.perf_counter()
        times.append((end - start) * 1000)

    mean_time = statistics.mean(times)
    std_time = statistics.stdev(times)

    print(f"  Output: {len(result)} edits")
    print(f"  Mean: {mean_time:.2f} ms (±{std_time:.2f} ms)")

    return {
        'name': name,
        'n_input': len(edits),
        'n_output': len(result),
        'mean_ms': mean_time,
        'std_ms': std_time,
        'times': times
    }


def main():
    print("="*70)
    print("WORST-CASE vs BEST-CASE BENCHMARK")
    print("="*70)
    print("\nThis benchmark compares:")
    print("  WORST: All edits same chr/pos, require full alignment")
    print("  BEST:  All edits different chr, early exit")
    print()

    # Test with different sizes
    sizes = [10, 25, 50, 100]

    for n in sizes:
        print("\n" + "="*70)
        print(f"Testing with {n} edits")
        print("="*70)

        worst_edits = create_worst_case_edits(n)
        best_edits = create_best_case_edits(n)

        worst_result = benchmark(worst_edits, "WORST CASE (requires alignment)", n_runs=5)
        best_result = benchmark(best_edits, "BEST CASE (early exit)", n_runs=5)

        # Calculate comparisons
        n_comparisons = (n * (n - 1)) // 2

        print(f"\n📊 Analysis for {n} edits:")
        print(f"  Comparisons:      {n_comparisons:6d}")
        print(f"  Worst case:       {worst_result['mean_ms']:8.2f} ms")
        print(f"  Best case:        {best_result['mean_ms']:8.2f} ms")
        print(f"  Slowdown:         {worst_result['mean_ms']/best_result['mean_ms']:8.1f}x")
        print(f"  Time per comp:    {worst_result['mean_ms']/n_comparisons*1000:8.2f} µs (worst)")
        print(f"                    {best_result['mean_ms']/n_comparisons*1000:8.2f} µs (best)")

        # Extrapolate to larger sizes
        if n == 50:
            # Estimate for real workload
            print(f"\n🔮 Extrapolation:")
            for real_n in [100, 200, 500]:
                real_comparisons = (real_n * (real_n - 1)) // 2
                estimated_worst = (worst_result['mean_ms'] / n_comparisons) * real_comparisons
                estimated_best = (best_result['mean_ms'] / n_comparisons) * real_comparisons
                print(f"  {real_n} edits (worst): ~{estimated_worst:8.0f} ms = {estimated_worst/1000:.1f}s")

    print("\n" + "="*70)
    print("KEY INSIGHTS")
    print("="*70)
    print("""
1. BEST CASE (different chromosomes): ~0.2-0.5ms
   - Early exit works perfectly
   - Already optimal

2. WORST CASE (similar edits): ??ms
   - Requires full Smith-Waterman alignment
   - This is where optimization matters!

3. OPTIMIZATION TARGETS:
   - Add pre-filters BEFORE expensive alignment
   - Cache expensive computations
   - Early exit based on cheap checks (length, position, etc.)
    """)

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
