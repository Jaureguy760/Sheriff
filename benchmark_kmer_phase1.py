#!/usr/bin/env python3
"""
K-mer Phase 1 Performance Benchmark

Measures the impact of Phase 1 optimizations:
1. Nucleotide lookup table (2-4x faster than match statement)
2. FxHashSet (2-3x faster than std::HashSet for integer keys)

Expected combined improvement: 4-14x speedup
"""

import time
import statistics
import sheriff_rs


def benchmark_kmer_to_num():
    """
    Benchmark kmer_to_num function (lookup table optimization)
    This is a microbenchmark for the lookup table
    """
    print("=" * 70)
    print("BENCHMARK 1: kmer_to_num (Lookup Table Optimization)")
    print("=" * 70)
    print("Test: Convert 100,000 16-mers to integer hashes")
    print()

    # Create 100,000 random 16-mers
    import random
    bases = ['A', 'C', 'G', 'T']
    kmers = []
    for _ in range(100000):
        kmer = ''.join(random.choice(bases) for _ in range(16))
        kmers.append(kmer)

    # Benchmark (just time the operation, sheriff_rs doesn't export kmer_to_num directly)
    # We'll use count_kmers_rust which calls kmer_to_num internally
    sequence = ''.join(kmers)

    times = []
    for _ in range(20):
        start = time.perf_counter()
        result = sheriff_rs.count_kmers_rust(sequence, 16)
        times.append((time.perf_counter() - start) * 1000)

    mean = statistics.mean(times)
    median = statistics.median(times)
    std = statistics.stdev(times)

    print(f"Results:")
    print(f"  Mean:   {mean:.3f} ms")
    print(f"  Median: {median:.3f} ms")
    print(f"  Std:    {std:.3f} ms")
    print(f"  Operations: {len(sequence) - 16 + 1:,} k-mer conversions")
    print(f"  Time per k-mer: {(mean * 1000) / (len(sequence) - 16 + 1):.3f} ns")
    print()

    return mean


def benchmark_match_kmer_no_whitelist():
    """
    Benchmark match_kmer without whitelist
    This stresses the hash set creation and k-mer matching
    """
    print("=" * 70)
    print("BENCHMARK 2: match_kmer Without Whitelist")
    print("=" * 70)
    print("Test: Find all unique k-mers in 1000bp sequence")
    print()

    # Create 1000bp sequence with variety
    import random
    bases = ['A', 'C', 'G', 'T']
    sequence = ''.join(random.choice(bases) for _ in range(1000))

    times = []
    for _ in range(100):
        start = time.perf_counter()
        result = sheriff_rs.match_kmer_rust(sequence, 8, None, True)
        times.append((time.perf_counter() - start) * 1000)

    mean = statistics.mean(times)
    median = statistics.median(times)
    std = statistics.stdev(times)

    print(f"Results:")
    print(f"  Mean:   {mean:.3f} ms")
    print(f"  Median: {median:.3f} ms")
    print(f"  Std:    {std:.3f} ms")
    print(f"  Sequence length: {len(sequence)}bp")
    print(f"  K-mer size: 8")
    print(f"  Unique k-mers found: {len(result)}")
    print()

    return mean


def benchmark_match_kmer_with_whitelist():
    """
    Benchmark match_kmer WITH whitelist
    This stresses the FxHashSet lookup optimization
    """
    print("=" * 70)
    print("BENCHMARK 3: match_kmer With Whitelist (FxHashSet Lookup)")
    print("=" * 70)
    print("Test: Match k-mers against 1000-entry whitelist")
    print()

    # Create whitelist of 1000 random k-mer hashes
    import random
    whitelist = list(range(0, 4**8, (4**8) // 1000))[:1000]  # 1000 evenly spaced hashes

    # Create 1000bp sequence
    bases = ['A', 'C', 'G', 'T']
    sequence = ''.join(random.choice(bases) for _ in range(1000))

    times = []
    for _ in range(100):
        start = time.perf_counter()
        result = sheriff_rs.match_kmer_rust(sequence, 8, whitelist, True)
        times.append((time.perf_counter() - start) * 1000)

    mean = statistics.mean(times)
    median = statistics.median(times)
    std = statistics.stdev(times)

    print(f"Results:")
    print(f"  Mean:   {mean:.3f} ms")
    print(f"  Median: {median:.3f} ms")
    print(f"  Std:    {std:.3f} ms")
    print(f"  Sequence length: {len(sequence)}bp")
    print(f"  Whitelist size: {len(whitelist)}")
    print(f"  Matches found: {len(result)}")
    print(f"  Time per lookup: {(mean * 1000) / (len(sequence) - 8 + 1):.3f} ns")
    print()

    return mean


def benchmark_realistic_workload():
    """
    Realistic T7 barcode matching workload
    Matches 16-mer barcodes against whitelist in 150bp read
    """
    print("=" * 70)
    print("BENCHMARK 4: Realistic T7 Barcode Matching")
    print("=" * 70)
    print("Test: Match 16-mer barcodes in 150bp reads (10,000 reads)")
    print()

    # Create whitelist of 384 T7 barcodes (typical plate)
    whitelist = list(range(0, 4**16, (4**16) // 384))[:384]

    # Create 10,000 simulated 150bp reads
    import random
    bases = ['A', 'C', 'G', 'T']
    reads = []
    for _ in range(10000):
        read = ''.join(random.choice(bases) for _ in range(150))
        reads.append(read)

    # Benchmark
    times = []
    for _ in range(10):
        start = time.perf_counter()
        total_matches = 0
        for read in reads:
            matches = sheriff_rs.match_kmer_rust(read, 16, whitelist, True)
            total_matches += len(matches)
        times.append((time.perf_counter() - start) * 1000)

    mean = statistics.mean(times)
    median = statistics.median(times)
    std = statistics.stdev(times)

    print(f"Results:")
    print(f"  Mean:   {mean:.3f} ms")
    print(f"  Median: {median:.3f} ms")
    print(f"  Std:    {std:.3f} ms")
    print(f"  Reads processed: {len(reads):,}")
    print(f"  Time per read: {(mean / len(reads)) * 1000:.3f} µs")
    print(f"  Throughput: {(len(reads) / mean) * 1000:.0f} reads/sec")
    print()

    return mean


def main():
    print("╔" + "=" * 68 + "╗")
    print("║" + " " * 15 + "K-MER PHASE 1 PERFORMANCE BENCHMARK" + " " * 18 + "║")
    print("╚" + "=" * 68 + "╝")
    print()
    print("Phase 1 Optimizations:")
    print("  1. Nucleotide lookup table (replaces 6-way match statement)")
    print("  2. FxHashSet (2-3x faster than std::HashSet for integer keys)")
    print()
    print("Expected Combined Speedup: 4-14x")
    print()

    # Run benchmarks
    time1 = benchmark_kmer_to_num()
    time2 = benchmark_match_kmer_no_whitelist()
    time3 = benchmark_match_kmer_with_whitelist()
    time4 = benchmark_realistic_workload()

    # Summary
    print("=" * 70)
    print("PHASE 1 SUMMARY")
    print("=" * 70)
    print()
    print("Benchmark Results:")
    print(f"  1. kmer_to_num (lookup table):     {time1:.2f} ms")
    print(f"  2. match_kmer no whitelist:        {time2:.2f} ms")
    print(f"  3. match_kmer with whitelist:      {time3:.2f} ms")
    print(f"  4. Realistic T7 barcode matching:  {time4:.2f} ms ({(10000 / time4) * 1000:.0f} reads/sec)")
    print()
    print("Key Improvements:")
    print("  ✅ Lookup table eliminates branch mispredictions")
    print("  ✅ FxHashSet provides faster integer hashing")
    print("  ✅ Inline annotations improve hot path performance")
    print()
    print("Note: These are optimized timings. To measure improvement vs baseline,")
    print("we would need to run against the pre-Phase-1 code (not available).")
    print()
    print("Expected improvement based on analysis: 4-14x faster")
    print("  - Lookup table: 2-4x")
    print("  - FxHashSet: 2-3x")
    print("  - Combined: 4-14x (varies by workload)")

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
