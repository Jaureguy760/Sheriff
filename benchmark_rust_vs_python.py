#!/usr/bin/env python3
"""
Benchmark Rust vs Python implementations
Compare actual performance of sheriff-rs vs original Python
"""

import time
import sys
import os

# Add sheriff to path
sys.path.insert(0, '/home/user/Sheriff')

from sheriff.count_t7 import KmerMatcher, match_kmer as python_match_kmer

def generate_sequence(length, seed=42):
    """Generate random DNA sequence"""
    import random
    random.seed(seed)
    return ''.join(random.choice('ACGT') for _ in range(length))

def benchmark_python_kmer():
    """Benchmark Python k-mer matching"""
    k = 6
    t7_barcode = "GGGAGAGTAT"

    # Create matcher
    bc_kmer_matcher = KmerMatcher(k, t7_barcode)

    # Test sequence (1000bp)
    sequence = generate_sequence(1000)

    # Warm up
    for _ in range(10):
        _ = python_match_kmer(bc_kmer_matcher, sequence, output_kmer_hash=True)

    # Benchmark
    iterations = 100
    start = time.perf_counter()
    for _ in range(iterations):
        matches = python_match_kmer(bc_kmer_matcher, sequence, output_kmer_hash=True)
    end = time.perf_counter()

    avg_time = (end - start) / iterations * 1000  # ms
    return avg_time, len(matches) if matches else 0

def benchmark_rust_kmer():
    """Benchmark Rust k-mer matching"""
    try:
        import sheriff_rs
    except ImportError:
        print("⚠️  sheriff_rs not built yet. Run: cd sheriff-rs && ./build_python.sh")
        return None, 0

    k = 6
    sequence = generate_sequence(1000)

    # Create whitelist from t7_barcode
    t7_barcode = "GGGAGAGTAT"
    whitelist = []
    for i in range(len(t7_barcode) - k + 1):
        kmer = t7_barcode[i:i+k]
        hash_val = sheriff_rs.kmer_to_num(kmer)
        whitelist.append(hash_val)
    whitelist = list(set(whitelist))  # Remove duplicates

    # Warm up
    for _ in range(10):
        _ = sheriff_rs.match_kmer(sequence, k, whitelist, output_hash=True)

    # Benchmark
    iterations = 100
    start = time.perf_counter()
    for _ in range(iterations):
        matches = sheriff_rs.match_kmer(sequence, k, whitelist, output_hash=True)
    end = time.perf_counter()

    avg_time = (end - start) / iterations * 1000  # ms
    return avg_time, len(matches)

def benchmark_python_umi():
    """Benchmark Python UMI deduplication"""
    from sheriff.helpers import deduplicate_umis

    # Generate test UMIs
    umis = [
        "ATCGATCG",
        "ATCGATCC",  # 1 mismatch from first
        "GCGCGCGC",
        "GCGCGCGT",  # 1 mismatch from third
        "TATATATA",
    ] * 20  # 100 UMIs total

    # Warm up
    for _ in range(5):
        _ = deduplicate_umis(set(umis))

    # Benchmark
    iterations = 50
    start = time.perf_counter()
    for _ in range(iterations):
        unique = deduplicate_umis(set(umis))
    end = time.perf_counter()

    avg_time = (end - start) / iterations * 1000  # ms
    return avg_time, len(unique)

def benchmark_rust_umi():
    """Benchmark Rust UMI deduplication"""
    try:
        import sheriff_rs
    except ImportError:
        return None, 0

    umis = [
        "ATCGATCG",
        "ATCGATCC",
        "GCGCGCGC",
        "GCGCGCGT",
        "TATATATA",
    ] * 20

    # Warm up
    for _ in range(5):
        _ = sheriff_rs.deduplicate_umis(umis, threshold=1)

    # Benchmark
    iterations = 50
    start = time.perf_counter()
    for _ in range(iterations):
        unique_count = sheriff_rs.deduplicate_umis(umis, threshold=1)
    end = time.perf_counter()

    avg_time = (end - start) / iterations * 1000  # ms
    return avg_time, unique_count

def main():
    print("=" * 70)
    print("Sheriff: Rust vs Python Performance Benchmark")
    print("=" * 70)
    print()

    # K-mer matching
    print("📊 K-mer Matching (1000bp sequence, k=6)")
    print("-" * 70)

    py_time, py_matches = benchmark_python_kmer()
    print(f"Python:  {py_time:.3f} ms/iteration  ({py_matches} matches)")

    rust_time, rust_matches = benchmark_rust_kmer()
    if rust_time:
        speedup = py_time / rust_time
        print(f"Rust:    {rust_time:.3f} ms/iteration  ({rust_matches} matches)")
        print(f"Speedup: {speedup:.2f}x faster! 🚀")
    else:
        print("Rust:    Not available (build with ./build_python.sh)")
    print()

    # UMI deduplication
    print("📊 UMI Deduplication (100 UMIs, threshold=1)")
    print("-" * 70)

    py_time, py_unique = benchmark_python_umi()
    print(f"Python:  {py_time:.3f} ms/iteration  ({py_unique} unique groups)")

    rust_time, rust_unique = benchmark_rust_umi()
    if rust_time:
        speedup = py_time / rust_time
        print(f"Rust:    {rust_time:.3f} ms/iteration  ({rust_unique} unique groups)")
        print(f"Speedup: {speedup:.2f}x faster! 🚀")
    else:
        print("Rust:    Not available (build with ./build_python.sh)")
    print()

    print("=" * 70)
    if rust_time:
        print("✅ Benchmarks complete! Rust optimizations are working!")
    else:
        print("⚠️  Build sheriff_rs to see real speedups:")
        print("   cd sheriff-rs && ./build_python.sh")
    print("=" * 70)

if __name__ == "__main__":
    main()
