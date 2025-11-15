#!/usr/bin/env python3
"""
K-mer matching benchmark: Python vs Rust

Compares performance of:
1. Python (numpy + KmerMatcher) - baseline
2. Rust (sheriff_rs.match_kmer_rust) - optimized

Tests with realistic T7 barcode detection workload.
"""

import time
import numpy as np
import sys
from pathlib import Path

# Add Sheriff to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from sheriff.count_t7 import KmerMatcher, match_kmer

# Try to import Rust
try:
    import sheriff_rs
    HAS_RUST = True
except ImportError:
    HAS_RUST = False
    print("⚠️  Rust not available, skipping Rust benchmarks")


def generate_test_sequence(length, with_n=False):
    """Generate random DNA sequence for testing."""
    bases = ['A', 'C', 'G', 'T']
    seq = ''.join(np.random.choice(bases, size=length))

    if with_n:
        # Add some N bases randomly (10%)
        seq_list = list(seq)
        n_positions = np.random.choice(len(seq_list), size=length//10, replace=False)
        for pos in n_positions:
            seq_list[pos] = 'N'
        seq = ''.join(seq_list)

    return seq


def benchmark_python_kmer(sequences, k, whitelist_kmers, num_iterations):
    """Benchmark Python k-mer matching."""
    # Create KmerMatcher with whitelist
    matcher = KmerMatcher(k=k, sequences=whitelist_kmers)

    start = time.time()

    for _ in range(num_iterations):
        for seq in sequences:
            # Match k-mers with whitelist filtering
            matches = match_kmer(matcher, seq, output_kmer_hash=True)

    duration = time.time() - start

    return duration


def benchmark_rust_kmer(sequences, k, whitelist_hashes, num_iterations):
    """Benchmark Rust k-mer matching."""
    start = time.time()

    for _ in range(num_iterations):
        for seq in sequences:
            # Match k-mers with whitelist filtering
            matches = sheriff_rs.match_kmer_rust(
                seq, k, whitelist=whitelist_hashes, output_hash=True
            )

    duration = time.time() - start

    return duration


def verify_outputs_match(sequences, k, whitelist_kmers):
    """Verify that Python and Rust produce identical results."""
    # Python implementation
    matcher = KmerMatcher(k=k, sequences=whitelist_kmers)

    # Convert whitelist k-mers to hashes for Rust
    whitelist_hashes = [matcher.kmer_to_num(kmer) for kmer in whitelist_kmers]

    mismatches = []

    for i, seq in enumerate(sequences):
        # Python result
        py_matches = match_kmer(matcher, seq, output_kmer_hash=True)
        if py_matches is None:
            py_matches = tuple()
        py_set = set(py_matches)

        # Rust result
        rust_matches = sheriff_rs.match_kmer_rust(
            seq, k, whitelist=whitelist_hashes, output_hash=True
        )
        rust_set = set(rust_matches)

        if py_set != rust_set:
            mismatches.append({
                'seq_idx': i,
                'seq': seq[:50] + '...' if len(seq) > 50 else seq,
                'python': sorted(py_set),
                'rust': sorted(rust_set),
                'diff': (py_set - rust_set, rust_set - py_set)
            })

    return mismatches


def run_benchmark():
    """Run comprehensive k-mer benchmark."""
    print("=" * 80)
    print("K-mer Matching Benchmark: Python vs Rust")
    print("=" * 80)
    print()

    # Test configuration
    k = 13  # Typical k-mer size for T7 barcodes
    num_sequences = 1000
    seq_length = 50  # Typical indel sequence length
    num_iterations = 10

    # Generate whitelist (simulate T7 barcodes)
    # Typical: 20 T7 barcodes
    whitelist_kmers = []
    for _ in range(20):
        whitelist_kmers.append(generate_test_sequence(k))

    # Generate test sequences
    print(f"Configuration:")
    print(f"  K-mer size: {k}")
    print(f"  Sequences: {num_sequences}")
    print(f"  Sequence length: {seq_length}")
    print(f"  Whitelist size: {len(whitelist_kmers)} k-mers")
    print(f"  Iterations: {num_iterations}")
    print()

    sequences = [generate_test_sequence(seq_length, with_n=True) for _ in range(num_sequences)]

    # Verify outputs match
    if HAS_RUST:
        print("Verifying output consistency...")
        mismatches = verify_outputs_match(sequences[:100], k, whitelist_kmers)

        if mismatches:
            print(f"  ❌ Found {len(mismatches)} mismatches!")
            for mm in mismatches[:3]:  # Show first 3
                print(f"    Seq {mm['seq_idx']}: {mm['seq']}")
                print(f"      Python: {mm['python']}")
                print(f"      Rust: {mm['rust']}")
        else:
            print("  ✅ Python and Rust produce identical results!")
        print()

    # Benchmark Python
    print("[1/2] Benchmarking Python (numpy + KmerMatcher)...")
    print("-" * 80)

    py_duration = benchmark_python_kmer(sequences, k, whitelist_kmers, num_iterations)

    total_matches = num_sequences * num_iterations
    py_throughput = total_matches / py_duration

    print(f"  Total matches: {total_matches:,}")
    print(f"  Duration: {py_duration:.3f}s")
    print(f"  Throughput: {py_throughput:,.0f} matches/sec")
    print()

    if not HAS_RUST:
        print("Rust not available - skipping Rust benchmark")
        return

    # Benchmark Rust
    print("[2/2] Benchmarking Rust (match_kmer_rust)...")
    print("-" * 80)

    # Create KmerMatcher to convert whitelist to hashes
    matcher = KmerMatcher(k=k, sequences=whitelist_kmers)
    whitelist_hashes = [matcher.kmer_to_num(kmer) for kmer in whitelist_kmers]

    rust_duration = benchmark_rust_kmer(sequences, k, whitelist_hashes, num_iterations)

    rust_throughput = total_matches / rust_duration

    print(f"  Total matches: {total_matches:,}")
    print(f"  Duration: {rust_duration:.3f}s")
    print(f"  Throughput: {rust_throughput:,.0f} matches/sec")
    print()

    # Summary
    print("=" * 80)
    print("PERFORMANCE SUMMARY")
    print("=" * 80)
    print()

    speedup = py_duration / rust_duration

    print(f"{'Implementation':<20} {'Time (s)':<15} {'Throughput (matches/s)':<30} {'Speedup':<15}")
    print("-" * 80)
    print(f"{'Python':<20} {py_duration:<15.3f} {py_throughput:<30,.0f} {'1.0x (baseline)':<15}")
    print(f"{'Rust':<20} {rust_duration:<15.3f} {rust_throughput:<30,.0f} {f'{speedup:.1f}x':<15}")
    print()

    print(f"🚀 Rust is {speedup:.1f}x faster than Python for k-mer matching!")
    print()

    # Realistic workload projection
    print("Projected impact on Sheriff pipeline (937M reads):")
    print(f"  If k-mer matching takes ~45 min in Python:")
    print(f"  With Rust: ~{45/speedup:.0f} min ({speedup:.1f}x speedup)")
    print(f"  Time saved: ~{45 - 45/speedup:.0f} minutes")
    print()


if __name__ == '__main__':
    np.random.seed(42)  # Reproducible results
    run_benchmark()
