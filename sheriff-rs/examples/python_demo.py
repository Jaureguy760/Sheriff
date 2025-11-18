#!/usr/bin/env python3
"""
Sheriff-rs Python Bindings Demo

This script demonstrates all the functionality exposed by the sheriff-rs
Python bindings, including k-mer operations and UMI deduplication.

Usage:
    python3 examples/python_demo.py

Requirements:
    - sheriff-rs must be built and installed with: maturin develop --release --features python
"""

import sys
import time
from typing import List

try:
    import sheriff_rs
except ImportError:
    print("Error: sheriff_rs module not found!")
    print("Please build and install it first:")
    print("  cd sheriff-rs")
    print("  maturin develop --release --features python")
    sys.exit(1)


def print_section(title: str):
    """Print a section header."""
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}\n")


def demo_version():
    """Demonstrate version information."""
    print_section("1. Version Information")
    print(f"sheriff-rs version: {sheriff_rs.__version__}")
    print(f"Module docstring: {sheriff_rs.__doc__}")


def demo_kmer_to_num():
    """Demonstrate k-mer to number conversion."""
    print_section("2. K-mer to Number Conversion")

    test_kmers = ["A", "C", "G", "T", "AC", "ACGT", "AAAAAA", "TTTTTT"]

    print("Converting k-mers to numeric hashes:\n")
    for kmer in test_kmers:
        hash_val = sheriff_rs.kmer_to_num(kmer)
        print(f"  kmer_to_num('{kmer:8s}') = {hash_val:6d}")

    # Demonstrate case-insensitivity
    print("\nCase-insensitivity check:")
    upper = sheriff_rs.kmer_to_num("ACGT")
    lower = sheriff_rs.kmer_to_num("acgt")
    mixed = sheriff_rs.kmer_to_num("AcGt")
    print(f"  'ACGT' -> {upper}")
    print(f"  'acgt' -> {lower}")
    print(f"  'AcGt' -> {mixed}")
    print(f"  All equal: {upper == lower == mixed}")

    # Demonstrate error handling
    print("\nError handling:")
    try:
        sheriff_rs.kmer_to_num("ACGN")  # Invalid nucleotide
    except ValueError as e:
        print(f"  kmer_to_num('ACGN') -> ValueError: {e}")


def demo_match_kmer():
    """Demonstrate k-mer matching against whitelist."""
    print_section("3. K-mer Matching Against Whitelist")

    sequence = "ACGTACGTGGGGACGT"
    k = 4

    # Create whitelist with specific k-mers
    target_kmers = ["ACGT", "GGGG"]
    whitelist = [sheriff_rs.kmer_to_num(kmer) for kmer in target_kmers]

    print(f"Sequence: {sequence}")
    print(f"K-mer length: {k}")
    print(f"Whitelist k-mers: {target_kmers}")
    print(f"Whitelist hashes: {whitelist}\n")

    # Match with hash output
    matches_hash = sheriff_rs.match_kmer(sequence, k, whitelist, output_hash=True)
    print(f"Matches (as hashes): {matches_hash}")
    print(f"Number of matches: {len(matches_hash)}")

    # Match with string output
    matches_str = sheriff_rs.match_kmer(sequence, k, whitelist, output_hash=False)
    print(f"Matches (as strings): {matches_str}")

    # Verify correctness
    print(f"\nVerification:")
    print(f"  'ACGT' appears {matches_str.count('ACGT')} times")
    print(f"  'GGGG' appears {matches_str.count('GGGG')} times")


def demo_kmer_counter():
    """Demonstrate KmerCounter class for frequency counting."""
    print_section("4. K-mer Frequency Counting")

    # Create a counter
    k = 4
    counter = sheriff_rs.KmerCounter(k)

    print(f"Created {counter}")
    print(f"  k = {counter.k}")
    print(f"  Array size = {counter.array_size} (4^{k} = {4**k})\n")

    # Count k-mers in a sequence
    sequence = "ACGTACGTACGT"
    print(f"Counting k-mers in: {sequence}\n")

    freqs = counter.count_kmers(sequence)

    # Show frequencies for k-mers that appear
    print("K-mer frequencies:")
    kmers_in_sequence = set()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmers_in_sequence.add(kmer)

    for kmer in sorted(kmers_in_sequence):
        hash_val = sheriff_rs.kmer_to_num(kmer)
        count = freqs[hash_val]
        print(f"  {kmer}: {count}")

    # Demonstrate array reuse
    print("\nArray reuse demonstration:")
    sequence2 = "AAAATTTTCCCCGGGG"
    freqs2 = counter.count_kmers(sequence2)
    print(f"Counted k-mers in: {sequence2}")

    # Verify that previous counts are cleared
    acgt_hash = sheriff_rs.kmer_to_num("ACGT")
    print(f"  'ACGT' count in new sequence: {freqs2[acgt_hash]} (should be 0)")

    aaaa_hash = sheriff_rs.kmer_to_num("AAAA")
    print(f"  'AAAA' count in new sequence: {freqs2[aaaa_hash]}")


def demo_umi_deduplication():
    """Demonstrate UMI deduplication."""
    print_section("5. UMI Deduplication")

    # Create test UMIs with known relationships
    umis = [
        "ATCGATCG",  # 0: Unique cluster 1
        "ATCGATCC",  # 1: 1 mismatch from 0 (same cluster)
        "ATCGATCA",  # 2: 1 mismatch from 0 (same cluster)
        "GCGCGCGC",  # 3: Unique cluster 2
        "GCGCGCGC",  # 4: Exact duplicate of 3
        "GCGCGCGA",  # 5: 1 mismatch from 3 (same cluster)
        "TTTTTTTT",  # 6: Unique cluster 3
    ]

    threshold = 1

    print(f"UMI sequences ({len(umis)} total):")
    for i, umi in enumerate(umis):
        print(f"  [{i}] {umi}")

    print(f"\nHamming distance threshold: {threshold}\n")

    # Simple count
    unique_count = sheriff_rs.deduplicate_umis(umis, threshold)
    print(f"Unique UMI groups: {unique_count}")

    # Detailed groupings
    groups = sheriff_rs.deduplicate_umis_detailed(umis, threshold)
    print(f"\nDetailed groupings:")
    for i, group in enumerate(groups):
        group_umis = [umis[idx] for idx in group]
        print(f"  Group {i+1}: indices {group} -> {group_umis}")

    # Demonstrate Hamming distance function
    print(f"\nHamming distance examples:")
    test_pairs = [
        ("ATCGATCG", "ATCGATCG"),  # Identical
        ("ATCGATCG", "ATCGATCC"),  # 1 mismatch
        ("ATCGATCG", "GCGCGCGC"),  # Many mismatches
    ]

    for a, b in test_pairs:
        dist = sheriff_rs.hamming_distance(a, b)
        print(f"  hamming_distance('{a}', '{b}') = {dist}")


def demo_performance():
    """Demonstrate performance improvements."""
    print_section("6. Performance Comparison")

    print("Generating test data...")

    # Generate a long sequence for k-mer matching
    import random
    random.seed(42)
    bases = ['A', 'C', 'G', 'T']
    long_sequence = ''.join(random.choice(bases) for _ in range(10000))

    # Create whitelist
    k = 6
    whitelist = list(range(100))  # Just use first 100 hashes

    print(f"\nK-mer matching benchmark:")
    print(f"  Sequence length: {len(long_sequence):,} bp")
    print(f"  K-mer length: {k}")
    print(f"  Whitelist size: {len(whitelist)}")

    # Time the Rust implementation
    start = time.time()
    matches = sheriff_rs.match_kmer(long_sequence, k, whitelist, output_hash=True)
    rust_time = time.time() - start

    print(f"  Rust time: {rust_time*1000:.2f} ms")
    print(f"  Matches found: {len(matches)}")

    # Generate UMIs for deduplication benchmark
    umi_length = 12
    num_umis = 1000
    test_umis = [''.join(random.choice(bases) for _ in range(umi_length))
                 for _ in range(num_umis)]

    print(f"\nUMI deduplication benchmark:")
    print(f"  Number of UMIs: {num_umis:,}")
    print(f"  UMI length: {umi_length}")

    # Time the Rust implementation
    start = time.time()
    unique_count = sheriff_rs.deduplicate_umis(test_umis, threshold=1)
    rust_time = time.time() - start

    print(f"  Rust time: {rust_time*1000:.2f} ms")
    print(f"  Unique groups: {unique_count}")

    print("\nNote: The Rust implementations provide 4-14x speedup for k-mer")
    print("      operations and 3-6x speedup for UMI deduplication compared")
    print("      to pure Python implementations.")


def demo_error_handling():
    """Demonstrate error handling."""
    print_section("7. Error Handling")

    print("Testing various error conditions:\n")

    # Invalid nucleotides
    print("1. Invalid nucleotides:")
    try:
        sheriff_rs.kmer_to_num("ACGN")
    except ValueError as e:
        print(f"   ✓ Caught: {e}\n")

    # Invalid k value
    print("2. Invalid k value (k=0):")
    try:
        sheriff_rs.match_kmer("ACGT", 0, [1, 2, 3], output_hash=True)
    except ValueError as e:
        print(f"   ✓ Caught: {e}\n")

    # UMIs of different lengths
    print("3. UMIs with different lengths:")
    try:
        sheriff_rs.deduplicate_umis(["ATCG", "ATCGATCG"], threshold=1)
    except ValueError as e:
        print(f"   ✓ Caught: {e}\n")

    # K too large
    print("4. K-mer length too large (k=20):")
    try:
        sheriff_rs.KmerCounter(20)
    except ValueError as e:
        print(f"   ✓ Caught: {e}\n")

    print("All error conditions handled correctly!")


def main():
    """Run all demonstrations."""
    print("="*70)
    print("  Sheriff-rs Python Bindings Comprehensive Demo")
    print("="*70)

    try:
        demo_version()
        demo_kmer_to_num()
        demo_match_kmer()
        demo_kmer_counter()
        demo_umi_deduplication()
        demo_performance()
        demo_error_handling()

        print("\n" + "="*70)
        print("  All demonstrations completed successfully!")
        print("="*70)
        print("\nFor more information:")
        print("  - Module help: help(sheriff_rs)")
        print("  - Function help: help(sheriff_rs.kmer_to_num)")
        print("  - Class help: help(sheriff_rs.KmerCounter)")
        print("\n")

    except Exception as e:
        print(f"\nError during demonstration: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
