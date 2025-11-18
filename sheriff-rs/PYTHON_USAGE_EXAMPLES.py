#!/usr/bin/env python3
"""
Sheriff-rs Python Bindings - Usage Examples

This file contains practical examples showing how to use the sheriff-rs
Python bindings in real-world scenarios.

Prerequisites:
    pip install maturin
    cd sheriff-rs
    maturin develop --release --features python
"""

import sheriff_rs

# ============================================================================
# Example 1: Basic K-mer Hashing
# ============================================================================

def example_kmer_hashing():
    """Convert k-mers to numeric hashes."""
    print("="*70)
    print("Example 1: K-mer Hashing")
    print("="*70)

    kmers = ["ACGT", "AAAA", "TTTT", "GCGC", "ATCGATCG"]

    for kmer in kmers:
        hash_val = sheriff_rs.kmer_to_num(kmer)
        print(f"  {kmer:10s} -> {hash_val:6d}")

    print()


# ============================================================================
# Example 2: K-mer Whitelist Matching
# ============================================================================

def example_kmer_matching():
    """Match k-mers against a whitelist."""
    print("="*70)
    print("Example 2: K-mer Whitelist Matching")
    print("="*70)

    # Define target k-mers
    target_kmers = ["ACGT", "CGTA", "GTAC", "TACG"]
    whitelist = [sheriff_rs.kmer_to_num(kmer) for kmer in target_kmers]

    # Sequence to scan
    sequence = "ACGTACGTACGT"

    print(f"Sequence: {sequence}")
    print(f"Targets:  {target_kmers}")
    print()

    # Find matches (as hashes)
    matches = sheriff_rs.match_kmer(sequence, 4, whitelist, output_hash=True)
    print(f"Matches (hashes): {matches}")

    # Find matches (as strings)
    matches_str = sheriff_rs.match_kmer(sequence, 4, whitelist, output_hash=False)
    print(f"Matches (strings): {matches_str}")
    print()


# ============================================================================
# Example 3: K-mer Frequency Counting
# ============================================================================

def example_kmer_counting():
    """Count k-mer frequencies efficiently."""
    print("="*70)
    print("Example 3: K-mer Frequency Counting")
    print("="*70)

    # Create counter (reusable!)
    counter = sheriff_rs.KmerCounter(k=4)

    sequences = [
        "ACGTACGTACGT",
        "AAAATTTTCCCCGGGG",
        "ACGTACGTGGGGACGT"
    ]

    for seq in sequences:
        # Count k-mers
        freqs = counter.count_kmers(seq)

        # Find non-zero frequencies
        kmers_found = {}
        for kmer in ["ACGT", "AAAA", "TTTT", "CCCC", "GGGG", "CGTA", "GTAC"]:
            hash_val = sheriff_rs.kmer_to_num(kmer)
            count = freqs[hash_val]
            if count > 0:
                kmers_found[kmer] = count

        print(f"\nSequence: {seq}")
        print(f"K-mers found:")
        for kmer, count in sorted(kmers_found.items()):
            print(f"  {kmer}: {count}")

    print()


# ============================================================================
# Example 4: UMI Deduplication (Simple)
# ============================================================================

def example_umi_simple():
    """Simple UMI deduplication - just count unique groups."""
    print("="*70)
    print("Example 4: UMI Deduplication (Simple)")
    print("="*70)

    umis = [
        "ATCGATCG",  # Group 1
        "ATCGATCC",  # Group 1 (1 mismatch from above)
        "ATCGATCA",  # Group 1 (1 mismatch from first)
        "GCGCGCGC",  # Group 2
        "GCGCGCGC",  # Group 2 (duplicate)
        "GCGCGCGA",  # Group 2 (1 mismatch)
        "TTTTTTTT",  # Group 3 (unique)
    ]

    print(f"Total UMIs: {len(umis)}")
    print("UMIs:")
    for i, umi in enumerate(umis):
        print(f"  [{i}] {umi}")

    threshold = 1
    unique_count = sheriff_rs.deduplicate_umis(umis, threshold)

    print(f"\nThreshold: {threshold}")
    print(f"Unique UMI groups: {unique_count}")
    print()


# ============================================================================
# Example 5: UMI Deduplication (Detailed)
# ============================================================================

def example_umi_detailed():
    """Detailed UMI deduplication - get grouping information."""
    print("="*70)
    print("Example 5: UMI Deduplication (Detailed)")
    print("="*70)

    umis = [
        "ATCGATCG",
        "ATCGATCC",
        "GCGCGCGC",
        "GCGCGCGA",
        "TTTTTTTT",
    ]

    print(f"Total UMIs: {len(umis)}")
    print("UMIs:")
    for i, umi in enumerate(umis):
        print(f"  [{i}] {umi}")

    threshold = 1
    groups = sheriff_rs.deduplicate_umis_detailed(umis, threshold)

    print(f"\nThreshold: {threshold}")
    print(f"Found {len(groups)} groups:")
    for i, group in enumerate(groups):
        group_umis = [umis[idx] for idx in group]
        print(f"  Group {i+1}: indices {group}")
        for idx in group:
            print(f"    [{idx}] {umis[idx]}")

    print()


# ============================================================================
# Example 6: Hamming Distance Calculations
# ============================================================================

def example_hamming_distance():
    """Calculate Hamming distances between sequences."""
    print("="*70)
    print("Example 6: Hamming Distance")
    print("="*70)

    pairs = [
        ("ATCGATCG", "ATCGATCG"),  # Identical
        ("ATCGATCG", "ATCGATCC"),  # 1 mismatch
        ("ATCGATCG", "ATCGTTCC"),  # 2 mismatches
        ("ATCGATCG", "GCGCGCGC"),  # Many mismatches
        ("AAAAAAAA", "TTTTTTTT"),  # All different
    ]

    for seq_a, seq_b in pairs:
        dist = sheriff_rs.hamming_distance(seq_a, seq_b)
        match = "identical" if dist == 0 else f"{dist} mismatch{'es' if dist > 1 else ''}"
        print(f"  {seq_a} vs {seq_b}: {match}")

    print()


# ============================================================================
# Example 7: Performance Comparison (K-mers)
# ============================================================================

def example_performance_kmers():
    """Compare performance of k-mer operations."""
    import time
    import random

    print("="*70)
    print("Example 7: K-mer Performance")
    print("="*70)

    # Generate test data
    random.seed(42)
    bases = ['A', 'C', 'G', 'T']
    sequence = ''.join(random.choice(bases) for _ in range(100000))

    k = 6
    whitelist = list(range(500))  # First 500 possible k-mer hashes

    print(f"Sequence length: {len(sequence):,} bp")
    print(f"K-mer length: {k}")
    print(f"Whitelist size: {len(whitelist)}")

    # Benchmark matching
    start = time.time()
    matches = sheriff_rs.match_kmer(sequence, k, whitelist, output_hash=True)
    elapsed = time.time() - start

    print(f"\nResults:")
    print(f"  Matches found: {len(matches):,}")
    print(f"  Time: {elapsed*1000:.2f} ms")
    print(f"  Speed: {len(sequence)/elapsed/1e6:.2f} Mbp/s")
    print()


# ============================================================================
# Example 8: Performance Comparison (UMIs)
# ============================================================================

def example_performance_umis():
    """Compare performance of UMI deduplication."""
    import time
    import random

    print("="*70)
    print("Example 8: UMI Deduplication Performance")
    print("="*70)

    # Generate test UMIs
    random.seed(42)
    bases = ['A', 'C', 'G', 'T']
    umi_length = 12
    num_umis = 5000

    umis = [''.join(random.choice(bases) for _ in range(umi_length))
            for _ in range(num_umis)]

    print(f"Number of UMIs: {num_umis:,}")
    print(f"UMI length: {umi_length}")

    # Benchmark deduplication
    start = time.time()
    unique_count = sheriff_rs.deduplicate_umis(umis, threshold=1)
    elapsed = time.time() - start

    print(f"\nResults:")
    print(f"  Unique groups: {unique_count:,}")
    print(f"  Time: {elapsed*1000:.2f} ms")
    print(f"  Speed: {num_umis/elapsed:.0f} UMIs/s")
    print()


# ============================================================================
# Example 9: Real-World Pipeline
# ============================================================================

def example_pipeline():
    """Realistic analysis pipeline using multiple functions."""
    print("="*70)
    print("Example 9: Complete Analysis Pipeline")
    print("="*70)

    # Step 1: Define barcode whitelist
    print("Step 1: Building k-mer whitelist...")
    barcode_sequences = ["ACGTAC", "GTACGT", "TACGTA", "CGTACG"]
    barcode_whitelist = [sheriff_rs.kmer_to_num(bc) for bc in barcode_sequences]
    print(f"  Whitelist size: {len(barcode_whitelist)} barcodes")

    # Step 2: Scan reads for barcodes
    print("\nStep 2: Scanning reads for barcodes...")
    reads = [
        "ACGTACGGGGTACGTATTTCGTACG",
        "GGGACGTACCCCTACGTAAAA",
        "TTTTACGTACGGGGGGG",
    ]

    k = 6
    for i, read in enumerate(reads):
        matches = sheriff_rs.match_kmer(read, k, barcode_whitelist, output_hash=False)
        print(f"  Read {i+1}: found {len(matches)} barcode(s): {matches}")

    # Step 3: Count k-mer frequencies
    print("\nStep 3: Analyzing k-mer composition...")
    counter = sheriff_rs.KmerCounter(k)

    for i, read in enumerate(reads):
        freqs = counter.count_kmers(read)
        total_kmers = sum(freqs)
        unique_kmers = sum(1 for f in freqs if f > 0)
        print(f"  Read {i+1}: {total_kmers} total, {unique_kmers} unique k-mers")

    # Step 4: Deduplicate UMIs
    print("\nStep 4: Deduplicating UMIs...")
    umis = [
        "ATCGATCG",
        "ATCGATCC",
        "ATCGATCA",
        "GCGCGCGC",
        "GCGCGCGA",
        "TTTTTTTT",
        "TTTTTTCT",
    ]

    unique_count = sheriff_rs.deduplicate_umis(umis, threshold=1)
    print(f"  Input UMIs: {len(umis)}")
    print(f"  Unique molecules: {unique_count}")

    # Step 5: Summary
    print("\nStep 5: Summary")
    print(f"  Total reads processed: {len(reads)}")
    print(f"  Barcodes detected: Yes")
    print(f"  UMIs deduplicated: {len(umis)} -> {unique_count}")
    print()


# ============================================================================
# Main
# ============================================================================

def main():
    """Run all examples."""
    examples = [
        example_kmer_hashing,
        example_kmer_matching,
        example_kmer_counting,
        example_umi_simple,
        example_umi_detailed,
        example_hamming_distance,
        example_performance_kmers,
        example_performance_umis,
        example_pipeline,
    ]

    print("\n")
    print("*" * 70)
    print("  Sheriff-rs Python Bindings - Usage Examples")
    print("*" * 70)
    print("\n")

    for example in examples:
        example()

    print("*" * 70)
    print("  All examples completed!")
    print("*" * 70)
    print()


if __name__ == "__main__":
    try:
        main()
    except ImportError as e:
        print("Error: sheriff_rs module not found!")
        print("Please build it first:")
        print("  cd sheriff-rs")
        print("  maturin develop --release --features python")
        import sys
        sys.exit(1)
