#!/usr/bin/env python3
"""
Benchmark Rust vs Python on REAL Sheriff data
Extracts actual sequences and UMIs from the example BAM file
"""

import time
import sys
import pysam
from collections import defaultdict

sys.path.insert(0, '/home/user/Sheriff')

from sheriff.count_t7 import KmerMatcher, match_kmer as python_match_kmer
from sheriff.helpers import deduplicate_umis
import sheriff_rs

print("=" * 70)
print("REAL DATA Benchmark: Rust vs Python")
print("=" * 70)
print()

# ==============================================================================
# PART 1: Extract real data from BAM file
# ==============================================================================
print("📂 Loading real data from BAM file...")
bam_file = "/home/user/Sheriff/example_data/barcode_headAligned_anno.sorted.edit_regions_200kb.bam"

sequences = []
umis_by_cell = defaultdict(set)
read_count = 0
max_reads = 10000  # Sample first 10000 reads

with pysam.AlignmentFile(bam_file, "rb") as bam:
    for read in bam:
        if read_count >= max_reads:
            break

        # Extract sequence (skip if contains N)
        if read.query_sequence and 'N' not in read.query_sequence:
            sequences.append(read.query_sequence)

        # Extract UMI and cell barcode (if present)
        # Note: pN is the UMI tag, CB is the cell barcode tag
        if read.has_tag("pN") and read.has_tag("CB"):
            umi = read.get_tag("pN")
            cell_barcode = read.get_tag("CB")
            umis_by_cell[cell_barcode].add(umi)

        read_count += 1

print(f"  Extracted {len(sequences)} sequences from BAM")
print(f"  Extracted UMIs from {len(umis_by_cell)} cells")

# Count total UMIs across all cells
total_umis = sum(len(umis) for umis in umis_by_cell.values())
print(f"  Total unique UMIs: {total_umis}")
print()

# ==============================================================================
# PART 2: K-mer matching benchmark on real sequences
# ==============================================================================
print("📊 K-mer Matching Benchmark (Real BAM sequences)")
print("-" * 70)

k = 6
t7_barcode = "GGGAGAGTAT"

# Create matcher
bc_kmer_matcher = KmerMatcher(k, t7_barcode)

# Python k-mer matching
print("Testing Python k-mer matching...")
py_matches_all = []
py_start = time.perf_counter()

for seq in sequences:
    matches = python_match_kmer(bc_kmer_matcher, seq, output_kmer_hash=True)
    if matches:
        py_matches_all.append(matches)

py_time = (time.perf_counter() - py_start) * 1000
py_match_count = sum(len(m) for m in py_matches_all)

print(f"  Python: {py_time:.3f} ms  ({py_match_count} total k-mer matches)")

# Rust k-mer matching
print("Testing Rust k-mer matching...")

# Build whitelist from t7_barcode
whitelist = []
for i in range(len(t7_barcode) - k + 1):
    kmer = t7_barcode[i:i+k]
    hash_val = sheriff_rs.kmer_to_num(kmer)
    whitelist.append(hash_val)
whitelist = list(set(whitelist))

rust_matches_all = []
rust_start = time.perf_counter()

for seq in sequences:
    matches = sheriff_rs.match_kmer(seq, k, whitelist, output_hash=True)
    if matches:
        rust_matches_all.append(matches)

rust_time = (time.perf_counter() - rust_start) * 1000
rust_match_count = sum(len(m) for m in rust_matches_all)

print(f"  Rust:   {rust_time:.3f} ms  ({rust_match_count} total k-mer matches)")

kmer_speedup = py_time / rust_time if rust_time > 0 else 0
print(f"  Speedup: {kmer_speedup:.2f}x faster! 🚀")
print()

# ==============================================================================
# PART 3: UMI deduplication benchmark on real UMIs
# ==============================================================================
print("📊 UMI Deduplication Benchmark (Real BAM UMIs)")
print("-" * 70)

# Test on the largest cell (most UMIs)
largest_cell = max(umis_by_cell.items(), key=lambda x: len(x[1]))
cell_barcode, cell_umis = largest_cell
cell_umis_list = list(cell_umis)

print(f"  Testing cell: {cell_barcode}")
print(f"  UMIs in cell: {len(cell_umis_list)}")

# Python UMI deduplication
print("Testing Python UMI deduplication...")
py_start = time.perf_counter()
py_unique = deduplicate_umis(set(cell_umis_list))
py_time = (time.perf_counter() - py_start) * 1000

print(f"  Python: {py_time:.3f} ms  ({len(py_unique)} unique groups)")

# Rust UMI deduplication
print("Testing Rust UMI deduplication...")
rust_start = time.perf_counter()
rust_unique = sheriff_rs.deduplicate_umis(cell_umis_list, threshold=1)
rust_time = (time.perf_counter() - rust_start) * 1000

print(f"  Rust:   {rust_time:.3f} ms  ({rust_unique} unique groups)")

umi_speedup = py_time / rust_time if rust_time > 0 else 0
print(f"  Speedup: {umi_speedup:.2f}x faster! 🚀")
print()

# ==============================================================================
# PART 4: Multi-cell UMI benchmark
# ==============================================================================
print("📊 Multi-Cell UMI Deduplication Benchmark")
print("-" * 70)

# Test on multiple cells
test_cells = list(umis_by_cell.items())[:10]  # First 10 cells
total_test_umis = sum(len(umis) for _, umis in test_cells)

print(f"  Testing {len(test_cells)} cells")
print(f"  Total UMIs: {total_test_umis}")

# Python multi-cell
py_start = time.perf_counter()
py_results = []
for cell_bc, umis in test_cells:
    unique = deduplicate_umis(umis)
    py_results.append(len(unique))
py_time = (time.perf_counter() - py_start) * 1000

print(f"  Python: {py_time:.3f} ms  ({sum(py_results)} total unique)")

# Rust multi-cell
rust_start = time.perf_counter()
rust_results = []
for cell_bc, umis in test_cells:
    unique = sheriff_rs.deduplicate_umis(list(umis), threshold=1)
    rust_results.append(unique)
rust_time = (time.perf_counter() - rust_start) * 1000

print(f"  Rust:   {rust_time:.3f} ms  ({sum(rust_results)} total unique)")

multi_speedup = py_time / rust_time if rust_time > 0 else 0
print(f"  Speedup: {multi_speedup:.2f}x faster! 🚀")
print()

# ==============================================================================
# Summary
# ==============================================================================
print("=" * 70)
print("SUMMARY: Real Data Performance")
print("=" * 70)
print(f"K-mer Matching:       {kmer_speedup:6.2f}x speedup")
print(f"UMI Dedup (1 cell):   {umi_speedup:6.2f}x speedup")
print(f"UMI Dedup (10 cells): {multi_speedup:6.2f}x speedup")
print()
print("✅ These are the REAL speedups on actual Sheriff data!")
print("=" * 70)
