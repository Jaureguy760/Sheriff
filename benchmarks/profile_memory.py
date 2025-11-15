#!/usr/bin/env python3
"""
Memory profiling for Sheriff using memory_profiler

Usage:
    python -m memory_profiler benchmarks/profile_memory.py

Or with line-by-line profiling:
    mprof run benchmarks/profile_memory.py
    mprof plot
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    from memory_profiler import profile
except ImportError:
    print("Installing memory_profiler...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "memory_profiler"])
    from memory_profiler import profile

from sheriff.count_t7 import KmerMatcher, match_kmer
from sheriff.helpers import deduplicate_umis, cell_umi_counts_FAST
import numpy as np


@profile
def profile_kmer_matcher():
    """Profile k-mer matcher memory usage"""
    k = 6
    t7_barcode = "GGGAGAGTAT"
    matcher = KmerMatcher(k, t7_barcode)

    # Simulate processing many sequences
    test_sequence = "ATCGATCGATCGGGGAGAGTATATCGATCGATCG" * 100

    results = []
    for i in range(1000):
        freq_array = np.zeros((4 ** k), dtype=np.uint8)
        kmers = [test_sequence[j:j+k] for j in range(len(test_sequence) - k + 1)]
        for kmer in kmers:
            try:
                num = matcher.kmer_to_num(kmer)
                results.append(num)
            except KeyError:
                pass

    return len(results)


@profile
def profile_umi_deduplication():
    """Profile UMI deduplication memory usage"""
    # Generate test data
    umi_length = 10
    n_cells = 100
    umis_per_cell = 100

    test_umis = []
    for _ in range(n_cells):
        umis = set()
        for _ in range(umis_per_cell):
            umi = ''.join(np.random.choice(['A', 'C', 'G', 'T'], umi_length))
            umis.add(umi)
        test_umis.append(umis)

    # Profile slow deduplication
    total = 0
    for umi_set in test_umis:
        deduped = deduplicate_umis(umi_set)
        total += len(deduped)

    return total


@profile
def profile_umi_deduplication_fast():
    """Profile fast UMI deduplication memory usage"""
    # Generate test data
    umi_length = 10
    n_cells = 100
    umis_per_cell = 100

    cell_bc_indexes = np.arange(n_cells, dtype=np.uint32)
    cell_umis = []

    for _ in range(n_cells):
        umis = []
        for _ in range(umis_per_cell):
            umi = ''.join(np.random.choice(['A', 'C', 'G', 'T'], umi_length))
            umis.append(umi)
        cell_umis.append(np.array(umis, dtype=f'<U{umi_length}'))

    # Profile fast deduplication
    counts = cell_umi_counts_FAST(n_cells, cell_bc_indexes, cell_umis)

    return counts.sum()


if __name__ == '__main__':
    print("Profiling k-mer matcher...")
    result1 = profile_kmer_matcher()
    print(f"Processed {result1} k-mers\n")

    print("Profiling UMI deduplication (slow)...")
    result2 = profile_umi_deduplication()
    print(f"Deduplicated {result2} UMI sets\n")

    print("Profiling UMI deduplication (fast)...")
    result3 = profile_umi_deduplication_fast()
    print(f"Deduplicated {result3} UMIs\n")

    print("\nFor line-by-line profiling, run:")
    print("  mprof run benchmarks/profile_memory.py")
    print("  mprof plot")
