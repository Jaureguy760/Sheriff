#!/usr/bin/env python3
"""
CORRECTED Benchmark: Rust vs Python UMI Deduplication
Fixed to compare apples-to-apples
"""

import time
import sys
sys.path.insert(0, '/home/user/Sheriff')

from sheriff.helpers import deduplicate_umis
import sheriff_rs

def benchmark_fair(n_umis):
    """Fair benchmark with truly unique UMIs"""
    # Generate unique UMIs
    umis = [f'{"ACGT"[i%4]}{"ACGT"[(i//4)%4]}{"ACGT"[(i//16)%4]}{"ACGT"[(i//64)%4]}{"ACGT"[(i//256)%4]}{"ACGT"[(i//1024)%4]}{"ACGT"[(i//4096)%4]}{"ACGT"[(i//16384)%4]}'[:8]
            for i in range(n_umis)]

    # Ensure we have n_umis unique strings
    umis_set = set(umis)
    print(f'  Generated {len(umis_set)} unique UMIs')

    # Python
    start = time.perf_counter()
    py_result = deduplicate_umis(umis_set)
    py_time = (time.perf_counter() - start) * 1000

    # Rust
    start = time.perf_counter()
    rust_result = sheriff_rs.deduplicate_umis(list(umis_set), threshold=1)
    rust_time = (time.perf_counter() - start) * 1000

    speedup = py_time / rust_time if rust_time > 0 else 0

    return {
        'n': len(umis_set),
        'py_time': py_time,
        'rust_time': rust_time,
        'py_result': len(py_result),
        'rust_result': rust_result,
        'speedup': speedup,
    }

print("=" * 70)
print("CORRECTED Benchmark: UMI Deduplication (Apples-to-Apples)")
print("=" * 70)
print()

for n in [10, 50, 100, 200, 500]:
    print(f"Testing {n} unique UMIs:")
    result = benchmark_fair(n)

    winner = "Rust 🚀" if result['speedup'] > 1 else "Python 🐍"
    print(f"  Python: {result['py_time']:8.3f} ms -> {result['py_result']} groups")
    print(f"  Rust:   {result['rust_time']:8.3f} ms -> {result['rust_result']} groups")
    print(f"  Winner: {winner} ({result['speedup']:.2f}x)")
    print()

print("=" * 70)
print("Conclusion: Both implementations have O(n²) complexity.")
print("Python uses pure Python sets (slow), Rust uses Union-Find (fast).")
print("For Sheriff's real usage (small batches per cell), both are adequate.")
print("=" * 70)
