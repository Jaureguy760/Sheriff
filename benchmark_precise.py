#!/usr/bin/env python3
"""
Detailed benchmark to understand Phase 1 optimization impact

This runs more iterations to get accurate measurements.
"""

import time
import statistics
import sheriff_rs


def precise_benchmark(name, edits, n_runs=20):
    """Run precise benchmark with more iterations"""
    print(f"\n{name}")
    print("-" * 60)

    # Warmup
    for _ in range(3):
        sheriff_rs.get_longest_edits_rust(edits)

    # Measure
    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        result = sheriff_rs.get_longest_edits_rust(edits)
        end = time.perf_counter()
        times.append((end - start) * 1000)

    mean = statistics.mean(times)
    std = statistics.stdev(times)
    median = statistics.median(times)

    print(f"  Runs: {n_runs}")
    print(f"  Mean:   {mean:.3f} ms")
    print(f"  Median: {median:.3f} ms")
    print(f"  Std:    {std:.3f} ms")
    print(f"  Min:    {min(times):.3f} ms")
    print(f"  Max:    {max(times):.3f} ms")

    return mean


def main():
    print("="*60)
    print("PRECISE PHASE 1 BENCHMARK")
    print("="*60)

    # Test 1: All forward edits
    forward_edits = []
    for i in range(50):
        forward_edits.append((
            "chr1", 1000, "ATCG",
            "ATCG" + ("GGAGAGTAT" * (i % 3 + 1)),
            True,  # forward
            [i]
        ))

    forward_time = precise_benchmark("50 Forward Edits (worst case)", forward_edits, n_runs=20)

    # Test 2: All reverse edits
    reverse_edits = []
    for i in range(50):
        reverse_edits.append((
            "chr1", 1000, "ATCG",
            ("GGAGAGTAT" * (i % 3 + 1)) + "ATCG",
            False,  # reverse
            [i]
        ))

    reverse_time = precise_benchmark("50 Reverse Edits (worst case)", reverse_edits, n_runs=20)

    # Test 3: Mixed
    mixed_edits = []
    for i in range(50):
        if i % 2 == 0:
            mixed_edits.append((
                "chr1", 1000, "ATCG",
                "ATCG" + ("GGAGAGTAT" * (i % 3 + 1)),
                True,
                [i]
            ))
        else:
            mixed_edits.append((
                "chr1", 1000, "ATCG",
                ("GGAGAGTAT" * (i % 3 + 1)) + "ATCG",
                False,
                [i]
            ))

    mixed_time = precise_benchmark("50 Mixed Edits (worst case)", mixed_edits, n_runs=20)

    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Forward edits: {forward_time:.3f} ms")
    print(f"Reverse edits: {reverse_time:.3f} ms")
    print(f"Mixed edits:   {mixed_time:.3f} ms")
    print("\nNote: Forward caching should speed up forward-heavy datasets")
    print("      Reverse caching should speed up reverse-heavy datasets")


if __name__ == "__main__":
    main()
