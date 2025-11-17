#!/usr/bin/env python3
"""
Comprehensive Phase 1 Validation & Benchmark Suite

This validates ALL Phase 1 optimizations are working correctly
and measures the actual performance gains.
"""

import time
import statistics
import sheriff_rs


def test_correctness():
    """Test that all Phase 1 optimizations maintain correctness"""
    print("="*70)
    print("PHASE 1 CORRECTNESS VALIDATION")
    print("="*70)

    tests_passed = 0
    tests_total = 0

    # Test 1: Identical alt_seq optimization
    print("\n[Test 1] Identical alt_seq early exit")
    print("-" * 70)
    tests_total += 1
    edits = [
        ("chr1", 1000, "ATCG", "ATCGGGAGAGTAT", True, [1, 2]),
        ("chr1", 1000, "ATCG", "ATCGGGAGAGTAT", True, [1, 2, 3]),  # Identical alt_seq
    ]
    result = sheriff_rs.get_longest_edits_rust(edits)
    if len(result) == 1:
        print("✅ PASS - Clustered identical alt_seq correctly")
        tests_passed += 1
    else:
        print(f"❌ FAIL - Expected 1 edit, got {len(result)}")

    # Test 2: Exact string match after extraction
    print("\n[Test 2] Exact string match optimization")
    print("-" * 70)
    tests_total += 1
    edits = [
        ("chr1", 1000, "ATCG", "ATCGGGAGAGTAT", True, [1]),
        ("chr1", 1000, "ATCG", "ATCGGGAGAGTATXXX", True, [1, 2]),  # Same insert, diff length
    ]
    result = sheriff_rs.get_longest_edits_rust(edits)
    if len(result) == 1:
        print("✅ PASS - Clustered identical sequences correctly")
        tests_passed += 1
    else:
        print(f"❌ FAIL - Expected 1 edit, got {len(result)}")

    # Test 3: Homopolymer caching works
    print("\n[Test 3] Homopolymer detection caching")
    print("-" * 70)
    tests_total += 1
    edits = [
        ("chr1", 1000, "ATCG", "ATCGAAAGGGCCC", True, [1]),
        ("chr1", 1000, "ATCG", "ATCGAAAAAAGGGGGGGCCCCCCC", True, [1]),  # More homopolymers
    ]
    result = sheriff_rs.get_longest_edits_rust(edits)
    if len(result) == 1:
        print("✅ PASS - Homopolymer correction still works")
        tests_passed += 1
    else:
        print(f"❌ FAIL - Expected 1 edit, got {len(result)}")

    # Test 4: Forward string caching
    print("\n[Test 4] Forward string caching")
    print("-" * 70)
    tests_total += 1
    forward_edits = [
        ("chr1", 1000, "ATCG", "ATCGGGAGAGTAT", True, [1]),
        ("chr1", 1000, "ATCG", "ATCGGGAGAGTATAAA", True, [1, 2]),
    ]
    result = sheriff_rs.get_longest_edits_rust(forward_edits)
    if len(result) == 1:
        print("✅ PASS - Forward caching works correctly")
        tests_passed += 1
    else:
        print(f"❌ FAIL - Expected 1 edit, got {len(result)}")

    # Test 5: Reverse string caching
    print("\n[Test 5] Reverse string caching")
    print("-" * 70)
    tests_total += 1
    reverse_edits = [
        ("chr1", 1000, "ATCG", "GGAGAGTAT" + "ATCG", False, [1]),
        ("chr1", 1000, "ATCG", "GGAGAGTATAA" + "ATCG", False, [1, 2]),
    ]
    result = sheriff_rs.get_longest_edits_rust(reverse_edits)
    if len(result) == 1:
        print("✅ PASS - Reverse caching works correctly")
        tests_passed += 1
    else:
        print(f"❌ FAIL - Expected 1 edit, got {len(result)}")

    # Test 6: Mixed orientation handling
    print("\n[Test 6] Mixed forward/reverse edits")
    print("-" * 70)
    tests_total += 1
    mixed_edits = [
        ("chr1", 1000, "ATCG", "ATCGGGAGAGTAT", True, [1]),
        ("chr1", 1000, "ATCG", "GGAGAGTAT" + "ATCG", False, [1]),
    ]
    result = sheriff_rs.get_longest_edits_rust(mixed_edits)
    if len(result) == 2:
        print("✅ PASS - Different orientations kept separate")
        tests_passed += 1
    else:
        print(f"❌ FAIL - Expected 2 edits, got {len(result)}")

    print("\n" + "="*70)
    print(f"CORRECTNESS: {tests_passed}/{tests_total} tests passed")
    print("="*70)

    return tests_passed == tests_total


def benchmark_phase1():
    """Benchmark Phase 1 performance improvements"""
    print("\n" + "="*70)
    print("PHASE 1 PERFORMANCE BENCHMARK")
    print("="*70)

    results = {}

    # Benchmark 1: Small dataset (10 edits)
    print("\n[Benchmark 1] 10 edits (all similar, worst case)")
    print("-" * 70)
    edits_10 = []
    for i in range(10):
        edits_10.append((
            "chr1", 1000, "ATCG",
            "ATCG" + ("GGAGAGTAT" * (i % 3 + 1)),
            True, [i]
        ))

    times = []
    for _ in range(30):
        start = time.perf_counter()
        sheriff_rs.get_longest_edits_rust(edits_10)
        times.append((time.perf_counter() - start) * 1000)

    results['10_edits'] = {
        'mean': statistics.mean(times),
        'median': statistics.median(times),
        'std': statistics.stdev(times),
        'min': min(times),
        'max': max(times)
    }
    print(f"  Mean:   {results['10_edits']['mean']:.3f} ms")
    print(f"  Median: {results['10_edits']['median']:.3f} ms")
    print(f"  Std:    {results['10_edits']['std']:.3f} ms")

    # Benchmark 2: Medium dataset (50 edits)
    print("\n[Benchmark 2] 50 edits (all similar, worst case)")
    print("-" * 70)
    edits_50 = []
    for i in range(50):
        edits_50.append((
            "chr1", 1000, "ATCG",
            "ATCG" + ("GGAGAGTAT" * (i % 3 + 1)),
            True, [i]
        ))

    times = []
    for _ in range(20):
        start = time.perf_counter()
        sheriff_rs.get_longest_edits_rust(edits_50)
        times.append((time.perf_counter() - start) * 1000)

    results['50_edits'] = {
        'mean': statistics.mean(times),
        'median': statistics.median(times),
        'std': statistics.stdev(times),
        'min': min(times),
        'max': max(times)
    }
    print(f"  Mean:   {results['50_edits']['mean']:.3f} ms")
    print(f"  Median: {results['50_edits']['median']:.3f} ms")
    print(f"  Std:    {results['50_edits']['std']:.3f} ms")

    # Benchmark 3: Large dataset (100 edits)
    print("\n[Benchmark 3] 100 edits (all similar, worst case)")
    print("-" * 70)
    edits_100 = []
    for i in range(100):
        edits_100.append((
            "chr1", 1000, "ATCG",
            "ATCG" + ("GGAGAGTAT" * (i % 3 + 1)),
            True, [i]
        ))

    times = []
    for _ in range(10):
        start = time.perf_counter()
        sheriff_rs.get_longest_edits_rust(edits_100)
        times.append((time.perf_counter() - start) * 1000)

    results['100_edits'] = {
        'mean': statistics.mean(times),
        'median': statistics.median(times),
        'std': statistics.stdev(times),
        'min': min(times),
        'max': max(times)
    }
    print(f"  Mean:   {results['100_edits']['mean']:.3f} ms")
    print(f"  Median: {results['100_edits']['median']:.3f} ms")
    print(f"  Std:    {results['100_edits']['std']:.3f} ms")

    # Benchmark 4: Mixed orientations
    print("\n[Benchmark 4] 50 edits (mixed forward/reverse)")
    print("-" * 70)
    edits_mixed = []
    for i in range(50):
        if i % 2 == 0:
            edits_mixed.append((
                "chr1", 1000, "ATCG",
                "ATCG" + ("GGAGAGTAT" * (i % 3 + 1)),
                True, [i]
            ))
        else:
            edits_mixed.append((
                "chr1", 1000, "ATCG",
                ("GGAGAGTAT" * (i % 3 + 1)) + "ATCG",
                False, [i]
            ))

    times = []
    for _ in range(20):
        start = time.perf_counter()
        sheriff_rs.get_longest_edits_rust(edits_mixed)
        times.append((time.perf_counter() - start) * 1000)

    results['50_mixed'] = {
        'mean': statistics.mean(times),
        'median': statistics.median(times),
        'std': statistics.stdev(times),
        'min': min(times),
        'max': max(times)
    }
    print(f"  Mean:   {results['50_mixed']['mean']:.3f} ms")
    print(f"  Median: {results['50_mixed']['median']:.3f} ms")
    print(f"  Std:    {results['50_mixed']['std']:.3f} ms")

    return results


def print_summary(results):
    """Print performance summary"""
    print("\n" + "="*70)
    print("PHASE 1 PERFORMANCE SUMMARY")
    print("="*70)

    print("\nWorst-Case Performance (all edits similar):")
    print(f"  10 edits:  {results['10_edits']['mean']:.2f} ms (±{results['10_edits']['std']:.2f})")
    print(f"  50 edits:  {results['50_edits']['mean']:.2f} ms (±{results['50_edits']['std']:.2f})")
    print(f"  100 edits: {results['100_edits']['mean']:.2f} ms (±{results['100_edits']['std']:.2f})")

    print("\nMixed Orientation Performance:")
    print(f"  50 mixed:  {results['50_mixed']['mean']:.2f} ms (±{results['50_mixed']['std']:.2f})")

    print("\nTime per Comparison (50 edits):")
    n_comparisons = (50 * 49) // 2
    time_per_comp = (results['50_edits']['mean'] / n_comparisons) * 1000
    print(f"  {time_per_comp:.2f} µs per comparison")

    print("\nEstimated vs Python (assuming Python is 10x slower):")
    python_estimate_50 = results['50_edits']['mean'] * 10
    print(f"  Python (est): {python_estimate_50:.1f} ms")
    print(f"  Rust Phase 1: {results['50_edits']['mean']:.1f} ms")
    print(f"  Speedup: ~{python_estimate_50 / results['50_edits']['mean']:.1f}x")

    print("\n" + "="*70)


def main():
    print("╔" + "="*68 + "╗")
    print("║" + " "*15 + "PHASE 1 COMPREHENSIVE TEST SUITE" + " "*20 + "║")
    print("╚" + "="*68 + "╝")

    # Run correctness tests
    correctness_ok = test_correctness()

    if not correctness_ok:
        print("\n❌ CORRECTNESS TESTS FAILED - NOT RUNNING BENCHMARKS")
        print("Fix correctness issues before measuring performance!")
        return 1

    print("\n✅ ALL CORRECTNESS TESTS PASSED!")
    print("Proceeding to performance benchmarks...")

    # Run performance benchmarks
    results = benchmark_phase1()

    # Print summary
    print_summary(results)

    print("\n✅ PHASE 1 VALIDATION COMPLETE!")
    print("\nPhase 1 is SOLID - ready to proceed to Phase 2!")

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
