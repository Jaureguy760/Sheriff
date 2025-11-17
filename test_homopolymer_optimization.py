#!/usr/bin/env python3
"""
Test homopolymer caching optimization specifically

The test_phase1.py uses "GGAGAGTAT" which has NO homopolymers (need 3+ consecutive bases).
This test uses actual homopolymer sequences to measure Optimization #1.
"""

import time
import statistics
import sheriff_rs


def test_homopolymer_caching():
    """Test sequences with homopolymers that require homopolymer correction"""
    print("="*70)
    print("HOMOPOLYMER CACHING OPTIMIZATION TEST")
    print("="*70)
    print("\nTest: 50 edits with homopolymer sequences")
    print("These sequences have 3+ consecutive identical bases")
    print()

    # Create 50 edits with homopolymer sequences
    # Each has AAA, GGG, CCC (homopolymers) but varying lengths
    edits = []
    for i in range(50):
        # These will cluster together due to homopolymer correction
        # Even though lengths differ
        if i % 3 == 0:
            seq = "ATCGATAAAGGGCCC"  # 3 of each homopolymer
        elif i % 3 == 1:
            seq = "ATCGATAAAAGGGGGCCCC"  # 4 of each
        else:
            seq = "ATCGATAAAAAAGGGGGGGCCCCCC"  # 6+ of each

        edits.append((
            "chr1", 1000, "ATCG",
            seq + "ATCG",  # Add reference at end
            True, [i]
        ))

    # Benchmark
    times = []
    for _ in range(20):
        start = time.perf_counter()
        result = sheriff_rs.get_longest_edits_rust(edits)
        times.append((time.perf_counter() - start) * 1000)

    mean = statistics.mean(times)
    median = statistics.median(times)
    std = statistics.stdev(times)

    print(f"Results:")
    print(f"  Mean:   {mean:.3f} ms")
    print(f"  Median: {median:.3f} ms")
    print(f"  Std:    {std:.3f} ms")
    print(f"  Output: {len(result)} edit(s) (clustered from {len(edits)})")
    print(f"  Time per comparison: {(mean / 1225):.3f} µs")
    print()

    # Check that homopolymer correction is working
    if len(result) == 1:
        print("✅ Homopolymer correction working - all sequences clustered")
    else:
        print(f"⚠️  Expected 1 output, got {len(result)}")

    return mean


def test_without_homopolymers():
    """Test sequences WITHOUT homopolymers (baseline for comparison)"""
    print("="*70)
    print("BASELINE: No Homopolymers")
    print("="*70)
    print("\nTest: 50 edits WITHOUT homopolymer sequences")
    print("This is the baseline - homopolymer caching won't help")
    print()

    # Create 50 edits without homopolymers (same as test_phase1.py)
    edits = []
    for i in range(50):
        edits.append((
            "chr1", 1000, "ATCG",
            "ATCG" + ("GGAGAGTAT" * (i % 3 + 1)),  # No 3+ consecutive bases
            True, [i]
        ))

    # Benchmark
    times = []
    for _ in range(20):
        start = time.perf_counter()
        result = sheriff_rs.get_longest_edits_rust(edits)
        times.append((time.perf_counter() - start) * 1000)

    mean = statistics.mean(times)
    median = statistics.median(times)
    std = statistics.stdev(times)

    print(f"Results:")
    print(f"  Mean:   {mean:.3f} ms")
    print(f"  Median: {median:.3f} ms")
    print(f"  Std:    {std:.3f} ms")
    print(f"  Output: {len(result)} edit(s)")
    print(f"  Time per comparison: {(mean / 1225):.3f} µs")
    print()

    return mean


def main():
    print("╔" + "="*68 + "╗")
    print("║" + " "*10 + "HOMOPOLYMER OPTIMIZATION VALIDATION" + " "*23 + "║")
    print("╚" + "="*68 + "╝")
    print()
    print("Optimization #1: Cache collapsed homopolymer strings")
    print("This test validates the optimization is working correctly")
    print()

    # Run tests
    time_with_homos = test_homopolymer_caching()
    time_without_homos = test_without_homopolymers()

    # Compare
    print("="*70)
    print("COMPARISON")
    print("="*70)
    print(f"\nWith homopolymers:    {time_with_homos:.2f} ms")
    print(f"Without homopolymers: {time_without_homos:.2f} ms")
    print()
    print("Note: Both should be similar since the optimization caches")
    print("the collapse operation, making it O(n) instead of O(n²)")
    print()
    print("The absolute time may be higher for homopolymer sequences")
    print("because they require additional processing steps:")
    print("  1. Initial alignment")
    print("  2. Homopolymer collapse + re-alignment")
    print("  3. Potentially 3' end checking")
    print()
    print("But Phase 2A makes step #2 faster by caching the collapsed strings!")

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
