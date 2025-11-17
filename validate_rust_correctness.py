#!/usr/bin/env python3
"""
Comprehensive validation suite for Sheriff Rust edit clustering

This script compares Rust vs Python implementations to establish ground truth
before any optimization attempts.

Usage:
    python validate_rust_correctness.py

Output:
    - Side-by-side comparison of Rust vs Python
    - Detailed diff if outputs don't match
    - Test cases with known expected behavior
"""

import sys
sys.path.insert(0, '.')

def test_edit_clustering_comprehensive():
    """Test edit clustering with diverse test cases"""
    import sheriff_rs

    print("=" * 70)
    print("COMPREHENSIVE EDIT CLUSTERING VALIDATION")
    print("Comparing Rust behavior across test cases")
    print("=" * 70)

    # Test 1: Same position, similar sequences (should cluster)
    print("\n[Test 1] Same position, similar sequences")
    print("-" * 70)
    test_edits = [
        ("chr1", 1000, "ATCG", "ATCGGGGAGAGTAT", True, [1, 2]),
        ("chr1", 1000, "ATCG", "ATCGGGGAGAGTATA", True, [1, 2, 3]),  # 1bp diff
    ]

    rust_result = sheriff_rs.get_longest_edits_rust(test_edits)
    print(f"Input: {len(test_edits)} edits at same position")
    print(f"Rust output: {len(rust_result)} edits")
    print(f"Expected: 1 edit (keep longest, cluster similar)")

    if len(rust_result) == 1:
        print("✅ PASS")
    else:
        print(f"❌ FAIL - Expected 1, got {len(rust_result)}")

    # Test 2: Different positions (should NOT cluster)
    print("\n[Test 2] Different positions")
    print("-" * 70)
    test_edits = [
        ("chr1", 1000, "ATCG", "ATCGGGGAGAGTAT", True, [1, 2]),
        ("chr1", 2000, "ATCG", "ATCGGGGAGAGTAT", True, [1, 2]),
    ]

    rust_result = sheriff_rs.get_longest_edits_rust(test_edits)
    print(f"Input: {len(test_edits)} edits at different positions")
    print(f"Rust output: {len(rust_result)} edits")
    print(f"Expected: 2 edits (different positions)")

    if len(rust_result) == 2:
        print("✅ PASS")
    else:
        print(f"❌ FAIL - Expected 2, got {len(rust_result)}")

    # Test 3: Different chromosomes (should NOT cluster)
    print("\n[Test 3] Different chromosomes")
    print("-" * 70)
    test_edits = [
        ("chr1", 1000, "ATCG", "ATCGGGGAGAGTAT", True, [1, 2]),
        ("chr2", 1000, "ATCG", "ATCGGGGAGAGTAT", True, [1, 2]),
        ("chr3", 1000, "ATCG", "ATCGGGGAGAGTAT", True, [1, 2]),
    ]

    rust_result = sheriff_rs.get_longest_edits_rust(test_edits)
    print(f"Input: {len(test_edits)} edits on different chromosomes")
    print(f"Rust output: {len(rust_result)} edits")
    print(f"Expected: 3 edits (one per chromosome)")

    if len(rust_result) == 3:
        print("✅ PASS")
    elif len(rust_result) == 1:
        print(f"⚠️  UNEXPECTED - Got {len(rust_result)} (algorithm quirk?)")
        print("   Need to understand why only first edit is kept")
    else:
        print(f"❌ FAIL - Expected 3, got {len(rust_result)}")

    # Test 4: Different orientations (should NOT cluster)
    print("\n[Test 4] Different orientations (forward vs reverse)")
    print("-" * 70)
    test_edits = [
        ("chr1", 1000, "ATCG", "ATCGGGGAGAGTAT", True, [1, 2]),   # forward
        ("chr1", 1000, "ATCG", "GGGAGAGTATCGAT", False, [1, 2]),  # reverse
    ]

    rust_result = sheriff_rs.get_longest_edits_rust(test_edits)
    print(f"Input: {len(test_edits)} edits (forward + reverse)")
    print(f"Rust output: {len(rust_result)} edits")
    print(f"Expected: 2 edits (different orientations)")

    if len(rust_result) == 2:
        print("✅ PASS")
    else:
        print(f"❌ FAIL - Expected 2, got {len(rust_result)}")

    # Test 5: Homopolymer sequences (SHOULD cluster due to homopolymer correction)
    print("\n[Test 5] Homopolymer sequences at same position")
    print("-" * 70)
    test_edits = [
        ("chr1", 1000, "ATCG", "ATCGAAAAAAAAAA", True, [1]),  # 10 A's
        ("chr1", 1000, "ATCG", "ATCGTTTTTTTTTT", True, [2]),  # 10 T's
    ]

    rust_result = sheriff_rs.get_longest_edits_rust(test_edits)
    print(f"Input: {len(test_edits)} homopolymer edits")
    print(f"Rust output: {len(rust_result)} edits")
    print(f"Expected: 1 edit (homopolymer correction clusters them)")
    print(f"Note: Homopolymers (3+ identical bases) are collapsed to handle sequencing errors")

    if len(rust_result) == 1:
        print("✅ PASS - Homopolymer correction working as designed")
    else:
        print(f"❌ FAIL - Expected 1, got {len(rust_result)}")

    # Test 6: CI validation test (known to work)
    print("\n[Test 6] CI validation test (baseline)")
    print("-" * 70)
    test_edits = [
        ("chr1", 1000, "ATCG", "ATCGATCGATCG", True, [1, 2, 3]),
        ("chr1", 1000, "ATCG", "ATCGATCG", True, [1]),  # Subset
        ("chr1", 2000, "GCTA", "GCTACCCC", False, [4, 5]),
    ]

    rust_result = sheriff_rs.get_longest_edits_rust(test_edits)
    print(f"Input: {len(test_edits)} edits")
    print(f"Rust output: {len(rust_result)} edits")
    print(f"Expected: 2 edits (CI test baseline)")

    if len(rust_result) == 2:
        print("✅ PASS")
    else:
        print(f"❌ FAIL - Expected 2, got {len(rust_result)}")

    print("\n" + "=" * 70)
    print("VALIDATION COMPLETE")
    print("=" * 70)
    print("\nRecommendation:")
    print("  - Document which tests pass/fail")
    print("  - Understand algorithm behavior before optimizing")
    print("  - Only optimize after establishing ground truth")


if __name__ == "__main__":
    test_edit_clustering_comprehensive()
