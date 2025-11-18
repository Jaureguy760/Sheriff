#!/usr/bin/env python3
"""
Quick test script for sheriff-rs Python bindings

This script runs basic sanity checks to verify the bindings work correctly.
Run this after building with: maturin develop --release --features python

Usage:
    python3 test_python_bindings.py
"""

import sys

def test_import():
    """Test that the module can be imported."""
    print("Test 1: Module import...", end=" ")
    try:
        import sheriff_rs
        print("✓ PASS")
        return sheriff_rs
    except ImportError as e:
        print(f"✗ FAIL: {e}")
        return None


def test_version(sheriff_rs):
    """Test version info."""
    print("Test 2: Version info...", end=" ")
    try:
        version = sheriff_rs.__version__
        if version:
            print(f"✓ PASS (v{version})")
            return True
        else:
            print("✗ FAIL: No version")
            return False
    except Exception as e:
        print(f"✗ FAIL: {e}")
        return False


def test_kmer_to_num(sheriff_rs):
    """Test k-mer to number conversion."""
    print("Test 3: kmer_to_num...", end=" ")
    try:
        result = sheriff_rs.kmer_to_num("ACGT")
        expected = 27
        if result == expected:
            print(f"✓ PASS ({result})")
            return True
        else:
            print(f"✗ FAIL: Expected {expected}, got {result}")
            return False
    except Exception as e:
        print(f"✗ FAIL: {e}")
        return False


def test_match_kmer(sheriff_rs):
    """Test k-mer matching."""
    print("Test 4: match_kmer...", end=" ")
    try:
        whitelist = [sheriff_rs.kmer_to_num("ACGT")]
        matches = sheriff_rs.match_kmer("ACGTACGT", 4, whitelist, output_hash=True)
        if len(matches) == 2 and all(m == 27 for m in matches):
            print(f"✓ PASS (found {len(matches)} matches)")
            return True
        else:
            print(f"✗ FAIL: Expected [27, 27], got {matches}")
            return False
    except Exception as e:
        print(f"✗ FAIL: {e}")
        return False


def test_kmer_counter(sheriff_rs):
    """Test KmerCounter class."""
    print("Test 5: KmerCounter...", end=" ")
    try:
        counter = sheriff_rs.KmerCounter(4)
        freqs = counter.count_kmers("ACGTACGT")
        acgt_hash = sheriff_rs.kmer_to_num("ACGT")
        count = freqs[acgt_hash]
        if count == 2:
            print(f"✓ PASS (count={count})")
            return True
        else:
            print(f"✗ FAIL: Expected count=2, got {count}")
            return False
    except Exception as e:
        print(f"✗ FAIL: {e}")
        return False


def test_deduplicate_umis(sheriff_rs):
    """Test UMI deduplication."""
    print("Test 6: deduplicate_umis...", end=" ")
    try:
        umis = ["ATCGATCG", "ATCGATCC", "GCGCGCGC"]
        unique_count = sheriff_rs.deduplicate_umis(umis, threshold=1)
        if unique_count == 2:
            print(f"✓ PASS (unique_count={unique_count})")
            return True
        else:
            print(f"✗ FAIL: Expected 2 unique groups, got {unique_count}")
            return False
    except Exception as e:
        print(f"✗ FAIL: {e}")
        return False


def test_deduplicate_umis_detailed(sheriff_rs):
    """Test detailed UMI deduplication."""
    print("Test 7: deduplicate_umis_detailed...", end=" ")
    try:
        umis = ["ATCGATCG", "ATCGATCC", "GCGCGCGC"]
        groups = sheriff_rs.deduplicate_umis_detailed(umis, threshold=1)
        if len(groups) == 2:
            print(f"✓ PASS ({len(groups)} groups)")
            return True
        else:
            print(f"✗ FAIL: Expected 2 groups, got {len(groups)}")
            return False
    except Exception as e:
        print(f"✗ FAIL: {e}")
        return False


def test_hamming_distance(sheriff_rs):
    """Test Hamming distance function."""
    print("Test 8: hamming_distance...", end=" ")
    try:
        dist = sheriff_rs.hamming_distance("ATCG", "ATGG")
        if dist == 1:
            print(f"✓ PASS (distance={dist})")
            return True
        else:
            print(f"✗ FAIL: Expected distance=1, got {dist}")
            return False
    except Exception as e:
        print(f"✗ FAIL: {e}")
        return False


def test_error_handling(sheriff_rs):
    """Test error handling."""
    print("Test 9: Error handling...", end=" ")
    try:
        # This should raise ValueError
        sheriff_rs.kmer_to_num("ACGN")
        print("✗ FAIL: Should have raised ValueError")
        return False
    except ValueError:
        print("✓ PASS")
        return True
    except Exception as e:
        print(f"✗ FAIL: Wrong exception: {e}")
        return False


def main():
    """Run all tests."""
    print("="*60)
    print("  Sheriff-rs Python Bindings Test Suite")
    print("="*60)
    print()

    # Import module
    sheriff_rs = test_import()
    if not sheriff_rs:
        print("\n✗ Module import failed. Please build with:")
        print("  cd sheriff-rs")
        print("  maturin develop --release --features python")
        sys.exit(1)

    print()

    # Run tests
    tests = [
        test_version,
        test_kmer_to_num,
        test_match_kmer,
        test_kmer_counter,
        test_deduplicate_umis,
        test_deduplicate_umis_detailed,
        test_hamming_distance,
        test_error_handling,
    ]

    results = [test(sheriff_rs) for test in tests]

    # Summary
    print()
    print("="*60)
    passed = sum(results)
    total = len(results)
    if passed == total:
        print(f"  ✓ All {total} tests PASSED!")
    else:
        print(f"  {passed}/{total} tests passed, {total-passed} failed")
    print("="*60)
    print()

    if passed == total:
        print("✓ Python bindings are working correctly!")
        print()
        print("Next steps:")
        print("  - Run the full demo: python3 examples/python_demo.py")
        print("  - Read the docs: help(sheriff_rs)")
        print()
        sys.exit(0)
    else:
        print("✗ Some tests failed. Please check the output above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
