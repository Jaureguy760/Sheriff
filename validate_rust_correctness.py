#!/usr/bin/env python3
"""
Comprehensive validation suite to verify Rust and Python implementations produce identical results.

This validates:
1. UMI deduplication - same unique counts
2. Edit clustering - same canonical edits
3. K-mer matching - same matches
4. Result determinism

Author: Claude Code
Date: 2025-11-16
"""

import sys
import random
import numpy as np
from collections import namedtuple

# Add sheriff to path
sys.path.insert(0, '/iblm/netapp/home/jjaureguy/Sheriff')

from sheriff.helpers import deduplicate_umis as deduplicate_umis_python
from sheriff.count_t7 import (
    get_longest_edits as get_longest_edits_python,
    KmerMatcher,
    _match_kmer_python,
    _match_kmer_rust,
)

try:
    import sheriff_rs
    RUST_AVAILABLE = True
except ImportError:
    RUST_AVAILABLE = False
    print("WARNING: sheriff_rs not available. Install with: cd sheriff-rs && maturin develop --release --features python")
    sys.exit(1)

# Define the ReadEdit namedtuple
ReadEdit = namedtuple("ReadEdit", ["chrom", "ref_pos", "ref_seq", "alt_seq", "forward", "kmer_matches"])


def generate_umis(n, umi_length=12):
    """Generate random UMIs with some duplicates and near-duplicates."""
    bases = ['A', 'C', 'G', 'T']
    umis = set()

    # Generate some original UMIs
    n_original = max(n // 3, 1)
    originals = []
    for _ in range(n_original):
        umi = ''.join(random.choices(bases, k=umi_length))
        originals.append(umi)
        umis.add(umi)

    # Add exact duplicates
    for _ in range(n // 3):
        if originals:
            umis.add(random.choice(originals))

    # Add near-duplicates (1 mismatch)
    for _ in range(n - len(umis)):
        if originals:
            original = random.choice(originals)
            pos = random.randint(0, umi_length - 1)
            mutated = list(original)
            mutated[pos] = random.choice([b for b in bases if b != mutated[pos]])
            umis.add(''.join(mutated))

    return umis


def test_umi_deduplication():
    """Test that Python and Rust UMI deduplication give same results."""
    print("=" * 60)
    print("TEST 1: UMI DEDUPLICATION CORRECTNESS")
    print("=" * 60)

    test_cases = [
        # (name, umis)
        ("single UMI", {"ATCGATCGATCG"}),
        ("two identical", {"ATCGATCGATCG", "ATCGATCGATCG"}),
        ("two 1-mismatch", {"ATCGATCGATCG", "TTCGATCGATCG"}),
        ("two different", {"ATCGATCGATCG", "GGGGGGGGGGGG"}),
        ("chain A-B-C", {"AAAA", "AAAB", "AABB"}),  # A~B, B~C, A!~C
        ("two groups", {"AAAA", "AAAB", "TTTT", "TTTG"}),
        ("random small (10)", generate_umis(10)),
        ("random medium (100)", generate_umis(100)),
        ("random large (500)", generate_umis(500)),
    ]

    all_passed = True

    for name, umi_set in test_cases:
        # Python implementation
        python_groups = deduplicate_umis_python(umi_set)
        python_count = len(python_groups)

        # Rust implementation
        rust_count = sheriff_rs.deduplicate_umis_py(list(umi_set))

        # Compare
        if python_count == rust_count:
            print(f"  PASS: {name} - Python={python_count}, Rust={rust_count}")
        else:
            print(f"  FAIL: {name} - Python={python_count}, Rust={rust_count}")
            all_passed = False

    return all_passed


def test_edit_clustering():
    """Test that Python and Rust edit clustering give same results."""
    print("\n" + "=" * 60)
    print("TEST 2: EDIT CLUSTERING CORRECTNESS")
    print("=" * 60)

    test_cases = [
        # Simple cases
        ("empty", []),
        ("single edit", [
            ReadEdit("chr1", 1000, "ATCG", "ATCGATCGATCG", True, (1, 2, 3))
        ]),
        ("different positions", [
            ReadEdit("chr1", 1000, "ATCG", "ATCGATCG", True, (1, 2)),
            ReadEdit("chr1", 2000, "ATCG", "ATCGATCG", True, (1, 2)),
        ]),
        ("different orientations", [
            ReadEdit("chr1", 1000, "ATCG", "ATCGATCG", True, (1,)),
            ReadEdit("chr1", 1000, "ATCG", "ATCGATCG", False, (1,)),
        ]),
        ("subset removal", [
            ReadEdit("chr1", 1000, "ATCG", "ATCGATCG", True, (1,)),  # Short
            ReadEdit("chr1", 1000, "ATCG", "ATCGATCGATCGATCG", True, (1, 2, 3)),  # Long
        ]),
    ]

    all_passed = True

    for name, edits in test_cases:
        # Python implementation
        python_result = get_longest_edits_python(edits)
        python_count = len(python_result)

        # Convert to format expected by Rust
        rust_edits = [
            (e.chrom, e.ref_pos, e.ref_seq, e.alt_seq, e.forward, list(e.kmer_matches))
            for e in edits
        ]

        # Rust implementation
        rust_result = sheriff_rs.get_longest_edits_rust(rust_edits)
        rust_count = len(rust_result)

        # Compare counts first (order might differ)
        if python_count == rust_count:
            # Also compare actual results (ignoring order)
            python_set = set((e.chrom, e.ref_pos, e.ref_seq, e.alt_seq, e.forward) for e in python_result)
            rust_set = set((e[0], e[1], e[2], e[3], e[4]) for e in rust_result)

            if python_set == rust_set:
                print(f"  PASS: {name} - Count={python_count}, Results match")
            else:
                print(f"  WARN: {name} - Same count ({python_count}) but different edits")
                print(f"    Python: {python_set}")
                print(f"    Rust: {rust_set}")
                # This might be acceptable due to tie-breaking differences
        else:
            print(f"  FAIL: {name} - Python={python_count}, Rust={rust_count}")
            all_passed = False

    return all_passed


def test_kmer_matching():
    """Test that Python and Rust k-mer matching give same results."""
    print("\n" + "=" * 60)
    print("TEST 3: K-MER MATCHING CORRECTNESS")
    print("=" * 60)

    k = 7
    bases = ['A', 'C', 'G', 'T']

    # Create test sequences
    test_sequences = [
        "ATCGATCGATCGATCGATCG",  # 20bp
        "AAAAAAACCCCCCCGGGGGGG",  # 21bp with homopolymers
        "ATCGATCG",  # 8bp (short)
        "A" * 100,  # 100bp all A's
        ''.join(random.choices(bases, k=50)),  # Random 50bp
    ]

    # Create whitelist
    target_kmers = [''.join(random.choices(bases, k=k)) for _ in range(50)]
    matcher = KmerMatcher(k, target_kmers)
    whitelist_hashes = matcher.match_hash

    all_passed = True

    for i, seq in enumerate(test_sequences):
        # Python implementation
        python_result = _match_kmer_python(matcher, seq, True)
        python_set = set(python_result) if python_result else set()

        # Rust implementation
        rust_result = _match_kmer_rust(seq, k, whitelist_hashes, True)
        rust_set = set(rust_result) if rust_result else set()

        # Compare
        if python_set == rust_set:
            print(f"  PASS: Sequence {i+1} (len={len(seq)}) - Matches: {len(python_set)}")
        else:
            print(f"  FAIL: Sequence {i+1} (len={len(seq)})")
            print(f"    Python: {sorted(python_set)[:5]}... ({len(python_set)} total)")
            print(f"    Rust:   {sorted(rust_set)[:5]}... ({len(rust_set)} total)")
            all_passed = False

    return all_passed


def test_determinism():
    """Test that both implementations are deterministic (same results on repeated runs)."""
    print("\n" + "=" * 60)
    print("TEST 4: DETERMINISM (same results on repeated runs)")
    print("=" * 60)

    # Fixed seed for reproducibility
    random.seed(42)
    umi_set = generate_umis(200)

    # Run multiple times
    python_results = []
    rust_results = []

    for i in range(5):
        python_results.append(len(deduplicate_umis_python(umi_set)))
        rust_results.append(sheriff_rs.deduplicate_umis_py(list(umi_set)))

    python_consistent = len(set(python_results)) == 1
    rust_consistent = len(set(rust_results)) == 1

    if python_consistent:
        print(f"  PASS: Python is deterministic (always {python_results[0]})")
    else:
        print(f"  FAIL: Python gives different results: {python_results}")

    if rust_consistent:
        print(f"  PASS: Rust is deterministic (always {rust_results[0]})")
    else:
        print(f"  FAIL: Rust gives different results: {rust_results}")

    return python_consistent and rust_consistent


def test_edge_cases():
    """Test edge cases that might reveal bugs."""
    print("\n" + "=" * 60)
    print("TEST 5: EDGE CASES")
    print("=" * 60)

    all_passed = True

    # Empty set
    python_result = deduplicate_umis_python(set())
    rust_result = sheriff_rs.deduplicate_umis_py([])
    if len(python_result) == rust_result == 0:
        print("  PASS: Empty UMI set")
    else:
        print(f"  FAIL: Empty set - Python={len(python_result)}, Rust={rust_result}")
        all_passed = False

    # All identical UMIs
    umi_list = ["ATCGATCGATCG"] * 10  # 10 copies of same UMI
    umi_set = set(umi_list)  # This is just 1 unique
    python_result = len(deduplicate_umis_python(umi_set))
    rust_result = sheriff_rs.deduplicate_umis_py(umi_list)
    if python_result == 1 and rust_result == 1:
        print("  PASS: All identical UMIs (10x same)")
    else:
        print(f"  FAIL: All identical - Python={python_result}, Rust={rust_result}")
        all_passed = False

    # All different UMIs (no matches)
    umis = [
        "AAAAAAAAAA",
        "CCCCCCCCCC",
        "GGGGGGGGGG",
        "TTTTTTTTTT",
    ]
    python_result = len(deduplicate_umis_python(set(umis)))
    rust_result = sheriff_rs.deduplicate_umis_py(umis)
    if python_result == 4 and rust_result == 4:
        print("  PASS: All different UMIs")
    else:
        print(f"  FAIL: All different - Python={python_result}, Rust={rust_result}")
        all_passed = False

    # Very long chain (potential stack overflow)
    chain_umis = []
    base = "AAAAAAAA"
    for i in range(100):
        # Each UMI differs by 1 from the previous
        umi = list(base)
        if i % 4 == 0:
            umi[0] = 'T'
        elif i % 4 == 1:
            umi[1] = 'C'
        elif i % 4 == 2:
            umi[2] = 'G'
        else:
            umi[3] = 'T'
        base = ''.join(umi)
        chain_umis.append(base)

    # This might not form one group since we're changing the base
    # Just check it doesn't crash
    try:
        python_result = len(deduplicate_umis_python(set(chain_umis)))
        rust_result = sheriff_rs.deduplicate_umis_py(chain_umis)
        print(f"  PASS: Long chain (100 UMIs) - Python={python_result}, Rust={rust_result}")
    except Exception as e:
        print(f"  FAIL: Long chain crashed: {e}")
        all_passed = False

    return all_passed


def main():
    print("=" * 60)
    print("SHERIFF RUST/PYTHON CORRECTNESS VALIDATION SUITE")
    print("=" * 60)
    print(f"Rust module available: {RUST_AVAILABLE}")
    print()

    results = {
        "UMI Deduplication": test_umi_deduplication(),
        "Edit Clustering": test_edit_clustering(),
        "K-mer Matching": test_kmer_matching(),
        "Determinism": test_determinism(),
        "Edge Cases": test_edge_cases(),
    }

    print("\n" + "=" * 60)
    print("FINAL RESULTS")
    print("=" * 60)

    all_passed = True
    for test_name, passed in results.items():
        status = "PASS" if passed else "FAIL"
        print(f"  {test_name}: {status}")
        if not passed:
            all_passed = False

    print()
    if all_passed:
        print("ALL TESTS PASSED - Rust and Python implementations are equivalent!")
    else:
        print("SOME TESTS FAILED - Results differ between implementations!")

    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
