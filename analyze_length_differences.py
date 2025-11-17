#!/usr/bin/env python3
"""
Analyze length differences in edits that cluster together

This script helps us understand what length differences we actually see
in practice for edits that SHOULD cluster together, so we can determine
a safe threshold for length-based filtering.
"""

import sheriff_rs

def analyze_length_filter_safety():
    """
    Test different length filter thresholds to see what we'd miss
    """

    print("="*70)
    print("LENGTH FILTER SAFETY ANALYSIS")
    print("="*70)

    # Create test cases that SHOULD cluster but differ in length
    test_cases = [
        {
            'name': 'Truncated read (substring)',
            'edits': [
                ("chr1", 1000, "ATCG", "ATCGATCGATCG", True, [1, 2]),        # 8bp insert
                ("chr1", 1000, "ATCG", "ATCGATCGATCGATCG", True, [1, 2, 3]), # 12bp insert
            ],
            'len_diff': 4,
            'should_cluster': True
        },
        {
            'name': 'Homopolymer difference',
            'edits': [
                ("chr1", 1000, "ATCG", "ATCGATAATACTCTCCCTATTCACTCTGCGT", True, [1]),     # CCC
                ("chr1", 1000, "ATCG", "ATCGATAATACTCTCCCCCCATTCACT", True, [1]),         # CCCCCC
            ],
            'len_diff': 3,
            'should_cluster': True
        },
        {
            'name': 'Very different lengths (15bp diff)',
            'edits': [
                ("chr1", 1000, "ATCG", "ATCGGGAGAGTAT", True, [1]),                    # 9bp insert
                ("chr1", 1000, "ATCG", "ATCGGGAGAGTATAGAATGGAGCTTTT", True, [1]),     # 24bp insert
            ],
            'len_diff': 15,
            'should_cluster': '???'  # Need to test!
        },
        {
            'name': 'Completely different sequences (same length)',
            'edits': [
                ("chr1", 1000, "ATCG", "ATCGAAAAAAAAAA", True, [1]),  # 10 A's
                ("chr1", 1000, "ATCG", "ATCGTTTTTTTTTT", True, [2]),  # 10 T's
            ],
            'len_diff': 0,
            'should_cluster': False  # But homopolymer correction might cluster them!
        },
        {
            'name': 'Very different AND different lengths (20bp diff)',
            'edits': [
                ("chr1", 1000, "ATCG", "ATCGAAAAAA", True, [1]),                        # 6bp insert
                ("chr1", 1000, "ATCG", "ATCGTTTTTTTTTTTTTTTTTTTTTT", True, [1]),      # 22bp insert
            ],
            'len_diff': 16,
            'should_cluster': False
        },
    ]

    print("\nTesting various length differences:\n")

    for test in test_cases:
        edits = test['edits']
        result = sheriff_rs.get_longest_edits_rust(edits)

        clustered = len(result) == 1
        expected = test['should_cluster']

        print(f"Test: {test['name']}")
        print(f"  Length difference: {test['len_diff']}bp")
        print(f"  Input edits: {len(edits)}")
        print(f"  Output edits: {len(result)}")
        print(f"  Clustered: {clustered}")
        print(f"  Expected to cluster: {expected}")

        if expected == '???':
            print(f"  Result: Unknown - actual behavior is {clustered}")
        elif clustered == expected:
            print(f"  ✅ PASS")
        else:
            print(f"  ❌ FAIL - Expected {expected}, got {clustered}")
        print()

    print("="*70)
    print("CONCLUSIONS:")
    print("="*70)
    print("""
Based on these tests, we can determine:

1. What length difference is SAFE to filter?
   - If we see clustering at length_diff=X, we cannot filter at threshold X
   - We need to be VERY conservative

2. Does homopolymer correction affect length-based filtering?
   - YES - homopolymers can create apparent length differences

3. Recommendation:
   - Either use VERY high threshold (>20bp) and verify with real data
   - Or DON'T use length-based filtering at all (too risky)

The safest approach: Only use optimizations that cannot break correctness!
    """)

if __name__ == "__main__":
    analyze_length_filter_safety()
