#!/usr/bin/env python3
"""
Debug why certain sequences cluster together

This helps us understand if the algorithm is working correctly or if
there's a bug causing over-clustering.
"""

import sheriff_rs

def debug_sequence_pair(name, edit1, edit2):
    """Debug a specific pair of edits"""
    print(f"\n{'='*70}")
    print(f"Test: {name}")
    print(f"{'='*70}")

    chr1, pos1, ref1, alt1, fwd1, kmers1 = edit1
    chr2, pos2, ref2, alt2, fwd2, kmers2 = edit2

    # Extract the insert sequences (remove reference portion)
    reflen = len(ref1)
    if fwd1:
        seq1 = alt1[:-reflen] if reflen > 0 else alt1
    else:
        seq1 = alt1[reflen:]

    if fwd2:
        seq2 = alt2[:-reflen] if reflen > 0 else alt2
    else:
        seq2 = alt2[reflen:]

    print(f"\nEdit 1:")
    print(f"  Full: {alt1}")
    print(f"  Insert: {seq1} ({len(seq1)}bp)")
    print(f"  Chr: {chr1}, Pos: {pos1}, Forward: {fwd1}")

    print(f"\nEdit 2:")
    print(f"  Full: {alt2}")
    print(f"  Insert: {seq2} ({len(seq2)}bp)")
    print(f"  Chr: {chr2}, Pos: {pos2}, Forward: {fwd2}")

    print(f"\nProperties:")
    print(f"  Same chromosome: {chr1 == chr2}")
    print(f"  Same position: {pos1 == pos2}")
    print(f"  Same orientation: {fwd1 == fwd2}")
    print(f"  Length difference: {abs(len(seq1) - len(seq2))}bp")

    # Test clustering
    result = sheriff_rs.get_longest_edits_rust([edit1, edit2])
    clustered = len(result) == 1

    print(f"\nResult:")
    print(f"  Input: 2 edits")
    print(f"  Output: {len(result)} edit(s)")
    print(f"  Clustered: {clustered}")

    if clustered:
        print(f"  Kept edit: {result[0][3]} ({len(result[0][3])}bp total)")

    return clustered


def main():
    print("="*70)
    print("CLUSTERING LOGIC DEBUG")
    print("="*70)

    # Test 1: Homopolymers (should cluster due to homopolymer correction)
    debug_sequence_pair(
        "Homopolymers - all A's vs all T's",
        ("chr1", 1000, "ATCG", "ATCGAAAAAAAAAA", True, [1]),  # 10 A's
        ("chr1", 1000, "ATCG", "ATCGTTTTTTTTTT", True, [2]),  # 10 T's
    )

    # Test 2: Very different sequences, same length
    debug_sequence_pair(
        "Completely different sequences (no homopolymers)",
        ("chr1", 1000, "ATCG", "ATCGAGCTACGTAG", True, [1]),  # Mixed sequence
        ("chr1", 1000, "ATCG", "ATCGTCGATCGACG", True, [2]),  # Different mixed
    )

    # Test 3: Different length, different sequence
    debug_sequence_pair(
        "Different sequences AND different lengths",
        ("chr1", 1000, "ATCG", "ATCGAGCTACGTAG", True, [1]),              # 10bp insert
        ("chr1", 1000, "ATCG", "ATCGTCGATCGACGTCGATCGA", True, [1]),     # 18bp insert
    )

    # Test 4: Substring relationship
    debug_sequence_pair(
        "One is substring of other (should cluster)",
        ("chr1", 1000, "ATCG", "ATCGGGAGAGTAT", True, [1]),              # 9bp
        ("chr1", 1000, "ATCG", "ATCGGGAGAGTATAGAA", True, [1]),          # 13bp (extends the 9bp)
    )

    # Test 5: Very long homopolymer difference
    debug_sequence_pair(
        "Extreme homopolymer difference",
        ("chr1", 1000, "ATCG", "ATCGAAAGGGCCC", True, [1]),              # AAA, GGG, CCC
        ("chr1", 1000, "ATCG", "ATCGAAAAAAGGGGGGGCCCCCCC", True, [1]),  # Many A's, G's, C's
    )

    print("\n" + "="*70)
    print("KEY FINDINGS:")
    print("="*70)
    print("""
If homopolymers cause very different sequences to cluster:
  → Homopolymer correction is working as designed
  → Length-based filtering would break this correction
  → We CANNOT use simple length filtering

If non-homopolymer sequences cluster when they shouldn't:
  → There may be a bug in the clustering logic
  → Need to investigate the edit distance calculation

Recommendation: Focus on optimizations that don't depend on length!
    """)


if __name__ == "__main__":
    main()
