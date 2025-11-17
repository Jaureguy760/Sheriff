#!/usr/bin/env python3
"""
CI-friendly validation for Sheriff Rust acceleration.

This validates Rust functions work correctly without needing the full pipeline
or 3GB reference genome. Suitable for GitHub Actions.
"""

import sys
import json
import hashlib
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.absolute()
SHERIFF_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SHERIFF_ROOT))

CHECKSUM_FILE = SCRIPT_DIR / "expected_checksums.json"


def md5_text(text):
    """Calculate MD5 of text."""
    return hashlib.md5(text.encode('utf-8')).hexdigest()


def md5_file(filepath):
    """Calculate MD5 of binary file."""
    hash_md5 = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def md5_text_file(filepath):
    """Calculate MD5 of text file with normalized line endings."""
    hash_md5 = hashlib.md5()
    with open(filepath, "r") as f:
        content = f.read()
        content = content.replace('\r\n', '\n').rstrip('\n')
        hash_md5.update(content.encode('utf-8'))
    return hash_md5.hexdigest()


def verify_input_checksums():
    """Verify test input files match expected checksums."""
    print("=== Verifying Test Input Files ===")

    if not CHECKSUM_FILE.exists():
        print(f"WARNING: {CHECKSUM_FILE} not found, skipping checksum verification")
        return True

    with open(CHECKSUM_FILE) as f:
        expected = json.load(f)

    all_pass = True
    for fname, expected_md5 in expected.get("test_files", {}).items():
        fpath = SCRIPT_DIR / fname
        if not fpath.exists():
            print(f"  {fname}: MISSING")
            all_pass = False
            continue

        # Use binary MD5 for binary files, text MD5 for text files
        if fname.endswith(('.bam', '.bai')):
            actual_md5 = md5_file(fpath)
        else:
            actual_md5 = md5_text_file(fpath)

        if actual_md5 == expected_md5:
            print(f"  {fname}: OK ({expected_md5[:8]}...)")
        else:
            print(f"  {fname}: MISMATCH!")
            print(f"    Expected: {expected_md5}")
            print(f"    Got:      {actual_md5}")
            all_pass = False

    return all_pass


def test_rust_imports():
    """Test that Rust module imports correctly."""
    print("\n=== Testing Rust Module Imports ===")

    try:
        import sheriff_rs
        print(f"  sheriff_rs: OK")
        return True
    except ImportError as e:
        print(f"  sheriff_rs: FAIL - {e}")
        return False


def test_kmer_counting():
    """Test k-mer counting produces deterministic results."""
    print("\n=== Testing K-mer Counting ===")

    import sheriff_rs

    # Known test case
    seq = "ATCGATCGATCGATCGATCGATCG"
    k = 4
    result = sheriff_rs.count_kmers_rust(seq, k)

    # Should return 4^4 = 256 bins
    assert len(result) == 256, f"Expected 256 bins, got {len(result)}"

    # Calculate checksum of result for determinism
    result_str = ",".join(map(str, result))
    result_md5 = md5_text(result_str)

    # Known expected checksum (deterministic)
    expected_md5 = "03fd3b75164217eb64b8cf8a1bc2eb20"

    if result_md5 == expected_md5:
        print(f"  4-mer counting: OK (checksum {result_md5[:8]}...)")
        return True
    else:
        print(f"  4-mer counting: CHECKSUM MISMATCH")
        print(f"    Expected: {expected_md5}")
        print(f"    Got:      {result_md5}")
        # Show actual counts for debugging
        total = sum(result)
        print(f"    Total k-mers: {total}, bins: {len(result)}")
        return False


def test_umi_deduplication():
    """Test UMI deduplication produces correct results."""
    print("\n=== Testing UMI Deduplication ===")

    import sheriff_rs

    # Test case: 4 UMIs forming 2 groups
    # Group 1: ATCGATCG, ATCGATCC (differ by 1)
    # Group 2: TTTTTTTT, TTTTTTCT (differ by 1)
    umis = ["ATCGATCG", "ATCGATCC", "TTTTTTTT", "TTTTTTCT"]
    unique = sheriff_rs.deduplicate_umis_py(umis)

    if unique == 2:
        print(f"  UMI dedup (4 UMIs → 2 groups): OK")
        return True
    else:
        print(f"  UMI dedup: FAIL (expected 2, got {unique})")
        return False


def test_edit_clustering():
    """Test edit clustering produces correct results."""
    print("\n=== Testing Edit Clustering ===")

    import sheriff_rs

    # Test case: 3 edits where second is subset of first
    edits = [
        ("chr1", 1000, "ATCG", "ATCGATCGATCG", True, [1, 2, 3]),  # Longest
        ("chr1", 1000, "ATCG", "ATCGATCG", True, [1]),  # Subset - should be removed
        ("chr1", 2000, "GCTA", "GCTACCCC", False, [4, 5]),  # Different position
    ]

    longest = sheriff_rs.get_longest_edits_rust(edits)

    # Should return 2 edits (removing the subset)
    if len(longest) == 2:
        # Verify correct edits returned
        positions = {e[1] for e in longest}
        if positions == {1000, 2000}:
            # Verify the longest alt_seq was kept for position 1000
            for edit in longest:
                if edit[1] == 1000:
                    if edit[3] == "ATCGATCGATCG":
                        print(f"  Edit clustering (3→2 edits): OK")
                        return True
                    else:
                        print(f"  Edit clustering: WRONG alt_seq kept")
                        return False
        else:
            print(f"  Edit clustering: WRONG positions {positions}")
            return False
    else:
        print(f"  Edit clustering: FAIL (expected 2, got {len(longest)})")
        return False


def test_cell_umi_counting():
    """Test cell-level UMI counting."""
    print("\n=== Testing Cell UMI Counting ===")

    import sheriff_rs

    # Test case: 2 cells, cell 0 has 2 reads, cell 1 has 1 read
    cell_bc_indexes = [0, 1, 0]
    cell_umis = [["ATCGATCG"], ["TTTTTTTT", "CCCCCCCC"], ["GGGGGGGG"]]
    total_cells = 2

    counts = sheriff_rs.cell_umi_counts_py(cell_bc_indexes, cell_umis, total_cells)

    # Cell 0: ATCGATCG + GGGGGGGG = 2 unique
    # Cell 1: TTTTTTTT + CCCCCCCC = 2 unique
    if counts == [2, 2]:
        print(f"  Cell UMI counting: OK (counts={counts})")
        return True
    else:
        print(f"  Cell UMI counting: FAIL (expected [2, 2], got {counts})")
        return False


def main():
    print("Sheriff CI Validation")
    print("=" * 50)

    results = {}

    # Verify input checksums
    results["input_checksums"] = verify_input_checksums()

    # Test Rust imports
    if not test_rust_imports():
        results["rust_import"] = False
        print("\n❌ Rust module not available - cannot continue tests")
        all_pass = False
    else:
        results["rust_import"] = True

        # Run functional tests
        results["kmer_counting"] = test_kmer_counting()
        results["umi_dedup"] = test_umi_deduplication()
        results["edit_clustering"] = test_edit_clustering()
        results["cell_umi_counting"] = test_cell_umi_counting()

    # Summary
    print("\n" + "=" * 50)
    print("Summary")
    print("=" * 50)

    all_pass = all(results.values())
    for test_name, passed in results.items():
        status = "PASS" if passed else "FAIL"
        print(f"  {test_name}: {status}")

    if all_pass:
        print("\n✅ All CI validation tests passed!")
        sys.exit(0)
    else:
        print("\n❌ Some tests failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()
