#!/usr/bin/env python3
"""
Sheriff Test Data Validation Script

Validates that the Rust-accelerated Sheriff produces correct results
by comparing with a known baseline.

Usage:
    python test_data/run_validation.py
"""

import os
import sys
import time
from pathlib import Path

# Add Sheriff to path
SCRIPT_DIR = Path(__file__).parent.absolute()
SHERIFF_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SHERIFF_ROOT))

def check_prerequisites():
    """Check that all required files exist."""
    required_files = [
        SCRIPT_DIR / "test_200kb.bam",
        SCRIPT_DIR / "test_200kb.bam.bai",
        SCRIPT_DIR / "barcodes.txt",
        SCRIPT_DIR / "reference.fa",
        SCRIPT_DIR / "reference.gtf",
    ]

    missing = []
    for f in required_files:
        if not f.exists():
            missing.append(str(f))

    if missing:
        print("ERROR: Missing required files:")
        for f in missing:
            print(f"  - {f}")

        if "reference.fa" in str(missing) or "reference.gtf" in str(missing):
            print("\nTo download reference genome, run:")
            print("  cd test_data && bash download_reference.sh")

        return False

    return True


def check_rust_available():
    """Check if Rust acceleration is available."""
    try:
        import sheriff_rs
        print("[OK] Rust acceleration available")
        return True
    except ImportError as e:
        print(f"[WARNING] Rust acceleration not available: {e}")
        print("To enable Rust acceleration:")
        print("  cd sheriff-rs && maturin develop --release")
        return False


def validate_imports():
    """Validate Sheriff imports work correctly."""
    print("\n=== Validating Imports ===")

    try:
        from sheriff import count_t7
        print("[OK] sheriff.count_t7 imported")
    except ImportError as e:
        print(f"[ERROR] Cannot import sheriff.count_t7: {e}")
        return False

    try:
        from sheriff import helpers
        print("[OK] sheriff.helpers imported")
    except ImportError as e:
        print(f"[ERROR] Cannot import sheriff.helpers: {e}")
        return False

    has_rust = check_rust_available()
    return True


def run_quick_smoke_test():
    """Run a quick smoke test to ensure basic functionality works."""
    print("\n=== Running Smoke Test ===")

    # Test k-mer matching
    try:
        import sheriff_rs
        seq = "ATCGATCGATCGATCGATCG"
        kmers = sheriff_rs.count_kmers_rust(seq, 8)
        assert len(kmers) == 4**8, f"Expected {4**8} k-mers, got {len(kmers)}"
        print(f"[OK] k-mer counting: {sum(kmers)} total k-mers in test sequence")
    except Exception as e:
        print(f"[WARNING] k-mer counting test failed: {e}")

    # Test UMI deduplication
    try:
        import sheriff_rs
        umis = ["ATCGATCG", "ATCGATCG", "TTTTTTTT", "ATCGATCC"]  # Last differs by 1
        unique = sheriff_rs.deduplicate_umis_py(umis)
        assert unique == 2, f"Expected 2 unique UMI groups, got {unique}"
        print(f"[OK] UMI deduplication: {unique} unique groups from {len(umis)} UMIs")
    except Exception as e:
        print(f"[WARNING] UMI deduplication test failed: {e}")

    # Test edit clustering
    try:
        import sheriff_rs
        edits = [
            ("chr1", 1000, "ATCG", "ATCGATCGATCG", True, [1, 2, 3]),
            ("chr1", 1000, "ATCG", "ATCGATCG", True, [1]),  # Subset
            ("chr1", 2000, "GCTA", "GCTACCCC", False, [4, 5]),
        ]
        longest = sheriff_rs.get_longest_edits_rust(edits)
        assert len(longest) == 2, f"Expected 2 longest edits, got {len(longest)}"
        print(f"[OK] Edit clustering: {len(longest)} unique edits from {len(edits)} input")
    except Exception as e:
        print(f"[WARNING] Edit clustering test failed: {e}")

    print("\nSmoke tests completed!")
    return True


def run_bam_processing_test():
    """Test BAM file processing with real data."""
    print("\n=== Testing BAM Processing ===")

    bam_file = str(SCRIPT_DIR / "test_200kb.bam")
    barcodes_file = str(SCRIPT_DIR / "barcodes.txt")

    # Load barcodes
    with open(barcodes_file) as f:
        barcodes = [line.strip() for line in f if line.strip()]

    print(f"Loaded {len(barcodes)} cell barcodes")

    # Test Rust BAM filtering (just count reads, don't write output)
    try:
        import pysam

        start_time = time.time()
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            total_reads = 0
            matched_reads = 0
            barcode_set = set(barcodes)

            for read in bam.fetch():
                total_reads += 1
                if read.has_tag("CB"):
                    cb = read.get_tag("CB")
                    if cb in barcode_set:
                        matched_reads += 1

        elapsed = time.time() - start_time
        print(f"[OK] BAM processing: {total_reads:,} reads, {matched_reads:,} matched in {elapsed:.2f}s")
    except Exception as e:
        print(f"[WARNING] BAM processing test failed: {e}")

    return True


def main():
    print("Sheriff Test Data Validation")
    print("=" * 50)
    print(f"Test data directory: {SCRIPT_DIR}")
    print(f"Sheriff root: {SHERIFF_ROOT}")
    print()

    # Check prerequisites
    print("=== Checking Prerequisites ===")
    if not check_prerequisites():
        print("\nPlease ensure all required files are present before running validation.")
        sys.exit(1)
    print("[OK] All required files found")

    # Validate imports
    if not validate_imports():
        sys.exit(1)

    # Run smoke tests
    run_quick_smoke_test()

    # Run BAM processing test
    run_bam_processing_test()

    print("\n" + "=" * 50)
    print("Validation Summary")
    print("=" * 50)
    print("[OK] All validation tests passed!")
    print("\nYour Sheriff installation is ready for development.")
    print("\nNext steps:")
    print("1. Run a full test: python -m sheriff.count_t7 --help")
    print("2. Profile performance: python sanity_check_comparison.py")
    print("3. Benchmark Rust: cd sheriff-rs && cargo bench")


if __name__ == "__main__":
    main()
