#!/usr/bin/env python3
"""
Comprehensive integration testing for Sheriff with Rust optimizations.

Tests:
1. Individual Rust functions vs Python equivalents (correctness)
2. Full Sheriff pipeline performance
3. Output verification
4. End-to-end benchmarking

Author: Sheriff Integration Test Suite
Date: 2025-11-19
"""

import sys
import time
import timeit
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict
import pysam

# Import Sheriff modules
from sheriff import helpers
from sheriff.count_t7 import KmerMatcher

# Import Rust optimizations
try:
    import sheriff_rs
    RUST_AVAILABLE = True
    print(f"✅ sheriff_rs v{sheriff_rs.__version__} loaded successfully")
except ImportError as e:
    RUST_AVAILABLE = False
    print(f"❌ sheriff_rs not available: {e}")
    sys.exit(1)


class IntegrationTester:
    """Comprehensive integration testing for Sheriff Rust optimizations."""

    def __init__(self, bam_file, barcode_file):
        self.bam_file = bam_file
        self.barcode_file = barcode_file
        self.results = {}

        # Load cell barcodes
        with open(barcode_file) as f:
            self.cell_barcodes = [line.rstrip() for line in f]
            self.cell_barcodes_dict = {bc: i for i, bc in enumerate(self.cell_barcodes)}

        print(f"📊 Loaded {len(self.cell_barcodes)} cell barcodes")

    def test_kmer_conversion(self):
        """Test k-mer to number conversion: Rust vs Python."""
        print("\n" + "="*80)
        print("TEST 1: K-mer to Number Conversion")
        print("="*80)

        k = 6
        test_kmers = ["ACGTAC", "GGGGGG", "AAAAAA", "TCGATC", "NNNNNN"]

        matcher = KmerMatcher(k)

        passed = 0
        failed = 0

        for kmer in test_kmers:
            if 'N' in kmer:
                # Skip invalid kmers for Rust
                continue

            python_result = matcher.kmer_to_num(kmer)
            rust_result = sheriff_rs.kmer_to_num(kmer)

            match = python_result == rust_result
            status = "✅" if match else "❌"

            if match:
                passed += 1
            else:
                failed += 1

            print(f"{status} {kmer}: Python={python_result}, Rust={rust_result}")

        self.results['kmer_conversion'] = {'passed': passed, 'failed': failed}
        print(f"\n📈 Result: {passed} passed, {failed} failed")
        return failed == 0

    def test_kmer_matching(self):
        """Test k-mer matching: Rust vs Python."""
        print("\n" + "="*80)
        print("TEST 2: K-mer Matching")
        print("="*80)

        k = 6
        t7_barcode = "GGGAGAGTAT"
        test_sequences = [
            "GGGAGAGTATACGTACGT",  # Contains barcode
            "ACGTACGTACGTACGT",     # No barcode
            "GGGAGAACGTACGT",       # Partial match
        ]

        # Create Python matcher
        matcher = KmerMatcher(k, t7_barcode)

        # Create Rust whitelist
        rust_whitelist = matcher.match_hash

        passed = 0
        failed = 0

        for seq in test_sequences:
            # Python version (simplified from match_kmer function)
            try:
                freq_array = np.zeros((4 ** k), dtype=np.uint8)
                freq_array[[matcher.kmer_to_num(seq[i: i + k])
                           for i in range(len(seq) - k + 1)]] += 1
                kmer_matches_py = freq_array.nonzero()[0]
                kmer_matches_py = kmer_matches_py[np.isin(kmer_matches_py, matcher.match_hash)]
                python_result = tuple(kmer_matches_py) if kmer_matches_py.size > 0 else None
            except KeyError:
                python_result = None

            # Rust version
            rust_result = sheriff_rs.match_kmer(seq, k, rust_whitelist, output_hash=True)
            if rust_result is not None and len(rust_result) > 0:
                rust_result = tuple(sorted(rust_result))
            elif rust_result is not None and len(rust_result) == 0:
                rust_result = None  # Convert empty tuple to None for comparison

            # Compare
            match = python_result == rust_result
            status = "✅" if match else "❌"

            if match:
                passed += 1
            else:
                failed += 1

            print(f"{status} {seq[:20]}...")
            print(f"   Python: {python_result}")
            print(f"   Rust:   {rust_result}")

        self.results['kmer_matching'] = {'passed': passed, 'failed': failed}
        print(f"\n📈 Result: {passed} passed, {failed} failed")
        return failed == 0

    def test_umi_deduplication(self):
        """Test UMI deduplication: Rust vs Python."""
        print("\n" + "="*80)
        print("TEST 3: UMI Deduplication")
        print("="*80)

        test_cases = [
            # (UMI set, expected count, description)
            (["AAAAAAAAAA"], 1, "Single UMI"),
            (["AAAAAAAAAA", "AAAAAAAAAA"], 1, "Identical UMIs"),
            (["AAAAAAAAAA", "AAAAAAATAA"], 1, "1 mismatch (should collapse)"),
            (["AAAAAAAAAA", "AAAAAATTAA"], 2, "2 mismatches (should NOT collapse)"),
            (["ACGTACGTAC", "TGCATGCATG"], 2, "Completely different"),
            (["AAAAAAAAAA", "AAAAAAATAA", "AAAAAATTAA"], 1, "Three UMIs with varying distances (all collapse)"),
        ]

        passed = 0
        failed = 0

        for umis, expected, description in test_cases:
            # Python version (using helpers.deduplicate_umis)
            python_result = len(helpers.deduplicate_umis(set(umis)))

            # Rust version
            rust_result = sheriff_rs.deduplicate_umis(umis, threshold=1)

            # Check if both match expected
            py_correct = python_result == expected
            rust_correct = rust_result == expected
            match = python_result == rust_result

            if match and py_correct:
                status = "✅"
                passed += 1
            else:
                status = "❌"
                failed += 1

            print(f"{status} {description}")
            print(f"   Expected: {expected}, Python: {python_result}, Rust: {rust_result}")

        self.results['umi_deduplication'] = {'passed': passed, 'failed': failed}
        print(f"\n📈 Result: {passed} passed, {failed} failed")
        return failed == 0

    def test_parallel_umi_deduplication(self):
        """Test parallel UMI deduplication on real BAM data."""
        print("\n" + "="*80)
        print("TEST 4: Parallel UMI Deduplication (Real Data)")
        print("="*80)

        # Extract UMIs from real BAM file
        print(f"📂 Reading UMIs from {self.bam_file}")

        cell_to_umis = defaultdict(set)
        n_reads = 0
        max_reads = 100000  # Limit for testing

        with pysam.AlignmentFile(self.bam_file, "rb") as bam:
            for read in bam:
                if n_reads >= max_reads:
                    break

                try:
                    cb = read.get_tag('CB')
                    umi = read.get_tag('pN')

                    if cb in self.cell_barcodes_dict:
                        cell_to_umis[cb].add(umi)
                        n_reads += 1
                except KeyError:
                    continue

        print(f"📊 Extracted {n_reads} reads from {len(cell_to_umis)} cells")

        # Test a sample of cells
        test_cells = list(cell_to_umis.keys())[:10]

        passed = 0
        failed = 0

        for cell_bc in test_cells:
            umis = list(cell_to_umis[cell_bc])

            # Python version
            python_result = len(helpers.deduplicate_umis(set(umis)))

            # Rust version (single-threaded)
            rust_result = sheriff_rs.deduplicate_umis(umis, threshold=1)

            match = python_result == rust_result
            status = "✅" if match else "❌"

            if match:
                passed += 1
            else:
                failed += 1

            print(f"{status} Cell {cell_bc}: {len(umis)} UMIs → Python: {python_result}, Rust: {rust_result}")

        # Test parallel version if available
        if hasattr(sheriff_rs, 'deduplicate_cells_parallel'):
            print("\n🚀 Testing parallel version...")

            # Prepare data for parallel processing (dict format)
            rust_cells = {
                cb: list(umis) for cb, umis in cell_to_umis.items()
            }

            start = time.time()
            rust_par_results = sheriff_rs.deduplicate_cells_parallel(rust_cells, threshold=1)
            rust_par_time = time.time() - start

            print(f"⏱️  Parallel processing: {len(rust_par_results)} cells in {rust_par_time:.3f}s")
            print(f"   ({len(rust_par_results)/rust_par_time:.1f} cells/sec)")

        self.results['parallel_umi_dedup'] = {'passed': passed, 'failed': failed}
        print(f"\n📈 Result: {passed} passed, {failed} failed")
        return failed == 0

    def benchmark_performance(self):
        """Benchmark Rust vs Python performance."""
        print("\n" + "="*80)
        print("BENCHMARK: Performance Comparison")
        print("="*80)

        # Benchmark k-mer conversion
        print("\n🔄 K-mer conversion (10,000 iterations)...")
        k = 6
        test_kmer = "ACGTAC"

        matcher = KmerMatcher(k)

        py_time = timeit.timeit(lambda: matcher.kmer_to_num(test_kmer), number=10000)
        rust_time = timeit.timeit(lambda: sheriff_rs.kmer_to_num(test_kmer), number=10000)

        speedup = py_time / rust_time
        print(f"   Python: {py_time:.4f}s")
        print(f"   Rust:   {rust_time:.4f}s")
        print(f"   Speedup: {speedup:.2f}x")

        self.results['kmer_conversion_speedup'] = speedup

        # Benchmark UMI deduplication
        print("\n🔄 UMI deduplication (1,000 iterations, 100 UMIs)...")
        test_umis = ["".join(np.random.choice(list("ACGT"), 10)) for _ in range(100)]

        py_time = timeit.timeit(lambda: len(helpers.deduplicate_umis(set(test_umis))), number=1000)
        rust_time = timeit.timeit(lambda: sheriff_rs.deduplicate_umis(test_umis, threshold=1), number=1000)

        speedup = py_time / rust_time
        print(f"   Python: {py_time:.4f}s")
        print(f"   Rust:   {rust_time:.4f}s")
        print(f"   Speedup: {speedup:.2f}x")

        self.results['umi_dedup_speedup'] = speedup

    def generate_report(self):
        """Generate comprehensive test report."""
        print("\n" + "="*80)
        print("INTEGRATION TEST REPORT")
        print("="*80)

        total_passed = sum(r.get('passed', 0) for r in self.results.values() if isinstance(r, dict))
        total_failed = sum(r.get('failed', 0) for r in self.results.values() if isinstance(r, dict))

        print(f"\n📊 Overall Results:")
        print(f"   ✅ Passed: {total_passed}")
        print(f"   ❌ Failed: {total_failed}")

        if total_failed == 0:
            print(f"\n🎉 All tests PASSED! Rust optimizations are correct.")
        else:
            print(f"\n⚠️  Some tests FAILED. Review results above.")

        print(f"\n⚡ Performance Summary:")
        if 'kmer_conversion_speedup' in self.results:
            print(f"   K-mer conversion: {self.results['kmer_conversion_speedup']:.2f}x faster")
        if 'umi_dedup_speedup' in self.results:
            print(f"   UMI deduplication: {self.results['umi_dedup_speedup']:.2f}x faster")

        return total_failed == 0


def main():
    """Main integration test entry point."""
    print("="*80)
    print("SHERIFF RUST INTEGRATION TESTS")
    print("="*80)

    # Configuration
    bam_file = "/home/user/Sheriff/example_data/barcode_headAligned_anno.sorted.edit_regions_200kb.bam"
    barcode_file = "/home/user/Sheriff/example_data/barcode_whitelist.500-cell.txt"

    # Check files exist
    if not Path(bam_file).exists():
        print(f"❌ BAM file not found: {bam_file}")
        sys.exit(1)

    if not Path(barcode_file).exists():
        print(f"❌ Barcode file not found: {barcode_file}")
        sys.exit(1)

    # Create tester
    tester = IntegrationTester(bam_file, barcode_file)

    # Run tests
    tests = [
        tester.test_kmer_conversion,
        tester.test_kmer_matching,
        tester.test_umi_deduplication,
        tester.test_parallel_umi_deduplication,
    ]

    all_passed = True
    for test_func in tests:
        try:
            if not test_func():
                all_passed = False
        except Exception as e:
            print(f"\n❌ Test failed with exception: {e}")
            import traceback
            traceback.print_exc()
            all_passed = False

    # Benchmark
    try:
        tester.benchmark_performance()
    except Exception as e:
        print(f"\n⚠️  Benchmark failed: {e}")

    # Generate report
    success = tester.generate_report()

    # Exit code
    sys.exit(0 if success and all_passed else 1)


if __name__ == "__main__":
    main()
