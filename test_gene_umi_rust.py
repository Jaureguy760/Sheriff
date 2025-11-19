#!/usr/bin/env python3
"""
Test and benchmark script for Rust gene UMI counting optimization.

This script:
1. Tests correctness of Rust implementation vs Numba
2. Benchmarks performance improvement
3. Validates integration with Sheriff pipeline
"""

import sys
import time
import numpy as np
from collections import defaultdict
from numba.typed import List

# Import Sheriff helpers
sys.path.insert(0, '/home/user/Sheriff')
from sheriff.helpers import get_cell_by_gene_umi_counts, deduplicated_umi_count_FAST

def generate_test_data(n_genes=100, n_cells=50, avg_umis_per_cell=10, umi_length=12):
    """Generate realistic test data for gene UMI counting."""
    import random
    import string

    def random_umi(length):
        return ''.join(random.choice('ACGT') for _ in range(length))

    gene_indices = []
    gene_cell_indices = []
    gene_cell_umis = []

    for gene_i in range(n_genes):
        # Randomly select cells for this gene
        n_cells_for_gene = random.randint(5, min(20, n_cells))
        cells = sorted(random.sample(range(n_cells), n_cells_for_gene))

        gene_indices.append(gene_i)
        gene_cell_indices.append(np.array(cells, dtype=np.int64))

        # Generate UMIs for each cell
        cell_umis_list = []
        for cell in cells:
            n_umis = random.randint(1, avg_umis_per_cell * 2)
            umis = [random_umi(umi_length) for _ in range(n_umis)]

            # Add some duplicates to test deduplication
            if n_umis > 1 and random.random() < 0.3:
                umis.append(umis[0])  # Exact duplicate

            if n_umis > 1 and random.random() < 0.2:
                # Add UMI with 1 mismatch
                umi_with_error = list(umis[0])
                pos = random.randint(0, umi_length - 1)
                umi_with_error[pos] = random.choice([b for b in 'ACGT' if b != umi_with_error[pos]])
                umis.append(''.join(umi_with_error))

            cell_umis_list.append(np.array(umis, dtype=f'<U{umi_length}'))

        gene_cell_umis.append(List(cell_umis_list))

    gene_indices = np.array(gene_indices, dtype=np.int64)

    return n_cells, gene_indices, gene_cell_indices, gene_cell_umis


def test_correctness():
    """Test that Rust and Numba implementations produce identical results."""
    print("=" * 80)
    print("CORRECTNESS TEST")
    print("=" * 80)

    # Generate test data
    n_cells, gene_indices, gene_cell_indices, gene_cell_umis = generate_test_data(
        n_genes=10, n_cells=20, avg_umis_per_cell=5
    )

    print(f"Test data: {len(gene_indices)} genes, {n_cells} cells")

    # Test Numba implementation
    print("\nTesting Numba implementation...")
    numba_result = get_cell_by_gene_umi_counts(
        n_cells, gene_indices, gene_cell_indices, gene_cell_umis
    )
    print(f"Numba result shape: {numba_result.shape}")
    print(f"Numba non-zero entries: {numba_result.shape[0]}")

    # Test Rust implementation
    print("\nTesting Rust implementation...")
    try:
        import sheriff_rs

        # Convert to Rust format
        gene_indices_list = gene_indices.astype(np.uint32).tolist()
        gene_cell_indices_list = [arr.astype(np.uint32).tolist() for arr in gene_cell_indices]
        gene_cell_umis_list = [[umi_arr.tolist() for umi_arr in cell_umis] for cell_umis in gene_cell_umis]

        rust_result = sheriff_rs.count_gene_umis_rust(
            n_cells, gene_indices_list, gene_cell_indices_list, gene_cell_umis_list, 1
        )
        rust_result = np.array(rust_result, dtype=np.uint32)

        print(f"Rust result shape: {rust_result.shape}")
        print(f"Rust non-zero entries: {rust_result.shape[0]}")

        # Compare results
        print("\nComparing results...")
        if numba_result.shape != rust_result.shape:
            print(f"❌ Shape mismatch: Numba {numba_result.shape} vs Rust {rust_result.shape}")
            return False

        # Sort both results for comparison (order may differ)
        numba_sorted = numba_result[np.lexsort((numba_result[:, 2], numba_result[:, 1], numba_result[:, 0]))]
        rust_sorted = rust_result[np.lexsort((rust_result[:, 2], rust_result[:, 1], rust_result[:, 0]))]

        if not np.array_equal(numba_sorted, rust_sorted):
            print("❌ Results differ!")
            print("\nFirst 10 mismatches:")
            mask = ~np.all(numba_sorted == rust_sorted, axis=1)
            mismatches = np.where(mask)[0][:10]
            for i in mismatches:
                print(f"  Row {i}: Numba {numba_sorted[i]} vs Rust {rust_sorted[i]}")
            return False

        print("✅ Results match perfectly!")
        return True

    except ImportError as e:
        print(f"❌ Rust module not available: {e}")
        return False


def benchmark_performance():
    """Benchmark Rust vs Numba performance."""
    print("\n" + "=" * 80)
    print("PERFORMANCE BENCHMARK")
    print("=" * 80)

    test_sizes = [
        (50, 100, 10, "Small"),
        (200, 500, 20, "Medium"),
        (500, 1000, 30, "Large"),
    ]

    results = []

    for n_genes, n_cells, avg_umis, label in test_sizes:
        print(f"\n{label} dataset: {n_genes} genes, {n_cells} cells, ~{avg_umis} UMIs/cell")
        print("-" * 80)

        # Generate data
        total_cells, gene_indices, gene_cell_indices, gene_cell_umis = generate_test_data(
            n_genes=n_genes, n_cells=n_cells, avg_umis_per_cell=avg_umis
        )

        # Benchmark Numba
        print("Benchmarking Numba...")
        numba_times = []
        for i in range(3):
            start = time.perf_counter()
            numba_result = get_cell_by_gene_umi_counts(
                total_cells, gene_indices, gene_cell_indices, gene_cell_umis
            )
            elapsed = time.perf_counter() - start
            numba_times.append(elapsed)
            print(f"  Run {i+1}: {elapsed*1000:.2f} ms")

        numba_avg = np.mean(numba_times)
        numba_std = np.std(numba_times)
        print(f"  Average: {numba_avg*1000:.2f} ± {numba_std*1000:.2f} ms")

        # Benchmark Rust
        try:
            import sheriff_rs

            print("\nBenchmarking Rust...")
            # Convert data once (not included in benchmark)
            gene_indices_list = gene_indices.astype(np.uint32).tolist()
            gene_cell_indices_list = [arr.astype(np.uint32).tolist() for arr in gene_cell_indices]
            gene_cell_umis_list = [[umi_arr.tolist() for umi_arr in cell_umis] for cell_umis in gene_cell_umis]

            rust_times = []
            for i in range(3):
                start = time.perf_counter()
                rust_result = sheriff_rs.count_gene_umis_rust(
                    total_cells, gene_indices_list, gene_cell_indices_list, gene_cell_umis_list, 1
                )
                rust_result_array = np.array(rust_result, dtype=np.uint32)
                elapsed = time.perf_counter() - start
                rust_times.append(elapsed)
                print(f"  Run {i+1}: {elapsed*1000:.2f} ms")

            rust_avg = np.mean(rust_times)
            rust_std = np.std(rust_times)
            print(f"  Average: {rust_avg*1000:.2f} ± {rust_std*1000:.2f} ms")

            speedup = numba_avg / rust_avg
            print(f"\n📊 Speedup: {speedup:.2f}x")

            results.append({
                'label': label,
                'n_genes': n_genes,
                'n_cells': n_cells,
                'numba_ms': numba_avg * 1000,
                'rust_ms': rust_avg * 1000,
                'speedup': speedup
            })

        except ImportError as e:
            print(f"\n❌ Rust module not available: {e}")

    # Print summary
    if results:
        print("\n" + "=" * 80)
        print("BENCHMARK SUMMARY")
        print("=" * 80)
        print(f"{'Dataset':<10} {'Numba (ms)':<12} {'Rust (ms)':<12} {'Speedup':<10}")
        print("-" * 80)
        for r in results:
            print(f"{r['label']:<10} {r['numba_ms']:>10.2f}   {r['rust_ms']:>10.2f}   {r['speedup']:>8.2f}x")

        avg_speedup = np.mean([r['speedup'] for r in results])
        print("-" * 80)
        print(f"Average speedup: {avg_speedup:.2f}x")

        # Estimate impact on full pipeline
        # Gene UMI counting is ~5.25% of total runtime
        pipeline_impact = (1 - 1/avg_speedup) * 5.25
        print(f"Estimated pipeline speedup: +{pipeline_impact:.2f}% (if gene UMI counting is 5.25% of runtime)")


if __name__ == "__main__":
    print("Gene UMI Counting Rust Optimization - Test & Benchmark")
    print("=" * 80)

    # Test correctness
    if test_correctness():
        print("\n✅ Correctness test passed!")
    else:
        print("\n❌ Correctness test failed!")
        sys.exit(1)

    # Benchmark performance
    benchmark_performance()

    print("\n" + "=" * 80)
    print("Testing complete!")
