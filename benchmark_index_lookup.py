#!/usr/bin/env python3
"""
Benchmark script to measure the performance improvement from fixing .index() lookups
in the allelic calling code.
"""

import timeit
import numpy as np

def benchmark_index_lookup():
    """Benchmark the difference between list.index() and dict lookup"""

    # Simulate the data structure from the allelic calling code
    # Create a list of edit sites (simulated as namedtuples)
    from collections import namedtuple
    EditSite = namedtuple("EditSite", ["chrom", "ref_pos"])

    # Generate synthetic data
    num_edit_sites = 1000
    num_cells = 10000
    edits_per_cell = 5  # Average number of edits per cell

    # Create called_edit_sites list
    called_edit_sites = [
        EditSite(f"chr{i//100}", i % 100)
        for i in range(num_edit_sites)
    ]

    # Create a mapping of cell to edit sites (simulating cells_to_canonical_and_edits)
    cells_to_edits = {}
    for cell_idx in range(num_cells):
        cell_barcode = f"CELL_{cell_idx:06d}"
        num_edits = np.random.randint(1, edits_per_cell + 1)
        edit_indices = np.random.choice(num_edit_sites, num_edits, replace=False)
        cells_to_edits[cell_barcode] = {
            called_edit_sites[idx]: [f"edit_{idx}_var1"]
            for idx in edit_indices
        }

    print("=" * 70)
    print(f"Benchmark: .index() lookup vs dict lookup")
    print(f"  Edit sites: {len(called_edit_sites)}")
    print(f"  Cells: {len(cells_to_edits)}")
    print(f"  Total lookups: {sum(len(edits) for edits in cells_to_edits.values())}")
    print("=" * 70)

    # OLD METHOD: Using list.index()
    def old_method():
        cell_count = 0
        for cell_barcode, edit_sites_to_edits in cells_to_edits.items():
            cell_count += 1
            for edit_site in edit_sites_to_edits.keys():
                # This is O(n) per lookup
                edit_sitei = called_edit_sites.index(edit_site)
        return cell_count

    # NEW METHOD: Using pre-computed dict
    def new_method():
        # Pre-compute the index mapping (O(n) one-time cost)
        edit_site_to_index = {
            edit_site: i for i, edit_site in enumerate(called_edit_sites)
        }

        cell_count = 0
        for cell_barcode, edit_sites_to_edits in cells_to_edits.items():
            cell_count += 1
            for edit_site in edit_sites_to_edits.keys():
                # This is O(1) lookup
                edit_sitei = edit_site_to_index[edit_site]
        return cell_count

    # Run benchmarks
    num_iterations = 5

    print(f"\nRunning OLD method (list.index()) {num_iterations} times...")
    old_time = timeit.timeit(old_method, number=num_iterations)
    old_avg = old_time / num_iterations
    print(f"  Total time: {old_time:.4f}s")
    print(f"  Average: {old_avg:.4f}s per iteration")

    print(f"\nRunning NEW method (dict lookup) {num_iterations} times...")
    new_time = timeit.timeit(new_method, number=num_iterations)
    new_avg = new_time / num_iterations
    print(f"  Total time: {new_time:.4f}s")
    print(f"  Average: {new_avg:.4f}s per iteration")

    print("\n" + "=" * 70)
    speedup = old_avg / new_avg
    improvement = (old_avg - new_avg) / old_avg * 100
    print(f"Speedup: {speedup:.2f}x")
    print(f"Time saved: {improvement:.2f}%")
    print("=" * 70)

if __name__ == "__main__":
    benchmark_index_lookup()
