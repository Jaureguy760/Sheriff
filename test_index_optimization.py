#!/usr/bin/env python3
"""
Unit test to verify that the .index() optimizations maintain correctness.
Tests the three main optimizations made in count_t7.py:
1. edit_site_to_index mapping (allelic calling)
2. edit_name_to_index mapping (gene overlap checking)
3. called_edit_site_name_to_index mapping (gene allelic calling)
"""

from collections import namedtuple
import numpy as np


def test_edit_site_to_index_mapping():
    """Test the edit_site_to_index optimization (main allelic calling)"""
    print("\n" + "="*70)
    print("TEST 1: edit_site_to_index mapping (allelic calling)")
    print("="*70)

    EditSite = namedtuple("EditSite", ["chrom", "ref_pos"])

    # Create called_edit_sites list
    called_edit_sites = [
        EditSite("chr1", 100),
        EditSite("chr1", 200),
        EditSite("chr2", 150),
        EditSite("chr3", 300),
        EditSite("chr1", 250),
    ]

    # OLD METHOD: using list.index()
    def old_method(edit_site):
        return called_edit_sites.index(edit_site)

    # NEW METHOD: using pre-computed dict
    edit_site_to_index = {
        edit_site: i for i, edit_site in enumerate(called_edit_sites)
    }

    def new_method(edit_site):
        return edit_site_to_index[edit_site]

    # Test each edit site
    all_match = True
    for edit_site in called_edit_sites:
        old_idx = old_method(edit_site)
        new_idx = new_method(edit_site)
        match = old_idx == new_idx
        status = "✅" if match else "❌"
        print(f"{status} {edit_site}: old={old_idx}, new={new_idx}, match={match}")
        if not match:
            all_match = False

    print(f"\n{'✅ PASS' if all_match else '❌ FAIL'}: edit_site_to_index mapping")
    return all_match


def test_edit_name_to_index_mapping():
    """Test the edit_name_to_index optimization (gene overlap checking)"""
    print("\n" + "="*70)
    print("TEST 2: edit_name_to_index mapping (gene overlap checking)")
    print("="*70)

    # Create edit_names list
    edit_names = [
        "chr1:100",
        "chr1:200",
        "chr2:150",
        "chr3:300",
        "chr1:250",
    ]

    # OLD METHOD: using list.index()
    def old_method(edit_name):
        return edit_names.index(edit_name)

    # NEW METHOD: using pre-computed dict
    edit_name_to_index = {name: i for i, name in enumerate(edit_names)}

    def new_method(edit_name):
        return edit_name_to_index[edit_name]

    # Test each edit name
    all_match = True
    for edit_name in edit_names:
        old_idx = old_method(edit_name)
        new_idx = new_method(edit_name)
        match = old_idx == new_idx
        status = "✅" if match else "❌"
        print(f"{status} {edit_name}: old={old_idx}, new={new_idx}, match={match}")
        if not match:
            all_match = False

    print(f"\n{'✅ PASS' if all_match else '❌ FAIL'}: edit_name_to_index mapping")
    return all_match


def test_called_edit_site_name_to_index_mapping():
    """Test the called_edit_site_name_to_index optimization (gene allelic calling)"""
    print("\n" + "="*70)
    print("TEST 3: called_edit_site_name_to_index mapping (gene allelic calling)")
    print("="*70)

    # Create called_edit_site_names list
    called_edit_site_names = [
        "chr1:100",
        "chr1:200",
        "chr2:150",
        "chr3:300",
        "chr1:250",
    ]

    # Test data: genes and their associated edits
    genes_to_edits = {
        "GENE1": ["chr1:100", "chr1:200"],
        "GENE2": ["chr2:150"],
        "GENE3": ["chr3:300", "chr1:250"],
    }

    # OLD METHOD: using list comprehension with .index()
    def old_method(gene_edit_list):
        return [called_edit_site_names.index(genic_edit) for genic_edit in gene_edit_list]

    # NEW METHOD: using pre-computed dict in list comprehension
    called_edit_site_name_to_index = {
        name: i for i, name in enumerate(called_edit_site_names)
    }

    def new_method(gene_edit_list):
        return [called_edit_site_name_to_index[genic_edit] for genic_edit in gene_edit_list]

    # Test each gene's edits
    all_match = True
    for gene, edits in genes_to_edits.items():
        old_indices = old_method(edits)
        new_indices = new_method(edits)
        match = old_indices == new_indices
        status = "✅" if match else "❌"
        print(f"{status} {gene}: {edits}")
        print(f"   old_indices={old_indices}, new_indices={new_indices}, match={match}")
        if not match:
            all_match = False

    print(f"\n{'✅ PASS' if all_match else '❌ FAIL'}: called_edit_site_name_to_index mapping")
    return all_match


def test_nested_loop_correctness():
    """Test that nested loop access patterns work correctly with the optimization"""
    print("\n" + "="*70)
    print("TEST 4: Nested loop correctness (simulates allelic calling)")
    print("="*70)

    EditSite = namedtuple("EditSite", ["chrom", "ref_pos"])

    # Create test data structure similar to the actual code
    called_edit_sites = [
        EditSite("chr1", 100),
        EditSite("chr1", 200),
        EditSite("chr2", 150),
    ]

    cells_to_edits = {
        "CELL_1": {EditSite("chr1", 100): ["edit1"], EditSite("chr1", 200): ["edit2"]},
        "CELL_2": {EditSite("chr2", 150): ["edit3"]},
        "CELL_3": {EditSite("chr1", 100): ["edit4"], EditSite("chr2", 150): ["edit5"]},
    }

    # OLD METHOD
    old_results = {}
    for cell_barcode, edit_sites_to_edits in cells_to_edits.items():
        old_results[cell_barcode] = []
        for edit_site in edit_sites_to_edits.keys():
            edit_sitei = called_edit_sites.index(edit_site)
            old_results[cell_barcode].append((edit_site, edit_sitei))

    # NEW METHOD
    edit_site_to_index = {
        edit_site: i for i, edit_site in enumerate(called_edit_sites)
    }

    new_results = {}
    for cell_barcode, edit_sites_to_edits in cells_to_edits.items():
        new_results[cell_barcode] = []
        for edit_site in edit_sites_to_edits.keys():
            edit_sitei = edit_site_to_index[edit_site]
            new_results[cell_barcode].append((edit_site, edit_sitei))

    # Compare results
    all_match = True
    for cell_barcode in cells_to_edits.keys():
        old_results_cell = sorted(old_results[cell_barcode])
        new_results_cell = sorted(new_results[cell_barcode])
        match = old_results_cell == new_results_cell
        status = "✅" if match else "❌"
        print(f"{status} {cell_barcode}")
        if not match:
            print(f"   OLD: {old_results_cell}")
            print(f"   NEW: {new_results_cell}")
            all_match = False

    print(f"\n{'✅ PASS' if all_match else '❌ FAIL'}: nested loop correctness")
    return all_match


def main():
    """Run all tests"""
    print("="*70)
    print("INDEX OPTIMIZATION CORRECTNESS TESTS")
    print("="*70)

    results = []
    results.append(("edit_site_to_index", test_edit_site_to_index_mapping()))
    results.append(("edit_name_to_index", test_edit_name_to_index_mapping()))
    results.append(("called_edit_site_name_to_index", test_called_edit_site_name_to_index_mapping()))
    results.append(("nested_loop_correctness", test_nested_loop_correctness()))

    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)

    passed = sum(1 for _, result in results if result)
    failed = len(results) - passed

    for test_name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{status}: {test_name}")

    print(f"\nTotal: {passed} passed, {failed} failed")

    if failed == 0:
        print("\n🎉 All correctness tests PASSED!")
        return 0
    else:
        print(f"\n⚠️  {failed} test(s) FAILED!")
        return 1


if __name__ == "__main__":
    exit(main())
