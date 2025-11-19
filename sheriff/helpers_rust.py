"""Rust-optimized helper functions for Sheriff.

This module provides drop-in replacements for compute-intensive functions
in helpers.py using Rust implementations via sheriff_rs.
"""

import numpy as np

try:
    import sheriff_rs
    RUST_AVAILABLE = True
except ImportError:
    RUST_AVAILABLE = False
    print("Warning: sheriff_rs not available, falling back to Python implementations")


def deduplicate_umis_rust(umi_set, threshold=1):
    """Rust-optimized UMI deduplication.

    Args:
        umi_set: Set or list of UMI strings
        threshold: Hamming distance threshold for collapsing (default: 1)

    Returns:
        Number of unique UMIs after deduplication
    """
    if not RUST_AVAILABLE:
        from . import helpers
        return len(helpers.deduplicate_umis(umi_set))

    # Convert set to list for Rust
    umi_list = list(umi_set) if isinstance(umi_set, set) else umi_set

    return sheriff_rs.deduplicate_umis(umi_list, threshold=threshold)


def get_cell_counts_from_umi_dict_rust(cell_bc_to_umis, cell_barcodes_dict, threshold=1):
    """Rust-optimized version of get_cell_counts_from_umi_dict.

    Uses parallel Rust implementation for UMI deduplication across cells.

    Args:
        cell_bc_to_umis: Dict mapping cell barcodes to sets of UMIs
        cell_barcodes_dict: Dict mapping cell barcodes to indices
        threshold: Hamming distance threshold for UMI collapsing

    Returns:
        Numpy array of UMI counts per cell
    """
    if not RUST_AVAILABLE:
        from . import helpers
        return helpers.get_cell_counts_from_umi_dict(cell_bc_to_umis, cell_barcodes_dict)

    # Filter to whitelisted cells only
    cell_bc_to_umis = {
        cell_bc: umis
        for cell_bc, umis in cell_bc_to_umis.items()
        if cell_bc in cell_barcodes_dict
    }

    if len(cell_bc_to_umis) == 0:
        return np.zeros(len(cell_barcodes_dict), dtype=np.uint32)

    # Convert sets to lists for Rust
    rust_cells = {
        cell_bc: list(umis)
        for cell_bc, umis in cell_bc_to_umis.items()
    }

    # Use parallel Rust implementation
    if hasattr(sheriff_rs, 'deduplicate_cells_parallel'):
        rust_results = sheriff_rs.deduplicate_cells_parallel(rust_cells, threshold=threshold)
    else:
        # Fall back to sequential Rust
        rust_results = {
            cell_bc: sheriff_rs.deduplicate_umis(umis, threshold=threshold)
            for cell_bc, umis in rust_cells.items()
        }

    # Convert results to numpy array
    cell_umi_counts = np.zeros(len(cell_barcodes_dict), dtype=np.uint32)
    for cell_bc, count in rust_results.items():
        cell_index = cell_barcodes_dict[cell_bc]
        cell_umi_counts[cell_index] = count

    return cell_umi_counts


def match_kmer_rust(sequence, k, whitelist, output_hash=True):
    """Rust-optimized k-mer matching.

    Args:
        sequence: DNA sequence to search
        k: K-mer size
        whitelist: List of k-mer hashes to match against
        output_hash: If True, return hash indices; if False, return k-mer strings

    Returns:
        Tuple of matching k-mer hashes or None if no matches
    """
    if not RUST_AVAILABLE:
        # Fall back to Python implementation
        return None  # Would need to import and use Python version

    matches = sheriff_rs.match_kmer(sequence, k, whitelist, output_hash=output_hash)

    if matches is None or len(matches) == 0:
        return None

    return tuple(matches)


def kmer_to_num_rust(kmer):
    """Rust-optimized k-mer to number conversion.

    Args:
        kmer: DNA sequence string

    Returns:
        Integer hash value for the k-mer
    """
    if not RUST_AVAILABLE:
        # Fall back to Python implementation
        raise NotImplementedError("Rust not available and Python fallback not implemented")

    return sheriff_rs.kmer_to_num(kmer)
