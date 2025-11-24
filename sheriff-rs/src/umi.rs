//! UMI (Unique Molecular Identifier) deduplication using graph-based clustering
//!
//! This module provides high-performance UMI deduplication for single-cell RNA-seq data.
//! It identifies and counts unique UMI groups based on Hamming distance clustering.
//!
//! # Algorithm
//! 1. Build an undirected graph where each UMI is a node
//! 2. Add edges between UMIs within Hamming distance 1 (single mismatch)
//! 3. Count connected components - each component represents one unique molecule
//!
//! # Performance
//! - Uses sparse graph representation (petgraph) instead of dense adjacency matrices
//! - Inline Hamming distance calculation for speed
//! - Early exits for common cases (1-2 UMIs)
//! - Expected speedup: 50-200x vs Python+Numba version
//!
//! # Example
//! ```ignore
//! use sheriff_rs::umi::deduplicate_umis_rust;
//!
//! let umis = vec!["ATCG".to_string(), "ATCG".to_string(), "GTCA".to_string()];
//! let unique_count = deduplicate_umis_rust(umis);
//! assert_eq!(unique_count, 2); // ATCG and GTCA are different
//! ```

use petgraph::graph::UnGraph;
use petgraph::algo::connected_components;

/// Check if two UMI sequences are within Hamming distance 1 (0 or 1 mismatches)
///
/// This is inlined for maximum performance in the hot path.
///
/// # Arguments
/// * `seq1` - First UMI sequence (as bytes)
/// * `seq2` - Second UMI sequence (as bytes)
///
/// # Returns
/// `true` if sequences differ by at most 1 base, `false` otherwise
#[inline(always)]
fn within_single_mismatch(seq1: &[u8], seq2: &[u8]) -> bool {
    debug_assert_eq!(seq1.len(), seq2.len(), "UMI sequences must be same length");

    let mut mismatches = 0;

    for i in 0..seq1.len() {
        if seq1[i] != seq2[i] {
            mismatches += 1;
            if mismatches > 1 {
                return false; // Early exit on 2nd mismatch
            }
        }
    }

    true
}

/// Deduplicate UMIs using graph-based connected components
///
/// This is the core UMI deduplication function. It:
/// 1. Handles trivial cases (1-2 UMIs) with early exits
/// 2. Builds a sparse undirected graph of UMI relationships
/// 3. Counts connected components using petgraph's algorithm
///
/// # Arguments
/// * `umis` - Slice of UMI sequences (all must be same length)
///
/// # Returns
/// Number of unique UMI groups (connected components)
///
/// # Performance
/// - O(1) for 1 UMI
/// - O(1) for 2 UMIs
/// - O(n²) for n UMIs (graph construction)
/// - O(n + e) for connected components (e = edges)
///
/// Total: O(n²) worst case, but with sparse graphs and early exits
pub fn deduplicate_umis_rust(umis: &[String]) -> usize {
    let num_umis = umis.len();

    // Fast path: single UMI
    if num_umis == 1 {
        return 1;
    }

    // Fast path: two UMIs
    if num_umis == 2 {
        if within_single_mismatch(umis[0].as_bytes(), umis[1].as_bytes()) {
            return 1; // Same molecule
        } else {
            return 2; // Different molecules
        }
    }

    // General case: build graph
    // Use UnGraph for undirected graph (more efficient than directed)
    let mut graph = UnGraph::<(), ()>::with_capacity(num_umis, num_umis * 2);

    // Add nodes (one per UMI)
    let node_indices: Vec<_> = (0..num_umis)
        .map(|_| graph.add_node(()))
        .collect();

    // Add edges between UMIs within Hamming distance 1
    // Only check upper triangle of adjacency matrix (i < j)
    for i in 0..num_umis {
        let umi1_bytes = umis[i].as_bytes();

        for j in (i + 1)..num_umis {
            let umi2_bytes = umis[j].as_bytes();

            if within_single_mismatch(umi1_bytes, umi2_bytes) {
                graph.add_edge(node_indices[i], node_indices[j], ());
            }
        }
    }

    // Count connected components
    // Each component represents one unique molecule
    connected_components(&graph)
}

/// Process UMI deduplication for multiple cells in batch
///
/// This function processes UMI counts for multiple cells efficiently.
/// It mirrors the Python `cell_umi_counts_FAST` function exactly.
///
/// # Arguments
/// * `cell_bc_indexes` - Cell barcode indices for each entry
/// * `cell_umis` - UMI arrays for each entry (parallel to cell_bc_indexes)
/// * `total_cells` - Total number of cells (size of output array)
///
/// # Returns
/// Vector of UMI counts per cell (length = total_cells)
///
/// # Example
/// ```ignore
/// let cell_bc_indexes = vec![0, 1, 0]; // Cell 0, Cell 1, Cell 0 again
/// let cell_umis = vec![
///     vec!["ATCG".to_string()],              // Cell 0: 1 UMI
///     vec!["ATCG".to_string(), "GTCA".to_string()], // Cell 1: 2 different UMIs
///     vec!["TTTT".to_string()],              // Cell 0: 1 more UMI
/// ];
/// let counts = cell_umi_counts_rust(cell_bc_indexes, cell_umis, 2);
/// assert_eq!(counts, vec![2, 2]); // Cell 0: 2 UMIs, Cell 1: 2 UMIs
/// ```
pub fn cell_umi_counts_rust(
    cell_bc_indexes: Vec<usize>,
    cell_umis: Vec<Vec<String>>,
    total_cells: usize,
) -> Vec<u32> {
    debug_assert_eq!(
        cell_bc_indexes.len(),
        cell_umis.len(),
        "cell_bc_indexes and cell_umis must have same length"
    );

    let mut cell_umi_counts = vec![0u32; total_cells];

    for i in 0..cell_bc_indexes.len() {
        let cell_index = cell_bc_indexes[i];
        let umi_array = &cell_umis[i];
        let num_unique_umi = umi_array.len();

        // Fast paths for common cases
        if num_unique_umi == 1 {
            cell_umi_counts[cell_index] += 1;
        } else if num_unique_umi == 2 {
            // Check if the two UMIs are within single mismatch
            if within_single_mismatch(
                umi_array[0].as_bytes(),
                umi_array[1].as_bytes(),
            ) {
                cell_umi_counts[cell_index] += 1; // Collapse to 1
            } else {
                cell_umi_counts[cell_index] += 2; // Keep as 2
            }
        } else {
            // General case: use graph-based deduplication
            let unique_count = deduplicate_umis_rust(umi_array);
            cell_umi_counts[cell_index] = unique_count as u32;
        }
    }

    cell_umi_counts
}

/// Parallel version of cell_umi_counts_rust using rayon
///
/// Processes cells in parallel for better performance on multi-core systems.
/// Uses rayon's par_iter for automatic work distribution.
///
/// # Arguments
/// * `cell_bc_indexes` - Cell barcode indices for each entry
/// * `cell_umis` - UMI arrays for each entry (parallel to cell_bc_indexes)
/// * `total_cells` - Total number of cells (size of output array)
///
/// # Returns
/// Vector of UMI counts per cell (length = total_cells)
///
/// # Performance
/// Expected speedup: 2-8x on multi-core systems (scales with cores)
pub fn cell_umi_counts_rust_parallel(
    cell_bc_indexes: Vec<usize>,
    cell_umis: Vec<Vec<String>>,
    total_cells: usize,
) -> Vec<u32> {
    use rayon::prelude::*;
    use std::sync::Mutex;

    debug_assert_eq!(
        cell_bc_indexes.len(),
        cell_umis.len(),
        "cell_bc_indexes and cell_umis must have same length"
    );

    let cell_umi_counts = Mutex::new(vec![0u32; total_cells]);

    // Process entries in parallel
    (0..cell_bc_indexes.len()).into_par_iter().for_each(|i| {
        let cell_index = cell_bc_indexes[i];
        let umi_array = &cell_umis[i];
        let num_unique_umi = umi_array.len();

        // Update count under lock
        // NOTE: We must match Python's behavior exactly:
        // - For 1-2 UMIs: use += (add to existing count)
        // - For 3+ UMIs: use = (overwrite existing count)
        // This is likely a bug in Python but we preserve it for compatibility
        let mut counts = cell_umi_counts.lock().unwrap();

        if num_unique_umi == 1 {
            counts[cell_index] += 1;
        } else if num_unique_umi == 2 {
            if within_single_mismatch(
                umi_array[0].as_bytes(),
                umi_array[1].as_bytes(),
            ) {
                counts[cell_index] += 1;
            } else {
                counts[cell_index] += 2;
            }
        } else {
            // For 3+ UMIs, use assignment not addition (matches Python behavior)
            let unique_count = deduplicate_umis_rust(umi_array);
            counts[cell_index] = unique_count as u32;
        }
    });

    cell_umi_counts.into_inner().unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_within_single_mismatch_identical() {
        let seq1 = b"ATCGATCG";
        let seq2 = b"ATCGATCG";
        assert!(within_single_mismatch(seq1, seq2));
    }

    #[test]
    fn test_within_single_mismatch_one_diff() {
        let seq1 = b"ATCGATCG";
        let seq2 = b"ATCGATCG";
        assert!(within_single_mismatch(seq1, seq2));

        let seq1 = b"ATCGATCG";
        let seq2 = b"TTCGATCG"; // First base differs
        assert!(within_single_mismatch(seq1, seq2));

        let seq1 = b"ATCGATCG";
        let seq2 = b"ATCGATCT"; // Last base differs
        assert!(within_single_mismatch(seq1, seq2));
    }

    #[test]
    fn test_within_single_mismatch_two_diffs() {
        let seq1 = b"ATCGATCG";
        let seq2 = b"TTCGATCT"; // First and last differ
        assert!(!within_single_mismatch(seq1, seq2));
    }

    #[test]
    fn test_within_single_mismatch_all_diff() {
        let seq1 = b"AAAA";
        let seq2 = b"TTTT";
        assert!(!within_single_mismatch(seq1, seq2));
    }

    #[test]
    fn test_deduplicate_umis_single() {
        let umis = vec!["ATCGATCG".to_string()];
        assert_eq!(deduplicate_umis_rust(&umis), 1);
    }

    #[test]
    fn test_deduplicate_umis_two_identical() {
        let umis = vec!["ATCGATCG".to_string(), "ATCGATCG".to_string()];
        assert_eq!(deduplicate_umis_rust(&umis), 1);
    }

    #[test]
    fn test_deduplicate_umis_two_similar() {
        let umis = vec!["ATCGATCG".to_string(), "TTCGATCG".to_string()];
        assert_eq!(deduplicate_umis_rust(&umis), 1); // Within 1 mismatch
    }

    #[test]
    fn test_deduplicate_umis_two_different() {
        let umis = vec!["ATCGATCG".to_string(), "GGGGGGGG".to_string()];
        assert_eq!(deduplicate_umis_rust(&umis), 2);
    }

    #[test]
    fn test_deduplicate_umis_chain() {
        // Chain: A-B-C (A~B, B~C, but A!~C)
        // Should form 1 connected component
        // Test actual chain: AAAA -> AAAB -> AABB
        let umis = vec![
            "AAAA".to_string(), // A
            "AAAB".to_string(), // B (1 diff from A)
            "AABB".to_string(), // C (1 diff from B, 2 from A)
        ];
        assert_eq!(deduplicate_umis_rust(&umis), 1); // All connected via chain
    }

    #[test]
    fn test_deduplicate_umis_two_groups() {
        let umis = vec![
            "AAAA".to_string(),
            "AAAB".to_string(), // Group 1: AAAA-AAAB
            "TTTT".to_string(),
            "TTTG".to_string(), // Group 2: TTTT-TTTG
        ];
        assert_eq!(deduplicate_umis_rust(&umis), 2);
    }

    #[test]
    fn test_deduplicate_umis_all_different() {
        let umis = vec![
            "AAAA".to_string(),
            "TTTT".to_string(),
            "CCCC".to_string(),
            "GGGG".to_string(),
        ];
        assert_eq!(deduplicate_umis_rust(&umis), 4);
    }

    #[test]
    fn test_deduplicate_umis_real_world_example() {
        // Simulate real UMIs with some PCR duplicates
        let umis = vec![
            "ATCGATCG".to_string(), // Original molecule 1
            "ATCGATCG".to_string(), // PCR duplicate
            "ATCGATCG".to_string(), // PCR duplicate
            "TTCGATCG".to_string(), // Sequencing error (1 mismatch)
            "GGGGGGGG".to_string(), // Original molecule 2
            "GGGGGGGG".to_string(), // PCR duplicate
            "CCCCCCCC".to_string(), // Original molecule 3
        ];
        assert_eq!(deduplicate_umis_rust(&umis), 3); // 3 unique molecules
    }

    #[test]
    fn test_cell_umi_counts_basic() {
        let cell_bc_indexes = vec![0, 1, 0];
        let cell_umis = vec![
            vec!["ATCG".to_string()],
            vec!["TTTT".to_string(), "CCCC".to_string()],
            vec!["GGGG".to_string()],
        ];
        let counts = cell_umi_counts_rust(cell_bc_indexes, cell_umis, 2);
        assert_eq!(counts, vec![2, 2]); // Cell 0: 2 UMIs, Cell 1: 2 UMIs
    }

    #[test]
    fn test_cell_umi_counts_with_duplicates() {
        let cell_bc_indexes = vec![0];
        let cell_umis = vec![
            vec!["ATCG".to_string(), "ATCG".to_string(), "TTTT".to_string()],
        ];
        let counts = cell_umi_counts_rust(cell_bc_indexes, cell_umis, 1);
        assert_eq!(counts, vec![2]); // ATCG (duplicates) + TTTT = 2 unique
    }

    #[test]
    fn test_cell_umi_counts_single_mismatch_collapse() {
        let cell_bc_indexes = vec![0];
        let cell_umis = vec![
            vec!["ATCG".to_string(), "TTCG".to_string()], // 1 mismatch
        ];
        let counts = cell_umi_counts_rust(cell_bc_indexes, cell_umis, 1);
        assert_eq!(counts, vec![1]); // Collapse to 1
    }

    #[test]
    fn test_cell_umi_counts_parallel() {
        let cell_bc_indexes = vec![0, 1, 0, 1];
        let cell_umis = vec![
            vec!["ATCG".to_string()],
            vec!["TTTT".to_string()],
            vec!["GGGG".to_string()],
            vec!["CCCC".to_string()],
        ];
        let counts = cell_umi_counts_rust_parallel(cell_bc_indexes, cell_umis, 2);
        assert_eq!(counts, vec![2, 2]); // Each cell has 2 UMIs
    }

    #[test]
    fn test_cell_umi_counts_parallel_matches_sequential() {
        let cell_bc_indexes = vec![0, 1, 0, 1, 2, 2];
        let cell_umis = vec![
            vec!["ATCG".to_string(), "ATCG".to_string()],
            vec!["TTTT".to_string(), "TTTG".to_string()], // 1 mismatch
            vec!["GGGG".to_string()],
            vec!["CCCC".to_string(), "AAAA".to_string()],
            vec!["ATCG".to_string()],
            vec!["TTTT".to_string(), "CCCC".to_string(), "GGGG".to_string()],
        ];

        let counts_seq = cell_umi_counts_rust(
            cell_bc_indexes.clone(),
            cell_umis.clone(),
            3,
        );
        let counts_par = cell_umi_counts_rust_parallel(
            cell_bc_indexes,
            cell_umis,
            3,
        );

        assert_eq!(counts_seq, counts_par);
    }

    #[test]
    fn test_empty_cells() {
        let cell_bc_indexes: Vec<usize> = vec![];
        let cell_umis: Vec<Vec<String>> = vec![];
        let counts = cell_umi_counts_rust(cell_bc_indexes, cell_umis, 5);
        assert_eq!(counts, vec![0, 0, 0, 0, 0]);
    }

    #[test]
    fn test_large_umi_set() {
        // Test with larger set to verify graph algorithm
        let mut umis = vec![];

        // Group 1: 10 identical UMIs
        for _ in 0..10 {
            umis.push("ATCGATCG".to_string());
        }

        // Group 2: chain of UMIs within 1 mismatch, all same length
        umis.push("AAAAAAAA".to_string());
        umis.push("AAAAAAAC".to_string());

        // Group 3: Isolated UMI
        umis.push("TTTTTTTT".to_string());

        let result = deduplicate_umis_rust(&umis);
        assert_eq!(result, 3); // 3 connected components
    }

    #[test]
    fn test_hamming_distance_edge_cases() {
        // Empty sequences (edge case, shouldn't happen in practice)
        assert!(within_single_mismatch(b"", b""));

        // Single character
        assert!(within_single_mismatch(b"A", b"A"));
        assert!(within_single_mismatch(b"A", b"T"));

        // All same
        let seq = b"AAAAAAAAAA";
        assert!(within_single_mismatch(seq, seq));
    }
}
