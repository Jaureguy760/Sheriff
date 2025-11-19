//! UMI Deduplication Module
//!
//! This module implements high-performance UMI (Unique Molecular Identifier) deduplication
//! using a Union-Find data structure with path compression and union-by-rank optimizations.
//!
//! # Algorithm Overview
//!
//! The UMI deduplication problem: given a set of UMIs (DNA sequences), group together
//! all UMIs that are within a Hamming distance threshold (typically 1 mismatch).
//!
//! ## Complexity Analysis
//!
//! - **Pairwise comparisons:** O(n²) where n is the number of UMIs
//! - **Hamming distance per pair:** O(L) where L is UMI length (~12bp)
//! - **Union-Find operations:** O(α(n)) ≈ O(1) amortized, where α is inverse Ackermann
//! - **Total complexity:** O(n² × L × α(n)) ≈ O(n² × L)
//!
//! ## Performance Improvements vs Python
//!
//! 1. **Union-Find:** O(α(n)) find/union vs O(n) set operations in Python
//! 2. **Early exit:** Hamming distance stops at threshold+1 mismatches
//! 3. **Zero allocations:** Uses Vec for parent/rank arrays, no set copying
//! 4. **FxHashMap:** Fast integer hashing for grouping results
//!
//! Expected speedup: 3-6x over Python implementation

use rustc_hash::FxHashMap;
use rayon::prelude::*;

/// Union-Find data structure with path compression and union-by-rank
///
/// This data structure maintains a collection of disjoint sets and supports
/// two operations:
/// - `find(x)`: Find which set x belongs to (returns the root representative)
/// - `union(x, y)`: Merge the sets containing x and y
///
/// # Time Complexity
///
/// Both operations run in O(α(n)) amortized time, where α is the inverse
/// Ackermann function. For all practical purposes, α(n) ≤ 4, making these
/// operations effectively O(1).
///
/// # Space Complexity
///
/// O(n) where n is the number of elements
///
/// # Example
///
/// ```
/// use sheriff_rs::UnionFind;
///
/// let mut uf = UnionFind::new(5);
/// uf.union(0, 1);
/// uf.union(1, 2);
/// assert_eq!(uf.find(0), uf.find(2)); // 0 and 2 are in the same set
/// assert_ne!(uf.find(0), uf.find(3)); // 0 and 3 are in different sets
/// ```
#[derive(Debug, Clone)]
pub struct UnionFind {
    /// parent[i] = parent of element i
    /// If parent[i] == i, then i is a root (representative of its set)
    parent: Vec<usize>,
    
    /// rank[i] = upper bound on height of tree rooted at i
    /// Used for union-by-rank optimization
    rank: Vec<usize>,
}

impl UnionFind {
    /// Creates a new Union-Find structure with n elements
    ///
    /// Initially, each element is in its own set (parent[i] = i)
    ///
    /// # Arguments
    ///
    /// * `size` - Number of elements (0..size)
    ///
    /// # Example
    ///
    /// ```
    /// use sheriff_rs::UnionFind;
    /// let uf = UnionFind::new(10);
    /// ```
    pub fn new(size: usize) -> Self {
        Self {
            parent: (0..size).collect(),
            rank: vec![0; size],
        }
    }

    /// Find the root (representative) of the set containing x
    ///
    /// Uses path halving optimization: as we traverse to the root,
    /// we make each node point to its grandparent, compressing the path.
    ///
    /// # Path Compression
    ///
    /// Path compression ensures that subsequent find operations on the
    /// same path are faster. Over many operations, this gives amortized
    /// O(α(n)) time complexity.
    ///
    /// # Arguments
    ///
    /// * `x` - Element to find the root of
    ///
    /// # Returns
    ///
    /// The root (representative) of the set containing x
    ///
    /// # Example
    ///
    /// ```
    /// use sheriff_rs::UnionFind;
    ///
    /// let mut uf = UnionFind::new(5);
    /// uf.union(0, 1);
    /// uf.union(1, 2);
    /// // All of 0, 1, 2 should have the same root
    /// let root = uf.find(0);
    /// assert_eq!(root, uf.find(1));
    /// assert_eq!(root, uf.find(2));
    /// ```
    #[inline]
    pub fn find(&mut self, mut x: usize) -> usize {
        // Path halving: make each node point to its grandparent
        while self.parent[x] != x {
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }

    /// Union the sets containing x and y
    ///
    /// Uses union-by-rank optimization: always attach the shorter tree
    /// under the taller tree to keep the overall tree height small.
    ///
    /// # Arguments
    ///
    /// * `x` - Element in first set
    /// * `y` - Element in second set
    ///
    /// # Returns
    ///
    /// `true` if the sets were merged (they were different),
    /// `false` if they were already in the same set
    ///
    /// # Example
    ///
    /// ```
    /// use sheriff_rs::UnionFind;
    ///
    /// let mut uf = UnionFind::new(5);
    /// assert!(uf.union(0, 1));  // Sets merged, returns true
    /// assert!(!uf.union(0, 1)); // Already in same set, returns false
    /// ```
    #[inline]
    pub fn union(&mut self, x: usize, y: usize) -> bool {
        let root_x = self.find(x);
        let root_y = self.find(y);

        if root_x == root_y {
            return false; // Already in the same set
        }

        // Union by rank: attach smaller tree under larger tree
        if self.rank[root_x] < self.rank[root_y] {
            self.parent[root_x] = root_y;
        } else if self.rank[root_x] > self.rank[root_y] {
            self.parent[root_y] = root_x;
        } else {
            // Equal rank: choose one arbitrarily and increment its rank
            self.parent[root_y] = root_x;
            self.rank[root_x] += 1;
        }

        true
    }
}

/// Compute the Hamming distance between two byte sequences
///
/// The Hamming distance is the number of positions at which the
/// corresponding symbols differ.
///
/// # Time Complexity
///
/// O(L) where L = min(a.len(), b.len())
///
/// # Arguments
///
/// * `a` - First sequence
/// * `b` - Second sequence
///
/// # Returns
///
/// Number of positions where a and b differ
///
/// # Example
///
/// ```
/// use sheriff_rs::hamming_distance;
///
/// let dist = hamming_distance(b"ATCG", b"ATGG");
/// assert_eq!(dist, 1); // Differ at position 2 (C vs G)
/// ```
#[inline]
pub fn hamming_distance(a: &[u8], b: &[u8]) -> usize {
    a.iter()
        .zip(b.iter())
        .filter(|(x, y)| x != y)
        .count()
}

/// Check if two sequences are within a Hamming distance threshold
///
/// This function uses early termination: it stops counting as soon as
/// the number of mismatches exceeds the threshold, making it much faster
/// than computing the full Hamming distance.
///
/// # Time Complexity
///
/// - Best case: O(threshold + 1) when sequences differ significantly
/// - Worst case: O(L) where L is sequence length
///
/// # Arguments
///
/// * `a` - First sequence
/// * `b` - Second sequence
/// * `threshold` - Maximum allowed Hamming distance
///
/// # Returns
///
/// `true` if hamming_distance(a, b) <= threshold, `false` otherwise
///
/// # Example
///
/// ```
/// use sheriff_rs::within_hamming_threshold;
///
/// // Differ by 1 position
/// assert!(within_hamming_threshold(b"ATCG", b"ATGG", 1));
/// assert!(!within_hamming_threshold(b"ATCG", b"ATGG", 0));
///
/// // Differ by 3 positions
/// assert!(within_hamming_threshold(b"AAAA", b"TTTT", 3));
/// assert!(!within_hamming_threshold(b"AAAA", b"TTTT", 2));
/// ```
#[inline]
pub fn within_hamming_threshold(a: &[u8], b: &[u8], threshold: usize) -> bool {
    let mut mismatches = 0;
    
    for (x, y) in a.iter().zip(b.iter()) {
        if x != y {
            mismatches += 1;
            if mismatches > threshold {
                return false; // Early exit!
            }
        }
    }
    
    true
}

/// Deduplicate UMIs using Union-Find algorithm
///
/// Groups UMIs that are within a Hamming distance threshold into clusters.
/// Each cluster represents a set of UMIs that should be considered duplicates.
///
/// # Algorithm
///
/// 1. Initialize Union-Find with n elements (one per UMI)
/// 2. For each pair of UMIs:
///    - Compute Hamming distance with early exit
///    - If within threshold, union their sets
/// 3. Group UMIs by their root representative
///
/// # Time Complexity
///
/// O(n² × L × α(n)) where:
/// - n = number of UMIs
/// - L = UMI length
/// - α(n) = inverse Ackermann function ≈ O(1)
///
/// # Space Complexity
///
/// O(n) for Union-Find + O(k × n) for output groups where k is number of clusters
///
/// # Arguments
///
/// * `umis` - Slice of UMI byte sequences
/// * `threshold` - Maximum Hamming distance to consider UMIs as duplicates
///
/// # Returns
///
/// Vector of groups, where each group is a vector of UMI indices that are
/// within the threshold distance of each other
///
/// # Example
///
/// ```
/// use sheriff_rs::deduplicate_umis_unionfind;
///
/// let umis = vec![
///     b"ATCGATCG".as_slice(),
///     b"ATCGATCC".as_slice(), // 1 mismatch from umis[0]
///     b"GCGCGCGC".as_slice(), // Different cluster
/// ];
///
/// let groups = deduplicate_umis_unionfind(&umis, 1);
/// assert_eq!(groups.len(), 2); // Two clusters
///
/// // Find which cluster contains each UMI
/// let cluster_0 = groups.iter().find(|g| g.contains(&0)).unwrap();
/// let cluster_2 = groups.iter().find(|g| g.contains(&2)).unwrap();
///
/// assert!(cluster_0.contains(&1)); // UMIs 0 and 1 in same cluster
/// assert!(!cluster_2.contains(&0)); // UMI 2 in different cluster
/// ```
pub fn deduplicate_umis_unionfind(
    umis: &[&[u8]],
    threshold: usize,
) -> Vec<Vec<usize>> {
    let n = umis.len();
    
    if n == 0 {
        return vec![];
    }
    
    if n == 1 {
        return vec![vec![0]];
    }
    
    let mut uf = UnionFind::new(n);

    // Build union-find structure by comparing all pairs
    // Time: O(n² × L × α(n))
    for i in 0..n {
        for j in (i + 1)..n {
            if within_hamming_threshold(umis[i], umis[j], threshold) {
                uf.union(i, j);
            }
        }
    }

    // Group UMIs by their root representative
    // Time: O(n × α(n))
    let mut groups: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
    for i in 0..n {
        let root = uf.find(i);
        groups.entry(root).or_default().push(i);
    }

    // Convert to vector of groups
    groups.into_values().collect()
}

/// Deduplicate UMIs for multiple cells in parallel
///
/// This function processes multiple cells in parallel using Rayon, where each cell
/// has a set of UMIs that need to be deduplicated independently.
///
/// # Performance
///
/// This is the **key parallelization opportunity** in Sheriff:
/// - BAM file reading is inherently sequential (compressed format)
/// - Per-cell processing is embarrassingly parallel (cells are independent)
/// - Expected speedup: 6-8x on 8-core machines
///
/// # Time Complexity
///
/// For C cells with average U UMIs per cell:
/// - Sequential: O(C × U² × L)
/// - Parallel (P cores): O((C/P) × U² × L)
///
/// # Arguments
///
/// * `cells` - HashMap mapping cell barcodes to vectors of UMI sequences
/// * `threshold` - Maximum Hamming distance to consider UMIs as duplicates
///
/// # Returns
///
/// HashMap mapping cell barcodes to the number of unique UMI groups in that cell
///
/// # Example
///
/// ```
/// use sheriff_rs::deduplicate_cells_parallel;
/// use std::collections::HashMap;
///
/// let mut cells = HashMap::new();
/// cells.insert(
///     b"CELL001".to_vec(),
///     vec![b"ATCGATCG".to_vec(), b"ATCGATCC".to_vec()],
/// );
/// cells.insert(
///     b"CELL002".to_vec(),
///     vec![b"GCGCGCGC".to_vec(), b"GCGCGCGA".to_vec()],
/// );
///
/// let results = deduplicate_cells_parallel(cells, 1);
/// assert_eq!(results.len(), 2);
/// assert_eq!(results[&b"CELL001".to_vec()], 1); // Both UMIs in same group
/// assert_eq!(results[&b"CELL002".to_vec()], 1); // Both UMIs in same group
/// ```
pub fn deduplicate_cells_parallel(
    cells: std::collections::HashMap<Vec<u8>, Vec<Vec<u8>>>,
    threshold: usize,
) -> std::collections::HashMap<Vec<u8>, usize> {
    cells
        .par_iter()
        .map(|(cell_barcode, umis)| {
            // Convert Vec<Vec<u8>> to Vec<&[u8]> for deduplicate_umis_unionfind
            let umi_refs: Vec<&[u8]> = umis.iter().map(|u| u.as_slice()).collect();
            let groups = deduplicate_umis_unionfind(&umi_refs, threshold);
            (cell_barcode.clone(), groups.len())
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    // Test 1: Union-Find Basic Operations
    #[test]
    fn test_union_find_basic() {
        let mut uf = UnionFind::new(5);
        
        // Initially, each element is in its own set
        assert_eq!(uf.find(0), 0);
        assert_eq!(uf.find(1), 1);
        assert_eq!(uf.find(2), 2);
        
        // Union 0 and 1
        assert!(uf.union(0, 1));
        assert_eq!(uf.find(0), uf.find(1));
        
        // Union 1 and 2 (should merge with 0's set)
        assert!(uf.union(1, 2));
        assert_eq!(uf.find(0), uf.find(2));
        assert_eq!(uf.find(1), uf.find(2));
        
        // Try to union already-connected elements
        assert!(!uf.union(0, 2));
        
        // Element 3 should still be in its own set
        assert_ne!(uf.find(0), uf.find(3));
        
        // Union 3 with the main group
        assert!(uf.union(2, 3));
        assert_eq!(uf.find(0), uf.find(3));
        
        // Element 4 should still be in its own set
        assert_ne!(uf.find(0), uf.find(4));
    }

    // Test 2: Path Compression Correctness
    #[test]
    fn test_path_compression() {
        let mut uf = UnionFind::new(10);
        
        // Create a long chain: 0 -> 1 -> 2 -> 3 -> 4
        uf.union(0, 1);
        uf.union(1, 2);
        uf.union(2, 3);
        uf.union(3, 4);
        
        // All should have the same root
        let root = uf.find(0);
        assert_eq!(root, uf.find(1));
        assert_eq!(root, uf.find(2));
        assert_eq!(root, uf.find(3));
        assert_eq!(root, uf.find(4));
        
        // After find operations, path should be compressed
        // (hard to test directly, but we can verify correctness)
        assert_eq!(uf.find(0), uf.find(4));
    }

    // Test 3: Hamming Distance Edge Cases
    #[test]
    fn test_hamming_distance_edge_cases() {
        // Identical sequences
        assert_eq!(hamming_distance(b"ATCG", b"ATCG"), 0);
        
        // Completely different
        assert_eq!(hamming_distance(b"AAAA", b"TTTT"), 4);
        
        // Single mismatch
        assert_eq!(hamming_distance(b"ATCG", b"ATCC"), 1);
        
        // Multiple mismatches
        assert_eq!(hamming_distance(b"ATCGATCG", b"ATGGATCC"), 2);
        
        // Empty sequences
        assert_eq!(hamming_distance(b"", b""), 0);
        
        // Different lengths (compares minimum length)
        assert_eq!(hamming_distance(b"ATCG", b"AT"), 0);
        assert_eq!(hamming_distance(b"ATCG", b"TTCG"), 1);
    }

    // Test 4: Within Hamming Threshold
    #[test]
    fn test_within_hamming_threshold() {
        let umi1 = b"ATCGATCG";
        let umi2 = b"ATCGATCC"; // 1 mismatch
        let umi3 = b"ATCGTTCC"; // 2 mismatches from umi1
        let umi4 = b"GCGCGCGC"; // 8 mismatches from umi1 (all positions differ)

        // Test threshold 0 (exact match only)
        assert!(within_hamming_threshold(umi1, umi1, 0));
        assert!(!within_hamming_threshold(umi1, umi2, 0));

        // Test threshold 1
        assert!(within_hamming_threshold(umi1, umi1, 1));
        assert!(within_hamming_threshold(umi1, umi2, 1));
        assert!(!within_hamming_threshold(umi1, umi3, 1));
        assert!(!within_hamming_threshold(umi1, umi4, 1));

        // Test threshold 2
        assert!(within_hamming_threshold(umi1, umi2, 2));
        assert!(within_hamming_threshold(umi1, umi3, 2));
        assert!(!within_hamming_threshold(umi1, umi4, 2));

        // Test large threshold (umi1 and umi4 differ by 8, so threshold must be >= 8)
        assert!(!within_hamming_threshold(umi1, umi4, 7));
        assert!(within_hamming_threshold(umi1, umi4, 8));
    }

    // Test 5: Full Deduplication with Known Inputs
    #[test]
    fn test_deduplicate_known_inputs() {
        let umis = vec![
            b"ATCGATCG".as_slice(),
            b"ATCGATCC".as_slice(), // 1 mismatch from umis[0]
            b"GCGCGCGC".as_slice(), // Different cluster
            b"GCGCGCGC".as_slice(), // Exact duplicate of umis[2]
            b"GCGCGCGA".as_slice(), // 1 mismatch from umis[2]
        ];
        
        let groups = deduplicate_umis_unionfind(&umis, 1);
        
        // Should have 2 groups
        assert_eq!(groups.len(), 2);
        
        // Find which group contains each UMI
        let group_0 = groups.iter().find(|g| g.contains(&0)).unwrap();
        let group_2 = groups.iter().find(|g| g.contains(&2)).unwrap();
        
        // Group 0 should contain UMIs 0 and 1 (1 mismatch apart)
        assert_eq!(group_0.len(), 2);
        assert!(group_0.contains(&0));
        assert!(group_0.contains(&1));
        
        // Group 2 should contain UMIs 2, 3, and 4
        assert_eq!(group_2.len(), 3);
        assert!(group_2.contains(&2));
        assert!(group_2.contains(&3));
        assert!(group_2.contains(&4));
    }

    // Test 6: Empty and Single UMI Cases
    #[test]
    fn test_deduplicate_edge_cases() {
        // Empty input
        let groups = deduplicate_umis_unionfind(&[], 1);
        assert_eq!(groups.len(), 0);
        
        // Single UMI
        let umis = vec![b"ATCGATCG".as_slice()];
        let groups = deduplicate_umis_unionfind(&umis, 1);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0], vec![0]);
    }

    // Test 7: No Duplicates (All Different)
    #[test]
    fn test_no_duplicates() {
        let umis = vec![
            b"AAAAAAAA".as_slice(),
            b"CCCCCCCC".as_slice(),
            b"GGGGGGGG".as_slice(),
            b"TTTTTTTT".as_slice(),
        ];
        
        let groups = deduplicate_umis_unionfind(&umis, 1);
        
        // Each UMI should be in its own group
        assert_eq!(groups.len(), 4);
        for group in groups {
            assert_eq!(group.len(), 1);
        }
    }

    // Test 8: All Duplicates (All Same)
    #[test]
    fn test_all_duplicates() {
        let umis = vec![
            b"ATCGATCG".as_slice(),
            b"ATCGATCG".as_slice(),
            b"ATCGATCG".as_slice(),
            b"ATCGATCG".as_slice(),
        ];
        
        let groups = deduplicate_umis_unionfind(&umis, 1);
        
        // All should be in one group
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].len(), 4);
    }

    // Test 9: Transitive Clustering
    #[test]
    fn test_transitive_clustering() {
        // Test that if A ~= B and B ~= C, then A, B, C are in same cluster
        // even if A and C are not within threshold
        let umis = vec![
            b"AAAAAAAA".as_slice(), // 0
            b"AAAAATAA".as_slice(), // 1 - 1 mismatch from 0
            b"AAAAATTA".as_slice(), // 2 - 1 mismatch from 1, 2 mismatches from 0
        ];
        
        let groups = deduplicate_umis_unionfind(&umis, 1);
        
        // All three should be in the same cluster due to transitivity
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].len(), 3);
    }

    // Test 10: Compare Against Python-like Results
    #[test]
    fn test_python_equivalence() {
        // This test mimics the behavior of the Python implementation
        // in sheriff/helpers.py:225-253
        
        let umis = vec![
            b"ATCGATCG".as_slice(),
            b"ATCGATCC".as_slice(),
            b"GCGCGCGC".as_slice(),
        ];
        
        // Python would create these groups with threshold=1:
        // Group 1: {ATCGATCG, ATCGATCC} (1 mismatch)
        // Group 2: {GCGCGCGC}
        
        let groups = deduplicate_umis_unionfind(&umis, 1);
        assert_eq!(groups.len(), 2);
        
        // Verify the grouping matches Python output
        let mut group_sizes: Vec<usize> = groups.iter().map(|g| g.len()).collect();
        group_sizes.sort();
        assert_eq!(group_sizes, vec![1, 2]);
    }
}
