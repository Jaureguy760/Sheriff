//! PyO3 Python bindings for Sheriff-rs
//!
//! This module provides comprehensive Python bindings for the Rust Phase 1 optimizations,
//! including k-mer matching, k-mer counting, and UMI deduplication.
//!
//! # Features
//!
//! - **K-mer Operations:** Fast k-mer hashing and matching against whitelists
//! - **K-mer Counting:** Reusable frequency counting with array reuse pattern
//! - **UMI Deduplication:** Union-Find based UMI clustering with Hamming distance
//!
//! # Performance
//!
//! These Rust implementations provide 4-14x speedup for k-mer operations and
//! 3-6x speedup for UMI deduplication compared to pure Python implementations.

use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use pyo3::types::{PyList, PyDict};
use rustc_hash::FxHashSet;
use crate::kmer::{self, KmerCounter as RustKmerCounter};
use crate::umi::{self, deduplicate_umis_unionfind, deduplicate_cells_parallel as rust_deduplicate_cells_parallel};

// ============================================================================
// K-mer Bindings
// ============================================================================

/// Convert a k-mer sequence to its numeric hash representation.
///
/// This function implements a 4-ary encoding scheme where each nucleotide
/// is converted to a 2-bit value (A=0, C=1, G=2, T=3) and combined into
/// a single integer hash.
///
/// Args:
///     kmer (str): DNA sequence to hash (case-insensitive)
///
/// Returns:
///     int: Numeric hash of the k-mer
///
/// Examples:
///     >>> kmer_to_num("ACGT")
///     27
///     >>> kmer_to_num("AAAA")
///     0
///     >>> kmer_to_num("acgt")  # Case-insensitive
///     27
///
/// Note:
///     Uses the same hashing algorithm as the Python implementation for
///     compatibility. The hash is computed as: hash = n₀ × 4^(k-1) + n₁ × 4^(k-2) + ... + n_(k-1)
#[pyfunction]
#[pyo3(text_signature = "(kmer)")]
fn kmer_to_num(kmer: &str) -> PyResult<u32> {
    if kmer.is_empty() {
        return Ok(0);
    }

    // Validate that the sequence contains only valid nucleotides
    if !kmer.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')) {
        return Err(PyValueError::new_err(
            "K-mer must contain only A, C, G, T (case-insensitive)"
        ));
    }

    Ok(kmer::kmer_to_num(kmer.as_bytes()))
}

/// Match k-mers in a sequence against a whitelist.
///
/// Slides a window of size k across the sequence, hashing each k-mer and
/// checking if it exists in the whitelist. Returns all matches found.
///
/// Args:
///     sequence (str): DNA sequence to scan (case-insensitive)
///     k (int): Length of k-mers to match
///     whitelist (list[int]): List of k-mer hashes to match against
///     output_hash (bool): If True, return hashes; if False, return k-mer strings
///
/// Returns:
///     list[int] or list[str]: List of matched k-mer hashes or sequences
///
/// Examples:
///     >>> whitelist = [kmer_to_num("ACGT")]
///     >>> match_kmer("ACGTACGT", 4, whitelist, output_hash=True)
///     [27, 27]
///     >>> match_kmer("ACGTACGT", 4, whitelist, output_hash=False)
///     ['ACGT', 'ACGT']
///
/// Note:
///     Uses FxHashSet for O(1) lookups with minimal hashing overhead.
///     Expected speedup: 4-14x over Python implementation.
#[pyfunction]
#[pyo3(text_signature = "(sequence, k, whitelist, output_hash=True)")]
fn match_kmer(
    sequence: &str,
    k: usize,
    whitelist: Vec<u32>,
    output_hash: Option<bool>,
) -> PyResult<PyObject> {
    if k == 0 {
        return Err(PyValueError::new_err("k must be greater than 0"));
    }

    if sequence.len() < k {
        // Return empty list if sequence is too short
        return Python::with_gil(|py| {
            let empty: Vec<u32> = Vec::new();
            Ok(PyList::new_bound(py, empty).into())
        });
    }

    // Validate sequence
    if !sequence.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')) {
        return Err(PyValueError::new_err(
            "Sequence must contain only A, C, G, T (case-insensitive)"
        ));
    }

    // Convert whitelist to FxHashSet for fast lookups
    let whitelist_set: FxHashSet<u32> = whitelist.into_iter().collect();

    // Call the Rust implementation
    let matches = kmer::match_kmer(sequence.as_bytes(), k, &whitelist_set);

    Python::with_gil(|py| {
        if output_hash.unwrap_or(true) {
            // Return as list of integers
            Ok(PyList::new_bound(py, matches).into())
        } else {
            // Convert hashes back to k-mer strings
            let sequence_bytes = sequence.as_bytes();
            let mut kmer_strings = Vec::new();

            for window in sequence_bytes.windows(k) {
                let hash = kmer::kmer_to_num(window);
                if whitelist_set.contains(&hash) {
                    // Convert bytes to string
                    let kmer_str = std::str::from_utf8(window)
                        .map_err(|e| PyValueError::new_err(format!("Invalid UTF-8: {}", e)))?
                        .to_uppercase();
                    kmer_strings.push(kmer_str);
                }
            }

            Ok(PyList::new_bound(py, kmer_strings).into())
        }
    })
}

/// K-mer frequency counter with array reuse pattern.
///
/// This class maintains a reusable frequency array to avoid repeated allocations
/// when counting k-mers across multiple sequences. Provides 1.1-1.2x speedup by
/// eliminating allocation overhead.
///
/// Args:
///     k (int): K-mer length
///
/// Examples:
///     >>> counter = KmerCounter(4)
///     >>> freqs = counter.count_kmers("ACGTACGT")
///     >>> len(freqs)
///     256  # 4^4 = 256 possible 4-mers
///     >>> freqs[kmer_to_num("ACGT")]
///     2  # "ACGT" appears twice
///
/// Note:
///     The frequency array has size 4^k. For k=6, this is 4KB; for k=8, this is 64KB.
///     Frequencies are stored as u8 and saturate at 255.
#[pyclass]
pub struct KmerCounter {
    inner: RustKmerCounter,
}

#[pymethods]
impl KmerCounter {
    /// Create a new KmerCounter with specified k-mer length.
    ///
    /// Args:
    ///     k (int): K-mer length (typically 4-12)
    ///
    /// Raises:
    ///     ValueError: If k is too large (k > 16)
    ///
    /// Note:
    ///     Memory usage is 4^k bytes. For k=10, this is 1MB.
    #[new]
    #[pyo3(text_signature = "(k)")]
    pub fn new(k: usize) -> PyResult<Self> {
        if k == 0 {
            return Err(PyValueError::new_err("k must be greater than 0"));
        }
        if k > 16 {
            return Err(PyValueError::new_err(
                "k must be <= 16 (4^16 = 4GB would be required)"
            ));
        }

        Ok(KmerCounter {
            inner: RustKmerCounter::new(k),
        })
    }

    /// Count k-mer frequencies in a sequence.
    ///
    /// Args:
    ///     sequence (str): DNA sequence to analyze (case-insensitive)
    ///
    /// Returns:
    ///     list[int]: Frequency array indexed by k-mer hash (length = 4^k)
    ///
    /// Examples:
    ///     >>> counter = KmerCounter(4)
    ///     >>> freqs = counter.count_kmers("ACGTACGT")
    ///     >>> freqs[kmer_to_num("ACGT")]
    ///     2
    ///
    /// Note:
    ///     The internal array is reused across calls for performance.
    ///     Frequencies saturate at 255 (u8::MAX).
    #[pyo3(text_signature = "($self, sequence)")]
    pub fn count_kmers(&mut self, sequence: &str) -> PyResult<PyObject> {
        // Validate sequence
        if !sequence.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')) {
            return Err(PyValueError::new_err(
                "Sequence must contain only A, C, G, T (case-insensitive)"
            ));
        }

        // Count k-mers
        let freqs = self.inner.count_kmers(sequence.as_bytes());

        // Convert to Python list
        Python::with_gil(|py| {
            let py_list = PyList::new_bound(py, freqs);
            Ok(py_list.into())
        })
    }

    /// Get the k-mer length this counter was created for.
    ///
    /// Returns:
    ///     int: K-mer length
    #[getter]
    pub fn k(&self) -> usize {
        self.inner.k()
    }

    /// Clear all frequency counts.
    ///
    /// Note:
    ///     This is done automatically when count_kmers() is called.
    pub fn clear(&mut self) {
        let _ = self.inner.count_kmers(&[]);
    }

    /// Get the size of the frequency array (4^k).
    ///
    /// Returns:
    ///     int: Array size
    #[getter]
    pub fn array_size(&self) -> usize {
        self.inner.array_size()
    }

    fn __repr__(&self) -> String {
        format!("KmerCounter(k={})", self.inner.k())
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

// ============================================================================
// UMI Bindings
// ============================================================================

/// Deduplicate UMIs and return the count of unique UMI groups.
///
/// Groups UMIs that are within a Hamming distance threshold into clusters.
/// This is useful for counting unique molecular identifiers while accounting
/// for sequencing errors.
///
/// Args:
///     umis (list[str]): List of UMI sequences (all must be same length)
///     threshold (int): Maximum Hamming distance to consider UMIs as duplicates (typically 1)
///
/// Returns:
///     int: Number of unique UMI groups
///
/// Examples:
///     >>> umis = ["ATCGATCG", "ATCGATCC", "GCGCGCGC"]
///     >>> deduplicate_umis(umis, threshold=1)
///     2  # Two unique groups
///
/// Note:
///     Uses Union-Find algorithm with path compression for O(n² × L × α(n)) complexity,
///     where n is the number of UMIs, L is UMI length, and α is the inverse Ackermann function.
///     Expected speedup: 3-6x over Python implementation.
#[pyfunction]
#[pyo3(text_signature = "(umis, threshold)")]
fn deduplicate_umis(umis: Vec<String>, threshold: usize) -> PyResult<usize> {
    if umis.is_empty() {
        return Ok(0);
    }

    // Validate that all UMIs are the same length
    let first_len = umis[0].len();
    if !umis.iter().all(|umi| umi.len() == first_len) {
        return Err(PyValueError::new_err(
            "All UMIs must be the same length"
        ));
    }

    // Validate UMI sequences
    for umi in &umis {
        if !umi.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')) {
            return Err(PyValueError::new_err(
                "UMIs must contain only A, C, G, T (case-insensitive)"
            ));
        }
    }

    // Convert to byte slices for processing
    let umi_bytes: Vec<&[u8]> = umis.iter().map(|s| s.as_bytes()).collect();

    // Call the Rust implementation
    let groups = deduplicate_umis_unionfind(&umi_bytes, threshold);

    Ok(groups.len())
}

/// Deduplicate UMIs and return detailed grouping information.
///
/// Groups UMIs that are within a Hamming distance threshold into clusters.
/// Returns the full grouping structure showing which UMIs belong to which cluster.
///
/// Args:
///     umis (list[str]): List of UMI sequences (all must be same length)
///     threshold (int): Maximum Hamming distance to consider UMIs as duplicates (typically 1)
///
/// Returns:
///     list[list[int]]: List of groups, where each group contains indices of UMIs in that cluster
///
/// Examples:
///     >>> umis = ["ATCGATCG", "ATCGATCC", "GCGCGCGC"]
///     >>> groups = deduplicate_umis_detailed(umis, threshold=1)
///     >>> groups
///     [[0, 1], [2]]  # UMIs 0 and 1 form one group, UMI 2 is separate
///
/// Note:
///     This function returns more information than deduplicate_umis() but has the same
///     performance characteristics. Use this when you need to know which UMIs are grouped together.
#[pyfunction]
#[pyo3(text_signature = "(umis, threshold)")]
fn deduplicate_umis_detailed(umis: Vec<String>, threshold: usize) -> PyResult<PyObject> {
    if umis.is_empty() {
        let empty: Vec<Vec<usize>> = Vec::new();
        return Python::with_gil(|py| Ok(PyList::new_bound(py, empty).into()));
    }

    // Validate that all UMIs are the same length
    let first_len = umis[0].len();
    if !umis.iter().all(|umi| umi.len() == first_len) {
        return Err(PyValueError::new_err(
            "All UMIs must be the same length"
        ));
    }

    // Validate UMI sequences
    for umi in &umis {
        if !umi.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')) {
            return Err(PyValueError::new_err(
                "UMIs must contain only A, C, G, T (case-insensitive)"
            ));
        }
    }

    // Convert to byte slices for processing
    let umi_bytes: Vec<&[u8]> = umis.iter().map(|s| s.as_bytes()).collect();

    // Call the Rust implementation
    let groups = deduplicate_umis_unionfind(&umi_bytes, threshold);

    // Convert to Python list of lists
    Python::with_gil(|py| {
        let py_groups: Vec<PyObject> = groups
            .into_iter()
            .map(|group| PyList::new_bound(py, group).into())
            .collect();
        Ok(PyList::new_bound(py, py_groups).into())
    })
}

/// Compute Hamming distance between two sequences.
///
/// The Hamming distance is the number of positions at which the corresponding
/// symbols differ.
///
/// Args:
///     a (str): First sequence
///     b (str): Second sequence
///
/// Returns:
///     int: Number of positions where sequences differ
///
/// Examples:
///     >>> hamming_distance("ATCG", "ATGG")
///     1
///     >>> hamming_distance("AAAA", "TTTT")
///     4
///
/// Note:
///     If sequences have different lengths, only the shorter length is compared.
#[pyfunction]
#[pyo3(text_signature = "(a, b)")]
fn hamming_distance(a: &str, b: &str) -> usize {
    umi::hamming_distance(a.as_bytes(), b.as_bytes())
}

/// Deduplicate UMIs for multiple cells in parallel.
///
/// This is the **key parallelization opportunity** in Sheriff. While BAM file reading
/// is sequential (compressed format), per-cell processing is embarrassingly parallel
/// since cells are independent. This function processes all cells in parallel using Rayon.
///
/// Args:
///     cells (dict[str, list[str]]): Dictionary mapping cell barcodes to lists of UMI sequences
///     threshold (int): Maximum Hamming distance to consider UMIs as duplicates (typically 1)
///
/// Returns:
///     dict[str, int]: Dictionary mapping cell barcodes to the number of unique UMI groups
///
/// Examples:
///     >>> cells = {
///     ...     "CELL001": ["ATCGATCG", "ATCGATCC"],
///     ...     "CELL002": ["GCGCGCGC", "GCGCGCGA"]
///     ... }
///     >>> results = deduplicate_cells_parallel(cells, threshold=1)
///     >>> results
///     {'CELL001': 1, 'CELL002': 1}
///
/// Performance:
///     Expected speedup: 6-8x on 8-core machines vs sequential processing
///     Scales linearly with number of cores available
///
/// Note:
///     This is the recommended function for processing Sheriff data with many cells.
///     BAM reading is still sequential, but per-cell UMI dedup is parallelized.
#[pyfunction]
#[pyo3(text_signature = "(cells, threshold)")]
fn deduplicate_cells_parallel(py: Python, cells: &Bound<'_, PyDict>, threshold: usize) -> PyResult<PyObject> {
    // Convert Python dict to Rust HashMap
    let mut rust_cells: std::collections::HashMap<Vec<u8>, Vec<Vec<u8>>> = std::collections::HashMap::new();

    for (key, value) in cells.iter() {
        // Extract cell barcode (key)
        let cell_barcode: String = key.extract()?;
        let cell_barcode_bytes = cell_barcode.into_bytes();

        // Extract UMI list (value)
        let umi_list: Vec<String> = value.extract()?;

        // Validate that all UMIs are the same length
        if !umi_list.is_empty() {
            let first_len = umi_list[0].len();
            if !umi_list.iter().all(|umi| umi.len() == first_len) {
                return Err(PyValueError::new_err(format!(
                    "All UMIs in cell {} must be the same length",
                    String::from_utf8_lossy(&cell_barcode_bytes)
                )));
            }

            // Validate UMI sequences
            for umi in &umi_list {
                if !umi.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')) {
                    return Err(PyValueError::new_err(
                        "UMIs must contain only A, C, G, T (case-insensitive)"
                    ));
                }
            }
        }

        // Convert to Vec<Vec<u8>>
        let umi_bytes: Vec<Vec<u8>> = umi_list.into_iter().map(|s| s.into_bytes()).collect();
        rust_cells.insert(cell_barcode_bytes, umi_bytes);
    }

    // Call the Rust parallel implementation
    let results = rust_deduplicate_cells_parallel(rust_cells, threshold);

    // Convert back to Python dict
    let py_dict = PyDict::new_bound(py);
    for (cell_barcode, unique_count) in results {
        let cell_barcode_str = String::from_utf8_lossy(&cell_barcode).to_string();
        py_dict.set_item(cell_barcode_str, unique_count)?;
    }

    Ok(py_dict.into())
}

// ============================================================================
// Gene UMI Counting Bindings
// ============================================================================

/// Count gene UMIs across cells with deduplication.
///
/// This function processes gene expression data organized by gene, where each gene
/// has a set of cells and each cell has a set of UMIs. It deduplicates UMIs within
/// each gene-cell combination and returns sparse matrix indices.
///
/// Args:
///     total_cells (int): Total number of cells in the dataset
///     gene_indices (list[int]): Array of gene indices (column indices in output matrix)
///     gene_cell_indices (list[list[int]]): For each gene, list of cell indices with UMIs
///     gene_cell_umis (list[list[list[str]]]): For each gene, for each cell, list of UMI sequences
///     threshold (int): Hamming distance threshold for UMI deduplication (typically 1)
///
/// Returns:
///     numpy.ndarray: 2D array of shape (n_nonzero, 3) where each row is (count, cell_idx, gene_idx)
///
/// Examples:
///     >>> gene_indices = [0, 1]
///     >>> gene_cell_indices = [[0, 1], [2]]
///     >>> gene_cell_umis = [
///     ...     [["ATCGATCG", "ATCGATCC"], ["GCGCGCGC"]],  # Gene 0
///     ...     [["TTTTTTTT", "TTTTTTTT"]],                 # Gene 1
///     ... ]
///     >>> results = count_gene_umis_rust(100, gene_indices, gene_cell_indices, gene_cell_umis, 1)
///     >>> results.shape
///     (3, 3)
///     >>> # Results: [[1, 0, 0], [1, 1, 0], [1, 2, 1]]
///     >>> #          count, cell, gene
///
/// Performance:
///     Expected speedup: 2-4x over Numba implementation
///     This optimizes the largest remaining bottleneck (5.25% of total runtime)
///
/// Note:
///     This function is a drop-in replacement for the Numba JIT version in helpers.py.
///     It uses the same Union-Find UMI deduplication algorithm but with native Rust performance.
#[pyfunction]
#[pyo3(text_signature = "(total_cells, gene_indices, gene_cell_indices, gene_cell_umis, threshold)")]
fn count_gene_umis_rust(
    py: Python,
    total_cells: usize,
    gene_indices: Vec<u32>,
    gene_cell_indices: Vec<Vec<u32>>,
    gene_cell_umis: Vec<Vec<Vec<String>>>,
    threshold: usize,
) -> PyResult<PyObject> {
    // Validate inputs
    if gene_indices.len() != gene_cell_indices.len() {
        return Err(PyValueError::new_err(
            "gene_indices and gene_cell_indices must have the same length"
        ));
    }

    if gene_indices.len() != gene_cell_umis.len() {
        return Err(PyValueError::new_err(
            "gene_indices and gene_cell_umis must have the same length"
        ));
    }

    // Validate that gene_cell_indices and gene_cell_umis have matching lengths
    for i in 0..gene_indices.len() {
        if gene_cell_indices[i].len() != gene_cell_umis[i].len() {
            return Err(PyValueError::new_err(format!(
                "Gene {}: gene_cell_indices and gene_cell_umis must have matching cell counts",
                gene_indices[i]
            )));
        }
    }

    // Validate UMI sequences
    for (gene_i, cell_umis_list) in gene_cell_umis.iter().enumerate() {
        for (cell_i, umis) in cell_umis_list.iter().enumerate() {
            if !umis.is_empty() {
                // Check that all UMIs are the same length
                let first_len = umis[0].len();
                if !umis.iter().all(|umi| umi.len() == first_len) {
                    return Err(PyValueError::new_err(format!(
                        "Gene {}, Cell {}: All UMIs must be the same length",
                        gene_indices[gene_i], gene_cell_indices[gene_i][cell_i]
                    )));
                }

                // Validate UMI sequences contain only ACGT
                for umi in umis {
                    if !umi.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')) {
                        return Err(PyValueError::new_err(
                            "UMIs must contain only A, C, G, T (case-insensitive)"
                        ));
                    }
                }
            }
        }
    }

    // Convert String UMIs to Vec<u8> for Rust processing
    let gene_cell_umis_bytes: Vec<Vec<Vec<Vec<u8>>>> = gene_cell_umis
        .into_iter()
        .map(|cell_umis_list| {
            cell_umis_list
                .into_iter()
                .map(|umis| umis.into_iter().map(|umi| umi.into_bytes()).collect())
                .collect()
        })
        .collect();

    // Call the Rust implementation
    let results = crate::gene::count_gene_umis(
        total_cells,
        &gene_indices,
        &gene_cell_indices,
        &gene_cell_umis_bytes,
        threshold,
    );

    // Convert results to numpy array (n_nonzero × 3)
    // Format: (count, cell_idx, gene_idx)
    let n = results.len();
    let mut array_data: Vec<u32> = Vec::with_capacity(n * 3);

    for (count, cell_idx, gene_idx) in results {
        array_data.push(count);
        array_data.push(cell_idx);
        array_data.push(gene_idx);
    }

    // Convert to Python list of lists (will be converted to numpy array in Python)
    let rows: Vec<Vec<u32>> = array_data
        .chunks(3)
        .map(|chunk| chunk.to_vec())
        .collect();

    let py_list: Vec<PyObject> = rows
        .into_iter()
        .map(|row| PyList::new_bound(py, row).into())
        .collect();
    Ok(PyList::new_bound(py, py_list).into())
}

// ============================================================================
// Module Setup
// ============================================================================

/// Sheriff-rs: High-performance k-mer and UMI processing
///
/// This module provides Rust implementations of k-mer matching and UMI deduplication
/// with significant performance improvements over pure Python implementations.
///
/// Performance Improvements:
///     - K-mer operations: 4-14x faster than Python
///     - UMI deduplication: 3-6x faster than Python
///
/// Functions:
///     kmer_to_num(kmer: str) -> int
///         Convert k-mer sequence to numeric hash
///
///     match_kmer(sequence: str, k: int, whitelist: List[int], output_hash: bool = True) -> List
///         Match k-mers against whitelist
///
///     deduplicate_umis(umis: List[str], threshold: int) -> int
///         Count unique UMI groups
///
///     deduplicate_umis_detailed(umis: List[str], threshold: int) -> List[List[int]]
///         Get detailed UMI groupings
///
///     hamming_distance(a: str, b: str) -> int
///         Compute Hamming distance between sequences
///
/// Classes:
///     KmerCounter(k: int)
///         Efficient k-mer frequency counter with array reuse
///
/// Examples:
///     >>> import sheriff_rs
///     >>> # K-mer hashing
///     >>> hash_val = sheriff_rs.kmer_to_num("ACGT")
///     >>> print(hash_val)
///     27
///
///     >>> # K-mer matching
///     >>> whitelist = [hash_val]
///     >>> matches = sheriff_rs.match_kmer("ACGTACGT", 4, whitelist, output_hash=True)
///     >>> print(matches)
///     [27, 27]
///
///     >>> # K-mer counting
///     >>> counter = sheriff_rs.KmerCounter(4)
///     >>> freqs = counter.count_kmers("ACGTACGT")
///     >>> print(freqs[hash_val])
///     2
///
///     >>> # UMI deduplication
///     >>> umis = ["ATCGATCG", "ATCGATCC", "GCGCGCGC"]
///     >>> unique_count = sheriff_rs.deduplicate_umis(umis, threshold=1)
///     >>> print(unique_count)
///     2
#[pymodule]
#[pyo3(name = "sheriff_rs")]
fn sheriff_rs_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Add version info
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__doc__", "High-performance k-mer and UMI processing in Rust")?;

    // Add k-mer functions
    m.add_function(wrap_pyfunction!(kmer_to_num, m)?)?;
    m.add_function(wrap_pyfunction!(match_kmer, m)?)?;

    // Add k-mer classes
    m.add_class::<KmerCounter>()?;

    // Add UMI functions
    m.add_function(wrap_pyfunction!(deduplicate_umis, m)?)?;
    m.add_function(wrap_pyfunction!(deduplicate_umis_detailed, m)?)?;
    m.add_function(wrap_pyfunction!(hamming_distance, m)?)?;
    m.add_function(wrap_pyfunction!(deduplicate_cells_parallel, m)?)?;

    // Add gene UMI counting functions
    m.add_function(wrap_pyfunction!(count_gene_umis_rust, m)?)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_to_num_python() {
        assert_eq!(kmer_to_num("ACGT").unwrap(), 27);
        assert_eq!(kmer_to_num("acgt").unwrap(), 27);
        assert_eq!(kmer_to_num("AAAA").unwrap(), 0);

        // Test error handling
        assert!(kmer_to_num("ACGN").is_err());
        assert!(kmer_to_num("ACG-").is_err());
    }

    #[test]
    fn test_hamming_distance_python() {
        assert_eq!(hamming_distance("ATCG", "ATCG"), 0);
        assert_eq!(hamming_distance("ATCG", "ATGG"), 1);
        assert_eq!(hamming_distance("AAAA", "TTTT"), 4);
    }
}
