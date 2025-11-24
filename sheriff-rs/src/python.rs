use pyo3::prelude::*;
use pyo3::types::PyDict;
use crate::bam_filter::{filter_bam_by_barcodes, filter_bam_by_barcodes_parallel, filter_bam_by_barcodes_chromosome_parallel, load_whitelist};
use crate::umi::{deduplicate_umis_rust, cell_umi_counts_rust, cell_umi_counts_rust_parallel};
use crate::edit_clustering::{Edit, get_longest_edits};
use crate::gene_counts::gene_counts_per_cell;
use std::collections::HashSet;

/// Filter BAM file by cell barcode whitelist (Python wrapper - file-based)
///
/// # Arguments
/// * `input_path` - Path to input BAM file
/// * `output_path` - Path to output BAM file
/// * `whitelist_path` - Path to whitelist text file
///
/// # Returns
/// Dict with keys: reads_processed, reads_kept, reads_rejected, duration_seconds
#[pyfunction]
fn filter_bam_rust(
    input_path: String,
    output_path: String,
    whitelist_path: String,
) -> PyResult<PyObject> {
    use std::time::Instant;

    let start = Instant::now();

    // Load whitelist
    let whitelist = load_whitelist(&whitelist_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    // Filter BAM
    let result = filter_bam_by_barcodes(&input_path, &output_path, &whitelist)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    let duration = start.elapsed();

    // Return dict
    Python::with_gil(|py| {
        let dict = PyDict::new(py);
        dict.set_item("reads_processed", result.reads_processed)?;
        dict.set_item("reads_kept", result.reads_kept)?;
        dict.set_item("reads_rejected", result.reads_rejected)?;
        dict.set_item("duration_seconds", duration.as_secs_f64())?;
        Ok(dict.into())
    })
}

/// Filter BAM file by cell barcode whitelist (Python wrapper - direct barcodes)
///
/// This version accepts barcodes directly as a list, avoiding temporary file overhead.
///
/// # Arguments
/// * `input_path` - Path to input BAM file
/// * `output_path` - Path to output BAM file
/// * `barcodes` - List of cell barcodes
///
/// # Returns
/// Dict with keys: reads_processed, reads_kept, reads_rejected, duration_seconds
#[pyfunction]
fn filter_bam_by_barcodes_rust(
    input_path: String,
    output_path: String,
    barcodes: Vec<String>,
) -> PyResult<PyObject> {
    use std::time::Instant;

    let start = Instant::now();

    // Convert to HashSet
    let whitelist: HashSet<String> = barcodes.into_iter().collect();

    // Filter BAM
    let result = filter_bam_by_barcodes(&input_path, &output_path, &whitelist)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    let duration = start.elapsed();

    // Return dict
    Python::with_gil(|py| {
        let dict = PyDict::new(py);
        dict.set_item("reads_processed", result.reads_processed)?;
        dict.set_item("reads_kept", result.reads_kept)?;
        dict.set_item("reads_rejected", result.reads_rejected)?;
        dict.set_item("duration_seconds", duration.as_secs_f64())?;
        Ok(dict.into())
    })
}

/// Count k-mers in DNA sequence (Python wrapper)
///
/// # Arguments
/// * `sequence` - DNA sequence string
/// * `k` - K-mer length
///
/// # Returns
/// List of k-mer counts (length 4^k)
#[pyfunction]
fn count_kmers_rust(sequence: String, k: usize) -> Vec<u8> {
    crate::kmer::count_kmers(sequence.as_bytes(), k)
}

/// Match k-mers in DNA sequence against whitelist (Python wrapper)
///
/// Finds all k-mers in sequence that appear in the whitelist.
/// Handles 'N' bases by skipping k-mers containing them.
///
/// # Arguments
/// * `sequence` - DNA sequence string
/// * `k` - K-mer length
/// * `whitelist` - Optional list of k-mer hashes to match against (None = all k-mers)
/// * `output_hash` - Return hashes (true) or k-mer strings (false)
///
/// # Returns
/// List of matching k-mer hashes (if output_hash=True) or strings (if output_hash=False)
#[pyfunction]
#[pyo3(signature = (sequence, k, whitelist=None, output_hash=true))]
fn match_kmer_rust(
    sequence: String,
    k: usize,
    whitelist: Option<Vec<usize>>,
    output_hash: bool,
) -> PyResult<Vec<PyObject>> {
    use std::collections::HashSet;

    // Convert whitelist to FxHashSet (PHASE 1 OPTIMIZATION: 2-3x faster than std::HashSet)
    let whitelist_set = whitelist.map(|v| v.into_iter().collect::<HashSet<usize>>());

    // Call Rust implementation
    let matches = crate::kmer::match_kmer(
        &sequence,
        k,
        whitelist_set.as_ref(),
        output_hash,
    );

    // Convert results to Python objects
    Python::with_gil(|py| {
        let result: Vec<PyObject> = matches.iter().map(|m| {
            match m {
                crate::kmer::KmerMatch::Hash(h) => h.to_object(py),
                crate::kmer::KmerMatch::String(s) => s.to_object(py),
            }
        }).collect();
        Ok(result)
    })
}

/// Filter BAM file by cell barcode whitelist (Python wrapper - parallel version)
///
/// This version uses rayon's par_bridge for parallel filtering.
/// Expected speedup: 2-5x on medium/large files.
///
/// # Arguments
/// * `input_path` - Path to input BAM file
/// * `output_path` - Path to output BAM file
/// * `barcodes` - List of cell barcodes
///
/// # Returns
/// Dict with keys: reads_processed, reads_kept, reads_rejected, duration_seconds
#[pyfunction]
fn filter_bam_by_barcodes_rust_parallel(
    input_path: String,
    output_path: String,
    barcodes: Vec<String>,
) -> PyResult<PyObject> {
    use std::time::Instant;

    let start = Instant::now();

    // Convert to HashSet
    let whitelist: HashSet<String> = barcodes.into_iter().collect();

    // Filter BAM with parallel processing
    let result = filter_bam_by_barcodes_parallel(&input_path, &output_path, &whitelist)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    let duration = start.elapsed();

    // Return dict
    Python::with_gil(|py| {
        let dict = PyDict::new(py);
        dict.set_item("reads_processed", result.reads_processed)?;
        dict.set_item("reads_kept", result.reads_kept)?;
        dict.set_item("reads_rejected", result.reads_rejected)?;
        dict.set_item("duration_seconds", duration.as_secs_f64())?;
        Ok(dict.into())
    })
}

/// Filter BAM file by cell barcode whitelist (Python wrapper - chromosome-based parallel)
///
/// This version uses chromosome-based parallelism for maximum speedup.
/// Expected speedup: 10-50x on large files.
/// Requires indexed BAM file (*.bam.bai).
///
/// # Arguments
/// * `input_path` - Path to input BAM file (must be indexed)
/// * `output_path` - Path to output BAM file
/// * `barcodes` - List of cell barcodes
/// * `num_threads` - Number of parallel threads (None = auto)
///
/// # Returns
/// Dict with keys: reads_processed, reads_kept, reads_rejected, duration_seconds
#[pyfunction]
#[pyo3(signature = (input_path, output_path, barcodes, num_threads=None))]
fn filter_bam_by_barcodes_rust_chromosome(
    input_path: String,
    output_path: String,
    barcodes: Vec<String>,
    num_threads: Option<usize>,
) -> PyResult<PyObject> {
    use std::time::Instant;

    let start = Instant::now();

    // Convert to HashSet
    let whitelist: HashSet<String> = barcodes.into_iter().collect();

    // Filter BAM with chromosome-based parallel processing
    let result = filter_bam_by_barcodes_chromosome_parallel(
        &input_path,
        &output_path,
        &whitelist,
        num_threads,
    )
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    let duration = start.elapsed();

    // Return dict
    Python::with_gil(|py| {
        let dict = PyDict::new(py);
        dict.set_item("reads_processed", result.reads_processed)?;
        dict.set_item("reads_kept", result.reads_kept)?;
        dict.set_item("reads_rejected", result.reads_rejected)?;
        dict.set_item("duration_seconds", duration.as_secs_f64())?;
        Ok(dict.into())
    })
}

/// Deduplicate UMIs using graph-based connected components (Python wrapper)
///
/// Fast Rust implementation of UMI deduplication using petgraph.
/// Identifies unique molecules by clustering UMIs within Hamming distance 1.
///
/// # Arguments
/// * `umis` - List of UMI sequences (all must be same length)
///
/// # Returns
/// Number of unique UMI groups (connected components)
///
/// # Example (Python)
/// ```python
/// from sheriff_rs import deduplicate_umis_rust
/// umis = ["ATCG", "ATCG", "TTTT"]
/// unique_count = deduplicate_umis_rust(umis)
/// print(unique_count)  # Output: 2
/// ```
#[pyfunction]
fn deduplicate_umis_py(umis: Vec<String>) -> usize {
    deduplicate_umis_rust(&umis)  // Pass by reference, no clone needed
}

/// Process UMI deduplication for multiple cells (Python wrapper)
///
/// Processes UMI counts for multiple cells efficiently, matching Python's
/// `cell_umi_counts_FAST` function exactly.
///
/// # Arguments
/// * `cell_bc_indexes` - List of cell barcode indices
/// * `cell_umis` - List of UMI arrays (parallel to cell_bc_indexes)
/// * `total_cells` - Total number of cells
///
/// # Returns
/// List of UMI counts per cell (length = total_cells)
///
/// # Example (Python)
/// ```python
/// from sheriff_rs import cell_umi_counts_py
/// cell_bc_indexes = [0, 1, 0]
/// cell_umis = [["ATCG"], ["TTTT", "CCCC"], ["GGGG"]]
/// counts = cell_umi_counts_py(cell_bc_indexes, cell_umis, 2)
/// print(counts)  # Output: [2, 2]
/// ```
#[pyfunction]
fn cell_umi_counts_py(
    cell_bc_indexes: Vec<usize>,
    cell_umis: Vec<Vec<String>>,
    total_cells: usize,
) -> Vec<u32> {
    cell_umi_counts_rust(cell_bc_indexes, cell_umis, total_cells)
}

/// Process UMI deduplication for multiple cells in parallel (Python wrapper)
///
/// Parallel version using rayon for multi-core speedup.
///
/// # Arguments
/// * `cell_bc_indexes` - List of cell barcode indices
/// * `cell_umis` - List of UMI arrays (parallel to cell_bc_indexes)
/// * `total_cells` - Total number of cells
///
/// # Returns
/// List of UMI counts per cell (length = total_cells)
///
/// # Performance
/// Expected speedup: 2-8x over sequential version on multi-core systems
#[pyfunction]
fn cell_umi_counts_py_parallel(
    cell_bc_indexes: Vec<usize>,
    cell_umis: Vec<Vec<String>>,
    total_cells: usize,
) -> Vec<u32> {
    cell_umi_counts_rust_parallel(cell_bc_indexes, cell_umis, total_cells)
}

/// Get longest edits from edit clustering (Python wrapper)
///
/// Clusters similar edits and returns only the longest/canonical edit from each cluster.
/// This is a high-performance Rust implementation of Python's `get_longest_edits` function.
///
/// # Arguments
/// * `edits` - List of tuples (chrom, ref_pos, ref_seq, alt_seq, forward, kmer_matches)
///             where kmer_matches is a list of integers
///
/// # Returns
/// List of tuples representing the longest/canonical edits in the same format
///
/// # Performance
/// Expected speedup: 20-100x over Python implementation
/// - Python: 8.9s for 352k reads, 48 minutes for 114M reads
/// - Rust: Sub-second for 352k reads, minutes for 114M reads
///
/// # Example (Python)
/// ```python
/// from sheriff_rs import get_longest_edits_rust
///
/// # Each edit is a tuple: (chrom, ref_pos, ref_seq, alt_seq, forward, kmer_matches)
/// edits = [
///     ("chr1", 1000, "ATCG", "ATCGATCGATCG", True, [1, 2, 3]),
///     ("chr1", 1000, "ATCG", "ATCGATCG", True, [1]),  # Subset, will be removed
///     ("chr1", 2000, "GCTA", "GCTACCCC", False, [4, 5]),
/// ]
///
/// longest = get_longest_edits_rust(edits)
/// # Returns: [("chr1", 1000, "ATCG", "ATCGATCGATCG", True, [1, 2, 3]),
/// #           ("chr1", 2000, "GCTA", "GCTACCCC", False, [4, 5])]
/// ```
#[pyfunction]
fn get_longest_edits_rust(
    edits: Vec<(String, i64, String, String, bool, Vec<usize>)>,
) -> Vec<(String, i64, String, String, bool, Vec<usize>)> {
    // Convert Python tuples to Rust Edit structs
    let rust_edits: Vec<Edit> = edits
        .into_iter()
        .map(|(chrom, ref_pos, ref_seq, alt_seq, forward, kmer_matches)| {
            Edit::new(chrom, ref_pos, ref_seq, alt_seq, forward, kmer_matches)
        })
        .collect();

    // Call Rust implementation
    let longest_edits = get_longest_edits(rust_edits);

    // Convert back to Python tuples
    longest_edits
        .into_iter()
        .map(|edit| {
            (
                edit.chrom,
                edit.ref_pos,
                edit.ref_seq,
                edit.alt_seq,
                edit.forward,
                edit.kmer_matches,
            )
        })
        .collect()
}

/// Gene counts per cell using Rust (Python wrapper)
///
/// Computes a gene x cell UMI count matrix from a BAM file.
/// Matches Python's `python_gene_counts` behavior exactly.
///
/// # Arguments
/// * `bam_path` - Path to BAM file
/// * `barcodes` - List of cell barcodes
/// * `gene_ids` - List of gene IDs (GX/GN/gn tag) in desired column order
/// * `max_reads` - Maximum reads to process (0 = unlimited, default)
///
/// # Returns
/// List of lists (genes x cells) of u32 counts
///
/// # Example (Python)
/// ```python
/// from sheriff_rs import gene_counts_py
/// counts = gene_counts_py("data.bam", ["AAACCCAAG", "AAACCCAAT"], ["GENE1", "GENE2"], 100000)
/// # Returns: [[5, 3], [2, 8]]  # 2 genes x 2 cells
/// ```
#[pyfunction]
#[pyo3(signature = (bam_path, barcodes, gene_ids, max_reads=0))]
fn gene_counts_py(
    bam_path: String,
    barcodes: Vec<String>,
    gene_ids: Vec<String>,
    max_reads: usize,
) -> PyResult<Vec<Vec<u32>>> {
    let counts = gene_counts_per_cell(&bam_path, &barcodes, &gene_ids, max_reads)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    Ok(counts)
}

/// Sheriff-rs Python module
#[pymodule]
fn sheriff_rs(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(filter_bam_rust, m)?)?;
    m.add_function(wrap_pyfunction!(filter_bam_by_barcodes_rust, m)?)?;
    m.add_function(wrap_pyfunction!(filter_bam_by_barcodes_rust_parallel, m)?)?;
    m.add_function(wrap_pyfunction!(filter_bam_by_barcodes_rust_chromosome, m)?)?;
    m.add_function(wrap_pyfunction!(count_kmers_rust, m)?)?;
    m.add_function(wrap_pyfunction!(match_kmer_rust, m)?)?;
    m.add_function(wrap_pyfunction!(deduplicate_umis_py, m)?)?;
    m.add_function(wrap_pyfunction!(cell_umi_counts_py, m)?)?;
    m.add_function(wrap_pyfunction!(cell_umi_counts_py_parallel, m)?)?;
    m.add_function(wrap_pyfunction!(get_longest_edits_rust, m)?)?;
    m.add_function(wrap_pyfunction!(gene_counts_py, m)?)?;
    Ok(())
}
