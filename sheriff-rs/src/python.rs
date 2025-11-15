use pyo3::prelude::*;
use pyo3::types::PyDict;
use crate::bam_filter::{filter_bam_by_barcodes, filter_bam_by_barcodes_parallel, filter_bam_by_barcodes_chromosome_parallel, load_whitelist};
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

    // Convert whitelist to HashSet
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

/// Sheriff-rs Python module
#[pymodule]
fn sheriff_rs(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(filter_bam_rust, m)?)?;
    m.add_function(wrap_pyfunction!(filter_bam_by_barcodes_rust, m)?)?;
    m.add_function(wrap_pyfunction!(filter_bam_by_barcodes_rust_parallel, m)?)?;
    m.add_function(wrap_pyfunction!(filter_bam_by_barcodes_rust_chromosome, m)?)?;
    m.add_function(wrap_pyfunction!(count_kmers_rust, m)?)?;
    m.add_function(wrap_pyfunction!(match_kmer_rust, m)?)?;
    Ok(())
}
