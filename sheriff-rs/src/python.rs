use pyo3::prelude::*;
use pyo3::types::PyDict;
use crate::bam_filter::{filter_bam_by_barcodes, load_whitelist};
use std::collections::HashSet;

/// Filter BAM file by cell barcode whitelist (Python wrapper)
///
/// # Arguments
/// * `input_path` - Path to input BAM file
/// * `output_path` - Path to output BAM file
/// * `whitelist_path` - Path to whitelist text file
///
/// # Returns
/// Dict with keys: reads_processed, reads_kept, reads_rejected
#[pyfunction]
fn filter_bam_rust(
    input_path: String,
    output_path: String,
    whitelist_path: String,
) -> PyResult<PyObject> {
    // Load whitelist
    let whitelist = load_whitelist(&whitelist_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    // Filter BAM
    let result = filter_bam_by_barcodes(&input_path, &output_path, &whitelist)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    // Return dict
    Python::with_gil(|py| {
        let dict = PyDict::new(py);
        dict.set_item("reads_processed", result.reads_processed)?;
        dict.set_item("reads_kept", result.reads_kept)?;
        dict.set_item("reads_rejected", result.reads_rejected)?;
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

/// Sheriff-rs Python module
#[pymodule]
fn sheriff_rs(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(filter_bam_rust, m)?)?;
    m.add_function(wrap_pyfunction!(count_kmers_rust, m)?)?;
    Ok(())
}
