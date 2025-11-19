//! BAM file processing module with rust-htslib
//!
//! This module provides high-performance utilities for reading and processing
//! BAM (Binary Alignment Map) files using rust-htslib for zero-copy operations.

use rust_htslib::bam::{Read, Reader, Record};
use std::path::Path;
use rayon::prelude::*;

/// Error type for BAM processing
#[derive(Debug)]
pub enum BamError {
    /// IO error
    IoError(std::io::Error),
    /// BAM parsing error
    ParseError(String),
    /// Missing required tag
    MissingTag(String),
}

impl std::fmt::Display for BamError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BamError::IoError(e) => write!(f, "IO error: {}", e),
            BamError::ParseError(s) => write!(f, "Parse error: {}", s),
            BamError::MissingTag(s) => write!(f, "Missing tag: {}", s),
        }
    }
}

impl std::error::Error for BamError {}

impl From<std::io::Error> for BamError {
    fn from(err: std::io::Error) -> Self {
        BamError::IoError(err)
    }
}

impl From<rust_htslib::errors::Error> for BamError {
    fn from(err: rust_htslib::errors::Error) -> Self {
        BamError::ParseError(format!("{:?}", err))
    }
}

pub type Result<T> = std::result::Result<T, BamError>;

/// Zero-copy BAM record processor
pub struct BamProcessor {
    reader: Reader,
}

impl BamProcessor {
    /// Create a new BAM processor from a file path
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        let reader = Reader::from_path(path)?;
        Ok(BamProcessor { reader })
    }

    /// Enable multi-threaded BAM decompression
    pub fn set_threads(&mut self, n_threads: usize) -> Result<()> {
        self.reader
            .set_threads(n_threads)
            .map_err(|e| BamError::ParseError(format!("Failed to set threads: {:?}", e)))
    }

    /// Get the underlying reader
    pub fn reader(&self) -> &Reader {
        &self.reader
    }

    /// Get mutable reference to the underlying reader
    pub fn reader_mut(&mut self) -> &mut Reader {
        &mut self.reader
    }
}

/// Extract UMI tag from BAM record (zero-copy)
///
/// Sheriff uses the "pN" tag for UMI sequences.
/// Returns a byte slice pointing directly into the BAM record (no allocation).
#[inline]
pub fn get_umi_tag(record: &Record) -> Option<&[u8]> {
    use rust_htslib::bam::record::Aux;
    match record.aux(b"pN").ok()? {
        Aux::String(s) => Some(s.as_bytes()),
        _ => None,
    }
}

/// Extract cell barcode tag from BAM record (zero-copy)
///
/// Sheriff uses the "CB" tag for cell barcodes.
/// Returns a byte slice pointing directly into the BAM record (no allocation).
#[inline]
pub fn get_cell_barcode_tag(record: &Record) -> Option<&[u8]> {
    use rust_htslib::bam::record::Aux;
    match record.aux(b"CB").ok()? {
        Aux::String(s) => Some(s.as_bytes()),
        _ => None,
    }
}

/// Extract both UMI and cell barcode tags (zero-copy)
///
/// Returns (umi, cell_barcode) as byte slices.
#[inline]
pub fn get_umi_and_barcode(record: &Record) -> Option<(&[u8], &[u8])> {
    let umi = get_umi_tag(record)?;
    let cell_barcode = get_cell_barcode_tag(record)?;
    Some((umi, cell_barcode))
}

/// Extract sequence from BAM record (zero-copy via indexing)
///
/// Returns the DNA sequence as a Vec<u8>.
/// Note: rust-htslib doesn't provide zero-copy sequence access due to
/// BAM's 4-bit encoding, but this is still much faster than Python.
#[inline]
pub fn get_sequence(record: &Record) -> Vec<u8> {
    record.seq().as_bytes()
}

/// Check if record has required tags for Sheriff processing
#[inline]
pub fn has_required_tags(record: &Record) -> bool {
    get_umi_tag(record).is_some() && get_cell_barcode_tag(record).is_some()
}

/// Process BAM records with a callback function
///
/// This provides an iterator-style interface for processing BAM records
/// with zero-copy tag extraction.
///
/// # Example
/// ```no_run
/// use sheriff_rs::bam::{BamProcessor, process_records, get_umi_and_barcode};
///
/// let mut processor = BamProcessor::new("data.bam").unwrap();
/// let mut count = 0;
///
/// process_records(&mut processor, |record| {
///     if let Some((umi, cb)) = get_umi_and_barcode(record) {
///         count += 1;
///     }
///     Ok(())
/// }).unwrap();
/// ```
pub fn process_records<F>(processor: &mut BamProcessor, mut callback: F) -> Result<()>
where
    F: FnMut(&Record) -> Result<()>,
{
    for result in processor.reader_mut().records() {
        let record = result?;
        callback(&record)?;
    }
    Ok(())
}

/// Statistics from BAM processing
#[derive(Debug, Default, Clone)]
pub struct BamStats {
    /// Total reads processed
    pub total_reads: usize,
    /// Reads with UMI tag
    pub reads_with_umi: usize,
    /// Reads with cell barcode tag
    pub reads_with_cb: usize,
    /// Reads with both UMI and CB tags
    pub reads_with_both: usize,
    /// Total bases sequenced
    pub total_bases: usize,
}

impl BamStats {
    /// Create new empty stats
    pub fn new() -> Self {
        Self::default()
    }

    /// Update stats with a BAM record
    #[inline]
    pub fn update(&mut self, record: &Record) {
        self.total_reads += 1;
        self.total_bases += record.seq_len();

        let has_umi = get_umi_tag(record).is_some();
        let has_cb = get_cell_barcode_tag(record).is_some();

        if has_umi {
            self.reads_with_umi += 1;
        }
        if has_cb {
            self.reads_with_cb += 1;
        }
        if has_umi && has_cb {
            self.reads_with_both += 1;
        }
    }
}

/// Collect statistics from a BAM file
pub fn collect_stats<P: AsRef<Path>>(path: P) -> Result<BamStats> {
    let mut processor = BamProcessor::new(path)?;
    let mut stats = BamStats::new();

    process_records(&mut processor, |record| {
        stats.update(record);
        Ok(())
    })?;

    Ok(stats)
}

/// Process BAM records in parallel using Rayon
///
/// This function reads all records into memory first, then processes them in parallel.
/// This is efficient for medium-sized BAM files but may use significant memory for large files.
///
/// # Example
/// ```no_run
/// use sheriff_rs::bam::{process_records_parallel, get_umi_and_barcode};
/// use std::sync::atomic::{AtomicUsize, Ordering};
///
/// let count = AtomicUsize::new(0);
/// process_records_parallel("data.bam", |record| {
///     if let Some((_umi, _cb)) = get_umi_and_barcode(record) {
///         count.fetch_add(1, Ordering::Relaxed);
///     }
/// }).unwrap();
/// ```
pub fn process_records_parallel<P, F>(path: P, callback: F) -> Result<()>
where
    P: AsRef<Path>,
    F: Fn(&Record) + Sync + Send,
{
    let mut processor = BamProcessor::new(path)?;

    // Collect all records into a Vec (this uses memory!)
    let records: Vec<Record> = processor
        .reader_mut()
        .records()
        .collect::<std::result::Result<Vec<_>, _>>()?;

    // Process in parallel using Rayon
    records.par_iter().for_each(|record| {
        callback(record);
    });

    Ok(())
}

/// Collect statistics from a BAM file in parallel
///
/// This is faster than `collect_stats` for large BAM files with many cores.
/// Uses atomic operations to safely accumulate statistics across threads.
pub fn collect_stats_parallel<P: AsRef<Path>>(path: P) -> Result<BamStats> {
    use std::sync::atomic::{AtomicUsize, Ordering};

    let total_reads = AtomicUsize::new(0);
    let reads_with_umi = AtomicUsize::new(0);
    let reads_with_cb = AtomicUsize::new(0);
    let reads_with_both = AtomicUsize::new(0);
    let total_bases = AtomicUsize::new(0);

    process_records_parallel(path, |record| {
        total_reads.fetch_add(1, Ordering::Relaxed);
        total_bases.fetch_add(record.seq_len(), Ordering::Relaxed);

        let has_umi = get_umi_tag(record).is_some();
        let has_cb = get_cell_barcode_tag(record).is_some();

        if has_umi {
            reads_with_umi.fetch_add(1, Ordering::Relaxed);
        }
        if has_cb {
            reads_with_cb.fetch_add(1, Ordering::Relaxed);
        }
        if has_umi && has_cb {
            reads_with_both.fetch_add(1, Ordering::Relaxed);
        }
    })?;

    Ok(BamStats {
        total_reads: total_reads.load(Ordering::Relaxed),
        reads_with_umi: reads_with_umi.load(Ordering::Relaxed),
        reads_with_cb: reads_with_cb.load(Ordering::Relaxed),
        reads_with_both: reads_with_both.load(Ordering::Relaxed),
        total_bases: total_bases.load(Ordering::Relaxed),
    })
}

/// Group records by cell barcode in parallel
///
/// Returns a HashMap where keys are cell barcodes and values are vectors of (UMI, sequence) tuples.
/// This is useful for per-cell processing in Sheriff.
///
/// Note: This loads all data into memory, so use with caution on very large BAM files.
pub fn group_by_cell_parallel<P: AsRef<Path>>(
    path: P,
) -> Result<std::collections::HashMap<Vec<u8>, Vec<(Vec<u8>, Vec<u8>)>>> {
    use std::sync::Mutex;

    let mut processor = BamProcessor::new(path)?;

    // Collect all records
    let records: Vec<Record> = processor
        .reader_mut()
        .records()
        .collect::<std::result::Result<Vec<_>, _>>()?;

    // Use DashMap for concurrent HashMap access (or Mutex<HashMap>)
    let cell_map: Mutex<std::collections::HashMap<Vec<u8>, Vec<(Vec<u8>, Vec<u8>)>>> =
        Mutex::new(std::collections::HashMap::new());

    records.par_iter().for_each(|record| {
        if let Some((umi, cb)) = get_umi_and_barcode(record) {
            let umi_owned = umi.to_vec();
            let cb_owned = cb.to_vec();
            let seq = get_sequence(record);

            let mut map = cell_map.lock().unwrap();
            map.entry(cb_owned)
                .or_insert_with(Vec::new)
                .push((umi_owned, seq));
        }
    });

    Ok(cell_map.into_inner().unwrap())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bam_error_display() {
        let err = BamError::MissingTag("pN".to_string());
        assert_eq!(format!("{}", err), "Missing tag: pN");
    }

    #[test]
    fn test_bam_stats_creation() {
        let stats = BamStats::new();
        assert_eq!(stats.total_reads, 0);
        assert_eq!(stats.reads_with_umi, 0);
    }

    // Integration test with real BAM file (requires example data)
    #[test]
    #[ignore] // Only run when example_data is available
    fn test_real_bam_processing() {
        let bam_path = "../example_data/barcode_headAligned_anno.sorted.edit_regions_200kb.bam";

        if !std::path::Path::new(bam_path).exists() {
            return; // Skip if file doesn't exist
        }

        let stats = collect_stats(bam_path).expect("Failed to collect stats");

        assert!(stats.total_reads > 0, "Should have processed some reads");
        assert!(stats.reads_with_both > 0, "Should have reads with UMI and CB tags");

        println!("BAM Stats:");
        println!("  Total reads: {}", stats.total_reads);
        println!("  Reads with UMI: {}", stats.reads_with_umi);
        println!("  Reads with CB: {}", stats.reads_with_cb);
        println!("  Reads with both: {}", stats.reads_with_both);
        println!("  Total bases: {}", stats.total_bases);
    }

    #[test]
    #[ignore]
    fn test_zero_copy_tag_extraction() {
        let bam_path = "../example_data/barcode_headAligned_anno.sorted.edit_regions_200kb.bam";

        if !std::path::Path::new(bam_path).exists() {
            return;
        }

        let mut processor = BamProcessor::new(bam_path).expect("Failed to open BAM");
        let mut found_tags = false;

        process_records(&mut processor, |record| {
            if let Some((umi, cb)) = get_umi_and_barcode(record) {
                // These are zero-copy byte slices
                assert!(umi.len() > 0, "UMI should not be empty");
                assert!(cb.len() > 0, "Cell barcode should not be empty");
                found_tags = true;
            }
            Ok(())
        }).expect("Failed to process records");

        assert!(found_tags, "Should have found at least one record with tags");
    }

    #[test]
    #[ignore]
    fn test_parallel_stats_collection() {
        let bam_path = "../example_data/barcode_headAligned_anno.sorted.edit_regions_200kb.bam";

        if !std::path::Path::new(bam_path).exists() {
            return;
        }

        // Get sequential stats
        let stats_seq = collect_stats(bam_path).expect("Failed to collect stats");

        // Get parallel stats
        let stats_par = collect_stats_parallel(bam_path).expect("Failed to collect parallel stats");

        // Results should be identical
        assert_eq!(stats_seq.total_reads, stats_par.total_reads);
        assert_eq!(stats_seq.reads_with_umi, stats_par.reads_with_umi);
        assert_eq!(stats_seq.reads_with_cb, stats_par.reads_with_cb);
        assert_eq!(stats_seq.reads_with_both, stats_par.reads_with_both);
        assert_eq!(stats_seq.total_bases, stats_par.total_bases);

        println!("Parallel BAM Stats Match!");
        println!("  Total reads: {}", stats_par.total_reads);
    }

    #[test]
    #[ignore]
    fn test_group_by_cell_parallel() {
        let bam_path = "../example_data/barcode_headAligned_anno.sorted.edit_regions_200kb.bam";

        if !std::path::Path::new(bam_path).exists() {
            return;
        }

        let cell_map = group_by_cell_parallel(bam_path).expect("Failed to group by cell");

        assert!(!cell_map.is_empty(), "Should have found cells");

        let total_reads: usize = cell_map.values().map(|v| v.len()).sum();
        println!("Grouped {} reads into {} cells", total_reads, cell_map.len());

        // Check a random cell has valid data
        if let Some((cb, reads)) = cell_map.iter().next() {
            assert!(!cb.is_empty(), "Cell barcode should not be empty");
            assert!(!reads.is_empty(), "Cell should have reads");
            assert!(!reads[0].0.is_empty(), "UMI should not be empty");
            assert!(!reads[0].1.is_empty(), "Sequence should not be empty");
        }
    }
}
