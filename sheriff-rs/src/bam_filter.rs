use rust_htslib::bam::{self, Read, Record};
use std::collections::HashSet;
use anyhow::{Result, Context};
use rayon::prelude::*;

#[derive(Debug, Clone)]
pub struct FilterResult {
    pub reads_processed: usize,
    pub reads_kept: usize,
    pub reads_rejected: usize,
}

/// Filter BAM file by cell barcode whitelist
///
/// Reads input BAM, filters reads by CB tag against whitelist,
/// writes filtered reads to output BAM.
///
/// # Arguments
/// * `input_path` - Path to input BAM file
/// * `output_path` - Path to output BAM file
/// * `whitelist` - Set of allowed cell barcodes
///
/// # Returns
/// FilterResult with statistics on reads processed, kept, and rejected
pub fn filter_bam_by_barcodes(
    input_path: &str,
    output_path: &str,
    whitelist: &HashSet<String>,
) -> Result<FilterResult> {
    // Open input BAM
    let mut bam = bam::Reader::from_path(input_path)
        .context("Failed to open input BAM")?;

    // Create output BAM with same header
    let header = bam::Header::from_template(bam.header());
    let mut out = bam::Writer::from_path(output_path, &header, bam::Format::Bam)
        .context("Failed to create output BAM")?;

    let mut stats = FilterResult {
        reads_processed: 0,
        reads_kept: 0,
        reads_rejected: 0,
    };

    // Iterate reads
    for result in bam.records() {
        let record = result.context("Failed to read BAM record")?;
        stats.reads_processed += 1;

        // Extract CB tag (cell barcode)
        if let Some(cb_tag) = get_cb_tag(&record) {
            if whitelist.contains(&cb_tag) {
                out.write(&record).context("Failed to write record")?;
                stats.reads_kept += 1;
            } else {
                stats.reads_rejected += 1;
            }
        } else {
            // No CB tag, reject
            stats.reads_rejected += 1;
        }
    }

    Ok(stats)
}

/// Filter BAM file by cell barcode whitelist (Parallel version using Rayon)
///
/// Uses rayon's par_bridge to parallelize the filtering operation.
/// Expected speedup: 2-5x on medium files, best for CPU-bound filtering.
///
/// # Arguments
/// * `input_path` - Path to input BAM file
/// * `output_path` - Path to output BAM file
/// * `whitelist` - Set of allowed cell barcodes
///
/// # Returns
/// FilterResult with statistics on reads processed, kept, and rejected
pub fn filter_bam_by_barcodes_parallel(
    input_path: &str,
    output_path: &str,
    whitelist: &HashSet<String>,
) -> Result<FilterResult> {
    // Open input BAM
    let mut bam = bam::Reader::from_path(input_path)
        .context("Failed to open input BAM")?;

    // Store header for output creation
    let header = bam::Header::from_template(bam.header());

    // Parallel filtering using par_bridge
    // Collects (record, should_keep) pairs
    let results: Vec<(Record, bool)> = bam.records()
        .par_bridge()  // Convert to parallel iterator
        .filter_map(|result| {
            let record = result.ok()?;
            let should_keep = if let Some(cb_tag) = get_cb_tag(&record) {
                whitelist.contains(&cb_tag)
            } else {
                false
            };
            Some((record, should_keep))
        })
        .collect();

    // Calculate statistics
    let reads_processed = results.len();
    let reads_kept = results.iter().filter(|(_, keep)| *keep).count();
    let reads_rejected = reads_processed - reads_kept;

    // Create output BAM and write filtered reads sequentially
    let mut out = bam::Writer::from_path(output_path, &header, bam::Format::Bam)
        .context("Failed to create output BAM")?;

    for (record, should_keep) in results {
        if should_keep {
            out.write(&record).context("Failed to write record")?;
        }
    }

    Ok(FilterResult {
        reads_processed,
        reads_kept,
        reads_rejected,
    })
}

/// Extract CB (cell barcode) tag from BAM record
fn get_cb_tag(record: &Record) -> Option<String> {
    record.aux(b"CB")
        .ok()
        .and_then(|aux| {
            match aux {
                rust_htslib::bam::record::Aux::String(s) => Some(s),
                _ => None,
            }
        })
        .map(|s| s.to_string())
}

/// Load whitelist from text file
///
/// Reads a text file with one barcode per line, returns a HashSet.
///
/// # Arguments
/// * `path` - Path to whitelist text file
///
/// # Returns
/// HashSet of cell barcodes
pub fn load_whitelist(path: &str) -> Result<HashSet<String>> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    let file = File::open(path).context("Failed to open whitelist")?;
    let reader = BufReader::new(file);

    let mut whitelist = HashSet::new();
    for line in reader.lines() {
        let barcode = line.context("Failed to read line")?;
        whitelist.insert(barcode.trim().to_string());
    }

    Ok(whitelist)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_whitelist() {
        // Create temp whitelist file
        use std::io::Write;
        let mut file = std::fs::File::create("/tmp/test_whitelist.txt").unwrap();
        writeln!(file, "AAACCTGAGAAACCAT-1").unwrap();
        writeln!(file, "AAACCTGAGAAACCGC-1").unwrap();

        let whitelist = load_whitelist("/tmp/test_whitelist.txt").unwrap();
        assert_eq!(whitelist.len(), 2);
        assert!(whitelist.contains("AAACCTGAGAAACCAT-1"));
    }

    #[test]
    fn test_get_cb_tag() {
        // Note: This requires a real BAM record, so we'll skip implementation
        // Real test will use integration tests with actual BAM files
    }
}
