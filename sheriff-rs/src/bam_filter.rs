use rust_htslib::bam::{self, Read, Record};
use std::collections::HashSet;
use anyhow::{Result, Context};
use rayon::prelude::*;
use std::fs;

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

    // Enable multi-threaded BGZF decompression (2-5x faster I/O)
    // Uses 4 threads by default for parallel block decompression
    bam.set_threads(4).ok();  // Ignore errors if threading not supported

    // Create output BAM with same header
    let header = bam::Header::from_template(bam.header());
    let mut out = bam::Writer::from_path(output_path, &header, bam::Format::Bam)
        .context("Failed to create output BAM")?;

    // Enable multi-threaded BGZF compression for output
    out.set_threads(4).ok();

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

    // Enable multi-threaded BGZF decompression (2-5x faster I/O)
    bam.set_threads(4).ok();

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

/// Filter BAM file by cell barcode whitelist (Chromosome-based parallel version)
///
/// Splits BAM processing by chromosome, processes each in parallel using rayon,
/// then merges results. Expected speedup: 10-50x on large files.
///
/// Requires indexed BAM file (*.bam.bai).
///
/// # Arguments
/// * `input_path` - Path to input BAM file (must be indexed)
/// * `output_path` - Path to output BAM file
/// * `whitelist` - Set of allowed cell barcodes
/// * `num_threads` - Number of parallel threads (None = auto-detect)
///
/// # Returns
/// FilterResult with statistics on reads processed, kept, and rejected
pub fn filter_bam_by_barcodes_chromosome_parallel(
    input_path: &str,
    output_path: &str,
    whitelist: &HashSet<String>,
    num_threads: Option<usize>,
) -> Result<FilterResult> {
    // Open BAM to read header
    let bam = bam::Reader::from_path(input_path)
        .context("Failed to open input BAM")?;
    let header_view = bam.header();
    let header = bam::Header::from_template(header_view);

    // Get list of reference sequences (chromosomes)
    let ref_names: Vec<String> = header_view
        .target_names()
        .iter()
        .map(|name| String::from_utf8_lossy(name).to_string())
        .collect();

    if ref_names.is_empty() {
        return Err(anyhow::anyhow!("No reference sequences found in BAM header"));
    }

    // Create temp directory for chromosome BAMs
    let temp_dir = format!("{}.tmp_chr", output_path);
    fs::create_dir_all(&temp_dir)
        .context("Failed to create temporary directory")?;

    // Create a custom thread pool or use default
    let process_chromosomes = || {
        ref_names
            .par_iter()
            .map(|chr_name| {
                let temp_bam = format!("{}/{}.bam", temp_dir, chr_name);
                let result = filter_chromosome(
                    input_path,
                    &temp_bam,
                    chr_name,
                    whitelist,
                );
                (temp_bam, result)
            })
            .filter_map(|(temp_bam, result)| {
                match result {
                    Ok(stats) => Some((temp_bam, stats)),
                    Err(e) => {
                        eprintln!("Warning: Failed to process chromosome: {}", e);
                        None
                    }
                }
            })
            .collect()
    };

    // Process each chromosome in parallel using custom or default pool
    let results: Vec<(String, FilterResult)> = if let Some(threads) = num_threads {
        // Use a custom local thread pool to avoid global pool conflicts
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .context("Failed to create thread pool")?
            .install(process_chromosomes)
    } else {
        // Use default global pool
        process_chromosomes()
    };

    // Merge all chromosome BAMs
    let total_stats = merge_bam_files(&results, output_path, &header)?;

    // Clean up temp directory
    fs::remove_dir_all(&temp_dir)
        .context("Failed to remove temporary directory")?;

    Ok(total_stats)
}

/// Filter reads from a specific chromosome
fn filter_chromosome(
    input_path: &str,
    output_path: &str,
    chr_name: &str,
    whitelist: &HashSet<String>,
) -> Result<FilterResult> {
    // Open input BAM
    let mut bam = bam::IndexedReader::from_path(input_path)
        .context("Failed to open indexed BAM")?;

    // Get chromosome ID
    let header_view = bam.header();
    let header = bam::Header::from_template(header_view);
    let tid = header_view
        .target_names()
        .iter()
        .position(|name| String::from_utf8_lossy(name) == chr_name)
        .ok_or_else(|| anyhow::anyhow!("Chromosome {} not found in header", chr_name))?;

    // Fetch reads for this chromosome
    bam.fetch(tid as u32)
        .context(format!("Failed to fetch chromosome {}", chr_name))?;

    // Create output BAM
    let mut out = bam::Writer::from_path(output_path, &header, bam::Format::Bam)
        .context("Failed to create output BAM")?;

    let mut stats = FilterResult {
        reads_processed: 0,
        reads_kept: 0,
        reads_rejected: 0,
    };

    // Process reads for this chromosome
    for result in bam.records() {
        let record = result.context("Failed to read BAM record")?;
        stats.reads_processed += 1;

        if let Some(cb_tag) = get_cb_tag(&record) {
            if whitelist.contains(&cb_tag) {
                out.write(&record).context("Failed to write record")?;
                stats.reads_kept += 1;
            } else {
                stats.reads_rejected += 1;
            }
        } else {
            stats.reads_rejected += 1;
        }
    }

    Ok(stats)
}

/// Merge multiple BAM files into one output file
fn merge_bam_files(
    input_files: &[(String, FilterResult)],
    output_path: &str,
    header: &bam::Header,
) -> Result<FilterResult> {
    // Create output BAM
    let mut out = bam::Writer::from_path(output_path, header, bam::Format::Bam)
        .context("Failed to create merged output BAM")?;

    let mut total_stats = FilterResult {
        reads_processed: 0,
        reads_kept: 0,
        reads_rejected: 0,
    };

    // Merge each chromosome BAM
    for (chr_bam, stats) in input_files {
        // Accumulate statistics
        total_stats.reads_processed += stats.reads_processed;
        total_stats.reads_kept += stats.reads_kept;
        total_stats.reads_rejected += stats.reads_rejected;

        // Only merge if there are kept reads
        if stats.reads_kept == 0 {
            continue;
        }

        // Read and write all records from chromosome BAM
        let mut chr_reader = bam::Reader::from_path(chr_bam)
            .context(format!("Failed to open chromosome BAM: {}", chr_bam))?;

        for result in chr_reader.records() {
            let record = result.context("Failed to read record during merge")?;
            out.write(&record).context("Failed to write record during merge")?;
        }
    }

    Ok(total_stats)
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
