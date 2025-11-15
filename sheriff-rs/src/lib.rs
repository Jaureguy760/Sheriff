pub mod bam_filter;
pub mod kmer;

// Re-export main functions
pub use bam_filter::{filter_bam_by_barcodes, load_whitelist, FilterResult};
pub use kmer::{kmer_to_num, count_kmers};

// Python bindings (Phase 2)
#[cfg(feature = "python")]
pub mod python;
