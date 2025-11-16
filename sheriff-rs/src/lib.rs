pub mod bam_filter;
pub mod kmer;
pub mod umi;
pub mod edit_clustering;

// Re-export main functions
pub use bam_filter::{filter_bam_by_barcodes, load_whitelist, FilterResult};
pub use kmer::{kmer_to_num, count_kmers};
pub use umi::{deduplicate_umis_rust, cell_umi_counts_rust, cell_umi_counts_rust_parallel};
pub use edit_clustering::{Edit, get_longest_edits, bio_edit_distance};

// Python bindings (Phase 2)
#[cfg(feature = "python")]
pub mod python;
