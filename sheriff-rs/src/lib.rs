//! Sheriff-rs: High-performance Rust implementation for k-mer and UMI processing
//!
//! This crate provides efficient implementations of k-mer analysis and UMI
//! (Unique Molecular Identifier) deduplication with optional PyO3 bindings
//! for Python integration.

pub mod kmer;
pub mod umi;
pub mod bam;

#[cfg(feature = "python")]
pub mod python;

// Re-export commonly used items
pub use kmer::*;
pub use umi::*;
pub use bam::*;

/// Get library version
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_version() {
        assert!(!version().is_empty());
    }
}
