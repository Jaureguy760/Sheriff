//! BAM file processing module
//!
//! This module provides utilities for reading and processing BAM (Binary Alignment Map) files.

/// A simple BAM record representation
#[derive(Debug, Clone)]
pub struct BamRecord {
    /// Query name
    pub qname: String,
    /// Sequence
    pub seq: String,
    /// Quality scores
    pub qual: String,
}

impl BamRecord {
    /// Create a new BamRecord
    pub fn new(qname: String, seq: String, qual: String) -> Self {
        BamRecord { qname, seq, qual }
    }

    /// Get the sequence length
    pub fn seq_len(&self) -> usize {
        self.seq.len()
    }
}

/// A BAM file reader
#[derive(Debug)]
pub struct BamReader {
    /// File path
    path: String,
}

impl BamReader {
    /// Create a new BamReader
    pub fn new(path: String) -> Self {
        BamReader { path }
    }

    /// Get the file path
    pub fn path(&self) -> &str {
        &self.path
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bam_record_creation() {
        let record = BamRecord::new(
            "read1".to_string(),
            "ACGTACGTACGT".to_string(),
            "IIIIIIIIIIII".to_string(),
        );
        assert_eq!(record.seq_len(), 12);
    }

    #[test]
    fn test_bam_reader_creation() {
        let reader = BamReader::new("/path/to/file.bam".to_string());
        assert_eq!(reader.path(), "/path/to/file.bam");
    }
}
