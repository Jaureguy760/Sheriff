/// Fast k-mer to integer hashing using bit shifts
///
/// Converts DNA sequence k-mer to a unique integer:
/// A = 0, C = 1, G = 2, T = 3
///
/// # Arguments
/// * `kmer` - Byte slice of DNA sequence
///
/// # Returns
/// u32 hash value
pub fn kmer_to_num(kmer: &[u8]) -> u32 {
    let mut result = 0u32;
    for &nucleotide in kmer {
        result = result.wrapping_shl(2);
        result += match nucleotide {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 0,  // Treat unknown as A
        };
    }
    result
}

/// Count k-mers in sequence
///
/// Creates a frequency array where index is k-mer hash,
/// value is count (saturating at u8::MAX).
///
/// # Arguments
/// * `sequence` - Byte slice of DNA sequence
/// * `k` - K-mer length
///
/// # Returns
/// Vector of k-mer counts (length 4^k)
pub fn count_kmers(sequence: &[u8], k: usize) -> Vec<u8> {
    let array_size = 4usize.pow(k as u32);
    let mut counts = vec![0u8; array_size];

    if sequence.len() < k {
        return counts;
    }

    for window in sequence.windows(k) {
        let hash = kmer_to_num(window) as usize;
        counts[hash] = counts[hash].saturating_add(1);
    }

    counts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_to_num_basic() {
        assert_eq!(kmer_to_num(b"AA"), 0);  // 00
        assert_eq!(kmer_to_num(b"AC"), 1);  // 01
        assert_eq!(kmer_to_num(b"AG"), 2);  // 02
        assert_eq!(kmer_to_num(b"AT"), 3);  // 03
        assert_eq!(kmer_to_num(b"CA"), 4);  // 10
    }

    #[test]
    fn test_kmer_to_num_case_insensitive() {
        assert_eq!(kmer_to_num(b"aa"), kmer_to_num(b"AA"));
        assert_eq!(kmer_to_num(b"CgT"), kmer_to_num(b"CGT"));
    }

    #[test]
    fn test_count_kmers_simple() {
        let seq = b"AACGTT";
        let counts = count_kmers(seq, 2);

        // AA appears once at position 0
        assert_eq!(counts[kmer_to_num(b"AA") as usize], 1);
        // AC appears once at position 1
        assert_eq!(counts[kmer_to_num(b"AC") as usize], 1);
        // CG appears once at position 2
        assert_eq!(counts[kmer_to_num(b"CG") as usize], 1);
        // GT appears once at position 3
        assert_eq!(counts[kmer_to_num(b"GT") as usize], 1);
        // TT appears once at position 4
        assert_eq!(counts[kmer_to_num(b"TT") as usize], 1);
    }

    #[test]
    fn test_count_kmers_empty() {
        let seq = b"A";
        let counts = count_kmers(seq, 2);

        // All zeros (sequence too short)
        assert_eq!(counts.iter().sum::<u8>(), 0);
    }

    #[test]
    fn test_count_kmers_repeats() {
        let seq = b"AAAA";
        let counts = count_kmers(seq, 2);

        // AA appears 3 times (positions 0, 1, 2)
        assert_eq!(counts[kmer_to_num(b"AA") as usize], 3);
    }
}
