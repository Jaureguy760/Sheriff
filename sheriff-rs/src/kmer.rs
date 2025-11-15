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

/// Convert k-mer hash back to DNA sequence string
///
/// # Arguments
/// * `num` - K-mer hash value
/// * `k` - K-mer length
///
/// # Returns
/// K-mer as String
pub fn num_to_kmer(mut num: usize, k: usize) -> String {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut kmer = vec![0u8; k];

    for i in (0..k).rev() {
        kmer[i] = bases[num & 0b11];  // Get last 2 bits
        num >>= 2;                     // Shift right by 2
    }

    String::from_utf8(kmer).unwrap()
}

/// Match k-mers in DNA sequence against whitelist
///
/// Finds all k-mers in sequence that appear in the whitelist.
/// Handles 'N' bases by skipping k-mers containing them.
///
/// # Arguments
/// * `sequence` - DNA sequence string
/// * `k` - K-mer length
/// * `whitelist` - Optional set of k-mer hashes to match against
/// * `output_hash` - Return hashes (true) or k-mer strings (false)
///
/// # Returns
/// Vector of matching k-mer hashes or strings, empty if no matches
pub fn match_kmer(
    sequence: &str,
    k: usize,
    whitelist: Option<&std::collections::HashSet<usize>>,
    output_hash: bool,
) -> Vec<KmerMatch> {
    if sequence.len() < k {
        return Vec::new();
    }

    // Create frequency array
    let array_size = 4usize.pow(k as u32);
    let mut freq = vec![0u8; array_size];

    // Count k-mer occurrences, skipping k-mers with 'N'
    let seq_bytes = sequence.as_bytes();
    for i in 0..=(sequence.len() - k) {
        let kmer = &seq_bytes[i..i + k];

        // Skip k-mers containing 'N' or 'n'
        if kmer.iter().any(|&b| b == b'N' || b == b'n') {
            continue;
        }

        let hash = kmer_to_num(kmer) as usize;
        freq[hash] = freq[hash].saturating_add(1);
    }

    // Find non-zero k-mers
    let mut matches: Vec<usize> = freq.iter()
        .enumerate()
        .filter(|(_, &count)| count > 0)
        .map(|(idx, _)| idx)
        .collect();

    // Filter by whitelist if provided
    if let Some(whitelist_set) = whitelist {
        matches.retain(|&hash| whitelist_set.contains(&hash));
    }

    // Return in requested format
    if output_hash {
        matches.into_iter().map(KmerMatch::Hash).collect()
    } else {
        matches.into_iter()
            .map(|hash| KmerMatch::String(num_to_kmer(hash, k)))
            .collect()
    }
}

/// K-mer match result - either hash or string
#[derive(Debug, Clone, PartialEq)]
pub enum KmerMatch {
    Hash(usize),
    String(String),
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

    #[test]
    fn test_num_to_kmer() {
        assert_eq!(num_to_kmer(0, 2), "AA");  // 00
        assert_eq!(num_to_kmer(1, 2), "AC");  // 01
        assert_eq!(num_to_kmer(2, 2), "AG");  // 02
        assert_eq!(num_to_kmer(3, 2), "AT");  // 03
        assert_eq!(num_to_kmer(4, 2), "CA");  // 10
        assert_eq!(num_to_kmer(27, 4), "ACGT"); // 0*64 + 1*16 + 2*4 + 3
    }

    #[test]
    fn test_match_kmer_no_whitelist() {
        let seq = "AAACGTTT";
        let matches = match_kmer(seq, 4, None, true);

        // Should find: AAAC, AACG, ACGT, CGTT, GTTT (5 unique k-mers)
        assert_eq!(matches.len(), 5);

        // Verify they're hashes
        for m in matches {
            assert!(matches!(m, KmerMatch::Hash(_)));
        }
    }

    #[test]
    fn test_match_kmer_with_whitelist() {
        let seq = "AAACGTTT";
        let whitelist: std::collections::HashSet<usize> = [27].iter().copied().collect(); // "ACGT" = 27
        let matches = match_kmer(seq, 4, Some(&whitelist), true);

        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0], KmerMatch::Hash(27));
    }

    #[test]
    fn test_match_kmer_output_strings() {
        let seq = "ACGT";
        let matches = match_kmer(seq, 4, None, false);

        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0], KmerMatch::String("ACGT".to_string()));
    }

    #[test]
    fn test_match_kmer_with_n() {
        let seq = "AAANGTTT"; // Contains N
        let matches = match_kmer(seq, 4, None, true);

        // Should skip: AAAN, NANG, ANGT, NGTT
        // Should find: GTTT
        assert!(matches.len() < 5); // Less than if all were counted
    }
}
