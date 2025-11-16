/// Edit clustering for Sheriff bioinformatics pipeline
///
/// This module implements high-performance edit clustering to identify canonical
/// T7 edits from sequencing data. The algorithm clusters similar edits together
/// and selects the longest representative edit from each cluster.
///
/// # Performance
/// - Target: 20-100x speedup over Python implementation
/// - Python baseline: 8.9s for 352k reads, 48 minutes for 114M reads
/// - Uses rust-bio for sequence alignment
/// - Optimized with aggressive compiler flags (LTO, opt-level=3)

use bio::alignment::pairwise::*;
use std::cmp::Ordering;
use std::collections::HashSet;

/// ReadEdit represents a single edit detected in a sequencing read
///
/// This corresponds to Python's namedtuple:
/// ReadEdit = namedtuple("ReadEdit", ["chrom", "ref_pos", "ref_seq", "alt_seq", "forward", "kmer_matches"])
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Edit {
    pub chrom: String,
    pub ref_pos: i64,
    pub ref_seq: String,
    pub alt_seq: String,
    pub forward: bool,
    pub kmer_matches: Vec<usize>,
}

impl Edit {
    pub fn new(
        chrom: String,
        ref_pos: i64,
        ref_seq: String,
        alt_seq: String,
        forward: bool,
        kmer_matches: Vec<usize>,
    ) -> Self {
        Edit {
            chrom,
            ref_pos,
            ref_seq,
            alt_seq,
            forward,
            kmer_matches,
        }
    }
}

/// Detect if a sequence has homopolymers (3+ consecutive identical bases)
fn has_homopolymer(seq: &str) -> bool {
    if seq.len() < 3 {
        return false;
    }
    let bytes = seq.as_bytes();
    let mut count = 1;
    for i in 1..bytes.len() {
        if bytes[i] == bytes[i - 1] {
            count += 1;
            if count >= 3 {
                return true;
            }
        } else {
            count = 1;
        }
    }
    false
}

/// Replace homopolymer runs (3+ consecutive identical bases) with single base
fn collapse_homopolymers(seq: &str) -> String {
    if seq.is_empty() {
        return String::new();
    }
    let bytes = seq.as_bytes();
    let mut result = Vec::with_capacity(seq.len());
    let mut prev_char = bytes[0];
    let mut count = 1;

    for i in 1..bytes.len() {
        if bytes[i] == prev_char {
            count += 1;
        } else {
            // Only keep if run length < 3, otherwise collapse to single
            if count >= 3 {
                result.push(prev_char);
            } else {
                for _ in 0..count {
                    result.push(prev_char);
                }
            }
            prev_char = bytes[i];
            count = 1;
        }
    }
    // Handle last run
    if count >= 3 {
        result.push(prev_char);
    } else {
        for _ in 0..count {
            result.push(prev_char);
        }
    }
    String::from_utf8(result).unwrap_or_default()
}

/// Calculate biologically-relevant edit distance between two sequences
///
/// This function implements the `bio_edit_distance` logic from Python:
/// - Uses local alignment (Smith-Waterman via rust-bio)
/// - Counts mismatches only until shorter sequence is fully aligned
/// - Returns 0 if one sequence is a perfect substring of the other
///
/// # Arguments
/// * `seq_a` - First sequence
/// * `seq_b` - Second sequence
/// * `start_from_first_smallest_seq_aln` - If true, start counting from first char of shorter seq
/// * `alns_to_compare` - Optional limit on alignment length to compare
///
/// # Returns
/// Edit distance between sequences
pub fn bio_edit_distance(
    seq_a: &str,
    seq_b: &str,
    start_from_first_smallest_seq_aln: bool,
    alns_to_compare: Option<usize>,
) -> usize {
    // Create scoring scheme matching Python's pairwise2.align.localms parameters:
    // match=1, mismatch=-1, gap_open=-0.5, gap_extend=-0.5
    let mut scoring = Scoring::new(-1, -1, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
    scoring.gap_open = -1;  // gap_open is a field, not a method!
    scoring.gap_extend = -1;

    let mut aligner = Aligner::with_scoring(scoring);
    let alignment = aligner.local(seq_a.as_bytes(), seq_b.as_bytes());

    // Determine which sequence is shorter
    let seq_lens = [seq_a.len(), seq_b.len()];
    let shorter_seq_len = seq_lens.iter().min().unwrap();
    let shorter_seq_index = if seq_a.len() <= seq_b.len() { 0 } else { 1 };

    // Reconstruct alignment strings from operations
    let (aligned_a, aligned_b) = reconstruct_alignment(seq_a, seq_b, &alignment);

    // Count mismatches until all of shorter sequence is aligned
    let mut n_seen = 0;
    let mut edit_dist = 0;
    let mut start_count = !start_from_first_smallest_seq_aln;

    let max_compare = alns_to_compare.unwrap_or(aligned_a.len()).min(aligned_a.len());

    let aln_seqs = [aligned_a.as_bytes(), aligned_b.as_bytes()];

    for i in 0..max_compare {
        if start_count && aligned_a.as_bytes()[i] != aligned_b.as_bytes()[i] {
            edit_dist += 1;
        }

        if aln_seqs[shorter_seq_index][i] != b'-' {
            n_seen += 1;
            start_count = true;

            if n_seen == *shorter_seq_len {
                break;
            }
        }
    }

    edit_dist
}

/// Reconstruct aligned sequences from alignment operations
///
/// This matches Python's pairwise2.align.localms output format:
/// - Leading dashes for unaligned prefix of the OTHER sequence
/// - Internal gaps from alignment operations
/// - Trailing dashes for unaligned suffix of the OTHER sequence
///
/// Python example:
/// seq1 = 'AGGGGCCCTTACCATACCCC' (len=20)
/// seq2 = 'CGCCACTGGAGACTGGGTTAG' (len=21)
/// Result:
/// seqA: '---------AG---GGGCCCTTACCATACCCC'  (9 leading dashes)
/// seqB: 'CGCCACTGGAGACTGGG---TTAG--------'  (8 trailing dashes)
///
/// If alignment.xstart = 9, alignment.ystart = 0:
/// - seq2[0:9] appears first, with dashes in seqA
/// - Then the aligned portion
/// - Then remaining seq1, with dashes in seqB
fn reconstruct_alignment(seq_a: &str, seq_b: &str, alignment: &bio::alignment::Alignment) -> (String, String) {
    let mut aligned_a = String::new();
    let mut aligned_b = String::new();

    let seq_a_bytes = seq_a.as_bytes();
    let seq_b_bytes = seq_b.as_bytes();

    // Add leading padding for unaligned prefix
    // In bio-rs aligner.local(x, y):
    //   x = first arg (seq_a), y = second arg (seq_b)
    //   xstart = position where alignment starts in seq_a
    //   ystart = position where alignment starts in seq_b

    // Safety check: ensure indices are within bounds
    let xstart = alignment.xstart.min(seq_a_bytes.len());
    let ystart = alignment.ystart.min(seq_b_bytes.len());

    // If ystart > 0, seq_b has leading unaligned chars - add dashes to aligned_a
    for idx in 0..ystart {
        aligned_a.push('-');
        aligned_b.push(seq_b_bytes[idx] as char);
    }

    // If xstart > 0, seq_a has leading unaligned chars - add dashes to aligned_b
    for idx in 0..xstart {
        aligned_a.push(seq_a_bytes[idx] as char);
        aligned_b.push('-');
    }

    let mut i = xstart;  // Current position in seq_a
    let mut j = ystart;  // Current position in seq_b

    // Process alignment operations
    for op in &alignment.operations {
        match op {
            bio::alignment::AlignmentOperation::Match => {
                if i < seq_a_bytes.len() && j < seq_b_bytes.len() {
                    aligned_a.push(seq_a_bytes[i] as char);
                    aligned_b.push(seq_b_bytes[j] as char);
                    i += 1;
                    j += 1;
                }
            }
            bio::alignment::AlignmentOperation::Subst => {
                if i < seq_a_bytes.len() && j < seq_b_bytes.len() {
                    aligned_a.push(seq_a_bytes[i] as char);
                    aligned_b.push(seq_b_bytes[j] as char);
                    i += 1;
                    j += 1;
                }
            }
            bio::alignment::AlignmentOperation::Del => {
                if i < seq_a_bytes.len() {
                    aligned_a.push(seq_a_bytes[i] as char);
                    aligned_b.push('-');
                    i += 1;
                }
            }
            bio::alignment::AlignmentOperation::Ins => {
                if j < seq_b_bytes.len() {
                    aligned_a.push('-');
                    aligned_b.push(seq_b_bytes[j] as char);
                    j += 1;
                }
            }
            _ => {}
        }
    }

    // Add trailing padding for unaligned suffix
    // Add remaining characters from seq_a and dashes to seq_b
    for idx in i..seq_a_bytes.len() {
        aligned_a.push(seq_a_bytes[idx] as char);
        aligned_b.push('-');
    }

    // Add remaining characters from seq_b and dashes to seq_a
    for idx in j..seq_b_bytes.len() {
        aligned_b.push(seq_b_bytes[idx] as char);
        aligned_a.push('-');
    }

    (aligned_a, aligned_b)
}

/// Get longest edits from a set of edits, removing sub-edits
///
/// This is the Rust implementation of Python's `get_longest_edits` function.
/// The algorithm:
/// 1. Sorts edits by alt_seq length
/// 2. Compares each pair (shorter vs longer)
/// 3. Calculates edit distance with homopolymer correction
/// 4. Groups similar edits, keeping only the longest
///
/// # Biological Logic
/// - Edits with different orientations (forward/reverse) are distinct
/// - Edits at different positions are distinct
/// - Handles homopolymer sequencing errors (consecutive identical bases)
/// - Allows small differences at 3' end (sequencing adapters/TSO)
/// - Edit distance threshold ≤ 2 for clustering
///
/// # Arguments
/// * `edits` - Vector of Edit structs
///
/// # Returns
/// Vector of canonical (longest) edits
pub fn get_longest_edits(mut edits: Vec<Edit>) -> Vec<Edit> {
    if edits.is_empty() {
        return vec![];
    }

    // Sort by alt_seq length (ascending)
    edits.sort_by_key(|e| e.alt_seq.len());

    // Use HashSet for O(1) membership checks instead of O(n) Vec::contains()
    let mut longest_edits_set: HashSet<Edit> = HashSet::new();
    let mut sub_edits_set: HashSet<Edit> = HashSet::new();

    for i in 0..edits.len() {
        let edit_1 = &edits[i];

        for j in (i + 1)..edits.len() {
            let edit_2 = &edits[j];

            // Different orientation or position -> definitely different alleles
            if edit_1.forward != edit_2.forward || edit_1.ref_pos != edit_2.ref_pos {
                // Add both as longest edits if not already sub-edits
                if !sub_edits_set.contains(edit_1) {
                    longest_edits_set.insert(edit_1.clone());
                }
                if !sub_edits_set.contains(edit_2) {
                    longest_edits_set.insert(edit_2.clone());
                }
                continue;
            }

            // Extract sequences for comparison (remove reference portion)
            let reflen = edit_1.ref_seq.len();
            let (edit_1_seq, edit_2_seq) = if !edit_1.forward && !edit_2.forward {
                // Reverse direction: cut off reference and reverse
                let seq1 = &edit_1.alt_seq[reflen..];
                let seq2 = &edit_2.alt_seq[reflen..];
                let rev1: String = seq1.chars().rev().collect();
                let rev2: String = seq2.chars().rev().collect();
                (rev1, rev2)
            } else if edit_1.forward && edit_2.forward {
                // Forward direction: already in correct orientation
                let seq1 = &edit_1.alt_seq[..edit_1.alt_seq.len() - reflen];
                let seq2 = &edit_2.alt_seq[..edit_2.alt_seq.len() - reflen];
                (seq1.to_string(), seq2.to_string())
            } else {
                continue; // Should not happen given earlier check
            };

            // Calculate initial edit distance
            let mut dist_between_seqs = bio_edit_distance(&edit_1_seq, &edit_2_seq, true, None);

            // Homopolymer correction if distance > 1
            if dist_between_seqs > 1 {
                // Check for homopolymers (3+ consecutive identical bases)
                let has_homopolymer_1 = has_homopolymer(&edit_1_seq);
                let has_homopolymer_2 = has_homopolymer(&edit_2_seq);

                if has_homopolymer_1 || has_homopolymer_2 {
                    // Replace homopolymers with single base
                    let edit_1_seq_homo_fix = collapse_homopolymers(&edit_1_seq);
                    let edit_2_seq_homo_fix = collapse_homopolymers(&edit_2_seq);

                    // Recalculate distance with correction
                    dist_between_seqs = bio_edit_distance(
                        &edit_1_seq_homo_fix,
                        &edit_2_seq_homo_fix,
                        true,
                        None
                    ).saturating_sub(1); // -1 because homopolymer adds extra error
                }
            }

            // Check 3' end if still different
            if dist_between_seqs > 2 {
                // Compare last 10 bases at 3' end (reversed to start from 3' end)
                let rev1: String = edit_1_seq.chars().rev().collect();
                let rev2: String = edit_2_seq.chars().rev().collect();

                let three_prime_edit_dist = bio_edit_distance(&rev1, &rev2, false, Some(10));

                // Only 1 mismatch at 3' end -> consider same sequence
                if three_prime_edit_dist <= 1 {
                    dist_between_seqs = three_prime_edit_dist;
                }
            }

            // Determine if edits should be clustered together
            if dist_between_seqs <= 2 {
                // Determine which is the sub-edit (shorter one, or fewer kmer matches if same length)
                let subedit = if edit_1.alt_seq.len() != edit_2.alt_seq.len() {
                    if edit_1.alt_seq.len() < edit_2.alt_seq.len() {
                        edit_1
                    } else {
                        edit_2
                    }
                } else {
                    // Same length, choose based on kmer matches
                    if edit_1.kmer_matches.len() >= edit_2.kmer_matches.len() {
                        edit_2
                    } else {
                        edit_1
                    }
                };

                // Remove subedit from longest_edits if present (O(1) operation)
                longest_edits_set.remove(subedit);

                // Add subedit to sub_edits set (O(1) operation)
                sub_edits_set.insert(subedit.clone());

                // Add the long_edit to longest_edits if not already a subedit
                let actual_long_edit = if edit_1 == subedit { edit_2 } else { edit_1 };
                if !sub_edits_set.contains(actual_long_edit) {
                    longest_edits_set.insert(actual_long_edit.clone());
                }
            } else {
                // They are different - add both as longest edits if not sub-edits
                if !sub_edits_set.contains(edit_1) {
                    longest_edits_set.insert(edit_1.clone());
                }
                if !sub_edits_set.contains(edit_2) {
                    longest_edits_set.insert(edit_2.clone());
                }
            }
        }
    }

    // Convert HashSet to sorted Vec for deterministic output
    let mut longest_edits: Vec<Edit> = longest_edits_set.into_iter().collect();
    longest_edits.sort_by(|a, b| {
        match a.chrom.cmp(&b.chrom) {
            Ordering::Equal => a.ref_pos.cmp(&b.ref_pos),
            other => other,
        }
    });

    longest_edits
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bio_edit_distance_identical() {
        let seq_a = "ATCGATCG";
        let seq_b = "ATCGATCG";
        let dist = bio_edit_distance(seq_a, seq_b, true, None);
        assert_eq!(dist, 0);
    }

    #[test]
    fn test_bio_edit_distance_substring() {
        // One sequence is substring of another -> should be 0
        let seq_a = "ATCG";
        let seq_b = "ATCGATCG";
        let dist = bio_edit_distance(seq_a, seq_b, true, None);
        assert_eq!(dist, 0);
    }

    #[test]
    fn test_bio_edit_distance_one_mismatch() {
        let seq_a = "ATCGATCG";
        let seq_b = "ATCGATCC"; // last base different
        let dist = bio_edit_distance(seq_a, seq_b, true, None);
        assert!(dist >= 1);
    }

    #[test]
    fn test_edit_clustering_empty() {
        let edits = vec![];
        let result = get_longest_edits(edits);
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_edit_clustering_single() {
        let edits = vec![Edit::new(
            "chr1".to_string(),
            1000,
            "ATCG".to_string(),
            "ATCGATCGATCG".to_string(),
            true,
            vec![1, 2, 3],
        )];
        let result = get_longest_edits(edits);
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_edit_clustering_different_positions() {
        // Edits at different positions should be kept separate
        let edits = vec![
            Edit::new(
                "chr1".to_string(),
                1000,
                "ATCG".to_string(),
                "ATCGATCGATCG".to_string(),
                true,
                vec![1, 2],
            ),
            Edit::new(
                "chr1".to_string(),
                2000, // Different position
                "ATCG".to_string(),
                "ATCGATCGATCG".to_string(),
                true,
                vec![1, 2],
            ),
        ];
        let result = get_longest_edits(edits);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_edit_clustering_different_orientations() {
        // Edits with different orientations should be kept separate
        let edits = vec![
            Edit::new(
                "chr1".to_string(),
                1000,
                "ATCG".to_string(),
                "ATCGATCGATCG".to_string(),
                true, // forward
                vec![1, 2],
            ),
            Edit::new(
                "chr1".to_string(),
                1000,
                "ATCG".to_string(),
                "ATCGATCGATCG".to_string(),
                false, // reverse
                vec![1, 2],
            ),
        ];
        let result = get_longest_edits(edits);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_edit_clustering_subset_removal() {
        // Shorter edit should be removed when it's a subset of longer edit
        let edits = vec![
            Edit::new(
                "chr1".to_string(),
                1000,
                "ATCG".to_string(),
                "ATCGATCG".to_string(), // Short
                true,
                vec![1],
            ),
            Edit::new(
                "chr1".to_string(),
                1000,
                "ATCG".to_string(),
                "ATCGATCGATCGATCG".to_string(), // Long
                true,
                vec![1, 2, 3],
            ),
        ];
        let result = get_longest_edits(edits);
        // Should keep only the longer edit
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].alt_seq.len(), "ATCGATCGATCGATCG".len());
    }

    #[test]
    fn test_homopolymer_correction() {
        // Test that homopolymer sequences are correctly clustered
        // These differ only in the number of consecutive 'C's
        let edits = vec![
            Edit::new(
                "chr1".to_string(),
                1000,
                "ATCG".to_string(),
                "ATCGATAATACTCTCCCTATTCACTCTGCGT".to_string(),
                true,
                vec![1],
            ),
            Edit::new(
                "chr1".to_string(),
                1000,
                "ATCG".to_string(),
                "ATCGATAATACTCTCCCCCCATTCACT".to_string(), // More C's
                true,
                vec![1],
            ),
        ];
        let result = get_longest_edits(edits);
        // Should cluster these together and keep the longer one
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].alt_seq, "ATCGATAATACTCTCCCTATTCACTCTGCGT");
    }
}
