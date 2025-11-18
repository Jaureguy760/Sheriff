//! K-mer Processing Demo
//!
//! Demonstrates the Phase 1 optimizations for k-mer matching.

use rustc_hash::FxHashSet;
use sheriff_rs::kmer::{nucleotide_to_bits, kmer_to_num, match_kmer, KmerCounter};

fn main() {
    println!("=== Sheriff K-mer Phase 1 Optimizations Demo ===\n");

    // Demo 1: Nucleotide to bits conversion
    println!("1. Nucleotide Lookup Table:");
    println!("   A -> {}", nucleotide_to_bits(b'A'));
    println!("   C -> {}", nucleotide_to_bits(b'C'));
    println!("   G -> {}", nucleotide_to_bits(b'G'));
    println!("   T -> {}", nucleotide_to_bits(b'T'));
    println!("   (Case-insensitive: a -> {})\n", nucleotide_to_bits(b'a'));

    // Demo 2: K-mer to numeric hash
    println!("2. K-mer to Numeric Hash:");
    let kmers: Vec<&[u8]> = vec![b"A", b"AC", b"ACG", b"ACGT"];
    for kmer in &kmers {
        let hash = kmer_to_num(kmer);
        println!("   {} -> {}", String::from_utf8_lossy(kmer), hash);
    }
    println!();

    // Demo 3: K-mer matching with whitelist
    println!("3. K-mer Matching (FxHashSet):");
    let sequence = b"ACGTACGTGGGGACGT";
    let k = 4;
    
    // Build whitelist with specific k-mers
    let mut whitelist = FxHashSet::default();
    whitelist.insert(kmer_to_num(b"ACGT")); // 27
    whitelist.insert(kmer_to_num(b"GGGG")); // 170
    
    println!("   Sequence: {}", String::from_utf8_lossy(sequence));
    println!("   Whitelist: ACGT (hash={}), GGGG (hash={})", 
             kmer_to_num(b"ACGT"), kmer_to_num(b"GGGG"));
    
    let matches = match_kmer(sequence, k, &whitelist);
    println!("   Found {} matches: {:?}\n", matches.len(), matches);

    // Demo 4: K-mer frequency counting with array reuse
    println!("4. K-mer Frequency Counter (Array Reuse):");
    let mut counter = KmerCounter::new(4);
    
    let seq1 = b"ACGTACGTACGT";
    println!("   Sequence 1: {}", String::from_utf8_lossy(seq1));
    let freqs1 = counter.count_kmers(seq1);
    let acgt_hash = kmer_to_num(b"ACGT");
    println!("   ACGT appears {} times", freqs1[acgt_hash as usize]);
    
    let seq2 = b"AAAAAAAA";
    println!("\n   Sequence 2: {}", String::from_utf8_lossy(seq2));
    let freqs2 = counter.count_kmers(seq2);
    let aaaa_hash = kmer_to_num(b"AAAA");
    println!("   AAAA appears {} times", freqs2[aaaa_hash as usize]);
    println!("   ACGT appears {} times (array was reused/zeroed)", 
             freqs2[acgt_hash as usize]);

    println!("\n=== Performance Characteristics ===");
    println!("✓ Nucleotide lookup: O(1) const array access");
    println!("✓ K-mer hashing: Iterative, zero allocations");
    println!("✓ FxHashSet: ~0.4ns per hash (vs ~1-2ns for SipHash)");
    println!("✓ Array reuse: Eliminates repeated allocations");
    println!("\nExpected speedup: 4-14x over Python baseline");
}
