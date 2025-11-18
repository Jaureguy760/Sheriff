//! UMI Deduplication Demo
//!
//! Demonstrates the Union-Find based UMI deduplication algorithm

use sheriff_rs::{UnionFind, deduplicate_umis_unionfind, hamming_distance, within_hamming_threshold};

fn main() {
    println!("=== UMI Deduplication Demo ===\n");
    
    // Example 1: Union-Find basics
    println!("1. Union-Find Data Structure:");
    let mut uf = UnionFind::new(5);
    println!("   Created UnionFind with 5 elements");
    
    uf.union(0, 1);
    uf.union(1, 2);
    println!("   After union(0,1) and union(1,2):");
    println!("   - Elements 0, 1, 2 are in the same set");
    println!("   - find(0) == find(2): {}", uf.find(0) == uf.find(2));
    println!();
    
    // Example 2: Hamming distance
    println!("2. Hamming Distance:");
    let umi1 = b"ATCGATCG";
    let umi2 = b"ATCGATCC";
    let dist = hamming_distance(umi1, umi2);
    println!("   Distance between {} and {}: {}", 
             std::str::from_utf8(umi1).unwrap(),
             std::str::from_utf8(umi2).unwrap(),
             dist);
    println!("   Within threshold 1: {}", within_hamming_threshold(umi1, umi2, 1));
    println!();
    
    // Example 3: UMI deduplication
    println!("3. UMI Deduplication (threshold = 1):");
    let umis = vec![
        b"ATCGATCG".as_slice(),
        b"ATCGATCC".as_slice(), // 1 mismatch from umis[0]
        b"ATCGATCA".as_slice(), // 1 mismatch from umis[0]
        b"GCGCGCGC".as_slice(), // Different cluster
        b"GCGCGCGC".as_slice(), // Exact duplicate
        b"GCGCGCGA".as_slice(), // 1 mismatch from umis[3]
    ];
    
    println!("   Input UMIs:");
    for (i, umi) in umis.iter().enumerate() {
        println!("     [{}] {}", i, std::str::from_utf8(umi).unwrap());
    }
    
    let groups = deduplicate_umis_unionfind(&umis, 1);
    
    println!("\n   Output groups: {} clusters", groups.len());
    for (i, group) in groups.iter().enumerate() {
        print!("     Cluster {}: [", i + 1);
        for (j, &idx) in group.iter().enumerate() {
            if j > 0 { print!(", "); }
            print!("{}", idx);
        }
        println!("]");
        
        // Show the UMIs in this cluster
        print!("       UMIs: ");
        for (j, &idx) in group.iter().enumerate() {
            if j > 0 { print!(", "); }
            print!("{}", std::str::from_utf8(umis[idx]).unwrap());
        }
        println!();
    }
    
    println!("\n=== Performance Characteristics ===");
    println!("   - Union-Find operations: O(α(n)) ≈ O(1) amortized");
    println!("   - Hamming distance: O(L) with early exit optimization");
    println!("   - Total deduplication: O(n² × L × α(n))");
    println!("   - Space complexity: O(n)");
    println!("\n   Expected speedup over Python: 3-6x");
}
