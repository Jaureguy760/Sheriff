# Sheriff Rust Optimization: Comprehensive Implementation Plan

**Author:** Claude Code (AI Agent)
**Date:** 2025-11-18
**Goal:** Implement 10-100x performance improvements for Sheriff's core bottlenecks

---

## Executive Summary

After deep analysis of Sheriff's Python codebase and research into Rust performance optimizations, I've identified **three critical bottlenecks** with combined **10-100x speedup potential**:

1. **K-mer Matching** (sheriff/count_t7.py:46-158): 10-50x potential
2. **UMI Deduplication** (sheriff/helpers.py:225-253): 5-15x potential
3. **BAM Processing** (sheriff/helpers.py:264-400): 5-10x potential

This plan provides a **phased implementation strategy** starting with quick wins and progressing to advanced optimizations.

---

## Part 1: K-mer Matching Analysis & Optimization

### Current Python Implementation Analysis

**Location:** `sheriff/count_t7.py` lines 46-158

#### Algorithm Overview
```python
class KmerMatcher:
    def __init__(self, k, sequences=None):
        self.k = k
        self.hash_symbol = {"A":0, "C":1, "G":2, "T":3}  # ❌ DICT LOOKUP
        self.match_hash = []  # ❌ PYTHON LIST

    def kmer_to_num(self, kmer):
        # ❌ RECURSIVE (lines 84-89)
        if len(kmer) < 1:
            return 0
        return (4*self.kmer_to_num(kmer[:-1:])) + self.hash_symbol[kmer[-1]]

def match_kmer(bc_kmer_matcher, indel_seq, output_kmer_hash):
    freq_array = np.zeros((4 ** k), dtype=np.uint8)  # ❌ ALLOCATES 4^k ARRAY

    # ❌ NO ROLLING HASH - RECOMPUTES ENTIRE K-MER EVERY WINDOW
    freq_array[[bc_kmer_matcher.kmer_to_num(indel_seq[i: i + k])
                for i in range(len(indel_seq) - k + 1)]] += 1
```

#### Critical Bottlenecks Identified

| Line | Issue | Impact | Fix |
|------|-------|--------|-----|
| 51 | Dict lookup for nucleotides | 2-4x slower | Lookup table (const array) |
| 54 | Python list for hashes | 1.5-2x slower | Vec/FxHashSet |
| 84-89 | Recursive kmer_to_num | 10-20x slower | Iterative + inline |
| 124 | Allocates 4^k array (k=6 → 4096 bytes) | 1.2-1.5x | Array reuse pattern |
| 127-128 | No rolling hash | 3-5x slower | ntHash rolling hash |
| N/A | No SIMD | 2-4x potential | SIMD for batch processing |

### Rust Optimization Strategy

#### Phase 1: Quick Wins (4-14x speedup, <4 hours implementation)

**1.1 Nucleotide Lookup Table (2-4x improvement)**

```rust
// BEFORE (Python equivalent):
// hash_symbol = {"A":0, "C":1, "G":2, "T":3}

// AFTER (Rust):
#[inline(always)]
const fn nucleotide_to_bits(nuc: u8) -> u8 {
    // Lookup table using ASCII values
    // A=65, C=67, G=71, T=84, a=97, c=99, g=103, t=116
    const LOOKUP: [u8; 256] = {
        let mut table = [0u8; 256];
        table[b'A' as usize] = 0;
        table[b'C' as usize] = 1;
        table[b'G' as usize] = 2;
        table[b'T' as usize] = 3;
        table[b'a' as usize] = 0;
        table[b'c' as usize] = 1;
        table[b'g' as usize] = 2;
        table[b't' as usize] = 3;
        table
    };
    LOOKUP[nuc as usize]
}
```

**Why this works:**
- O(1) array indexing vs dict lookup
- `#[inline(always)]` eliminates function call overhead
- Compiler can optimize to single instruction (e.g., `AND` + `SHR`)
- CPU L1 cache friendly (256 bytes fits in single cache line)

**1.2 FxHashSet for Integer Keys (2-3x improvement)**

```rust
use rustc_hash::FxHashSet;  // Import from rustc-hash crate

// BEFORE: std::collections::HashSet (uses SipHash - cryptographically secure)
// AFTER: FxHashSet (uses FxHash - optimized for integer keys)

pub struct KmerMatcher {
    k: usize,
    match_hash: FxHashSet<u32>,  // 2-3x faster lookups for u32 keys
}
```

**Research findings (from web search):**
- FxHash: ~0.4ns per hash (1 CPU cycle) for integers
- SipHash: ~1-2ns per hash (3-5 CPU cycles)
- **Trade-off:** FxHash not DoS-resistant (fine for bioinformatics)
- **Best for:** u32/u64 keys (perfect for k-mer hashes)

**1.3 Iterative kmer_to_num with Inline (3-6x improvement)**

```rust
#[inline]
pub fn kmer_to_num(kmer: &[u8]) -> u32 {
    let mut result = 0u32;
    for &nucleotide in kmer {
        result = result.wrapping_shl(2);  // Multiply by 4
        result += nucleotide_to_bits(nucleotide) as u32;
    }
    result
}
```

**Why faster than Python:**
- No recursion overhead (Python: ~100ns per call)
- `wrapping_shl(2)` compiles to single CPU instruction
- `#[inline]` allows compiler to unroll loop for fixed k
- Zero allocations (Python creates intermediate strings)

**1.4 Array Reuse Pattern (1.1-1.2x improvement)**

```rust
pub struct KmerCounter {
    freq_array: Vec<u8>,  // Reuse across calls
    k: usize,
}

impl KmerCounter {
    pub fn new(k: usize) -> Self {
        Self {
            freq_array: vec![0; 4usize.pow(k as u32)],
            k,
        }
    }

    pub fn count_kmers(&mut self, sequence: &[u8]) -> &[u8] {
        // Zero the array (faster than re-allocating)
        self.freq_array.fill(0);

        for window in sequence.windows(self.k) {
            let hash = kmer_to_num(window);
            self.freq_array[hash as usize] += 1;
        }

        &self.freq_array
    }
}
```

**Combined Phase 1 Expected Speedup: 4-14x**

---

#### Phase 2: Rolling Hash (13-72x total speedup)

**2.1 ntHash Implementation (3-5x additional)**

Based on [SimdMinimizers research](https://curiouscoding.nl/posts/simd-minimizers/):

```rust
/// ntHash: Rolling hash for DNA sequences
/// Algorithm: Mohamadi et al. 2016 - Bioinformatics
///
/// Key insight: Update hash incrementally using only:
/// - Outgoing base (left side)
/// - Incoming base (right side)
///
/// Time complexity: O(1) per window (vs O(k) for naive)
pub struct NtHash {
    k: usize,
    hash: u64,
    mask: u64,
}

impl NtHash {
    const SEED_TAB: [u64; 256] = Self::init_seed_table();

    const fn init_seed_table() -> [u64; 256] {
        let mut table = [0u64; 256];
        // Based on ntHash paper - specific values for A,C,G,T
        table[b'A' as usize] = 0x3c8bfbb395c60474;
        table[b'C' as usize] = 0x3193c18562a02b4c;
        table[b'G' as usize] = 0x20323ed082572324;
        table[b'T' as usize] = 0x295549f54be24456;
        table[b'a' as usize] = 0x3c8bfbb395c60474;
        table[b'c' as usize] = 0x3193c18562a02b4c;
        table[b'g' as usize] = 0x20323ed082572324;
        table[b't' as usize] = 0x295549f54be24456;
        table
    }

    pub fn new(k: usize, initial_kmer: &[u8]) -> Self {
        assert_eq!(initial_kmer.len(), k);

        // Compute initial hash
        let mut hash = 0u64;
        for &nuc in initial_kmer {
            hash ^= Self::SEED_TAB[nuc as usize].rotate_left(k as u32);
        }

        let mask = (1u64 << (2 * k)) - 1;

        Self { k, hash, mask }
    }

    #[inline]
    pub fn roll(&mut self, out_nuc: u8, in_nuc: u8) -> u64 {
        // Remove outgoing base contribution
        self.hash ^= Self::SEED_TAB[out_nuc as usize].rotate_left(self.k as u32);

        // Rotate hash
        self.hash = self.hash.rotate_left(1);

        // Add incoming base contribution
        self.hash ^= Self::SEED_TAB[in_nuc as usize];

        self.hash & self.mask
    }
}
```

**Usage Example:**

```rust
pub fn match_kmer_rolling(
    sequence: &[u8],
    k: usize,
    whitelist: &FxHashSet<u64>,
) -> Vec<u64> {
    let mut matches = Vec::new();

    if sequence.len() < k {
        return matches;
    }

    // Initialize with first k-mer
    let mut nt_hash = NtHash::new(k, &sequence[0..k]);

    if whitelist.contains(&nt_hash.hash) {
        matches.push(nt_hash.hash);
    }

    // Roll through remaining windows - O(1) per window!
    for i in k..sequence.len() {
        let hash = nt_hash.roll(sequence[i - k], sequence[i]);

        if whitelist.contains(&hash) {
            matches.push(hash);
        }
    }

    matches
}
```

**Performance Analysis:**

| Method | Time per k-mer | Sequence 1000bp, k=6 |
|--------|----------------|----------------------|
| Naive (Python recursive) | ~500ns | ~500µs |
| Iterative (Phase 1) | ~50ns | ~50µs |
| Rolling hash | ~1ns | ~1µs |

**Expected Rolling Hash Speedup: 3-5x over Phase 1**

**2.2 SIMD Vectorization (2-4x additional)**

For processing multiple sequences in parallel:

```rust
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

/// SIMD-accelerated batch k-mer matching
/// Processes 4 sequences simultaneously using AVX2
pub unsafe fn match_kmer_simd_batch(
    sequences: &[&[u8]; 4],
    k: usize,
    whitelist: &FxHashSet<u64>,
) -> [Vec<u64>; 4] {
    // Implementation based on SimdMinimizers approach
    // Uses AVX2 256-bit registers to process 4 lanes
    todo!("Advanced optimization - Phase 3")
}
```

**Combined Phase 2 Expected Speedup: 13-72x over Python baseline**

---

## Part 2: UMI Deduplication Analysis & Optimization

### Current Python Implementation Analysis

**Location:** `sheriff/helpers.py` lines 225-253

#### Algorithm Overview

```python
def deduplicate_umis(umi_set):
    """ O(n²) pairwise comparison - CRITICAL BOTTLENECK """
    umis_to_match_umis = {umi: {umi} for umi in umi_set}

    # ❌ O(n²) - compares ALL pairs
    for umi_1, umi_2 in itertools.combinations(umi_set, 2):

        if within_single_mismatch(umi_1, umi_2):  # ❌ Hamming distance in Python
            # ❌ O(n) set union operations in loop
            umis_to_match_umis[umi_1] = umis_to_match_umis[umi_1].union(umis_to_match_umis[umi_2])
            umis_to_match_umis[umi_2] = umis_to_match_umis[umi_2].union(umis_to_match_umis[umi_1])

            # ❌ Nested loop for syncing - O(n³) worst case!
            for umi in umis_to_match_umis[umi_1]:
                umis_to_match_umis[umi_1] = umis_to_match_umis[umi_1].union(umis_to_match_umis[umi])
                # ...

    # ❌ O(n) deduplication at end
    unique_umi_sets = []
    for umi_set_ in umis_to_match_umis.values():
        if umi_set_ not in unique_umi_sets:
            unique_umi_sets.append(umi_set_)

    return unique_umi_sets
```

#### Performance Analysis

For n UMIs:
- **Pairwise comparisons:** n×(n-1)/2 = **O(n²)**
- **Hamming distance per pair:** O(L) where L = UMI length (~12bp)
- **Set unions:** O(n) per union × O(n²) comparisons = **O(n³) worst case**
- **Deduplication:** O(n²) list membership checks

**Total: O(n³) worst case, O(n²) typical**

**Real-world impact:**
- 100 UMIs: ~5,000 comparisons
- 1,000 UMIs: ~500,000 comparisons
- 10,000 UMIs: ~50,000,000 comparisons

### Rust Optimization Strategy

#### Phase 1: Union-Find Algorithm (3-6x speedup)

**Algorithm:** Path-compressed Union-Find with union-by-rank

```rust
use rustc_hash::FxHashMap;

/// Union-Find data structure with path compression
/// Time complexity: O(α(n)) ≈ O(1) amortized per operation
/// where α is inverse Ackermann function
pub struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<usize>,
}

impl UnionFind {
    pub fn new(size: usize) -> Self {
        Self {
            parent: (0..size).collect(),
            rank: vec![0; size],
        }
    }

    /// Find with path compression
    #[inline]
    pub fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            // Path compression: make grandparent the parent
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }

    /// Union by rank
    #[inline]
    pub fn union(&mut self, x: usize, y: usize) -> bool {
        let root_x = self.find(x);
        let root_y = self.find(y);

        if root_x == root_y {
            return false;  // Already in same set
        }

        // Union by rank - attach smaller tree under larger
        if self.rank[root_x] < self.rank[root_y] {
            self.parent[root_x] = root_y;
        } else if self.rank[root_x] > self.rank[root_y] {
            self.parent[root_y] = root_x;
        } else {
            self.parent[root_y] = root_x;
            self.rank[root_x] += 1;
        }

        true
    }
}

/// Hamming distance - optimized with early exit
#[inline]
fn hamming_distance(a: &[u8], b: &[u8]) -> usize {
    a.iter()
        .zip(b.iter())
        .filter(|(x, y)| x != y)
        .count()
}

/// Fast Hamming distance with threshold
#[inline]
fn within_hamming_threshold(a: &[u8], b: &[u8], threshold: usize) -> bool {
    let mut mismatches = 0;
    for (x, y) in a.iter().zip(b.iter()) {
        if x != y {
            mismatches += 1;
            if mismatches > threshold {
                return false;  // Early exit!
            }
        }
    }
    true
}

pub fn deduplicate_umis_unionfind(
    umis: &[&[u8]],
    threshold: usize,
) -> Vec<Vec<usize>> {
    let n = umis.len();
    let mut uf = UnionFind::new(n);

    // Build union-find structure - O(n² × α(n))
    for i in 0..n {
        for j in (i + 1)..n {
            if within_hamming_threshold(umis[i], umis[j], threshold) {
                uf.union(i, j);
            }
        }
    }

    // Group UMIs by connected component - O(n × α(n))
    let mut groups: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
    for i in 0..n {
        let root = uf.find(i);
        groups.entry(root).or_default().push(i);
    }

    groups.into_values().collect()
}
```

**Why Union-Find is faster:**

| Operation | Python (set unions) | Rust Union-Find |
|-----------|---------------------|-----------------|
| Find root | O(n) list traversal | O(α(n)) ≈ O(1) |
| Union | O(n) set copy | O(α(n)) ≈ O(1) |
| Space | O(n²) worst case | O(n) always |

**Expected Phase 1 Speedup: 3-6x**

---

#### Phase 2: Hash-Based Exact Deduplication (2-10x additional)

**Optimization:** Use hash map for exact matches before Hamming distance

```rust
pub fn deduplicate_umis_hybrid(
    umis: &[&[u8]],
    threshold: usize,
) -> Vec<Vec<usize>> {
    let n = umis.len();

    // Phase 1: Group exact duplicates with hash map - O(n)
    let mut exact_groups: FxHashMap<&[u8], Vec<usize>> = FxHashMap::default();
    for (i, &umi) in umis.iter().enumerate() {
        exact_groups.entry(umi).or_default().push(i);
    }

    // Phase 2: Only compare representatives - O(m²) where m << n
    let representatives: Vec<_> = exact_groups.keys().copied().collect();
    let m = representatives.len();

    let mut uf = UnionFind::new(m);

    for i in 0..m {
        for j in (i + 1)..m {
            if within_hamming_threshold(representatives[i], representatives[j], threshold) {
                uf.union(i, j);
            }
        }
    }

    // Phase 3: Merge groups - O(n)
    let mut final_groups: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
    for (rep_idx, &rep) in representatives.iter().enumerate() {
        let root = uf.find(rep_idx);
        let original_indices = &exact_groups[rep];
        final_groups.entry(root).or_default().extend(original_indices);
    }

    final_groups.into_values().collect()
}
```

**Performance:**

For data with 50% exact duplicates:
- n = 10,000 UMIs
- m = 5,000 unique UMIs
- Comparisons: 50M → 12.5M (4x reduction)

**Expected Phase 2 Additional Speedup: 2-10x (data-dependent)**

---

#### Phase 3: BK-Tree for Approximate Matching (5-8x additional)

**Advanced optimization** for ultra-large UMI sets:

```rust
/// BK-Tree: Metric tree for approximate string matching
/// Allows querying all strings within distance d in O(log n) average case
pub struct BKTree {
    root: Option<Box<BKNode>>,
}

struct BKNode {
    umi: Vec<u8>,
    children: FxHashMap<usize, Box<BKNode>>,
}

impl BKTree {
    pub fn new() -> Self {
        Self { root: None }
    }

    pub fn insert(&mut self, umi: Vec<u8>) {
        // Insert into tree based on Hamming distance
        todo!("Advanced optimization - Phase 3")
    }

    pub fn query(&self, umi: &[u8], threshold: usize) -> Vec<&[u8]> {
        // Find all UMIs within threshold distance
        todo!("Advanced optimization - Phase 3")
    }
}
```

**Expected Phase 3 Speedup: 15-60x total over Python baseline**

---

## Part 3: BAM Processing Analysis & Optimization

### Current Python Implementation

**Location:** `sheriff/helpers.py` lines 264-400

#### Bottlenecks Identified

1. **String allocation on every read** (line 320):
   ```python
   cell_barcode = read.get_tag('CB')  # Allocates Python string
   ```

2. **Par-bridge buffers ALL records** (line 311):
   ```python
   with ProcessPoolExecutor(max_workers=n_cpus) as executor:
       # Buffers all futures in memory
   ```

3. **No BGZF threading** for BAM compression

### Rust Optimization Strategy

#### Phase 1: rust-htslib + Rayon (2-3x speedup)

```rust
use rust_htslib::bam::{Read, Reader, Record};
use rayon::prelude::*;

/// Zero-copy cell barcode extraction
#[inline]
fn get_cb_tag_bytes(record: &Record) -> Option<&[u8]> {
    record.aux(b"CB").ok()?.string()
}

/// Parallel BAM processing with streaming
pub fn process_bam_parallel(
    bam_path: &str,
    n_threads: usize,
) -> Result<ProcessedData> {
    let mut reader = Reader::from_path(bam_path)?;
    reader.set_threads(n_threads)?;  // Enable BGZF threading

    // Streaming par-bridge (doesn't buffer all records)
    reader
        .records()
        .par_bridge()
        .filter_map(|result| {
            let record = result.ok()?;
            let cb = get_cb_tag_bytes(&record)?;

            // Process record...
            Some(processed_result)
        })
        .collect()
}
```

**Expected Phase 1 Speedup: 2-3x**

---

## Implementation Roadmap

### Priority Order (ROI-Optimized)

| Phase | Module | Complexity | Time | Expected Speedup | Priority |
|-------|--------|------------|------|------------------|----------|
| 1A | K-mer Phase 1 | Low | 4h | 4-14x | **HIGHEST** |
| 1B | UMI Phase 1 | Medium | 6h | 3-6x | **HIGH** |
| 1C | BAM Phase 1 | Low | 4h | 2-3x | MEDIUM |
| 2A | K-mer Phase 2 | Medium | 8h | 13-72x total | HIGH |
| 2B | UMI Phase 2 | Medium | 6h | 6-60x total | MEDIUM |

### Project Structure

```
Sheriff/
├── sheriff-rs/           # NEW Rust implementation
│   ├── Cargo.toml
│   ├── src/
│   │   ├── lib.rs
│   │   ├── kmer.rs      # K-mer matching optimizations
│   │   ├── umi.rs       # UMI deduplication
│   │   ├── bam.rs       # BAM processing
│   │   └── python.rs    # PyO3 Python bindings
│   └── benches/         # Criterion benchmarks
│       ├── kmer_bench.rs
│       └── umi_bench.rs
├── sheriff/             # Existing Python code
│   ├── count_t7.py
│   └── helpers.py
└── tests/
    └── test_rust_python_equivalence.py
```

### Validation Strategy

1. **Unit tests:** Ensure Rust implementations match Python output exactly
2. **Property tests:** Use `proptest` crate for fuzzing
3. **Benchmarks:** Criterion.rs for performance measurement
4. **Integration tests:** Process real Sheriff data, compare results

---

## API Design

### Python Interface (PyO3)

```python
# sheriff-rs Python API (drop-in replacement)
import sheriff_rs

# K-mer matching (10-50x faster)
matches = sheriff_rs.match_kmer(
    sequence="ATCGATCGATCG",
    k=6,
    whitelist=[27, 54, 108],  # k-mer hashes
    output_hash=True,
    use_rolling_hash=True,  # Enable Phase 2 optimization
)

# UMI deduplication (5-15x faster)
n_unique = sheriff_rs.deduplicate_umis(
    umis=["ATCGATCG", "ATCGATCC", "GCGCGCGC"],
    threshold=1,  # Hamming distance threshold
    method="unionfind",  # or "hybrid", "bktree"
)

# BAM processing (5-10x faster)
results = sheriff_rs.process_bam(
    bam_path="data.bam",
    cell_barcodes=["ACGT-1", "CGTA-1"],
    n_threads=16,
)
```

---

## Benchmarking Plan

### Benchmark Suite

```rust
// benches/kmer_bench.rs
use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};

fn bench_kmer_matching(c: &mut Criterion) {
    let mut group = c.benchmark_group("kmer_matching");

    for seq_len in [100, 1000, 10000] {
        let sequence = generate_random_sequence(seq_len);

        group.bench_with_input(
            BenchmarkId::new("naive", seq_len),
            &sequence,
            |b, seq| b.iter(|| match_kmer_naive(black_box(seq))),
        );

        group.bench_with_input(
            BenchmarkId::new("phase1", seq_len),
            &sequence,
            |b, seq| b.iter(|| match_kmer_phase1(black_box(seq))),
        );

        group.bench_with_input(
            BenchmarkId::new("rolling_hash", seq_len),
            &sequence,
            |b, seq| b.iter(|| match_kmer_rolling(black_box(seq))),
        );
    }

    group.finish();
}

criterion_group!(benches, bench_kmer_matching);
criterion_main!(benches);
```

### Expected Benchmark Results

```
kmer_matching/naive/100        time:   [50.2 µs 50.5 µs 50.8 µs]
kmer_matching/phase1/100       time:   [3.8 µs 3.9 µs 4.0 µs]
                               change: [-92.8% -92.3% -91.9%] (12.9x faster)
kmer_matching/rolling_hash/100 time:   [0.96 µs 0.98 µs 1.0 µs]
                               change: [-98.1% -98.0% -97.9%] (51.5x faster)
```

---

## Risk Mitigation

### Potential Issues & Solutions

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Rust-Python type conversion overhead | Medium | Low | Use zero-copy PyO3 conversions |
| FxHash collision issues | Low | Medium | Add tests with pathological inputs |
| Rolling hash numerical stability | Low | High | Use well-tested ntHash implementation |
| Union-Find correctness | Medium | Critical | Extensive property testing |
| BAM parallel processing deadlocks | Low | High | Use Rayon's proven parallel iterators |

---

## Next Steps

1. **Implement K-mer Phase 1** (HIGHEST priority)
   - [ ] Create `sheriff-rs` Rust crate
   - [ ] Implement nucleotide lookup table
   - [ ] Implement FxHashSet integration
   - [ ] Add PyO3 Python bindings
   - [ ] Write tests comparing Python vs Rust output
   - [ ] Benchmark and validate 4-14x speedup

2. **Implement UMI Phase 1** (HIGH priority)
   - [ ] Implement Union-Find with path compression
   - [ ] Implement Hamming distance with early exit
   - [ ] Add Python bindings
   - [ ] Validate correctness
   - [ ] Benchmark and validate 3-6x speedup

3. **Phase 2 optimizations** (based on Phase 1 results)

---

## References

### Research Papers
- **ntHash:** Mohamadi et al. 2016, "ntHash: recursive nucleotide hashing" - Bioinformatics
- **SimdMinimizers:** 2025, "Computing random minimizers, fast" - bioRxiv
- **UMI-tools:** Smith et al. 2017, "UMI-tools: modeling sequencing errors in Unique Molecular Identifiers" - Genome Research

### Rust Resources
- **Rust Performance Book:** https://nnethercote.github.io/perf-book/
- **FxHash:** https://github.com/rust-lang/rustc-hash
- **rust-htslib:** https://github.com/rust-bio/rust-htslib
- **PyO3:** https://pyo3.rs/

### Sheriff Codebase
- **K-mer matching:** `sheriff/count_t7.py` lines 46-158
- **UMI deduplication:** `sheriff/helpers.py` lines 225-253
- **BAM processing:** `sheriff/helpers.py` lines 264-400

---

**End of Plan**

This comprehensive plan provides both:
1. **Quick wins** (Phase 1) for immediate 4-14x speedups with low risk
2. **Advanced optimizations** (Phase 2-3) for 10-100x total speedups

All optimizations are grounded in:
- Deep analysis of the existing Python code
- Research into state-of-the-art Rust performance techniques
- Proven algorithms from bioinformatics literature

Ready for implementation! 🚀
