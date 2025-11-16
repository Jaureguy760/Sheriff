# Sheriff-RS Optimization Engineering Plan

**Date:** 2025-11-16
**Estimated Total Effort:** 4-8 weeks
**Expected Performance Gain:** Additional 2-10x on top of current 10-44x

---

## Current State (Verified Correct)

| Component | Current Speedup | Correctness |
|-----------|----------------|-------------|
| UMI Deduplication | 33-44x | PASS |
| Edit Clustering | 16-25x | PASS |
| K-mer Matching | 13.7x | PASS |

**Key Finding:** Edit struct already has `#[derive(Hash)]` on line 20 - the O(n³) fix is a quick win!

---

## Phase 1: CRITICAL - O(n³) → O(n²) Fix

**Timeline:** 1-2 days
**Impact:** 10-100x faster edit clustering for large datasets
**Risk:** LOW

### The Problem

```rust
// Current: O(n) lookup done O(n²) times = O(n³)
let mut longest_edits: Vec<Edit> = Vec::new();    // BAD
let mut sub_edits: Vec<Edit> = Vec::new();        // BAD

if !sub_edits.contains(edit_1) { ... }  // O(n) scan
if !longest_edits.contains(edit_1) { ... }  // O(n) scan
longest_edits.retain(|e| e != subedit);  // O(n) scan
```

### The Fix

```rust
use std::collections::HashSet;

// New: O(1) lookup
let mut longest_edits: HashSet<Edit> = HashSet::new();  // GOOD
let mut sub_edits: HashSet<Edit> = HashSet::new();      // GOOD

if !sub_edits.contains(edit_1) { ... }  // O(1) lookup
if !longest_edits.contains(edit_1) { ... }  // O(1) lookup
longest_edits.remove(subedit);  // O(1) remove
```

### Files to Modify

1. `sheriff-rs/src/edit_clustering.rs`
   - Line 251: `let mut longest_edits: HashSet<Edit> = HashSet::new();`
   - Line 252: `let mut sub_edits: HashSet<Edit> = HashSet::new();`
   - Line 264: `longest_edits.insert(edit_1.clone());`
   - Line 267: `longest_edits.insert(edit_2.clone());`
   - Line 362: `longest_edits.remove(subedit);`
   - Lines 366, 372, 377, 380: Convert `.push()` to `.insert()`
   - Lines 387-394: Convert HashSet back to Vec for return, then sort

### Breaking Changes

**None** - Internal optimization only, same API and results.

### Validation

```bash
python validate_rust_correctness.py  # Must still pass
python benchmark_comprehensive.py --dataset medium  # Measure improvement
```

---

## Phase 2: Aligner Reuse (HIGH PRIORITY)

**Timeline:** 3-5 days
**Impact:** Eliminate n² allocations
**Risk:** MEDIUM

### The Problem

```rust
// Creates NEW Aligner for every comparison (up to 3n² times!)
pub fn bio_edit_distance(...) -> usize {
    let mut aligner = Aligner::with_scoring(scoring);  // ALLOCATION
    let alignment = aligner.local(...);
    ...
}
```

### The Fix - Thread-Local Aligner

```rust
use std::cell::RefCell;
use bio::alignment::pairwise::*;

thread_local! {
    static ALIGNER: RefCell<Aligner<MatchParams>> = RefCell::new({
        let mut scoring = Scoring::new(-1, -1, |a: u8, b: u8| if a == b { 1 } else { -1 });
        scoring.gap_open = -1;
        scoring.gap_extend = -1;
        Aligner::with_capacity_and_scoring(500, 500, scoring)  // Pre-allocate
    });
}

pub fn bio_edit_distance(seq_a: &str, seq_b: &str, ...) -> usize {
    ALIGNER.with(|aligner| {
        let mut aligner = aligner.borrow_mut();
        let alignment = aligner.local(seq_a.as_bytes(), seq_b.as_bytes());
        reconstruct_alignment(seq_a, seq_b, &alignment)
        // ... rest of logic
    })
}
```

### Alternative - Pass Aligner as Parameter

```rust
pub fn get_longest_edits(mut edits: Vec<Edit>) -> Vec<Edit> {
    // Create once, reuse for all comparisons
    let mut scoring = Scoring::new(...);
    let mut aligner = Aligner::with_capacity_and_scoring(500, 500, scoring);

    for i in 0..edits.len() {
        for j in (i+1)..edits.len() {
            let dist = bio_edit_distance_with_aligner(&seq1, &seq2, &mut aligner, ...);
        }
    }
}
```

### Files to Modify

1. `sheriff-rs/src/edit_clustering.rs`
   - Lines 122-170: Refactor `bio_edit_distance()` to accept aligner
   - Lines 291, 305, 320: Update call sites

### Validation

- Benchmark to determine optimal capacity (500x500 is a guess)
- Profile memory allocations before/after

---

## Phase 3: String Allocation Reduction

**Timeline:** 3-5 days
**Impact:** 2-4x fewer heap allocations
**Risk:** LOW

### The Problem

```rust
// Inside O(n²) loop - allocates 2-4 strings per iteration
let rev1: String = seq1.chars().rev().collect();  // MALLOC
let rev2: String = seq2.chars().rev().collect();  // MALLOC
// ... later ...
let rev1: String = edit_1_seq.chars().rev().collect();  // MALLOC AGAIN
let rev2: String = edit_2_seq.chars().rev().collect();  // MALLOC AGAIN
```

### The Fix - Pre-compute Once

```rust
pub fn get_longest_edits(mut edits: Vec<Edit>) -> Vec<Edit> {
    edits.sort_by_key(|e| e.alt_seq.len());

    // Pre-compute all reversed sequences ONCE (O(n) instead of O(n²))
    let reversed_alts: Vec<String> = edits.iter()
        .map(|e| e.alt_seq.chars().rev().collect())
        .collect();

    for i in 0..edits.len() {
        let edit_1 = &edits[i];
        let rev_1 = &reversed_alts[i];  // No allocation

        for j in (i+1)..edits.len() {
            let edit_2 = &edits[j];
            let rev_2 = &reversed_alts[j];  // No allocation
            // ... use rev_1, rev_2 directly
        }
    }
}
```

### Alternative - In-Place Byte Comparison

```rust
// DNA is ASCII, so byte operations are safe
fn compare_reversed_suffixes(s1: &str, s2: &str, len: usize) -> usize {
    let b1 = s1.as_bytes();
    let b2 = s2.as_bytes();
    let mut mismatches = 0;

    for i in 0..len.min(b1.len()).min(b2.len()) {
        if b1[b1.len() - 1 - i] != b2[b2.len() - 1 - i] {
            mismatches += 1;
        }
    }
    mismatches
}
```

### Files to Modify

1. `sheriff-rs/src/edit_clustering.rs`
   - Lines 274-288: Remove inline reversals
   - Lines 317-318: Use pre-computed or byte-based comparison

---

## Phase 4: BAM Streaming I/O

**Timeline:** 5-7 days
**Impact:** Memory: GB → MB
**Risk:** MEDIUM

### The Problem

```rust
// Loads ENTIRE BAM file into memory
let results: Vec<(Record, bool)> = bam.records()
    .par_bridge()
    .filter_map(...)
    .collect();  // 100M reads = 20-50 GB RAM
```

### The Fix - Streaming with Record Reuse

```rust
pub fn filter_bam_streaming(
    input_path: &str,
    output_path: &str,
    whitelist: &HashSet<String>,
) -> Result<FilterResult> {
    let mut bam = bam::Reader::from_path(input_path)?;
    bam.set_threads(4).ok();

    let header = bam::Header::from_template(bam.header());
    let mut out = bam::Writer::from_path(output_path, &header, bam::Format::Bam)?;
    out.set_threads(4).ok();

    let mut record = Record::new();  // REUSE single Record
    let mut stats = FilterResult::default();

    while bam.read(&mut record).is_ok() {
        stats.reads_processed += 1;

        if let Some(cb) = get_cb_tag(&record) {
            if whitelist.contains(&cb) {
                out.write(&record)?;
                stats.reads_kept += 1;
            } else {
                stats.reads_rejected += 1;
            }
        } else {
            stats.reads_rejected += 1;
        }
    }

    Ok(stats)
}
```

### Files to Modify

1. `sheriff-rs/src/bam_filter.rs`
   - Add new `filter_bam_streaming()` function
   - Deprecate `filter_bam_by_barcodes_parallel()` (memory hog)

2. `sheriff-rs/src/python.rs`
   - Add PyO3 wrapper for streaming version

---

## Phase 5: PyO3 Zero-Copy (OPTIONAL)

**Timeline:** 7-14 days
**Impact:** 10-30% reduction in FFI overhead
**Risk:** HIGH (breaking changes)

### The Problem

```rust
// Every string crosses the FFI boundary via copy
fn deduplicate_umis_py(umis: Vec<String>) -> usize {
    // Python → Rust: N string copies
    deduplicate_umis_rust(&umis)
}
```

### The Fix - Use Python Native Types

```rust
use pyo3::types::{PyList, PyString};

#[pyfunction]
fn deduplicate_umis_py_zero_copy(umis: &PyList) -> PyResult<usize> {
    let umi_vec: Vec<String> = umis.iter()
        .filter_map(|item| {
            item.downcast::<PyString>()
                .ok()
                .and_then(|s| s.to_str().ok())
                .map(|s| s.to_string())
        })
        .collect();
    Ok(deduplicate_umis_rust(&umi_vec))
}
```

### Consideration

- Actual FFI overhead may be negligible compared to algorithm work
- Benchmark before committing to this refactor
- May not be worth the breaking changes

---

## Testing Strategy

### Unit Tests (Rust)

```bash
cd sheriff-rs
cargo test --release
```

### Integration Tests (Python)

```bash
python validate_rust_correctness.py
```

### Performance Benchmarks

```bash
# Before each phase
python benchmark_comprehensive.py --dataset medium > before.log

# After each phase
python benchmark_comprehensive.py --dataset medium > after.log

# Compare
diff before.log after.log
```

### Memory Profiling

```bash
# Install heaptrack
cargo build --release
heaptrack ./target/release/sheriff-rs [args]

# Or use valgrind
valgrind --tool=massif ./target/release/sheriff-rs [args]
```

---

## Success Metrics

### Phase 1 (O(n³) Fix)
- 200 edits: 25x → 50-100x speedup
- Memory allocations: 50% reduction
- All correctness tests pass

### Phase 2 (Aligner Reuse)
- Edit distance calls: 0 allocations per call
- Overall edit clustering: 30-50x speedup

### Phase 3 (String Allocations)
- Heap allocations: 90% reduction in edit clustering
- Cache hit rate improvement

### Phase 4 (BAM Streaming)
- Memory usage: <100MB for any BAM size
- Throughput maintained at 3.7x

---

## Implementation Order

1. **Week 1-2: Phase 1**
   - HashSet replacement (QUICK WIN)
   - Validate correctness
   - Benchmark improvement

2. **Week 2-3: Phase 2**
   - Aligner reuse implementation
   - Capacity tuning
   - Memory profiling

3. **Week 3-4: Phase 3**
   - String allocation optimization
   - Byte-level operations
   - Final benchmarking

4. **Week 4-6: Phase 4**
   - BAM streaming implementation
   - Large file testing
   - Memory validation

5. **Week 6-8: Phase 5 (Optional)**
   - PyO3 optimization prototype
   - Benchmark FFI overhead
   - Decide if worth pursuing

---

## Deliverables

1. **Phase 1 Complete**
   - O(n²) edit clustering algorithm
   - Benchmark results showing improvement
   - Updated RUST_AUDIT_REPORT.md

2. **Phase 2-3 Complete**
   - Memory-optimized Rust code
   - Profiling data showing allocation reduction
   - Performance comparison graphs

3. **Phase 4 Complete**
   - Streaming BAM filter
   - Memory usage graphs
   - Documentation for new API

4. **Final**
   - Updated README with performance claims
   - Comprehensive benchmark suite
   - CI/CD pipeline for regression testing

---

## Risk Mitigation

1. **Correctness Regression**
   - Run `validate_rust_correctness.py` after every change
   - Compare Python vs Rust results on real datasets
   - Maintain backwards compatibility

2. **Performance Regression**
   - Benchmark before/after each change
   - Profile hotspots with `cargo flamegraph`
   - A/B test in production-like environment

3. **Breaking Changes**
   - Keep old API alongside new (deprecation warnings)
   - Version bump for major changes
   - Document migration path

---

## Questions to Resolve

1. **Phase 2:** Thread-local vs. parameter passing for Aligner?
2. **Phase 3:** Pre-compute all reversals vs. lazy computation?
3. **Phase 4:** Streaming vs. chunked parallel for BAM?
4. **Phase 5:** Is PyO3 overhead actually measurable?

---

## Next Step

**Start Phase 1 immediately** - it's the highest impact, lowest risk change. The Edit struct already has Hash implemented, so this is a 2-4 hour fix that could yield 10-100x improvement for large edit sets.

```bash
# Execute Phase 1
cd /iblm/netapp/home/jjaureguy/Sheriff/sheriff-rs
# Edit src/edit_clustering.rs lines 251-252, 264, 267, 362, 366, 372, 377, 380, 387-394
maturin develop --release --features python
python ../validate_rust_correctness.py
python ../benchmark_comprehensive.py --dataset medium
```
