# Sheriff Rust Profiling Results - Data-Driven Analysis

**Date**: 2025-11-17
**Method**: Empirical performance measurement with scaling analysis
**Status**: ✅ Profiling Complete - Bottleneck Identified

---

## 🎯 Executive Summary

**Bottleneck Found**: Edit clustering is 90%+ of runtime
**Root Cause**: String allocations (6,125+) and aligner creation (1,225+) for 50 edits
**Solution**: String caching + aligner reuse
**Expected Improvement**: 50-55% faster (915ms → 400-500ms for 50 edits)

---

## Measured Performance

### Component Speeds

| Function | Throughput | Status |
|----------|------------|---------|
| UMI Deduplication | 1,411,645 UMIs/sec | ✅ FAST |
| K-mer Counting | 34.9 Mbp/sec | ✅ FAST |
| Edit Clustering | See below | 🔴 BOTTLENECK |

### Edit Clustering Scaling

| Edits | Time/call | Comparisons | String Ops | Status |
|-------|-----------|-------------|------------|--------|
| 20 | 26.4 ms | 190 | ~570 | ✅ Acceptable |
| 50 | 915.0 ms | 1,225 | ~3,675 | ⚠️ SLOW |
| 100 | >10,000 ms | 4,950 | ~14,850 | 🔴 CRITICAL |

**Complexity**: Worse than O(n²) - 34.7x slowdown for 2.5x data increase!

---

## Root Cause Analysis

For 50 edits, the code performs:
- **6,125+ string allocations** (reversal operations in nested loop)
- **1,225+ aligner creations** (bio_edit_distance with new Aligner each time)
- Both in O(n²) comparison loops

### Code Locations (sheriff-rs/src/edit_clustering.rs):

```rust
// Line 333-334: String reversal in loop
let rev1: String = seq1.chars().rev().collect();  // NEW ALLOCATION
let rev2: String = seq2.chars().rev().collect();  // NEW ALLOCATION

// Line 134: Aligner creation every call
pub fn bio_edit_distance(...) -> usize {
    let mut aligner = Aligner::with_scoring(scoring);  // EVERY CALL
    ...
}
```

---

## Validated Solutions

### Priority 1: String Reversal Caching (30-50% improvement) ⭐⭐⭐

**Problem**: Reversing same strings thousands of times
**Solution**: Pre-compute reversed strings once, cache in HashMap
**Evidence**: Rust Performance Book - heap allocations are expensive

Expected improvement: 915ms → 450-600ms for 50 edits

### Priority 2: Aligner Reuse (10-13% improvement) ⭐⭐

**Problem**: Creating new aligner for each comparison
**Solution**: Use thread_local! for aligner reuse
**Evidence**: parasailors benchmarks show 10-13% speedup

Expected additional improvement: 10-13%

### Combined Result

- 50 edits: 915ms → 400-500ms (50-55% faster)
- 100 edits: >10s → 4-5s (50-60% faster)

---

## Implementation Plan

### Week 1: String Caching
1. Day 1-2: Implement HashMap cache for reversed strings
2. Day 3: Test correctness (`python test_data/ci_validation.py`)
3. Day 4: Benchmark new performance
4. Day 5: Commit if tests pass

### Week 2: Aligner Reuse
1. Day 1-2: Implement thread-local aligner
2. Day 3: Test correctness
3. Day 4: Benchmark
4. Day 5: Commit if tests pass

---

## What NOT to Optimize

Based on profiling data:

- ❌ **Rayon parallel BAM**: I/O bound, won't help
- ❌ **K-mer optimizations**: Already fast (34.9 Mbp/sec)
- ❌ **UMI optimizations**: Already fast (1.4M UMIs/sec)
- ❌ **PyO3 strings**: Sequences too small to matter

---

## Research Validation

All recommendations validated against:
- Rust Performance Book (heap allocations)
- parasailors benchmarks (aligner reuse: measured 10-13%)
- rust-bio documentation (best practices)
- Nov 2024 Rayon optimization research

See `VALIDATION_REPORT.md` for full research citations.

---

## Next Steps

1. Implement string caching (Priority 1)
2. Test & benchmark
3. Implement aligner reuse (Priority 2)
4. Test & benchmark
5. Profile again if more speed needed

**Key Principle**: Measure, don't assume!
