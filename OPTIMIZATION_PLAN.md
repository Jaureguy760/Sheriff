# Rust Edit Clustering Optimization Plan
**Status**: Baseline validated and correct ✅
**Date**: 2025-11-17
**Branch**: `claude/rust-optimization-upgrade-01WHgvkxgGsPTW6SMJBeKZDa`

---

## 🎯 Objective
Optimize edit clustering performance while **maintaining 100% correctness** on all validation tests.

**CRITICAL RULE**: Speed without correctness is worthless. Validate after EVERY change.

---

## ✅ Current Status: Baseline Validated

### Bugs Fixed
1. **Missing chromosome check** - Edits on different chromosomes were being compared
2. **Missing final edits** - Edits never compared weren't added to output

### Validation Results
```
✅ All 6 comprehensive correctness tests PASS
✅ CI validation with 352,535 real reads PASS
✅ All edge cases handled correctly
```

**Baseline correctness is GUARANTEED** - this is our ground truth.

---

## 📊 Step 1: Performance Baseline (IN PROGRESS)

### Questions to Answer:
1. What is the **actual runtime** of edit clustering on 352k reads?
2. Where exactly is time spent? (alignment vs string ops vs logic)
3. How many edits are typically processed? (10? 50? 500?)
4. What's the distribution of chromosomes/positions?

### Benchmarking Strategy:
```python
# Test with REAL genomic data from test_200kb.bam
import time
import sheriff_rs

# Extract real edits from pipeline
edits = extract_edits_from_bam('test_data/test_200kb.bam')

# Benchmark multiple runs for statistical significance
runs = 10
times = []
for _ in range(runs):
    start = time.perf_counter()
    result = sheriff_rs.get_longest_edits_rust(edits)
    end = time.perf_counter()
    times.append(end - start)

print(f"Mean: {sum(times)/len(times)*1000:.2f}ms")
print(f"Std: {statistics.stdev(times)*1000:.2f}ms")
print(f"Number of edits: {len(edits)} → {len(result)}")
```

### Profiling with py-spy:
```bash
# Profile actual Sheriff run on real data
py-spy record -o profile.svg -- python -c "
from sheriff.pipeline import process_bam
process_bam('test_data/test_200kb.bam', ...)
"
```

---

## 🔍 Step 2: Bottleneck Analysis

### Algorithm Complexity Analysis

**Current algorithm**: O(n²) pairwise comparisons
- For n edits: n*(n-1)/2 comparisons
- Each comparison involves:
  1. ✅ **Cheap**: Chromosome/position/orientation check (O(1))
  2. ✅ **Cheap**: String slicing (O(k) where k=ref_len, small)
  3. ❌ **EXPENSIVE**: `bio_edit_distance()` - Smith-Waterman alignment
  4. ❌ **EXPENSIVE**: Homopolymer detection (O(m) where m=seq_len)
  5. ❌ **EXPENSIVE**: Homopolymer collapse + re-alignment
  6. ❌ **EXPENSIVE**: 3' end check (reverse + alignment)

**Key insight**: If we can skip expensive operations for dissimilar sequences, we win big!

### Potential Bottlenecks (by likelihood):
1. **Smith-Waterman alignment** (bio_edit_distance) - O(m*n) per call
2. **Repeated string operations** - reverse, substring, homopolymer collapse
3. **Hash operations** - AHashSet lookups (though already optimized)
4. **Cloning edits** - inserting into HashSet requires clone()

---

## 🎨 Step 3: Optimization Candidates

### LOW RISK (✅ Safe, measurable impact)

#### A. Pre-compute string slices
**Idea**: Cache `edit.alt_seq[..len-reflen]` for forward edits
**Risk**: Low - we already do this for reverse
**Expected gain**: 5-10% (eliminate repeated slicing)
**Validation**: Should produce identical output

```rust
// Pre-compute forward sequences too
let mut forward_cache: AHashMap<usize, &str> = AHashMap::new();
for (idx, edit) in edits.iter().enumerate() {
    if edit.forward {
        let reflen = edit.ref_seq.len();
        forward_cache.insert(idx, &edit.alt_seq[..edit.alt_seq.len() - reflen]);
    }
}
```

#### B. Length-based pre-filter
**Idea**: Skip alignment if length difference > threshold
**Risk**: Low - biologically justified
**Expected gain**: 20-30% (skip many comparisons)
**Validation**: Check that similar-length edits still cluster correctly

```rust
// If sequences differ by >10bp, they're definitely different
let len_diff = (edit_1_seq.len() as i64 - edit_2_seq.len() as i64).abs();
if len_diff > 10 {
    // Add both, continue - they're too different
    continue;
}
```

#### C. Early exit on position distance
**Idea**: Edits >50bp apart are likely different
**Risk**: Low - biologically justified
**Expected gain**: 10-20% for datasets with sparse edits
**Validation**: Ensure nearby edits still cluster

```rust
// If ref_pos differs by >50, they're likely different edits
if (edit_1.ref_pos - edit_2.ref_pos).abs() > 50 {
    continue;
}
```

---

### MEDIUM RISK (⚠️ Requires careful validation)

#### D. Hamming distance pre-filter
**Idea**: Quick Hamming distance check before expensive alignment
**Risk**: MEDIUM - We tried this before and it BROKE CORRECTNESS
**Expected gain**: 30-50% if done correctly
**Validation**: MUST test on all edge cases, homopolymers, etc.

```rust
// Only for same-length sequences
if edit_1_seq.len() == edit_2_seq.len() {
    let hamming = edit_1_seq.bytes().zip(edit_2_seq.bytes())
        .filter(|(a, b)| a != b).count();

    // If Hamming distance > 5, definitely different
    if hamming > 5 {
        continue;
    }
}
```

**WARNING**: Last time we tried this, we broke correctness!
**Lesson learned**: Must validate against ALL test cases, not just CI baseline.

#### E. Parallel processing
**Idea**: Use rayon to parallelize the O(n²) loop
**Risk**: MEDIUM - thread safety, HashSet concurrent access
**Expected gain**: 2-4x on multi-core
**Validation**: Ensure output is deterministic

---

### HIGH RISK (🚨 Could break correctness)

#### F. Change clustering threshold
**Idea**: Tighten threshold from ≤2 to ≤1
**Risk**: HIGH - changes biological interpretation
**Expected gain**: Unknown - might cluster less, reducing downstream processing
**Validation**: Would need biologist approval

❌ **DO NOT ATTEMPT** without explicit approval

#### G. Different alignment algorithm
**Idea**: Replace Smith-Waterman with faster approximate algorithm
**Risk**: HIGH - fundamental change to accuracy
**Expected gain**: 50-80% potentially
**Validation**: Would need extensive validation on real data

❌ **DO NOT ATTEMPT** without extensive testing

---

## 🧪 Step 4: Testing Strategy

### For EACH optimization attempt:

#### 1. Pre-optimization checklist:
- [ ] Document exact change being made
- [ ] Predict expected speedup
- [ ] Identify what could go wrong
- [ ] Save baseline performance numbers

#### 2. Implementation:
- [ ] Make ONE change at a time
- [ ] Add clear comments explaining optimization
- [ ] Keep old code in comments for reference

#### 3. Validation (MUST ALL PASS):
```bash
# Correctness tests
python validate_rust_correctness.py  # Must show 6/6 PASS

# CI validation with real data
python test_data/ci_validation.py    # Must pass all checks

# Rust unit tests
cd sheriff-rs && cargo test          # Must pass
```

#### 4. Performance measurement:
```python
# Benchmark before/after with same data
python benchmark_edit_clustering.py
# Expected output:
#   Baseline: X.XX ms
#   Optimized: Y.YY ms
#   Speedup: Z.Zx
```

#### 5. Regression testing:
- [ ] Test on multiple datasets (not just test_200kb.bam)
- [ ] Test edge cases (1 edit, 1000 edits, identical sequences)
- [ ] Test worst case (all edits similar, requiring full alignment)

#### 6. Commit strategy:
```bash
# Good commit message format:
git commit -m "Optimize: [specific change]

Performance:
- Before: X.X ms
- After: Y.Y ms
- Speedup: Z.Zx

Validation: All tests pass
- validate_rust_correctness.py: 6/6 PASS
- ci_validation.py: PASS
- cargo test: PASS
"
```

---

## 📈 Step 5: Incremental Optimization Process

### Workflow:
1. **Profile** → Identify bottleneck
2. **Plan** → Choose ONE optimization from above
3. **Predict** → Estimate impact and risk
4. **Implement** → Make minimal change
5. **Validate** → Run ALL tests
6. **Benchmark** → Measure actual speedup
7. **Document** → Record results
8. **Commit** → Save working state
9. **Repeat** → Go to step 1

### Red flags that mean STOP:
- ❌ Any validation test fails
- ❌ Output differs from baseline on any test
- ❌ Speedup is suspiciously large (>10x from one change)
- ❌ Code becomes unclear or hard to understand
- ❌ "I think this should work..." - NO, prove it!

### Green flags that mean continue:
- ✅ All validation tests pass
- ✅ Output matches baseline exactly
- ✅ Speedup is measurable and explained
- ✅ Code is clear with good comments
- ✅ Performance gain is reproducible

---

## 🎯 Success Criteria

### Minimum viable optimization:
- [ ] 2-5x speedup on edit clustering
- [ ] 100% correctness maintained
- [ ] All tests pass
- [ ] Code remains readable

### Stretch goals:
- [ ] 10x speedup on edit clustering
- [ ] Comprehensive benchmark suite
- [ ] Performance regression tests in CI/CD
- [ ] Profiling results documented

---

## 📝 Next Actions

1. ✅ **DONE**: Fix baseline bugs
2. 🔄 **IN PROGRESS**: Benchmark baseline performance
3. ⏭️ **NEXT**: Profile with py-spy to find hotspots
4. ⏭️ **THEN**: Implement LOW RISK optimizations first
5. ⏭️ **FINALLY**: Consider MEDIUM RISK if needed

---

## 💭 Lessons Learned

From previous optimization attempt:

> **User feedback**: "this sounds suspicious as fuck? did you validate against our test dataset?"

**Lesson**: NEVER claim speedup without validation. The user will catch it.

**Approach**: Measure first, optimize second, validate always.

**Mantra**: Speed without correctness is worthless. 🎯
