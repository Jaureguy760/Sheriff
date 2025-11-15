# Sheriff Optimization Analysis: Bioinformatics Perspective

## Research-Driven Insights (2024-2025)

### Key Findings

**Numba vs Cython for Genomics:**
- Numba: 90%+ speed improvements, scales better with large datasets
- Cython: Better for <1000 elements, saturates at 100-150x speedup
- **Numba wins for genomics** due to scaling with dataset size
- **Critical limitation**: Numba doesn't work with pandas/scipy (common in bioinformatics)

**scRNA-seq Performance Bottlenecks:**
- **I/O dominates**: BAM file reading is typically the bottleneck, not computation
- STARsolo is 10x faster than Cellranger primarily due to I/O optimization
- kallisto/bustools uses constant memory with streaming (no intermediate BAM files)
- **Memory requirements**: 32-64 GB recommended for alignment

---

## Where is Sheriff's REAL Bottleneck?

### Current Architecture Analysis

```python
# Sheriff's processing pipeline:
1. Iterate ALL reads in BAM (I/O intensive) ← LIKELY BOTTLENECK
2. Filter by cell barcode (Python dict lookup) ← Fast
3. Match k-mer barcodes (some Python, some optimized) ← Medium
4. FAISS nearest neighbor search ← Already optimized (C++)
5. UMI deduplication (mixed: Python + Numba) ← Partially optimized
6. Gene counting (parallel, pysam) ← Reasonably optimized
```

### Profiling Hypothesis

Based on similar tools and Sheriff's architecture:

| Component | Estimated % Time | Optimization Status |
|-----------|-----------------|-------------------|
| **BAM iteration** | **40-60%** | ⚠️ NOT OPTIMIZED |
| K-mer matching | 15-25% | ⚠️ Partially optimized |
| FAISS indexing | 5-10% | ✅ Optimized (C++) |
| UMI deduplication | 10-20% | ✅/⚠️ Mixed (Numba available) |
| Gene counting | 10-15% | ✅ Parallelized |

**Conclusion:** K-mer hashing optimization may only improve 15-25% of runtime!

---

## Alternative Optimization Strategies

### Strategy 1: I/O Optimization (Highest Impact)

**Problem:** Sheriff reads entire BAM sequentially

**Current (count_t7.py:318):**
```python
for i, read in enumerate(bam):  # Reads ALL reads
    cell_barcode = read.get_tag('CB')
    if cell_barcode not in cell_barcodes:
        continue  # Wastes 50-90% of iterations!
```

**Optimization Options:**

#### Option A: Pre-filter BAM by Cell Barcode
```bash
# One-time preprocessing
samtools view -h input.bam | \
    awk -F'\t' 'BEGIN{while(getline<"whitelist.txt")w[$1]}/^@/{print}/CB:Z:/{split($0,a,"CB:Z:");split(a[2],b,"\t");if(b[1] in w)print}' | \
    samtools view -b > filtered.bam
samtools index filtered.bam
```

**Impact:** 50-90% reduction in reads processed
**Trade-off:** Extra preprocessing step, more disk space

#### Option B: Streaming with Indexed Fetch
```python
# Use BAM cell barcode index if available (requires CB:Z tag indexing)
# Or process by genomic region in parallel
for chromosome in chromosomes:
    for read in bam.fetch(chromosome):
        # Process in parallel by region
```

**Impact:** Better parallelization, reduced memory
**Trade-off:** More complex code

#### Option C: Use Rust/C++ for BAM Parsing
```python
# Use rust-bio or HTSlib directly
# Similar to what kallisto/bustools does
```

**Impact:** 5-10x faster BAM parsing
**Trade-off:** Requires Rust/C++ integration

---

### Strategy 2: Vectorization over Pure Speed

**Insight:** NumPy vectorization often beats Cython/Numba for array operations

**Current k-mer matching (count_t7.py:124-128):**
```python
freq_array = np.zeros((4 ** k), dtype=np.uint8)
freq_array[[bc_kmer_matcher.kmer_to_num(indel_seq[i: i + k])
            for i in range(len(indel_seq) - k + 1)]] += 1
```

**Vectorized alternative:**
```python
# Use numpy rolling window
from numpy.lib.stride_tricks import sliding_window_view

def kmer_hash_vectorized(sequence, k):
    # Convert sequence to numeric array
    seq_numeric = np.array([BASE_MAP[b] for b in sequence])

    # Create sliding windows
    windows = sliding_window_view(seq_numeric, k)

    # Vectorized hash computation
    hash_vals = (windows @ (4 ** np.arange(k)[::-1])).astype(np.uint32)

    # Bincount is faster than manual indexing
    return np.bincount(hash_vals, minlength=4**k).astype(np.uint8)
```

**Impact:** 2-5x faster than current Python loops
**Benefit:** Pure NumPy, no compilation needed

---

### Strategy 3: Cython for String-Heavy Operations

**When to use Cython over Numba:**
- String manipulation (Numba weak here)
- Integration with C libraries (pysam internals)
- Memory views for efficient array access
- Need exact control over memory layout

**Where Cython makes sense in Sheriff:**

#### Optimized K-mer Hashing (Cython)
```cython
# sheriff/kmer_cython.pyx
cimport numpy as np
import numpy as np

cdef dict BASE_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

cpdef np.ndarray[np.uint32_t, ndim=1] kmer_hash_cython(str sequence, int k):
    cdef int i, j, hash_val
    cdef int n = len(sequence) - k + 1
    cdef np.ndarray[np.uint32_t, ndim=1] hashes = np.zeros(n, dtype=np.uint32)

    for i in range(n):
        hash_val = 0
        for j in range(k):
            hash_val = hash_val * 4 + BASE_MAP[sequence[i + j]]
        hashes[i] = hash_val

    return hashes
```

**Build setup (setup.py addition):**
```python
from Cython.Build import cythonize
import numpy

ext_modules = cythonize("sheriff/kmer_cython.pyx")

setup(
    ...,
    ext_modules=ext_modules,
    include_dirs=[numpy.get_include()],
)
```

**Pros:**
- 5-10x faster for string operations
- Standard in bioinformatics (pysam, HTSeq use it)
- Better string handling than Numba

**Cons:**
- Requires C compiler (complicates installation)
- More complex build process
- conda/pip wheels needed for distribution

---

### Strategy 4: Hybrid Approach (Recommended)

**Combine best of each tool:**

```python
# Fast path: NumPy vectorization (no dependencies)
def kmer_hash_numpy(sequence, k):
    # Pure NumPy implementation
    ...

# Faster path: Cython (optional, if C compiler available)
try:
    from .kmer_cython import kmer_hash_cython as kmer_hash
except ImportError:
    kmer_hash = kmer_hash_numpy  # Fallback

# Existing: Numba for UMI deduplication (works well)
@jit(nopython=True)
def deduplicate_umis_fast(...):
    # Keep existing Numba code
```

**Benefits:**
- Works without C compiler (NumPy fallback)
- Fast when Cython available
- Maintains Numba for what it's good at (UMI dedup)

---

## Rust Integration (Future Consideration)

**Why genomics tools are moving to Rust:**
- minimap2, salmon: C/Rust for core algorithms
- Safety + Performance
- Better parallelism than Python

**Rust for Sheriff (long-term):**
```rust
// BAM iteration + k-mer matching in Rust
// Use rust-bio crate

use rust_htslib::{bam, bam::Read};

pub fn process_bam(bam_path: &str) -> Vec<Edit> {
    let mut bam = bam::Reader::from_path(bam_path).unwrap();
    // Ultra-fast parallel BAM processing
    ...
}
```

**PyO3 binding:**
```python
# Python calls Rust
from sheriff_rs import process_bam_fast

edits = process_bam_fast("data.bam", cell_barcodes)
```

**Trade-off:**
- 10-100x speedup for BAM processing
- Requires Rust toolchain
- More complex distribution

---

## Recommendation: Phased Approach

### Phase 1: Low-Hanging Fruit (This Week)
**No new dependencies, pure Python/NumPy optimizations**

1. ✅ Vectorize k-mer hashing with NumPy
2. ✅ Fix namedtuple recreation
3. ✅ Array reuse for freq_array
4. ✅ Ensure fast UMI dedup everywhere
5. ✅ Add simple BAM caching/filtering

**Expected:** 1.5-2x speedup
**Risk:** Very low
**Time:** 2-3 hours

### Phase 2: Optional Cython (Next Sprint)
**If users have C compiler**

1. Cython k-mer hashing
2. Cython BAM tag extraction
3. Provide pip wheels for common platforms

**Expected:** 2-3x speedup (combined with Phase 1)
**Risk:** Low (has NumPy fallback)
**Time:** 1-2 days

### Phase 3: I/O Optimization (When Needed)
**For users processing 10k+ cell libraries**

1. BAM pre-filtering utilities
2. Parallel region-based processing
3. Consider Rust for hot paths

**Expected:** 5-10x speedup for large datasets
**Risk:** Medium (more complex)
**Time:** 1 week

---

## Profiling Action Items

**Before optimizing, let's MEASURE:**

```bash
# 1. Profile with cProfile (find real bottleneck)
python -m cProfile -o sheriff.prof -m sheriff <args>
snakeviz sheriff.prof

# 2. Line-by-line profiling
pip install line_profiler
kernprof -l -v sheriff/count_t7.py

# 3. I/O profiling
strace -c sheriff <args>  # See system calls
```

**Expected findings:**
- 40-60% time in pysam read iteration
- 15-25% in k-mer matching
- 10-20% in UMI deduplication

---

## Decision Matrix

| Optimization | Speedup | Complexity | Dependencies | Maintenance |
|--------------|---------|------------|--------------|-------------|
| NumPy vectorization | 2-3x | Low | None | Low |
| Numba (current) | 10-100x | Low | Numba | Low |
| **Cython** | **5-10x** | **Medium** | **C compiler** | **Medium** |
| Rust | 10-100x | High | Rust | High |
| BAM pre-filtering | 2-5x | Low | None | Low |
| Parallel I/O | 3-8x | Medium | None | Medium |

**Verdict:**
1. Start with NumPy vectorization (easy win)
2. Add Cython as optional enhancement (conda builds)
3. Profile before pursuing Rust

---

## Questions for You

1. **What dataset sizes do users typically process?**
   - <1k cells: Cython overkill
   - 10k+ cells: Cython/Rust worth it

2. **Installation constraints?**
   - Conda-only: Can ship Cython pre-compiled
   - pip on varied systems: NumPy safer

3. **Performance target?**
   - Current: ~30s for example data
   - Target: <15s? <5s?

4. **Most important: Speed or Memory?**
   - Speed: Cython + parallel I/O
   - Memory: Streaming + constant memory algorithms

Let me know your priorities and I'll implement accordingly!
