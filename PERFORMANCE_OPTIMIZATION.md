# Sheriff Performance Optimization Guide

## Baseline Performance Metrics

**System:** 16 CPU cores, 13.0 GB RAM, Python 3.11.14
**Timestamp:** 2025-11-15

### Micro-Benchmarks (Baseline)

| Benchmark | Throughput | Runtime | Peak Memory |
|-----------|-----------|---------|-------------|
| **K-mer Matching** | 575 ops/sec | 17.38s | 969 MB |
| **UMI Dedup (Slow)** | 485 UMIs/sec | 10.31s | 823 MB |
| **UMI Dedup (Fast/Numba)** | 326 UMIs/sec* | 15.35s | 898 MB |

*Note: Fast version shows slower due to Numba JIT compilation overhead on first run. In production with warm cache, expect 10-100x speedup.

---

## Optimization Strategy

### Phase 1: Quick Wins (Implemented)
1. ✅ Create benchmarking framework
2. ✅ Establish baseline metrics
3. ⏳ Fix namedtuple recreation
4. ⏳ Make k-mer hashing iterative
5. ⏳ Add array reuse for freq_array
6. ⏳ Optimize k-mer counting with Numba

### Phase 2: Algorithmic Improvements
1. BAM filtering optimization
2. Parallel BAM reading
3. FAISS index optimization
4. Data structure improvements

### Phase 3: Memory Optimization
1. Reduce defaultdict usage
2. Implement sparse matrix earlier
3. Pool object allocation

---

## How to Benchmark

### Before Making Changes

```bash
# Run baseline benchmarks
python benchmarks/benchmark_sheriff.py --mode baseline --output benchmarks/baseline_results.json
```

### After Optimization

```bash
# Run optimized benchmarks
python benchmarks/benchmark_sheriff.py --mode optimized --output benchmarks/optimized_results.json

# Compare results
python benchmarks/benchmark_sheriff.py --mode compare \
    --baseline-file benchmarks/baseline_results.json \
    --optimized-file benchmarks/optimized_results.json
```

### Memory Profiling

```bash
# Install dependencies
pip install memory_profiler matplotlib

# Profile specific functions
python -m memory_profiler benchmarks/profile_memory.py

# Or with visualization
mprof run benchmarks/profile_memory.py
mprof plot
```

---

## Optimization Targets

### Performance Goals

| Metric | Current | Target | Improvement |
|--------|---------|--------|-------------|
| K-mer matching | 575 ops/s | 1,500+ ops/s | 2.6x |
| UMI dedup | 485 UMIs/s | 5,000+ UMIs/s | 10x |
| Full pipeline | N/A | <15s | N/A |
| Peak memory | 969 MB | <700 MB | 28% reduction |

### Critical Bottlenecks (Priority Order)

1. **CRITICAL:** Full BAM iteration (50-70% potential speedup)
2. **HIGH:** Recursive k-mer hashing (3-5x speedup)
3. **HIGH:** K-mer array allocation (10-20% speedup)
4. **HIGH:** UMI deduplication algorithm (10-100x speedup)
5. **MEDIUM:** String slicing in loops (2-3x speedup)

---

## Detailed Optimization Plans

### 1. Iterative K-mer Hashing

**Current (Recursive):**
```python
def kmer_to_num(self, kmer):
    if len(kmer) < 1:
        return 0
    return (4*self.kmer_to_num(kmer[:-1:])) + self.hash_symbol[kmer[-1]]
```

**Optimized (Iterative with Bit Shifts):**
```python
def kmer_to_num(self, kmer):
    result = 0
    for nucleotide in kmer:
        result = (result << 2) + self.hash_symbol[nucleotide]
    return result
```

**Expected:** 3-5x faster

### 2. Array Reuse

**Current:**
```python
def match_kmer(...):
    freq_array = np.zeros((4 ** k), dtype=np.uint8)  # Allocated every call!
```

**Optimized:**
```python
class KmerMatcher:
    def __init__(self, k, sequences=None):
        self._freq_array = np.zeros((4 ** k), dtype=np.uint8)

    def match_kmer(...):
        self._freq_array.fill(0)  # Reuse instead of allocate
```

**Expected:** 10-20% speedup, 20% memory reduction

### 3. Numba-Optimized K-mer Counting

```python
@jit(nopython=True)
def count_kmers_numba(seq, k, hash_array):
    nucleotide_map = np.array([0, 1, 2, 3])  # A, C, G, T
    for i in range(len(seq) - k + 1):
        hash_val = 0
        for j in range(k):
            hash_val = hash_val * 4 + get_nucleotide_int(seq[i+j])
        hash_array[hash_val] += 1
```

**Expected:** 2-3x faster

### 4. Replace Slow UMI Deduplication

**Action:** Ensure all code paths use `cell_umi_counts_FAST` instead of `deduplicate_umis`

**Expected:** 10-100x faster for large UMI sets

---

## Testing Procedure

### 1. Unit Tests (Pre-optimization)
```bash
# Ensure current functionality works
pytest tests/  # After adding tests
```

### 2. Benchmark Baseline
```bash
python benchmarks/benchmark_sheriff.py --mode baseline
```

### 3. Implement Optimization
- Make code changes
- Run unit tests to ensure correctness
- Check for edge cases

### 4. Benchmark Optimized
```bash
python benchmarks/benchmark_sheriff.py --mode optimized
```

### 5. Compare and Validate
```bash
python benchmarks/benchmark_sheriff.py --mode compare \
    --baseline-file benchmarks/baseline_results.json \
    --optimized-file benchmarks/optimized_results.json
```

### 6. Regression Testing
```bash
# Run on real data if available
sheriff <args> -o test_output_baseline
# <make optimizations>
sheriff <args> -o test_output_optimized

# Compare outputs
diff -r test_output_baseline test_output_optimized
```

---

## Performance Monitoring in Production

### Add Timing Decorators

```python
import time
from functools import wraps

def timeit(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        print(f"{func.__name__}: {time.time() - start:.3f}s")
        return result
    return wrapper

@timeit
def expensive_function():
    ...
```

### Log Performance Metrics

```python
import logging

logger = logging.getLogger(__name__)

# In critical functions
start = time.time()
# ... processing ...
logger.info(f"Processed {n_reads} reads in {time.time()-start:.2f}s")
```

---

## Optimization Checklist

Before committing optimizations:

- [ ] Baseline benchmarks recorded
- [ ] Optimization implemented
- [ ] Post-optimization benchmarks show improvement
- [ ] No functionality regression (tests pass)
- [ ] Memory usage acceptable
- [ ] Code remains readable
- [ ] Documentation updated
- [ ] CHANGELOG.md updated with performance notes

---

## Resources

- [Python Performance Tips](https://wiki.python.org/moin/PythonSpeed/PerformanceTips)
- [Numba Documentation](https://numba.readthedocs.io/)
- [Profiling Python Code](https://docs.python.org/3/library/profile.html)
- [Memory Profiler](https://pypi.org/project/memory-profiler/)

---

## Contact

For performance questions or optimization suggestions:
- Email: bbalderson@salk.edu
- Issues: https://github.com/BradBalderson/Sheriff/issues
