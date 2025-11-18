# Sheriff-rs Python Bindings

Comprehensive Python bindings for the Sheriff-rs Rust library, exposing all Phase 1 optimizations for high-performance k-mer matching and UMI deduplication.

## Overview

The Python bindings provide significant performance improvements over pure Python implementations:
- **K-mer operations**: 4-14x faster
- **UMI deduplication**: 3-6x faster

All while maintaining a simple, Pythonic API that's compatible with the existing Sheriff codebase.

## Quick Start

### 1. Build and Install

```bash
cd sheriff-rs
./build_python.sh
```

Or manually:
```bash
pip install maturin
maturin develop --release --features python
```

### 2. Verify Installation

```bash
python3 test_python_bindings.py
```

### 3. Try the Demo

```bash
python3 examples/python_demo.py
```

## API Reference

### K-mer Functions

#### `kmer_to_num(kmer: str) -> int`

Convert a k-mer sequence to its numeric hash representation.

```python
import sheriff_rs

hash_val = sheriff_rs.kmer_to_num("ACGT")
print(hash_val)  # 27
```

**Parameters:**
- `kmer` (str): DNA sequence (A, C, G, T, case-insensitive)

**Returns:**
- `int`: Numeric hash of the k-mer

**Raises:**
- `ValueError`: If sequence contains invalid characters

---

#### `match_kmer(sequence: str, k: int, whitelist: List[int], output_hash: bool = True) -> List`

Match k-mers in a sequence against a whitelist.

```python
import sheriff_rs

sequence = "ACGTACGT"
whitelist = [sheriff_rs.kmer_to_num("ACGT")]

# Return hashes
matches = sheriff_rs.match_kmer(sequence, 4, whitelist, output_hash=True)
print(matches)  # [27, 27]

# Return k-mer strings
matches_str = sheriff_rs.match_kmer(sequence, 4, whitelist, output_hash=False)
print(matches_str)  # ['ACGT', 'ACGT']
```

**Parameters:**
- `sequence` (str): DNA sequence to scan
- `k` (int): K-mer length
- `whitelist` (List[int]): List of k-mer hashes to match
- `output_hash` (bool, optional): Return hashes if True, strings if False. Default: True

**Returns:**
- `List[int]` or `List[str]`: Matched k-mer hashes or sequences

**Raises:**
- `ValueError`: If k <= 0 or sequence contains invalid characters

---

### KmerCounter Class

Efficient k-mer frequency counter with array reuse pattern.

```python
import sheriff_rs

# Create counter for k=4
counter = sheriff_rs.KmerCounter(4)

# Count k-mers in a sequence
freqs = counter.count_kmers("ACGTACGT")

# Access frequency for a specific k-mer
acgt_hash = sheriff_rs.kmer_to_num("ACGT")
print(freqs[acgt_hash])  # 2

# Reuse counter for another sequence (efficient!)
freqs2 = counter.count_kmers("AAAATTTT")
```

**Methods:**

##### `__init__(k: int)`

Create a new KmerCounter.

**Parameters:**
- `k` (int): K-mer length (must be 1-16)

**Raises:**
- `ValueError`: If k is invalid

##### `count_kmers(sequence: str) -> List[int]`

Count k-mer frequencies in a sequence.

**Parameters:**
- `sequence` (str): DNA sequence to analyze

**Returns:**
- `List[int]`: Frequency array of length 4^k, indexed by k-mer hash

**Raises:**
- `ValueError`: If sequence contains invalid characters

##### `clear()`

Clear all frequency counts (automatically done by `count_kmers`).

**Properties:**
- `k` (int): K-mer length
- `array_size` (int): Size of frequency array (4^k)

---

### UMI Functions

#### `deduplicate_umis(umis: List[str], threshold: int) -> int`

Count unique UMI groups after deduplication.

```python
import sheriff_rs

umis = ["ATCGATCG", "ATCGATCC", "GCGCGCGC"]
unique_count = sheriff_rs.deduplicate_umis(umis, threshold=1)
print(unique_count)  # 2
```

**Parameters:**
- `umis` (List[str]): List of UMI sequences (all must be same length)
- `threshold` (int): Maximum Hamming distance for grouping (typically 1)

**Returns:**
- `int`: Number of unique UMI groups

**Raises:**
- `ValueError`: If UMIs have different lengths or contain invalid characters

---

#### `deduplicate_umis_detailed(umis: List[str], threshold: int) -> List[List[int]]`

Get detailed UMI grouping information.

```python
import sheriff_rs

umis = ["ATCGATCG", "ATCGATCC", "GCGCGCGC"]
groups = sheriff_rs.deduplicate_umis_detailed(umis, threshold=1)
print(groups)  # [[0, 1], [2]]
```

**Parameters:**
- `umis` (List[str]): List of UMI sequences
- `threshold` (int): Maximum Hamming distance for grouping

**Returns:**
- `List[List[int]]`: List of groups, each containing UMI indices

**Raises:**
- `ValueError`: If UMIs have different lengths or contain invalid characters

---

#### `hamming_distance(a: str, b: str) -> int`

Compute Hamming distance between two sequences.

```python
import sheriff_rs

dist = sheriff_rs.hamming_distance("ATCG", "ATGG")
print(dist)  # 1
```

**Parameters:**
- `a` (str): First sequence
- `b` (str): Second sequence

**Returns:**
- `int`: Number of positions where sequences differ

---

## Performance Tips

### 1. Reuse KmerCounter Objects

```python
# Good - reuses the internal array
counter = sheriff_rs.KmerCounter(k)
for sequence in sequences:
    freqs = counter.count_kmers(sequence)
    # process frequencies...

# Bad - creates new counter each time
for sequence in sequences:
    counter = sheriff_rs.KmerCounter(k)
    freqs = counter.count_kmers(sequence)
```

### 2. Use Hash Output for Matching

```python
# Faster - returns integers
matches = sheriff_rs.match_kmer(seq, k, whitelist, output_hash=True)

# Slower - converts to strings
matches = sheriff_rs.match_kmer(seq, k, whitelist, output_hash=False)
```

### 3. Choose the Right UMI Function

```python
# Fast - only returns count
count = sheriff_rs.deduplicate_umis(umis, threshold)

# Slower - returns full grouping details
groups = sheriff_rs.deduplicate_umis_detailed(umis, threshold)
```

## Integration with Existing Code

The bindings are designed to be drop-in replacements for Python implementations:

```python
# Old Python code:
class KmerAnalyzer:
    def kmer_to_num(self, kmer):
        # Python implementation...
        pass

# New Rust-accelerated code:
import sheriff_rs

class KmerAnalyzer:
    def kmer_to_num(self, kmer):
        return sheriff_rs.kmer_to_num(kmer)
```

## Examples

### Complete K-mer Analysis

```python
import sheriff_rs

# Create k-mer counter
counter = sheriff_rs.KmerCounter(k=6)

# Build whitelist of target k-mers
target_kmers = ["ACGTAC", "GTACGT", "TACGTA"]
whitelist = [sheriff_rs.kmer_to_num(kmer) for kmer in target_kmers]

# Process sequences
sequences = ["ACGTACGTACGTA", "GTACGTACGTAC"]

for seq in sequences:
    # Count all k-mers
    freqs = counter.count_kmers(seq)

    # Match against whitelist
    matches = sheriff_rs.match_kmer(seq, 6, whitelist, output_hash=True)

    print(f"Sequence: {seq}")
    print(f"  Total k-mers: {sum(freqs)}")
    print(f"  Whitelist matches: {len(matches)}")
```

### UMI Deduplication Pipeline

```python
import sheriff_rs

def process_umis(umi_list, threshold=1):
    """Process UMIs and return unique count and groupings."""

    # Get unique count (fast)
    unique_count = sheriff_rs.deduplicate_umis(umi_list, threshold)

    # Get detailed groupings if needed
    if need_details:
        groups = sheriff_rs.deduplicate_umis_detailed(umi_list, threshold)
        return unique_count, groups
    else:
        return unique_count, None

# Process
umis = ["ATCGATCG", "ATCGATCC", "GCGCGCGC", "GCGCGCGA"]
count, groups = process_umis(umis, threshold=1)

print(f"Unique UMIs: {count}")
for i, group in enumerate(groups):
    print(f"  Group {i+1}: {[umis[idx] for idx in group]}")
```

## Building for Distribution

### Create a Wheel

```bash
./build_python.sh wheel
```

The wheel will be in `target/wheels/` and can be distributed:

```bash
pip install target/wheels/sheriff_rs-*.whl
```

### Building for Multiple Python Versions

```bash
maturin build --release --features python \
    --interpreter python3.8 python3.9 python3.10 python3.11
```

## Troubleshooting

### Module Not Found

```python
>>> import sheriff_rs
ModuleNotFoundError: No module named 'sheriff_rs'
```

**Solution:** Build and install the module:
```bash
cd sheriff-rs
maturin develop --release --features python
```

### Functions Missing or Wrong Results

**Solution:** Rebuild with the python feature enabled:
```bash
maturin develop --release --features python --force
```

### Performance Not as Expected

1. Make sure you built with `--release` flag
2. Reuse `KmerCounter` objects
3. Use `output_hash=True` when possible
4. Profile your Python code to find actual bottlenecks

## Documentation

Get help in Python:

```python
import sheriff_rs

# Module documentation
help(sheriff_rs)

# Function documentation
help(sheriff_rs.kmer_to_num)
help(sheriff_rs.match_kmer)

# Class documentation
help(sheriff_rs.KmerCounter)
```

## Files

- `src/python.rs` - PyO3 bindings implementation
- `BUILD_PYTHON.md` - Detailed build instructions
- `build_python.sh` - Automated build script
- `test_python_bindings.py` - Test suite
- `examples/python_demo.py` - Comprehensive examples

## Support

For issues or questions:
1. Check the troubleshooting section
2. Review the build guide: `BUILD_PYTHON.md`
3. Run the test suite: `python3 test_python_bindings.py`
4. Try the demo: `python3 examples/python_demo.py`

## License

Same as the Sheriff project.
