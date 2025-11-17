# Sheriff - Claude Cloud Development Setup

Complete guide for setting up Sheriff development in Claude Cloud environments.

## Quick Start (5-10 minutes)

```bash
# 1. Clone repository
git clone https://github.com/BradBalderson/Sheriff.git
cd Sheriff

# 2. Download self-contained test data (23MB - includes reference!)
wget -O- -q https://github.com/BradBalderson/Sheriff/releases/latest/download/test_data_chr19.tar.gz | tar xvzf -
gunzip test_data/chr19.fa.gz test_data/chr19.gtf.gz

# 3. Set up Python environment
python3 -m venv venv
source venv/bin/activate
pip install -e .
pip install maturin pytest

# 4. Build Rust acceleration (if cargo available)
cd sheriff-rs && maturin develop --release && cd ..

# 5. Validate everything works
python test_data/ci_validation.py          # Tests Rust functions
python test_data/validate_chr19_test.py    # Tests data integrity
```

## What Claude Cloud Gets

### Test Data (23MB total)
- **BAM file**: 42k real reads on chr19
- **Reference genome**: GRCh38 chr19 (58MB uncompressed)
- **GTF annotations**: Ensembl 110 gene annotations
- **Cell barcodes**: 4,680 barcodes from real dataset
- **Edit sites**: 4 known CRISPR edit positions
- **MD5 checksums**: Deterministic validation

### Validated Rust Functions
- K-mer counting with checksum verification
- UMI deduplication (graph-based clustering)
- Edit clustering (O(n²) optimized algorithm)
- Cell UMI counting
- BAM filtering (parallel processing)

### Full Development Capabilities
- Edit Python code → test immediately
- Edit Rust code → `maturin develop --release` → test
- Run benchmarks
- Compare Python vs Rust implementations
- Profile performance bottlenecks

## Dev Environment Structure

```
Sheriff/
├── venv/                      # Python virtual environment
├── sheriff/                   # Python source code
│   ├── count_t7.py           # Main pipeline (edit here)
│   ├── helpers.py            # Core utilities
│   └── rust_accelerated.py   # Rust wrappers
├── sheriff-rs/               # Rust acceleration
│   ├── Cargo.toml            # Dependencies
│   └── src/
│       ├── edit_clustering.rs  # Edit deduplication
│       ├── umi.rs              # UMI deduplication
│       └── python.rs           # PyO3 bindings
├── test_data/                # Self-contained test dataset
│   ├── test_chr19.bam        # 42k reads
│   ├── chr19.fa              # Reference genome (58MB)
│   ├── chr19.gtf             # Gene annotations
│   ├── ci_validation.py      # Quick validation
│   └── expected_checksums.json  # MD5 hashes
└── sanity_check_comparison.py  # Python vs Rust benchmark
```

## Development Workflows

### 1. Quick Code Validation (1 second)
```bash
python test_data/ci_validation.py
```
Tests Rust functions with MD5 checksums - catches regressions instantly.

### 2. Data Integrity Check (2 seconds)
```bash
python test_data/validate_chr19_test.py
```
Validates BAM, reference, GTF, chromosome mapping all work together.

### 3. Full Pipeline Test (30-60 seconds)
```python
# Run Sheriff on single-chromosome test data
from sheriff import count_t7

# You can now run the actual pipeline!
# (Implementation depends on your specific needs)
```

### 4. Rust Development Cycle
```bash
# 1. Edit Rust code
vim sheriff-rs/src/edit_clustering.rs

# 2. Rebuild (30 seconds)
cd sheriff-rs && maturin develop --release && cd ..

# 3. Test
python test_data/ci_validation.py
```

## What You Can Work On

### Python Code (instant iteration)
- `sheriff/count_t7.py` - Main pipeline logic
- `sheriff/helpers.py` - Cell/UMI processing
- `sheriff/rust_accelerated.py` - Rust wrappers

### Rust Code (needs rebuild)
- `sheriff-rs/src/edit_clustering.rs` - Edit deduplication (current bottleneck)
- `sheriff-rs/src/umi.rs` - UMI deduplication
- `sheriff-rs/src/kmer.rs` - K-mer operations
- `sheriff-rs/src/python.rs` - PyO3 bindings

### Optimization Opportunities
1. **mRNA counting** - Currently Python, 80% of runtime
2. **Edit clustering** - Rust but O(n²), can be O(n log n)
3. **Parallel processing** - Many functions exist but unused
4. **Memory optimization** - AHash enabled but more possible

## Testing Without Full Pipeline

Even without running the complete Sheriff pipeline, you can:

```python
# Test individual Rust functions
import sheriff_rs

# K-mer counting
kmers = sheriff_rs.count_kmers_rust("ATCGATCGATCG", 4)

# UMI deduplication
unique = sheriff_rs.deduplicate_umis_py(["ATCG", "ATCC", "TTTT"])

# Edit clustering
edits = [("chr1", 1000, "A", "ATCG", True, [1,2])]
longest = sheriff_rs.get_longest_edits_rust(edits)

# BAM filtering (with test data)
import pysam
with pysam.AlignmentFile("test_data/test_chr19.bam", "rb") as bam:
    for read in bam.fetch():
        # Process reads...
        pass
```

## Cloud Environment Requirements

### Minimal (Python only)
- Python 3.10+
- pip/venv
- ~200MB disk (Sheriff + venv + test data)
- ~500MB RAM

### Full (Python + Rust)
- Python 3.10+
- Rust toolchain (cargo, rustc)
- ~500MB disk (+ Rust toolchain)
- ~1GB RAM (for Rust compilation)

### Install Rust (if not available)
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source ~/.cargo/env
```

## Checksum Validation

All test data is verified by MD5 checksums:

```json
{
  "test_chr19.bam": "4feec0ac4b8707fcd64744f4386be1a0",
  "chr19.fa.gz": "ca4f1a9927913fe5a62585a492d13c1a",
  "chr19.gtf.gz": "a77168a66e7fadf40bcfa171dbd6c537",
  "barcodes_chr19.txt": "ad89c4d2b2e5c78e0e213f71b693f641",
  "edit_sites_chr19.txt": "8de1610a3f293815e9281275acda7426"
}
```

## Syncing Back to HPC

After making changes in Claude Cloud:

```bash
# On Claude Cloud
git add -A
git commit -m "Optimizations from cloud dev"
git push origin rust-optimization-upgrade

# On HPC
cd /iblm/netapp/data3/jjaureguy/software/Sheriff
git pull origin rust-optimization-upgrade
# Rebuild Rust if changed
cd sheriff-rs && maturin develop --release && cd ..
# Run full benchmark on large data
python sanity_check_comparison.py
```

## Limitations in Cloud

1. **No 3GB reference** - Use single-chrom test (23MB) instead
2. **No HPC compute** - For final benchmarks, sync back to cluster
3. **No large datasets** - Test logic on small data, scale on HPC

## What Works Great in Cloud

1. **Code editing** - Python and Rust changes
2. **Logic validation** - MD5 checksums catch bugs
3. **Unit testing** - Fast iteration on algorithms
4. **Documentation** - Update README, docs
5. **Git operations** - Commit, branch, merge
6. **Small benchmarks** - Test on 42k reads vs 114M

## Example Development Session

```bash
# Setup (once)
git clone https://github.com/BradBalderson/Sheriff.git && cd Sheriff
wget -O- -q .../test_data_chr19.tar.gz | tar xvzf -
gunzip test_data/chr19.fa.gz test_data/chr19.gtf.gz
python3 -m venv venv && source venv/bin/activate
pip install -e . && pip install maturin
cd sheriff-rs && maturin develop --release && cd ..

# Daily workflow
source venv/bin/activate
python test_data/ci_validation.py  # Sanity check

# Make changes...
vim sheriff/count_t7.py

# Test
python test_data/ci_validation.py

# Commit
git add -A && git commit -m "Improved edit clustering"
git push
```

## Summary

**You can now:**
- ✅ Develop Sheriff in Claude Cloud
- ✅ Test with real biological data (42k reads, 4 edit sites)
- ✅ Validate Rust functions with MD5 checksums
- ✅ Edit both Python and Rust code
- ✅ Run fast iteration cycles (1-2 second tests)
- ✅ Sync changes back to HPC for full-scale testing

**Package size:** 23MB (vs 3GB+ for full reference)
**Setup time:** 5-10 minutes
**Test runtime:** 1-2 seconds for validation
