# Sheriff - Claude AI Agent Setup

Quick setup for Claude AI agents to work on Sheriff Rust optimization.

## One-Command Setup

```bash
# Clone the optimization branch
git clone -b rust-optimization-upgrade https://github.com/BradBalderson/Sheriff.git
cd Sheriff

# Run setup (installs deps, downloads test data)
bash setup_cloud_dev.sh

# Activate environment
source venv/bin/activate

# Verify
python test_data/ci_validation.py
```

## If Test Data Download Fails

Upload `test_data_chr19.tar.gz` manually (23MB), then:
```bash
tar xzf test_data_chr19.tar.gz
gunzip test_data/chr19.fa.gz test_data/chr19.gtf.gz
```

The package contains:
- BAM file with 42k real reads
- GRCh38 chr19 reference genome
- Gene annotations
- 4,680 cell barcodes
- 4 known edit sites
- MD5 checksums for validation

## Key Files to Work On

### Python (instant iteration)
- `sheriff/count_t7.py` - Main pipeline logic
- `sheriff/helpers.py` - Cell/UMI processing
- `sheriff/rust_accelerated.py` - Rust wrappers

### Rust (rebuild after changes)
- `sheriff-rs/src/edit_clustering.rs` - Edit deduplication
- `sheriff-rs/src/umi.rs` - UMI deduplication
- `sheriff-rs/Cargo.toml` - Dependencies

## Current State (rust-optimization-upgrade branch)

### Completed
- PyO3 0.21, bio 2.3, rayon 1.11 updates
- AHash enabled for fast hashing
- target-cpu=native optimizations
- Python fallbacks removed (Rust-only)
- 1.45x speedup achieved

### Optimization Opportunities
1. **mRNA counting** - 80% of runtime, still in Python
2. **Edit clustering** - O(n²), can be O(n log n)
3. **Unused Rust functions** - `cell_umi_counts_py_parallel()` not wired up
4. **SIMD operations** - Not yet implemented

## Validation Commands

```bash
# Quick Rust function test (1 sec)
python test_data/ci_validation.py

# Data integrity check (2 sec)
python test_data/validate_chr19_test.py

# Check current git status
git status

# View recent commits
git log --oneline -10
```

## Development Cycle

```bash
# 1. Make changes
vim sheriff/count_t7.py

# 2. Test
python test_data/ci_validation.py

# 3. If editing Rust, rebuild
cd sheriff-rs && maturin develop --release && cd ..

# 4. Commit
git add -A
git commit -m "Description of changes"

# 5. Push to branch
git push origin rust-optimization-upgrade
```

## What Works Without HPC

✅ Edit Python/Rust code
✅ Test Rust functions with MD5 checksums
✅ Validate data integrity
✅ Small-scale benchmarking (42k reads)
✅ Git operations

## What Needs HPC

❌ Full benchmark (114M reads)
❌ Large-scale performance testing
❌ Multi-GB datasets

## MD5 Checksums

Verify test data integrity:
```
test_chr19.bam: 4feec0ac4b8707fcd64744f4386be1a0
chr19.fa.gz: ca4f1a9927913fe5a62585a492d13c1a
chr19.gtf.gz: a77168a66e7fadf40bcfa171dbd6c537
```

## Current Performance

From last benchmark (200kb test, 352k reads):
- Python-only: 68.4s
- Rust-accelerated: 47.4s
- **Speedup: 1.45x**

Key bottleneck: mRNA counting (80% of runtime, not in Rust)

## Files Created This Session

- `test_data/` - Self-contained test dataset (23MB)
- `setup_cloud_dev.sh` - One-command setup
- `CLOUD_DEV_SETUP.md` - Full dev environment guide
- `ci_validation.py` - MD5 checksum validation
- `validate_chr19_test.py` - Data integrity tests
- `create_single_chrom_release.sh` - Package builder
- `expected_checksums.json` - Deterministic validation
