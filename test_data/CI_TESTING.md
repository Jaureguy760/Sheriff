# Sheriff CI/CD Testing with Real Genomic Data

This directory contains **real experimental data** from a CRISPR editing experiment, used for continuous integration and validation.

---

## Test Data Files

### Core Test Files (Self-Contained)

| File | Size | Description |
|------|------|-------------|
| `test_200kb.bam` | 30 MB | 352,535 reads from edited cells |
| `test_200kb.bam.bai` | 97 KB | BAM index |
| `barcodes.txt` | 6.5 KB | 583 valid cell barcodes |
| `edit_sites.txt` | 1.4 KB | 43 known CRISPR edit sites |
| `blacklist.bed` | 50 KB | 1,362 noise regions to exclude |
| `blacklist_seqs.txt` | 170 B | 6 TSO sequences to filter |

**Total:** ~31 MB (fast CI testing!)

### What This Data Represents

**Biological Context:**
- Superb-seq experiment: CRISPR editing + RNA-seq
- T7 donor template (GGGAGAGTAT) marks edit sites
- Split-pipe barcoded: Cell (CB tag) + UMI (UB tag)
- Aligned to GRCh38/hg38 chromosome 19 region

**Data Stats:**
- 352,535 reads processed
- 583 valid cells (from 500-cell library)
- 304,213 reads match whitelisted barcodes (86%)
- 43 ground-truth edit sites for validation

---

## CI Tests Run

### 1. `ci_validation.py` - Rust Function Validation

**What it tests:**
```bash
python test_data/ci_validation.py
```

✅ **Input file checksums** - Data integrity
✅ **Rust module import** - sheriff_rs loads
✅ **K-mer counting** - 4-mer hashing with MD5 checksum
✅ **UMI deduplication** - Hamming distance clustering
✅ **Edit clustering** - Longest edit selection algorithm
✅ **Cell UMI counting** - Per-cell molecule counting

**Purpose:** Validates core Rust functions with checksummed outputs

---

### 2. `test_integration_smoke.py` - API Integration

**What it tests:**
```bash
python test_integration_smoke.py
```

✅ **Function signatures** - New parameters exist
✅ **Object creation** - PipelineProgress, CheckpointManager, PipelineResults
✅ **Backward compatibility** - Old code still works

**Purpose:** Ensures Priority 3 features integrate correctly

---

### 3. `test_priority3_integration.py` - Feature Validation

**What it tests:**
```bash
python test_priority3_integration.py
```

✅ **Checkpoint save/load** - State persistence
✅ **Progress tracking** - Progress bars render
✅ **Results display** - Rich table formatting
✅ **CLI flags** - --resume, --enable-checkpoints, --checkpoint

**Purpose:** Validates user-facing Priority 3 features

---

## CI/CD Pipeline (.github/workflows/test.yml)

### On Every Push/PR:

1. **Install Dependencies**
   - Python 3.10 & 3.11 matrix
   - samtools, bcftools, bedtools
   - Rust toolchain (stable)

2. **Build Rust Module**
   - `maturin build --release --features python`
   - Install sheriff_rs wheel

3. **Validate Test Data**
   - Check files exist (352k reads)
   - Verify checksums match

4. **Run All Tests**
   - CI validation (Rust functions)
   - Integration smoke tests
   - Full Priority 3 validation
   - Rust acceleration smoke test

5. **Report Results**
   - ✅ All tests pass → PR can merge
   - ❌ Any test fails → PR blocked

---

## Pre-commit Hooks (.pre-commit-config.yaml)

**Runs locally before commit:**

✅ **Code formatting** - Black, isort
✅ **Rust formatting** - rustfmt
✅ **File hygiene** - trailing whitespace, YAML syntax
✅ **Test data validation** - Checksums on push

**Install:**
```bash
pip install pre-commit
pre-commit install
```

Now every commit auto-checks code quality!

---

## Why Use Real Data for CI?

### Traditional Approach (Bad)
```python
def test_edit_clustering():
    # Mock data
    edits = [("chr1", 100, "A", "T", True, [])]
    result = cluster(edits)
    assert len(result) == 1  # Passes but meaningless!
```

### Sheriff Approach (Good) ✅
```python
def test_edit_clustering():
    # REAL 352k reads from actual experiment
    # REAL cell barcodes from real cells
    # REAL edit sites we know exist
    result = run_sheriff_on_real_data()
    assert MD5(result) == "known_good_checksum"
```

**Benefits:**
- Catches real-world bugs
- Validates against known biological truth
- Tests actual I/O (BAM files, not mocks)
- Checksums ensure bit-exact reproducibility

---

## Adding New Tests

### For New Features:

1. **Add to ci_validation.py**
   ```python
   print("=== Testing New Feature ===")
   result = sheriff_rs.new_function(test_data)
   expected_checksum = hashlib.md5(str(result).encode()).hexdigest()
   print(f"  Checksum: {expected_checksum}")
   assert result == expected_output
   ```

2. **Run locally first**
   ```bash
   python test_data/ci_validation.py
   ```

3. **Commit changes**
   - GitHub Actions will run on push
   - PR blocked if tests fail

### For Optimization Validation:

**ALWAYS compare outputs before claiming speedup!**

```python
# Get baseline output
baseline = run_baseline_version(test_data)

# Get optimized output
optimized = run_optimized_version(test_data)

# MUST BE IDENTICAL
assert baseline == optimized  # ← THIS IS CRITICAL!

# Only then measure speed
time_baseline = benchmark(run_baseline_version)
time_optimized = benchmark(run_optimized_version)
speedup = time_baseline / time_optimized
```

**Lesson learned:** Speed without correctness is worthless!

---

## Test Data Provenance

**Source:** BradBalderson/Sheriff example_data
**Experiment:** Superb-seq CRISPR editing with T7 donor
**Assembly:** GRCh38/hg38
**Region:** Chromosome 19 (±200kb of edit sites)
**Citation:** [Sheriff bioRxiv preprint]

**Checksums (MD5):**
```
test_200kb.bam:     d6f65018...
barcodes.txt:       478d09b5...
edit_sites.txt:     080071d6...
blacklist.bed:      6a923066...
blacklist_seqs.txt: ff3dffed...
```

Run `python ci_validation.py` to verify integrity!

---

## Summary

✅ **352k reads** tested on every commit
✅ **Rust functions** validated with checksums
✅ **Real genomic data** catches real bugs
✅ **Auto-blocks** broken PRs
✅ **Fast CI** (~2-3 min per run)

**This is production-grade testing for bioinformatics!** 🧬
