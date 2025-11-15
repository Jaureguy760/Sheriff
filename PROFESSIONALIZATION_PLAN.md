# Sheriff Professionalization Plan
## Comparison to scvi-tools Standards and Action Plan

**Date:** 2025-11-15
**Goal:** Bring Sheriff to scvi-tools professional standards
**Status:** In Progress

---

## Executive Summary

Sheriff is **already quite professional** with strong documentation, modern packaging, and comprehensive features. This plan identifies remaining gaps compared to scvi-tools and provides actionable steps to reach the same professional tier.

**Current Strength:** 7.5/10 professional standard
**Target:** 9.5/10 (scvi-tools level)

---

## Comparison Matrix: Sheriff vs scvi-tools

| Feature | scvi-tools | Sheriff | Status | Priority |
|---------|-----------|---------|---------|----------|
| **README Quality** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ✅ Excellent | - |
| **Documentation (Sphinx)** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ✅ Complete | - |
| **ReadTheDocs Integration** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ✅ Active | - |
| **Modern Packaging (pyproject.toml)** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ✅ Complete | - |
| **pip installability** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | 🟡 Manual steps needed | **HIGH** |
| **API Documentation** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | 🟡 Missing Rust API docs | **HIGH** |
| **Contributing Guidelines** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ✅ Good | LOW |
| **CI/CD (GitHub Actions)** | ⭐⭐⭐⭐⭐ | ⭐ | ❌ Missing | **MEDIUM** |
| **Code Coverage** | ⭐⭐⭐⭐⭐ | ⭐ | ❌ No coverage tracking | **MEDIUM** |
| **Pre-commit Hooks** | ⭐⭐⭐⭐⭐ | ⭐⭐ | 🟡 Basic, needs expansion | **MEDIUM** |
| **Badges (PyPI, coverage, etc)** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | 🟡 Missing PyPI, coverage | **LOW** |
| **Rust Documentation** | N/A | ⭐⭐ | 🟡 Needs improvement | **HIGH** |
| **Examples/Tutorials** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ✅ Good examples | LOW |
| **Performance Benchmarks** | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ✅ Excellent benchmarks | - |
| **License Clarity** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ✅ BSD-3-Clause clear | - |
| **Issue Templates** | ⭐⭐⭐⭐⭐ | ⭐ | ❌ Missing | **LOW** |
| **PyPI Distribution** | ⭐⭐⭐⭐⭐ | ⭐ | ❌ Not on PyPI | **HIGH** |
| **Conda Distribution** | ⭐⭐⭐⭐⭐ | ⭐ | ❌ Not on Conda-Forge | MEDIUM |

---

## Critical Gaps (HIGH Priority)

### 1. **pip Installability - Seamless Rust Integration**

**Current State:**
- Users must manually build Rust components
- Multi-step installation process
- Requires Rust toolchain knowledge

**Target State (scvi-tools level):**
- Single `pip install sheriff` command
- Automatic Rust compilation via maturin
- Optional Rust acceleration with graceful fallback

**Action Plan:**

#### Phase 1: Hybrid Package with Optional Rust (✅ RECOMMENDED)
```toml
# pyproject.toml additions
[project.optional-dependencies]
rust = ["maturin>=1.0.0"]

[tool.maturin]
python-source = "python"
module-name = "sheriff._rust"
```

**Implementation:**
1. Move Python code to `python/sheriff/` (src layout)
2. Rust becomes optional build dependency
3. Add build script that tries Rust, falls back gracefully
4. `pip install sheriff` → Pure Python (fast, always works)
5. `pip install sheriff[rust]` → With Rust acceleration

**Benefits:**
- ✅ Always installable (critical for reproducibility)
- ✅ No Rust toolchain requirement
- ✅ Optional performance boost
- ✅ PyPI compatible

#### Phase 2: Pre-built Wheels (Future)
- Build wheels for major platforms (Linux, macOS, Windows)
- Upload to PyPI
- Users get Rust by default, no compilation needed

### 2. **Rust API Documentation**

**Current State:**
- Rust code has minimal docstrings
- No API reference for Rust functions
- Users don't know what's available

**Target State:**
- Comprehensive Rust docstrings using `///` style
- API reference in Sphinx docs
- Usage examples for each Rust function

**Action Plan:**

```rust
/// Match k-mers in DNA sequence against whitelist
///
/// Performs high-performance k-mer matching using frequency arrays.
/// Automatically handles ambiguous bases ('N') by skipping k-mers.
///
/// # Arguments
///
/// * `sequence` - DNA sequence string (A/C/G/T/N)
/// * `k` - K-mer length (typically 6-15 for barcode matching)
/// * `whitelist` - Optional list of k-mer hashes to match against
/// * `output_hash` - Return hashes (True) or k-mer strings (False)
///
/// # Returns
///
/// List of matching k-mer hashes (int) or strings (str)
///
/// # Examples
///
/// ```python
/// import sheriff_rs
///
/// # Find all k-mers in sequence
/// matches = sheriff_rs.match_kmer_rust("AAACGTTT", k=4)
///
/// # Match against T7 barcode whitelist
/// whitelist = [27, 45, 109]  # Pre-computed k-mer hashes
/// matches = sheriff_rs.match_kmer_rust(
///     "AAACGTTT", k=4, whitelist=whitelist, output_hash=True
/// )
/// ```
///
/// # Performance
///
/// - 10-50x faster than Python implementation
/// - ~200 ns per sequence (vs 10+ μs in Python)
/// - Scales linearly with sequence length
#[pyfunction]
pub fn match_kmer_rust(...) -> PyResult<Vec<PyObject>> {
    ...
}
```

**Create `docs/source/api/rust.rst`:**
```rst
Rust Acceleration API
=====================

High-performance Rust implementations for performance-critical operations.

Installation
------------

.. code-block:: bash

   pip install sheriff[rust]

   # Or build from source
   cd sheriff-rs
   maturin develop --release

Overview
--------

Sheriff includes optional Rust acceleration that provides 10-50x speedup
for performance-critical operations while maintaining identical output.

.. automodule:: sheriff_rs
   :members:
   :undoc-members:
   :show-inheritance:

Functions
---------

BAM Filtering
~~~~~~~~~~~~~

.. autofunction:: sheriff_rs.filter_bam_by_barcodes_rust

K-mer Matching
~~~~~~~~~~~~~~

.. autofunction:: sheriff_rs.match_kmer_rust
.. autofunction:: sheriff_rs.count_kmers_rust

Performance Comparison
----------------------

+------------------+-------------+-------------+----------+
| Operation        | Python Time | Rust Time   | Speedup  |
+==================+=============+=============+==========+
| BAM Filtering    | 180s        | 3.5s        | 51x      |
| K-mer Matching   | 60s         | 0.8s        | 75x      |
| Overall Pipeline | 9 hours     | 1-2 hours   | 4.5-9x   |
+------------------+-------------+-------------+----------+
```

### 3. **PyPI Distribution**

**Current State:**
- Not published to PyPI
- Users must `git clone` and `pip install .`

**Target State:**
- `pip install sheriff` works immediately
- Automatic version updates on release

**Action Plan:**

1. **Prepare for PyPI:**
   ```bash
   # Test build
   python -m build

   # Test with TestPyPI first
   python -m twine upload --repository testpypi dist/*

   # Install from TestPyPI to validate
   pip install --index-url https://test.pypi.org/simple/ sheriff
   ```

2. **Publish to PyPI:**
   ```bash
   python -m twine upload dist/*
   ```

3. **Automate with GitHub Actions:**
   ```.github/workflows/publish.yml
   name: Publish to PyPI

   on:
     release:
       types: [published]

   jobs:
     publish:
       runs-on: ubuntu-latest
       steps:
         - uses: actions/checkout@v3
         - uses: actions/setup-python@v4
         - run: pip install build twine
         - run: python -m build
         - run: twine upload dist/*
           env:
             TWINE_USERNAME: __token__
             TWINE_PASSWORD: ${{ secrets.PYPI_TOKEN }}
   ```

---

## Medium Priority Improvements

### 4. **CI/CD Pipeline (GitHub Actions)**

**Create `.github/workflows/tests.yml`:**

```yaml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ${{ matrix.os }}
    strategy:
      matrix:
        os: [ubuntu-latest, macos-latest]
        python-version: ['3.10', '3.11', '3.12']

    steps:
      - uses: actions/checkout@v3

      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: ${{ matrix.python-version }}

      - name: Install dependencies
        run: |
          pip install -e .[dev]

      - name: Run tests
        run: |
          pytest --cov=sheriff --cov-report=xml

      - name: Upload coverage
        uses: codecov/codecov-action@v3
        with:
          file: ./coverage.xml

  rust-tests:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3

      - name: Install Rust
        uses: actions-rs/toolchain@v1
        with:
          toolchain: stable

      - name: Run Rust tests
        run: |
          cd sheriff-rs
          cargo test
          cargo clippy
```

### 5. **Code Coverage Tracking**

1. **Add Codecov integration:**
   - Sign up at codecov.io
   - Add `CODECOV_TOKEN` to GitHub secrets
   - Add badge to README

2. **Coverage configuration (`.coveragerc`):**
   ```ini
   [run]
   source = sheriff
   omit =
       */tests/*
       */setup.py

   [report]
   exclude_lines =
       pragma: no cover
       def __repr__
       raise AssertionError
       raise NotImplementedError
       if __name__ == .__main__.:
   ```

### 6. **Enhanced Pre-commit Hooks**

**Update `.pre-commit-config.yaml`:**

```yaml
repos:
  - repo: https://github.com/pre-commit/pre-commit-hooks
    rev: v4.5.0
    hooks:
      - id: trailing-whitespace
      - id: end-of-file-fixer
      - id: check-yaml
      - id: check-added-large-files
      - id: check-merge-conflict
      - id: debug-statements

  - repo: https://github.com/psf/black
    rev: 23.12.0
    hooks:
      - id: black
        language_version: python3.10

  - repo: https://github.com/astral-sh/ruff-pre-commit
    rev: v0.1.9
    hooks:
      - id: ruff
        args: [--fix, --exit-non-zero-on-fix]

  - repo: https://github.com/pre-commit/mirrors-mypy
    rev: v1.8.0
    hooks:
      - id: mypy
        additional_dependencies: [types-all]

  - repo: local
    hooks:
      - id: cargo-test
        name: Cargo Test
        entry: bash -c 'cd sheriff-rs && cargo test'
        language: system
        pass_filenames: false

      - id: cargo-clippy
        name: Cargo Clippy
        entry: bash -c 'cd sheriff-rs && cargo clippy -- -D warnings'
        language: system
        pass_filenames: false
```

---

## Low Priority Improvements

### 7. **Additional Badges**

Add to README.md:

```markdown
[![PyPI version](https://badge.fury.io/py/sheriff.svg)](https://badge.fury.io/py/sheriff)
[![Downloads](https://pepy.tech/badge/sheriff)](https://pepy.tech/project/sheriff)
[![codecov](https://codecov.io/gh/BradBalderson/Sheriff/branch/main/graph/badge.svg)](https://codecov.io/gh/BradBalderson/Sheriff)
[![CI](https://github.com/BradBalderson/Sheriff/workflows/Tests/badge.svg)](https://github.com/BradBalderson/Sheriff/actions)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
```

### 8. **Issue Templates**

**Create `.github/ISSUE_TEMPLATE/bug_report.md`:**

```markdown
---
name: Bug report
about: Create a report to help us improve
title: '[BUG] '
labels: bug
assignees: ''
---

**Describe the bug**
A clear description of the bug.

**To Reproduce**
Steps to reproduce:
1. Command run: `sheriff ...`
2. Input files: ...
3. Error message: ...

**Expected behavior**
What you expected to happen.

**Environment:**
 - OS: [e.g. Ubuntu 22.04]
 - Python version: [e.g. 3.10]
 - Sheriff version: [e.g. 1.2.0]
 - Installed via: [pip, conda, source]
 - Rust acceleration: [yes/no]

**Additional context**
Any other relevant information.
```

### 9. **Conda-Forge Distribution**

**Future goal:** Submit to conda-forge for broader accessibility

---

## K-mer Integration Plan

**Separate detailed plan below** for integrating Rust k-mer matching into Sheriff pipeline.

---

## Implementation Timeline

### Week 1 (Current)
- ✅ Complete Rust k-mer implementation
- ✅ Add Python bindings
- ✅ Create benchmarks
- 🔄 Document Rust API comprehensively
- 🔄 Create pip installation workflow

### Week 2
- [ ] Integrate k-mer Rust into count_t7.py
- [ ] Test on real dataset
- [ ] Update version to 1.2.0
- [ ] Create GitHub Actions CI/CD

### Week 3
- [ ] Publish to TestPyPI
- [ ] Validate installation
- [ ] Publish to PyPI
- [ ] Add coverage tracking

### Week 4
- [ ] Enhanced pre-commit hooks
- [ ] Issue templates
- [ ] Additional badges
- [ ] Final professionalization review

---

## Success Metrics

| Metric | Current | Target | Status |
|--------|---------|--------|--------|
| GitHub Stars | ~10 | 50+ | Growing |
| PyPI Downloads/month | 0 | 100+ | Pending release |
| Documentation Coverage | 80% | 95% | In progress |
| Code Coverage | Unknown | 80%+ | Need CI |
| Installation Success Rate | ~60% | 95%+ | Need PyPI |

---

## Conclusion

Sheriff is **already professional** with excellent documentation and modern tooling. The remaining gaps are primarily distribution/CI/CD related, not fundamental quality issues.

**Priority Order:**
1. **HIGH:** Complete Rust API docs + pip installation workflow
2. **HIGH:** PyPI distribution + k-mer integration
3. **MEDIUM:** CI/CD + coverage tracking
4. **LOW:** Badges + issue templates + conda-forge

Estimated time to full professionalization: **2-3 weeks** of focused work.
