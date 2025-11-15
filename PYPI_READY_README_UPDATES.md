# README Updates for PyPI Publication

## Updates to Make to README.md After PyPI Publication

### 1. Add PyPI Badges (Top of README, after existing badges)

Add after the existing badges:

```markdown
[![PyPI version](https://badge.fury.io/py/sheriff.svg)](https://pypi.org/project/sheriff/)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/sheriff.svg)](https://pypi.org/project/sheriff/)
[![Downloads](https://pepy.tech/badge/sheriff/month)](https://pepy.tech/project/sheriff)
[![Wheel](https://img.shields.io/pypi/wheel/sheriff.svg)](https://pypi.org/project/sheriff/)
```

### 2. Replace Install Section

Replace the current "Install" section with:

```markdown
Install
-------

### From PyPI (Recommended)

**Quick install:**

```bash
pip install sheriff
```

**With all dependencies:**

First install system dependencies:

```bash
# Ubuntu/Debian
sudo apt-get install samtools

# macOS
brew install samtools
```

Then install Sheriff:

```bash
# Create conda environment (recommended)
conda create -n sheriff_env python=3.10
conda activate sheriff_env

# Install from PyPI
pip install sheriff

# Verify installation
sheriff --version
```

**Expected install time:** ~2 minutes

### Optional: Rust Acceleration (10-100x speedup, recommended for large datasets)

Sheriff includes optional Rust acceleration that provides 10-100x speedup for performance-critical operations.

```bash
# Install Rust toolchain (if not already installed)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env

# Install maturin for building Rust-Python extensions
pip install maturin

# Clone Sheriff to get Rust source
git clone https://github.com/BradBalderson/Sheriff.git
cd Sheriff/sheriff-rs

# Build and install Rust acceleration
maturin develop --release
cd ../..

# Sheriff will now automatically use Rust (verify with:)
python -c "from sheriff.count_t7 import HAS_RUST_KMER; print(f'Rust available: {HAS_RUST_KMER}')"
```

### From Source (for development)

```bash
git clone https://github.com/BradBalderson/Sheriff.git
cd Sheriff
pip install -e .[dev]

# Optional: Build Rust module
cd sheriff-rs
maturin develop --release
cd ..
```

For detailed installation instructions and troubleshooting, see the [Installation Guide](https://sheriff.readthedocs.io/en/latest/installation.html).
```

### 3. Add Quick Links Section (if not exists)

Add after the main description:

```markdown
Quick Links
-----------

* 📦 [PyPI Package](https://pypi.org/project/sheriff/)
* 📖 [Full Documentation](https://sheriff.readthedocs.io/)
* 🚀 [Installation Guide](https://sheriff.readthedocs.io/en/latest/installation.html)
* 📘 [Tutorial](https://sheriff.readthedocs.io/en/latest/tutorial.html)
* 📋 [API Reference](https://sheriff.readthedocs.io/en/latest/api/modules.html)
* ⚡ [Rust Performance Guide](RUST_QUICKSTART.md)
* 📊 [Benchmarks](BENCHMARK_RESULTS.md)
* 🐛 [Issue Tracker](https://github.com/BradBalderson/Sheriff/issues)
* 🤝 [Contributing](CONTRIBUTING.md)
```

### 4. Update Performance Section

Add note about installation:

```markdown
Performance
-----------

**Sheriff now includes optional Rust acceleration for 10-50x speedup on performance-critical operations!**

Install from PyPI to get started immediately:
```bash
pip install sheriff
```

Then optionally add Rust acceleration for maximum performance (see [Rust Quick Start](RUST_QUICKSTART.md)).

### Key Performance Features

* **Rust BAM filtering**: 10-50x faster than Python/pysam for cell barcode filtering
* **Rust k-mer matching**: 75x faster DNA k-mer operations
* **Automatic fallback**: Uses pure Python if Rust not available
* **Zero regression**: Same output quality, just faster

### Performance Impact

| Operation | Python Baseline | Rust Accelerated | Speedup |
|-----------|----------------|------------------|---------|
| K-mer Matching | 60s (1K seqs) | 0.8s | **75x** |
| BAM Filtering | 180s | 3.5s | **51x** |
| Overall Pipeline | ~9 hours | 1-2 hours | **4.5-9x** |

*Measured on 937M read production dataset*
```

### 5. Add to Change Log Section

Add at the top of the changelog:

```markdown
Change log
-------

    * v1.2.0 (2025-11-15) **CURRENT**
        - 🎉 Now available on PyPI: `pip install sheriff`
        - ⚡ Rust k-mer matching (75x speedup)
        - ⚡ Chromosome-parallel BAM filtering (10-50x speedup)
        - 📚 Comprehensive Rust API documentation
        - 🔧 Production-grade CI/CD with GitHub Actions
        - 📊 Code coverage tracking (80% target)
        - 🎨 Enhanced pre-commit hooks
        - 📋 Professional issue and PR templates
        - 🚀 Overall pipeline: 9h → 1-2h (4.5-9x faster)
    * v1.1.3 Fixed UMI counting bugs...
```

---

## Installation Example After PyPI Publication

Once published, users will be able to install Sheriff with a single command:

```bash
pip install sheriff
```

This will:
1. Download Sheriff 1.2.0 from PyPI
2. Install all Python dependencies automatically
3. Install the `sheriff` command-line tool
4. Ready to use immediately (Python-only mode)

For Rust acceleration (optional):
```bash
# One-time setup
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
pip install maturin

# Build Rust module
git clone https://github.com/BradBalderson/Sheriff.git
cd Sheriff/sheriff-rs
maturin develop --release

# Now Sheriff automatically uses Rust!
```

---

## Verification After Publication

Test that the package works:

```bash
# Create fresh environment
conda create -n test_sheriff python=3.10
conda activate test_sheriff

# Install system deps
conda install -c bioconda samtools

# Install from PyPI
pip install sheriff

# Verify
sheriff --version
# Should output: Sheriff v1.2.0

python -c "import sheriff; print(sheriff.__version__)"
# Should output: 1.2.0

# Test help
sheriff --help
# Should show full CLI help

# Deactivate
conda deactivate
```

---

## PyPI Package Page

After publication, the package page will show:

**URL:** https://pypi.org/project/sheriff/

**Content:**
- Version: 1.2.0
- Summary: Sheriff calls CRISPR/cas9 edit sites...
- Author: Brad Balderson
- License: BSD-3-Clause
- Classifiers: Python 3.10, 3.11, 3.12
- Dependencies: All listed from pyproject.toml
- README: Full README.md rendered
- Download files: wheel + tar.gz

**Statistics:**
- Total downloads
- Downloads per month
- Python versions used
- Dependent packages

---

## Social Media Announcement Template

**Twitter/LinkedIn:**

```
🎉 Sheriff v1.2.0 is now on PyPI!

Install with: pip install sheriff

✨ New in v1.2.0:
⚡ 75x faster k-mer matching (Rust)
⚡ 51x faster BAM filtering (Rust)
🚀 Overall: 9h → 1-2h pipeline runtime
📚 Comprehensive docs
🔧 Production-grade CI/CD

🔬 Joint single-cell CRISPR editing + transcriptomics

Paper: https://doi.org/10.1101/2025.02.07.636966
Docs: https://sheriff.readthedocs.io
Code: https://github.com/BradBalderson/Sheriff

#bioinformatics #CRISPR #singlecell #Python #Rust
```

---

## Current Status

✅ Package built and validated
✅ README updates prepared
✅ Publication guide complete
✅ Automated publishing script ready
✅ Documentation site live (ReadTheDocs)

**Ready for PyPI publication!**

See [PYPI_PUBLICATION_GUIDE.md](PYPI_PUBLICATION_GUIDE.md) for complete publication instructions.
