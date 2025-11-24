# Sheriff Development Guide

This guide helps you set up a development environment for contributing to Sheriff.

## Table of Contents

- [Setting Up Development Environment](#setting-up-development-environment)
- [Building Rust Modules](#building-rust-modules)
- [Running Tests](#running-tests)
- [Code Style Guidelines](#code-style-guidelines)
- [Submitting Pull Requests](#submitting-pull-requests)
- [Release Process](#release-process)

---

## Setting Up Development Environment

### Prerequisites

- Python 3.10 or later
- Rust toolchain (for Rust acceleration)
- Git
- conda or mamba (recommended)

### Initial Setup

```bash
# Clone the repository
git clone https://github.com/BradBalderson/Sheriff.git
cd Sheriff

# Create conda environment
conda create -n sheriff-dev python=3.11
conda activate sheriff-dev

# Install system dependencies
conda install -c bioconda samtools pysam
conda install -c conda-forge scipy numpy gtfparse faiss-cpu numba biopython

# Install Sheriff in editable mode
pip install -e .[dev]

# Install Rust toolchain
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env

# Install maturin for Rust-Python bindings
pip install maturin

# Build Rust module
cd sheriff-rs
maturin develop --release
cd ..
```

### Verify Installation

```bash
# Test Sheriff CLI
sheriff --help

# Test Python import
python -c "import sheriff; from sheriff.count_t7 import HAS_RUST_KMER; print(f'Rust available: {HAS_RUST_KMER}')"

# Run test suite
pytest tests/ -v
```

---

## Building Rust Modules

### Development Build (Fast, Debug Symbols)

```bash
cd sheriff-rs
maturin develop
cd ..
```

### Release Build (Optimized, Slower Compile)

```bash
cd sheriff-rs
maturin develop --release
cd ..
```

### Run Rust Tests

```bash
cd sheriff-rs
cargo test --verbose
cargo clippy -- -D warnings
cargo fmt -- --check
cd ..
```

### Build Wheels for Distribution

```bash
cd sheriff-rs
maturin build --release --out ../dist-rust/
cd ..
```

---

## Running Tests

### Quick Test Suite (Python Only)

```bash
pytest tests/ -v -m "not rust"
```

### Full Test Suite (Python + Rust)

```bash
# Ensure Rust module is built
cd sheriff-rs && maturin develop --release && cd ..

# Run all tests
pytest tests/ -v
```

### With Coverage

```bash
pytest tests/ -v --cov=sheriff --cov-report=html --cov-report=term
```

View coverage report:
```bash
open htmlcov/index.html  # macOS
xdg-open htmlcov/index.html  # Linux
```

### Run Specific Test File

```bash
pytest tests/test_rust_integration.py -v
```

### Run Specific Test Function

```bash
pytest tests/test_rust_integration.py::test_kmer_matching -v
```

### Run Benchmarks

```bash
# K-mer benchmark
python benchmarks/benchmark_kmer.py

# BAM filtering benchmark
python benchmarks/compare_all_modes.py

# Full benchmark suite
bash scripts/run_benchmarks.sh
```

---

## Code Style Guidelines

Sheriff follows modern Python best practices and Rust conventions.

### Python Style

**Formatting**: Black with 120 character line length
```bash
black sheriff/ tests/ benchmarks/ --line-length 120
```

**Linting**: Ruff
```bash
ruff check sheriff/ tests/ benchmarks/ --fix
```

**Type Checking**: MyPy
```bash
mypy sheriff/ --ignore-missing-imports
```

**Docstrings**: NumPy style
```python
def my_function(param1: str, param2: int) -> bool:
    """
    Short description of function.

    Longer description explaining what the function does,
    its purpose, and any important details.

    Parameters
    ----------
    param1 : str
        Description of param1
    param2 : int
        Description of param2

    Returns
    -------
    bool
        Description of return value

    Examples
    --------
    >>> my_function("hello", 42)
    True
    """
    ...
```

### Rust Style

**Formatting**: rustfmt
```bash
cd sheriff-rs
cargo fmt
```

**Linting**: clippy
```bash
cd sheriff-rs
cargo clippy -- -D warnings
```

**Documentation**: Standard Rust doc comments
```rust
/// Match k-mers in DNA sequence against whitelist
///
/// # Arguments
///
/// * `sequence` - DNA sequence string
/// * `k` - K-mer length
///
/// # Returns
///
/// List of matching k-mer hashes
///
/// # Examples
///
/// ```
/// let matches = match_kmer("AAACGTTT", 4, None, true);
/// ```
pub fn match_kmer(...) -> Vec<KmerMatch> {
    ...
}
```

### Pre-commit Hooks

Install pre-commit hooks to automatically check code style:

```bash
pip install pre-commit
pre-commit install

# Run manually on all files
pre-commit run --all-files
```

---

## Submitting Pull Requests

### Before Submitting

1. **Create a feature branch**:
   ```bash
   git checkout -b feature/my-new-feature
   ```

2. **Make your changes**:
   - Write code
   - Add tests
   - Update documentation

3. **Run the test suite**:
   ```bash
   pytest tests/ -v
   ```

4. **Run code quality checks**:
   ```bash
   black sheriff/ tests/
   ruff check sheriff/ tests/ --fix
   pre-commit run --all-files
   ```

5. **Update documentation**:
   - Update docstrings
   - Update README if needed
   - Add entry to CHANGELOG.md

6. **Commit your changes**:
   ```bash
   git add .
   git commit -m "Add feature: description"
   ```

### Pull Request Checklist

Use the [Pull Request Template](.github/PULL_REQUEST_TEMPLATE.md) and ensure:

- [ ] Code follows Sheriff's style guidelines
- [ ] All tests pass
- [ ] New tests added for new functionality
- [ ] Documentation updated
- [ ] CHANGELOG.md updated
- [ ] No new warnings introduced
- [ ] Performance impact documented (if applicable)

### Submitting the PR

```bash
git push origin feature/my-new-feature
```

Then create a pull request on GitHub.

---

## Release Process

### Version Numbering

Sheriff follows [Semantic Versioning](https://semver.org/):

- **Major** (1.0.0): Breaking changes
- **Minor** (1.1.0): New features, backward compatible
- **Patch** (1.1.1): Bug fixes

### Creating a Release

1. **Update version in `pyproject.toml`**:
   ```toml
   version = "1.2.1"
   ```

2. **Update CHANGELOG.md**:
   - Add release date
   - Summarize changes
   - Link to issues/PRs

3. **Commit version bump**:
   ```bash
   git add pyproject.toml CHANGELOG.md
   git commit -m "Bump version to 1.2.1"
   git push
   ```

4. **Create and push tag**:
   ```bash
   git tag v1.2.1
   git push origin v1.2.1
   ```

5. **GitHub Actions automatically**:
   - Runs full test suite
   - Builds distribution packages
   - Publishes to PyPI
   - Creates GitHub release

6. **Verify release**:
   ```bash
   pip install sheriff==1.2.1
   sheriff --version
   ```

---

## Project Structure

```
Sheriff/
├── .github/
│   ├── workflows/          # GitHub Actions CI/CD
│   └── ISSUE_TEMPLATE/     # Issue templates
├── sheriff/
│   ├── __init__.py
│   ├── __main__.py         # CLI entry point
│   ├── count_t7.py         # Main pipeline
│   ├── helpers.py          # Helper functions
│   └── bam_utils.py        # BAM utilities
├── sheriff-rs/
│   ├── src/
│   │   ├── lib.rs          # Rust library root
│   │   ├── bam_filter.rs   # BAM filtering
│   │   ├── kmer.rs         # K-mer matching
│   │   └── python.rs       # Python bindings
│   ├── Cargo.toml
│   └── tests/              # Rust tests
├── tests/
│   ├── test_rust_integration.py
│   ├── test_production_workflow.py
│   └── conftest.py         # Pytest fixtures
├── benchmarks/
│   ├── benchmark_kmer.py
│   └── compare_all_modes.py
├── docs/
│   └── source/
│       ├── api/            # API documentation
│       └── ...
├── example_data/           # Example datasets
├── pyproject.toml          # Package configuration
├── CHANGELOG.md            # Release history
├── DEVELOPMENT.md          # This file
└── README.md               # Project README
```

---

## Common Tasks

### Add a New Python Function

1. Add function to `sheriff/count_t7.py` or `sheriff/helpers.py`
2. Add docstring (NumPy style)
3. Add unit test in `tests/`
4. Update API documentation in `docs/source/api/`

### Add a New Rust Function

1. Add function to appropriate module in `sheriff-rs/src/`
2. Add Rust unit tests in the same file
3. Add Python binding in `sheriff-rs/src/python.rs`
4. Register in `sheriff_rs` module
5. Rebuild: `cd sheriff-rs && maturin develop --release`
6. Add Python integration test in `tests/`
7. Update Rust API docs in `docs/source/api/rust.rst`

### Update Documentation

```bash
# Edit .rst files in docs/source/
cd docs
make clean
make html

# View locally
open build/html/index.html  # macOS
xdg-open build/html/index.html  # Linux
```

### Profile Performance

```python
import cProfile
import pstats

cProfile.run('sheriff_main()', 'profile_stats')
p = pstats.Stats('profile_stats')
p.sort_stats('cumulative').print_stats(20)
```

---

## Troubleshooting

### Rust Build Fails

```bash
# Clean Rust build
cd sheriff-rs
cargo clean
cargo build --release

# Reinstall maturin
pip uninstall maturin
pip install maturin

# Rebuild
maturin develop --release
```

### Tests Fail

```bash
# Ensure Rust module is built
cd sheriff-rs && maturin develop --release && cd ..

# Run tests with verbose output
pytest tests/ -vv

# Run single test for debugging
pytest tests/test_rust_integration.py::test_specific_function -vv -s
```

### Pre-commit Hooks Fail

```bash
# Run individual hooks
pre-commit run black --all-files
pre-commit run ruff --all-files

# Skip hooks temporarily (not recommended)
git commit --no-verify
```

---

## Getting Help

- **Documentation**: https://sheriff.readthedocs.io/
- **Issues**: https://github.com/BradBalderson/Sheriff/issues
- **Discussions**: https://github.com/BradBalderson/Sheriff/discussions
- **Email**: bbalderson@salk.edu

---

## License

Sheriff is licensed under the BSD 3-Clause License. See [LICENSE.md](LICENSE.md) for details.

## Citation

If you use Sheriff in your research, please cite:

> Lorenzini et al. (2025). Joint single-cell profiling of CRISPR-Cas9 edits
> and transcriptomes reveals widespread off-target events and their effects on
> gene expression. bioRxiv 2025.02.07.636966.
