# Contributing to Sheriff

Thank you for your interest in contributing to Sheriff! This document provides guidelines and instructions for contributing.

## Getting Started

### Prerequisites

- Python 3.10 or higher
- Conda or Mamba (recommended)
- Git

### Development Setup

1. **Fork and clone the repository**
   ```bash
   git clone https://github.com/YOUR_USERNAME/Sheriff.git
   cd Sheriff
   ```

2. **Create a development environment**
   ```bash
   mamba create -n sheriff_dev python=3.10
   mamba activate sheriff_dev
   ```

3. **Install dependencies**
   ```bash
   # Install main dependencies
   pip install -r requirements.txt

   # Install development dependencies
   pip install pytest pytest-cov black ruff mypy

   # Install Sheriff in editable mode
   pip install -e .
   ```

4. **Verify installation**
   ```bash
   sheriff --help
   ```

## Development Workflow

### Code Style

We follow PEP 8 style guidelines with some modifications:

- **Line length**: Maximum 120 characters
- **Formatter**: Black (configured in `pyproject.toml`)
- **Linter**: Ruff and Flake8
- **Type checking**: MyPy (gradual adoption)

#### Format your code

```bash
# Format with Black
black sheriff/

# Lint with Ruff
ruff check sheriff/

# Type check with MyPy
mypy sheriff/
```

### Making Changes

1. **Create a new branch**
   ```bash
   git checkout -b feature/your-feature-name
   # or
   git checkout -b fix/your-bugfix-name
   ```

2. **Make your changes**
   - Write clear, concise code
   - Add comments for complex logic
   - Update documentation if needed

3. **Test your changes**
   ```bash
   # Run tests (when available)
   pytest tests/

   # Test on example data
   ./run_example.sh  # If you create an example script
   ```

4. **Commit your changes**
   ```bash
   git add .
   git commit -m "Clear description of your changes"
   ```

   **Commit message guidelines:**
   - Use present tense ("Add feature" not "Added feature")
   - Be descriptive but concise
   - Reference issues if applicable (e.g., "Fix #123: Description")

### Pull Request Process

1. **Update your branch**
   ```bash
   git fetch upstream
   git rebase upstream/main
   ```

2. **Push to your fork**
   ```bash
   git push origin feature/your-feature-name
   ```

3. **Create a Pull Request**
   - Go to the [Sheriff repository](https://github.com/BradBalderson/Sheriff)
   - Click "New Pull Request"
   - Select your fork and branch
   - Fill in the PR template with:
     - Description of changes
     - Related issues
     - Testing performed
     - Any breaking changes

4. **Address review feedback**
   - Make requested changes
   - Push updates to your branch
   - Respond to reviewer comments

## Contribution Guidelines

### Code Quality

- **No bare exceptions**: Use specific exception types
  ```python
  # Bad
  try:
      risky_operation()
  except:
      pass

  # Good
  try:
      risky_operation()
  except (ValueError, KeyError) as e:
      logger.warning(f"Operation failed: {e}")
  ```

- **Use logging instead of print**: For library code, use Python's logging module
  ```python
  import logging
  logger = logging.getLogger(__name__)
  logger.info("Processing started")
  ```

- **Type hints**: Add type hints for new functions (gradual adoption)
  ```python
  def process_reads(bam_file: str, verbosity: int = 1) -> pd.DataFrame:
      ...
  ```

### Documentation

- **Docstrings**: Add docstrings for all public functions and classes
  ```python
  def match_kmer(bc_kmer_matcher, indel_seq, output_kmer_hash):
      """
      Match k-mers in a sequence against a reference set.

      Args:
          bc_kmer_matcher: KmerMatcher object with reference k-mers
          indel_seq: Sequence to search for k-mer matches
          output_kmer_hash: Whether to return hashes or k-mer strings

      Returns:
          Tuple of matched k-mers or None if no matches
      """
  ```

- **Update README**: If you add features, update the README.md
- **Changelog**: Add your changes to the changelog section in README.md

### Testing

While we're building our test suite, please:

- Test your changes with the example data
- Verify edge cases manually
- Document expected behavior
- Consider writing tests for new features

**Future test structure:**
```
tests/
├── test_kmer_matching.py
├── test_edit_calling.py
├── test_umi_counting.py
└── fixtures/
    └── test_data.bam
```

### Performance Considerations

Sheriff processes large genomic datasets, so performance matters:

- **Profile new code**: Use `cProfile` or `timeit` for bottlenecks
- **Consider Numba**: For tight loops, consider Numba JIT compilation
- **Memory efficiency**: Use appropriate data types (e.g., `np.uint32` for counts)
- **Parallel processing**: Leverage the existing multiprocessing framework

## Reporting Issues

### Bug Reports

Include:
- Sheriff version (`sheriff --version` or check `sheriff/__init__.py`)
- Python version and OS
- Steps to reproduce
- Expected vs. actual behavior
- Error messages (full traceback)
- Sample data (if shareable)

### Feature Requests

Include:
- Use case description
- Proposed solution
- Potential impact
- Willingness to contribute implementation

## Questions?

- **Documentation**: Check the [README.md](README.md)
- **Issues**: Search [existing issues](https://github.com/BradBalderson/Sheriff/issues)
- **Contact**: Email bbalderson@salk.edu

## Citation

If you contribute to Sheriff, you may be acknowledged in future publications. Please include citation information in your contributions:

```
Joint single-cell profiling of CRISPR-Cas9 edits and transcriptomes reveals widespread off-target events and their effects on gene expression
Michael H. Lorenzini, Brad Balderson, Karthyayani Sajeev, Aaron J. Ho, Graham McVicker
bioRxiv 2025.02.07.636966; doi: https://doi.org/10.1101/2025.02.07.636966
```

## License

By contributing, you agree that your contributions will be licensed under the BSD License used by this project.

---

Thank you for contributing to Sheriff! 🎯
