# Building Sheriff-rs Python Bindings

This guide explains how to build and use the Python bindings for the Sheriff-rs Rust library.

## Prerequisites

1. **Rust toolchain** (install from https://rustup.rs/):
   ```bash
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   source $HOME/.cargo/env
   ```

2. **Python 3.7+** with pip

3. **Maturin** (Rust/Python build tool):
   ```bash
   pip install maturin
   ```

## Building the Python Module

### Option 1: Development Build (Recommended for Testing)

Build and install in development mode with editable install:

```bash
cd sheriff-rs
maturin develop --release --features python
```

This will:
- Compile the Rust code with optimizations (`--release`)
- Enable Python bindings (`--features python`)
- Install the module in your current Python environment
- Allow you to immediately import and use `sheriff_rs`

### Option 2: Production Wheel

Build a redistributable wheel file:

```bash
cd sheriff-rs
maturin build --release --features python
```

The wheel will be created in `target/wheels/` and can be installed with:

```bash
pip install target/wheels/sheriff_rs-*.whl
```

### Option 3: Direct Installation

Install directly without keeping the wheel:

```bash
cd sheriff-rs
maturin build --release --features python --interpreter python3
pip install target/wheels/sheriff_rs-*.whl
```

## Verifying the Installation

After building, verify the installation:

```bash
python3 -c "import sheriff_rs; print(sheriff_rs.__version__)"
```

You should see the version number printed (e.g., `0.1.0`).

## Troubleshooting

### Issue: "feature `python` is not enabled"

**Solution:** Make sure to include `--features python` in your build command:
```bash
maturin develop --release --features python
```

### Issue: "maturin: command not found"

**Solution:** Install maturin:
```bash
pip install maturin
```

### Issue: Rust compiler not found

**Solution:** Install Rust:
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
```

### Issue: Build fails with linking errors

**Solution:** Make sure you have development tools installed:

**Ubuntu/Debian:**
```bash
sudo apt-get install build-essential python3-dev
```

**MacOS:**
```bash
xcode-select --install
```

### Issue: Module imports but functions are missing

**Solution:** This usually means the bindings weren't compiled. Rebuild with:
```bash
maturin develop --release --features python --force
```

## Quick Start Example

After installation, try this quick test:

```python
import sheriff_rs

# Test k-mer hashing
hash_val = sheriff_rs.kmer_to_num("ACGT")
print(f"Hash of 'ACGT': {hash_val}")  # Should print: 27

# Test k-mer matching
whitelist = [hash_val]
matches = sheriff_rs.match_kmer("ACGTACGT", 4, whitelist, output_hash=True)
print(f"Matches found: {matches}")  # Should print: [27, 27]

# Test UMI deduplication
umis = ["ATCGATCG", "ATCGATCC", "GCGCGCGC"]
unique_count = sheriff_rs.deduplicate_umis(umis, threshold=1)
print(f"Unique UMI groups: {unique_count}")  # Should print: 2
```

## Performance Notes

The Rust implementations provide significant speedups:
- **K-mer operations:** 4-14x faster than pure Python
- **UMI deduplication:** 3-6x faster than pure Python

For best performance:
1. Always build with `--release` flag
2. Reuse `KmerCounter` objects across multiple sequences
3. Use `deduplicate_umis()` instead of `deduplicate_umis_detailed()` when you only need counts

## Development Workflow

For iterative development:

1. Make changes to `src/python.rs` or other Rust files
2. Rebuild with `maturin develop --release --features python`
3. Test your changes in Python without reinstalling

## API Documentation

For detailed API documentation, see:
- Module docstring: `help(sheriff_rs)`
- Function docstrings: `help(sheriff_rs.kmer_to_num)`
- Class docstrings: `help(sheriff_rs.KmerCounter)`

Or run the example script:
```bash
python3 examples/python_demo.py
```

## Integration with Existing Python Code

To integrate into your existing Sheriff Python codebase:

1. Import the module:
   ```python
   import sheriff_rs
   ```

2. Replace bottleneck functions:
   ```python
   # Old Python code:
   # hash_val = self.kmer_to_num(kmer)

   # New Rust-accelerated code:
   hash_val = sheriff_rs.kmer_to_num(kmer)
   ```

3. Use `KmerCounter` for repeated counting:
   ```python
   counter = sheriff_rs.KmerCounter(k=6)
   for sequence in sequences:
       freqs = counter.count_kmers(sequence)
       # Process frequencies...
   ```

## Building for Different Python Versions

Maturin automatically detects your Python version. To build for multiple versions:

```bash
# Build for all Python versions found on your system
maturin build --release --features python --interpreter python3.8 python3.9 python3.10 python3.11
```

## Cross-Platform Notes

### Linux
- Most distributions work out of the box
- Make sure `python3-dev` is installed

### macOS
- Requires Xcode command line tools
- May need to set `MACOSX_DEPLOYMENT_TARGET`:
  ```bash
  export MACOSX_DEPLOYMENT_TARGET=10.9
  maturin build --release --features python
  ```

### Windows
- Requires Visual Studio C++ Build Tools
- Use PowerShell or cmd.exe (not WSL) for building
- May need to specify Python path explicitly

## Next Steps

- See `examples/python_demo.py` for comprehensive usage examples
- Run benchmarks to verify performance improvements
- Integrate into your existing Sheriff Python codebase

## Support

For issues or questions:
1. Check the troubleshooting section above
2. Review PyO3 documentation: https://pyo3.rs
3. Review Maturin documentation: https://maturin.rs
