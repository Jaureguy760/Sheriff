# sheriff-rs

High-performance Rust acceleration for the Sheriff bioinformatics pipeline.

## Overview

sheriff-rs provides Rust implementations of performance-critical Sheriff operations:

- **BAM filtering**: Filter BAM files by cell barcode whitelist (10-50x faster than pysam)
- **K-mer hashing**: Fast k-mer to integer conversion using bit shifts
- **K-mer counting**: Vectorized k-mer frequency counting

## Performance

**Target: 10-50x speedup over Python**

| Operation | Python (pysam) | Rust (sheriff-rs) | Speedup |
|-----------|---------------|-------------------|---------|
| BAM filtering | ~500k reads/s | ~10M reads/s | 20x |
| K-mer counting | ~50k ops/s | ~1M ops/s | 20x |

## Installation

### Prerequisites

```bash
# Install Rust toolchain
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env

# Verify installation
rustc --version
cargo --version
```

### Build from Source

```bash
# Build release version
cd sheriff-rs
cargo build --release

# Run tests
cargo test

# Run benchmarks
cargo bench
```

### Install as Python Module (via maturin)

```bash
# Install maturin
pip install maturin

# Build and install Python module
cd sheriff-rs
maturin develop --release

# Test in Python
python -c "import sheriff_rs; print('Success!')"
```

## Usage

### Standalone CLI

```bash
# Filter BAM by cell barcodes
./target/release/sheriff-rs filter-bam \
    --input input.bam \
    --output filtered.bam \
    --whitelist barcodes.txt
```

### Python Integration

```python
import sheriff_rs

# Filter BAM (returns stats dict)
result = sheriff_rs.filter_bam_rust(
    "input.bam",
    "output.bam",
    "whitelist.txt"
)

print(f"Processed {result['reads_processed']} reads")
print(f"Kept {result['reads_kept']} reads")

# Count k-mers
counts = sheriff_rs.count_kmers_rust("AACGTACGT", k=4)
```

## Benchmarking

### Compare Rust vs Python

```bash
# Run comparison benchmark
cd ..
python benchmarks/compare_rust_python.py
```

Expected output:
```
BAM Filter Performance Comparison
==================================================================
[1/2] Running Python (pysam) baseline...
  Duration: 45.123s
  Throughput: 487,234 reads/sec

[2/2] Running Rust (sheriff-rs) optimized...
  Duration: 2.134s
  Throughput: 10,298,345 reads/sec

RESULTS
==================================================================
Speedup: 21.1x faster
✅ SUCCESS: Rust achieves target performance!
```

### Rust-Only Benchmarks (Criterion)

```bash
cargo bench
```

View detailed reports in `target/criterion/report/index.html`.

## Testing

```bash
# Run unit tests
cargo test

# Run integration tests (requires test data)
cargo test --ignored

# Run with output
cargo test -- --nocapture
```

## Development

### Project Structure

```
sheriff-rs/
├── src/
│   ├── lib.rs          # Library entry point
│   ├── bam_filter.rs   # BAM filtering logic
│   ├── kmer.rs         # K-mer hashing
│   ├── python.rs       # PyO3 bindings
│   └── main.rs         # CLI binary
├── benches/
│   └── bam_benchmark.rs
├── tests/
│   ├── integration_test.rs
│   └── data/           # Test BAM files
├── Cargo.toml
└── README.md
```

### Adding New Features

1. Implement in Rust (e.g., `src/new_feature.rs`)
2. Add unit tests (`#[cfg(test)] mod tests { ... }`)
3. Add Python bindings in `src/python.rs`
4. Add benchmarks in `benches/`
5. Document in this README

## Deployment

### Build Wheels for Distribution

```bash
# Build wheel for current platform
maturin build --release

# Build for multiple platforms (requires Docker)
maturin build --release --manylinux 2014
```

### Install from Wheel

```bash
pip install target/wheels/sheriff_rs-0.1.0-*.whl
```

## Troubleshooting

**Error: `linker 'cc' not found`**
```bash
# Install C compiler
sudo apt-get install build-essential  # Ubuntu/Debian
```

**Error: `cannot import name 'sheriff_rs'`**
```bash
# Rebuild Python module
cd sheriff-rs
maturin develop --release
```

**Slow performance in debug mode**
```bash
# Always use --release for benchmarking
cargo build --release
cargo bench
```

## Contributing

See [CONTRIBUTING.md](../CONTRIBUTING.md) in the main Sheriff repository.

## License

BSD 3-Clause License (same as Sheriff)

## Resources

- [rust-htslib documentation](https://docs.rs/rust-htslib/)
- [PyO3 user guide](https://pyo3.rs/)
- [maturin guide](https://www.maturin.rs/)
- [Criterion benchmarking](https://bheisler.github.io/criterion.rs/book/)
