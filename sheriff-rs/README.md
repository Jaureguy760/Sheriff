# Sheriff-rs

A high-performance Rust implementation of bioinformatics algorithms for k-mer analysis, UMI (Unique Molecular Identifier) deduplication, and BAM file processing.

## Features

- **K-mer Analysis**: Efficient k-mer extraction and counting using FxHashMap for fast hashing
- **UMI Deduplication**: Fast-hash-based unique identifier tracking and deduplication
- **BAM Processing**: Binary alignment map record representation and processing utilities
- **Parallel Processing**: Leverages Rayon for multi-threaded operations
- **Python Bindings**: Optional PyO3 bindings for seamless Python integration
- **Performance**: Optimized Rust implementation with benchmarks using Criterion

## Building

### Standard Build

```bash
cargo build --release
```

### With Python Support

```bash
cargo build --release --features python
```

### Running Tests

```bash
cargo test
```

### Running Benchmarks

```bash
cargo bench
```

## Project Structure

- `src/lib.rs` - Main library interface
- `src/kmer.rs` - K-mer counter and related utilities
- `src/umi.rs` - UMI deduplicator and related utilities
- `src/bam.rs` - BAM record handling
- `src/python.rs` - PyO3 Python bindings (feature-gated)
- `benches/` - Criterion benchmarks

## Dependencies

- **rustc-hash 2.0** - High-performance FxHashMap/FxHashSet
- **rayon 1.8** - Data-level parallelism
- **pyo3 0.21** - Python bindings (optional)
- **criterion 0.5** - Benchmarking framework (dev-dependency)

## Features

- `python` - Enable PyO3 Python bindings for the library

## Usage Examples

### K-mer Counter

```rust
use sheriff_rs::KmerCounter;

let mut counter = KmerCounter::new(21);  // 21-mer counter
assert_eq!(counter.k(), 21);
assert!(counter.is_empty());
```

### UMI Deduplicator

```rust
use sheriff_rs::UmiDeduplicator;

let mut dedup = UmiDeduplicator::new(12);  // 12 base UMI
dedup.add("ACGTACGTACGT".to_string());
assert_eq!(dedup.unique_count(), 1);
assert!(dedup.contains("ACGTACGTACGT"));
```

### BAM Record

```rust
use sheriff_rs::BamRecord;

let record = BamRecord::new(
    "read1".to_string(),
    "ACGTACGTACGT".to_string(),
    "IIIIIIIIIIII".to_string(),
);
assert_eq!(record.seq_len(), 12);
```

## Python Usage

When compiled with the `python` feature, the library can be used from Python:

```python
import sheriff_rs

# K-mer counter
kmer_counter = sheriff_rs.PyKmerCounter(21)
print(f"K-mer length: {kmer_counter.k()}")

# UMI deduplicator
umi_dedup = sheriff_rs.PyUmiDeduplicator(12)
umi_dedup.add("ACGTACGTACGT")
print(f"Unique UMIs: {umi_dedup.unique_count()}")
```

## Performance

The library is optimized for performance with:

- Fast hashing using FxHashMap/FxHashSet
- Zero-copy operations where possible
- Multi-threaded parallelization support via Rayon
- Criterion-based benchmarking infrastructure

## License

See LICENSE.md in the parent directory.
