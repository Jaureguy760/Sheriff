# Sheriff Rust Implementation Plan

## Test-Driven Development Approach for High-Performance BAM Processing

**Goal:** Achieve 10-50x speedup for BAM iteration and filtering using Rust, with comprehensive benchmarking to validate each phase.

**Timeline:** 3-4 weeks (1 week per phase + testing)

---

## Phase 0: Prerequisites & Setup (Day 1)

### Install Rust Toolchain

```bash
# Install Rust (if not already installed)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env

# Verify installation
rustc --version
cargo --version

# Install maturin for Python-Rust integration
pip install maturin
```

### Create Rust Project Structure

```bash
# Create Rust library crate
cargo new --lib sheriff-rs
cd sheriff-rs

# Project structure:
sheriff-rs/
├── Cargo.toml              # Rust dependencies
├── src/
│   ├── lib.rs             # Main library entry
│   ├── bam_filter.rs      # BAM filtering logic
│   └── kmer.rs            # K-mer hashing (optional)
├── benches/
│   └── bam_benchmark.rs   # Criterion benchmarks
├── tests/
│   └── integration_test.rs
└── README.md
```

### Key Dependencies (Cargo.toml)

```toml
[package]
name = "sheriff-rs"
version = "0.1.0"
edition = "2021"

[dependencies]
rust-htslib = "0.46"        # BAM/SAM parsing (HTSlib wrapper)
anyhow = "1.0"              # Error handling
rayon = "1.8"               # Parallelism
ahash = "0.8"               # Fast hashing
pyo3 = { version = "0.20", features = ["extension-module"], optional = true }

[dev-dependencies]
criterion = { version = "0.5", features = ["html_reports"] }

[lib]
name = "sheriff_rs"
crate-type = ["cdylib", "rlib"]  # cdylib for Python, rlib for Rust

[features]
default = []
python = ["pyo3"]

[[bench]]
name = "bam_benchmark"
harness = false
```

---

## Phase 1: Standalone Rust BAM Filter (Week 1)

### Objective
Create a high-performance BAM filter utility that:
- Reads BAM file
- Filters reads by cell barcode whitelist
- Outputs filtered BAM
- **10-50x faster than Python pysam iteration**

### 1.1: Write Tests FIRST (Test-Driven Development)

**File:** `sheriff-rs/tests/integration_test.rs`

```rust
use sheriff_rs::bam_filter::filter_bam_by_barcodes;
use std::collections::HashSet;

#[test]
fn test_filter_bam_basic() {
    // Create test data
    let input_bam = "tests/data/test_input.bam";
    let output_bam = "tests/data/test_output.bam";

    // Whitelist barcodes
    let mut whitelist = HashSet::new();
    whitelist.insert("AAACCTGAGAAACCAT-1".to_string());
    whitelist.insert("AAACCTGAGAAACCGC-1".to_string());

    // Filter BAM
    let result = filter_bam_by_barcodes(input_bam, output_bam, &whitelist);

    assert!(result.is_ok());
    assert_eq!(result.unwrap().reads_processed, 1000);
    assert_eq!(result.unwrap().reads_kept, 50);
}

#[test]
fn test_filter_bam_performance() {
    // Performance regression test
    let input_bam = "tests/data/large_test.bam";
    let output_bam = "tests/data/large_output.bam";

    let whitelist = load_whitelist("tests/data/whitelist.txt");

    let start = std::time::Instant::now();
    let result = filter_bam_by_barcodes(input_bam, output_bam, &whitelist).unwrap();
    let duration = start.elapsed();

    println!("Processed {} reads in {:?}", result.reads_processed, duration);
    println!("Throughput: {} reads/sec", result.reads_processed as f64 / duration.as_secs_f64());

    // Regression check: should process at least 100k reads/sec
    assert!(result.reads_processed as f64 / duration.as_secs_f64() > 100_000.0);
}
```

**Create test data:**
```bash
# Use samtools to create small test BAM from example data
cd sheriff-rs
mkdir -p tests/data
samtools view -h -s 0.01 ../example_data/barcode_headAligned_anno.sorted.edit_regions_200kb.bam | samtools view -b > tests/data/test_input.bam
samtools index tests/data/test_input.bam
```

### 1.2: Implement BAM Filter (TDD Style)

**File:** `sheriff-rs/src/bam_filter.rs`

```rust
use rust_htslib::bam::{self, Read, Record};
use std::collections::HashSet;
use anyhow::{Result, Context};

pub struct FilterResult {
    pub reads_processed: usize,
    pub reads_kept: usize,
    pub reads_rejected: usize,
}

/// Filter BAM file by cell barcode whitelist
pub fn filter_bam_by_barcodes(
    input_path: &str,
    output_path: &str,
    whitelist: &HashSet<String>,
) -> Result<FilterResult> {
    // Open input BAM
    let mut bam = bam::Reader::from_path(input_path)
        .context("Failed to open input BAM")?;

    // Create output BAM with same header
    let header = bam::Header::from_template(bam.header());
    let mut out = bam::Writer::from_path(output_path, &header, bam::Format::Bam)
        .context("Failed to create output BAM")?;

    let mut stats = FilterResult {
        reads_processed: 0,
        reads_kept: 0,
        reads_rejected: 0,
    };

    // Iterate reads
    for result in bam.records() {
        let record = result.context("Failed to read BAM record")?;
        stats.reads_processed += 1;

        // Extract CB tag (cell barcode)
        if let Some(cb_tag) = get_cb_tag(&record) {
            if whitelist.contains(&cb_tag) {
                out.write(&record).context("Failed to write record")?;
                stats.reads_kept += 1;
            } else {
                stats.reads_rejected += 1;
            }
        } else {
            // No CB tag, reject
            stats.reads_rejected += 1;
        }
    }

    Ok(stats)
}

/// Extract CB (cell barcode) tag from BAM record
fn get_cb_tag(record: &Record) -> Option<String> {
    record.aux(b"CB")
        .ok()
        .and_then(|aux| aux.string().ok())
        .map(|s| s.to_string())
}

/// Load whitelist from text file
pub fn load_whitelist(path: &str) -> Result<HashSet<String>> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    let file = File::open(path).context("Failed to open whitelist")?;
    let reader = BufReader::new(file);

    let mut whitelist = HashSet::new();
    for line in reader.lines() {
        let barcode = line.context("Failed to read line")?;
        whitelist.insert(barcode.trim().to_string());
    }

    Ok(whitelist)
}
```

**File:** `sheriff-rs/src/lib.rs`

```rust
pub mod bam_filter;
pub mod kmer;

// Re-export main functions
pub use bam_filter::{filter_bam_by_barcodes, load_whitelist, FilterResult};

// Python bindings (Phase 2)
#[cfg(feature = "python")]
pub mod python;
```

### 1.3: Benchmark Against Python (Criterion)

**File:** `sheriff-rs/benches/bam_benchmark.rs`

```rust
use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use sheriff_rs::bam_filter::{filter_bam_by_barcodes, load_whitelist};
use std::collections::HashSet;

fn bench_bam_filter(c: &mut Criterion) {
    let mut group = c.benchmark_group("bam_filter");

    // Load whitelist once
    let whitelist = load_whitelist("tests/data/whitelist.txt").unwrap();

    group.bench_function("filter_small_bam", |b| {
        b.iter(|| {
            filter_bam_by_barcodes(
                black_box("tests/data/test_input.bam"),
                black_box("/tmp/test_output.bam"),
                black_box(&whitelist)
            ).unwrap()
        });
    });

    group.finish();
}

criterion_group!(benches, bench_bam_filter);
criterion_main!(benches);
```

**Run benchmarks:**
```bash
cd sheriff-rs
cargo bench

# Output will show:
# bam_filter/filter_small_bam
#                         time:   [12.345 ms 12.567 ms 12.789 ms]
#                         thrpt: [78,234 reads/s 79,567 reads/s 80,890 reads/s]
```

### 1.4: Compare with Python Baseline

**File:** `benchmarks/compare_rust_python.py`

```python
#!/usr/bin/env env python3
"""
Compare Rust BAM filter vs Python pysam iteration performance
"""
import time
import pysam
import subprocess
import json

def benchmark_python_filter(bam_path, whitelist_path, output_path):
    """Python baseline using pysam"""
    # Load whitelist
    with open(whitelist_path) as f:
        whitelist = set(line.strip() for line in f)

    # Process BAM
    start = time.time()
    bam = pysam.AlignmentFile(bam_path, 'rb')
    out = pysam.AlignmentFile(output_path, 'wb', template=bam)

    reads_processed = 0
    reads_kept = 0

    for read in bam:
        reads_processed += 1
        try:
            cb = read.get_tag('CB')
            if cb in whitelist:
                out.write(read)
                reads_kept += 1
        except KeyError:
            pass

    duration = time.time() - start

    bam.close()
    out.close()

    return {
        'reads_processed': reads_processed,
        'reads_kept': reads_kept,
        'duration_seconds': duration,
        'throughput_reads_per_sec': reads_processed / duration
    }

def benchmark_rust_filter(bam_path, whitelist_path, output_path):
    """Rust implementation using sheriff-rs CLI"""
    cmd = [
        'cargo', 'run', '--release', '--',
        'filter-bam',
        '--input', bam_path,
        '--output', output_path,
        '--whitelist', whitelist_path
    ]

    start = time.time()
    result = subprocess.run(cmd, cwd='sheriff-rs', capture_output=True, text=True)
    duration = time.time() - start

    # Parse output (JSON from Rust)
    stats = json.loads(result.stdout)

    return {
        **stats,
        'duration_seconds': duration,
        'throughput_reads_per_sec': stats['reads_processed'] / duration
    }

if __name__ == '__main__':
    test_bam = 'sheriff-rs/tests/data/test_input.bam'
    whitelist = 'example_data/barcode_whitelist.500-cell.txt'

    print("=" * 60)
    print("BAM Filter Performance Comparison")
    print("=" * 60)

    # Python baseline
    print("\n[1/2] Running Python (pysam) baseline...")
    python_results = benchmark_python_filter(
        test_bam, whitelist, '/tmp/python_output.bam'
    )

    print(f"  Reads processed: {python_results['reads_processed']:,}")
    print(f"  Reads kept: {python_results['reads_kept']:,}")
    print(f"  Duration: {python_results['duration_seconds']:.3f}s")
    print(f"  Throughput: {python_results['throughput_reads_per_sec']:,.0f} reads/sec")

    # Rust optimized
    print("\n[2/2] Running Rust (sheriff-rs) optimized...")
    rust_results = benchmark_rust_filter(
        test_bam, whitelist, '/tmp/rust_output.bam'
    )

    print(f"  Reads processed: {rust_results['reads_processed']:,}")
    print(f"  Reads kept: {rust_results['reads_kept']:,}")
    print(f"  Duration: {rust_results['duration_seconds']:.3f}s")
    print(f"  Throughput: {rust_results['throughput_reads_per_sec']:,.0f} reads/sec")

    # Comparison
    speedup = python_results['duration_seconds'] / rust_results['duration_seconds']

    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)
    print(f"Speedup: {speedup:.1f}x faster")
    print(f"Target: 10-50x faster")

    if speedup >= 10:
        print("✅ SUCCESS: Rust achieves target performance!")
    else:
        print(f"⚠️  WARNING: Only {speedup:.1f}x speedup (target: 10x)")

    # Save results
    results = {
        'python': python_results,
        'rust': rust_results,
        'speedup': speedup
    }

    with open('benchmarks/rust_python_comparison.json', 'w') as f:
        json.dump(results, f, indent=2)

    print("\nResults saved to: benchmarks/rust_python_comparison.json")
```

### 1.5: Create Standalone CLI (for testing)

**File:** `sheriff-rs/src/main.rs`

```rust
use sheriff_rs::bam_filter::{filter_bam_by_barcodes, load_whitelist};
use clap::Parser;
use anyhow::Result;

#[derive(Parser)]
#[command(name = "sheriff-rs")]
#[command(about = "High-performance BAM filtering for Sheriff")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Parser)]
enum Commands {
    FilterBam {
        #[arg(short, long)]
        input: String,

        #[arg(short, long)]
        output: String,

        #[arg(short, long)]
        whitelist: String,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::FilterBam { input, output, whitelist } => {
            eprintln!("Loading whitelist from: {}", whitelist);
            let barcodes = load_whitelist(&whitelist)?;
            eprintln!("Loaded {} barcodes", barcodes.len());

            eprintln!("Filtering BAM: {} -> {}", input, output);
            let result = filter_bam_by_barcodes(&input, &output, &barcodes)?;

            // Output JSON for benchmarking
            let json = serde_json::json!({
                "reads_processed": result.reads_processed,
                "reads_kept": result.reads_kept,
                "reads_rejected": result.reads_rejected,
            });
            println!("{}", json);

            eprintln!("✅ Complete: {} reads processed, {} kept",
                     result.reads_processed, result.reads_kept);
        }
    }

    Ok(())
}
```

**Add to Cargo.toml:**
```toml
[dependencies]
clap = { version = "4.4", features = ["derive"] }
serde_json = "1.0"

[[bin]]
name = "sheriff-rs"
path = "src/main.rs"
```

**Test CLI:**
```bash
cargo build --release
./target/release/sheriff-rs filter-bam \
    --input tests/data/test_input.bam \
    --output /tmp/filtered.bam \
    --whitelist tests/data/whitelist.txt
```

---

## Phase 2: Python Integration via PyO3 (Week 2)

### Objective
Create Python bindings so Sheriff can call Rust functions directly, with automatic fallback to Python if Rust unavailable.

### 2.1: Add PyO3 Bindings

**File:** `sheriff-rs/src/python.rs`

```rust
use pyo3::prelude::*;
use pyo3::types::PyDict;
use crate::bam_filter::{filter_bam_by_barcodes, load_whitelist};
use std::collections::HashSet;

#[pyfunction]
fn filter_bam_rust(
    input_path: String,
    output_path: String,
    whitelist_path: String,
) -> PyResult<PyObject> {
    // Load whitelist
    let whitelist = load_whitelist(&whitelist_path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    // Filter BAM
    let result = filter_bam_by_barcodes(&input_path, &output_path, &whitelist)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    // Return dict
    Python::with_gil(|py| {
        let dict = PyDict::new(py);
        dict.set_item("reads_processed", result.reads_processed)?;
        dict.set_item("reads_kept", result.reads_kept)?;
        dict.set_item("reads_rejected", result.reads_rejected)?;
        Ok(dict.into())
    })
}

#[pymodule]
fn sheriff_rs(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(filter_bam_rust, m)?)?;
    Ok(())
}
```

**Build Python module:**
```bash
cd sheriff-rs
maturin develop --release
```

**Test in Python:**
```python
import sheriff_rs

result = sheriff_rs.filter_bam_rust(
    "tests/data/test_input.bam",
    "/tmp/filtered.bam",
    "tests/data/whitelist.txt"
)

print(f"Processed {result['reads_processed']} reads")
print(f"Kept {result['reads_kept']} reads")
```

### 2.2: Integrate into Sheriff with Fallback

**File:** `sheriff/bam_utils.py` (NEW)

```python
"""
High-performance BAM utilities with optional Rust acceleration
"""
import pysam
from typing import Set, Dict

# Try to import Rust acceleration
try:
    import sheriff_rs
    HAS_RUST = True
except ImportError:
    HAS_RUST = False
    import warnings
    warnings.warn(
        "Rust acceleration not available. "
        "Install with: cd sheriff-rs && maturin develop --release"
    )

def filter_bam_by_barcodes(
    input_bam: str,
    output_bam: str,
    cell_barcodes: Set[str],
    use_rust: bool = True
) -> Dict[str, int]:
    """
    Filter BAM file by cell barcode whitelist.

    Uses Rust implementation if available (10-50x faster),
    falls back to Python pysam if not.

    Args:
        input_bam: Input BAM file path
        output_bam: Output BAM file path
        cell_barcodes: Set of allowed cell barcodes
        use_rust: Prefer Rust if available (default True)

    Returns:
        Dict with keys: reads_processed, reads_kept, reads_rejected
    """
    if HAS_RUST and use_rust:
        # Write whitelist to temp file for Rust
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            for bc in cell_barcodes:
                f.write(f"{bc}\n")
            whitelist_path = f.name

        try:
            result = sheriff_rs.filter_bam_rust(
                input_bam, output_bam, whitelist_path
            )
            return result
        finally:
            import os
            os.unlink(whitelist_path)

    else:
        # Fallback to Python implementation
        return _filter_bam_python(input_bam, output_bam, cell_barcodes)

def _filter_bam_python(
    input_bam: str,
    output_bam: str,
    cell_barcodes: Set[str]
) -> Dict[str, int]:
    """Python fallback implementation using pysam"""
    bam = pysam.AlignmentFile(input_bam, 'rb')
    out = pysam.AlignmentFile(output_bam, 'wb', template=bam)

    stats = {
        'reads_processed': 0,
        'reads_kept': 0,
        'reads_rejected': 0
    }

    for read in bam:
        stats['reads_processed'] += 1

        try:
            cb = read.get_tag('CB')
            if cb in cell_barcodes:
                out.write(read)
                stats['reads_kept'] += 1
            else:
                stats['reads_rejected'] += 1
        except KeyError:
            stats['reads_rejected'] += 1

    bam.close()
    out.close()

    return stats
```

### 2.3: Update count_t7.py to Use BAM Filter

**File:** `sheriff/count_t7.py` (MODIFY)

```python
# Add at top
from .bam_utils import filter_bam_by_barcodes

# In run_count_t7 function, BEFORE main BAM iteration:
def run_count_t7(...):
    # ... existing setup code ...

    # Load cell barcodes
    cell_barcodes = set()
    with open(cell_barcode_whitelist, 'r') as f:
        for line in f:
            cell_barcodes.add(line.strip())

    # OPTIMIZATION: Pre-filter BAM by cell barcodes (10-50x faster with Rust)
    import tempfile
    filtered_bam_path = tempfile.mktemp(suffix='.filtered.bam')

    print(f"Pre-filtering BAM by {len(cell_barcodes)} cell barcodes...")
    filter_stats = filter_bam_by_barcodes(
        bam_filepath,
        filtered_bam_path,
        cell_barcodes,
        use_rust=True  # Will fall back to Python if Rust unavailable
    )

    print(f"  Kept {filter_stats['reads_kept']:,} / {filter_stats['reads_processed']:,} reads")
    print(f"  Reduction: {100 * (1 - filter_stats['reads_kept']/filter_stats['reads_processed']):.1f}%")

    # Now iterate filtered BAM (50-90% fewer reads!)
    bam = pysam.AlignmentFile(filtered_bam_path, 'rb')

    # ... rest of existing code ...

    # Cleanup
    import os
    os.unlink(filtered_bam_path)
```

### 2.4: Benchmark Full Pipeline

**File:** `benchmarks/benchmark_full_pipeline_rust.py`

```python
#!/usr/bin/env python3
"""
Benchmark full Sheriff pipeline with and without Rust acceleration
"""
import time
import subprocess
import json

def run_sheriff(use_rust: bool) -> dict:
    """Run Sheriff with or without Rust"""
    env = os.environ.copy()
    if not use_rust:
        env['SHERIFF_NO_RUST'] = '1'  # Disable Rust

    start = time.time()

    result = subprocess.run([
        'sheriff',
        'example_data/barcode_headAligned_anno.sorted.edit_regions_200kb.bam',
        'example_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa',
        'example_data/barcode_whitelist.500-cell.txt',
        'example_data/Homo_sapiens.GRCh38.110.gtf',
        '-o', f'/tmp/sheriff_output_rust_{use_rust}'
    ], env=env, capture_output=True, text=True)

    duration = time.time() - start

    return {
        'duration': duration,
        'stdout': result.stdout,
        'stderr': result.stderr
    }

if __name__ == '__main__':
    print("Benchmarking Full Sheriff Pipeline")
    print("=" * 60)

    print("\n[1/2] Running WITHOUT Rust...")
    python_result = run_sheriff(use_rust=False)
    print(f"  Duration: {python_result['duration']:.2f}s")

    print("\n[2/2] Running WITH Rust...")
    rust_result = run_sheriff(use_rust=True)
    print(f"  Duration: {rust_result['duration']:.2f}s")

    speedup = python_result['duration'] / rust_result['duration']

    print("\n" + "=" * 60)
    print(f"Speedup: {speedup:.1f}x faster with Rust")
    print("=" * 60)

    results = {
        'python_only': python_result['duration'],
        'with_rust': rust_result['duration'],
        'speedup': speedup
    }

    with open('benchmarks/full_pipeline_rust_comparison.json', 'w') as f:
        json.dump(results, f, indent=2)
```

---

## Phase 3: Additional Optimizations (Week 3-4, Optional)

### Objective
If Phase 1-2 doesn't achieve 10x speedup, optimize additional hot paths.

### 3.1: Rust K-mer Hashing (if needed)

**File:** `sheriff-rs/src/kmer.rs`

```rust
use std::collections::HashMap;
use ahash::AHashMap;  // Faster than std HashMap

/// Fast k-mer to integer hashing
pub fn kmer_to_num(kmer: &[u8]) -> u32 {
    let mut result = 0u32;
    for &nucleotide in kmer {
        result = result.wrapping_shl(2);
        result += match nucleotide {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 0,  // Treat unknown as A
        };
    }
    result
}

/// Count k-mers in sequence
pub fn count_kmers(sequence: &[u8], k: usize) -> Vec<u8> {
    let array_size = 4usize.pow(k as u32);
    let mut counts = vec![0u8; array_size];

    if sequence.len() < k {
        return counts;
    }

    for window in sequence.windows(k) {
        let hash = kmer_to_num(window) as usize;
        counts[hash] = counts[hash].saturating_add(1);
    }

    counts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_to_num() {
        assert_eq!(kmer_to_num(b"AA"), 0);
        assert_eq!(kmer_to_num(b"AC"), 1);
        assert_eq!(kmer_to_num(b"AT"), 3);
        assert_eq!(kmer_to_num(b"CA"), 4);
    }

    #[test]
    fn test_count_kmers() {
        let seq = b"AACGTT";
        let counts = count_kmers(seq, 2);

        // AA appears once
        assert_eq!(counts[kmer_to_num(b"AA") as usize], 1);
        // AC appears once
        assert_eq!(counts[kmer_to_num(b"AC") as usize], 1);
    }
}
```

**PyO3 binding:**
```rust
// In python.rs
#[pyfunction]
fn count_kmers_rust(sequence: String, k: usize) -> Vec<u8> {
    crate::kmer::count_kmers(sequence.as_bytes(), k)
}
```

### 3.2: Parallel BAM Processing with Rayon

**File:** `sheriff-rs/src/bam_filter.rs` (enhanced)

```rust
use rayon::prelude::*;

/// Parallel BAM filter (for very large BAMs)
pub fn filter_bam_parallel(
    input_path: &str,
    output_path: &str,
    whitelist: &HashSet<String>,
    threads: usize,
) -> Result<FilterResult> {
    // Split BAM by chromosome, process in parallel
    // Merge results
    // (Advanced - only if needed)
    todo!("Implement parallel processing if single-threaded isn't fast enough")
}
```

---

## Testing & Validation Strategy

### Unit Tests (Rust)

```bash
# Run all tests
cargo test

# Run with output
cargo test -- --nocapture

# Test specific module
cargo test bam_filter
```

### Integration Tests (Python)

**File:** `tests/test_rust_integration.py`

```python
import pytest
import os
import pysam

try:
    import sheriff_rs
    HAS_RUST = True
except ImportError:
    HAS_RUST = False

@pytest.mark.skipif(not HAS_RUST, reason="Rust not available")
def test_rust_bam_filter():
    """Test Rust BAM filter produces same output as Python"""
    from sheriff.bam_utils import filter_bam_by_barcodes

    input_bam = "example_data/barcode_headAligned_anno.sorted.edit_regions_200kb.bam"
    whitelist_path = "example_data/barcode_whitelist.500-cell.txt"

    # Load whitelist
    with open(whitelist_path) as f:
        cell_barcodes = set(line.strip() for line in f)

    # Run Rust version
    rust_output = "/tmp/test_rust_output.bam"
    rust_stats = filter_bam_by_barcodes(
        input_bam, rust_output, cell_barcodes, use_rust=True
    )

    # Run Python version
    python_output = "/tmp/test_python_output.bam"
    python_stats = filter_bam_by_barcodes(
        input_bam, python_output, cell_barcodes, use_rust=False
    )

    # Compare stats
    assert rust_stats['reads_processed'] == python_stats['reads_processed']
    assert rust_stats['reads_kept'] == python_stats['reads_kept']

    # Compare output BAMs
    rust_bam = pysam.AlignmentFile(rust_output, 'rb')
    python_bam = pysam.AlignmentFile(python_output, 'rb')

    rust_reads = list(rust_bam)
    python_reads = list(python_bam)

    assert len(rust_reads) == len(python_reads), "Different number of reads in output"

    # Cleanup
    os.unlink(rust_output)
    os.unlink(python_output)
```

### Benchmark Validation

**Checklist before declaring success:**
- [ ] Rust filter produces identical output to Python (byte-for-byte BAM comparison)
- [ ] Rust achieves ≥10x speedup on example data
- [ ] Memory usage doesn't increase significantly
- [ ] Works with conda/pip installation
- [ ] Gracefully falls back to Python if Rust unavailable
- [ ] Documentation updated

---

## Installation & Distribution

### Development Install

```bash
# Install Sheriff in editable mode
pip install -e .

# Build Rust extension
cd sheriff-rs
maturin develop --release
cd ..

# Verify Rust is available
python -c "import sheriff_rs; print('Rust available!')"
```

### Production Install (with Rust)

**Option 1: maturin wheels**
```bash
# Build wheel with Rust included
cd sheriff-rs
maturin build --release

# Install wheel
pip install target/wheels/sheriff_rs-0.1.0-*.whl
```

**Option 2: Add to setup.py**
```python
# setup.py
from setuptools import setup
from setuptools_rust import Binding, RustExtension

setup(
    # ... existing setup ...
    rust_extensions=[
        RustExtension(
            "sheriff.sheriff_rs",
            path="sheriff-rs/Cargo.toml",
            binding=Binding.PyO3,
            optional=True  # Don't fail if Rust not available
        )
    ],
)
```

### Conda Distribution

**File:** `.conda/recipe/meta.yaml`

```yaml
package:
  name: sheriff
  version: {{ GIT_DESCRIBE_TAG }}

build:
  number: 0
  script:
    - {{ PYTHON }} -m pip install . --no-deps -vv
    # Build Rust if available
    - cd sheriff-rs && maturin build --release || true

requirements:
  build:
    - rust >=1.70  # Optional
    - maturin      # Optional
  host:
    - python
    - pip
  run:
    - python >=3.10
    - pysam >=0.19.0
    # ... other deps ...
```

---

## Performance Targets & Success Criteria

### Phase 1 Success Criteria
- [ ] Rust BAM filter implemented and tested
- [ ] Achieves **10-50x speedup** vs Python pysam on example data
- [ ] Produces identical output (validated with `samtools diff`)
- [ ] Benchmark comparison saved to `benchmarks/rust_python_comparison.json`

### Phase 2 Success Criteria
- [ ] PyO3 bindings working from Python
- [ ] Automatic fallback to Python if Rust unavailable
- [ ] Full Sheriff pipeline uses Rust filter
- [ ] Overall pipeline **3-10x faster** than baseline
- [ ] No regression in output quality

### Phase 3 Success Criteria (if needed)
- [ ] Additional hot paths identified via profiling
- [ ] K-mer hashing in Rust (if needed)
- [ ] Parallel processing (if needed)
- [ ] Final speedup: **10-50x overall**

---

## Troubleshooting

### "error: linker `cc` not found"
```bash
# Install C compiler (required for Rust)
sudo apt-get install build-essential  # Ubuntu/Debian
```

### "ImportError: cannot import name 'sheriff_rs'"
```bash
# Rebuild Rust extension
cd sheriff-rs
maturin develop --release
```

### "Rust version slower than expected"
- Check release build: `cargo build --release` (not debug)
- Profile with `cargo flamegraph`
- Ensure using ahash instead of std HashMap

### "BAM outputs differ between Rust and Python"
- Check BAM header preservation
- Verify tag encoding (CB:Z: vs CB:A:)
- Use `samtools diff` to find differences

---

## Resources & References

### Rust for Bioinformatics
- [rust-bio cookbook](https://rust-bio.github.io/)
- [rust-htslib documentation](https://docs.rs/rust-htslib/)
- [PyO3 user guide](https://pyo3.rs/)
- [maturin guide](https://www.maturin.rs/)

### Performance Profiling
- [Criterion.rs benchmarking](https://bheisler.github.io/criterion.rs/book/)
- [cargo-flamegraph](https://github.com/flamegraph-rs/flamegraph)
- [perf profiling](https://perf.wiki.kernel.org/)

### Similar Projects
- [pyo3-pack examples](https://github.com/PyO3/pyo3-pack)
- [polars](https://github.com/pola-rs/polars) (Rust + Python)
- [cryptography](https://github.com/pyca/cryptography) (Rust + Python)

---

## Next Steps

1. **Day 1:** Set up Rust project, write tests
2. **Day 2-3:** Implement BAM filter in Rust
3. **Day 4:** Benchmark and validate
4. **Week 2:** PyO3 integration
5. **Week 3:** Full pipeline integration and testing
6. **Week 4:** Documentation, distribution, optimization

**Let's start with Phase 0!** 🚀

Would you like me to:
1. Create the initial Rust project structure?
2. Generate test data from example BAM?
3. Write the first integration test?
