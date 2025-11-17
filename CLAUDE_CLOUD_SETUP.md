# Sheriff - Claude Cloud Deployment Guide

This guide explains how to set up Sheriff for development in Claude Cloud environments.

## Quick Setup for Claude Cloud

```bash
# 1. Clone the repository
git clone https://github.com/BradBalderson/Sheriff.git
cd Sheriff

# 2. Download test data (31MB)
wget -O- -q https://github.com/BradBalderson/Sheriff/releases/latest/download/test_data.tar.gz | tar xvzf -

# 3. Set up dev environment (creates venv, installs dependencies)
bash setup_dev_environment.sh

# 4. Activate environment
source venv/bin/activate

# 5. [OPTIONAL] Download reference genome (3GB - skip for quick testing)
cd test_data && bash download_reference.sh && cd ..

# 6. Validate installation
python test_data/run_validation.py
```

## What's Included

### Test Data Package (31MB)
- `test_200kb.bam` - Small BAM with ~352k reads
- `test_200kb.bam.bai` - BAM index
- `barcodes.txt` - 500 cell barcode whitelist
- `edit_sites.txt` - Known edit positions
- `blacklist.bed` - Exclusion regions
- Configuration files for realistic testing

### Development Environment
- Python 3.10+ virtual environment
- All Sheriff dependencies (pysam, numpy, pandas, scanpy)
- Maturin for Rust builds
- Sheriff installed in editable mode

### Rust Acceleration
- PyO3 0.21 bindings
- rust-bio 2.3 sequence alignment
- AHash for fast lookups
- Rayon parallelism
- CPU-native optimizations

## Testing Without Reference Genome

For quick code validation without the 3GB reference genome:

```python
# Test Rust acceleration directly
import sheriff_rs

# k-mer counting
kmers = sheriff_rs.count_kmers_rust("ATCGATCGATCG", 4)
print(f"K-mer counts: {len(kmers)} bins")

# UMI deduplication
umis = ["ATCGATCG", "ATCGATCC", "TTTTTTTT"]  # 2 groups (first two are 1 edit apart)
unique = sheriff_rs.deduplicate_umis_py(umis)
print(f"Unique UMI groups: {unique}")

# Edit clustering
edits = [
    ("chr1", 1000, "A", "ATCG", True, [1, 2]),
    ("chr1", 1000, "A", "AT", True, [1]),  # subset of above
]
longest = sheriff_rs.get_longest_edits_rust(edits)
print(f"Longest edits: {len(longest)}")
```

## Full Pipeline Testing

With reference genome downloaded:

```python
from sheriff import count_t7

# Run edit detection on test BAM
results = count_t7.count_t7_edits(
    bam_file="test_data/test_200kb.bam",
    ref_fasta_file="test_data/reference.fa",
    white_list_file="test_data/barcodes.txt",
    edit_site_list="test_data/edit_sites.txt",
    output_prefix="test_output",
)
```

## Performance Comparison

Run the benchmarking script:

```bash
python sanity_check_comparison.py
```

Expected results (on test dataset):
- Python-only: ~60 seconds
- Rust-accelerated: ~45 seconds (1.5x faster)

## Key Files for Development

```
Sheriff/
├── sheriff/                    # Python package
│   ├── count_t7.py            # Main pipeline
│   ├── helpers.py             # Core utilities
│   └── rust_accelerated.py    # Rust wrappers
├── sheriff-rs/                 # Rust acceleration
│   ├── Cargo.toml             # Dependencies
│   └── src/
│       ├── edit_clustering.rs # Edit deduplication
│       ├── umi.rs             # UMI deduplication
│       └── python.rs          # PyO3 bindings
├── test_data/                  # Test dataset
├── setup_dev_environment.sh    # Dev setup script
└── sanity_check_comparison.py  # Benchmarking
```

## Rebuilding Rust After Changes

```bash
cd sheriff-rs
maturin develop --release
cd ..
```

## Common Issues

### Missing libclang
If updating rust-htslib, stick with version 0.46:
```toml
rust-htslib = "0.46"  # Pre-built bindings, no libclang needed
```

### Rust module not found
Ensure you've built with maturin and activated the venv:
```bash
source venv/bin/activate
cd sheriff-rs && maturin develop --release && cd ..
```

### Reference genome issues
The BAM uses `hg38_*` chromosome names. Ensure your reference matches.

## GitHub Release Process

To create a new release with test data:

1. Run `bash create_test_data_release.sh`
2. Go to GitHub releases page
3. Create new release with version tag
4. Attach `test_data.tar.gz` as release asset
5. Update download URL if different from latest

## Size Considerations

- Test data package: ~31MB (fits GitHub releases)
- Reference genome: ~3GB (downloaded separately)
- Full dev environment: ~500MB (with venv)
- Compiled Rust module: ~20MB
