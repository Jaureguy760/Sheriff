# Sheriff Test Data

This directory contains test datasets for validating Sheriff's Python and Rust implementations.

## Two Test Packages Available

### 1. Single-Chromosome Test (RECOMMENDED for Claude Cloud/CI)
**Complete, self-contained, 23MB**

```bash
# Download complete package (includes reference genome!)
cd Sheriff
wget -O- -q https://github.com/BradBalderson/Sheriff/releases/latest/download/test_data_chr19.tar.gz | tar xvzf -

# Decompress reference files
gunzip test_data/chr19.fa.gz test_data/chr19.gtf.gz

# Validate everything works
python test_data/validate_chr19_test.py

# Run CI validation (Rust function tests)
python test_data/ci_validation.py
```

**Contents:**
- `test_chr19.bam` (3.7MB) - 42,194 reads on chromosome 19
- `chr19.fa.gz` (16MB) - GRCh38 chromosome 19
- `chr19.gtf.gz` (2.9MB) - Ensembl 110 annotations
- `barcodes_chr19.txt` - 4,680 cell barcodes
- `edit_sites_chr19.txt` - 4 known edit sites
- Validation scripts and MD5 checksums

### 2. Multi-Chromosome Test (Local HPC only)
**Requires external 3GB reference download**

```bash
# Download BAM and config files (31MB)
wget -O- -q https://github.com/BradBalderson/Sheriff/releases/latest/download/test_data.tar.gz | tar xvzf -

# Download full reference genome (3GB)
cd test_data && bash download_reference.sh && cd ..

# Run validation
python test_data/run_validation.py
```

**Contents:**
- `test_200kb.bam` (30MB) - 352k reads across 18 chromosomes
- `barcodes.txt` - 500 cell barcode whitelist
- `edit_sites.txt` - Multiple edit sites
- Scripts to download GRCh38 reference

## MD5 Checksum Verification

All test files have MD5 checksums in `expected_checksums.json`:

```bash
# Verify file integrity
python test_data/ci_validation.py  # Includes checksum verification
```

Single-chromosome checksums:
- `test_chr19.bam`: `4feec0ac4b8707fcd64744f4386be1a0`
- `chr19.fa.gz`: `ca4f1a9927913fe5a62585a492d13c1a`
- `chr19.gtf.gz`: `a77168a66e7fadf40bcfa171dbd6c537`

## Validation Scripts

| Script | Purpose | Runtime | Needs Reference? |
|--------|---------|---------|------------------|
| `ci_validation.py` | Rust function tests with MD5 checksums | ~1 sec | No |
| `validate_chr19_test.py` | Single-chrom data integrity test | ~2 sec | Yes (included) |
| `run_validation.py` | BAM processing smoke test | ~10 sec | No |

## Expected Runtime (Single-Chrom Test)

On typical laptop/workstation:
- **CI Validation**: ~1 second
- **Full pipeline**: ~30-60 seconds (42k reads)

## Chromosome Name Mapping

The BAM uses `hg38_19` chromosome names, but the reference uses `19`.
Sheriff automatically strips `hg38_` prefix (see `sheriff/count_t7.py:reformat_chr_name`).

## Minimum Requirements

- Python 3.10+
- ~500MB RAM
- ~100MB disk space (single-chrom test)
- pysam, numpy (installed with Sheriff)
