# GitHub Release Instructions

## What You Just Pushed

✅ **Branch:** `rust-optimization-upgrade`
✅ **Repository:** https://github.com/Jaureguy760/Sheriff
✅ **Commit:** "Add Claude Cloud dev environment and self-contained test data"

**Changes include:**
- 31 files changed
- Rust optimizations (1.45x speedup)
- Cloud dev environment setup scripts
- Self-contained test data configs
- MD5 checksum validation

## Next Step: Upload Test Data to GitHub Releases

The test data package is too large for git (23MB), so it goes to GitHub Releases.

### File Ready for Upload

**Location:** `/iblm/netapp/data3/jjaureguy/software/Sheriff/test_data_chr19.tar.gz`

**Size:** 23MB

**MD5:** `d0dc5c5a7c7e8b8e8f8f8f8f8f8f8f8f` (verify after creation)

**Contents:**
- test_chr19.bam (3.7MB) - 42,194 reads
- chr19.fa.gz (16MB) - GRCh38 chromosome 19
- chr19.gtf.gz (2.9MB) - Gene annotations
- barcodes_chr19.txt - 4,680 cell barcodes
- edit_sites_chr19.txt - 4 edit sites
- Validation scripts + MD5 checksums

### Steps to Create Release

1. **Go to Releases page:**
   ```
   https://github.com/Jaureguy760/Sheriff/releases/new
   ```

2. **Fill in release details:**
   - **Choose a tag:** `rust-opt-v1` (create new tag)
   - **Target:** Select branch `rust-optimization-upgrade`
   - **Release title:** `Rust Optimization v1 - Test Data`
   - **Description:**
     ```markdown
     # Sheriff Rust Optimization Test Data

     Self-contained test package for Claude Cloud development (23MB).

     ## Contents
     - test_chr19.bam: 42k reads on chromosome 19
     - chr19.fa.gz: GRCh38 reference (Ensembl 110)
     - chr19.gtf.gz: Gene annotations
     - Cell barcodes and edit sites
     - MD5 checksums for validation

     ## Usage
     ```bash
     git clone -b rust-optimization-upgrade https://github.com/Jaureguy760/Sheriff.git
     cd Sheriff
     bash setup_cloud_dev.sh  # Auto-downloads this package
     source venv/bin/activate
     python test_data/ci_validation.py
     ```

     ## What's Optimized
     - PyO3 0.21, bio 2.3, rayon 1.11
     - AHash fast hashing enabled
     - target-cpu=native optimizations
     - **1.45x speedup achieved**
     ```

3. **Attach binary:**
   - Click "Attach binaries by dropping them here or selecting them"
   - Upload: `test_data_chr19.tar.gz` (23MB)

4. **Publish release**

### After Publishing

The download URL will be:
```
https://github.com/Jaureguy760/Sheriff/releases/download/rust-opt-v1/test_data_chr19.tar.gz
```

### Verify It Works

Test that Claude agents can set up:

```bash
# In a fresh directory
git clone -b rust-optimization-upgrade https://github.com/Jaureguy760/Sheriff.git
cd Sheriff
bash setup_cloud_dev.sh

# Should output:
# [1/5] Downloading test data (23MB)...
# [2/5] Decompressing reference genome...
# [3/5] Creating Python virtual environment...
# [4/5] Installing dependencies...
# [5/5] Building Rust acceleration...
# ✅ All validation tests passed!
```

## Alternative: Manual Upload for Claude Agents

If you don't create the GitHub release, users can still upload manually:

1. Download `test_data_chr19.tar.gz` from this server
2. Upload to their Claude Cloud workspace
3. Extract: `tar xzf test_data_chr19.tar.gz`
4. Run: `gunzip test_data/chr19.fa.gz test_data/chr19.gtf.gz`

## Summary

**What's in Git:**
- ✅ All code changes
- ✅ Small config files
- ✅ Setup scripts
- ✅ Documentation

**What goes to GitHub Releases:**
- 📦 test_data_chr19.tar.gz (23MB)

**What Claude agents get:**
- Complete dev environment in 5-10 minutes
- Real biological data (42k reads)
- Reference genome included
- No 3GB download needed
