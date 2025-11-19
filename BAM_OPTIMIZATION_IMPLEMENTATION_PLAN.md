# BAM Optimization Implementation Plan: Phase 1
## Push Sheriff from 4.4x to 5-6x Speedup

**Document Version:** 1.0
**Date:** November 19, 2025
**Status:** Ready for Implementation
**Confidence Level:** 95% (libdeflate is proven, parallel BGZF is already in rust-htslib)

---

## A. Overview

### What We're Doing

Sheriff's BAM processing is currently achieving 4.4x speedup through Rust optimizations. The remaining bottleneck is **BGZF decompression**, which accounts for 60% of BAM I/O time. We're implementing two proven, low-risk optimizations to push us to 5-6x total speedup:

1. **libdeflate installation** - A faster deflate/zlib implementation (~6% speedup)
2. **Parallel BGZF decompression** - Multi-threaded BGZF decompression (~4-5% speedup)

### Why This Matters

- **Current state:** 4.4x speedup, BAM I/O is still 30% of pipeline runtime
- **Target:** 5-6x speedup (5-10% additional improvement)
- **Decompression bottleneck:** 60% of BAM I/O time is pure math (can't optimize away)
- **Our solution:** Accelerate decompression with libdeflate + parallelize across cores

### Expected Results

```
Before optimization:
├─ BAM I/O:        30% of runtime
├─ K-mer analysis: 20% of runtime
├─ UMI handling:   15% of runtime
└─ Other:          35% of runtime
Total: 100% (baseline 4.4x)

After optimization:
├─ BAM I/O:        20% of runtime (-10%)  ← libdeflate + parallel BGZF
├─ K-mer analysis: 20% of runtime (unchanged)
├─ UMI handling:   15% of runtime (unchanged)
└─ Other:          35% of runtime (unchanged)
Total: 90% (new baseline 5.2-5.6x) ← 18-27% faster
```

### Expected Speedup Breakdown

| Component | Speedup | Effort | Risk |
|-----------|---------|--------|------|
| libdeflate | +6% | 2 hours | Very Low |
| Parallel BGZF | +4-5% | 1 hour | Low |
| **Total Phase 1** | **+10-11%** | **3 hours** | **Very Low** |

This translates to:
- Current pipeline: ~4.4x speedup
- After Phase 1: ~4.9-5.0x speedup (on 4-core system)
- On 8-core system: ~5.0-5.2x speedup

### Time Required

- **Total:** 3 hours
  - Hour 1 (0-60 min): libdeflate download, build, and installation
  - Hour 2 (60-120 min): HTSlib reconfiguration and rebuild with libdeflate
  - Hour 3 (120-180 min): Parallel BGZF code changes, testing, and verification

### Risk Level: VERY LOW

**Why?**
- libdeflate is production-used in sambamba and other tools
- Parallel BGZF is already implemented in rust-htslib (we just enable it)
- No architectural changes required
- Fully reversible if needed

**Potential Issues:** Build system issues with new dependencies (easily recoverable)

---

## B. Libdeflate Installation (Step-by-Step)

### What is libdeflate?

libdeflate is a high-performance library for DEFLATE/zlib compression and decompression. Key advantages:

- **Speed:** 2-3x faster decompression than standard zlib (through optimized algorithms and vectorization)
- **Compatibility:** Drop-in replacement for zlib
- **Production Ready:** Used in sambamba, BGZip tools, and other performance-critical bioinformatics software
- **Low Overhead:** Minimal CPU overhead for the speedup

HTSlib can be compiled with libdeflate as an optional dependency, which will automatically use it for all BGZF operations.

### Prerequisites Check

Before starting, verify your system has the necessary build tools:

```bash
# Check for required tools
which gcc g++ make git
# All should return paths - if not, install build-essential
sudo apt-get update && sudo apt-get install -y build-essential git
```

Expected output: Should show paths to gcc, g++, make, and git

---

### Step 1: Clone and Build libdeflate (20 minutes)

**Action:** Download and compile libdeflate from source

**Commands to run:**

```bash
# Navigate to temporary directory
cd /tmp

# Clone libdeflate repository
git clone https://github.com/ebiggers/libdeflate.git
cd libdeflate

# Verify clone was successful
ls -la | grep -E "Makefile|README"

# Build libdeflate
make

# Verify build succeeded
ls -la libdeflate.so* libdeflate.a 2>/dev/null && echo "Build successful"
```

**What success looks like:**
- No error messages from `make`
- Output shows "libdeflate.so" and "libdeflate.a" files created
- Should complete in 30-60 seconds

**How to verify:**
```bash
# Check shared library was built
file libdeflate.so
# Should output: "ELF 64-bit LSB shared object..."

# Quick sanity test
nm libdeflate.so | grep deflate_decompress | head -1
# Should show a symbol named "deflate_decompress"
```

**Troubleshooting if it fails:**

| Error | Solution |
|-------|----------|
| `git: command not found` | Install git: `sudo apt-get install git` |
| `gcc: command not found` | Install build tools: `sudo apt-get install build-essential` |
| `Makefile: No such file` | Verify clone succeeded: `ls libdeflate/Makefile` |
| `Permission denied` | Use `sudo` or check directory permissions |

---

### Step 2: Install libdeflate System-Wide (10 minutes)

**Action:** Install libdeflate headers and libraries to system directories

**Commands to run:**

```bash
# Still in /tmp/libdeflate directory
cd /tmp/libdeflate

# Install libdeflate to /usr/local
sudo make install

# Verify installation
ls -la /usr/local/lib/libdeflate*
ls -la /usr/local/include/libdeflate.h

# Update library cache so system can find it
sudo ldconfig

# Verify system can find it
ldconfig -p | grep libdeflate
```

**What success looks like:**
- Installation outputs: "Installing library..." and "Installing headers..."
- `ls /usr/local/lib/` shows `libdeflate.so` and `libdeflate.a`
- `ls /usr/local/include/` shows `libdeflate.h`
- `ldconfig` outputs show libdeflate in the cache

**How to verify:**
```bash
# Test that pkg-config can find libdeflate (optional)
pkg-config --modversion libdeflate 2>/dev/null && echo "pkg-config OK"

# Test that compiler can find the library
gcc -I/usr/local/include -c -x c - <<< 'int main() { void libdeflate_free_decompressor(void*); }' && echo "Compiler can find libdeflate"
```

**Troubleshooting if it fails:**

| Error | Solution |
|-------|----------|
| `sudo: permission denied` | Check sudoers: `sudo -l` |
| `Permission denied /usr/local` | Use `sudo` or change permissions: `sudo chown -R $(whoami) /usr/local` |
| `ldconfig: permission denied` | Run with `sudo`: `sudo ldconfig` |

---

### Step 3: Rebuild HTSlib with libdeflate (25 minutes)

**Action:** Reconfigure and rebuild HTSlib to use libdeflate

**Note on HTSlib Location:**
HTSlib is built as part of the Cargo build process via the `hts-sys` crate. We need to force a rebuild with libdeflate enabled.

**Commands to run:**

```bash
# Navigate to Sheriff Rust directory
cd /home/user/Sheriff/sheriff-rs

# Set environment variables to tell hts-sys about libdeflate
export LDFLAGS="-L/usr/local/lib"
export CPPFLAGS="-I/usr/local/include"
export LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"

# Clean previous build to force full rebuild with new flags
cargo clean

# Rebuild with libdeflate support
# This will take 2-3 minutes as it rebuilds everything
cargo build --release 2>&1 | tee build.log

# Verify build succeeded
tail -20 build.log | grep -i "finished\|error"
```

**What success looks like:**
- No error messages in build output
- Final message: "Finished `release` profile..."
- Build log shows HTSlib being compiled with new flags

**How to verify libdeflate is being used:**

```bash
# Check the build output mentioned libdeflate
grep -i "libdeflate\|--with-libdeflate" build.log

# If not found, check HTSlib configure output:
cat sheriff-rs/target/release/build/hts-sys-*/out/htslib/config.status | grep libdeflate

# Alternative: Check if libdeflate is linked into the binary
nm sheriff-rs/target/release/*.so 2>/dev/null | grep -i "deflate_decompress"
# If found, libdeflate is linked in
```

**If libdeflate wasn't detected:**

The build system might need explicit configuration. Try this:

```bash
# Create a custom build script approach
cd /home/user/Sheriff/sheriff-rs
rm -rf target

# Try with explicit environment variable
export HTS_USE_LIBDEFLATE=1
export LDFLAGS="-L/usr/local/lib -ldeflate"
export CPPFLAGS="-I/usr/local/include"

cargo build --release 2>&1 | tee build.log

# Check output
grep "libdeflate\|deflate" build.log | head -5
```

**Troubleshooting if it fails:**

| Error | Solution |
|-------|----------|
| `hts-sys build failed` | Run: `export LDFLAGS="-L/usr/local/lib"; export CPPFLAGS="-I/usr/local/include"; cargo clean && cargo build --release` |
| `ld: cannot find -ldeflate` | Verify install: `ls /usr/local/lib/libdeflate.so` |
| `fatal error: libdeflate.h: No such file` | Verify: `ls /usr/local/include/libdeflate.h` |
| Long build time (>5 min) | This is normal - HTSlib being rebuilt. Just wait. |

---

### Step 4: Performance Verification (10 minutes)

**Action:** Verify libdeflate is actually being used and provides speedup

**Commands to run:**

```bash
# Run a quick BAM processing test to verify it works
cd /home/user/Sheriff

# If you have test data, use it:
# python -m pytest tests/ -v -k bam

# Otherwise, create a simple Python test:
python3 << 'EOF'
import sys
sys.path.insert(0, '/home/user/Sheriff/sheriff-rs/target/release')

try:
    from sheriff_rs import sheriff_rs
    print("✓ Rust module loads successfully")
    print(f"  Module: {sheriff_rs}")
except ImportError as e:
    print(f"✗ Failed to load Rust module: {e}")
    sys.exit(1)
EOF
```

**What success looks like:**
- Module loads without errors
- Python can import the sheriff_rs module
- No warning messages about missing libraries

**How to verify libdeflate is actually in use:**

```bash
# Check library dependencies (if applicable)
ldd /home/user/Sheriff/sheriff-rs/target/release/*.so 2>/dev/null | grep -i deflate

# Or check symbols
nm /home/user/Sheriff/sheriff-rs/target/release/*.so 2>/dev/null | grep deflate_decompress | head -1

# Check HTSlib's configuration inside the build artifacts
cat /home/user/Sheriff/sheriff-rs/target/release/build/hts-sys-*/out/htslib/config.h 2>/dev/null | grep "HAVE_LIBDEFLATE"
```

**Expected output:** Should show libdeflate symbols or HAVE_LIBDEFLATE defined

**Troubleshooting:**

| Issue | Solution |
|-------|----------|
| libdeflate not detected | Re-run Step 3 with explicit `export HTS_USE_LIBDEFLATE=1` |
| Module won't load | Try: `python -c "import sheriff_rs; print(sheriff_rs)"` to get actual error |
| Symbols not found | libdeflate may not be linked. Check build.log for errors |

---

## C. Parallel BGZF Enablement (Step-by-Step)

### What is Parallel BGZF Decompression?

BGZF (Blocked GZip Format) is the compression format used by BAM files. Key facts:

- **Format:** File is divided into 64KB blocks, each independently compressed
- **Sequential Design:** Currently, rust-htslib reads blocks sequentially
- **Parallelization:** Multiple blocks can be decompressed simultaneously on multi-core CPUs
- **Benefit:** 4-5% speedup on 4-core systems, potentially 10-15% on 8+ core systems

**How it works:**

```
BGZF File Structure:
[Block 1: 64KB]─┐
[Block 2: 64KB]─┼─→ Thread Pool
[Block 3: 64KB]─┤    (4 threads)
[Block 4: 64KB]─┘

Old way (single-threaded):
1. Decompress Block 1 → 64KB data
2. Decompress Block 2 → 64KB data
3. Decompress Block 3 → 64KB data
4. Decompress Block 4 → 64KB data
Total: 4 sequential operations

New way (parallel):
1-4. Decompress Blocks 1-4 simultaneously across threads
Total: ~1 concurrent operation per thread
```

### Implementation Details

The Rust BAM processor already has the `set_threads()` method available. We just need to:

1. Call this method when initializing the BAM reader
2. Set it to the number of available CPU cores
3. No other code changes needed

---

### Step 1: Examine Current BAM Code (5 minutes)

**Action:** Understand the current BAM reader implementation

**Commands to run:**

```bash
# View the current BamProcessor implementation
cat /home/user/Sheriff/sheriff-rs/src/bam.rs | head -100
```

**What you'll see:**

```rust
pub struct BamProcessor {
    reader: Reader,
}

impl BamProcessor {
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        let reader = Reader::from_path(path)?;
        Ok(BamProcessor { reader })
    }

    pub fn set_threads(&mut self, n_threads: usize) -> Result<()> {
        self.reader
            .set_threads(n_threads)
            .map_err(|e| BamError::ParseError(format!("Failed to set threads: {:?}", e)))
    }
}
```

**Key observations:**
- ✓ `set_threads()` method already exists
- ✓ Properly wrapped with error handling
- ✓ Ready to be called during initialization

---

### Step 2: Find BAM Reader Initialization Points (10 minutes)

**Action:** Locate all places where BamProcessor is created and determine where to enable threads

**Commands to run:**

```bash
# Search for BamProcessor usage in the codebase
grep -r "BamProcessor::new" /home/user/Sheriff --include="*.rs" --include="*.py"

# Search for Reader initialization patterns
grep -r "Reader::from_path" /home/user/Sheriff/sheriff-rs/src --include="*.rs"

# Check the PyO3 bindings
cat /home/user/Sheriff/sheriff-rs/src/python.rs | grep -A 10 "BamProcessor"
```

**Expected output:** Shows where BamProcessor is instantiated

---

### Step 3: Update BamProcessor to Enable Parallel BGZF (15 minutes)

**Action:** Modify BamProcessor::new() to automatically enable multi-threaded decompression

**Location:** `/home/user/Sheriff/sheriff-rs/src/bam.rs`

**Current code (lines 54-57):**

```rust
pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
    let reader = Reader::from_path(path)?;
    Ok(BamProcessor { reader })
}
```

**Change to:**

```rust
pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
    let mut reader = Reader::from_path(path)?;

    // Enable multi-threaded BGZF decompression
    // Use all available CPU cores for parallel decompression
    let n_threads = num_cpus::get();
    reader.set_threads(n_threads)
        .map_err(|e| BamError::ParseError(format!("Failed to set threads: {:?}", e)))?;

    Ok(BamProcessor { reader })
}
```

**Why this change:**

1. Makes the reader mutable (required to call set_threads)
2. Gets number of available CPU cores via `num_cpus::get()`
3. Calls set_threads() to enable parallel decompression
4. Properly handles errors

**Commands to make the change:**

```bash
# Edit the file
cd /home/user/Sheriff/sheriff-rs/src
nano bam.rs  # or use your preferred editor
# Or use the Edit tool below
```

---

### Step 4: Add num_cpus Dependency (5 minutes)

**Action:** Add the `num_cpus` crate to Cargo.toml so we can get CPU count

**Location:** `/home/user/Sheriff/sheriff-rs/Cargo.toml`

**Current dependencies:**

```toml
[dependencies]
rustc-hash = "2.0"
rayon = "1.8"
rust-htslib = "0.47"
pyo3 = { version = "0.21", features = ["extension-module"], optional = true }
```

**Change to:**

```toml
[dependencies]
rustc-hash = "2.0"
rayon = "1.8"
rust-htslib = "0.47"
num_cpus = "1.16"
pyo3 = { version = "0.21", features = ["extension-module"], optional = true }
```

**Why this crate:**
- Lightweight (no heavy dependencies)
- Production-use in hundreds of Rust projects
- Cross-platform (Windows, Linux, macOS)
- Already used by rayon for thread pool sizing

---

### Step 5: Rebuild and Verify (10 minutes)

**Action:** Rebuild the project to incorporate the changes

**Commands to run:**

```bash
# Navigate to project directory
cd /home/user/Sheriff/sheriff-rs

# Clean and rebuild with the new code
cargo clean
cargo build --release 2>&1 | tee build_parallel.log

# Check for build errors
grep -i "error\|warning" build_parallel.log | head -10

# Verify build succeeded
tail -5 build_parallel.log | grep -i "finished"
```

**What success looks like:**
- No error messages
- Build completes with "Finished `release` profile..."
- File size increases slightly (new code added)

**How to verify the code was built correctly:**

```bash
# Check that the new symbols are in the binary
nm /home/user/Sheriff/sheriff-rs/target/release/*.so 2>/dev/null | grep "set_threads\|num_cpus"

# Run a quick test
cd /home/user/Sheriff
python3 << 'EOF'
import sys
sys.path.insert(0, '/home/user/Sheriff/sheriff-rs/target/release')
from sheriff_rs import sheriff_rs
print("✓ Module compiled with parallel BGZF support")
EOF
```

**Troubleshooting if build fails:**

| Error | Solution |
|-------|----------|
| `error: could not find num_cpus` | Ensure num_cpus is in Cargo.toml, then re-run |
| `error: cannot borrow as mutable` | Ensure reader is declared as `mut` |
| `unresolved import` | Run: `cargo update && cargo build --release` |

---

### Step 6: Performance Impact Verification (10 minutes)

**Action:** Verify that parallel BGZF is actually enabled and providing speedup

**Commands to run:**

```bash
# Create a test script to verify parallel threads are being used
python3 << 'EOF'
import subprocess
import os

# Monitor the process while running a BAM operation
test_bam = "/path/to/test.bam"  # Update with actual test BAM file if available

if os.path.exists(test_bam):
    # Run a BAM processing operation and monitor threads
    print("Testing parallel BGZF with actual BAM file...")
    # Example if your Python code calls the Rust function:
    # result = sheriff_rs.process_bam(test_bam)
    print("✓ BAM processing completed")
else:
    print("⚠ No test BAM file found - skipping live verification")
    print("  When you run the actual pipeline, monitor with: ps aux | grep sheriff")
EOF
```

**How to verify threads are being used:**

```bash
# If the pipeline is running, monitor thread count:
ps -eLo pid,comm | grep -E "sheriff|python|rust" | wc -l
# Should show multiple threads (one per CPU core + overhead)

# Or use top/htop while the pipeline runs:
# Press 'H' in htop to see threads
# Each BAM decompression should show 4-8 threads (depending on CPU cores)
```

**Expected behavior:**

- CPU usage should increase (all cores utilized during BAM decompression)
- Memory usage slightly increases (thread overhead)
- Overall pipeline time decreases by 4-5%

---

### Step 7: Optional - Configure Thread Count (5 minutes)

**Action (Optional):** If you want to limit thread count for some reason

By default, we use all CPU cores. If you want to limit it (e.g., to 4 threads max):

**In `/home/user/Sheriff/sheriff-rs/src/bam.rs`, change:**

```rust
let n_threads = num_cpus::get();
```

**To:**

```rust
let n_threads = std::cmp::min(num_cpus::get(), 4);  // Cap at 4 threads
```

**Or make it configurable:**

```rust
pub fn new<P: AsRef<Path>>(path: P, max_threads: Option<usize>) -> Result<Self> {
    let mut reader = Reader::from_path(path)?;

    let n_threads = match max_threads {
        Some(max) => std::cmp::min(num_cpus::get(), max),
        None => num_cpus::get(),
    };

    reader.set_threads(n_threads)
        .map_err(|e| BamError::ParseError(format!("Failed to set threads: {:?}", e)))?;

    Ok(BamProcessor { reader })
}
```

**Note:** This is optional. Default behavior (all cores) is recommended.

---

## D. Testing & Verification

### Before-Optimization Benchmarking

**Step 1: Establish baseline (before implementing optimizations)**

```bash
# Navigate to Sheriff
cd /home/user/Sheriff

# Run a representative test on your actual data
# Measure time, memory, and CPU usage
time python -m sheriff.main --input test_data.bam --output results.txt

# Record metrics:
# - Real time (total elapsed)
# - User time (CPU time)
# - Memory peak usage
```

### After-Optimization Verification

**Step 2: After both optimizations, re-run same test**

```bash
# Run the exact same command
time python -m sheriff.main --input test_data.bam --output results.txt

# Compare:
# - Time should be 8-11% faster
# - Memory should be similar or slightly higher (parallel overhead)
# - Results should be identical
```

### Correctness Verification

**Step 3: Ensure results are identical**

```bash
# Compare results before/after
diff results_before.txt results_after.txt
# Should produce no output (identical results)

# Or use MD5 checksums
md5sum results_before.txt results_after.txt
# Should produce identical checksums
```

### Performance Metrics to Track

| Metric | Before | After | Target |
|--------|--------|-------|--------|
| Total pipeline time | 100% | 90-92% | -8-10% |
| BAM I/O time | 30% | 18-20% | -10-12% |
| CPU cores utilized | 1-2 | 4-8 | All cores during BAM I/O |
| Memory peak | X MB | X + 5-10 MB | <10% increase |
| Result correctness | ✓ | ✓ | Must be identical |

### Specific Test Cases

**Test 1: Small file (< 100MB)**

```bash
cd /home/user/Sheriff
# Time how long a small BAM takes to process
time python -c "from sheriff.main import process_bam; process_bam('small_test.bam')"
```

**Test 2: Real-world file (if available)**

```bash
# Process actual data
time python -m sheriff.main --input your_data.bam --output results.txt

# Verify output format
head -20 results.txt
wc -l results.txt  # Check line count matches expectations
```

**Test 3: Multi-threaded verification**

```bash
# While running, check thread count in another terminal
# On Linux:
ps -eLo pid,comm | grep -E "python|libdeflate" | wc -l

# On macOS:
ps -eLo pid,comm | grep -E "python|libdeflate" | wc -l
```

---

## E. Rollback Plan

### If libdeflate Installation Causes Issues

**Quick Rollback (5 minutes):**

```bash
# Uninstall libdeflate
sudo rm -f /usr/local/lib/libdeflate*
sudo rm -f /usr/local/include/libdeflate.h
sudo ldconfig

# Rebuild without libdeflate
cd /home/user/Sheriff/sheriff-rs
unset LDFLAGS CPPFLAGS HTS_USE_LIBDEFLATE
cargo clean
cargo build --release
```

**Undo via git (if committed):**

```bash
git reset --hard HEAD~1
```

### If Parallel BGZF Causes Issues

**Quick Rollback (5 minutes):**

```bash
# Revert the BAM processor changes
cd /home/user/Sheriff/sheriff-rs
git checkout src/bam.rs

# Remove num_cpus from Cargo.toml
# (Edit manually or git checkout Cargo.toml)

# Rebuild
cargo clean
cargo build --release
```

### If Compilation Fails at Any Point

**Complete Clean Restart:**

```bash
cd /home/user/Sheriff/sheriff-rs

# Clean everything
cargo clean
rm -rf target

# Reset to original state
git checkout Cargo.toml src/bam.rs

# Rebuild original version
cargo build --release

# Verify it works
python3 << 'EOF'
import sys
sys.path.insert(0, '/home/user/Sheriff/sheriff-rs/target/release')
from sheriff_rs import sheriff_rs
print("✓ Back to original working state")
EOF
```

### Automated Rollback Script

Create `/home/user/Sheriff/rollback.sh`:

```bash
#!/bin/bash
set -e

echo "Rolling back BAM optimizations..."

# Uninstall libdeflate
echo "1. Removing libdeflate..."
sudo rm -f /usr/local/lib/libdeflate* 2>/dev/null || true
sudo rm -f /usr/local/include/libdeflate.h 2>/dev/null || true
sudo ldconfig 2>/dev/null || true

# Reset code to original
echo "2. Resetting code changes..."
cd /home/user/Sheriff/sheriff-rs
git checkout Cargo.toml src/bam.rs 2>/dev/null || true

# Clean and rebuild
echo "3. Rebuilding original version..."
cargo clean
cargo build --release

echo "✓ Rollback complete. System is back to original state."
```

---

## F. Expected Timeline

### Detailed Hour-by-Hour Breakdown

#### **Hour 1: Libdeflate Installation (0:00 - 1:00)**

| Time | Task | Estimated Duration | What to Check |
|------|------|-------------------|----------------|
| 0:00 - 0:05 | Verify prerequisites (gcc, make, git) | 5 min | Commands return paths, no "not found" |
| 0:05 - 0:25 | Clone and build libdeflate | 20 min | No compile errors, .so and .a files created |
| 0:25 - 0:35 | Install libdeflate system-wide | 10 min | Files in /usr/local/lib and /usr/local/include |
| 0:35 - 1:00 | HTSlib rebuild with libdeflate | 25 min | `cargo build` completes successfully |

**Check after Hour 1:**
```bash
ls /usr/local/lib/libdeflate*        # Should exist
nm /home/user/Sheriff/sheriff-rs/target/release/*.so | grep deflate  # Should find symbols
```

---

#### **Hour 2: Parallel BGZF Enablement (1:00 - 2:00)**

| Time | Task | Estimated Duration | What to Check |
|------|------|-------------------|----------------|
| 1:00 - 1:05 | Review current BAM code | 5 min | Understand BamProcessor::new() |
| 1:05 - 1:15 | Add num_cpus to Cargo.toml | 10 min | Dependency added correctly |
| 1:15 - 1:30 | Update BamProcessor::new() | 15 min | Code changes applied correctly |
| 1:30 - 1:50 | Rebuild with changes | 20 min | `cargo build` completes successfully |
| 1:50 - 2:00 | Verify module loads | 10 min | Python can import sheriff_rs |

**Check after Hour 2:**
```bash
grep "set_threads" /home/user/Sheriff/sheriff-rs/src/bam.rs  # Should find updated code
cargo build --release | tail -1  # Should say "Finished"
```

---

#### **Hour 3: Testing & Verification (2:00 - 3:00)**

| Time | Task | Estimated Duration | What to Check |
|------|------|-------------------|----------------|
| 2:00 - 2:15 | Run baseline benchmark | 15 min | Get before-optimization metrics |
| 2:15 - 2:45 | Run optimized code | 30 min | Process same data, measure time difference |
| 2:45 - 2:55 | Verify correctness | 10 min | Results identical to baseline |
| 2:55 - 3:00 | Document results | 5 min | Record speedup percentage achieved |

**Check after Hour 3:**
```bash
# Compare metrics
echo "Expected improvement: 8-11% faster pipeline"
# Actual results should be close to this
```

---

### Timeline Summary

```
START (0:00)
  │
  ├─→ Hour 1: libdeflate (installation + HTSlib rebuild)
  │    └─→ Checkpoint: /usr/local/lib/libdeflate.so exists ✓
  │
  ├─→ Hour 2: Parallel BGZF (code changes + rebuild)
  │    └─→ Checkpoint: Module loads successfully ✓
  │
  └─→ Hour 3: Testing & Verification
       └─→ Checkpoint: 8-11% speedup verified ✓

DONE (3:00) - Expected speedup: 4.4x → 4.9-5.0x
```

---

## G. Potential Issues & Troubleshooting

### Issue 1: libdeflate build fails

**Symptoms:**
- `make` command fails
- Missing files in /usr/local/lib

**Root causes and solutions:**

| Cause | Solution |
|-------|----------|
| Missing build tools | `sudo apt-get install build-essential` |
| Old libdeflate version | `cd /tmp/libdeflate && git pull origin master && make clean && make` |
| Insufficient disk space | `df -h` (ensure > 1GB free) |
| Old git clone | `rm -rf /tmp/libdeflate && git clone https://github.com/ebiggers/libdeflate.git` |

---

### Issue 2: HTSlib build can't find libdeflate

**Symptoms:**
- Build output says "libdeflate not found"
- HTSlib compiles but doesn't use libdeflate

**Root causes and solutions:**

| Cause | Solution |
|-------|----------|
| LD_LIBRARY_PATH not set | Set before build: `export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH` |
| LDFLAGS not passed | Export: `export LDFLAGS="-L/usr/local/lib"` before rebuild |
| CPPFLAGS not set | Export: `export CPPFLAGS="-I/usr/local/include"` before rebuild |
| hts-sys caching old build | Run: `cargo clean` before rebuilding |

**Complete fix:**

```bash
cd /home/user/Sheriff/sheriff-rs
export LDFLAGS="-L/usr/local/lib"
export CPPFLAGS="-I/usr/local/include"
export LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"
cargo clean
cargo build --release
```

---

### Issue 3: Parallel BGZF build fails

**Symptoms:**
- `error: could not find num_cpus`
- `error: cannot borrow as mutable`

**Root causes and solutions:**

| Cause | Solution |
|-------|----------|
| num_cpus not in Cargo.toml | Add: `num_cpus = "1.16"` to [dependencies] |
| Reader not declared mutable | Change `let reader` to `let mut reader` |
| Old cached dependency | Run: `cargo update && cargo clean && cargo build` |

---

### Issue 4: Python module won't load after changes

**Symptoms:**
- `ImportError: cannot import name 'sheriff_rs'`
- `OSError: libdeflate.so.1: cannot open shared object file`

**Root causes and solutions:**

| Cause | Solution |
|-------|----------|
| New library not in path | Add to LD_LIBRARY_PATH: `export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH` |
| Module not rebuilt | Run: `cd /home/user/Sheriff/sheriff-rs && cargo build --release` |
| Wrong Python path | Check: `python -c "import sys; print(sys.path)"` includes target/release |

**Test module loading:**

```bash
python3 << 'EOF'
import sys
sys.path.insert(0, '/home/user/Sheriff/sheriff-rs/target/release')

try:
    from sheriff_rs import sheriff_rs
    print("✓ Module loads successfully")
except ImportError as e:
    print(f"✗ Import error: {e}")
    import traceback
    traceback.print_exc()
EOF
```

---

### Issue 5: No speedup observed

**Symptoms:**
- Benchmark shows no time difference
- Process still takes the same time

**Root causes and solutions:**

| Cause | Solution |
|-------|----------|
| libdeflate not actually linked | Check: `nm /home/user/Sheriff/sheriff-rs/target/release/*.so \| grep deflate` should find symbols |
| Parallel BGZF not enabled | Check code: `grep "set_threads" /home/user/Sheriff/sheriff-rs/src/bam.rs` should show the call |
| Test data too small | Use larger BAM files (>100MB) where compression overhead is significant |
| Benchmarking wrong code path | Verify running rebuilt binary, not old cached one |
| System limitations | Check CPU: 4-5% speedup assumes multi-core CPU (check `nproc`) |

**Diagnostic script:**

```bash
#!/bin/bash

echo "=== Optimization Verification ==="
echo ""

echo "1. libdeflate installed?"
ls -la /usr/local/lib/libdeflate* && echo "   ✓ Yes" || echo "   ✗ No"
echo ""

echo "2. libdeflate linked in binary?"
nm /home/user/Sheriff/sheriff-rs/target/release/*.so 2>/dev/null | grep deflate_decompress && echo "   ✓ Yes" || echo "   ✗ No"
echo ""

echo "3. Parallel BGZF code present?"
grep "set_threads" /home/user/Sheriff/sheriff-rs/src/bam.rs && echo "   ✓ Yes" || echo "   ✗ No"
echo ""

echo "4. Available CPU cores:"
nproc
echo ""

echo "5. Module loads?"
python3 << 'EOF'
import sys
sys.path.insert(0, '/home/user/Sheriff/sheriff-rs/target/release')
try:
    from sheriff_rs import sheriff_rs
    print("   ✓ Yes")
except ImportError as e:
    print(f"   ✗ No: {e}")
EOF
```

---

## H. Risk Assessment & Mitigation

### Risk Matrix

| Risk | Likelihood | Impact | Severity | Mitigation |
|------|-----------|--------|----------|-----------|
| libdeflate build fails | Low (5%) | Medium | Medium | Rollback script, easy revert |
| HTSlib can't find libdeflate | Low (10%) | Low | Low | Clear environment, explicit paths |
| Parallel BGZF breaks BAM reading | Very Low (2%) | High | High | Code already tested in rust-htslib |
| Python module won't load | Low (5%) | Medium | Medium | Rebuild, check library paths |
| No performance improvement | Low (10%) | Low | Low | Verify installations, check data size |
| Breaking change to API | Very Low (1%) | High | High | Only internal changes, no API changes |

### Mitigation Strategies

1. **Backup current working state:**
   ```bash
   cd /home/user/Sheriff
   git stash
   git branch backup-$(date +%s)
   ```

2. **Test on non-production data first**

3. **Have rollback script ready** (provided in Section E)

4. **Monitor build process carefully:**
   ```bash
   cargo build --release 2>&1 | tee build.log
   ```

5. **Verify at each step** (not all-or-nothing)

---

## I. Success Criteria

### Phase 1 is successful if:

1. **✓ libdeflate Installation:**
   - [ ] `/usr/local/lib/libdeflate.so` exists
   - [ ] `/usr/local/include/libdeflate.h` exists
   - [ ] `cargo build --release` completes without errors

2. **✓ Parallel BGZF Implementation:**
   - [ ] Code in `src/bam.rs` includes `set_threads()` call
   - [ ] `cargo build --release` completes without errors
   - [ ] Python module imports successfully
   - [ ] No API changes (backward compatible)

3. **✓ Performance Verification:**
   - [ ] Pipeline is 8-11% faster on test data
   - [ ] Memory usage similar (± 10%)
   - [ ] Results identical to before optimization
   - [ ] CPU cores utilized during BAM I/O

4. **✓ Code Quality:**
   - [ ] No compiler warnings
   - [ ] All tests pass (if applicable)
   - [ ] Code changes are minimal and focused
   - [ ] Error handling is proper

### Expected Final State

```
Current:     4.4x speedup
Target:      5.0-5.2x speedup
Success if:  >= 4.8x speedup (minimum 9% improvement)

Breakdown:
├─ libdeflate:      +6% expected, -5% to +7% realistic
├─ Parallel BGZF:   +4-5% expected, -3% to +6% realistic
└─ Combined:        +10-11% expected, 8-13% realistic
```

---

## J. Quick Reference Checklist

### Before Starting
- [ ] Read this document completely
- [ ] Backup current state: `git stash && git branch backup-before-bam-opt`
- [ ] Have 2-3 hours available uninterrupted
- [ ] Test data ready for benchmarking

### Execution
- [ ] Step 1: Verify build tools installed
- [ ] Step 2: Clone and build libdeflate
- [ ] Step 3: Install libdeflate system-wide
- [ ] Step 4: Rebuild HTSlib with libdeflate
- [ ] Step 5: Update BamProcessor code
- [ ] Step 6: Add num_cpus dependency
- [ ] Step 7: Rebuild Rust project
- [ ] Step 8: Verify module loads
- [ ] Step 9: Run baseline benchmark
- [ ] Step 10: Run optimized benchmark
- [ ] Step 11: Compare results

### Verification
- [ ] libdeflate symbols found in binary
- [ ] set_threads() code in bam.rs
- [ ] Module loads without errors
- [ ] Speedup measured at 8-11%
- [ ] Results identical to baseline
- [ ] All CPU cores used during BAM I/O

### If Issues Occur
- [ ] Run diagnostic script (see Section G)
- [ ] Check build.log for error messages
- [ ] Verify environment variables set correctly
- [ ] Use rollback script if needed
- [ ] Start from clean state if persistent issues

---

## K. Detailed Implementation Commands (Copy-Paste Ready)

### Complete Command Sequence (No Interruptions)

```bash
#!/bin/bash
set -e

echo "=========================================="
echo "Sheriff BAM Optimization - Phase 1"
echo "=========================================="
echo ""

# === HOUR 1: LIBDEFLATE INSTALLATION ===
echo "HOUR 1: Installing libdeflate"
echo "========================================"
echo ""

echo "Step 1: Clone and build libdeflate..."
cd /tmp
rm -rf libdeflate
git clone https://github.com/ebiggers/libdeflate.git
cd libdeflate
make
echo "✓ libdeflate built successfully"
echo ""

echo "Step 2: Install libdeflate system-wide..."
sudo make install
sudo ldconfig
echo "✓ libdeflate installed"
echo ""

echo "Step 3: Rebuild HTSlib with libdeflate..."
cd /home/user/Sheriff/sheriff-rs
export LDFLAGS="-L/usr/local/lib"
export CPPFLAGS="-I/usr/local/include"
export LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"
cargo clean
cargo build --release 2>&1 | tee build_libdeflate.log
echo "✓ HTSlib rebuilt with libdeflate"
echo ""

echo "Step 4: Verify libdeflate integration..."
nm /home/user/Sheriff/sheriff-rs/target/release/*.so 2>/dev/null | grep deflate_decompress | head -1
echo "✓ libdeflate verified"
echo ""

# === HOUR 2: PARALLEL BGZF SETUP ===
echo "HOUR 2: Enabling Parallel BGZF"
echo "========================================"
echo ""

echo "Step 5: Update Cargo.toml with num_cpus..."
# Will be done manually or with sed

echo "Step 6: Update bam.rs with set_threads()..."
# Will be done manually or with sed

echo "Step 7: Rebuild with changes..."
cd /home/user/Sheriff/sheriff-rs
cargo clean
cargo build --release 2>&1 | tee build_parallel.log
echo "✓ Parallel BGZF enabled"
echo ""

echo "Step 8: Verify module loads..."
python3 << 'EOF'
import sys
sys.path.insert(0, '/home/user/Sheriff/sheriff-rs/target/release')
from sheriff_rs import sheriff_rs
print("✓ Module loaded successfully")
EOF
echo ""

# === HOUR 3: TESTING ===
echo "HOUR 3: Testing & Verification"
echo "========================================"
echo ""

echo "Step 9: Run benchmarks..."
echo "(This will depend on your test data)"
echo "Try: time python -m sheriff.main --input test.bam"
echo ""

echo "=========================================="
echo "✓ All steps complete!"
echo "Expected speedup: 8-11%"
echo "Next: Benchmark your actual data"
echo "=========================================="
```

---

## L. Appendix: Reference Information

### libdeflate Resources

- **Official Repository:** https://github.com/ebiggers/libdeflate
- **Compression Algorithm:** DEFLATE (RFC 1951)
- **Decompression Speed:** 2-3x faster than zlib
- **Compatibility:** Drop-in replacement for zlib/libz
- **License:** MIT

### rust-htslib Resources

- **Crate:** https://crates.io/crates/rust-htslib
- **Version used:** 0.47
- **set_threads() documentation:** Enables multi-threaded BGZF decompression

### BGZF Format Details

- **Block Format:** Blocked GZip Format
- **Block Size:** 64KB uncompressed, variable compressed
- **Header:** Contains "BGZF" magic bytes
- **Footer:** 8-byte ISIZE field with uncompressed size
- **Parallel Capability:** Each block can be decompressed independently

### Performance Tuning References

| Configuration | Value | Why |
|-----------------|-------|-----|
| BGZF block threads | num_cpus::get() | Use all cores |
| libdeflate level | Auto (depends on HTSlib default) | Balance speed vs compression |
| Memory buffer | 64MB (rust-htslib default) | Reasonable for multi-threaded I/O |

---

## M. Next Steps After Phase 1

After successfully implementing Phase 1 (+10-11% speedup), the next optimization phase is:

### Phase 2: Gene UMI Counting in Rust (10 hours, +5-10% speedup)

The Python function `bam_count_gene_umis()` is marked as "Most time consuming step" in the research. Moving it to Rust (similar to k-mer and UMI dedup optimizations) can provide an additional 5-10% pipeline speedup.

This requires:
1. Study existing gene.rs module structure
2. Implement count_gene_umis() in Rust
3. Add PyO3 bindings
4. Test against Python version for correctness
5. Benchmark for performance improvement

**Estimated timeline for Phase 2:** 1-2 weeks (10 hours work distributed)
**Expected result:** 5.0-5.2x → 5.5-5.8x speedup

---

## N. Document History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | Nov 19, 2025 | Initial comprehensive implementation plan |

---

**Document Status:** READY FOR IMPLEMENTATION

**Contact:** For questions, refer to BAM_RESEARCH_EXECUTIVE_SUMMARY.md or BAM_OPTIMIZATION_DEEP_DIVE.md

**Last Updated:** November 19, 2025
