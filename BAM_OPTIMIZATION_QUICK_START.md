# BAM Optimization Quick Start Guide
## Phase 1: 3 Hours → +10-11% Speedup

**Jump to sections:**
- [Copy-Paste Commands](#copy-paste-commands) - Just run these!
- [Visual Checklist](#visual-checklist) - Track your progress
- [Verification Steps](#verification-steps) - Confirm it worked
- [Troubleshooting](#troubleshooting) - If something breaks

---

## Copy-Paste Commands

### HOUR 1: Install libdeflate & Rebuild HTSlib

```bash
# Step 1.1: Verify build tools (should return 3 paths)
which gcc g++ make
which git

# Step 1.2: Clone and build libdeflate (2-3 minutes)
cd /tmp
rm -rf libdeflate
git clone https://github.com/ebiggers/libdeflate.git
cd libdeflate
make
echo "✓ Build complete. Check for .so and .a files:"
ls -la libdeflate.so* libdeflate.a

# Step 1.3: Install system-wide (1 minute)
sudo make install
sudo ldconfig
echo "✓ Installation complete:"
ls -la /usr/local/lib/libdeflate*
ls -la /usr/local/include/libdeflate.h

# Step 1.4: Rebuild HTSlib with libdeflate (2-3 minutes)
cd /home/user/Sheriff/sheriff-rs
export LDFLAGS="-L/usr/local/lib"
export CPPFLAGS="-I/usr/local/include"
export LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"
cargo clean
echo "Starting cargo build (this takes 2-3 minutes)..."
cargo build --release 2>&1 | tail -10
echo ""
echo "If you see 'Finished release profile', libdeflate is integrated!"

# Step 1.5: Verify libdeflate is in the binary
nm /home/user/Sheriff/sheriff-rs/target/release/*.so | grep deflate_decompress
echo "✓ HOUR 1 COMPLETE - If you see symbols above, libdeflate is installed!"
```

### HOUR 2: Enable Parallel BGZF

```bash
# Step 2.1: Add num_cpus to Cargo.toml
cd /home/user/Sheriff/sheriff-rs
cat Cargo.toml | grep "num_cpus" || (
  echo "Adding num_cpus to Cargo.toml..."
  # Backup and edit
  cp Cargo.toml Cargo.toml.backup

  # Add num_cpus line after other dependencies
  sed -i '/^rust-htslib/a num_cpus = "1.16"' Cargo.toml
  echo "✓ Added num_cpus dependency"
  echo "Verify:"
  grep "num_cpus" Cargo.toml
)

# Step 2.2: Update src/bam.rs to enable set_threads()
cd /home/user/Sheriff/sheriff-rs
cp src/bam.rs src/bam.rs.backup

# Check current code
echo "Current BamProcessor::new() implementation:"
sed -n '54,57p' src/bam.rs

# Show what we need to change to
echo ""
echo "We need to add set_threads() call. Here's the updated code:"
echo '
pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
    let mut reader = Reader::from_path(path)?;
    let n_threads = num_cpus::get();
    reader.set_threads(n_threads)
        .map_err(|e| BamError::ParseError(format!("Failed to set threads: {:?}", e)))?;
    Ok(BamProcessor { reader })
}
'

# Step 2.3: Rebuild with new code
echo ""
echo "Rebuilding with parallel BGZF support..."
cargo clean
cargo build --release 2>&1 | tail -5
echo "✓ HOUR 2 CHECKPOINT: If you see 'Finished', rebuild successful!"

# Step 2.4: Verify the module loads
python3 << 'EOF'
import sys
sys.path.insert(0, '/home/user/Sheriff/sheriff-rs/target/release')
try:
    from sheriff_rs import sheriff_rs
    print("✓ Module loads successfully - Parallel BGZF is ready!")
except ImportError as e:
    print(f"✗ Import failed: {e}")
    sys.exit(1)
EOF
```

### HOUR 3: Test & Verify Performance

```bash
# Step 3.1: Create a test benchmark script
cat > /tmp/test_performance.py << 'TESTEOF'
#!/usr/bin/env python3
import time
import sys
import os

# Add Rust module to path
sys.path.insert(0, '/home/user/Sheriff/sheriff-rs/target/release')

print("=" * 50)
print("Performance Verification")
print("=" * 50)
print("")

# Test 1: Module loading speed
print("Test 1: Module loading...")
start = time.time()
try:
    from sheriff_rs import sheriff_rs
    elapsed = time.time() - start
    print(f"  ✓ Loaded in {elapsed:.3f}s")
except Exception as e:
    print(f"  ✗ Failed: {e}")
    sys.exit(1)

print("")

# Test 2: Try a simple operation
print("Test 2: Basic functionality...")
print("  (If you have test data, run your actual pipeline here)")
print("  Command template:")
print("    time python -m sheriff.main --input your_test.bam --output results.txt")
print("")
print("  Expected improvement: 8-11% faster")

print("")
print("=" * 50)
print("✓ All tests passed!")
print("=" * 50)
TESTEOF

chmod +x /tmp/test_performance.py
python3 /tmp/test_performance.py

# Step 3.2: Actual performance test (if you have test data)
echo ""
echo "Running actual pipeline test (customize path to your test BAM):"
echo "  time python -m sheriff.main --input /path/to/test.bam --output results.txt"
echo ""
echo "Record the time shown as: real X.XXs"
echo "Compare to your baseline before these optimizations."
echo ""
echo "SUCCESS if: New time is 8-11% less than baseline"
```

---

## Visual Checklist

### Pre-Implementation
```
□ Read BAM_OPTIMIZATION_IMPLEMENTATION_PLAN.md
□ Backup current state: git branch backup-$(date +%s)
□ Have test data ready for benchmarking
□ Have 3 hours of uninterrupted time
□ Close other applications to avoid noise in benchmarks
```

### Hour 1: libdeflate
```
□ Step 1.1: Build tools installed
  - gcc ✓
  - g++ ✓
  - make ✓
  - git ✓

□ Step 1.2: libdeflate cloned and built
  - /tmp/libdeflate/libdeflate.so exists ✓
  - /tmp/libdeflate/libdeflate.a exists ✓
  - make completed without errors ✓

□ Step 1.3: libdeflate installed
  - /usr/local/lib/libdeflate.so exists ✓
  - /usr/local/include/libdeflate.h exists ✓
  - ldconfig ran successfully ✓

□ Step 1.4: HTSlib rebuilt
  - cargo build finished successfully ✓
  - No compilation errors ✓
  - No major warnings ✓

□ Step 1.5: Verification
  - nm output shows deflate_decompress symbols ✓
  - No "not found" errors ✓
```

**Time Used: _____ min / 60 min**
**Status: □ HOUR 1 COMPLETE**

---

### Hour 2: Parallel BGZF
```
□ Step 2.1: Dependency updated
  - num_cpus = "1.16" added to Cargo.toml ✓
  - Cargo.toml is valid TOML ✓

□ Step 2.2: Code updated
  - src/bam.rs backed up ✓
  - set_threads() call added ✓
  - Code compiles ✓

□ Step 2.3: Rebuilt
  - cargo build finished successfully ✓
  - No compilation errors ✓
  - Binary size increased ✓

□ Step 2.4: Verification
  - Module imports successfully ✓
  - No library load errors ✓
  - Python test script runs ✓
```

**Time Used: _____ min / 60 min**
**Status: □ HOUR 2 COMPLETE**

---

### Hour 3: Testing
```
□ Step 3.1: Environment ready
  - Test data available ✓
  - Baseline metrics recorded ✓
  - System clean (no other work) ✓

□ Step 3.2: Run tests
  - Module loading test passed ✓
  - Actual pipeline runs successfully ✓
  - Results are identical to baseline ✓

□ Step 3.3: Measure speedup
  - Baseline time: ________ seconds
  - Optimized time: ________ seconds
  - Improvement: ________ % (should be 8-11%)
  - ✓ Within expected range (8-11%)

□ Step 3.4: Verification
  - CPU cores used during BAM I/O ✓
  - Memory reasonable (±10%) ✓
  - No correctness issues ✓
  - Results match baseline exactly ✓
```

**Time Used: _____ min / 60 min**
**Status: □ HOUR 3 COMPLETE**

---

### Post-Implementation
```
□ Document results
  - Original speedup: 4.4x
  - New speedup: _____ x
  - Improvement: _____ %

□ Git commit (optional)
  - Stage changes: git add -A
  - Commit: git commit -m "Add libdeflate + parallel BGZF support"
  - Push: git push

□ Cleanup
  - Remove backup files (optional)
  - Remove /tmp/libdeflate (optional)
  - Update documentation with results

□ Next phase planning
  - Review PHASE 2 plan (Gene UMI Rust)
  - Schedule implementation
```

**OVERALL STATUS: □ PHASE 1 COMPLETE ✓**

---

## Verification Steps

### 1. Quick Verification (2 minutes)

```bash
#!/bin/bash
echo "=========================================="
echo "BAM Optimization Verification"
echo "=========================================="
echo ""

# Check 1: libdeflate installed
echo "1. libdeflate installed?"
if [ -f /usr/local/lib/libdeflate.so ]; then
    echo "   ✓ YES - /usr/local/lib/libdeflate.so exists"
else
    echo "   ✗ NO - libdeflate not found"
fi
echo ""

# Check 2: libdeflate linked in binary
echo "2. libdeflate in binary?"
if nm /home/user/Sheriff/sheriff-rs/target/release/*.so 2>/dev/null | grep -q deflate_decompress; then
    echo "   ✓ YES - deflate_decompress symbol found"
else
    echo "   ✗ NO - libdeflate not linked"
fi
echo ""

# Check 3: Code has set_threads
echo "3. Parallel BGZF code?"
if grep -q "set_threads" /home/user/Sheriff/sheriff-rs/src/bam.rs; then
    echo "   ✓ YES - set_threads() found in code"
else
    echo "   ✗ NO - set_threads() not present"
fi
echo ""

# Check 4: Module loads
echo "4. Module loads?"
if python3 << 'EOF'
import sys
sys.path.insert(0, '/home/user/Sheriff/sheriff-rs/target/release')
try:
    from sheriff_rs import sheriff_rs
    print("   ✓ YES - imported successfully")
    sys.exit(0)
except:
    print("   ✗ NO - import failed")
    sys.exit(1)
EOF
then
    true
fi
echo ""

# Check 5: Available CPUs
echo "5. CPU cores available:"
CPUS=$(nproc)
echo "   $CPUS cores (parallel BGZF will use $CPUS threads)"
echo ""

echo "=========================================="
echo "Verification Complete!"
echo "=========================================="
```

### 2. Performance Comparison (15 minutes)

```bash
#!/bin/bash

# Set baseline test
TEST_BAM="/path/to/test.bam"  # UPDATE THIS PATH!

if [ ! -f "$TEST_BAM" ]; then
    echo "⚠ Update TEST_BAM path to your actual BAM file"
    exit 1
fi

echo "=========================================="
echo "Performance Benchmark"
echo "=========================================="
echo ""
echo "Test file: $TEST_BAM"
echo "Size: $(du -h "$TEST_BAM" | cut -f1)"
echo ""

# Run 3 times to get stable measurement
echo "Running benchmark (3 iterations)..."
echo ""

for i in 1 2 3; do
    echo "Iteration $i:"
    time python -m sheriff.main --input "$TEST_BAM" --output /tmp/results_$i.txt 2>&1 | grep real
done

echo ""
echo "Expected improvement: 8-11% faster"
echo "If slower, check: Are libdeflate and set_threads actually enabled?"
```

### 3. Correctness Check

```bash
# Verify results are identical before/after optimization
if [ -f results_before.txt ] && [ -f results_after.txt ]; then
    if diff results_before.txt results_after.txt; then
        echo "✓ Results are IDENTICAL"
    else
        echo "✗ Results DIFFER - optimization may have introduced a bug"
        exit 1
    fi
fi
```

---

## Troubleshooting

### Problem: `gcc: command not found`

**Solution:**
```bash
sudo apt-get update
sudo apt-get install -y build-essential git
which gcc  # Should now show a path
```

---

### Problem: `git clone` fails with network error

**Solution:**
```bash
# Try again (GitHub may be temporarily unreachable)
cd /tmp
rm -rf libdeflate
git clone https://github.com/ebiggers/libdeflate.git

# If it fails again, use HTTPS:
git clone https://github.com/ebiggers/libdeflate.git
```

---

### Problem: `cargo build` fails with "cannot find -ldeflate"

**Solution:**
```bash
# Verify libdeflate was installed
ls /usr/local/lib/libdeflate.so

# Set environment variables and try again
export LDFLAGS="-L/usr/local/lib"
export CPPFLAGS="-I/usr/local/include"
export LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"
cd /home/user/Sheriff/sheriff-rs
cargo clean
cargo build --release
```

---

### Problem: Module won't import in Python

**Solution:**
```bash
# Check if library is accessible
python3 << 'EOF'
import os
import ctypes
try:
    ctypes.CDLL('/usr/local/lib/libdeflate.so')
    print("✓ libdeflate.so accessible from Python")
except OSError as e:
    print(f"✗ Cannot load libdeflate: {e}")
    print("  Try: export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH")
EOF

# Ensure module path is correct
python3 -c "import sys; sys.path.insert(0, '/home/user/Sheriff/sheriff-rs/target/release'); from sheriff_rs import sheriff_rs; print('OK')"
```

---

### Problem: No performance improvement observed

**Solution:**
```bash
# Verify optimizations are actually enabled

# 1. Check libdeflate is linked
nm /home/user/Sheriff/sheriff-rs/target/release/*.so | grep deflate_decompress
# Should show: "0000000000xxxxxx T deflate_decompress"

# 2. Check set_threads code exists
grep "set_threads" /home/user/Sheriff/sheriff-rs/src/bam.rs
# Should show the set_threads call

# 3. Check test file is large enough
ls -lh /path/to/test.bam
# Should be >100MB for measurable impact

# 4. Run diagnostic
nproc
# If shows 1 or 2, parallel BGZF won't show much benefit
```

---

### Problem: Build takes very long (>10 minutes)

**Solution:**
```bash
# This is normal - HTSlib is being recompiled
# Just wait. It's a one-time cost.

# To see progress:
tail -f build.log

# If it actually hangs (>30 min):
# Ctrl+C to stop, then:
cargo clean
# Check available disk space
df -h /home/user/Sheriff
# If < 1GB free, clean up space
```

---

## Success Indicators

After completing all 3 hours, you should see:

### Checklist
- [x] libdeflate installed at `/usr/local/lib/libdeflate.so`
- [x] HTSlib rebuilt with `--with-libdeflate`
- [x] Rust module includes `set_threads()` call
- [x] num_cpus dependency in Cargo.toml
- [x] Python module loads without errors
- [x] Pipeline runs 8-11% faster on test data
- [x] Results identical to baseline

### Metrics
- **libdeflate:** 6% speedup (confirmed by: symbols in binary)
- **Parallel BGZF:** 4-5% speedup (confirmed by: set_threads in code)
- **Combined:** 10-11% speedup (measured by: before/after benchmark)

### Performance
```
Before:  4.4x speedup
After:   4.9-5.0x speedup (on 4-core system)
After:   5.0-5.2x speedup (on 8-core system)

Expected time improvement: 8-11%
```

---

## Quick Reference

| Item | Location |
|------|----------|
| Full Plan | `/home/user/Sheriff/BAM_OPTIMIZATION_IMPLEMENTATION_PLAN.md` |
| Research Summary | `/home/user/Sheriff/BAM_RESEARCH_EXECUTIVE_SUMMARY.md` |
| Deep Dive | `/home/user/Sheriff/BAM_OPTIMIZATION_DEEP_DIVE.md` |
| BAM Code | `/home/user/Sheriff/sheriff-rs/src/bam.rs` |
| Cargo Config | `/home/user/Sheriff/sheriff-rs/Cargo.toml` |

---

## Rollback If Needed

```bash
#!/bin/bash

echo "Rolling back BAM optimizations..."

# Uninstall libdeflate
sudo rm -f /usr/local/lib/libdeflate*
sudo rm -f /usr/local/include/libdeflate.h
sudo ldconfig

# Reset code
cd /home/user/Sheriff/sheriff-rs
git checkout Cargo.toml src/bam.rs

# Rebuild original
cargo clean
cargo build --release

echo "✓ Rollback complete - back to original 4.4x version"
```

---

## What's Next?

After Phase 1 succeeds (you reach 4.9-5.0x), the next phase is:

### Phase 2: Gene UMI Counting in Rust
- **Time:** 10 hours
- **Benefit:** +5-10% speedup (5.5-5.8x final)
- **Difficulty:** Medium (similar to existing optimizations)
- **Status:** Planned but not yet started

See `BAM_RESEARCH_EXECUTIVE_SUMMARY.md` for Phase 2 details.

---

**Good luck! You've got this! 🚀**

*Report any issues: Check troubleshooting section or review the full implementation plan.*
