# BAM Optimization - Visual Reference Guide
## Quick Visual Overview of Phase 1 Implementation

---

## 1. Overall Architecture

### Current State (4.4x speedup)
```
Sheriff Pipeline
├─ K-mer Analysis (Rust)         [20% time]
├─ UMI Deduplication (Rust)      [15% time]
├─ BAM Processing (HTSlib)       [30% time] ← BOTTLENECK!
│  ├─ BGZF Decompression        [60% of BAM I/O]
│  ├─ BAM Parsing               [20% of BAM I/O]
│  └─ Tag Extraction            [20% of BAM I/O]
├─ Gene UMI Counting (Python)    [10% time]
└─ Other Processing              [25% time]

Total Pipeline: 100% → 4.4x baseline
```

### After Phase 1 (4.9-5.0x speedup)
```
Sheriff Pipeline
├─ K-mer Analysis (Rust)         [20% time]
├─ UMI Deduplication (Rust)      [15% time]
├─ BAM Processing (HTSlib)       [20% time] ← OPTIMIZED!
│  ├─ BGZF Decompression        [36% of BAM I/O] ← libdeflate
│  ├─ BGZF Parallelism          [multi-threaded] ← parallel BGZF
│  ├─ BAM Parsing               [20% of BAM I/O]
│  └─ Tag Extraction            [20% of BAM I/O]
├─ Gene UMI Counting (Python)    [10% time]
└─ Other Processing              [25% time]

Total Pipeline: 90% → 4.9-5.0x speedup
Improvement: +11% faster
```

---

## 2. Three-Hour Implementation Timeline

```
START (Hour 0:00)
│
├─────────────────────────────────────────────────────────┐
│                                                          │
│  HOUR 1 (0:00 - 1:00): Install libdeflate             │
│  ════════════════════════════════════════════          │
│                                                          │
│  ┌─ 0:00-0:05: Verify build tools                      │
│  │  └─ Check: gcc, g++, make, git installed            │
│  │                                                      │
│  ├─ 0:05-0:25: Clone & build libdeflate (20 min)       │
│  │  └─ cd /tmp && git clone libdeflate && make          │
│  │     Result: /tmp/libdeflate/libdeflate.so           │
│  │                                                      │
│  ├─ 0:25-0:35: Install system-wide (10 min)            │
│  │  └─ sudo make install && sudo ldconfig              │
│  │     Result: /usr/local/lib/libdeflate.so            │
│  │                                                      │
│  └─ 0:35-1:00: Rebuild HTSlib (25 min)                 │
│     └─ cargo clean && cargo build --release            │
│        Result: Binary linked with libdeflate            │
│                                                          │
│  CHECKPOINT ✓:                                          │
│  - /usr/local/lib/libdeflate.so exists                 │
│  - Binary has deflate_decompress symbols              │
│                                                          │
└─────────────────────────────────────────────────────────┘
          ↓
┌─────────────────────────────────────────────────────────┐
│                                                          │
│  HOUR 2 (1:00 - 2:00): Enable parallel BGZF           │
│  ═══════════════════════════════════════════           │
│                                                          │
│  ┌─ 1:00-1:10: Add num_cpus dependency                  │
│  │  └─ Edit Cargo.toml: add num_cpus = "1.16"          │
│  │     File: /home/user/Sheriff/sheriff-rs/Cargo.toml  │
│  │                                                      │
│  ├─ 1:10-1:25: Update bam.rs                           │
│  │  └─ Edit: src/bam.rs BamProcessor::new()            │
│  │     Add: set_threads(num_cpus::get()) call          │
│  │     File: /home/user/Sheriff/sheriff-rs/src/bam.rs │
│  │                                                      │
│  ├─ 1:25-1:50: Rebuild (25 min)                        │
│  │  └─ cargo clean && cargo build --release            │
│  │     Result: Binary with parallel BGZF               │
│  │                                                      │
│  └─ 1:50-2:00: Verify module (10 min)                  │
│     └─ python3: from sheriff_rs import sheriff_rs       │
│        Result: Module loads successfully                │
│                                                          │
│  CHECKPOINT ✓:                                          │
│  - Code has set_threads() call                         │
│  - Module imports without errors                       │
│                                                          │
└─────────────────────────────────────────────────────────┘
          ↓
┌─────────────────────────────────────────────────────────┐
│                                                          │
│  HOUR 3 (2:00 - 3:00): Test & Verify                  │
│  ══════════════════════════════════════════            │
│                                                          │
│  ┌─ 2:00-2:20: Run baseline benchmark                  │
│  │  └─ time python -m sheriff ... < baseline.bam       │
│  │     Record: Time T1                                 │
│  │                                                      │
│  ├─ 2:20-2:45: Run optimized version                   │
│  │  └─ time python -m sheriff ... < baseline.bam       │
│  │     Record: Time T2                                 │
│  │                                                      │
│  ├─ 2:45-2:55: Calculate improvement                   │
│  │  └─ Improvement = (T1 - T2) / T1 * 100%             │
│  │     Expected: 8-11%                                 │
│  │                                                      │
│  └─ 2:55-3:00: Verify correctness                      │
│     └─ diff results_before.txt results_after.txt       │
│        Expected: No output (identical)                 │
│                                                          │
│  CHECKPOINT ✓:                                          │
│  - Speedup measured at 8-11%                          │
│  - Results are identical to baseline                   │
│                                                          │
└─────────────────────────────────────────────────────────┘
          ↓
       DONE! (3:00)

   Current: 4.4x speedup
   After:   4.9-5.0x speedup
   Gain:    +11% improvement
```

---

## 3. Code Change Summary

### Change 1: Cargo.toml (Add dependency)

**Location:** `/home/user/Sheriff/sheriff-rs/Cargo.toml`

```diff
[dependencies]
rustc-hash = "2.0"
rayon = "1.8"
rust-htslib = "0.47"
+ num_cpus = "1.16"
pyo3 = { version = "0.21", features = ["extension-module"], optional = true }
```

**Impact:** 1 line addition

---

### Change 2: src/bam.rs (Enable parallel decompression)

**Location:** `/home/user/Sheriff/sheriff-rs/src/bam.rs`

**Before:**
```rust
pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
    let reader = Reader::from_path(path)?;
    Ok(BamProcessor { reader })
}
```

**After:**
```rust
pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
    let mut reader = Reader::from_path(path)?;

    // Enable multi-threaded BGZF decompression
    let n_threads = num_cpus::get();
    reader.set_threads(n_threads)
        .map_err(|e| BamError::ParseError(format!("Failed to set threads: {:?}", e)))?;

    Ok(BamProcessor { reader })
}
```

**Impact:** 7 lines changed/added

---

### Change 3: Build environment (Set during cargo build)

**When running:** `cargo build --release`

```bash
export LDFLAGS="-L/usr/local/lib"
export CPPFLAGS="-I/usr/local/include"
export LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"
cargo clean
cargo build --release
```

**Impact:** 4 environment variables set

---

## 4. File Modification Map

```
/home/user/Sheriff/
├── sheriff-rs/
│   ├── Cargo.toml                    ← MODIFY (add num_cpus)
│   ├── src/
│   │   ├── bam.rs                    ← MODIFY (add set_threads)
│   │   └── [other files]             [NO CHANGES]
│   └── target/
│       └── release/
│           └── *.so                  ← REBUILD (will link libdeflate)
│
└── /usr/local/                       ← NEW FILES (libdeflate)
    ├── lib/
    │   ├── libdeflate.so             ← INSTALLED
    │   └── libdeflate.a              ← INSTALLED
    └── include/
        └── libdeflate.h              ← INSTALLED
```

---

## 5. Performance Impact Diagram

### BAM Processing Timeline (Before)
```
Time: 0.0s ─────────────────────────────────────────── 1.0s
            ├─ Decompress Block 1  [0.2s]
            ├─ Decompress Block 2  [0.2s]
            ├─ Decompress Block 3  [0.2s]
            ├─ Decompress Block 4  [0.2s]
            └─ Parse & Extract     [0.2s]

Single-threaded decompression
Total: 1.0s (baseline)
```

### BAM Processing Timeline (After - With libdeflate)
```
Time: 0.0s ─────────────────────── 0.85s
            ├─ Decompress Blocks [0.17s] ← libdeflate faster
            └─ Parse & Extract   [0.2s]

Total: 0.85s (14% faster)
```

### BAM Processing Timeline (After - With parallel BGZF)
```
Time: 0.0s ──────────────────── 0.74s
            ├─ Thread 1: Blocks 1,5,9...   ┐
            ├─ Thread 2: Blocks 2,6,10...  │ Parallel
            ├─ Thread 3: Blocks 3,7,11...  │ decompression
            ├─ Thread 4: Blocks 4,8,12...  ┤ [0.17s → 0.05s each]
            └─ Parse & Extract             │
                       [0.2s]

Total: 0.74s (26% faster than original)
```

---

## 6. Dependency Hierarchy

```
Sheriff-rs
├─ Cargo.toml
│  ├─ rust-htslib 0.47
│  │  └─ (will be rebuilt to use libdeflate)
│  │
│  ├─ num_cpus 1.16 ← NEW DEPENDENCY
│  │  └─ Standard library only (zero heavy deps)
│  │
│  └─ rayon 1.8
│     └─ (parallel iterators, not used for BAM)
│
System
├─ /usr/local/lib/libdeflate.so ← NEW SYSTEM LIBRARY
│  └─ Compiled from source
│
└─ /home/user/Sheriff/sheriff-rs/target/release/
   └─ Binary (rebuilt with both optimizations)
```

---

## 7. Verification Checklist (Visual)

### Verification 1: System Installation
```
✓ Step 1: libdeflate built
  /tmp/libdeflate/
  ├─ libdeflate.so     ✓ EXISTS
  └─ libdeflate.a      ✓ EXISTS

✓ Step 2: libdeflate installed
  /usr/local/lib/
  ├─ libdeflate.so     ✓ EXISTS
  ├─ libdeflate.so.1   ✓ EXISTS
  └─ libdeflate.a      ✓ EXISTS

  /usr/local/include/
  └─ libdeflate.h      ✓ EXISTS

✓ Step 3: HTSlib rebuilt
  Cargo build output: "Finished release profile"
```

### Verification 2: Code Changes
```
✓ Step 1: Cargo.toml updated
  [dependencies]
  num_cpus = "1.16"    ✓ PRESENT

✓ Step 2: bam.rs updated
  pub fn new(...) {
    let mut reader = Reader::from_path(path)?;
    let n_threads = num_cpus::get();     ✓ PRESENT
    reader.set_threads(n_threads)?;      ✓ PRESENT
  }
```

### Verification 3: Binary Verification
```
✓ Step 1: Binary linked with libdeflate
  nm sheriff-rs/target/release/*.so | grep deflate_decompress
  → Shows symbol                         ✓ FOUND

✓ Step 2: Module loads
  python: from sheriff_rs import sheriff_rs
  → No ImportError                       ✓ LOADS

✓ Step 3: CPU cores recognized
  num_cpus::get() → Returns 4 (or your CPU count)
  → set_threads(4) is called            ✓ WORKS
```

### Verification 4: Performance
```
✓ Speedup measured
  Before: 100% time (1.0x baseline)
  After:  90% time (1.11x faster)
  Improvement: +11% ✓ EXPECTED RANGE

✓ Results identical
  diff results_before.txt results_after.txt
  → No output                            ✓ IDENTICAL

✓ CPU utilization
  During BAM I/O: All cores active       ✓ CONFIRMED
```

---

## 8. Rollback Decision Tree

```
Something went wrong?
│
├─→ Did HTSlib build fail?
│   ├─→ Set environment variables & retry
│   │   export LDFLAGS="-L/usr/local/lib"
│   │   cargo clean && cargo build
│   │
│   └─→ Still fails? Uninstall libdeflate
│       sudo rm /usr/local/lib/libdeflate*
│       git checkout Cargo.toml src/bam.rs
│       cargo clean && cargo build
│
├─→ Does module not load?
│   ├─→ Check LD_LIBRARY_PATH
│   │   export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
│   │   python -c "from sheriff_rs import sheriff_rs"
│   │
│   └─→ Still fails? See above rollback
│
├─→ Is performance same as before?
│   ├─→ Check if set_threads is actually called
│   │   grep "set_threads" /home/user/Sheriff/sheriff-rs/src/bam.rs
│   │
│   ├─→ Check if test data is large enough
│   │   ls -lh test.bam  (should be >100MB)
│   │
│   └─→ If still wrong, libdeflate may not be used
│       Check: nm *.so | grep deflate_decompress
│
└─→ If unsure, full rollback
    bash /home/user/Sheriff/rollback.sh
    (Or: git checkout Cargo.toml src/bam.rs + uninstall libdeflate)
```

---

## 9. Success Criteria Summary

### ✓ Success if all present:

```
Hardware/System:
  □ Linux system (or macOS/Windows - tested on Linux)
  □ Multi-core CPU (2+ cores, ideally 4+)
  □ At least 500MB free disk space
  □ 30 minutes uninterrupted time

Installation:
  □ libdeflate installed at /usr/local/lib/libdeflate.so
  □ libdeflate.h at /usr/local/include/libdeflate.h
  □ num_cpus in Cargo.toml
  □ set_threads() in src/bam.rs

Build:
  □ cargo build --release completes without errors
  □ No major compiler warnings
  □ Binary size increases (new code included)

Verification:
  □ nm output shows deflate_decompress symbols
  □ Python module imports successfully
  □ Pipeline runs 8-11% faster
  □ Results identical to baseline

Overall:
  □ 4.4x → 4.9-5.0x speedup achieved
  □ Total time spent: ~3 hours
  □ Risk realized: Very low (no issues)
```

---

## 10. Next Steps After Phase 1

```
Phase 1 Complete ✓
    ↓
Current speedup: 4.9-5.0x ✓
    ↓
    ├─→ Goal achieved for now?
    │   └─→ Celebrate! 🎉 You saved everyone 11% of pipeline time
    │
    └─→ Want to push further?
        └─→ Phase 2: Gene UMI Counting in Rust
            ├─ Effort: 10 hours
            ├─ Benefit: +5-10% more speedup (5.5-5.8x final)
            ├─ Complexity: Medium (similar pattern to Phase 1)
            └─ When: Next sprint
```

---

## 11. File Reference Quick Links

| File | Purpose | Edit? |
|------|---------|-------|
| `/home/user/Sheriff/BAM_OPTIMIZATION_IMPLEMENTATION_PLAN.md` | Complete detailed plan | No - Reference |
| `/home/user/Sheriff/BAM_OPTIMIZATION_QUICK_START.md` | Copy-paste commands | No - Reference |
| `/home/user/Sheriff/BAM_RESEARCH_EXECUTIVE_SUMMARY.md` | Research findings | No - Reference |
| `/home/user/Sheriff/sheriff-rs/Cargo.toml` | Dependencies | **Yes** |
| `/home/user/Sheriff/sheriff-rs/src/bam.rs` | BAM processing | **Yes** |

---

## 12. Time Estimate Breakdown

```
Activity                          Duration    Typical Range
────────────────────────────────  ──────────  ─────────────
Hour 1: libdeflate
  Prerequisites check              5 min       3-7 min
  Clone & build libdeflate         20 min      15-25 min
  Install system-wide              10 min      8-12 min
  Rebuild HTSlib                   25 min      20-30 min

Hour 2: Parallel BGZF
  Update dependencies              10 min      8-12 min
  Update code                      15 min      10-20 min
  Rebuild                          20 min      15-25 min
  Verify module                    10 min      5-15 min

Hour 3: Testing
  Baseline benchmark               15 min      10-20 min
  Optimized benchmark              30 min      25-35 min
  Verify correctness               10 min      8-12 min

TOTAL                              180 min     160-200 min
                                   (3 hrs)     (2.7-3.3 hrs)
```

---

**Good luck with your implementation! 🚀**

*Visual Reference Complete - See other documentation for detailed steps*
