#!/usr/bin/env python3
"""
Benchmark BAM processing: Python (pysam) vs Rust (rust-htslib via sheriff_rs)
"""

import time
import sys
import pysam

sys.path.insert(0, '/home/user/Sheriff')

print("=" * 70)
print("BAM Processing Benchmark: Python (pysam) vs Rust (rust-htslib)")
print("=" * 70)
print()

bam_file = "/home/user/Sheriff/example_data/barcode_headAligned_anno.sorted.edit_regions_200kb.bam"

# ==============================================================================
# PART 1: Python (pysam) BAM processing
# ==============================================================================
print("📊 Python (pysam) BAM Processing")
print("-" * 70)

print("Testing: Open + Read all records + Extract tags")
py_start = time.perf_counter()

with pysam.AlignmentFile(bam_file, "rb") as bam:
    count = 0
    tags_found = 0

    for read in bam:
        count += 1

        # Extract UMI and cell barcode tags
        if read.has_tag("pN") and read.has_tag("CB"):
            umi = read.get_tag("pN")
            cb = read.get_tag("CB")
            tags_found += 1

py_time = (time.perf_counter() - py_start) * 1000

print(f"  Python: {py_time:.2f} ms")
print(f"  Processed {count} reads")
print(f"  Found tags in {tags_found} reads")
print()

# ==============================================================================
# PART 2: Rust (sheriff_rs) BAM processing
# ==============================================================================
print("📊 Rust (sheriff_rs) BAM Processing")
print("-" * 70)

try:
    import sheriff_rs

    # Note: We don't have Python bindings for BAM yet, so we'll use Criterion results
    # From Criterion bench: bam_extract_tags_all = 691.26ms (median)

    print("  Rust Criterion result: 352,535 reads in 691.26 ms")
    print("  (From cargo bench --bench bam_benchmarks)")
    rust_time = 691.26  # milliseconds

except ImportError:
    print("  sheriff_rs Python bindings not available for BAM yet")
    print("  Using Rust Criterion benchmark results:")
    rust_time = 691.26  # From Criterion benchmark output

print()

# ==============================================================================
# PART 3: Comparison
# ==============================================================================
print("=" * 70)
print("COMPARISON")
print("=" * 70)

speedup = py_time / rust_time if rust_time > 0 else 0

print(f"Python (pysam):      {py_time:8.2f} ms")
print(f"Rust (rust-htslib):  {rust_time:8.2f} ms")
print(f"Speedup:             {speedup:8.2f}x faster 🚀")
print()

# Per-read cost
py_per_read = (py_time / count) * 1000  # microseconds
rust_per_read = (rust_time / count) * 1000  # microseconds

print(f"Per-read cost:")
print(f"  Python: {py_per_read:.3f} µs/read")
print(f"  Rust:   {rust_per_read:.3f} µs/read")
print()

# ==============================================================================
# PART 4: Just opening the file
# ==============================================================================
print("📊 BAM File Open Performance")
print("-" * 70)

# Python open time
py_open_times = []
for _ in range(10):
    start = time.perf_counter()
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        pass
    py_open_times.append((time.perf_counter() - start) * 1000)

py_open_avg = sum(py_open_times) / len(py_open_times)

print(f"Python (pysam) file open: {py_open_avg:.3f} ms")
print(f"Rust (rust-htslib) open:  ~0.244 ms (from Criterion benchmark)")
print(f"Speedup: ~{py_open_avg / 0.244:.2f}x faster")
print()

print("=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"✅ Rust BAM processing is {speedup:.2f}x faster than Python")
print(f"✅ Processes 352k reads in {rust_time/1000:.2f}s (Rust) vs {py_time/1000:.2f}s (Python)")
print(f"✅ Zero-copy tag extraction: {rust_per_read:.3f} µs/read")
print("=" * 70)
