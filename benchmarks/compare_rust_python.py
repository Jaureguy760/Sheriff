#!/usr/bin/env python3
"""
Compare Rust BAM filter vs Python pysam iteration performance

This script benchmarks:
1. Python baseline using pysam (current implementation)
2. Rust optimized using sheriff_rs Python module

Expected speedup: 10-50x with Rust
"""
import time
import pysam
import json
import os
import sys
from pathlib import Path

# Import Sheriff's bam_utils
try:
    from sheriff.bam_utils import filter_bam_by_barcodes, HAS_RUST
except ImportError:
    print("ERROR: Could not import sheriff.bam_utils")
    print("Please install Sheriff: pip install -e .")
    sys.exit(1)

def benchmark_python_filter(bam_path, whitelist_path, output_path):
    """Python baseline using pysam"""
    print("  Loading whitelist...")
    with open(whitelist_path) as f:
        whitelist = set(line.strip() for line in f)
    print(f"  Loaded {len(whitelist):,} barcodes")

    print("  Opening BAM...")
    bam = pysam.AlignmentFile(bam_path, 'rb')
    out = pysam.AlignmentFile(output_path, 'wb', template=bam)

    reads_processed = 0
    reads_kept = 0

    print("  Filtering reads...")
    start = time.time()

    for read in bam:
        reads_processed += 1
        try:
            cb = read.get_tag('CB')
            if cb in whitelist:
                out.write(read)
                reads_kept += 1
        except KeyError:
            pass

        # Progress indicator
        if reads_processed % 100000 == 0:
            elapsed = time.time() - start
            throughput = reads_processed / elapsed
            print(f"    {reads_processed:,} reads ({throughput:,.0f} reads/sec)", end='\r')

    duration = time.time() - start

    bam.close()
    out.close()

    print()  # New line after progress

    return {
        'reads_processed': reads_processed,
        'reads_kept': reads_kept,
        'reads_rejected': reads_processed - reads_kept,
        'duration_seconds': duration,
        'throughput_reads_per_sec': reads_processed / duration
    }

def benchmark_rust_filter(bam_path, whitelist_path, output_path):
    """Rust implementation using sheriff_rs Python module"""
    if not HAS_RUST:
        print(f"  ERROR: Rust module not available")
        print(f"  Please install: cd sheriff-rs && maturin build --release --features python")
        print(f"  Then: pip install target/wheels/sheriff_rs-*.whl")
        return None

    print("  Running Rust filter via Python module...")

    # Load whitelist
    with open(whitelist_path) as f:
        whitelist = set(line.strip() for line in f)

    start = time.time()

    # Use sheriff.bam_utils which automatically uses Rust
    result = filter_bam_by_barcodes(
        bam_path,
        output_path,
        whitelist,
        use_rust=True,
        verbose=False  # Don't print progress
    )

    duration = time.time() - start
    result['duration_seconds'] = duration
    result['throughput_reads_per_sec'] = result['reads_processed'] / duration

    return result

def format_number(n):
    """Format number with commas"""
    if isinstance(n, float):
        return f"{n:,.2f}"
    return f"{n:,}"

def print_results(results, label):
    """Pretty print benchmark results"""
    print(f"\n  Results:")
    print(f"    Reads processed: {format_number(results['reads_processed'])}")
    print(f"    Reads kept:      {format_number(results['reads_kept'])} ({100*results['reads_kept']/results['reads_processed']:.1f}%)")
    print(f"    Reads rejected:  {format_number(results['reads_rejected'])} ({100*results['reads_rejected']/results['reads_processed']:.1f}%)")
    print(f"    Duration:        {results['duration_seconds']:.3f}s")
    print(f"    Throughput:      {format_number(results['throughput_reads_per_sec'])} reads/sec")

if __name__ == '__main__':
    # Check for example data
    example_bam = 'example_data/barcode_headAligned_anno.sorted.edit_regions_200kb.bam'
    example_whitelist = 'example_data/barcode_whitelist.500-cell.txt'

    if not os.path.exists(example_bam):
        print(f"ERROR: Example BAM not found: {example_bam}")
        print("Please ensure example data is present")
        sys.exit(1)

    if not os.path.exists(example_whitelist):
        print(f"ERROR: Whitelist not found: {example_whitelist}")
        sys.exit(1)

    print("=" * 70)
    print("Sheriff BAM Filter Performance Comparison")
    print("=" * 70)
    print()
    print(f"Test data:")
    print(f"  BAM file:  {example_bam}")
    print(f"  Whitelist: {example_whitelist}")
    print()

    # Python baseline
    print("[1/2] Running Python (pysam) baseline...")
    python_output = '/tmp/python_output.bam'
    python_results = benchmark_python_filter(
        example_bam, example_whitelist, python_output
    )
    print_results(python_results, "Python")

    # Rust optimized
    print()
    print("[2/2] Running Rust (sheriff-rs) optimized...")
    rust_output = '/tmp/rust_output.bam'
    rust_results = benchmark_rust_filter(
        example_bam, example_whitelist, rust_output
    )

    if rust_results is None:
        print("\n⚠️  Rust benchmark failed. See errors above.")
        sys.exit(1)

    print_results(rust_results, "Rust")

    # Comparison
    speedup = python_results['duration_seconds'] / rust_results['duration_seconds']
    time_saved = python_results['duration_seconds'] - rust_results['duration_seconds']

    print()
    print("=" * 70)
    print("COMPARISON")
    print("=" * 70)
    print(f"  Speedup:     {speedup:.1f}x faster with Rust")
    print(f"  Time saved:  {time_saved:.2f}s ({100*time_saved/python_results['duration_seconds']:.1f}% reduction)")
    print()

    if speedup >= 10:
        print("  ✅ SUCCESS: Rust achieves target performance (≥10x speedup)!")
    elif speedup >= 5:
        print(f"  ⚠️  PARTIAL: {speedup:.1f}x speedup (target: 10x, but close!)")
    else:
        print(f"  ❌ MISS: Only {speedup:.1f}x speedup (target: 10x)")
        print("     Consider profiling to identify bottlenecks")

    print()

    # Verify outputs match
    print("Verifying outputs match...")
    python_count = sum(1 for _ in pysam.AlignmentFile(python_output, 'rb'))
    rust_count = sum(1 for _ in pysam.AlignmentFile(rust_output, 'rb'))

    if python_count == rust_count:
        print(f"  ✅ Outputs match: {python_count:,} reads in both files")
    else:
        print(f"  ❌ MISMATCH: Python={python_count:,}, Rust={rust_count:,}")

    # Save results
    results = {
        'test_data': {
            'bam': example_bam,
            'whitelist': example_whitelist,
        },
        'python': python_results,
        'rust': rust_results,
        'comparison': {
            'speedup': speedup,
            'time_saved_seconds': time_saved,
            'target_met': speedup >= 10,
        }
    }

    output_file = 'benchmarks/rust_python_comparison.json'
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)

    print()
    print(f"Results saved to: {output_file}")
    print()

    # Cleanup
    for f in [python_output, rust_output]:
        if os.path.exists(f):
            os.unlink(f)

    print("=" * 70)
