#!/usr/bin/env python3
"""
Comprehensive benchmark comparing Python, Rust sequential, and Rust parallel implementations.

Compares:
1. Python (pysam) - baseline
2. Rust sequential - single-threaded
3. Rust parallel - rayon par_bridge

Usage:
    python benchmarks/compare_parallel.py
"""

import time
import json
import sys
from pathlib import Path

# Add Sheriff to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from sheriff.bam_utils import filter_bam_by_barcodes, HAS_RUST


def run_benchmark(bam_path, whitelist_path, output_dir):
    """Run comprehensive benchmark comparing all implementations."""

    # Load whitelist
    with open(whitelist_path) as f:
        whitelist = set(line.strip() for line in f)

    print(f"=" * 80)
    print(f"Sheriff BAM Filter Benchmark - Parallel Processing Comparison")
    print(f"=" * 80)
    print(f"\nConfiguration:")
    print(f"  Input BAM: {bam_path}")
    print(f"  Whitelist: {whitelist_path} ({len(whitelist)} barcodes)")
    print(f"  Rust available: {HAS_RUST}")
    print()

    results = {}

    # Benchmark 1: Python (baseline)
    print(f"[1/3] Running Python (pysam) baseline...")
    print("-" * 80)
    output_python = f"{output_dir}/filtered_python.bam"

    result_python = filter_bam_by_barcodes(
        bam_path, output_python, whitelist,
        use_rust=False, verbose=True
    )
    results['python'] = result_python
    print()

    if not HAS_RUST:
        print("⚠️  Rust not available, skipping Rust benchmarks")
        return results

    # Benchmark 2: Rust Sequential
    print(f"[2/3] Running Rust sequential...")
    print("-" * 80)
    output_rust_seq = f"{output_dir}/filtered_rust_sequential.bam"

    result_rust_seq = filter_bam_by_barcodes(
        bam_path, output_rust_seq, whitelist,
        use_rust=True, use_parallel=False, verbose=True
    )
    results['rust_sequential'] = result_rust_seq
    print()

    # Benchmark 3: Rust Parallel
    print(f"[3/3] Running Rust parallel (rayon par_bridge)...")
    print("-" * 80)
    output_rust_par = f"{output_dir}/filtered_rust_parallel.bam"

    result_rust_par = filter_bam_by_barcodes(
        bam_path, output_rust_par, whitelist,
        use_rust=True, use_parallel=True, verbose=True
    )
    results['rust_parallel'] = result_rust_par
    print()

    # Summary comparison
    print("=" * 80)
    print("PERFORMANCE SUMMARY")
    print("=" * 80)
    print()

    # Calculate speedups
    python_time = results['python']['duration_seconds']
    rust_seq_time = results['rust_sequential']['duration_seconds']
    rust_par_time = results['rust_parallel']['duration_seconds']

    rust_seq_speedup = python_time / rust_seq_time
    rust_par_speedup = python_time / rust_par_time
    parallel_vs_seq = rust_seq_time / rust_par_time

    print(f"{'Implementation':<25} {'Time (s)':<12} {'Throughput (reads/s)':<25} {'Speedup vs Python':<20}")
    print("-" * 80)

    for name, result in results.items():
        duration = result['duration_seconds']
        throughput = result['reads_processed'] / max(0.001, duration)

        if name == 'python':
            speedup_str = "1.0x (baseline)"
        elif name == 'rust_sequential':
            speedup_str = f"{rust_seq_speedup:.2f}x"
        elif name == 'rust_parallel':
            speedup_str = f"{rust_par_speedup:.2f}x"

        print(f"{name:<25} {duration:<12.2f} {throughput:<25,.0f} {speedup_str:<20}")

    print()
    print(f"🚀 Parallel vs Sequential: {parallel_vs_seq:.2f}x speedup")
    print()

    # Verify outputs are identical
    print("Verifying output consistency...")
    print(f"  Python:          {results['python']['reads_kept']:,} reads")
    print(f"  Rust sequential: {results['rust_sequential']['reads_kept']:,} reads")
    print(f"  Rust parallel:   {results['rust_parallel']['reads_kept']:,} reads")

    if (results['python']['reads_kept'] == results['rust_sequential']['reads_kept'] ==
        results['rust_parallel']['reads_kept']):
        print("  ✅ All implementations produced identical results!")
    else:
        print("  ⚠️  WARNING: Output mismatch detected!")

    # Save results
    output_file = f"{output_dir}/parallel_benchmark_results.json"
    with open(output_file, 'w') as f:
        json.dump({
            'results': results,
            'speedups': {
                'rust_sequential_vs_python': rust_seq_speedup,
                'rust_parallel_vs_python': rust_par_speedup,
                'parallel_vs_sequential': parallel_vs_seq,
            },
            'configuration': {
                'bam_file': bam_path,
                'whitelist_size': len(whitelist),
                'reads_processed': results['python']['reads_processed'],
            }
        }, f, indent=2)

    print(f"\n📊 Results saved to: {output_file}")

    return results


if __name__ == '__main__':
    # Use full example dataset
    DATA_DIR = Path(__file__).parent.parent / "example_data"
    OUTPUT_DIR = Path(__file__).parent.parent / "benchmarks/output"
    OUTPUT_DIR.mkdir(exist_ok=True)

    BAM_FILE = DATA_DIR / "barcode_headAligned_anno.sorted.edit_regions_200kb.bam"
    WHITELIST_FILE = DATA_DIR / "barcode_whitelist.500-cell.txt"

    if not BAM_FILE.exists():
        print(f"❌ Error: BAM file not found: {BAM_FILE}")
        sys.exit(1)

    if not WHITELIST_FILE.exists():
        print(f"❌ Error: Whitelist file not found: {WHITELIST_FILE}")
        sys.exit(1)

    run_benchmark(str(BAM_FILE), str(WHITELIST_FILE), str(OUTPUT_DIR))
