#!/usr/bin/env python3
"""
Sheriff Performance Benchmarking Suite

Measures performance metrics before and after optimizations:
- Runtime (total and per-stage)
- Memory usage (peak and average)
- CPU utilization
- I/O operations

Usage:
    python benchmarks/benchmark_sheriff.py --mode baseline
    python benchmarks/benchmark_sheriff.py --mode optimized
    python benchmarks/benchmark_sheriff.py --compare baseline_results.json optimized_results.json
"""

import argparse
import json
import time
import tracemalloc
import psutil
import os
import sys
from pathlib import Path
from datetime import datetime
import subprocess

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from sheriff.count_t7 import run_count_t7, get_barcoded_edits, KmerMatcher
from sheriff.helpers import deduplicate_umis, cell_umi_counts_FAST
import numpy as np
import pandas as pd


class PerformanceMonitor:
    """Monitor system resources during execution"""

    def __init__(self):
        self.process = psutil.Process()
        self.start_time = None
        self.start_memory = None
        self.peak_memory = 0
        self.cpu_percentages = []

    def start(self):
        """Start monitoring"""
        tracemalloc.start()
        self.start_time = time.time()
        self.start_memory = self.process.memory_info().rss / 1024 / 1024  # MB
        self.peak_memory = self.start_memory

    def sample(self):
        """Sample current metrics"""
        current_memory = self.process.memory_info().rss / 1024 / 1024  # MB
        self.peak_memory = max(self.peak_memory, current_memory)
        cpu_percent = self.process.cpu_percent(interval=0.1)
        self.cpu_percentages.append(cpu_percent)

    def stop(self):
        """Stop monitoring and return metrics"""
        end_time = time.time()
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        return {
            'runtime_seconds': end_time - self.start_time,
            'memory_start_mb': self.start_memory,
            'memory_peak_mb': self.peak_memory,
            'memory_traced_peak_mb': peak / 1024 / 1024,
            'cpu_avg_percent': np.mean(self.cpu_percentages) if self.cpu_percentages else 0,
            'cpu_max_percent': max(self.cpu_percentages) if self.cpu_percentages else 0,
        }


def benchmark_kmer_matching(n_iterations=10000):
    """Benchmark k-mer matching operations"""
    print("Benchmarking k-mer matching...")

    # Setup
    k = 6
    t7_barcode = "GGGAGAGTAT"
    test_sequence = "ATCGATCGATCGGGGAGAGTATATCGATCGATCG" * 10

    matcher = KmerMatcher(k, t7_barcode)

    monitor = PerformanceMonitor()
    monitor.start()

    # Benchmark
    results = []
    for i in range(n_iterations):
        # Simulate k-mer to num conversion
        kmers = [test_sequence[j:j+k] for j in range(len(test_sequence) - k + 1)]
        for kmer in kmers:
            try:
                num = matcher.kmer_to_num(kmer)
                results.append(num)
            except KeyError:
                pass

        if i % 1000 == 0:
            monitor.sample()

    metrics = monitor.stop()
    metrics['operations'] = n_iterations
    metrics['ops_per_second'] = n_iterations / metrics['runtime_seconds']

    return metrics


def benchmark_umi_deduplication(n_cells=100, umis_per_cell=50):
    """Benchmark UMI deduplication"""
    print("Benchmarking UMI deduplication...")

    # Generate test UMI data
    umi_length = 10
    test_umis = []
    for _ in range(n_cells):
        umis = set()
        for _ in range(umis_per_cell):
            umi = ''.join(np.random.choice(['A', 'C', 'G', 'T'], umi_length))
            umis.add(umi)
        test_umis.append(umis)

    monitor = PerformanceMonitor()
    monitor.start()

    # Benchmark
    total_deduped = 0
    for umi_set in test_umis:
        # Test slow version
        deduped = deduplicate_umis(umi_set)
        total_deduped += len(deduped)
        monitor.sample()

    metrics = monitor.stop()
    metrics['cells_processed'] = n_cells
    metrics['total_umis'] = n_cells * umis_per_cell
    metrics['umis_per_second'] = (n_cells * umis_per_cell) / metrics['runtime_seconds']

    return metrics


def benchmark_umi_deduplication_fast(n_cells=100, umis_per_cell=50):
    """Benchmark fast UMI deduplication with Numba"""
    print("Benchmarking FAST UMI deduplication...")

    # Generate test UMI data
    umi_length = 10
    cell_bc_indexes = np.arange(n_cells, dtype=np.uint32)
    cell_umis = []

    for _ in range(n_cells):
        umis = []
        for _ in range(umis_per_cell):
            umi = ''.join(np.random.choice(['A', 'C', 'G', 'T'], umi_length))
            umis.append(umi)
        cell_umis.append(np.array(umis, dtype=f'<U{umi_length}'))

    monitor = PerformanceMonitor()
    monitor.start()

    # Benchmark
    counts = cell_umi_counts_FAST(n_cells, cell_bc_indexes, cell_umis)
    monitor.sample()

    metrics = monitor.stop()
    metrics['cells_processed'] = n_cells
    metrics['total_umis'] = n_cells * umis_per_cell
    metrics['umis_per_second'] = (n_cells * umis_per_cell) / metrics['runtime_seconds']

    return metrics


def benchmark_full_pipeline_example_data():
    """Benchmark full Sheriff pipeline on example data"""
    print("Benchmarking full pipeline on example data...")

    # Check if example data exists
    example_dir = Path("example_data")
    if not example_dir.exists():
        print("Warning: example_data directory not found, skipping full pipeline benchmark")
        return None

    # Setup paths
    bam_file = "example_data/barcode_headAligned_anno.sorted.edit_regions_200kb.bam"
    ref_file = "example_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    barcode_file = "example_data/barcode_whitelist.500-cell.txt"
    gtf_file = "example_data/Homo_sapiens.GRCh38.110.gtf"

    # Check if files exist
    required_files = [bam_file, barcode_file]
    for file in required_files:
        if not Path(file).exists():
            print(f"Warning: {file} not found, skipping full pipeline benchmark")
            return None

    output_dir = f"benchmark_output_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

    monitor = PerformanceMonitor()
    monitor.start()

    try:
        # Run Sheriff
        run_count_t7(
            bam_file=bam_file,
            ref_file=ref_file,
            barcode_file=barcode_file,
            gtf_file=gtf_file if Path(gtf_file).exists() else None,
            t7_barcode="GGGAGAGTAT",
            k=6,
            edit_dist=140,
            edit_site_min_cells=1,
            outdir=output_dir,
            verbosity=0,
            n_cpus=1,
        )

        monitor.sample()
        success = True
    except Exception as e:
        print(f"Error running full pipeline: {e}")
        success = False

    metrics = monitor.stop()
    metrics['success'] = success

    # Cleanup
    if Path(output_dir).exists():
        import shutil
        shutil.rmtree(output_dir)

    return metrics


def run_all_benchmarks(output_file=None):
    """Run all benchmarks and save results"""
    print("="*60)
    print("Sheriff Performance Benchmark Suite")
    print("="*60)

    results = {
        'timestamp': datetime.now().isoformat(),
        'system': {
            'cpu_count': psutil.cpu_count(),
            'total_memory_gb': psutil.virtual_memory().total / 1024**3,
            'python_version': sys.version,
        },
        'benchmarks': {}
    }

    # Run micro-benchmarks
    print("\n--- Micro-benchmarks ---")
    results['benchmarks']['kmer_matching'] = benchmark_kmer_matching(n_iterations=10000)
    print(f"  K-mer matching: {results['benchmarks']['kmer_matching']['ops_per_second']:.0f} ops/sec")

    results['benchmarks']['umi_dedup_slow'] = benchmark_umi_deduplication(n_cells=100, umis_per_cell=50)
    print(f"  UMI dedup (slow): {results['benchmarks']['umi_dedup_slow']['umis_per_second']:.0f} UMIs/sec")

    results['benchmarks']['umi_dedup_fast'] = benchmark_umi_deduplication_fast(n_cells=100, umis_per_cell=50)
    print(f"  UMI dedup (fast): {results['benchmarks']['umi_dedup_fast']['umis_per_second']:.0f} UMIs/sec")

    # Calculate speedup
    speedup = (results['benchmarks']['umi_dedup_fast']['umis_per_second'] /
               results['benchmarks']['umi_dedup_slow']['umis_per_second'])
    print(f"  UMI dedup speedup: {speedup:.1f}x")

    # Run full pipeline if data available
    print("\n--- Full Pipeline Benchmark ---")
    full_metrics = benchmark_full_pipeline_example_data()
    if full_metrics:
        results['benchmarks']['full_pipeline'] = full_metrics
        print(f"  Runtime: {full_metrics['runtime_seconds']:.1f} seconds")
        print(f"  Peak Memory: {full_metrics['memory_peak_mb']:.1f} MB")

    # Save results
    if output_file:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nResults saved to: {output_file}")

    print("\n" + "="*60)
    return results


def compare_results(baseline_file, optimized_file):
    """Compare baseline and optimized benchmark results"""
    print("="*60)
    print("Performance Comparison Report")
    print("="*60)

    with open(baseline_file) as f:
        baseline = json.load(f)

    with open(optimized_file) as f:
        optimized = json.load(f)

    print(f"\nBaseline: {baseline['timestamp']}")
    print(f"Optimized: {optimized['timestamp']}")
    print("\n" + "-"*60)

    # Compare each benchmark
    for benchmark_name in baseline['benchmarks']:
        if benchmark_name not in optimized['benchmarks']:
            continue

        base = baseline['benchmarks'][benchmark_name]
        opt = optimized['benchmarks'][benchmark_name]

        print(f"\n{benchmark_name.upper().replace('_', ' ')}:")
        print("-" * 40)

        # Runtime comparison
        if 'runtime_seconds' in base and 'runtime_seconds' in opt:
            speedup = base['runtime_seconds'] / opt['runtime_seconds']
            print(f"  Runtime:")
            print(f"    Baseline:  {base['runtime_seconds']:.3f}s")
            print(f"    Optimized: {opt['runtime_seconds']:.3f}s")
            print(f"    Speedup:   {speedup:.2f}x {'✓' if speedup > 1 else '✗'}")

        # Memory comparison
        if 'memory_peak_mb' in base and 'memory_peak_mb' in opt:
            reduction = (1 - opt['memory_peak_mb'] / base['memory_peak_mb']) * 100
            print(f"  Peak Memory:")
            print(f"    Baseline:  {base['memory_peak_mb']:.1f} MB")
            print(f"    Optimized: {opt['memory_peak_mb']:.1f} MB")
            print(f"    Reduction: {reduction:.1f}% {'✓' if reduction > 0 else '✗'}")

        # Throughput comparison
        if 'ops_per_second' in base and 'ops_per_second' in opt:
            improvement = (opt['ops_per_second'] / base['ops_per_second'] - 1) * 100
            print(f"  Throughput:")
            print(f"    Baseline:  {base['ops_per_second']:.0f} ops/sec")
            print(f"    Optimized: {opt['ops_per_second']:.0f} ops/sec")
            print(f"    Improvement: {improvement:+.1f}% {'✓' if improvement > 0 else '✗'}")

        if 'umis_per_second' in base and 'umis_per_second' in opt:
            improvement = (opt['umis_per_second'] / base['umis_per_second'] - 1) * 100
            print(f"  UMI Throughput:")
            print(f"    Baseline:  {base['umis_per_second']:.0f} UMIs/sec")
            print(f"    Optimized: {opt['umis_per_second']:.0f} UMIs/sec")
            print(f"    Improvement: {improvement:+.1f}% {'✓' if improvement > 0 else '✗'}")

    print("\n" + "="*60)


def main():
    parser = argparse.ArgumentParser(description='Sheriff Performance Benchmarking')
    parser.add_argument('--mode', choices=['baseline', 'optimized', 'compare'],
                       default='baseline', help='Benchmark mode')
    parser.add_argument('--output', help='Output JSON file for results')
    parser.add_argument('--baseline-file', help='Baseline results file for comparison')
    parser.add_argument('--optimized-file', help='Optimized results file for comparison')

    args = parser.parse_args()

    if args.mode == 'compare':
        if not args.baseline_file or not args.optimized_file:
            print("Error: --baseline-file and --optimized-file required for comparison")
            sys.exit(1)
        compare_results(args.baseline_file, args.optimized_file)
    else:
        output_file = args.output or f'benchmarks/{args.mode}_results.json'
        run_all_benchmarks(output_file)


if __name__ == '__main__':
    main()
