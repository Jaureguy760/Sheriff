# Sheriff Benchmarking Suite

Comprehensive performance benchmarking and profiling tools for Sheriff optimization.

## Quick Start

### 1. Run Baseline Benchmarks (Before Optimization)

```bash
# Run all benchmarks and save baseline
python benchmarks/benchmark_sheriff.py --mode baseline --output benchmarks/baseline_results.json

# This will test:
# - K-mer matching performance
# - UMI deduplication (slow vs fast)
# - Full pipeline on example data (if available)
```

### 2. Make Code Optimizations

Edit the source code with your optimizations.

### 3. Run Post-Optimization Benchmarks

```bash
# Run benchmarks on optimized code
python benchmarks/benchmark_sheriff.py --mode optimized --output benchmarks/optimized_results.json
```

### 4. Compare Results

```bash
# Generate comparison report
python benchmarks/benchmark_sheriff.py --mode compare \
    --baseline-file benchmarks/baseline_results.json \
    --optimized-file benchmarks/optimized_results.json
```

Example output:
```
Performance Comparison Report
============================================================

K-MER MATCHING:
----------------------------------------
  Runtime:
    Baseline:  2.453s
    Optimized: 0.821s
    Speedup:   2.99x ✓
  Peak Memory:
    Baseline:  145.2 MB
    Optimized: 98.7 MB
    Reduction: 32.0% ✓
```

## Detailed Profiling

### Memory Profiling

```bash
# Install memory_profiler
pip install memory_profiler matplotlib

# Profile memory usage
python -m memory_profiler benchmarks/profile_memory.py

# Or with plots
mprof run benchmarks/profile_memory.py
mprof plot
```

### CPU Profiling with cProfile

```bash
# Profile full pipeline
python -m cProfile -o sheriff.prof -m sheriff \
    example_data/barcode_headAligned_anno.sorted.edit_regions_200kb.bam \
    example_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    example_data/barcode_whitelist.500-cell.txt \
    example_data/Homo_sapiens.GRCh38.110.gtf \
    -o benchmark_output

# Visualize with snakeviz
pip install snakeviz
snakeviz sheriff.prof
```

### Line-by-Line Profiling

```bash
# Install line_profiler
pip install line_profiler

# Add @profile decorator to functions you want to profile
# Then run:
kernprof -l -v your_script.py
```

## Benchmark Components

### Micro-benchmarks

**K-mer Matching** (`benchmark_kmer_matching`)
- Tests k-mer hashing and matching speed
- Measures ops/second
- Default: 10,000 iterations

**UMI Deduplication** (`benchmark_umi_deduplication`)
- Tests slow (O(n³)) vs fast (Numba) implementations
- Measures UMIs/second
- Default: 100 cells × 50 UMIs

### Full Pipeline Benchmark

**Example Data Processing** (`benchmark_full_pipeline_example_data`)
- Runs complete Sheriff workflow
- Measures end-to-end performance
- Requires example data files

## Interpreting Results

### Good Optimization Signs
- ✓ Speedup > 1.0x (faster)
- ✓ Memory reduction > 0% (less memory)
- ✓ Higher throughput (ops/sec or UMIs/sec)

### Warning Signs
- ✗ Speedup < 1.0x (slower!)
- ✗ Memory increase
- ✗ Lower throughput

## CI/CD Integration

Add to your GitHub Actions workflow:

```yaml
- name: Run Performance Benchmarks
  run: |
    python benchmarks/benchmark_sheriff.py --mode baseline
    # Store results as artifacts

- name: Compare Performance
  run: |
    python benchmarks/benchmark_sheriff.py --mode compare \
      --baseline-file benchmarks/baseline_results.json \
      --optimized-file benchmarks/optimized_results.json
```

## Custom Benchmarks

Create your own benchmarks by extending the framework:

```python
from benchmarks.benchmark_sheriff import PerformanceMonitor

def benchmark_my_function():
    monitor = PerformanceMonitor()
    monitor.start()

    # Your code here
    for i in range(1000):
        my_function()
        monitor.sample()

    metrics = monitor.stop()
    return metrics
```

## Troubleshooting

**"example_data not found"**
- Download example data first or skip full pipeline benchmark

**"Out of memory"**
- Reduce benchmark iterations
- Use smaller test datasets
- Close other applications

**Import errors**
- Ensure Sheriff is installed: `pip install -e .`
- Install benchmark dependencies: `pip install psutil memory_profiler`

## Benchmark Results Archive

Store baseline results for version tracking:

```bash
benchmarks/
├── v1.1.3_baseline.json
├── v1.2.0_optimized.json
├── v1.2.1_optimized.json
└── ...
```

## Performance Targets

Based on example data (500 cells):

| Metric | Current | Target | Status |
|--------|---------|--------|--------|
| Full pipeline runtime | ~30s | <15s | 🎯 |
| Peak memory | ~500MB | <350MB | 🎯 |
| K-mer matching | 50k ops/s | 150k ops/s | 🎯 |
| UMI dedup | 1k UMIs/s | 10k UMIs/s | ✅ (Fast version) |

## Resources

- [Python profiling guide](https://docs.python.org/3/library/profile.html)
- [memory_profiler docs](https://pypi.org/project/memory-profiler/)
- [psutil documentation](https://psutil.readthedocs.io/)
