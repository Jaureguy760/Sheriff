#!/usr/bin/env python3
"""
Apples-to-Apples Benchmark: Original Sheriff vs Rust-Accelerated Sheriff

This script compares:
- Sheriff_ORIGINAL (pure Python, v1.1.3)
- Sheriff (Rust-accelerated, v1.2.0)

Running on the same dataset with memory/time profiling.
"""

import sys
import time
import subprocess
import tracemalloc
import psutil
import json
import csv
from pathlib import Path
from datetime import datetime

# Configuration
ORIGINAL_DIR = Path("/iblm/netapp/data3/jjaureguy/software/Sheriff_ORIGINAL")
RUST_DIR = Path("/iblm/netapp/data3/jjaureguy/software/Sheriff")
DATA_DIR = Path("/iblm/netapp/data3/jjaureguy/software/Sheriff/full_dataset")

# Test datasets (start small, work up)
DATASETS = {
    "chr21_1pct": DATA_DIR / "chr21_1pct.bam",
    "chr21_5pct": DATA_DIR / "chr21_5pct.bam",
    "chr21_20pct": DATA_DIR / "chr21_20pct.bam",
    # "full_114GB": DATA_DIR / "full_114GB.bam",  # uncomment when ready
}

# Reference files (shared)
REF_FASTA = Path("/iblm/netapp/data3/jjaureguy/software/bbalderson/reference/reference.fa")
REF_GTF = Path("/iblm/netapp/data3/jjaureguy/software/bbalderson/reference/gencode.v39.annotation.gtf")
BLACKLIST = Path("/iblm/netapp/data3/jjaureguy/software/bbalderson/references/hg38-blacklist.v2.bed")

# Output directory
OUTPUT_DIR = Path("/iblm/netapp/data3/jjaureguy/software/Sheriff/benchmark_results")
OUTPUT_DIR.mkdir(exist_ok=True)


def get_memory_mb():
    """Get current memory usage in MB."""
    process = psutil.Process()
    return process.memory_info().rss / 1024 / 1024


def run_sheriff_version(version_dir, bam_file, output_prefix, cpu=4):
    """
    Run a specific version of Sheriff.

    Returns dict with timing and output paths.
    """
    version_name = "ORIGINAL" if "ORIGINAL" in str(version_dir) else "RUST"
    print(f"\n{'='*60}")
    print(f"Running {version_name} Sheriff on {bam_file.name}")
    print(f"{'='*60}")

    # Prepare output directory
    out_dir = OUTPUT_DIR / f"{output_prefix}_{version_name}"
    out_dir.mkdir(exist_ok=True)

    # Build command - run sheriff directly
    cmd = [
        sys.executable,
        "-m", "sheriff",
        "count_t7",
        str(bam_file),
        str(REF_FASTA),
        str(REF_GTF),
        str(out_dir),
        "--cpu", str(cpu),
    ]

    # Add blacklist if exists
    if BLACKLIST.exists():
        cmd.extend(["--blacklist", str(BLACKLIST)])

    # Set PYTHONPATH to use specific version
    env = {
        **dict(list(subprocess.os.environ.items())),
        "PYTHONPATH": str(version_dir),
    }

    # If running ORIGINAL, disable Rust by not having the module
    if version_name == "ORIGINAL":
        # Original doesn't have rust_accelerated module, so it will use pure Python
        pass
    else:
        # For Rust version, ensure it uses Rust
        env["USE_RUST_UMI"] = "1"
        env["USE_RUST_EDIT"] = "1"

    print(f"Command: {' '.join(cmd)}")
    print(f"PYTHONPATH: {env['PYTHONPATH']}")

    # Track time and memory
    start_time = time.time()
    start_mem = get_memory_mb()
    peak_mem = start_mem

    # Run process and monitor
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=env,
        cwd=str(version_dir),
    )

    # Monitor output and memory
    output_lines = []
    while True:
        line = process.stdout.readline()
        if not line and process.poll() is not None:
            break
        if line:
            print(f"  [{version_name}] {line.rstrip()}")
            output_lines.append(line)
            # Update peak memory
            current_mem = get_memory_mb()
            peak_mem = max(peak_mem, current_mem)

    # Wait for completion
    return_code = process.wait()
    end_time = time.time()
    end_mem = get_memory_mb()

    elapsed = end_time - start_time

    result = {
        "version": version_name,
        "bam_file": str(bam_file),
        "output_dir": str(out_dir),
        "elapsed_seconds": elapsed,
        "start_mem_mb": start_mem,
        "end_mem_mb": end_mem,
        "peak_mem_mb": peak_mem,
        "return_code": return_code,
        "output_lines": output_lines,
    }

    print(f"\n{version_name} completed in {elapsed:.2f}s (return code: {return_code})")
    print(f"Memory: start={start_mem:.1f}MB, peak={peak_mem:.1f}MB, end={end_mem:.1f}MB")

    return result


def compare_outputs(original_dir, rust_dir):
    """
    Compare output files between original and Rust versions.

    Returns dict with comparison results.
    """
    print(f"\n{'='*60}")
    print("Comparing outputs for correctness")
    print(f"{'='*60}")

    original_path = Path(original_dir)
    rust_path = Path(rust_dir)

    comparisons = {}

    # Compare all output files
    for orig_file in original_path.glob("*"):
        if orig_file.is_file():
            rust_file = rust_path / orig_file.name

            if rust_file.exists():
                # Compare file sizes
                orig_size = orig_file.stat().st_size
                rust_size = rust_file.stat().st_size

                # Compare contents (for text files)
                if orig_file.suffix in [".csv", ".txt", ".tsv"]:
                    with open(orig_file) as f1, open(rust_file) as f2:
                        orig_content = f1.read()
                        rust_content = f2.read()
                    content_match = orig_content == rust_content
                else:
                    content_match = orig_size == rust_size

                comparisons[orig_file.name] = {
                    "orig_size": orig_size,
                    "rust_size": rust_size,
                    "size_match": orig_size == rust_size,
                    "content_match": content_match,
                }

                status = "✓ MATCH" if content_match else "✗ DIFFER"
                print(f"  {orig_file.name}: {status}")
                if not content_match:
                    print(f"    Original: {orig_size} bytes, Rust: {rust_size} bytes")
            else:
                comparisons[orig_file.name] = {"missing_in_rust": True}
                print(f"  {orig_file.name}: ✗ Missing in Rust output")

    # Check for extra files in Rust output
    for rust_file in rust_path.glob("*"):
        if rust_file.is_file() and rust_file.name not in comparisons:
            comparisons[rust_file.name] = {"extra_in_rust": True}
            print(f"  {rust_file.name}: ⚠ Extra file in Rust output")

    return comparisons


def run_benchmark(dataset_name, bam_file, cpu=4):
    """
    Run complete benchmark for a single dataset.
    """
    print(f"\n{'#'*70}")
    print(f"# BENCHMARK: {dataset_name}")
    print(f"# File: {bam_file}")
    print(f"# CPUs: {cpu}")
    print(f"{'#'*70}")

    # Get file size
    file_size_gb = bam_file.stat().st_size / (1024**3)
    print(f"File size: {file_size_gb:.2f} GB")

    results = {
        "dataset": dataset_name,
        "file_size_gb": file_size_gb,
        "cpu": cpu,
        "timestamp": datetime.now().isoformat(),
    }

    # Run original version
    try:
        original_result = run_sheriff_version(ORIGINAL_DIR, bam_file, dataset_name, cpu)
        results["original"] = original_result
    except Exception as e:
        print(f"ERROR running original: {e}")
        traceback.print_exc()
        results["original"] = {"error": str(e)}
        return results

    # Run Rust version
    try:
        rust_result = run_sheriff_version(RUST_DIR, bam_file, dataset_name, cpu)
        results["rust"] = rust_result
    except Exception as e:
        print(f"ERROR running Rust version: {e}")
        traceback.print_exc()
        results["rust"] = {"error": str(e)}
        return results

    # Compare outputs
    if "error" not in results["original"] and "error" not in results["rust"]:
        comparisons = compare_outputs(
            results["original"]["output_dir"],
            results["rust"]["output_dir"]
        )
        results["comparisons"] = comparisons

        # Calculate speedup
        if results["original"]["elapsed_seconds"] > 0:
            speedup = results["original"]["elapsed_seconds"] / results["rust"]["elapsed_seconds"]
            results["speedup"] = speedup
            print(f"\n🚀 SPEEDUP: {speedup:.2f}x faster with Rust")

    # Save results
    result_file = OUTPUT_DIR / f"{dataset_name}_benchmark.json"
    with open(result_file, "w") as f:
        # Remove output_lines for cleaner JSON
        clean_results = results.copy()
        if "original" in clean_results and "output_lines" in clean_results["original"]:
            del clean_results["original"]["output_lines"]
        if "rust" in clean_results and "output_lines" in clean_results["rust"]:
            del clean_results["rust"]["output_lines"]
        json.dump(clean_results, f, indent=2)
    print(f"\nResults saved to: {result_file}")

    return results


def main():
    print("="*70)
    print("APPLES-TO-APPLES BENCHMARK: Original vs Rust-Accelerated Sheriff")
    print("="*70)
    print(f"Timestamp: {datetime.now().isoformat()}")
    print(f"Original Sheriff: {ORIGINAL_DIR}")
    print(f"Rust Sheriff: {RUST_DIR}")
    print(f"Output dir: {OUTPUT_DIR}")

    # Verify setup
    print("\nVerifying setup...")
    assert ORIGINAL_DIR.exists(), f"Original dir not found: {ORIGINAL_DIR}"
    assert RUST_DIR.exists(), f"Rust dir not found: {RUST_DIR}"
    assert (ORIGINAL_DIR / "sheriff").exists(), "Original sheriff package not found"
    assert (RUST_DIR / "sheriff").exists(), "Rust sheriff package not found"

    # Check Rust availability
    sys.path.insert(0, str(RUST_DIR))
    from sheriff.rust_accelerated import RUST_AVAILABLE
    print(f"Rust available: {RUST_AVAILABLE}")
    assert RUST_AVAILABLE, "Rust not available!"

    # Run benchmarks
    all_results = []
    for name, bam_file in DATASETS.items():
        if bam_file.exists():
            results = run_benchmark(name, bam_file, cpu=8)
            all_results.append(results)
        else:
            print(f"WARNING: Dataset {name} not found at {bam_file}")

    # Summary
    print("\n" + "="*70)
    print("BENCHMARK SUMMARY")
    print("="*70)

    summary_csv = OUTPUT_DIR / "benchmark_summary.csv"
    with open(summary_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Dataset", "File Size (GB)", "Original Time (s)", "Rust Time (s)",
            "Speedup", "Original Memory (MB)", "Rust Memory (MB)", "Output Match"
        ])

        for r in all_results:
            if "error" not in r.get("original", {}) and "error" not in r.get("rust", {}):
                orig_time = r["original"]["elapsed_seconds"]
                rust_time = r["rust"]["elapsed_seconds"]
                speedup = r.get("speedup", 0)
                orig_mem = r["original"]["peak_mem_mb"]
                rust_mem = r["rust"]["peak_mem_mb"]

                # Check if outputs match
                all_match = all(
                    c.get("content_match", False)
                    for c in r.get("comparisons", {}).values()
                    if isinstance(c, dict) and "content_match" in c
                )

                writer.writerow([
                    r["dataset"], f"{r['file_size_gb']:.2f}",
                    f"{orig_time:.2f}", f"{rust_time:.2f}",
                    f"{speedup:.2f}x", f"{orig_mem:.1f}", f"{rust_mem:.1f}",
                    "✓" if all_match else "✗"
                ])

                print(f"{r['dataset']}: {speedup:.2f}x speedup ({orig_time:.2f}s -> {rust_time:.2f}s)")
            else:
                print(f"{r['dataset']}: ERROR")

    print(f"\nSummary saved to: {summary_csv}")
    print("\nBenchmark complete!")


if __name__ == "__main__":
    main()
