#!/usr/bin/env python3
"""
One-command runner for Sheriff performance + validation with labeled plots.

Steps:
1) Optional baseline (Rust disabled) perf suite
2) Rust-enabled perf suite (default)
3) Validation suite (correctness)
4) Optional plots (PDF/PNG)

Examples:
  # Default: bundled example data, Rust on, validation, plots
  python benchmarks/run_perf_and_validation.py

  # With your BAM + whitelist, all reads, 256 cells/128 genes
  python benchmarks/run_perf_and_validation.py \
    --dataset mid=/data/mid.bam \
    --whitelist /data/whitelist.txt \
    --max-cells 256 --max-genes 128 --max-reads 0

  # Also run a baseline (Rust disabled) for comparison
  python benchmarks/run_perf_and_validation.py --baseline
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

from benchmarks import run_perf_suite, plot_perf


def run_validation(skip: bool) -> str:
    if skip:
        return "SKIPPED"
    print("[validation] Running validate_rust_correctness.py ...")
    res = subprocess.run(
        [sys.executable, "validate_rust_correctness.py"],
        cwd=Path(__file__).resolve().parents[1],
        capture_output=True,
        text=True,
    )
    out_file = Path("benchmarks/validation_output.txt")
    out_file.write_text(res.stdout + "\n" + res.stderr)
    status = "PASS" if res.returncode == 0 else "FAIL"
    print(f"[validation] Status: {status} (log: {out_file})")
    return status


def run_suite(out_json: Path, out_png: Path, dataset_args, baseline: bool):
    # Set/clear env toggle for Rust
    if baseline:
        os.environ["SHERIFF_DISABLE_RUST"] = "1"
        print("[suite] Running with SHERIFF_DISABLE_RUST=1 (baseline).")
    else:
        os.environ.pop("SHERIFF_DISABLE_RUST", None)

    # Build arg list for run_perf_suite.main-compatible argv
    argv = []
    for item in dataset_args.dataset or []:
        argv.extend(["--dataset", item])
    if dataset_args.whitelist:
        argv.extend(["--whitelist", dataset_args.whitelist])
    argv.extend(
        [
            "--max-cells",
            str(dataset_args.max_cells),
            "--max-genes",
            str(dataset_args.max_genes),
            "--max-reads",
            str(dataset_args.max_reads),
            "--out-json",
            str(out_json),
            "--out-png",
            str(out_png),
        ]
    )
    if dataset_args.no_plot:
        argv.append("--no-plot")

    # Invoke run_perf_suite.main() with synthetic argv
    sys_argv_backup = sys.argv
    try:
        sys.argv = ["run_perf_suite.py"] + argv
        run_perf_suite.main()
    finally:
        sys.argv = sys_argv_backup


def main():
    parser = argparse.ArgumentParser(description="Run perf + validation + plots.")
    parser.add_argument(
        "--dataset",
        action="append",
        help="Dataset spec name=/path/to.bam (repeat for multiple). Defaults to bundled example.",
    )
    parser.add_argument("--whitelist", help="Whitelist path (applied to all datasets if provided).")
    parser.add_argument("--max-cells", type=int, default=64)
    parser.add_argument("--max-genes", type=int, default=32)
    parser.add_argument("--max-reads", type=int, default=0, help="0 = all reads")
    parser.add_argument("--baseline", action="store_true", help="Also run Rust-disabled baseline suite.")
    parser.add_argument("--no-plot", action="store_true", help="Skip plot generation.")
    parser.add_argument("--skip-validation", action="store_true", help="Skip correctness validation.")
    args = parser.parse_args()

    out_dir = Path("benchmarks")
    out_dir.mkdir(exist_ok=True)

    # 1) Validation
    validation_status = run_validation(skip=args.skip_validation)

    # 2) Baseline suite (optional)
    if args.baseline:
        run_suite(
            out_json=out_dir / "results_perf_suite_baseline.json",
            out_png=out_dir / "results_perf_suite_baseline.png",
            dataset_args=args,
            baseline=True,
        )

    # 3) Rust-enabled suite
    run_suite(
        out_json=out_dir / "results_perf_suite.json",
        out_png=out_dir / "results_perf_suite.png",
        dataset_args=args,
        baseline=False,
    )

    # 4) Plots (reuse plot_perf to generate PDF/PNG from the rust suite JSON)
    if not args.no_plot:
        # plot_perf reads benchmarks/results_perf_suite.json by default; ensure it’s present
        sys_argv_backup = sys.argv
        try:
            sys.argv = ["plot_perf.py"]
            plot_perf.main()
        finally:
            sys.argv = sys_argv_backup

    # 5) Summary print
    print("\n=== Summary ===")
    print(f"Validation: {validation_status}")
    if args.baseline:
        print("Baseline JSON: benchmarks/results_perf_suite_baseline.json")
        print("Baseline PNG:  benchmarks/results_perf_suite_baseline.png")
    print("Rust JSON:     benchmarks/results_perf_suite.json")
    print("Rust plots:    benchmarks/results_perf_suite.pdf / .png")


if __name__ == "__main__":
    main()
