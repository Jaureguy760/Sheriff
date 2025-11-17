#!/usr/bin/env python3
"""
Sanity check: Compare Pure Python vs Rust-only Sheriff on small dataset.
Verifies outputs match and shows timing.
"""

import subprocess
import time
import sys
from pathlib import Path
import tempfile

# Paths
RUST_SHERIFF = Path("/iblm/netapp/data3/jjaureguy/software/Sheriff")
PYTHON_SHERIFF = Path("/iblm/netapp/data3/jjaureguy/software/Sheriff_ORIGINAL")

# Small test dataset
BAM_FILE = RUST_SHERIFF / "example_data/barcode_headAligned_anno.sorted.edit_regions_200kb.bam"
REF_FASTA = Path("/iblm/netapp/data4/bbalderson/edit_capture/data/t7_indel_calling/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
BARCODE_FILE = RUST_SHERIFF / "example_data/barcode_whitelist.500-cell.txt"
REF_GTF = Path("/iblm/netapp/data4/bbalderson/edit_capture/data/t7_indel_calling/Homo_sapiens.GRCh38.110.gtf")


def run_sheriff(sheriff_dir, output_dir, label, use_run_cmd=False):
    """Run sheriff with specific version."""
    print(f"\n{'='*60}")
    print(f"Running {label} Sheriff")
    print(f"{'='*60}")

    # Build command based on version
    if use_run_cmd:
        # New version uses 'run' subcommand
        cmd = [
            sys.executable, "-m", "sheriff", "run",
            str(BAM_FILE), str(REF_FASTA), str(BARCODE_FILE), str(REF_GTF),
            "--cpu", "4",
            "--outdir", str(output_dir),
        ]
    else:
        # Old version: direct arguments
        cmd = [
            sys.executable, "-m", "sheriff",
            str(BAM_FILE), str(REF_FASTA), str(BARCODE_FILE), str(REF_GTF),
            "--cpu", "4",
            "--outdir", str(output_dir),
        ]

    env = dict(subprocess.os.environ)
    # Put our version FIRST in path to override any installed version
    existing_path = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = f"{sheriff_dir}:{existing_path}" if existing_path else str(sheriff_dir)

    print(f"PYTHONPATH: {sheriff_dir}")
    print(f"Output dir: {output_dir}")
    print(f"Command: {' '.join(cmd[:6])}...")

    start = time.time()
    result = subprocess.run(
        cmd,
        cwd=str(sheriff_dir),
        env=env,
        capture_output=True,
        text=True,
        timeout=1800,  # 30 min timeout
    )
    elapsed = time.time() - start

    if result.returncode != 0:
        print(f"ERROR (exit code {result.returncode}):")
        print("STDERR:", result.stderr[-3000:] if len(result.stderr) > 3000 else result.stderr)
        if result.stdout:
            print("STDOUT (last 500 chars):", result.stdout[-500:])
        return None, elapsed

    print(f"✅ Completed in {elapsed:.2f} seconds")

    # Show key output lines
    for line in result.stdout.split('\n'):
        if any(k in line.lower() for k in ['edit site', 'umi', 'cell', 'complete', 'process']):
            print(f"  {line.strip()}")

    return result.stdout, elapsed


def compare_outputs(dir1, dir2):
    """Compare output files from two runs."""
    print(f"\n{'='*60}")
    print("Comparing Outputs")
    print(f"{'='*60}")

    dir1 = Path(dir1)
    dir2 = Path(dir2)

    files1 = {f.name: f for f in dir1.glob("*") if f.is_file()}
    files2 = {f.name: f for f in dir2.glob("*") if f.is_file()}

    all_match = True

    for name in sorted(set(files1.keys()) | set(files2.keys())):
        if name not in files1:
            print(f"  {name}: ⚠️ Only in Rust output")
            continue
        if name not in files2:
            print(f"  {name}: ⚠️ Only in Python output")
            continue

        f1, f2 = files1[name], files2[name]
        size1, size2 = f1.stat().st_size, f2.stat().st_size

        # For CSV/text files, compare contents
        if name.endswith(('.csv', '.txt', '.tsv')):
            with open(f1) as a, open(f2) as b:
                content1 = a.read()
                content2 = b.read()

            if content1 == content2:
                print(f"  {name}: ✅ EXACT MATCH ({size1} bytes)")
            else:
                print(f"  {name}: ❌ DIFFER (Python: {size1}b, Rust: {size2}b)")
                all_match = False
                # Show first difference
                lines1 = content1.split('\n')
                lines2 = content2.split('\n')
                for i, (l1, l2) in enumerate(zip(lines1, lines2)):
                    if l1 != l2:
                        print(f"    First diff at line {i+1}:")
                        print(f"      Python: {l1[:80]}")
                        print(f"      Rust:   {l2[:80]}")
                        break
        else:
            if size1 == size2:
                print(f"  {name}: ✅ Same size ({size1} bytes)")
            else:
                print(f"  {name}: ⚠️ Size differs (Python: {size1}b, Rust: {size2}b)")

    return all_match


def main():
    print("="*60)
    print("SANITY CHECK: Python-only vs Rust-only Sheriff")
    print("="*60)
    print(f"BAM file: {BAM_FILE.name}")
    print(f"Size: {BAM_FILE.stat().st_size / 1024 / 1024:.1f} MB")
    print(f"Barcode file: {BARCODE_FILE.name}")

    # Create temp output directories
    with tempfile.TemporaryDirectory() as tmpdir:
        python_out = Path(tmpdir) / "python_output"
        rust_out = Path(tmpdir) / "rust_output"
        python_out.mkdir()
        rust_out.mkdir()

        # Run Python version (old CLI - no subcommand)
        python_stdout, python_time = run_sheriff(
            PYTHON_SHERIFF, python_out, "PURE PYTHON (original)", use_run_cmd=False
        )

        # Run Rust version (new CLI - uses 'run' subcommand)
        rust_stdout, rust_time = run_sheriff(
            RUST_SHERIFF, rust_out, "RUST-ONLY", use_run_cmd=True
        )

        # Compare outputs
        if python_stdout and rust_stdout:
            all_match = compare_outputs(python_out, rust_out)

            # Summary
            print(f"\n{'='*60}")
            print("SUMMARY")
            print(f"{'='*60}")
            print(f"Python time: {python_time:.2f}s")
            print(f"Rust time:   {rust_time:.2f}s")
            if python_time > 0 and rust_time > 0:
                speedup = python_time / rust_time
                print(f"Speedup:     {speedup:.2f}x")

            if all_match:
                print(f"\n✅ OUTPUTS MATCH - Rust produces identical results!")
            else:
                print(f"\n⚠️ OUTPUTS DIFFER - Check if differences are expected")

            # Show output files
            print(f"\nPython output: {python_out}")
            print(f"Rust output:   {rust_out}")
            print("\nOutput files:")
            for f in sorted(rust_out.glob("*")):
                print(f"  {f.name}: {f.stat().st_size} bytes")
        else:
            print("\n❌ One or both runs failed!")
            if python_stdout is None:
                print("  Python run failed")
            if rust_stdout is None:
                print("  Rust run failed")


if __name__ == "__main__":
    main()
