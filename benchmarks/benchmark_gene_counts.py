#!/usr/bin/env python3
"""
Benchmark Rust vs Python gene counting on a subset of a BAM.

- Derives a deterministic subset of barcodes/genes from the input BAM
- Runs Python reference counting and Rust (`sheriff_rs.gene_counts_py`)
- Reports wall time, RSS (best-effort), and output equality
"""

import argparse
import json
import os
import time
from pathlib import Path
from collections import Counter
from typing import List, Tuple

import numpy as np
import pysam

try:
    import psutil
except ImportError:
    psutil = None

try:
    import sheriff_rs
except ImportError as exc:
    raise SystemExit("sheriff_rs not available. Build with `cd sheriff-rs && maturin develop --release`.") from exc


DEFAULT_BAM = Path(__file__).resolve().parent.parent / "example_data" / "barcode_headAligned_anno.sorted.edit_regions_200kb.bam"


def memory_mb() -> float:
    """Best-effort RSS measurement in MB."""
    if psutil:
        return psutil.Process(os.getpid()).memory_info().rss / (1024 * 1024)
    import resource  # Lazy import to avoid platform differences
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # ru_maxrss is KB on Linux, bytes on macOS
    if os.name == "posix" and os.uname().sysname == "Darwin":
        return usage / (1024 * 1024)
    return usage / 1024


def sample_barcodes_and_genes(bam_path: Path, max_cells: int, max_genes: int, max_reads: int) -> Tuple[List[str], List[str]]:
    """Collect most frequent barcodes/genes to keep runs fast."""
    cb_counter = Counter()
    gx_counter = Counter()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for i, rec in enumerate(bam.fetch(until_eof=True)):
            if max_reads and max_reads > 0 and i >= max_reads:
                break
            try:
                cb = rec.get_tag("CB")
            except KeyError:
                continue
            cb_counter[cb] += 1

            gene_tag = None
            for tag in ("GX", "GN", "gn"):
                if rec.has_tag(tag):
                    gene_tag = rec.get_tag(tag)
                    break
            if gene_tag:
                gx_counter[gene_tag] += 1

    barcodes = [bc for bc, _ in cb_counter.most_common(max_cells)]
    genes = [g for g, _ in gx_counter.most_common() if g][:max_genes]
    return barcodes, genes


def python_gene_counts(bam_path: Path, barcodes: List[str], gene_ids: List[str], max_reads: int) -> np.ndarray:
    """Simple Python reference implementation (UMI-unique per cell/gene)."""
    bc_to_idx = {bc: i for i, bc in enumerate(barcodes)}
    gene_to_idx = {g: i for i, g in enumerate(gene_ids)}
    umi_map = {}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for i, rec in enumerate(bam.fetch(until_eof=True)):
            if max_reads and max_reads > 0 and i >= max_reads:
                break
            try:
                cb = rec.get_tag("CB")
            except KeyError:
                continue
            cell_idx = bc_to_idx.get(cb)
            if cell_idx is None:
                continue

            gene_tag = None
            for tag in ("GX", "GN", "gn"):
                if rec.has_tag(tag):
                    gene_tag = rec.get_tag(tag)
                    break
            if not gene_tag:
                continue

            gene_idx = gene_to_idx.get(gene_tag)
            if gene_idx is None:
                continue

            umi_tag = None
            for tag in ("pN", "UB"):
                if rec.has_tag(tag):
                    umi_tag = rec.get_tag(tag)
                    break
            if umi_tag is None:
                continue

            umi_map.setdefault((gene_idx, cell_idx), set()).add(umi_tag)

    counts = np.zeros((len(gene_ids), len(barcodes)), dtype=np.uint32)
    for (gene_idx, cell_idx), umi_set in umi_map.items():
        counts[gene_idx, cell_idx] = len(umi_set)
    return counts


def run_variant(name, func):
    rss_before = memory_mb()
    start = time.perf_counter()
    result = func()
    duration = time.perf_counter() - start
    rss_after = memory_mb()
    return result, {"wall_sec": duration, "rss_mb": rss_after, "rss_mb_delta": max(0.0, rss_after - rss_before)}


def main():
    parser = argparse.ArgumentParser(description="Benchmark Rust vs Python gene counting.")
    parser.add_argument("--bam", type=Path, default=DEFAULT_BAM, help="Input BAM (indexed).")
    parser.add_argument("--max-cells", type=int, default=64, help="Max barcodes to include.")
    parser.add_argument("--max-genes", type=int, default=32, help="Max genes to include.")
    parser.add_argument("--max-reads", type=int, default=0, help="Max reads to scan (0 = all).")
    parser.add_argument("--out", type=Path, default=Path("benchmarks/results_gene_counts.json"), help="Path to write JSON results.")
    args = parser.parse_args()

    if not args.bam.exists():
        raise SystemExit(f"BAM not found: {args.bam}")

    barcodes, gene_ids = sample_barcodes_and_genes(args.bam, args.max_cells, args.max_genes, args.max_reads)
    if not barcodes or not gene_ids:
        raise SystemExit("Could not derive barcodes/genes from BAM; aborting.")

    print(f"Benchmarking on {args.bam} with {len(barcodes)} barcodes, {len(gene_ids)} genes (<= {args.max_reads} reads)")

    py_counts, py_metrics = run_variant(
        "python",
        lambda: python_gene_counts(args.bam, barcodes, gene_ids, args.max_reads),
    )
    print(f"[python] {py_metrics['wall_sec']:.2f}s, rss~{py_metrics['rss_mb']:.1f} MB")

    rust_counts, rust_metrics = run_variant(
        "rust",
        lambda: np.array(sheriff_rs.gene_counts_py(str(args.bam), barcodes, gene_ids), dtype=np.uint32),
    )
    print(f"[rust]   {rust_metrics['wall_sec']:.2f}s, rss~{rust_metrics['rss_mb']:.1f} MB")

    parity = py_counts.shape == rust_counts.shape and np.array_equal(py_counts, rust_counts)
    diff = np.abs(py_counts - rust_counts)
    max_delta = diff.max() if diff.size else 0

    print(f"Parity: {parity}, max delta: {max_delta}")

    results = {
        "bam": str(args.bam),
        "barcodes": len(barcodes),
        "genes": len(gene_ids),
        "max_reads": args.max_reads,
        "python": py_metrics,
        "rust": rust_metrics,
        "parity": parity,
        "max_delta": int(max_delta),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        json.dump(results, fh, indent=2)
    print(f"Wrote results to {args.out}")


if __name__ == "__main__":
    main()
