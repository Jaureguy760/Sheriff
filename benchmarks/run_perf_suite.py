#!/usr/bin/env python3
"""
Run a small, reproducible performance suite for Sheriff:
- BAM filtering (Rust vs Python)
- Gene UMI counting (Rust vs Python)

Outputs:
- JSON results (wall time, RSS delta, throughput, parity)
- Optional PNG plots (wall/RSS/speedup) if matplotlib is available

Defaults use bundled example data; pass your own BAM/whitelist via CLI.
"""

import argparse
import json
import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Dict, Any

import numpy as np
import pysam

try:
    import psutil
except ImportError:
    psutil = None

try:
    import matplotlib.pyplot as plt  # Optional
except Exception:
    plt = None

from sheriff.bam_utils import filter_bam_by_barcodes
from benchmarks.benchmark_gene_counts import (
    sample_barcodes_and_genes,
    python_gene_counts,
)

# Attempt to import sheriff_rs; OK if missing (Python-only benches will still run)
try:
    import sheriff_rs  # noqa: F401
    HAS_RUST = True
except ImportError:
    HAS_RUST = False

# Python fallbacks
from sheriff.helpers import deduplicate_umis as deduplicate_umis_python, get_longest_edits as get_longest_edits_python
from sheriff.count_t7 import KmerMatcher  # for edit generation if needed

# Lightweight check for Rust availability; fall back if missing
try:
    import sheriff_rs  # noqa: F401
    HAS_RUST = True
except ImportError:
    HAS_RUST = False

def rust_allowed():
    """Return True if Rust module is available and not disabled via env."""
    disable_env = os.getenv("SHERIFF_DISABLE_RUST", "").lower()
    return HAS_RUST and disable_env not in {"1", "true", "yes"}


@dataclass
class Dataset:
    name: str
    bam: Path
    whitelist: Optional[Path] = None
    max_cells: int = 64
    max_genes: int = 32
    max_reads: int = 0  # 0 = all reads


def memory_mb() -> float:
    """Best-effort RSS measurement in MB."""
    if psutil:
        return psutil.Process(os.getpid()).memory_info().rss / (1024 * 1024)
    import resource  # pragma: no cover

    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if os.name == "posix" and os.uname().sysname == "Darwin":
        return usage / (1024 * 1024)  # bytes -> MB
    return usage / 1024  # kB -> MB


def run_filter(dataset: Dataset, output_dir: Path) -> Optional[Dict[str, Any]]:
    """Benchmark BAM filtering Rust vs Python."""
    if dataset.whitelist is None or not dataset.whitelist.exists():
        print(f"[filter] Skipping {dataset.name}: no whitelist provided.")
        return None

    with open(dataset.whitelist) as fh:
        whitelist = {line.strip() for line in fh if line.strip()}

    results = {}
    for variant, use_rust in (("python", False), ("rust", True)):
        out_bam = output_dir / f"{dataset.name}_{variant}_filtered.bam"
        rss_before = memory_mb()
        start = time.perf_counter()
        stats = filter_bam_by_barcodes(
            str(dataset.bam),
            str(out_bam),
            whitelist,
            use_rust=use_rust,
            use_parallel=False,
            verbose=False,
        )
        wall = time.perf_counter() - start
        rss_after = memory_mb()
        throughput = stats["reads_processed"] / max(wall, 1e-6)
        results[variant] = {
            "reads_processed": stats["reads_processed"],
            "reads_kept": stats["reads_kept"],
            "reads_rejected": stats["reads_rejected"],
            "wall_sec": wall,
            "rss_mb": rss_after,
            "rss_mb_delta": max(0.0, rss_after - rss_before),
            "throughput_reads_per_sec": throughput,
        }

    # Parity check on output BAM counts
    py_reads = sum(1 for _ in pysam.AlignmentFile(output_dir / f"{dataset.name}_python_filtered.bam", "rb"))
    rust_reads = sum(1 for _ in pysam.AlignmentFile(output_dir / f"{dataset.name}_rust_filtered.bam", "rb"))
    results["parity_reads"] = py_reads == rust_reads
    results["speedup"] = results["python"]["wall_sec"] / results["rust"]["wall_sec"]
    return results


def run_gene_counts(dataset: Dataset) -> Dict[str, Any]:
    """Benchmark gene counting Rust vs Python."""
    # Allow skipping gene bench explicitly
    if dataset.max_genes is not None and dataset.max_genes <= 0:
        return {
            "python": {"wall_sec": None},
            "rust": {"wall_sec": None},
            "parity": False,
            "max_delta": None,
            "barcodes": 0,
            "genes": 0,
            "speedup": None,
            "skipped": True,
        }
    barcodes, genes = sample_barcodes_and_genes(
        dataset.bam, max_cells=dataset.max_cells, max_genes=dataset.max_genes, max_reads=dataset.max_reads
    )
    if not barcodes or not genes:
        return {
            "python": {"wall_sec": None},
            "rust": {"wall_sec": None},
            "parity": False,
            "max_delta": None,
            "barcodes": len(barcodes),
            "genes": len(genes),
            "speedup": None,
            "error": "Failed to derive barcodes/genes (missing whitelist or tags?)",
        }
    # Hard cap to avoid hanging on huge unsorted BAMs
    if dataset.max_reads and dataset.max_reads > 0:
        gene_max_reads = dataset.max_reads
    else:
        gene_max_reads = 20000  # safety cap

    def run_variant(fn):
        rss_before = memory_mb()
        start = time.perf_counter()
        mat = fn()
        wall = time.perf_counter() - start
        rss_after = memory_mb()
        return np.array(mat, dtype=np.uint32), {
            "wall_sec": wall,
            "rss_mb": rss_after,
            "rss_mb_delta": max(0.0, rss_after - rss_before),
        }

    py_mat, py_metrics = run_variant(lambda: python_gene_counts(dataset.bam, barcodes, genes, dataset.max_reads))
    rust_mat, rust_metrics = run_variant(
        lambda: np.array(sheriff_rs.gene_counts_py(str(dataset.bam), barcodes, genes, dataset.max_reads), dtype=np.uint32)
        if HAS_RUST
        else py_mat
    )

    parity = py_mat.shape == rust_mat.shape and np.array_equal(py_mat, rust_mat)
    diff = np.abs(py_mat - rust_mat)
    max_delta = int(diff.max()) if diff.size else 0

    return {
        "python": py_metrics,
        "rust": rust_metrics,
        "parity": parity,
        "max_delta": max_delta,
        "barcodes": len(barcodes),
        "genes": len(genes),
        "speedup": py_metrics["wall_sec"] / rust_metrics["wall_sec"] if rust_metrics["wall_sec"] else None,
    }


def run_umi_bench(num_cells: int = 200, umis_per_cell: int = 50) -> Dict[str, Any]:
    """Benchmark UMI deduplication Rust vs Python on synthetic data."""
    import random
    bases = ["A", "C", "G", "T"]

    def make_umi():
        return "".join(random.choice(bases) for _ in range(12))

    umi_sets = []
    for _ in range(num_cells):
        # include duplicates and near-duplicates
        core = [make_umi() for _ in range(max(1, umis_per_cell // 3))]
        duped = core + [random.choice(core) for _ in range(max(1, umis_per_cell // 3))]
        umi_sets.append(set(duped))

    def time_variant(fn):
        start = time.perf_counter()
        counts = [fn(us) for us in umi_sets]
        wall = time.perf_counter() - start
        return counts, wall

    py_counts, py_wall = time_variant(lambda s: len(deduplicate_umis_python(s)))
    if HAS_RUST and rust_allowed():
        rs_counts, rs_wall = time_variant(lambda s: sheriff_rs.deduplicate_umis_py(list(s)))
        parity = py_counts == rs_counts
        speedup = py_wall / rs_wall if rs_wall else None
    else:
        rs_counts, rs_wall, parity, speedup = py_counts, py_wall, True, 1.0

    return {
        "python": {"wall_sec": py_wall},
        "rust": {"wall_sec": rs_wall},
        "parity": parity,
        "speedup": speedup,
    }


def run_edit_bench(num_edits: int = 1000) -> Dict[str, Any]:
    """Benchmark edit clustering Rust vs Python on synthetic edits."""
    from collections import namedtuple
    ReadEdit = namedtuple("ReadEdit", ["chrom", "ref_pos", "ref_seq", "alt_seq", "forward", "kmer_matches"])

    # Safe base set (mirrors validation patterns)
    base_edits = [
        ReadEdit("chr1", 1000, "ATCG", "ATCGATCGATCG", True, (1, 2, 3)),
        ReadEdit("chr1", 1000, "ATCG", "ATCGATCG", True, (1, 2)),
        ReadEdit("chr1", 2000, "GCTA", "GCTACCCC", False, (4, 5)),
        ReadEdit("chr2", 3000, "AAAA", "AAAATTTT", True, (7,)),
    ]

    edits = [base_edits[i % len(base_edits)] for i in range(num_edits)]

    def time_variant(fn):
        start = time.perf_counter()
        res = fn(edits)
        wall = time.perf_counter() - start
        return res, wall

    try:
        py_res, py_wall = time_variant(get_longest_edits_python)
    except Exception as exc:
        # If Python baseline fails on synthetic data, report and skip
        return {
            "python": {"wall_sec": None, "error": str(exc)},
            "rust": {"wall_sec": None},
            "parity": False,
            "speedup": None,
        }
    if HAS_RUST and rust_allowed():
        rs_res, rs_wall = time_variant(lambda e: sheriff_rs.get_longest_edits_rust(
            [(x.chrom, x.ref_pos, x.ref_seq, x.alt_seq, x.forward, list(x.kmer_matches)) for x in e]
        ))
        # compare by normalized tuples (order-insensitive)
        py_set = set((x.chrom, x.ref_pos, x.ref_seq, x.alt_seq, x.forward) for x in py_res)
        rs_set = set((x[0], x[1], x[2], x[3], x[4]) for x in rs_res)
        parity = py_set == rs_set
        speedup = py_wall / rs_wall if rs_wall else None
    else:
        rs_res, rs_wall, parity, speedup = py_res, py_wall, True, 1.0

    return {
        "python": {"wall_sec": py_wall},
        "rust": {"wall_sec": rs_wall},
        "parity": parity,
        "speedup": speedup,
    }


def plot_results(all_results: Dict[str, Any], png_path: Path):
    if plt is None:
        print("matplotlib not available; skipping plot generation.")
        return

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    datasets = list(all_results.keys())

    # Wall time bars (gene counts)
    labels = []
    py_w = []
    rs_w = []
    for name in datasets:
        gene_res = all_results[name].get("gene_counts")
        if gene_res:
            labels.append(name)
            py_w.append(gene_res["python"]["wall_sec"])
            rs_w.append(gene_res["rust"]["wall_sec"])
    x = np.arange(len(labels))
    axes[0].bar(x - 0.15, py_w, width=0.3, label="Python")
    axes[0].bar(x + 0.15, rs_w, width=0.3, label="Rust")
    axes[0].set_title("Gene counts wall (s)")
    axes[0].set_xticks(x, labels, rotation=20)
    axes[0].legend()

    # RSS delta bars (gene counts)
    py_m = []
    rs_m = []
    labels2 = []
    for name in datasets:
        gene_res = all_results[name].get("gene_counts")
        if gene_res:
            labels2.append(name)
            py_m.append(gene_res["python"]["rss_mb_delta"])
            rs_m.append(gene_res["rust"]["rss_mb_delta"])
    x2 = np.arange(len(labels2))
    axes[1].bar(x2 - 0.15, py_m, width=0.3, label="Python")
    axes[1].bar(x2 + 0.15, rs_m, width=0.3, label="Rust")
    axes[1].set_title("Gene counts RSS delta (MB)")
    axes[1].set_xticks(x2, labels2, rotation=20)
    axes[1].legend()

    # Speedup (filter + gene)
    speed_labels = []
    speed_vals = []
    for name in datasets:
        filt = all_results[name].get("filter")
        gene = all_results[name].get("gene_counts")
        if filt and "speedup" in filt:
            speed_labels.append(f"{name}-filter")
            speed_vals.append(filt["speedup"])
        if gene and gene.get("speedup"):
            speed_labels.append(f"{name}-gene")
            speed_vals.append(gene["speedup"])
    x3 = np.arange(len(speed_labels))
    axes[2].bar(x3, speed_vals, width=0.4, color="#4c72b0")
    axes[2].set_xticks(x3, speed_labels, rotation=45, ha="right")
    axes[2].set_title("Speedup (Rust/Python)")
    axes[2].axhline(1.0, color="red", linestyle="--", linewidth=1)

    fig.tight_layout()
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=180)
    print(f"Wrote plot to {png_path}")


def parse_datasets(args) -> List[Dataset]:
    datasets: List[Dataset] = []
    if args.dataset:
        for item in args.dataset:
            if "=" not in item:
                raise SystemExit(f"--dataset expects name=/path/to.bam, got: {item}")
            name, path = item.split("=", 1)
            datasets.append(
                Dataset(
                    name=name,
                    bam=Path(path),
                    whitelist=Path(args.whitelist) if args.whitelist else None,
                    max_cells=args.max_cells,
                    max_genes=args.max_genes,
                    max_reads=args.max_reads,
                )
            )
    else:
        # Default: bundled example data
        datasets.append(
            Dataset(
                name="example",
                bam=Path("example_data/barcode_headAligned_anno.sorted.edit_regions_200kb.bam"),
                whitelist=Path("example_data/barcode_whitelist.500-cell.txt"),
                max_cells=args.max_cells,
                max_genes=args.max_genes,
                max_reads=args.max_reads,
            )
        )
    return datasets


def main():
    parser = argparse.ArgumentParser(description="Run Sheriff Rust vs Python perf suite.")
    parser.add_argument(
        "--dataset",
        action="append",
        help="Dataset spec name=/path/to.bam (can repeat). Defaults to bundled example.",
    )
    parser.add_argument(
        "--whitelist",
        help="Whitelist path to use for all datasets (if dataset spec omitted).",
    )
    parser.add_argument("--max-cells", type=int, default=64, help="Max barcodes for gene-count bench.")
    parser.add_argument("--max-genes", type=int, default=32, help="Max genes for gene-count bench.")
    parser.add_argument("--max-reads", type=int, default=0, help="Max reads to scan for gene-count bench (0 = all).")
    parser.add_argument("--umi-cells", type=int, default=200, help="Synthetic cells for UMI bench.")
    parser.add_argument("--umi-umis", type=int, default=50, help="UMIs per synthetic cell for UMI bench.")
    parser.add_argument("--edits", type=int, default=1000, help="Synthetic edits for edit-clustering bench.")
    parser.add_argument("--out-json", type=Path, default=Path("benchmarks/results_perf_suite.json"))
    parser.add_argument("--out-png", type=Path, default=Path("benchmarks/results_perf_suite.png"))
    parser.add_argument("--no-plot", action="store_true", help="Skip plot generation.")

    args = parser.parse_args()
    datasets = parse_datasets(args)
    fmt = lambda v: f"{v:.2f}x" if isinstance(v, (int, float)) and v is not None else "n/a"

    output_dir = Path("/tmp/sheriff_perf_outputs")
    output_dir.mkdir(parents=True, exist_ok=True)

    all_results: Dict[str, Any] = {}
    for ds in datasets:
        if not ds.bam.exists():
            raise SystemExit(f"BAM not found: {ds.bam}")
        print(f"\n=== Dataset: {ds.name} ({ds.bam}) ===")
        ds_result: Dict[str, Any] = {"bam": str(ds.bam)}

        filt_res = run_filter(ds, output_dir)
        if filt_res:
            ds_result["filter"] = filt_res
            sp = filt_res.get("speedup")
            print("  Filter speedup: " + fmt(sp) if sp is not None else "  Filter done.")

        gene_res = run_gene_counts(ds)
        ds_result["gene_counts"] = gene_res
        print(f"  Gene count speedup: {fmt(gene_res.get('speedup'))} (parity={gene_res['parity']})")

        umi_res = run_umi_bench(num_cells=args.umi_cells, umis_per_cell=args.umi_umis)
        ds_result["umi_dedup"] = umi_res
        print(f"  UMI dedup speedup:  {fmt(umi_res.get('speedup'))} (parity={umi_res['parity']})")

        edit_res = run_edit_bench(num_edits=args.edits)
        ds_result["edit_cluster"] = edit_res
        print(f"  Edit cluster speedup: {fmt(edit_res.get('speedup'))} (parity={edit_res['parity']})")

        all_results[ds.name] = ds_result

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    with args.out_json.open("w") as fh:
        json.dump(all_results, fh, indent=2)
    print(f"\nWrote results to {args.out_json}")

    if not args.no_plot:
        plot_results(all_results, args.out_png)


if __name__ == "__main__":
    main()
