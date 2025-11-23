#!/usr/bin/env python3
"""
Create publication-ready plots from perf suite outputs.

Reads:
  - benchmarks/results_perf_suite.json (from run_perf_suite.py)

Outputs:
  - benchmarks/results_perf_suite_wall.(pdf|png)     # wall seconds (Python vs Rust)
  - benchmarks/results_perf_suite_rss.(pdf|png)      # RSS delta (MB, Python vs Rust)
  - benchmarks/results_perf_suite_speed.(pdf|png)    # speedup (fold, Rust/Python)
  - benchmarks/results_perf_suite_panels.(pdf|png)   # combined 3-panel figure
"""

import json
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import seaborn as sns

# Visual style tuned for manuscripts/presentations
sns.set_theme(style="whitegrid", font="DejaVu Serif", context="talk")

PY_COLOR = "#a6bddb"   # muted blue
RS_COLOR = "#045a8d"   # deep blue
SP_COLOR = "#41ab5d"   # green for speedup

OPS = [
    ("filter", "Filter"),
    ("gene_counts", "Gene counts"),
    ("umi_dedup", "UMI dedup"),
    ("edit_cluster", "Edit cluster"),
]


def load_suite(path: Path) -> Dict:
    with path.open() as fh:
        return json.load(fh)


def fmt(val: float) -> str:
    if val is None:
        return ""
    if abs(val) >= 100:
        return f"{val:.0f}"
    if abs(val) >= 10:
        return f"{val:.1f}"
    return f"{val:.2f}"


def extract_entries(data: Dict) -> List[Dict]:
    entries: List[Dict] = []
    for ds_name, ds in data.items():
        for key, label in OPS:
            comp = ds.get(key)
            if not comp:
                continue
            entries.append(
                {
                    "dataset": ds_name,
                    "op": key,
                    "label": label,
                    "python_wall": comp.get("python", {}).get("wall_sec"),
                    "rust_wall": comp.get("rust", {}).get("wall_sec"),
                    "python_rss": comp.get("python", {}).get("rss_mb_delta"),
                    "rust_rss": comp.get("rust", {}).get("rss_mb_delta"),
                    "speedup": comp.get("speedup"),
                }
            )
    return entries


def make_labels(entries: List[Dict]) -> List[str]:
    datasets = {e["dataset"] for e in entries}
    multi = len(datasets) > 1
    labels = []
    for e in entries:
        base = e["label"]
        labels.append(f"{e['dataset']} · {base}" if multi else base)
    return labels


def bar_pair(ax, labels, py_vals, rs_vals, title, ylabel):
    x = range(len(labels))
    ax.bar([i - 0.22 for i in x], py_vals, width=0.44, label="Python", color=PY_COLOR, edgecolor="#34495e")
    ax.bar([i + 0.22 for i in x], rs_vals, width=0.44, label="Rust", color=RS_COLOR, edgecolor="#1b2833")
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_title(title, pad=8)
    ax.set_ylabel(ylabel)
    ax.legend(frameon=False, loc="upper right")
    # Annotate bars
    for i, (p, r) in enumerate(zip(py_vals, rs_vals)):
        ax.text(i - 0.22, p, fmt(p), ha="center", va="bottom", fontsize=9, rotation=90)
        ax.text(i + 0.22, r, fmt(r), ha="center", va="bottom", fontsize=9, rotation=90)


def make_wall(entries: List[Dict], out_pdf: Path, out_png: Path):
    labels = make_labels(entries)
    py_vals = [e["python_wall"] for e in entries if e["python_wall"] is not None and e["rust_wall"] is not None]
    rs_vals = [e["rust_wall"] for e in entries if e["python_wall"] is not None and e["rust_wall"] is not None]
    labels = [lbl for e, lbl in zip(entries, labels) if e["python_wall"] is not None and e["rust_wall"] is not None]

    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    bar_pair(ax, labels, py_vals, rs_vals, "Wall time", "Seconds")
    fig.tight_layout()
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=300)


def make_rss(entries: List[Dict], out_pdf: Path, out_png: Path):
    mask = [e for e in entries if e["python_rss"] is not None and e["rust_rss"] is not None]
    if not mask:
        return
    labels = make_labels(mask)
    py_vals = [e["python_rss"] for e in mask]
    rs_vals = [e["rust_rss"] for e in mask]

    fig, ax = plt.subplots(figsize=(8.5, 4.0))
    bar_pair(ax, labels, py_vals, rs_vals, "RSS delta", "MB")
    fig.tight_layout()
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=300)


def make_speed(entries: List[Dict], out_pdf: Path, out_png: Path):
    mask = [e for e in entries if e["speedup"] is not None]
    labels = make_labels(mask)
    vals = [e["speedup"] for e in mask]

    fig, ax = plt.subplots(figsize=(8.5, 4.25))
    ax.bar(range(len(vals)), vals, color=SP_COLOR, edgecolor="#1b5e20", width=0.55)
    ax.axhline(1.0, color="#d95f02", linestyle="--", linewidth=1.2, label="Parity (1x)")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=25, ha="right")
    ax.set_title("Speedup (Rust / Python)", pad=8)
    ax.set_ylabel("Fold")
    ax.set_yscale("log")
    # Annotate speedups
    for i, v in enumerate(vals):
        ax.text(i, v, fmt(v), ha="center", va="bottom", fontsize=9, rotation=90)
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=300)


def make_panels(entries: List[Dict], out_pdf: Path, out_png: Path):
    labels = make_labels(entries)
    wall_mask = [e for e in entries if e["python_wall"] is not None and e["rust_wall"] is not None]
    rss_mask = [e for e in entries if e["python_rss"] is not None and e["rust_rss"] is not None]
    sp_mask = [e for e in entries if e["speedup"] is not None]

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

    # Wall
    if wall_mask:
        wl = make_labels(wall_mask)
        bar_pair(axes[0], wl, [e["python_wall"] for e in wall_mask], [e["rust_wall"] for e in wall_mask], "Wall time", "Seconds")
    else:
        axes[0].set_visible(False)

    # RSS
    if rss_mask:
        rl = make_labels(rss_mask)
        bar_pair(axes[1], rl, [e["python_rss"] for e in rss_mask], [e["rust_rss"] for e in rss_mask], "RSS delta", "MB")
    else:
        axes[1].set_visible(False)

    # Speedup
    if sp_mask:
        sl = make_labels(sp_mask)
        vals = [e["speedup"] for e in sp_mask]
        axes[2].bar(range(len(vals)), vals, color=SP_COLOR, edgecolor="#1b5e20", width=0.55)
        axes[2].axhline(1.0, color="#d95f02", linestyle="--", linewidth=1.2, label="Parity (1x)")
        axes[2].set_xticks(range(len(sl)))
        axes[2].set_xticklabels(sl, rotation=25, ha="right")
        axes[2].set_title("Speedup (Rust / Python)", pad=8)
        axes[2].set_ylabel("Fold")
        axes[2].set_yscale("log")
        for i, v in enumerate(vals):
            axes[2].text(i, v, fmt(v), ha="center", va="bottom", fontsize=9, rotation=90)
        axes[2].legend(frameon=False, loc="upper left")
    else:
        axes[2].set_visible(False)

    fig.tight_layout()
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=300)


def main():
    suite_path = Path("benchmarks/results_perf_suite.json")
    out_wall_pdf = Path("benchmarks/results_perf_suite_wall.pdf")
    out_wall_png = Path("benchmarks/results_perf_suite_wall.png")
    out_rss_pdf = Path("benchmarks/results_perf_suite_rss.pdf")
    out_rss_png = Path("benchmarks/results_perf_suite_rss.png")
    out_speed_pdf = Path("benchmarks/results_perf_suite_speed.pdf")
    out_speed_png = Path("benchmarks/results_perf_suite_speed.png")
    out_panels_pdf = Path("benchmarks/results_perf_suite_panels.pdf")
    out_panels_png = Path("benchmarks/results_perf_suite_panels.png")

    if not suite_path.exists():
        raise SystemExit(f"Missing {suite_path}; run run_perf_suite.py first.")

    data = load_suite(suite_path)
    entries = extract_entries(data)

    out_wall_pdf.parent.mkdir(parents=True, exist_ok=True)
    make_wall(entries, out_wall_pdf, out_wall_png)
    make_rss(entries, out_rss_pdf, out_rss_png)
    make_speed(entries, out_speed_pdf, out_speed_png)
    make_panels(entries, out_panels_pdf, out_panels_png)


if __name__ == "__main__":
    main()
