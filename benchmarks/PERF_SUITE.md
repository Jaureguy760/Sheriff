# Sheriff Perf Suite (Python vs Rust)

How to reproduce the end-to-end benchmark + plots on the bundled example data or the larger superbseq BAM.

## Quick run (large superbseq BAM, capped reads)

```bash
# From repo root, with Sheriff installed (pip install -e .) and rust wheel built/installed
python benchmarks/run_perf_suite.py \
  --dataset superb=benchmarks/large_superb.bam \
  --whitelist benchmarks/large_superb_whitelist.txt \
  --max-cells 30 \
  --max-genes 200 \
  --max-reads 100000 \
  --log benchmarks/run_perf_suite_superb.log

# Plot publication-ready figures
python benchmarks/plot_perf.py
```

Outputs:
- `benchmarks/results_perf_suite.json` (timings, RSS, speedups, parity)
- `benchmarks/results_perf_suite_wall.(pdf|png)`   – wall time Python vs Rust
- `benchmarks/results_perf_suite_rss.(pdf|png)`    – RSS delta Python vs Rust
- `benchmarks/results_perf_suite_speed.(pdf|png)`  – fold speedup
- `benchmarks/results_perf_suite_panels.(pdf|png)` – 3-panel composite figure
- `benchmarks/run_perf_suite_superb.log`           – log of the run

Notes:
- The superb BAM is a symlink to `/iblm/netapp/data4/mlorenzini/.../barcode_headAligned_anno.bam` and is *not* tracked in git.
- The BAM is unsorted/unindexed; use `--max-reads` (e.g., 100k) to bound runtime when comparing Python vs Rust fairly.
- `--max-cells` and `--max-genes` limit the sampled barcodes/genes; adjust to taste if you want a larger slice.

## Small example data

For a quick sanity pass on the bundled 30 MB BAM in `example_data/`, run without the superb arguments:

```bash
python benchmarks/run_perf_suite.py --log benchmarks/run_perf_suite_example.log
python benchmarks/plot_perf.py
```

## Plot styling

`plot_perf.py` produces separate panels (wall, RSS, speedup) and a combined 3-panel figure with:
- Separate colors for Python (muted blue) and Rust (deep blue); speedup in green.
- Log-scaled speedup panel (useful when edit clustering dwarfs other ops).
- Per-bar annotations with rounded values and a parity (1x) reference line on speedup.

## Reusing the data

`results_perf_suite.json` is stable across runs; delete it if you want a fresh capture. Generated plots and the large BAM symlink are ignored by git to keep the tree clean.
