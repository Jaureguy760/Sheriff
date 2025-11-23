# Benchmarks

Quick pointers to run and plot the Python vs Rust perf suite.

## Superb BAM slice (recommended comparison)

```bash
python benchmarks/run_perf_suite.py \
  --dataset superb=benchmarks/large_superb.bam \
  --whitelist benchmarks/large_superb_whitelist.txt \
  --max-cells 30 \
  --max-genes 200 \
  --max-reads 100000 \
  --log benchmarks/run_perf_suite_superb.log

python benchmarks/plot_perf.py
```

Outputs land in `benchmarks/`:
- `results_perf_suite.json` (timings, RSS, speedups, parity)
- `results_perf_suite_wall|rss|speed|panels.{pdf,png}` (publication-ready plots)
- `run_perf_suite_superb.log` (run log)

Notes:
- The superb BAM is a local symlink and is not tracked in git; adjust `--max-reads` if you want a larger slice.
- Speedup plot uses log scale with a 1x parity line.

## Small example data

For a fast sanity check on the bundled 30 MB BAM:

```bash
python benchmarks/run_perf_suite.py --log benchmarks/run_perf_suite_example.log
python benchmarks/plot_perf.py
```

## More details

See `benchmarks/PERF_SUITE.md` for full context, styling notes, and output descriptions.***
