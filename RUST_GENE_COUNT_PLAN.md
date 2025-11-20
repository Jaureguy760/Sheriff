Rust Gene Counting Optimization Plan
====================================

Objective
---------
Reduce bottlenecks in gene UMI counting and provide a fair Rust vs Python benchmark path with clear toggles, validation, and memory/runtime gains.

Current State
-------------
- Rust accelerations: BAM filtering, k-mer matching, edit clustering, UMI dedup/per-cell counts, *initial* gene-count binding.
- Toggle: `--no-rust` / `SHERIFF_DISABLE_RUST` wired into k-mers and gene counting; Python fallback remains.
- Contig filtering: gene counting limited to contigs with mapped reads.
- Gap: Rust gene counting currently inflates RSS (~2.9 GB on chr19) due to dense per-gene/per-cell String hash sets.

Plan (Steps)
------------
1) Refactor Rust gene counting for low memory and speed
   - Stream per contig; avoid per-gene/per-cell String sets.
   - Use byte slices and numeric IDs (gene/barcode indexing) to avoid allocs.
   - Produce a compact matrix (Vec<Vec<u32>> or a flat Vec<u32>) and transpose in Python.
   - Preserve Python fallback; parity required.

2) Propagate Rust toggle everywhere
   - Ensure `use_rust`/`SHERIFF_DISABLE_RUST` gates k-mer matching, edit clustering, UMI dedup/counting, gene counting.
   - Add tests to assert Python vs Rust outputs match when toggled.

3) Benchmark harness
   - Add `benchmarks/benchmark_pipeline.py` to time full runs with/without Rust on:
     * chr19 test bundle (sanity)
     * mid-size BAM (e.g., 300–500k reads)
   - Collect wall time, max RSS (`/usr/bin/time -v`), and output equality checks (counts/edit files).

4) Validation & tests
   - Extend `validate_rust_correctness.py` to include gene count parity (Python vs Rust).
   - Add regression test that runs a small dataset with `--no-rust` vs Rust and asserts identical outputs.
   - Keep Python fallback intact; CI should pass without Rust installed.

5) Profile remaining hotspots
   - Identify heavy pandas/dense DataFrame steps; consider sparse/streaming or Rust helpers if needed.
   - Optional: skip empty chunks based on BAM index stats to further reduce work.

Sanity Checks
-------------
- Functional: Python vs Rust outputs equal (counts/edit files).
- Performance: Rust gene counting RSS lower than current approach; improved wall time on mid-size BAM.
- Toggle: `--no-rust` and `SHERIFF_DISABLE_RUST=1` reliably force Python path without errors.

Deliverables
------------
- Refactored Rust gene counting (low-RSS design) with PyO3 binding.
- Toggle propagated across accelerated paths.
- Benchmark script and documented results (chr19 + mid-size BAM).
- Updated docs/tests covering toggles, parity, and expected speedups.
