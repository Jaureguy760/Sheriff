# Sheriff Deep Integration Engineering Plan - True 10/10

## Executive Summary

**Goal**: Instrument count_t7.py for real-time progress, state checkpoints, and live metrics
**Approach**: Minimal invasive instrumentation at strategic points
**Time**: 4-6 hours (optimized plan)
**Risk**: Low (backward compatible, optional features)
**Result**: True 10/10 CLI with production-ready resume/progress/metrics

---

## Phase 1: Architecture Analysis (30 min)

### 1.1 count_t7.py Structure Map

```
run_count_t7() [Line 548-1373]
├── Setup (Lines 574-586)
│   └── Load barcodes, create output dir
│
├── Step 1: Barcoded Edits (Lines 588-606) ← INSTRUMENT
│   └── get_barcoded_edits() [Lines 307-463]
│       ├── BAM filtering (chromosome iteration)
│       └── K-mer matching (read iteration)
│
├── Step 2: Edit Site Clustering (Lines 608-788) ← INSTRUMENT
│   ├── Define canonical edit sites
│   ├── Cluster edits by distance
│   └── Filter by criteria (bidirectional, min cells)
│
├── Step 3: Non-barcoded Edits (Lines 790-820) ← INSTRUMENT
│   └── get_nonbarcoded_edits() [Lines 465-546]
│       └── Assign non-barcoded reads to edit sites
│
├── Step 4: Gene Counting (Lines 822-1200) ← INSTRUMENT
│   ├── Build gene models
│   ├── UMI counting (parallel by chromosome)
│   └── Allele-specific counting
│
└── Step 5: Output Generation (Lines 1202-1373) ← INSTRUMENT
    ├── Write edit_site_info.txt
    ├── Write UMI matrices
    └── Write gene count matrices
```

### 1.2 Integration Points

**Checkpoint Opportunities** (6 total):
1. **After validation**: Inputs checked, parameters set
2. **After barcoded edits**: get_barcoded_edits() complete
3. **After edit clustering**: Canonical sites defined
4. **After non-barcoded**: All edit reads assigned
5. **After UMI counting**: Gene counts computed
6. **After outputs**: All files written

**Progress Tracking** (5 major tasks):
1. **Barcoded edits**: Track chromosomes processed
2. **Edit clustering**: Track edits processed
3. **Non-barcoded assignment**: Track reads processed
4. **UMI counting**: Track genes/chunks processed
5. **Output generation**: Track files written

**Metrics Collection** (throughout):
- Input: Total reads, BAM size, barcode count
- Filtering: Reads processed, barcoded reads found
- Edits: Candidate sites, canonical sites, cells with edits
- Counting: Genes quantified, UMIs counted
- Performance: Time per step

---

## Phase 2: Function Signature Extension (15 min)

### 2.1 Add Optional Parameters

```python
def run_count_t7(
    bam_file,
    ref_file,
    barcode_file,
    gtf_file=None,
    # ... existing 20+ parameters ...

    # NEW: Priority 3 Integration (all optional, default=None)
    progress_tracker=None,      # PipelineProgress instance
    checkpoint_manager=None,     # CheckpointManager instance
    results_collector=None,      # PipelineResults instance
    enable_instrumentation=True, # Master switch for all features
):
    """
    Run Sheriff CRISPR edit calling pipeline.

    New Parameters (Priority 3):
        progress_tracker: Optional PipelineProgress for real-time updates
        checkpoint_manager: Optional CheckpointManager for resume capability
        results_collector: Optional PipelineResults for metrics collection
        enable_instrumentation: If False, skip all Priority 3 features (for testing)
    """
```

**Backward Compatibility**:
- All new parameters default to None
- If None, features are skipped (zero overhead)
- Existing code unchanged

### 2.2 Initialize Instrumentation

```python
# At start of run_count_t7() [after line 572]

# Initialize instrumentation helpers
_has_progress = progress_tracker is not None and enable_instrumentation
_has_checkpoint = checkpoint_manager is not None and enable_instrumentation
_has_results = results_collector is not None and enable_instrumentation

# Helper function for safe updates
def _update_progress(task, advance=1, **kwargs):
    if _has_progress:
        progress_tracker.update(task, advance=advance, **kwargs)

def _save_checkpoint(step_name, outputs=None, metrics=None):
    if _has_checkpoint:
        checkpoint_manager.update_outputs(**(outputs or {}))
        checkpoint_manager.update_metrics(**(metrics or {}))
        checkpoint_manager.save()

def _record_metric(key, value):
    if _has_results:
        results_collector.set_metric(key, value)

# Start timing
import time
_pipeline_start = time.time()
```

---

## Phase 3: Step-by-Step Instrumentation (2 hours)

### 3.1 Step 1: Barcoded Edits [Lines 588-606]

**Current Code**:
```python
start_bc = timeit.default_timer()
print("Counting barcoded edits...", file=sys.stdout, flush=True) if verbosity >= 1 else None

edit_counts, edit_bc_cell_umis, edit_reads, chr_edit_loc_indices, chr_edits_order_added = get_barcoded_edits(
    bam_file,
    cell_barcodes,
    ref_file,
    k, t7_barcode,
    output_kmer_hash=False,
    blacklist_seqs=blacklist_seqs,
    print_freq=print_freq, verbosity=verbosity,
)

print(f"Processed barcoded edits in {(timeit.default_timer()-start_bc)/60:.3f} minutes\n",
      file=sys.stdout, flush=True) if verbosity>= 1 else None
```

**Instrumented Code**:
```python
# Checkpoint: Check if can skip this step
if _has_checkpoint and checkpoint_manager.can_skip_step("barcoded_edits"):
    logger.info("Resuming: Loading barcoded edits from checkpoint")
    # Load from checkpoint
    edit_counts = checkpoint_manager.current_checkpoint.outputs.edit_counts
    # ... load other outputs
else:
    start_bc = time.time()

    # Progress: Add task for this step
    if _has_progress:
        # Count chromosomes for progress tracking
        import pysam
        bam = pysam.AlignmentFile(bam_file, "rb")
        n_chromosomes = len(bam.references)
        bam.close()
        progress_tracker.add_task("barcoded_edits", "Finding barcoded edits", total=n_chromosomes)

    print("Counting barcoded edits...", file=sys.stdout, flush=True) if verbosity >= 1 else None

    # Call function (need to pass progress through)
    edit_counts, edit_bc_cell_umis, edit_reads, chr_edit_loc_indices, chr_edits_order_added = get_barcoded_edits(
        bam_file,
        cell_barcodes,
        ref_file,
        k, t7_barcode,
        output_kmer_hash=False,
        blacklist_seqs=blacklist_seqs,
        print_freq=print_freq,
        verbosity=verbosity,
        progress_callback=lambda: _update_progress("barcoded_edits", advance=1) if _has_progress else None,
    )

    elapsed_bc = time.time() - start_bc

    # Metrics: Record counts
    _record_metric("barcoded_reads", len(edit_reads))
    _record_metric("candidate_edits", len(edit_counts))

    # Performance: Record timing
    if _has_results:
        results_collector.set_performance("barcoded_edits", elapsed_bc)

    # Checkpoint: Save after this step
    _save_checkpoint(
        "barcoded_edits",
        outputs={
            "edit_counts": edit_counts,
            "edit_bc_cell_umis": edit_bc_cell_umis,
            # ... other outputs (or save to temp files)
        },
        metrics={
            "barcoded_reads": len(edit_reads),
            "runtime_seconds": elapsed_bc,
        }
    )

    print(f"Processed barcoded edits in {elapsed_bc/60:.3f} minutes\n",
          file=sys.stdout, flush=True) if verbosity>= 1 else None
```

### 3.2 Modify get_barcoded_edits() for Progress [Lines 307-463]

**Strategy**: Add optional progress_callback parameter

**Minimal Change**:
```python
def get_barcoded_edits(
    bam_file, cell_barcodes, ref_file,
    k, t7_barcode,
    output_kmer_hash=False,
    blacklist_seqs=None,
    print_freq=1000000,
    verbosity=1,
    progress_callback=None,  # NEW: Optional progress callback
):
    # ... existing code ...

    # In chromosome iteration loop [around line 340]
    for chrom_index, chrom in enumerate(chromosomes):
        # ... process chromosome ...

        # NEW: Call progress callback after each chromosome
        if progress_callback:
            progress_callback()

        # ... rest of chromosome processing ...
```

**Impact**: Minimal - just one function call per chromosome (low overhead)

### 3.3 Step 2: Edit Clustering [Lines 608-788]

```python
# After edit clustering completes [around line 788]

start_cluster = time.time()

# Progress: Track clustering
if _has_progress:
    progress_tracker.add_task("clustering", "Clustering edit sites", total=len(edit_order))

# Existing clustering code with progress updates
for orderi, editi in enumerate(edit_order):
    # ... existing clustering logic ...

    # NEW: Update progress every 100 edits (not every one - too slow)
    if _has_progress and orderi % 100 == 0:
        progress_tracker.update("clustering", completed=orderi)

elapsed_cluster = time.time() - start_cluster

# Metrics
_record_metric("canonical_edits", len(canonical_edit_sites))
_record_metric("cells_with_edits", sum(canonical_edit_cell_counts))

# Performance
if _has_results:
    results_collector.set_performance("edit_clustering", elapsed_cluster)

# Checkpoint
_save_checkpoint(
    "edit_clustering",
    outputs={"canonical_sites": canonical_edit_sites},
    metrics={"canonical_edits": len(canonical_edit_sites)}
)
```

### 3.4 Step 3: Non-barcoded Edits [Lines 790-820]

```python
# Similar pattern to Step 1
if _has_checkpoint and checkpoint_manager.can_skip_step("nonbarcoded_edits"):
    # Resume: Skip this step
    logger.info("Resuming: Skipping non-barcoded edits (already done)")
else:
    start_nonbc = time.time()

    if _has_progress:
        progress_tracker.add_task("nonbarcoded", "Processing non-barcoded reads", total=None)

    # Call existing function
    nonbc_reads_data = get_nonbarcoded_edits(...)

    elapsed_nonbc = time.time() - start_nonbc

    _record_metric("nonbarcoded_reads", len(nonbc_reads_data))
    if _has_results:
        results_collector.set_performance("nonbarcoded_edits", elapsed_nonbc)

    _save_checkpoint("nonbarcoded_edits",
                     metrics={"nonbarcoded_reads": len(nonbc_reads_data)})
```

### 3.5 Step 4: UMI Counting [Lines 822-1200]

**Complex section with parallel processing**

```python
# UMI counting has chromosome-parallel processing
# Add progress per gene or per chunk

start_umi = time.time()

if _has_progress:
    # Estimate total work (number of genes)
    n_genes = len(gene_models) if gene_models else 1000
    progress_tracker.add_task("umi_counting", "Counting UMIs", total=n_genes)

# Existing UMI counting code
# Add periodic updates in the counting loop
for gene_idx, gene in enumerate(genes):
    # ... existing counting logic ...

    # Update progress every 100 genes
    if _has_progress and gene_idx % 100 == 0:
        progress_tracker.update("umi_counting", completed=gene_idx)

elapsed_umi = time.time() - start_umi

# Metrics
_record_metric("genes_quantified", len(gene_umi_counts))
_record_metric("umis_counted", total_umi_count)

if _has_results:
    results_collector.set_performance("umi_counting", elapsed_umi)

_save_checkpoint("umi_counting",
                 outputs={"umi_counts_file": umi_output_path},
                 metrics={"genes_quantified": len(gene_umi_counts)})
```

### 3.6 Step 5: Output Generation [Lines 1202-1373]

```python
start_output = time.time()

if _has_progress:
    progress_tracker.add_task("outputs", "Writing output files", total=3)

# Write edit site info
# ... existing code ...
_update_progress("outputs", advance=1)

# Write UMI matrix
# ... existing code ...
_update_progress("outputs", advance=1)

# Write gene counts
# ... existing code ...
_update_progress("outputs", advance=1)

elapsed_output = time.time() - start_output

# Record output files
if _has_results:
    results_collector.set_output("Edit Sites", edit_site_file)
    results_collector.set_output("UMI Counts", umi_file)
    results_collector.set_output("Gene Counts", gene_file)
    results_collector.set_performance("output_generation", elapsed_output)

# Final checkpoint
_save_checkpoint("outputs", status="completed")
```

---

## Phase 4: CLI Integration (30 min)

### 4.1 Update cli/run.py

```python
# In sheriff/cli/run.py run() function

# After config loading and validation

# Initialize instrumentation (if enabled)
progress_tracker = None
checkpoint_manager = None
results_collector = None

if verbosity >= 1:  # Only if not silent
    from ..progress import PipelineProgress
    progress_tracker = PipelineProgress(verbosity=verbosity)
    progress_tracker.__enter__()  # Start progress

if enable_checkpoints or resume:
    from ..checkpoint import CheckpointManager
    checkpoint_manager = CheckpointManager(
        checkpoint_dir=checkpoint_dir,
        config=pipeline_config,
        enabled=True
    )

    if resume:
        loaded = checkpoint_manager.load(checkpoint_path)
        if loaded:
            checkpoint_manager.display_resume_info(loaded)

from ..results import PipelineResults
results_collector = PipelineResults()
results_collector.start_time = time.time()

try:
    # Call pipeline with instrumentation
    run_count_t7(
        bam_file=bam_file,
        ref_file=ref_file,
        # ... all existing parameters ...

        # NEW: Pass instrumentation objects
        progress_tracker=progress_tracker,
        checkpoint_manager=checkpoint_manager,
        results_collector=results_collector,
        enable_instrumentation=True,
    )

    # Success
    results_collector.end_time = time.time()

    if checkpoint_manager:
        checkpoint_manager.mark_completed(success=True)

    if verbosity >= 1:
        results_collector.display_summary()

except KeyboardInterrupt:
    logger.warning("Pipeline interrupted by user (Ctrl+C)")
    if checkpoint_manager:
        checkpoint_manager.save()
        console.print(f"\n[yellow]✓ Checkpoint saved[/yellow]")
        console.print(f"[dim]Resume with: sheriff run --config {config} --resume[/dim]")
    raise

except Exception as e:
    logger.error(f"Pipeline failed: {e}", exc_info=True)
    if checkpoint_manager:
        checkpoint_manager.mark_completed(success=False, error=str(e))
    raise

finally:
    # Cleanup progress
    if progress_tracker:
        progress_tracker.__exit__(None, None, None)
```

---

## Phase 5: Error Handling & Edge Cases (45 min)

### 5.1 Checkpoint Compatibility

```python
# In CheckpointManager.load()

def _is_compatible(self, checkpoint):
    # Check version
    if checkpoint.version != current_version:
        logger.warning(f"Version mismatch: checkpoint v{checkpoint.version}, current v{current_version}")
        # Allow if minor version matches
        if checkpoint.version.split('.')[0:2] != current_version.split('.')[0:2]:
            return False

    # Check config hash
    if checkpoint.config_hash != self.compute_config_hash(self.config):
        logger.warning("Config changed since checkpoint - may not resume correctly")
        return False

    # Check outputs exist
    if checkpoint.outputs.filtered_bam:
        if not os.path.exists(checkpoint.outputs.filtered_bam):
            logger.error(f"Checkpoint output missing: {checkpoint.outputs.filtered_bam}")
            return False

    return True
```

### 5.2 Progress Bar Performance

```python
# Update progress intelligently (not every iteration)

# Bad: Update every read (millions of calls)
for read in reads:
    process(read)
    progress.update(task, advance=1)  # TOO SLOW

# Good: Update every N reads
UPDATE_INTERVAL = 1000
for i, read in enumerate(reads):
    process(read)
    if i % UPDATE_INTERVAL == 0:
        progress.update(task, advance=UPDATE_INTERVAL)

# Better: Update per chromosome or chunk
for chrom in chromosomes:
    reads = get_reads(chrom)
    process_all(reads)
    progress.update(task, advance=1)  # One update per chromosome
```

### 5.3 Memory Management

```python
# Don't store all data in checkpoints (too large)
# Instead, save to temp files and store paths

# Bad: Store data in memory
checkpoint.outputs.edit_counts = edit_counts  # Could be GB of data

# Good: Save to temp file
import pickle
temp_file = f"{checkpoint_dir}/edit_counts.pkl"
with open(temp_file, 'wb') as f:
    pickle.dump(edit_counts, f)
checkpoint.outputs.edit_counts_file = temp_file
```

---

## Phase 6: Testing Strategy (1 hour)

### 6.1 Unit Tests

```python
# tests/test_instrumentation.py

def test_progress_optional():
    """Test pipeline works without progress tracker"""
    run_count_t7(
        bam_file=test_bam,
        # ... params ...
        progress_tracker=None,  # Should work
    )

def test_checkpoint_save_load():
    """Test checkpoint round-trip"""
    # Run with checkpoints
    run_count_t7(..., checkpoint_manager=manager)

    # Verify checkpoint saved
    assert checkpoint_file.exists()

    # Load and validate
    loaded = manager.load()
    assert loaded.current_step == expected_step

def test_resume_skips_work():
    """Test resume actually skips completed steps"""
    # Run to checkpoint 2
    # Resume
    # Verify step 1 was skipped
```

### 6.2 Integration Test with Real Data

```bash
#!/bin/bash
# integration_test.sh

echo "=== Testing Sheriff Deep Integration ==="

# Test 1: Full run with all features
echo "Test 1: Full run with progress, checkpoints, results"
python -m sheriff run \
  --config examples/example-analysis.yaml \
  --enable-checkpoints \
  --log-file test.log \
  --verbosity 2

# Test 2: Interrupt and resume
echo "Test 2: Interrupt and resume"
timeout 30s python -m sheriff run \
  --config examples/example-analysis.yaml \
  --enable-checkpoints || true  # Will timeout

python -m sheriff run \
  --config examples/example-analysis.yaml \
  --resume

# Test 3: Verify outputs
echo "Test 3: Verify outputs exist"
test -f results/edit_site_info.txt || (echo "FAIL: edit_site_info.txt missing" && exit 1)
test -f results/*umi*.mtx || (echo "FAIL: UMI matrix missing" && exit 1)

echo "✓ All integration tests passed"
```

---

## Phase 7: Documentation (30 min)

### 7.1 Update Docstrings

```python
def run_count_t7(
    bam_file,
    ref_file,
    # ... params ...
    progress_tracker=None,
    checkpoint_manager=None,
    results_collector=None,
):
    """
    Run Sheriff CRISPR edit calling pipeline.

    Priority 3 Integration (New):
        The pipeline now supports real-time progress tracking, checkpointing,
        and metrics collection. These features are optional and have zero
        overhead when not used.

        Progress Tracking:
            Pass a PipelineProgress instance to see real-time progress bars:

            >>> from sheriff.progress import PipelineProgress
            >>> with PipelineProgress(verbosity=1) as progress:
            ...     run_count_t7(..., progress_tracker=progress)

        Checkpointing:
            Enable checkpoints to resume interrupted runs:

            >>> from sheriff.checkpoint import CheckpointManager
            >>> manager = CheckpointManager(".checkpoints", config=cfg, enabled=True)
            >>> run_count_t7(..., checkpoint_manager=manager)
            >>> # If interrupted, resume:
            >>> manager.load()
            >>> run_count_t7(..., checkpoint_manager=manager)

        Results:
            Collect pipeline metrics automatically:

            >>> from sheriff.results import PipelineResults
            >>> results = PipelineResults()
            >>> run_count_t7(..., results_collector=results)
            >>> results.display_summary()

    Parameters:
        ... existing params ...

        progress_tracker : PipelineProgress, optional
            Progress tracker for real-time updates. If None, no progress shown.

        checkpoint_manager : CheckpointManager, optional
            Checkpoint manager for resume capability. If None, no checkpoints.

        results_collector : PipelineResults, optional
            Results collector for metrics. If None, no metrics collected.

    Examples:
        Basic usage (backward compatible):
        >>> run_count_t7("sample.bam", "genome.fa", "barcodes.txt", "genes.gtf")

        With all Priority 3 features:
        >>> from sheriff.progress import PipelineProgress
        >>> from sheriff.checkpoint import CheckpointManager
        >>> from sheriff.results import PipelineResults
        >>>
        >>> with PipelineProgress(verbosity=1) as progress:
        ...     results = PipelineResults()
        ...     manager = CheckpointManager(".checkpoints", enabled=True)
        ...
        ...     run_count_t7(
        ...         "sample.bam", "genome.fa", "barcodes.txt", "genes.gtf",
        ...         progress_tracker=progress,
        ...         checkpoint_manager=manager,
        ...         results_collector=results
        ...     )
        ...
        ...     results.display_summary()
    """
```

### 7.2 Create Integration Guide

```markdown
# Sheriff Deep Integration Guide

## Real-Time Progress

Progress bars appear automatically when using the CLI:

```bash
sheriff run --config analysis.yaml

⠹ Finding barcoded edits ████████████░░░░  chr12/22 • 2m • ~1.5m remaining
```

## Checkpoint & Resume

Enable checkpoints for long runs:

```bash
# Run with checkpoints
sheriff run --config analysis.yaml --enable-checkpoints

# If interrupted (Ctrl+C, network drop, etc.):
sheriff run --config analysis.yaml --resume

♻️ Resuming from Checkpoint
  Progress: 65% complete
  Resuming from: UMI counting
```

## Results Summary

After completion, see full summary:

```
╭──────────── 📊 Sheriff Results ────────────╮
│  Input Reads: 12,456,789                  │
│  T7 Barcoded: 45,678 (0.37%)              │
│  Edit Sites: 23                           │
│  Runtime: 1h 23m                          │
╰───────────────────────────────────────────╯
```

## Programmatic Usage

```python
from sheriff.count_t7 import run_count_t7
from sheriff.progress import PipelineProgress
from sheriff.checkpoint import CheckpointManager
from sheriff.results import PipelineResults

# Initialize instrumentation
with PipelineProgress(verbosity=1) as progress:
    checkpoint_mgr = CheckpointManager(".checkpoints", enabled=True)
    results = PipelineResults()

    # Run pipeline
    run_count_t7(
        bam_file="sample.bam",
        ref_file="genome.fa",
        barcode_file="barcodes.txt",
        gtf_file="genes.gtf",
        progress_tracker=progress,
        checkpoint_manager=checkpoint_mgr,
        results_collector=results,
    )

    # Display results
    results.display_summary()
```
```

---

## Implementation Timeline

### Session 1 (2 hours):
- ✅ Phase 1: Architecture analysis (30 min)
- ✅ Phase 2: Function signatures (15 min)
- ⏳ Phase 3: Steps 1-2 instrumentation (1h 15min)

### Session 2 (2 hours):
- ⏳ Phase 3: Steps 3-5 instrumentation (1h 30min)
- ⏳ Phase 4: CLI integration (30 min)

### Session 3 (2 hours):
- ⏳ Phase 5: Error handling (45 min)
- ⏳ Phase 6: Testing (1 hour)
- ⏳ Phase 7: Documentation (30 min)

**Total: 6 hours** (can be done in 3 focused sessions)

---

## Success Criteria

✅ **Progress bars show during actual pipeline**
✅ **Checkpoints save after each major step**
✅ **Resume skips completed work correctly**
✅ **Results display real pipeline metrics**
✅ **Works with example_data/ test files**
✅ **Ctrl+C saves checkpoint and allows resume**
✅ **All features optional (backward compatible)**
✅ **Performance impact < 5%**
✅ **Zero breaking changes**

---

## Risk Mitigation

**Risk**: Breaking existing functionality
**Mitigation**:
- All features optional (default=None)
- Extensive testing at each phase
- Git commits after each phase

**Risk**: Performance degradation
**Mitigation**:
- Progress updates batched (not per-read)
- Checkpoints only every N minutes
- Instrumentation has if-guards (zero cost when disabled)

**Risk**: Memory overhead
**Mitigation**:
- Checkpoint to files, not memory
- Stream metrics, don't accumulate
- Cleanup temp files

**Risk**: Resume corruption
**Mitigation**:
- Validate checkpoint compatibility
- Atomic file writes
- Keep last N checkpoints as backup

---

## Expected Result

**True 10/10 CLI** with:
- ✅ Real-time progress (live bars during execution)
- ✅ Full resume capability (interrupt and continue)
- ✅ Live metrics (real pipeline data)
- ✅ Production-ready (error handling, logging)
- ✅ Backward compatible (zero breaking changes)
- ✅ Tested (works with real data)

**Honest rating after implementation: 10/10** ✨

Let's do this!
