# Sheriff Priority 3 - Complete Integration Plan

## Current State (Honest Assessment)

### ✅ What Works
- Infrastructure: Checkpoint manager, progress tracker, results display
- CLI flags: --resume, --enable-checkpoints, --checkpoint
- Demo scripts: Standalone features work perfectly
- Tests: Component-level testing passes

### ❌ What's Missing
- **count_t7.py instrumentation**: Pipeline doesn't use any Priority 3 features
- **Real data testing**: Haven't run with actual BAM files
- **Error handling**: No error recovery in pipeline
- **Metrics collection**: No real pipeline metrics captured
- **Production validation**: Haven't tested resume after real interruption

## Integration Architecture

### 1. count_t7.py Structure Analysis

**Current major steps:**
```python
def run_count_t7(...):
    # 1. Input validation
    # 2. BAM filtering (chromosome-based)
    # 3. K-mer matching (find T7 barcodes)
    # 4. Edit site calling (cluster edits)
    # 5. UMI counting (count unique molecules)
    # 6. Gene counting (quantify expression)
    # 7. Output generation
```

**Each step needs:**
- Progress bar with read count
- Checkpoint save point
- Metric collection
- Error handling

### 2. Progress Tracking Integration

**Add to function signature:**
```python
def run_count_t7(
    ...,
    enable_progress: bool = True,
    progress_callback = None  # For testing
):
```

**Instrument key loops:**
```python
# Example: BAM filtering
if enable_progress:
    with PipelineProgress(verbosity=verbosity) as progress:
        task = progress.add_task("filtering", "Filtering BAM", total=total_reads)

        for read in bam_reader:
            # ... process read
            progress.update("filtering", advance=1)

        progress.complete_task("filtering")
```

**Steps to track:**
1. BAM Filtering (reads processed)
2. K-mer Matching (reads scanned)
3. Edit Site Calling (sites evaluated)
4. UMI Counting (genes processed)
5. Output Generation (files written)

### 3. Checkpoint Integration

**Checkpoint points:**
```python
checkpoints = {
    "validation": Before processing starts,
    "bam_filtering": After filtered BAM created,
    "kmer_matching": After T7 reads identified,
    "edit_calling": After edit sites called,
    "umi_counting": After UMI counts complete,
    "output": After all outputs written
}
```

**Resume logic:**
```python
if checkpoint_manager and checkpoint_manager.can_skip_step("bam_filtering"):
    logger.info("Skipping BAM filtering (already completed)")
    # Load outputs from checkpoint
    filtered_bam = checkpoint.outputs.filtered_bam
else:
    # Do BAM filtering
    filtered_bam = filter_bam(...)

    # Save checkpoint
    if checkpoint_manager:
        checkpoint_manager.update_outputs(filtered_bam=filtered_bam)
        checkpoint_manager.save()
```

**Checkpoint data to save:**
- Intermediate file paths
- Reads processed count
- Edit sites found
- Current step
- Runtime so far

### 4. Metrics Collection

**Metrics to collect:**
```python
metrics = {
    # Input
    "input_reads": total_reads,
    "input_bam_size": os.path.getsize(bam_file),

    # Filtering
    "filtered_reads": len(filtered_reads),
    "filter_rate": filtered / total,

    # K-mer matching
    "t7_barcoded_reads": len(barcoded_reads),
    "barcode_rate": barcoded / filtered,

    # Edit calling
    "candidate_edits": len(candidates),
    "canonical_edits": len(canonical),
    "cells_with_edits": len(unique_cells),

    # UMI counting
    "genes_quantified": len(genes),
    "umis_counted": total_umis,

    # Performance
    "bam_filtering_seconds": t_filter,
    "kmer_matching_seconds": t_kmer,
    "umi_counting_seconds": t_umi,
    "total_runtime_seconds": total_time
}
```

**Pass to results display:**
```python
results = PipelineResults()
results.start_time = start
# ... collect during pipeline
results.end_time = end
results.display_summary()
```

### 5. Error Handling

**Wrap pipeline in try/catch:**
```python
try:
    # Initialize
    checkpoint_manager = CheckpointManager(...) if enable_checkpoints else None
    results = PipelineResults()

    # Run pipeline
    with PipelineProgress(verbosity=verbosity) as progress:
        # ... all steps
        pass

    # Success
    if checkpoint_manager:
        checkpoint_manager.mark_completed(success=True)
    results.display_summary()

except KeyboardInterrupt:
    logger.warning("Pipeline interrupted by user")
    if checkpoint_manager:
        checkpoint_manager.save()  # Save current state
        logger.info(f"Checkpoint saved - resume with: sheriff run --resume")
    raise

except Exception as e:
    logger.error(f"Pipeline failed: {e}", exc_info=True)
    if checkpoint_manager:
        checkpoint_manager.mark_completed(success=False, error=str(e))
    raise
```

**Error scenarios to handle:**
1. File not found (after validation)
2. Out of memory
3. Corrupted BAM
4. Write permission errors
5. Keyboard interrupt (Ctrl+C)
6. Network errors (if reading from network storage)

### 6. Logging Integration

**Use structured logging throughout:**
```python
from sheriff.logging_config import get_logger

logger = get_logger("pipeline")

# Log major events
logger.info("Starting BAM filtering")
logger.info(f"Filtered {n_reads:,} reads in {elapsed:.1f}s")
logger.warning(f"Low barcode match rate: {rate:.2%}")
logger.error(f"Failed to process chromosome {chrom}: {e}")
```

**Log checkpoints:**
```python
logger.debug(f"Checkpoint saved: {checkpoint_path}")
logger.info(f"Resume progress: {percent}% complete")
```

**Log metrics:**
```python
logger.info(f"Pipeline metrics: {json.dumps(metrics, indent=2)}")
```

### 7. Testing Strategy

**Test with example data:**
```python
# Small test dataset
bam_file = "example_data/barcode_headAligned_anno.sorted.edit_regions_200kb.bam"
# This is ~30MB, should run in < 1 minute

# Test scenarios:
1. Full run with progress bars
2. Run with checkpoints enabled
3. Interrupt and resume
4. Run with --dry-run
5. Run with --log-file
6. Run with all features together
```

**Validation:**
- Progress bars appear and update
- Checkpoints created in .sheriff_checkpoints/
- Resume loads checkpoint and skips completed steps
- Results table shows real metrics
- Log file contains structured logs
- Errors are caught and logged

## Implementation Order

### Phase 1: Read and Understand count_t7.py (30 min)
- Map out all major steps
- Identify loop points for progress
- Find checkpoint opportunities
- Understand current error handling

### Phase 2: Add Progress Tracking (1 hour)
- Add progress parameter to function signature
- Wrap BAM reading loop with progress
- Add progress to k-mer loop
- Add progress to UMI counting
- Test that bars appear

### Phase 3: Add Checkpoint Logic (1 hour)
- Add checkpoint_manager parameter
- Add checkpoint after each major step
- Implement resume logic (skip completed steps)
- Test save/load cycle

### Phase 4: Metrics Collection (30 min)
- Add results parameter
- Collect counts during processing
- Record timings
- Test results display

### Phase 5: Error Handling (30 min)
- Wrap in try/catch
- Handle interrupts
- Save checkpoint on error
- Test error scenarios

### Phase 6: Integration Testing (1 hour)
- Test with example data
- Test all features together
- Test resume after interrupt
- Validate all metrics accurate

### Phase 7: Documentation (30 min)
- Update docstrings
- Add usage examples
- Document checkpoint format
- Create integration test

## Success Criteria

✅ **Progress bars show during actual pipeline run**
✅ **Checkpoints save after each step**
✅ **Resume skips completed work**
✅ **Results table shows real pipeline metrics**
✅ **Errors are caught and logged properly**
✅ **Works with example data end-to-end**
✅ **All features work together**
✅ **100% backward compatible (features are optional)**

## Risks and Mitigations

**Risk**: count_t7.py is complex, hard to instrument
**Mitigation**: Start with minimal instrumentation, expand gradually

**Risk**: Progress bars slow down pipeline
**Mitigation**: Make them optional, update every N reads not every read

**Risk**: Checkpoint overhead
**Mitigation**: Only save every 5 minutes, not every operation

**Risk**: Breaking existing functionality
**Mitigation**: All new features optional, default=off, extensive testing

**Risk**: Resume doesn't work correctly
**Mitigation**: Thorough checkpoint validation, test with interrupts

## Next Steps

1. Read count_t7.py completely
2. Create instrumented version
3. Test with real data
4. Iterate based on results
5. Commit when actually working end-to-end

This time, we do it right.
