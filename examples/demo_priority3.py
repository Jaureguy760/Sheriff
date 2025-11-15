#!/usr/bin/env python3
"""Demo script showing Priority 3 features: Checkpoint/Resume, Progress, Results."""

import time
from sheriff.checkpoint import CheckpointManager, Checkpoint, PipelineStep, CheckpointMetrics, CheckpointOutputs
from sheriff.progress import PipelineProgress
from sheriff.results import PipelineResults


def demo_checkpoint_system():
    """Demonstrate checkpoint save/load/resume functionality."""
    print("\n" + "=" * 60)
    print("DEMO 1: Checkpoint System")
    print("=" * 60 + "\n")

    # Simulate pipeline configuration
    config = {
        "bam_file": "sample.bam",
        "ref_file": "genome.fa",
        "parameters": {"k": 6, "edit_dist": 140},
    }

    # Create checkpoint manager
    manager = CheckpointManager(checkpoint_dir=".demo_checkpoints", config=config, enabled=True)

    # Create initial checkpoint (simulating pipeline in progress)
    checkpoint = manager.create_checkpoint(
        version="1.2.0",
        current_step=PipelineStep.KMER_MATCHING,
        steps_completed=[PipelineStep.VALIDATION, PipelineStep.BAM_FILTERING],
        steps_remaining=[PipelineStep.KMER_MATCHING, PipelineStep.EDIT_CALLING, PipelineStep.UMI_COUNTING],
        metrics=CheckpointMetrics(reads_processed=1_250_000, reads_filtered=850_000, runtime_seconds=125.3),
        outputs=CheckpointOutputs(filtered_bam="results/filtered.bam"),
    )

    print(f"✓ Created checkpoint at {checkpoint.get_progress_percent()}% completion")
    print(f"  Current step: {checkpoint.current_step}")
    print(f"  Steps completed: {', '.join(checkpoint.steps_completed)}")

    # Save checkpoint
    filepath = manager.save()
    print(f"\n✓ Saved checkpoint to: {filepath}")

    # Simulate loading checkpoint (resuming)
    print("\n--- Simulating Resume ---\n")
    loaded_checkpoint = manager.load()
    if loaded_checkpoint:
        manager.display_resume_info(loaded_checkpoint)
        print(f"✓ Successfully loaded checkpoint from {loaded_checkpoint.timestamp}")

    # Cleanup
    import shutil

    shutil.rmtree(".demo_checkpoints", ignore_errors=True)
    print("\n✓ Demo checkpoint directory cleaned up")


def demo_progress_tracking():
    """Demonstrate Rich progress bars and live status."""
    print("\n" + "=" * 60)
    print("DEMO 2: Progress Tracking")
    print("=" * 60 + "\n")

    with PipelineProgress(verbosity=1) as progress:
        # Simulate BAM filtering
        total_reads = 1_000_000
        task = progress.add_task("bam_filtering", "BAM Filtering", total=total_reads)

        print("Simulating BAM filtering with progress bar...")
        for i in range(0, total_reads, 50_000):
            time.sleep(0.1)  # Simulate work
            progress.update("bam_filtering", advance=50_000)

        progress.complete_task("bam_filtering")
        print("\n✓ BAM filtering complete")

        # Simulate k-mer matching with live status
        task2 = progress.add_task("kmer_matching", "K-mer Matching", total=total_reads)

        print("\nSimulating k-mer matching with live status...")
        for chrom in ["chr1", "chr2", "chr3", "chr4", "chr5"]:
            chunk_size = total_reads // 5
            for i in range(0, chunk_size, 25_000):
                time.sleep(0.05)  # Simulate work
                progress.update("kmer_matching", advance=25_000)

        progress.complete_task("kmer_matching")
        print("\n✓ K-mer matching complete")


def demo_results_summary():
    """Demonstrate Rich result summary tables."""
    print("\n" + "=" * 60)
    print("DEMO 3: Results Summary")
    print("=" * 60 + "\n")

    # Create results collector
    results = PipelineResults()
    results.start_time = time.time() - 4980  # Simulate 1h 23m runtime
    results.end_time = time.time()

    # Add pipeline metrics
    results.set_metric("input_reads", 12_456_789)
    results.set_metric("filtered_reads", 8_234_123)
    results.set_metric("barcoded_reads", 45_678)
    results.set_metric("edit_sites", 23)
    results.set_metric("cells_with_edits", 1_234)
    results.set_metric("genes_quantified", 18_452)

    # Add performance metrics
    results.set_performance("bam_filtering", 720.5)
    results.set_performance("kmer_matching", 2700.3)
    results.set_performance("umi_counting", 1560.0)

    # Add output files
    results.set_output("Edit Sites", "results/edit_site_info.txt")
    results.set_output("UMI Counts", "results/umi_counts.mtx")
    results.set_output("Gene Counts", "results/gene_counts.mtx")
    results.set_output("Log File", "results/sheriff.log")

    # Display summary
    print("Displaying comprehensive results summary:\n")
    results.display_summary(show_performance=True, show_outputs=True)


def main():
    """Run all demos."""
    print("\n")
    print("╔" + "═" * 58 + "╗")
    print("║" + " " * 10 + "Sheriff Priority 3 Features Demo" + " " * 15 + "║")
    print("║" + " " * 58 + "║")
    print("║  Demonstrating:                                          ║")
    print("║    • Checkpoint/Resume System                            ║")
    print("║    • Progress Bars & Live Status                         ║")
    print("║    • Results Summary Tables                              ║")
    print("╚" + "═" * 58 + "╝")

    try:
        # Demo 1: Checkpoint system
        demo_checkpoint_system()

        time.sleep(1)

        # Demo 2: Progress tracking
        demo_progress_tracking()

        time.sleep(1)

        # Demo 3: Results summary
        demo_results_summary()

        print("\n" + "=" * 60)
        print("✓ All demos completed successfully!")
        print("=" * 60 + "\n")

        print("Key Features Demonstrated:")
        print("  ✓ Checkpointing - Save pipeline state for resume")
        print("  ✓ Progress Bars - Live feedback during execution")
        print("  ✓ Live Status - Real-time operation details")
        print("  ✓ Results Tables - Beautiful summary at completion")
        print("\nThese features bring Sheriff to 10/10 CLI excellence!")
        print()

    except KeyboardInterrupt:
        print("\n\n⚠ Demo interrupted by user")
    except Exception as e:
        print(f"\n\n✗ Demo failed: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
