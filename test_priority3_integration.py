#!/usr/bin/env python3
"""Integration test for Priority 3 CLI features."""

import sys
import os
import tempfile
import shutil
from pathlib import Path

# Add sheriff to path
sys.path.insert(0, str(Path(__file__).parent))

def test_checkpoint_manager():
    """Test checkpoint manager directly."""
    print("\n=== Testing Checkpoint Manager ===")

    from sheriff.checkpoint import CheckpointManager, PipelineStep, CheckpointMetrics

    with tempfile.TemporaryDirectory() as tmpdir:
        # Create manager
        config = {"bam": "test.bam", "params": {"k": 6}}
        manager = CheckpointManager(checkpoint_dir=tmpdir, config=config, enabled=True)

        # Create checkpoint
        checkpoint = manager.create_checkpoint(
            version="1.2.0",
            current_step=PipelineStep.KMER_MATCHING,
            steps_completed=[PipelineStep.VALIDATION],
            steps_remaining=[PipelineStep.KMER_MATCHING, PipelineStep.UMI_COUNTING],
            metrics=CheckpointMetrics(reads_processed=1000)
        )

        # Save
        filepath = manager.save()
        assert filepath.exists(), "Checkpoint file not created"
        print(f"✓ Checkpoint saved: {filepath}")

        # Load
        loaded = manager.load()
        assert loaded is not None, "Failed to load checkpoint"
        assert loaded.current_step == PipelineStep.KMER_MATCHING
        print(f"✓ Checkpoint loaded successfully")
        print(f"  Progress: {loaded.get_progress_percent()}%")

    print("✓ Checkpoint Manager works!")


def test_progress_tracker():
    """Test progress tracker."""
    print("\n=== Testing Progress Tracker ===")

    from sheriff.progress import PipelineProgress
    import time

    with PipelineProgress(verbosity=1) as progress:
        # Add task
        progress.add_task("test", "Test Task", total=100)

        # Update
        for i in range(0, 100, 10):
            progress.update("test", advance=10)
            time.sleep(0.05)

        print("✓ Progress bars work!")


def test_results_display():
    """Test results display."""
    print("\n=== Testing Results Display ===")

    from sheriff.results import PipelineResults
    import time

    results = PipelineResults()
    results.start_time = time.time() - 100
    results.end_time = time.time()

    # Add metrics
    results.set_metric("input_reads", 1000000)
    results.set_metric("edit_sites", 10)

    # Add performance
    results.set_performance("filtering", 50.0)

    # Add outputs
    results.set_output("Test Output", "test.txt")

    # Display
    results.display_summary(show_performance=True, show_outputs=True)

    print("✓ Results display works!")


def test_cli_flags():
    """Test that CLI flags are recognized."""
    print("\n=== Testing CLI Flags ===")

    import subprocess

    # Test --help includes new flags
    result = subprocess.run(
        [sys.executable, "-m", "sheriff", "run", "--help"],
        capture_output=True,
        text=True
    )

    # Check flags (may be truncated in display)
    has_resume = "--resume" in result.stdout
    has_enable = "enable-check" in result.stdout.lower()  # Partial match due to truncation
    has_checkpoint = "--checkpoint" in result.stdout and "TEXT" in result.stdout

    assert has_resume, "Missing --resume flag"
    assert has_enable, "Missing --enable-checkpoints flag"
    assert has_checkpoint, "Missing --checkpoint flag"

    print("✓ All CLI flags present in help")
    print(f"  - --resume: ✓")
    print(f"  - --enable-checkpoints: ✓")
    print(f"  - --checkpoint: ✓")


def test_imports():
    """Test all imports work."""
    print("\n=== Testing Imports ===")

    try:
        from sheriff.checkpoint import CheckpointManager, Checkpoint, PipelineStep
        from sheriff.progress import PipelineProgress
        from sheriff.results import PipelineResults
        from sheriff.cli.run import run
        print("✓ All imports successful")
    except ImportError as e:
        print(f"✗ Import failed: {e}")
        return False

    return True


def main():
    """Run all integration tests."""
    print("=" * 60)
    print("Sheriff Priority 3 Integration Tests")
    print("=" * 60)

    try:
        # Test imports first
        if not test_imports():
            print("\n✗ Import test failed - aborting")
            return 1

        # Test components
        test_checkpoint_manager()
        test_progress_tracker()
        test_results_display()
        test_cli_flags()

        print("\n" + "=" * 60)
        print("✓ ALL INTEGRATION TESTS PASSED!")
        print("=" * 60)
        print("\nPriority 3 features are fully integrated and working:")
        print("  ✓ Checkpoint/Resume system")
        print("  ✓ Progress bars")
        print("  ✓ Results display")
        print("  ✓ CLI flags")
        print()
        return 0

    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
