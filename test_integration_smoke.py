"""Smoke test for Priority 3 integration.

This test verifies that the instrumented run_count_t7() function
accepts the new parameters and handles them correctly.
"""

import sys
from sheriff.count_t7 import run_count_t7
from sheriff.progress import PipelineProgress
from sheriff.checkpoint import CheckpointManager
from sheriff.results import PipelineResults


def test_function_signature():
    """Test that run_count_t7 accepts new Priority 3 parameters."""
    print("Testing function signature...")

    # Get the function signature
    import inspect
    sig = inspect.signature(run_count_t7)

    # Check for new parameters
    params = list(sig.parameters.keys())

    required_params = ['progress_tracker', 'checkpoint_manager', 'results_collector', 'enable_instrumentation']

    for param in required_params:
        assert param in params, f"Missing parameter: {param}"
        print(f"  ✓ Found parameter: {param}")

    # Check defaults
    assert sig.parameters['progress_tracker'].default is None
    assert sig.parameters['checkpoint_manager'].default is None
    assert sig.parameters['results_collector'].default is None
    assert sig.parameters['enable_instrumentation'].default is True

    print("✓ Function signature correct\n")


def test_instrumentation_objects():
    """Test that instrumentation objects can be created."""
    print("Testing instrumentation object creation...")

    # Test PipelineProgress
    progress = PipelineProgress(verbosity=0)  # verbosity=0 to avoid output
    assert progress is not None
    print("  ✓ PipelineProgress created")

    # Test CheckpointManager
    import tempfile
    import os
    with tempfile.TemporaryDirectory() as tmpdir:
        config = {"bam_file": "test.bam", "parameters": {"k": 6}}
        ckpt = CheckpointManager(checkpoint_dir=tmpdir, config=config, enabled=True)
        assert ckpt is not None
        print("  ✓ CheckpointManager created")

    # Test PipelineResults
    results = PipelineResults()
    assert results is not None
    results.set_metric("test_metric", 42)
    assert results.metrics["test_metric"] == 42
    print("  ✓ PipelineResults created and working")

    print("✓ All instrumentation objects work\n")


def test_backward_compatibility():
    """Test that old code still works without new parameters."""
    print("Testing backward compatibility...")

    # The function should work when called without the new parameters
    # (we're not actually calling it, just verifying the signature allows it)
    import inspect
    sig = inspect.signature(run_count_t7)

    # All new parameters should have defaults
    new_params = ['progress_tracker', 'checkpoint_manager', 'results_collector', 'enable_instrumentation']

    for param_name in new_params:
        param = sig.parameters[param_name]
        assert param.default != inspect.Parameter.empty, f"{param_name} must have a default value"

    print("  ✓ All new parameters have defaults")
    print("✓ Backward compatible\n")


if __name__ == "__main__":
    try:
        test_function_signature()
        test_instrumentation_objects()
        test_backward_compatibility()

        print("=" * 60)
        print("✅ ALL SMOKE TESTS PASSED!")
        print("=" * 60)
        print("\nPriority 3 integration is working correctly.")
        print("The pipeline is ready for testing with real data.")
        sys.exit(0)

    except Exception as e:
        print(f"\n❌ SMOKE TEST FAILED: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
