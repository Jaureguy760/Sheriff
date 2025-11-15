"""Checkpoint manager for saving and loading pipeline state."""

import json
from pathlib import Path
from typing import Optional, Dict, Any
from datetime import datetime
from rich.console import Console

from .format import Checkpoint, PipelineStep, CheckpointMetrics, CheckpointOutputs

console = Console()


class CheckpointManager:
    """Manages pipeline checkpoints for resume functionality."""

    def __init__(
        self,
        checkpoint_dir: str = ".sheriff_checkpoints",
        config: Optional[Dict[str, Any]] = None,
        enabled: bool = False,
    ):
        """Initialize checkpoint manager.

        Args:
            checkpoint_dir: Directory to store checkpoints
            config: Pipeline configuration for hash computation
            enabled: Whether checkpointing is enabled
        """
        self.checkpoint_dir = Path(checkpoint_dir)
        self.config = config or {}
        self.enabled = enabled
        self.current_checkpoint: Optional[Checkpoint] = None

        if self.enabled:
            self.checkpoint_dir.mkdir(parents=True, exist_ok=True)

    def create_checkpoint(
        self,
        version: str,
        current_step: str,
        steps_completed: list,
        steps_remaining: list,
        outputs: Optional[CheckpointOutputs] = None,
        metrics: Optional[CheckpointMetrics] = None,
    ) -> Checkpoint:
        """Create a new checkpoint.

        Args:
            version: Sheriff version
            current_step: Current pipeline step
            steps_completed: List of completed steps
            steps_remaining: List of remaining steps
            outputs: Output files
            metrics: Pipeline metrics

        Returns:
            Created checkpoint
        """
        checkpoint = Checkpoint(
            version=version,
            config_hash=Checkpoint.compute_config_hash(self.config),
            current_step=current_step,
            steps_completed=steps_completed,
            steps_remaining=steps_remaining,
            outputs=outputs or CheckpointOutputs(),
            metrics=metrics or CheckpointMetrics(),
        )

        self.current_checkpoint = checkpoint
        return checkpoint

    def save(self, checkpoint: Optional[Checkpoint] = None) -> Optional[Path]:
        """Save checkpoint to disk.

        Args:
            checkpoint: Checkpoint to save (uses current if None)

        Returns:
            Path to saved checkpoint file, or None if disabled
        """
        if not self.enabled:
            return None

        checkpoint = checkpoint or self.current_checkpoint
        if not checkpoint:
            return None

        # Create filename with timestamp
        timestamp = checkpoint.timestamp.strftime("%Y%m%d_%H%M%S")
        filename = f"checkpoint_{timestamp}.json"
        filepath = self.checkpoint_dir / filename

        # Write checkpoint
        with open(filepath, "w") as f:
            f.write(checkpoint.to_json())

        return filepath

    def load(self, checkpoint_path: Optional[str] = None) -> Optional[Checkpoint]:
        """Load checkpoint from disk.

        Args:
            checkpoint_path: Specific checkpoint file, or None to load latest

        Returns:
            Loaded checkpoint, or None if not found
        """
        if not self.enabled:
            return None

        # Find checkpoint file
        if checkpoint_path:
            filepath = Path(checkpoint_path)
        else:
            filepath = self._find_latest_checkpoint()

        if not filepath or not filepath.exists():
            return None

        # Load and parse
        try:
            with open(filepath, "r") as f:
                checkpoint = Checkpoint.from_json(f.read())

            # Validate compatibility
            if not self._is_compatible(checkpoint):
                console.print(
                    f"[yellow]⚠ Checkpoint incompatible:[/yellow] "
                    f"config has changed (hash mismatch)"
                )
                return None

            self.current_checkpoint = checkpoint
            return checkpoint

        except Exception as e:
            console.print(f"[red]✗ Failed to load checkpoint:[/red] {e}")
            return None

    def _find_latest_checkpoint(self) -> Optional[Path]:
        """Find the most recent checkpoint file.

        Returns:
            Path to latest checkpoint, or None
        """
        if not self.checkpoint_dir.exists():
            return None

        checkpoints = list(self.checkpoint_dir.glob("checkpoint_*.json"))
        if not checkpoints:
            return None

        # Sort by modification time (most recent first)
        checkpoints.sort(key=lambda p: p.stat().st_mtime, reverse=True)
        return checkpoints[0]

    def _is_compatible(self, checkpoint: Checkpoint) -> bool:
        """Check if checkpoint is compatible with current config.

        Args:
            checkpoint: Checkpoint to validate

        Returns:
            True if compatible
        """
        # Check version (exact match required)
        from .. import __version__

        if checkpoint.version != __version__:
            console.print(
                f"[yellow]⚠ Version mismatch:[/yellow] "
                f"checkpoint is v{checkpoint.version}, current is v{__version__}"
            )
            return False

        # Check config hash
        current_hash = Checkpoint.compute_config_hash(self.config)
        if checkpoint.config_hash != current_hash:
            return False

        return True

    def cleanup_old_checkpoints(self, keep_last: int = 5):
        """Remove old checkpoint files, keeping most recent.

        Args:
            keep_last: Number of recent checkpoints to keep
        """
        if not self.checkpoint_dir.exists():
            return

        checkpoints = list(self.checkpoint_dir.glob("checkpoint_*.json"))
        if len(checkpoints) <= keep_last:
            return

        # Sort by modification time (most recent first)
        checkpoints.sort(key=lambda p: p.stat().st_mtime, reverse=True)

        # Remove old checkpoints
        for old_checkpoint in checkpoints[keep_last:]:
            old_checkpoint.unlink()

    def mark_completed(self, success: bool = True, error: Optional[str] = None):
        """Mark checkpoint as completed.

        Args:
            success: Whether pipeline completed successfully
            error: Error message if failed
        """
        if not self.enabled or not self.current_checkpoint:
            return

        self.current_checkpoint.status = "completed" if success else "failed"
        if error:
            self.current_checkpoint.error = error

        self.save()

        # Cleanup old checkpoints on success
        if success:
            self.cleanup_old_checkpoints()

    def update_metrics(self, **kwargs):
        """Update checkpoint metrics.

        Args:
            **kwargs: Metric key-value pairs to update
        """
        if not self.enabled or not self.current_checkpoint:
            return

        for key, value in kwargs.items():
            if hasattr(self.current_checkpoint.metrics, key):
                setattr(self.current_checkpoint.metrics, key, value)

    def update_outputs(self, **kwargs):
        """Update checkpoint outputs.

        Args:
            **kwargs: Output key-value pairs to update
        """
        if not self.enabled or not self.current_checkpoint:
            return

        for key, value in kwargs.items():
            if hasattr(self.current_checkpoint.outputs, key):
                setattr(self.current_checkpoint.outputs, key, value)

    def can_skip_step(self, step: str) -> bool:
        """Check if a step can be skipped (already completed in checkpoint).

        Args:
            step: Step to check

        Returns:
            True if step is already completed
        """
        if not self.enabled or not self.current_checkpoint:
            return False

        return step in self.current_checkpoint.steps_completed

    def display_resume_info(self, checkpoint: Checkpoint):
        """Display information about resumed checkpoint.

        Args:
            checkpoint: Checkpoint being resumed
        """
        from rich.panel import Panel
        from rich.table import Table

        # Create summary table
        table = Table(show_header=False, box=None, padding=(0, 1))
        table.add_column("Key", style="cyan")
        table.add_column("Value")

        table.add_row("Checkpoint Time", checkpoint.timestamp.strftime("%Y-%m-%d %H:%M:%S"))
        table.add_row("Progress", f"{checkpoint.get_progress_percent()}%")
        table.add_row("Completed Steps", ", ".join(checkpoint.steps_completed) or "None")
        table.add_row("Next Step", checkpoint.current_step)
        table.add_row("Remaining Steps", ", ".join(checkpoint.steps_remaining) or "None")

        panel = Panel(
            table,
            title="[bold green]♻️ Resuming from Checkpoint[/bold green]",
            border_style="green",
        )

        console.print()
        console.print(panel)
        console.print()
