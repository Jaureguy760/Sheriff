"""Progress tracking with Rich for Sheriff pipeline."""

import time
from typing import Optional, Dict
from contextlib import contextmanager
from rich.console import Console
from rich.progress import (
    Progress,
    SpinnerColumn,
    BarColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)
from rich.panel import Panel
from rich.live import Live
from rich.table import Table

console = Console()


class PipelineProgress:
    """Manages progress bars and live status for Sheriff pipeline."""

    def __init__(self, verbosity: int = 1):
        """Initialize progress tracker.

        Args:
            verbosity: Verbosity level (0=quiet, 1=normal, 2=verbose)
        """
        self.verbosity = verbosity
        self.progress: Optional[Progress] = None
        self.tasks: Dict[str, int] = {}
        self.start_time = time.time()
        self.current_status = ""
        self.live: Optional[Live] = None

    def __enter__(self):
        """Enter context manager."""
        if self.verbosity >= 1:
            self.progress = Progress(
                SpinnerColumn(),
                TextColumn("[bold blue]{task.description}"),
                BarColumn(complete_style="green", finished_style="bold green"),
                TaskProgressColumn(),
                TextColumn("•"),
                TimeElapsedColumn(),
                TextColumn("•"),
                TimeRemainingColumn(),
                console=console,
            )
            self.progress.__enter__()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit context manager."""
        if self.progress:
            self.progress.__exit__(exc_type, exc_val, exc_tb)
        return False

    def add_task(self, name: str, description: str, total: Optional[int] = None) -> int:
        """Add a new progress task.

        Args:
            name: Task identifier
            description: Task description
            total: Total units of work (None for indeterminate)

        Returns:
            Task ID
        """
        if not self.progress:
            return -1

        task_id = self.progress.add_task(description, total=total)
        self.tasks[name] = task_id
        return task_id

    def update(self, name: str, advance: int = 1, **kwargs):
        """Update task progress.

        Args:
            name: Task identifier
            advance: Amount to advance
            **kwargs: Additional updates (description, total, etc.)
        """
        if not self.progress or name not in self.tasks:
            return

        self.progress.update(self.tasks[name], advance=advance, **kwargs)

    def complete_task(self, name: str):
        """Mark task as complete.

        Args:
            name: Task identifier
        """
        if not self.progress or name not in self.tasks:
            return

        task_id = self.tasks[name]
        self.progress.update(task_id, completed=self.progress.tasks[task_id].total or 100)

    @contextmanager
    def live_status(self, title: str = "Current Operation"):
        """Context manager for live updating status panel.

        Args:
            title: Panel title

        Example:
            with progress.live_status("Processing"):
                progress.update_status("Step 1", "Details...")
        """
        if self.verbosity < 1:
            yield
            return

        def generate_status_panel():
            """Generate status panel content."""
            table = Table(show_header=False, box=None, padding=(0, 1))
            table.add_column("Key", style="dim")
            table.add_column("Value", style="bold")

            if self.current_status:
                for line in self.current_status.split("\n"):
                    if ":" in line:
                        key, value = line.split(":", 1)
                        table.add_row(key.strip(), value.strip())

            return Panel(table, title=f"[bold cyan]{title}[/bold cyan]", border_style="cyan")

        # Use Live display
        with Live(generate_status_panel(), console=console, refresh_per_second=4) as live:
            self.live = live
            try:
                yield self
            finally:
                self.live = None

    def update_status(self, **kwargs):
        """Update live status display.

        Args:
            **kwargs: Key-value pairs to display
        """
        if not self.live:
            return

        # Format status lines
        lines = [f"{key}: {value}" for key, value in kwargs.items()]
        self.current_status = "\n".join(lines)

        # Update live display
        table = Table(show_header=False, box=None, padding=(0, 1))
        table.add_column("Key", style="dim")
        table.add_column("Value", style="bold")

        for key, value in kwargs.items():
            table.add_row(str(key), str(value))

        panel = Panel(table, title="[bold cyan]Current Operation[/bold cyan]", border_style="cyan")
        self.live.update(panel)

    def get_elapsed_time(self) -> str:
        """Get elapsed time since start.

        Returns:
            Formatted elapsed time string
        """
        elapsed = time.time() - self.start_time
        hours = int(elapsed // 3600)
        minutes = int((elapsed % 3600) // 60)
        seconds = int(elapsed % 60)

        if hours > 0:
            return f"{hours}h {minutes}m {seconds}s"
        elif minutes > 0:
            return f"{minutes}m {seconds}s"
        else:
            return f"{seconds}s"

    @staticmethod
    def format_number(num: int) -> str:
        """Format large numbers with comma separators.

        Args:
            num: Number to format

        Returns:
            Formatted string
        """
        return f"{num:,}"
