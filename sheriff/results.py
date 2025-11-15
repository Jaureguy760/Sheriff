"""Result summary tables for Sheriff pipeline."""

from typing import Dict, Any, Optional
from pathlib import Path
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.text import Text

console = Console()


class PipelineResults:
    """Collects and displays pipeline results."""

    def __init__(self):
        """Initialize results collector."""
        self.metrics: Dict[str, Any] = {}
        self.outputs: Dict[str, str] = {}
        self.performance: Dict[str, float] = {}
        self.start_time: Optional[float] = None
        self.end_time: Optional[float] = None

    def set_metric(self, key: str, value: Any):
        """Set a pipeline metric.

        Args:
            key: Metric name
            value: Metric value
        """
        self.metrics[key] = value

    def set_output(self, key: str, filepath: str):
        """Set an output file.

        Args:
            key: Output description
            filepath: Path to output file
        """
        self.outputs[key] = filepath

    def set_performance(self, key: str, seconds: float):
        """Set a performance timing.

        Args:
            key: Operation name
            seconds: Duration in seconds
        """
        self.performance[key] = seconds

    def display_summary(self, show_performance: bool = True, show_outputs: bool = True):
        """Display comprehensive results summary.

        Args:
            show_performance: Include performance metrics
            show_outputs: Include output file paths
        """
        # Calculate total runtime
        total_runtime = self._format_duration(
            self.end_time - self.start_time if (self.end_time and self.start_time) else 0
        )

        # Create main table
        table = Table(show_header=False, box=None, padding=(0, 2), expand=True)
        table.add_column("Content", justify="left")

        # Pipeline Summary Section
        summary_section = self._create_summary_section()
        table.add_row(summary_section)

        # Performance Section (if requested)
        if show_performance and self.performance:
            table.add_row("")  # Spacer
            perf_section = self._create_performance_section(total_runtime)
            table.add_row(perf_section)

        # Output Files Section (if requested)
        if show_outputs and self.outputs:
            table.add_row("")  # Spacer
            output_section = self._create_outputs_section()
            table.add_row(output_section)

        # Wrap in panel
        panel = Panel(
            table,
            title="[bold green]📊 Sheriff Results[/bold green]",
            border_style="green",
            padding=(1, 2),
        )

        console.print()
        console.print(panel)
        console.print()

    def _create_summary_section(self) -> Table:
        """Create pipeline summary section.

        Returns:
            Summary table
        """
        table = Table(show_header=True, header_style="bold cyan", box=None)
        table.add_column("📊 Pipeline Summary", style="bold")
        table.add_column("Value", justify="right", style="green")

        # Add metrics
        metric_order = [
            ("input_reads", "Input Reads"),
            ("filtered_reads", "Filtered Reads"),
            ("barcoded_reads", "T7 Barcoded Reads"),
            ("edit_sites", "Edit Sites Detected"),
            ("cells_with_edits", "Cells with Edits"),
            ("genes_quantified", "Genes Quantified"),
        ]

        for key, label in metric_order:
            if key in self.metrics:
                value = self.metrics[key]
                formatted_value = self._format_metric(key, value)
                table.add_row(label, formatted_value)

        return table

    def _create_performance_section(self, total_runtime: str) -> Table:
        """Create performance section.

        Args:
            total_runtime: Formatted total runtime string

        Returns:
            Performance table
        """
        table = Table(show_header=True, header_style="bold cyan", box=None)
        table.add_column("⏱️  Performance", style="bold")
        table.add_column("Duration", justify="right", style="yellow")
        table.add_column("Note", style="dim")

        # Total runtime
        table.add_row("Total Runtime", total_runtime, "")

        # Individual step timings
        for step_name, duration_seconds in self.performance.items():
            duration_str = self._format_duration(duration_seconds)
            note = self._get_performance_note(step_name)
            # Format step name nicely
            display_name = "  └─ " + step_name.replace("_", " ").title()
            table.add_row(display_name, duration_str, note)

        return table

    def _create_outputs_section(self) -> Table:
        """Create outputs section.

        Returns:
            Outputs table
        """
        table = Table(show_header=True, header_style="bold cyan", box=None)
        table.add_column("📁 Output Files", style="bold")
        table.add_column("Path", style="blue")

        for description, filepath in self.outputs.items():
            # Check if file exists and get size
            path = Path(filepath)
            if path.exists():
                size = self._format_file_size(path.stat().st_size)
                filepath_with_size = f"{filepath} ({size})"
            else:
                filepath_with_size = filepath

            table.add_row(description, filepath_with_size)

        return table

    def _format_metric(self, key: str, value: Any) -> str:
        """Format a metric value for display.

        Args:
            key: Metric key
            value: Metric value

        Returns:
            Formatted string
        """
        # Handle percentages
        if isinstance(value, tuple) and len(value) == 2:
            count, total = value
            percent = (count / total * 100) if total > 0 else 0
            return f"{count:,} ({percent:.2f}%)"

        # Handle regular numbers
        if isinstance(value, (int, float)):
            if key in ["barcoded_reads", "edit_sites"]:
                # These typically have percentages - check if parent metric exists
                if "filtered_reads" in self.metrics:
                    total = self.metrics["filtered_reads"]
                    percent = (value / total * 100) if total > 0 else 0
                    return f"{value:,} ({percent:.2f}%)"
            return f"{value:,}"

        return str(value)

    @staticmethod
    def _format_duration(seconds: float) -> str:
        """Format duration in seconds to human-readable string.

        Args:
            seconds: Duration in seconds

        Returns:
            Formatted duration string
        """
        if seconds < 60:
            return f"{int(seconds)}s"
        elif seconds < 3600:
            minutes = int(seconds // 60)
            secs = int(seconds % 60)
            return f"{minutes}m {secs}s" if secs > 0 else f"{minutes}m"
        else:
            hours = int(seconds // 3600)
            minutes = int((seconds % 3600) // 60)
            return f"{hours}h {minutes}m" if minutes > 0 else f"{hours}h"

    @staticmethod
    def _format_file_size(size_bytes: int) -> str:
        """Format file size in bytes to human-readable string.

        Args:
            size_bytes: Size in bytes

        Returns:
            Formatted size string
        """
        for unit in ["B", "KB", "MB", "GB", "TB"]:
            if size_bytes < 1024.0:
                return f"{size_bytes:.1f} {unit}"
            size_bytes /= 1024.0
        return f"{size_bytes:.1f} PB"

    @staticmethod
    def _get_performance_note(step_name: str) -> str:
        """Get performance note for a step.

        Args:
            step_name: Step name

        Returns:
            Note string
        """
        # Add notes about Rust acceleration
        if "bam_filtering" in step_name.lower():
            return "Rust 10x speedup"
        elif "kmer_matching" in step_name.lower():
            return "Rust 75x speedup"
        return ""

    def export_json(self, filepath: str):
        """Export results to JSON file.

        Args:
            filepath: Path to output JSON file
        """
        import json

        data = {
            "metrics": self.metrics,
            "outputs": self.outputs,
            "performance": self.performance,
            "total_runtime_seconds": self.end_time - self.start_time if (self.end_time and self.start_time) else 0,
        }

        with open(filepath, "w") as f:
            json.dump(data, f, indent=2)

        console.print(f"[dim]Results exported to: {filepath}[/dim]")
