"""Input validation subcommand."""

import os
from typing import Optional
from typing_extensions import Annotated
import typer
from rich.console import Console
from rich.table import Table

console = Console()


def validate_command(
    bam_file: Annotated[Optional[str], typer.Argument(help="BAM file")] = None,
    ref_file: Annotated[Optional[str], typer.Argument(help="Reference genome FASTA")] = None,
    barcode_file: Annotated[Optional[str], typer.Argument(help="Barcode whitelist")] = None,
    gtf_file: Annotated[Optional[str], typer.Argument(help="GTF annotation file")] = None,
    config: Annotated[Optional[str], typer.Option("--config", "-c", help="Config file")] = None,
    blacklist_file: Annotated[Optional[str], typer.Option("--blacklist", help="Blacklist BED file")] = None,
    whitelist_file: Annotated[Optional[str], typer.Option("--whitelist", help="Whitelist BED file")] = None,
    cnv_file: Annotated[Optional[str], typer.Option("--cnv", help="CNV file")] = None,
):
    """Validate input files without running the pipeline.

    Can validate from config file or individual files.

    Examples:
        sheriff validate --config sheriff-config.yaml
        sheriff validate sample.bam ref.fa barcodes.txt genes.gtf
    """
    from ..config.loader import ConfigLoader

    files_to_check = {}

    # Load from config if provided
    if config:
        try:
            loader = ConfigLoader(config)
            cfg = loader.config
            files_to_check = {
                "BAM file": cfg.input.bam_file,
                "Reference FASTA": cfg.input.ref_file,
                "Barcode file": cfg.input.barcode_file,
                "GTF file": cfg.input.gtf_file,
            }
            if cfg.filters.blacklist_file:
                files_to_check["Blacklist file"] = cfg.filters.blacklist_file
            if cfg.filters.whitelist_file:
                files_to_check["Whitelist file"] = cfg.filters.whitelist_file
            if cfg.filters.cnv_file:
                files_to_check["CNV file"] = cfg.filters.cnv_file
        except (FileNotFoundError, ValueError) as e:
            console.print(f"[red]Config error:[/red] {e}")
            raise typer.Exit(code=1)
    else:
        # Check individual files
        if bam_file:
            files_to_check["BAM file"] = bam_file
        if ref_file:
            files_to_check["Reference FASTA"] = ref_file
        if barcode_file:
            files_to_check["Barcode file"] = barcode_file
        if gtf_file:
            files_to_check["GTF file"] = gtf_file
        if blacklist_file:
            files_to_check["Blacklist file"] = blacklist_file
        if whitelist_file:
            files_to_check["Whitelist file"] = whitelist_file
        if cnv_file:
            files_to_check["CNV file"] = cnv_file

    if not files_to_check:
        console.print("[yellow]No files to validate. Provide --config or file arguments.[/yellow]")
        raise typer.Exit(code=1)

    # Validate files
    table = Table(title="Input Validation Results", show_header=True, header_style="bold cyan")
    table.add_column("File Type", style="dim")
    table.add_column("Path")
    table.add_column("Status")
    table.add_column("Size", justify="right")

    all_valid = True
    for file_type, file_path in files_to_check.items():
        if os.path.exists(file_path):
            size = os.path.getsize(file_path)
            size_str = _format_size(size)
            table.add_row(file_type, file_path, "[green]✓ Exists[/green]", size_str)
        else:
            table.add_row(file_type, file_path, "[red]✗ Not found[/red]", "-")
            all_valid = False

    console.print(table)

    if all_valid:
        console.print("\n[green bold]✓ All files valid![/green bold]")
    else:
        console.print("\n[red bold]✗ Some files are missing[/red bold]")
        raise typer.Exit(code=1)


def _format_size(size_bytes: int) -> str:
    """Format file size in human-readable format."""
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} PB"
