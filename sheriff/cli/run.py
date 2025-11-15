"""Main pipeline run subcommand."""

import os
from typing import Optional
from typing_extensions import Annotated
import typer
from rich.console import Console

from .. import __version__
from ..count_t7 import run_count_t7
from ..config.loader import ConfigLoader
from ..logging_config import setup_logging, get_logger

console = Console()
app = typer.Typer(help="Run the Sheriff CRISPR edit calling pipeline")


@app.command()
def run(
    bam_file: Annotated[Optional[str], typer.Argument(help="BAM file")] = None,
    ref_file: Annotated[Optional[str], typer.Argument(help="Fasta containing ref genome")] = None,
    barcode_file: Annotated[Optional[str], typer.Argument(help="Barcode whitelist file")] = None,
    gtf_file: Annotated[Optional[str], typer.Argument(help="GTF file containing relevant gene data")] = None,
    config: Annotated[
        Optional[str], typer.Option("--config", "-c", help="Load configuration from YAML file")
    ] = None,
    t7_barcode: Annotated[
        Optional[str], typer.Option("--t7", "--t7_barcode", help="Target/query barcode sequence")
    ] = "GGGAGAGTAT",
    blacklist_file: Annotated[Optional[str], typer.Option("--blacklist", help="Blacklist BED file")] = None,
    whitelist_file: Annotated[Optional[str], typer.Option("--whitelist", help="Whitelist BED file")] = None,
    k: Annotated[
        Optional[int], typer.Option("-k", "--kmer", help="Size of kmers for pattern matching")
    ] = 6,
    edit_dist: Annotated[
        Optional[int], typer.Option("--edit_dist", "--edist", help="+/- distance from edit site grouping")
    ] = 140,
    edit_site_rev_comp_filt: Annotated[
        Optional[bool],
        typer.Option(
            "--bidirectional_inserts/--no-bidirectional_inserts",
            help="Require bidirectional donor insertion evidence",
        ),
    ] = True,
    stranded_edit_dist: Annotated[
        Optional[int],
        typer.Option("--stranded_edit_dist", help="Max distance between forward and reverse edit sites"),
    ] = 15,
    edit_site_min_cells: Annotated[
        Optional[int], typer.Option("--edit_site_min_cells", help="Minimum cells per edit site")
    ] = 3,
    nonbc_edit_dist: Annotated[
        Optional[int],
        typer.Option("--nonbc_edit_dist", "--nonbc", help="+/- distance to assign non-barcoded reads"),
    ] = 1000,
    ploidy: Annotated[Optional[int], typer.Option("--ploidy", help="Genome ploidy")] = 2,
    cnv_file: Annotated[
        Optional[str], typer.Option("--cnv", "--cnv_file", help="Copy-number-variation file (bedGraph)")
    ] = None,
    blacklist_seqs: Annotated[
        Optional[str], typer.Option("--blacklist_seqs", help="File of blacklist sequences")
    ] = None,
    mrna_count_mode: Annotated[
        Optional[str], typer.Option("--mrna_count_mode", help="Mode for quantifying gene expression")
    ] = "all",
    uncorrected_gene_count: Annotated[
        Optional[bool],
        typer.Option("--uncorrected_count/--no-uncorrected_count", help="Output uncorrected gene counts"),
    ] = False,
    outdir: Annotated[
        Optional[str], typer.Option("-o", "--out", "--outdir", help="Output directory")
    ] = None,
    verbosity: Annotated[
        Optional[int], typer.Option("--v", "-v", "--verbosity", help="Verbosity level (0-2)")
    ] = 1,
    n_cpus: Annotated[Optional[int], typer.Option("--cpu", help="Number of CPUs")] = 1,
    chunk_size_mb: Annotated[
        Optional[int], typer.Option("--chunk", help="Chunk size in MB for processing")
    ] = 250,
    dry_run: Annotated[
        Optional[bool], typer.Option("--dry-run", "--dry_run", help="Preview without executing")
    ] = False,
    resume: Annotated[
        Optional[bool], typer.Option("--resume", help="Resume from last checkpoint if available")
    ] = False,
    enable_checkpoints: Annotated[
        Optional[bool], typer.Option("--enable-checkpoints", help="Enable automatic checkpointing")
    ] = False,
    checkpoint_path: Annotated[
        Optional[str], typer.Option("--checkpoint", help="Specific checkpoint file to resume from")
    ] = None,
    log_file: Annotated[Optional[str], typer.Option("--log-file", help="Log file path")] = None,
    log_level: Annotated[
        Optional[str], typer.Option("--log-level", help="Log level (DEBUG, INFO, WARNING, ERROR)")
    ] = "INFO",
    profile: Annotated[Optional[bool], typer.Option("--profile", help="Enable profiling mode")] = False,
):
    """Run the Sheriff CRISPR edit calling pipeline.

    Can be run with individual arguments or a config file (config takes precedence).

    Examples:
        # With config file
        sheriff run --config sheriff-config.yaml

        # With arguments
        sheriff run sample.bam ref.fa barcodes.txt genes.gtf

        # Mix config and arguments (args override config)
        sheriff run --config base.yaml --n_cpus 16

        # Enable checkpointing for long runs
        sheriff run --config analysis.yaml --enable-checkpoints

        # Resume from last checkpoint
        sheriff run --config analysis.yaml --resume

        # Resume from specific checkpoint
        sheriff run --config analysis.yaml --checkpoint .sheriff_checkpoints/checkpoint_20251115_102345.json
    """
    # Setup logging first
    logger = setup_logging(log_file=log_file, log_level=log_level)
    logger.info(f"Sheriff v{__version__} - Starting")

    # Handle config file if provided
    if config:
        try:
            loader = ConfigLoader(config)
            logger.info(f"Loaded configuration from: {config}")

            # Merge config with CLI arguments
            cli_kwargs = {
                "bam_file": bam_file,
                "ref_file": ref_file,
                "barcode_file": barcode_file,
                "gtf_file": gtf_file,
                "t7_barcode": t7_barcode,
                "blacklist_file": blacklist_file,
                "whitelist_file": whitelist_file,
                "k": k,
                "edit_dist": edit_dist,
                "edit_site_rev_comp_filt": edit_site_rev_comp_filt,
                "stranded_edit_dist": stranded_edit_dist,
                "edit_site_min_cells": edit_site_min_cells,
                "nonbc_edit_dist": nonbc_edit_dist,
                "ploidy": ploidy,
                "cnv_file": cnv_file,
                "blacklist_seqs": blacklist_seqs,
                "mrna_count_mode": mrna_count_mode,
                "uncorrected_gene_count": uncorrected_gene_count,
                "outdir": outdir,
                "verbosity": verbosity,
                "n_cpus": n_cpus,
                "chunk_size_mb": chunk_size_mb,
            }
            merged = loader.merge_with_cli(cli_kwargs)

            # Extract values from merged config
            bam_file = merged["bam_file"]
            ref_file = merged["ref_file"]
            barcode_file = merged["barcode_file"]
            gtf_file = merged["gtf_file"]
            t7_barcode = merged["t7_barcode"]
            blacklist_file = merged["blacklist_file"]
            whitelist_file = merged["whitelist_file"]
            k = merged["k"]
            edit_dist = merged["edit_dist"]
            edit_site_rev_comp_filt = merged["edit_site_rev_comp_filt"]
            stranded_edit_dist = merged["stranded_edit_dist"]
            edit_site_min_cells = merged["edit_site_min_cells"]
            nonbc_edit_dist = merged["nonbc_edit_dist"]
            ploidy = merged["ploidy"]
            cnv_file = merged["cnv_file"]
            blacklist_seqs = merged["blacklist_seqs"]
            mrna_count_mode = merged["mrna_count_mode"]
            uncorrected_gene_count = merged["uncorrected_gene_count"]
            outdir = merged["outdir"]
            verbosity = merged["verbosity"]
            n_cpus = merged["n_cpus"]
            chunk_size_mb = merged["chunk_size_mb"]

        except (FileNotFoundError, ValueError) as e:
            logger.error(f"Config error: {e}")
            console.print(f"[red]✗ Config error:[/red] {e}")
            raise typer.Exit(code=1)
    else:
        # Validate required arguments when not using config
        if not all([bam_file, ref_file, barcode_file, gtf_file]):
            console.print("[red]✗ Error:[/red] Required arguments missing")
            console.print("Provide either --config or all of: bam_file ref_file barcode_file gtf_file")
            raise typer.Exit(code=1)

    # Validate input files
    _validate_inputs(bam_file, ref_file, barcode_file, gtf_file, blacklist_file, whitelist_file, cnv_file, blacklist_seqs, logger)

    # Dry-run mode
    if dry_run:
        _show_dry_run(
            bam_file,
            ref_file,
            barcode_file,
            gtf_file,
            t7_barcode,
            k,
            edit_dist,
            edit_site_min_cells,
            n_cpus,
            outdir,
            blacklist_file,
            whitelist_file,
            cnv_file,
        )
        return

    # Show start message
    if verbosity >= 1:
        console.print(f"\n[green bold]🚀 Sheriff v{__version__} - Starting analysis[/green bold]")
        console.print(f"  BAM: {bam_file}")
        console.print(f"  Output: {outdir if outdir else 'Current directory'}\n")

    logger.info("Running Sheriff pipeline")
    logger.info(f"Input BAM: {bam_file}")
    logger.info(f"CPUs: {n_cpus}")

    # Initialize checkpoint manager (if enabled or resuming)
    checkpoint_manager = None
    resumed_checkpoint = None
    if enable_checkpoints or resume:
        from ..checkpoint import CheckpointManager

        # Build config dict for checkpointing
        pipeline_config = {
            "bam_file": bam_file,
            "ref_file": ref_file,
            "barcode_file": barcode_file,
            "gtf_file": gtf_file,
            "parameters": {
                "t7_barcode": t7_barcode,
                "k": k,
                "edit_dist": edit_dist,
                "n_cpus": n_cpus,
            },
        }

        checkpoint_dir = os.path.join(outdir if outdir else ".", ".sheriff_checkpoints")
        checkpoint_manager = CheckpointManager(
            checkpoint_dir=checkpoint_dir, config=pipeline_config, enabled=True
        )

        # Try to resume if requested
        if resume:
            resumed_checkpoint = checkpoint_manager.load(checkpoint_path)
            if resumed_checkpoint:
                checkpoint_manager.display_resume_info(resumed_checkpoint)
                logger.info(f"Resuming from checkpoint: {resumed_checkpoint.get_progress_percent()}% complete")
            elif checkpoint_path:
                # User specified a checkpoint file that doesn't exist
                console.print(f"[yellow]⚠ Checkpoint file not found: {checkpoint_path}[/yellow]")
                console.print("[yellow]Starting from beginning...[/yellow]")
            else:
                # No checkpoint found, start fresh
                console.print("[dim]No checkpoint found, starting from beginning...[/dim]")

    # Initialize results collector
    from ..results import PipelineResults
    import time

    results = PipelineResults()
    results.start_time = time.time()

    # Run the pipeline
    try:
        run_count_t7(
            bam_file=bam_file,
            ref_file=ref_file,
            barcode_file=barcode_file,
            gtf_file=gtf_file,
            t7_barcode=t7_barcode,
            blacklist_file=blacklist_file,
            whitelist_file=whitelist_file,
            k=k,
            edit_dist=edit_dist,
            stranded_edit_dist=stranded_edit_dist,
            edit_site_min_cells=edit_site_min_cells,
            nonbc_edit_dist=nonbc_edit_dist,
            ploidy=ploidy,
            copy_number_variant_file=cnv_file,
            blacklist_seqs=blacklist_seqs,
            mrna_count_mode=mrna_count_mode,
            outdir=outdir,
            edit_site_rev_comp_filt=edit_site_rev_comp_filt,
            max_gene_count_reads=None,
            uncorrected_gene_count=uncorrected_gene_count,
            constrain_allele_calls=False,
            verbosity=verbosity,
            n_cpus=n_cpus,
            chunk_size_mb=chunk_size_mb,
        )

        # Record end time
        results.end_time = time.time()

        # Collect pipeline metrics (mock data for now - would be populated by instrumented pipeline)
        # In production, count_t7.py would pass metrics back or we'd parse output files
        if outdir:
            import glob

            # Try to infer results from output files
            edit_sites_file = os.path.join(outdir, "edit_site_info.txt")
            if os.path.exists(edit_sites_file):
                results.set_output("Edit Sites", edit_sites_file)

            umi_file = os.path.join(outdir, "*umi*.mtx")
            umi_files = glob.glob(umi_file)
            if umi_files:
                results.set_output("UMI Counts", umi_files[0])

            gene_file = os.path.join(outdir, "*gene*.mtx")
            gene_files = glob.glob(gene_file)
            if gene_files:
                results.set_output("Gene Counts", gene_files[0])

        logger.info("Pipeline completed successfully")

        # Mark checkpoint as completed
        if checkpoint_manager:
            checkpoint_manager.mark_completed(success=True)
            logger.info("Checkpoint marked as completed")

        # Display results summary
        if verbosity >= 1:
            console.print("\n[green bold]✓ Pipeline completed successfully![/green bold]")
            results.display_summary(show_performance=False, show_outputs=True)
        else:
            console.print("\n[green bold]✓ Pipeline completed successfully![/green bold]")

    except Exception as e:
        # Record failure
        if checkpoint_manager:
            checkpoint_manager.mark_completed(success=False, error=str(e))

        logger.error(f"Pipeline failed: {e}", exc_info=True)
        console.print(f"\n[red bold]✗ Pipeline failed:[/red bold] {e}")
        raise typer.Exit(code=1)


def _validate_inputs(
    bam_file: str,
    ref_file: str,
    barcode_file: str,
    gtf_file: str,
    blacklist_file: Optional[str],
    whitelist_file: Optional[str],
    cnv_file: Optional[str],
    blacklist_seqs: Optional[str],
    logger,
):
    """Validate input files exist."""
    errors = []

    required_files = {
        "BAM file": bam_file,
        "Reference FASTA": ref_file,
        "Barcode file": barcode_file,
        "GTF file": gtf_file,
    }

    for name, filepath in required_files.items():
        if not os.path.exists(filepath):
            errors.append(f"  ✗ {name} not found: {filepath}")

    optional_files = {
        "Blacklist file": blacklist_file,
        "Whitelist file": whitelist_file,
        "CNV file": cnv_file,
        "Blacklist sequences": blacklist_seqs,
    }

    for name, filepath in optional_files.items():
        if filepath is not None and not os.path.exists(filepath):
            errors.append(f"  ✗ {name} not found: {filepath}")

    if errors:
        logger.error("Input validation failed")
        console.print("\n[red bold]❌ Input validation failed:[/red bold]")
        for error in errors:
            console.print(error)
        console.print("")
        raise typer.Exit(code=1)

    logger.info("Input validation passed")


def _show_dry_run(
    bam_file,
    ref_file,
    barcode_file,
    gtf_file,
    t7_barcode,
    k,
    edit_dist,
    edit_site_min_cells,
    n_cpus,
    outdir,
    blacklist_file,
    whitelist_file,
    cnv_file,
):
    """Show dry-run configuration."""
    console.print("\n[cyan bold]🔍 DRY RUN MODE - Preview Only[/cyan bold]")
    console.print("\n✓ Input validation passed")
    console.print("\n[bold]Configuration:[/bold]")
    console.print(f"  BAM file: {bam_file}")
    console.print(f"  Reference: {ref_file}")
    console.print(f"  Barcodes: {barcode_file}")
    console.print(f"  GTF: {gtf_file}")
    console.print(f"  T7 barcode: {t7_barcode}")
    console.print(f"  K-mer size: {k}")
    console.print(f"  Edit distance: {edit_dist}")
    console.print(f"  Min cells per edit: {edit_site_min_cells}")
    console.print(f"  CPUs: {n_cpus}")
    console.print(f"  Output dir: {outdir if outdir else 'Current directory'}")

    if blacklist_file:
        console.print(f"  Blacklist: {blacklist_file}")
    if whitelist_file:
        console.print(f"  Whitelist: {whitelist_file}")
    if cnv_file:
        console.print(f"  CNV file: {cnv_file}")

    console.print("\n[green]✓ Dry run complete - no files were modified[/green]")


# Export the run function directly for __main__.py
# Keep app for potential future use
def run_command(*args, **kwargs):
    """Wrapper to call run directly."""
    return run(*args, **kwargs)
