from pathlib import Path
from typing import List, Optional
from typing_extensions import Annotated

import typer
import sys
import os

# Local Imports
from .count_t7 import run_count_t7
from . import __version__
from .cli.run import run as run_pipeline
from .cli.config_cmd import config_command
from .cli.validate import validate_command

# Create main app with subcommands
app = typer.Typer(
    pretty_exceptions_short=False, help="Sheriff - CRISPR/Cas9 Edit Site Detection in Single Cells", no_args_is_help=True
)


def version_callback(value: bool):
    """Print version and exit."""
    if value:
        typer.echo(f"Sheriff version {__version__}")
        raise typer.Exit()


@app.callback()
def common_options(
    version: Annotated[
        Optional[bool],
        typer.Option("--version", "-V", callback=version_callback, is_eager=True, help="Show version and exit"),
    ] = None,
):
    """Sheriff - CRISPR/Cas9 Edit Site Detection in Single Cells"""
    pass


# Add config subcommand (has multiple sub-commands: generate, validate, show)
app.add_typer(config_command, name="config", help="Configuration file utilities")

# Add run command (direct command, not a group)
app.command(name="run", help="Run the Sheriff CRISPR edit calling pipeline")(run_pipeline)

# Add validate command (direct command, not a group)
app.command(name="validate", help="Validate inputs without running pipeline")(validate_command)


# BACKWARD COMPATIBILITY: Keep get-t7-edits as an alias to run
@app.command(name="get-t7-edits", hidden=True)
def get_t7_edits(
    bam_file: Annotated[str, typer.Argument(help="BAM file")],
    ref_file: Annotated[str, typer.Argument(help="Fasta containing ref genome")],
    barcode_file: Annotated[str, typer.Argument(help="Text file containing whitelisted barcode per line")],
    gtf_file: Annotated[str, typer.Argument(help="GTF file containing relevant gene data")],
    t7_barcode: Annotated[
        Optional[str], typer.Option("--t7", "--t7_barcode", help=("Target/query barcode sequence to denote t7 reads."))
    ] = "GGGAGAGTAT",
    blacklist_file: Annotated[
        Optional[str],
        typer.Option(
            "--blacklist",
            "--blacklist_file",
            help=(
                "Bed file that specifies the location of blacklist regions, these generate a lot of endogenous t7 reads "
                "that can lead to slow processing time and false-positive edit-site calling."
            ),
        ),
    ] = None,
    whitelist_file: Annotated[
        Optional[str],
        typer.Option(
            "--whitelist",
            "--whitelist_file",
            help=(
                "Bed file that specifies the location of whitelist regions, which are known edit sites and so will call "
                "any barcoded reads implying an edit site intersecting these regions as canonical edit sites."
            ),
        ),
    ] = None,
    k: Annotated[
        Optional[int],
        typer.Option(
            "-k",
            "--kmer",
            "--kmer_size",
            help=("Size of kmers used to pattern match read barcodes to the t7 barcode."),
        ),
    ] = 6,
    edit_dist: Annotated[
        Optional[int],
        typer.Option("--edit_dist", "--edist", "--dist", help=("+/- distance from edit site to be grouped as same edit.")),
    ] = 140,
    edit_site_rev_comp_filt: Annotated[
        Optional[bool],
        typer.Option(
            "--bidirectional_inserts/--no-bidirectional_inserts",
            help=(
                "Candidate edit site must have evidence of bi-directional donor insertion to be called as a canonical edit site."
                "Highly recommended criteria. If turned off, it also turns off --stranded_edit_dist criteria, but this "
                "information about the edit sites are still recorded in the output via the 'stranded_edit_dist' column"
                "in edit_site_info.txt"
            ),
        ),
    ] = True,
    stranded_edit_dist: Annotated[
        Optional[int],
        typer.Option(
            "--stranded_edit_dist",
            help=(
                "Maximum allowed distance between the nearest forward and reverse edit sites at a given canonical edit site to qualify as real edit."
            ),
        ),
    ] = 15,
    edit_site_min_cells: Annotated[
        Optional[int], typer.Option("--edit_site_min_cells", help=("Minimum cells in edit site to be considered true edit."))
    ] = 3,
    nonbc_edit_dist: Annotated[
        Optional[int],
        typer.Option(
            "--nonbc_edit_dist",
            "--nonbc_edist",
            "--nonbc_dist",
            "--nonbc",
            help=("+/- distance from edit to mop up the non-barcoded reads."),
        ),
    ] = 1000,
    ploidy: Annotated[Optional[int], typer.Option("--ploidy", help=("Ploidy/Number of chromosomes in the genome."))] = 2,
    cnv_file: Annotated[
        Optional[str],
        typer.Option(
            "--cnv",
            "--cnv_file",
            "--copy_number_variant_file",
            help=("A bedGraph file that specifies copy-number-variation sites, " "that deviate from the ploidy number."),
        ),
    ] = None,
    blacklist_seqs: Annotated[
        Optional[str],
        typer.Option(
            "--blacklist_seqs",
            help=(
                "Text file of sequences, with a new sequence on each line, that may be present in read soft-clip sequences"
                "can confound t7 barcoded read calls. Currently only the TSO, which is a common left-over artifact."
            ),
        ),
    ] = None,
    mrna_count_mode: Annotated[
        Optional[str],
        typer.Option(
            "--mrna_count_mode",
            help=(
                "Mode for quantifying gene expression,"
                "'all' is to count all reads associated with a gene, 'polyT' is to only count polyT reads, indicating mature mRNA transcripts."
            ),
        ),
    ] = "all",
    uncorrected_gene_count: Annotated[
        Optional[bool],
        typer.Option(
            "--uncorrected_count/--no-uncorrected_count",
            help=("Whether to also output the gene counts WITHOUT removing reads around inferred edit sites."),
        ),
    ] = False,
    outdir: Annotated[
        Optional[str],
        typer.Option(
            "-o", "--out", "--outdir", "--out_dir", help=("Write output files to this location. " "Defaults to Current Working Directory")
        ),
    ] = None,
    verbosity: Annotated[
        Optional[int],
        typer.Option(
            "--v",
            "-v",
            "--verbosity",
            "-verbosity",
            help=("Verbosity levels. 0 errors only, 1 prints processing progress, 2 prints debugging information."),
        ),
    ] = 1,
    n_cpus: Annotated[
        Optional[int], typer.Option("--cpu", "-cpu", help=("Number of CPUs to use for processing, necessary to increase this for fast UMI counting."))
    ] = 1,
    chunk_size_mb: Annotated[
        Optional[int],
        typer.Option(
            "--chunk",
            "-chunk",
            help=(
                "Number of mega-bases to process at a time for gene UMI counting. Set this lower if get memory"
                "issues, currently set to > hg38 chr1 size (249 Mb), so is parallelized by chromosome."
                "Trade-off is can cause minor double-counting of UMIs if a genome chunk intersects a gene."
                "Very small / negligible difference in counts."
            ),
        ),
    ] = 250,
    dry_run: Annotated[
        Optional[bool],
        typer.Option(
            "--dry-run", "--dry_run", help=("Preview what will be executed without running the pipeline. " "Validates inputs and shows configuration.")
        ),
    ] = False,
):
    """[DEPRECATED] Use 'sheriff run' instead.

    This command is kept for backward compatibility.
    """
    from .logging_config import quick_setup
    from rich.console import Console

    console = Console()

    # Setup basic logging
    logger = quick_setup(verbosity)
    logger.warning("'get-t7-edits' is deprecated. Use 'sheriff run' instead.")

    # Validate input files before processing
    def validate_input_files(
        bam_file: str,
        ref_file: str,
        barcode_file: str,
        gtf_file: str,
        blacklist_file: Optional[str] = None,
        whitelist_file: Optional[str] = None,
        cnv_file: Optional[str] = None,
        blacklist_seqs: Optional[str] = None,
    ) -> None:
        """Validate that required input files exist before processing."""
        errors = []

        # Check required files
        required_files = {
            "BAM file": bam_file,
            "Reference FASTA": ref_file,
            "Barcode file": barcode_file,
            "GTF file": gtf_file,
        }

        for name, filepath in required_files.items():
            if not os.path.exists(filepath):
                errors.append(f"  ✗ {name} not found: {filepath}")

        # Check optional files if provided
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
            console.print("\n[red bold]❌ Input validation failed:[/red bold]")
            for error in errors:
                console.print(error)
            console.print("")
            raise typer.Exit(code=1)

    validate_input_files(
        bam_file=bam_file,
        ref_file=ref_file,
        barcode_file=barcode_file,
        gtf_file=gtf_file,
        blacklist_file=blacklist_file,
        whitelist_file=whitelist_file,
        cnv_file=cnv_file,
        blacklist_seqs=blacklist_seqs,
    )

    # Dry-run mode: show configuration and exit
    if dry_run:
        console.print("\n[cyan bold]🔍 DRY RUN MODE - Preview Only[/cyan bold]")
        console.print("\n✓ Input validation passed")
        console.print("\nConfiguration:")
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

        # Optional files
        if blacklist_file:
            console.print(f"  Blacklist: {blacklist_file}")
        if whitelist_file:
            console.print(f"  Whitelist: {whitelist_file}")
        if cnv_file:
            console.print(f"  CNV file: {cnv_file}")

        console.print("\n[green]✓ Dry run complete - no files were modified[/green]")
        raise typer.Exit()

    # Show processing start message
    if verbosity >= 1:
        console.print(f"\n[green bold]🚀 Sheriff v{__version__} - Starting analysis[/green bold]")
        console.print(f"  BAM: {bam_file}")
        console.print(f"  Output: {outdir if outdir else 'Current directory'}\n")

    # Run
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


def main():
    root_dir = Path(__file__).parent
    sys.path.append(str(root_dir))
    app()


if __name__ == "__main__":
    main()
