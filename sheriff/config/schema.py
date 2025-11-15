"""Pydantic schema for Sheriff configuration."""

from pathlib import Path
from typing import Optional, Literal
from pydantic import BaseModel, Field, field_validator
import yaml


class InputConfig(BaseModel):
    """Input file configuration."""

    bam_file: str = Field(..., description="BAM file path")
    ref_file: str = Field(..., description="Reference genome FASTA file")
    barcode_file: str = Field(..., description="Barcode whitelist file")
    gtf_file: str = Field(..., description="GTF annotation file")

    @field_validator("bam_file", "ref_file", "barcode_file", "gtf_file")
    @classmethod
    def validate_file_exists(cls, v: str) -> str:
        """Validate that file paths exist (optional, can be disabled)."""
        # Validation happens at runtime, not at schema level
        return v


class FiltersConfig(BaseModel):
    """Optional filter file configuration."""

    blacklist_file: Optional[str] = Field(None, description="Blacklist regions BED file")
    whitelist_file: Optional[str] = Field(None, description="Whitelist regions BED file")
    blacklist_seqs: Optional[str] = Field(None, description="Blacklist sequences file")
    cnv_file: Optional[str] = Field(None, description="Copy number variation file")


class ParametersConfig(BaseModel):
    """Pipeline parameters configuration."""

    t7_barcode: str = Field("GGGAGAGTAT", description="T7 barcode sequence")
    k: int = Field(6, ge=3, le=15, description="K-mer size for barcode matching")
    edit_dist: int = Field(140, ge=0, description="Edit site grouping distance")
    stranded_edit_dist: int = Field(15, ge=0, description="Max distance between forward/reverse edits")
    edit_site_min_cells: int = Field(3, ge=1, description="Minimum cells per edit site")
    nonbc_edit_dist: int = Field(1000, ge=0, description="Distance for non-barcoded read assignment")
    ploidy: int = Field(2, ge=1, description="Genome ploidy")
    mrna_count_mode: Literal["all", "polyT"] = Field("all", description="mRNA counting mode")


class ProcessingConfig(BaseModel):
    """Processing options configuration."""

    n_cpus: int = Field(1, ge=1, description="Number of CPUs to use")
    chunk_size_mb: int = Field(250, ge=1, description="Chunk size in MB for processing")
    verbosity: int = Field(1, ge=0, le=2, description="Verbosity level (0=errors, 1=info, 2=debug)")


class OutputConfig(BaseModel):
    """Output options configuration."""

    outdir: Optional[str] = Field(None, description="Output directory")
    edit_site_rev_comp_filt: bool = Field(True, description="Require bidirectional insert evidence")
    uncorrected_gene_count: bool = Field(False, description="Output uncorrected gene counts")


class ResumeConfig(BaseModel):
    """Resume functionality configuration."""

    enabled: bool = Field(False, description="Enable automatic checkpointing")
    checkpoint_dir: str = Field(".sheriff_checkpoints", description="Checkpoint directory")
    checkpoint_interval: int = Field(5, ge=1, description="Checkpoint interval in minutes")


class LoggingConfig(BaseModel):
    """Logging configuration."""

    file: Optional[str] = Field(None, description="Log file path")
    level: Literal["DEBUG", "INFO", "WARNING", "ERROR"] = Field("INFO", description="Log level")
    format: Literal["simple", "detailed", "json"] = Field("detailed", description="Log format")


class ProfilingConfig(BaseModel):
    """Profiling configuration."""

    enabled: bool = Field(False, description="Enable profiling")
    output: str = Field("sheriff_profile.json", description="Profiling output file")
    memory_tracking: bool = Field(True, description="Enable memory tracking")


class SheriffConfig(BaseModel):
    """Complete Sheriff configuration schema."""

    version: str = Field("1.2.0", description="Config file version")
    input: InputConfig
    filters: FiltersConfig = Field(default_factory=FiltersConfig)
    parameters: ParametersConfig = Field(default_factory=ParametersConfig)
    processing: ProcessingConfig = Field(default_factory=ProcessingConfig)
    output: OutputConfig = Field(default_factory=OutputConfig)
    resume: ResumeConfig = Field(default_factory=ResumeConfig)
    logging: LoggingConfig = Field(default_factory=LoggingConfig)
    profiling: ProfilingConfig = Field(default_factory=ProfilingConfig)

    model_config = {"extra": "forbid"}  # Disallow extra fields

    def to_cli_kwargs(self) -> dict:
        """Convert config to CLI keyword arguments."""
        return {
            # Input files
            "bam_file": self.input.bam_file,
            "ref_file": self.input.ref_file,
            "barcode_file": self.input.barcode_file,
            "gtf_file": self.input.gtf_file,
            # Filters
            "blacklist_file": self.filters.blacklist_file,
            "whitelist_file": self.filters.whitelist_file,
            "blacklist_seqs": self.filters.blacklist_seqs,
            "cnv_file": self.filters.cnv_file,
            # Parameters
            "t7_barcode": self.parameters.t7_barcode,
            "k": self.parameters.k,
            "edit_dist": self.parameters.edit_dist,
            "stranded_edit_dist": self.parameters.stranded_edit_dist,
            "edit_site_min_cells": self.parameters.edit_site_min_cells,
            "nonbc_edit_dist": self.parameters.nonbc_edit_dist,
            "ploidy": self.parameters.ploidy,
            "mrna_count_mode": self.parameters.mrna_count_mode,
            # Processing
            "n_cpus": self.processing.n_cpus,
            "chunk_size_mb": self.processing.chunk_size_mb,
            "verbosity": self.processing.verbosity,
            # Output
            "outdir": self.output.outdir,
            "edit_site_rev_comp_filt": self.output.edit_site_rev_comp_filt,
            "uncorrected_gene_count": self.output.uncorrected_gene_count,
        }


def load_config(config_path: str | Path) -> SheriffConfig:
    """Load and validate Sheriff configuration from YAML file.

    Args:
        config_path: Path to YAML configuration file

    Returns:
        Validated SheriffConfig object

    Raises:
        FileNotFoundError: If config file doesn't exist
        ValueError: If config is invalid
    """
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as f:
        config_dict = yaml.safe_load(f)

    try:
        return SheriffConfig(**config_dict)
    except Exception as e:
        raise ValueError(f"Invalid config file: {e}")


def merge_config(config: SheriffConfig, cli_kwargs: dict) -> dict:
    """Merge configuration with CLI arguments.

    CLI arguments take precedence over config file values.

    Args:
        config: SheriffConfig object from file
        cli_kwargs: Dictionary of CLI arguments

    Returns:
        Merged configuration dictionary
    """
    # Start with config values
    merged = config.to_cli_kwargs()

    # Override with CLI values (if not None)
    for key, value in cli_kwargs.items():
        if value is not None:
            merged[key] = value

    return merged
