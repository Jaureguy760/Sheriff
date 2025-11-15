"""Configuration loader with CLI integration."""

from pathlib import Path
from typing import Optional, Dict, Any
import yaml
from rich.console import Console
from rich.syntax import Syntax

from .schema import SheriffConfig, load_config, merge_config

console = Console()


class ConfigLoader:
    """Handles loading and merging of configuration from files and CLI."""

    def __init__(self, config_path: Optional[str] = None):
        """Initialize config loader.

        Args:
            config_path: Optional path to YAML config file
        """
        self.config_path = config_path
        self.config: Optional[SheriffConfig] = None

        if config_path:
            self.load()

    def load(self) -> SheriffConfig:
        """Load configuration from file.

        Returns:
            Loaded and validated SheriffConfig

        Raises:
            FileNotFoundError: If config file doesn't exist
            ValueError: If config is invalid
        """
        if not self.config_path:
            raise ValueError("No config path provided")

        self.config = load_config(self.config_path)
        return self.config

    def merge_with_cli(self, cli_kwargs: Dict[str, Any]) -> Dict[str, Any]:
        """Merge config file with CLI arguments.

        Args:
            cli_kwargs: Dictionary of CLI arguments

        Returns:
            Merged configuration
        """
        if self.config is None:
            # No config file, return CLI kwargs
            return cli_kwargs

        return merge_config(self.config, cli_kwargs)

    def display(self):
        """Display loaded configuration with syntax highlighting."""
        if self.config is None:
            console.print("[yellow]No configuration loaded[/yellow]")
            return

        # Convert to YAML for display
        config_dict = self.config.model_dump()
        yaml_str = yaml.dump(config_dict, default_flow_style=False, sort_keys=False)

        syntax = Syntax(yaml_str, "yaml", theme="monokai", line_numbers=True)
        console.print("\n[bold cyan]Loaded Configuration:[/bold cyan]")
        console.print(syntax)

    @staticmethod
    def generate_template(output_path: str = "sheriff-config.yaml", minimal: bool = False):
        """Generate a template configuration file.

        Args:
            output_path: Where to save the template
            minimal: If True, generate minimal config with only required fields
        """
        if minimal:
            template = {
                "version": "1.2.0",
                "input": {
                    "bam_file": "path/to/sample.bam",
                    "ref_file": "path/to/genome.fa",
                    "barcode_file": "path/to/barcodes.txt",
                    "gtf_file": "path/to/genes.gtf",
                },
            }
        else:
            template = {
                "version": "1.2.0",
                "input": {
                    "bam_file": "path/to/sample.bam",
                    "ref_file": "path/to/genome.fa",
                    "barcode_file": "path/to/barcodes.txt",
                    "gtf_file": "path/to/genes.gtf",
                },
                "filters": {
                    "blacklist_file": None,
                    "whitelist_file": None,
                    "blacklist_seqs": None,
                    "cnv_file": None,
                },
                "parameters": {
                    "t7_barcode": "GGGAGAGTAT",
                    "k": 6,
                    "edit_dist": 140,
                    "stranded_edit_dist": 15,
                    "edit_site_min_cells": 3,
                    "nonbc_edit_dist": 1000,
                    "ploidy": 2,
                    "mrna_count_mode": "all",
                },
                "processing": {"n_cpus": 1, "chunk_size_mb": 250, "verbosity": 1},
                "output": {
                    "outdir": None,
                    "edit_site_rev_comp_filt": True,
                    "uncorrected_gene_count": False,
                },
                "resume": {
                    "enabled": False,
                    "checkpoint_dir": ".sheriff_checkpoints",
                    "checkpoint_interval": 5,
                },
                "logging": {"file": None, "level": "INFO", "format": "detailed"},
                "profiling": {"enabled": False, "output": "sheriff_profile.json", "memory_tracking": True},
            }

        output_path = Path(output_path)

        # Add helpful comments to YAML
        yaml_str = yaml.dump(template, default_flow_style=False, sort_keys=False)

        # Add header comment
        header = """# Sheriff Configuration File
# Generated template - customize for your analysis
#
# Documentation: https://github.com/BradBalderson/Sheriff

"""
        yaml_str = header + yaml_str

        output_path.write_text(yaml_str)
        console.print(f"[green]✓ Generated config template:[/green] {output_path}")

    @staticmethod
    def validate_file(config_path: str):
        """Validate a configuration file.

        Args:
            config_path: Path to config file to validate
        """
        try:
            config = load_config(config_path)
            console.print(f"[green]✓ Configuration valid:[/green] {config_path}")
            console.print(f"\n[dim]Sheriff version: {config.version}[/dim]")
            return True
        except FileNotFoundError as e:
            console.print(f"[red]✗ File not found:[/red] {e}")
            return False
        except ValueError as e:
            console.print(f"[red]✗ Invalid configuration:[/red]")
            console.print(f"  {e}")
            return False
