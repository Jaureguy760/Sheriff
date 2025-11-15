"""Configuration management subcommand."""

from typing import Optional
from typing_extensions import Annotated
import typer
from rich.console import Console

from ..config.loader import ConfigLoader

console = Console()
app = typer.Typer(help="Configuration file utilities")


@app.command("generate")
def generate_template(
    output: Annotated[str, typer.Option("--output", "-o", help="Output file path")] = "sheriff-config.yaml",
    minimal: Annotated[bool, typer.Option("--minimal", help="Generate minimal config (required fields only)")] = False,
):
    """Generate a configuration file template.

    Examples:
        sheriff config generate
        sheriff config generate --output my-config.yaml
        sheriff config generate --minimal
    """
    ConfigLoader.generate_template(output, minimal=minimal)


@app.command("validate")
def validate_config(
    config_file: Annotated[str, typer.Argument(help="Config file to validate")],
):
    """Validate a configuration file.

    Examples:
        sheriff config validate sheriff-config.yaml
    """
    success = ConfigLoader.validate_file(config_file)
    if not success:
        raise typer.Exit(code=1)


@app.command("show")
def show_config(
    config_file: Annotated[str, typer.Argument(help="Config file to display")],
):
    """Display configuration file with syntax highlighting.

    Examples:
        sheriff config show sheriff-config.yaml
    """
    try:
        loader = ConfigLoader(config_file)
        loader.display()
    except (FileNotFoundError, ValueError) as e:
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(code=1)


# Export the app as config_command
config_command = app
