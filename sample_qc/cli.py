"""
cli.py — Click CLI commands for Sample QC pipeline.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import click
from rich.console import Console

from sample_qc import __version__
from sample_qc.metrics import run_all_metrics
from sample_qc.parser import load_proteomics_data
from sample_qc.report import generate_report

console = Console()


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def main() -> None:
    """🧬 Sample QC: Quality control metrics and reporting for genomic and proteomic samples.

    Calculate metrics, flag outliers, detect duplicates, and generate HTML, JSON,
    or Plain Text reports easily from your genomic sample sheets or proteomic abundance matrices.
    """
    pass


@main.command("run")
@click.argument(
    "input_path",
    type=click.Path(exists=True, dir_okay=True, path_type=Path),
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(file_okay=False, writable=True, path_type=Path),
    default=Path("./qc_results"),
    help="Directory to save output report files.",
)
@click.option(
    "-f",
    "--format",
    "fmt",
    type=click.Choice(["json", "html", "text", "all"]),
    default="all",
    help="Output report formats to generate.",
)
@click.option(
    "-i",
    "--input-format",
    "input_format",
    type=click.Choice(["fragpipe", "spectronaut", "genomics", "auto"]),
    default="auto",
    help="Input dataset format type.",
)
@click.option(
    "-s",
    "--sep",
    type=str,
    default=None,
    help="Delimiter characters. Auto-inferred if omitted.",
)
def run_pipeline(
    input_path: Path,
    output_dir: Path,
    fmt: str,
    input_format: str,
    sep: str | None,
) -> None:
    """Run full QC analysis on an input dataset (Genomics/FragPipe/Spectronaut)."""
    console.print(f"[bold indigo]🧬 Sample QC Pipeline v{__version__}[/bold indigo]")
    console.print(f"Loading input: [cyan]{input_path}[/cyan] (Format: [yellow]{input_format.upper()}[/yellow])")

    try:
        # 1. Parse sheet/matrix
        format_type = None if input_format == "auto" else input_format
        df = load_proteomics_data(input_path, format_type=format_type, sep=sep)
        
        # Check whether genomics or proteomics matrix was loaded
        is_proteomics = "total_reads" not in df.columns
        if is_proteomics:
            console.print(f"  ✓ Successfully parsed Proteomics matrix: [green]{len(df)}[/green] proteins across [green]{len(df.columns)}[/green] samples.")
        else:
            console.print(f"  ✓ Successfully parsed Genomics sheet: [green]{len(df)}[/green] samples.")

        # 2. Run analysis metrics
        console.print("Running descriptive statistics, PCA clustering, and outlier detection...")
        results = run_all_metrics(df)
        summary = results["summary"]

        console.print(
            f"  ✓ Process complete: [green]{summary['n_pass']} Passed[/green], "
            f"[red]{summary['n_fail']} Failed[/red] QC checks. "
            f"({summary['pass_rate'] * 100:.1f}% rate)"
        )

        # 3. Save report files
        output_dir.mkdir(parents=True, exist_ok=True)
        formats_to_gen = ["json", "html", "text"] if fmt == "all" else [fmt]

        for current_fmt in formats_to_gen:
            filename = f"qc_report.{current_fmt}"
            out_path = output_dir / filename
            generate_report(results, output_path=out_path, fmt=current_fmt)
            console.print(f"  ✓ Saved [bold]{current_fmt.upper()}[/bold] report to [yellow]{out_path}[/yellow]")

        # Print standard text report directly to console stdout for immediate inspection
        console.print("\n[bold]Summary Console Table Output:[/bold]")
        text_summary = generate_report(results, fmt="text")
        print(text_summary)

    except Exception as e:
        console.print(f"[bold red]Pipeline Error:[/bold red] {e}", err=True)
        sys.exit(1)


@main.command("report")
@click.argument(
    "json_results",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=False, writable=True, path_type=Path),
    required=True,
    help="Output file path (e.g. report.html).",
)
@click.option(
    "-f",
    "--format",
    "fmt",
    type=click.Choice(["html", "text"]),
    default="html",
    help="Output report format to render.",
)
def render_saved(json_results: Path, output: Path, fmt: str) -> None:
    """Render HTML or Text reports from previously saved QC JSON files."""
    console.print(f"[bold indigo]🧬 Sample QC Renderer[/bold indigo]")
    console.print(f"Loading results: [cyan]{json_results}[/cyan]")

    try:
        with open(json_results, "r", encoding="utf-8") as f:
            results = json.load(f)

        generate_report(results, output_path=output, fmt=fmt)
        console.print(f"  ✓ Generated [bold]{fmt.upper()}[/bold] report at [green]{output}[/green]")

    except Exception as e:
        console.print(f"[bold red]Render Error:[/bold red] {e}", err=True)
        sys.exit(1)


@main.command("version")
def print_version() -> None:
    """Display software version and packaging information."""
    console.print(f"🧬 sample-qc CLI tool version: [bold green]{__version__}[/bold green]")


if __name__ == "__main__":
    main()
