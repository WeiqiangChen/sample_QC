"""
sample_qc — Quality Control for proteomics sample data.
"""

__version__ = "0.1.0"
__author__ = "WeiQiang"

from sample_qc.metrics import run_all_metrics
from sample_qc.parser import load_proteomics_data
from sample_qc.report import generate_report
from sample_qc.cli import main

__all__ = [
    "main",
    "load_proteomics_data",
    "run_all_metrics",
    "generate_report",
]


def run_qc(
    input_path: str,
    sep: str | None = None,
    format_type: str | None = None,
) -> dict:
    """
    High-level convenience function for proteomics QC.

    Parameters
    ----------
    input_path : str
        Path to a proteomics input file or FragPipe project directory.
    sep : str | None
        Optional delimiter for generic abundance matrices.
    format_type : {"fragpipe", "spectronaut", "auto", None}
        Optional format override.

    Returns
    -------
    dict
        QC results dictionary ready for report generation.
    """
    df = load_proteomics_data(input_path, format_type=format_type, sep=sep)
    return run_all_metrics(df)
