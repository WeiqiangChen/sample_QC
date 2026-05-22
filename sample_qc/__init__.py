"""
sample_qc — Quality Control for genomic/sequencing sample data.
"""

__version__ = "0.1.0"
__author__ = "WeiQiang"

from sample_qc.metrics import (
    compute_basic_stats,
    compute_missing_rate,
    compute_outliers,
    compute_duplication_rate,
    compute_gc_stats,
    compute_coverage_stats,
    run_all_metrics,
)
from sample_qc.parser import parse_sample_sheet, validate_schema, load_proteomics_data
from sample_qc.report import generate_report

__all__ = [
    "compute_basic_stats",
    "compute_missing_rate",
    "compute_outliers",
    "compute_duplication_rate",
    "compute_gc_stats",
    "compute_coverage_stats",
    "run_all_metrics",
    "parse_sample_sheet",
    "validate_schema",
    "load_proteomics_data",
    "generate_report",
    "run_qc",
]


def run_qc(
    input_path: str,
    sep: str | None = None,
    format_type: str | None = None,
) -> dict:
    """
    High-level convenience function: parse genomic/proteomics input and run all QC metrics.

    Parameters
    ----------
    input_path : str
        Path to file or output directory.
    sep : str
        Column delimiter for sample sheet (if applicable).
    format_type : {"fragpipe", "spectronaut", "genomics", None}
        Optional format override.

    Returns
    -------
    dict
        QC results dictionary ready to be passed to ``generate_report``.
    """
    df = load_proteomics_data(input_path, format_type=format_type, sep=sep)
    return run_all_metrics(df)
