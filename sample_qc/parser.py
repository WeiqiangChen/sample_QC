"""
parser.py — Input file parsers and schema validation for sample_qc (Genomics + Proteomics).
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

# Minimum required columns in a genomics sample sheet
REQUIRED_COLUMNS: list[str] = [
    "sample_id",
    "total_reads",
    "mapped_reads",
    "gc_content",
    "coverage",
    "dup_rate",
]


def load_proteomics_data(
    path: str | Path,
    format_type: str | None = None,
    sep: str | None = None,
) -> pd.DataFrame:
    """
    Main unified data loader that automatically detects or uses specified format
    to load a numeric abundance matrix (proteins/peptides as rows, samples as columns)
    or standard sample sheets.

    Parameters
    ----------
    path : str or Path
        Path to file or output directory.
    format_type : {"fragpipe", "spectronaut", "genomics", None}
        Optional format override.
    sep : str
        Optional column delimiter for genomic/custom sample sheets.

    Returns
    -------
    pd.DataFrame
        If genomics: standard sample sheet DataFrame.
        If proteomics: wide abundance matrix (index=proteins, columns=samples).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Path does not exist: {path}")

    # Resolve folder structure for FragPipe
    if path.is_dir():
        combined_file = path / "combined_protein.tsv"
        if combined_file.exists():
            path = combined_file
            format_type = "fragpipe"
        else:
            raise FileNotFoundError(
                f"Directory provided but no 'combined_protein.tsv' found inside: {path}"
            )

    # Inferred format detection
    if format_type is None:
        # Read small sample of headers
        try:
            sample_df = pd.read_csv(path, sep="\t", nrows=5, engine="python")
        except Exception:
            try:
                sample_df = pd.read_csv(path, sep=",", nrows=5, engine="python")
            except Exception as e:
                raise ValueError(f"Failed to read file headers: {e}")

        # Check for Spectronaut long format
        if "R.FileName" in sample_df.columns and "PG.Quantity" in sample_df.columns:
            format_type = "spectronaut"
        # Check for FragPipe combined protein
        elif "Protein ID" in sample_df.columns and any(
            "Intensity" in c or "Spectral Count" in c for c in sample_df.columns
        ):
            format_type = "fragpipe"
        # Else check if fits basic genomics layout
        elif "sample_id" in sample_df.columns:
            format_type = "genomics"
        else:
            # Fallback default: try to parse as genomics
            format_type = "genomics"

    logger.info("Detected format: %s for input: %s", format_type.upper(), path)

    if format_type == "fragpipe":
        return parse_fragpipe_combined(path)
    elif format_type == "spectronaut":
        return parse_spectronaut_report(path)
    else:
        return parse_sample_sheet(path, sep=sep)


def parse_fragpipe_combined(path: str | Path) -> pd.DataFrame:
    """
    Parses a FragPipe DDA `combined_protein.tsv` output file.

    Extracts MaxLFQ Intensity columns (or general Intensity columns if LFQ is not present)
    per sample, using Protein ID as the index.

    Returns
    -------
    pd.DataFrame
        Abundance matrix where index=Protein ID and columns=Sample Names.
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t", na_values=["", "NaN", "0", 0], engine="python")
    df.columns = [c.strip() for c in df.columns]

    if "Protein ID" not in df.columns:
        raise ValueError(f"FragPipe file missing 'Protein ID' column. Headers: {list(df.columns)}")

    # Set index to Protein ID
    df = df.set_index("Protein ID")

    # Find sample intensity columns. Prefer "MaxLFQ Intensity" columns.
    lfq_cols = [c for c in df.columns if c.endswith("MaxLFQ Intensity")]
    intensity_cols = [
        c for c in df.columns 
        if c.endswith("Intensity") and c != "Intensity" and not c.endswith("MaxLFQ Intensity")
    ]

    target_cols = lfq_cols if lfq_cols else intensity_cols
    if not target_cols:
        raise ValueError(
            "Could not find any sample-specific quantitative intensity columns "
            f"ending with 'Intensity' or 'MaxLFQ Intensity' in {path.name}"
        )

    # Slice the quantitative columns
    abundance_df = df[target_cols].copy()

    # Clean sample names (e.g. remove " MaxLFQ Intensity" suffix)
    suffix_to_remove = " MaxLFQ Intensity" if lfq_cols else " Intensity"
    abundance_df.columns = [c.replace(suffix_to_remove, "").strip() for c in abundance_df.columns]

    # Convert all columns to numeric, coercion
    for col in abundance_df.columns:
        abundance_df[col] = pd.to_numeric(abundance_df[col], errors="coerce")

    logger.info(
        "Parsed FragPipe matrix: %d proteins across %d samples",
        len(abundance_df),
        len(abundance_df.columns),
    )
    return abundance_df


def parse_spectronaut_report(path: str | Path) -> pd.DataFrame:
    """
    Parses a Spectronaut DIA exported report file.

    Handles both wide format (columns are sample names) and standard pivoted long format
    (containing columns: `R.FileName`, `PG.ProteinGroups`, `PG.Quantity`).

    Returns
    -------
    pd.DataFrame
        Abundance matrix where index=PG.ProteinGroups and columns=Sample Names.
    """
    path = Path(path)
    sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    df = pd.read_csv(path, sep=sep, na_values=["", "NaN", "0", 0], engine="python")
    df.columns = [c.strip() for c in df.columns]

    # Detect long format
    if "R.FileName" in df.columns and "PG.ProteinGroups" in df.columns:
        quantity_col = "PG.Quantity" if "PG.Quantity" in df.columns else df.columns[-1]
        logger.info(
            "Pivoting long-format Spectronaut report using %s as intensity metric",
            quantity_col,
        )
        
        # Pivot the table to create a wide abundance matrix
        abundance_df = df.pivot_table(
            index="PG.ProteinGroups",
            columns="R.FileName",
            values=quantity_col,
            aggfunc="first",
        )
    else:
        # If it's already a wide format, we identify the primary index column (e.g. protein group)
        id_candidates = ["PG.ProteinGroups", "Protein ID", "ProteinGroups", "Protein"]
        id_col = None
        for candidate in id_candidates:
            if candidate in df.columns:
                id_col = candidate
                break
        
        if id_col is None:
            id_col = df.columns[0]
            logger.warning("Could not find protein group column. Using first column '%s' as index", id_col)

        df = df.set_index(id_col)
        # Keep numeric columns
        abundance_df = df.select_dtypes(include="number").copy()

    # Convert all columns to float
    for col in abundance_df.columns:
        abundance_df[col] = pd.to_numeric(abundance_df[col], errors="coerce")

    logger.info(
        "Parsed Spectronaut matrix: %d proteins across %d samples",
        len(abundance_df),
        len(abundance_df.columns),
    )
    return abundance_df


def parse_sample_sheet(
    path: str | Path,
    sep: str | None = None,
    required_cols: list[str] | None = None,
) -> pd.DataFrame:
    """
    Load a standard TSV or CSV sample sheet into a DataFrame (Genomics).
    """
    path = Path(path)
    if sep is None:
        sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","

    na_values = ["NA", "N/A", "nan", "NaN", ".", "", "null", "NULL"]
    df = pd.read_csv(path, sep=sep, na_values=na_values, engine="python")
    df.columns = [c.strip() for c in df.columns]

    validate_schema(df, required_cols=required_cols)
    _coerce_numeric_columns(df)

    return df


def validate_schema(
    df: pd.DataFrame,
    required_cols: list[str] | None = None,
) -> None:
    """
    Validate that all required columns are present in ``df`` (Genomics).
    """
    if required_cols is None:
        required_cols = REQUIRED_COLUMNS

    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"Sample sheet is missing required column(s): {missing}. "
            f"Found columns: {list(df.columns)}"
        )


def _coerce_numeric_columns(df: pd.DataFrame) -> None:
    """Convert columns that look numeric (but were read as object) to float."""
    for col in df.columns:
        if df[col].dtype == object:
            try:
                df[col] = pd.to_numeric(df[col], errors="ignore")
            except Exception:
                pass

