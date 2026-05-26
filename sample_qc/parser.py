"""
parser.py — Input parsers for proteomics abundance matrices.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def load_proteomics_data(
    path: str | Path,
    format_type: str | None = None,
    sep: str | None = None,
) -> pd.DataFrame:
    """
    Load proteomics input data for analysis.

    Parameters
    ----------
    path : str or Path
        Path to a file or FragPipe directory.
    format_type : {"fragpipe", "spectronaut", "auto", None}
        Optional format override.
    sep : str | None
        Optional delimiter for generic matrix files.

    Returns
    -------
    pd.DataFrame
        Abundance matrix with proteins/peptides as rows and runs as columns.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Path does not exist: {path}")

    if path.is_dir():
        combined_file = path / "combined_protein.tsv"
        if combined_file.exists():
            path = combined_file
            format_type = "fragpipe"
        else:
            raise FileNotFoundError(
                f"Directory provided but no 'combined_protein.tsv' was found in {path}"
            )

    if format_type is None or format_type == "auto":
        sample_df = _read_header(path)
        if "R.FileName" in sample_df.columns and "PG.Quantity" in sample_df.columns:
            format_type = "spectronaut"
        elif "Protein ID" in sample_df.columns and any(
            c.endswith("MaxLFQ Intensity") or (c.endswith("Intensity") and c != "Intensity")
            for c in sample_df.columns
        ):
            format_type = "fragpipe"
        else:
            format_type = "matrix"

    if format_type == "fragpipe":
        return parse_fragpipe_combined(path)
    if format_type == "spectronaut":
        return parse_spectronaut_report(path)
    if format_type == "matrix":
        return parse_abundance_matrix(path, sep=sep)

    raise ValueError(
        f"Unknown proteomics format '{format_type}'. Use 'fragpipe', 'spectronaut', or 'auto'."
    )


def _read_header(path: Path) -> pd.DataFrame:
    """Read the first lines of a file to detect column headers."""
    for sep in ["\t", ","]:
        try:
            return pd.read_csv(path, sep=sep, nrows=5, engine="python")
        except Exception:
            continue
    raise ValueError(f"Unable to read input file headers: {path}")


def parse_fragpipe_combined(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    df = pd.read_csv(path, sep="\t", na_values=["", "NaN", "0", 0], engine="python")
    df.columns = [c.strip() for c in df.columns]

    if "Protein ID" not in df.columns:
        raise ValueError(
            f"FragPipe file missing 'Protein ID' column. Headers: {list(df.columns)}"
        )

    df = df.set_index("Protein ID")
    lfq_cols = [c for c in df.columns if c.endswith("MaxLFQ Intensity")]
    intensity_cols = [
        c for c in df.columns
        if c.endswith("Intensity") and c != "Intensity" and not c.endswith("MaxLFQ Intensity")
    ]
    target_cols = lfq_cols if lfq_cols else intensity_cols
    if not target_cols:
        raise ValueError(
            "Could not find any intensity columns ending with 'Intensity' or 'MaxLFQ Intensity'."
        )

    abundance_df = df[target_cols].copy()
    suffix_to_remove = " MaxLFQ Intensity" if lfq_cols else " Intensity"
    abundance_df.columns = [c.replace(suffix_to_remove, "").strip() for c in abundance_df.columns]

    for col in abundance_df.columns:
        abundance_df[col] = pd.to_numeric(abundance_df[col], errors="coerce")

    logger.info(
        "Parsed FragPipe matrix: %d proteins across %d runs",
        len(abundance_df),
        len(abundance_df.columns),
    )
    return abundance_df


def parse_spectronaut_report(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    df = pd.read_csv(path, sep=sep, na_values=["", "NaN", "0", 0], engine="python")
    df.columns = [c.strip() for c in df.columns]

    if "R.FileName" in df.columns and "PG.ProteinGroups" in df.columns:
        quantity_col = "PG.Quantity" if "PG.Quantity" in df.columns else df.columns[-1]
        abundance_df = df.pivot_table(
            index="PG.ProteinGroups",
            columns="R.FileName",
            values=quantity_col,
            aggfunc="first",
        )
    else:
        id_candidates = ["PG.ProteinGroups", "Protein ID", "ProteinGroups", "Protein"]
        id_col = next((c for c in id_candidates if c in df.columns), None)
        if id_col is None:
            id_col = df.columns[0]
        df = df.set_index(id_col)
        abundance_df = df.select_dtypes(include="number").copy()

    for col in abundance_df.columns:
        abundance_df[col] = pd.to_numeric(abundance_df[col], errors="coerce")

    logger.info(
        "Parsed Spectronaut matrix: %d proteins across %d runs",
        len(abundance_df),
        len(abundance_df.columns),
    )
    return abundance_df


def parse_abundance_matrix(path: str | Path, sep: str | None = None) -> pd.DataFrame:
    path = Path(path)
    if sep is None:
        sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","

    df = pd.read_csv(path, sep=sep, na_values=["", "NaN"], engine="python")
    df.columns = [c.strip() for c in df.columns]

    if df.shape[1] < 2:
        raise ValueError("Proteomics matrix must contain at least one identifier column and one sample column.")

    id_col = df.columns[0]
    df = df.set_index(id_col)
    abundance_df = df.select_dtypes(include="number").copy()

    if abundance_df.empty:
        raise ValueError("No numeric abundance columns were found in the input file.")

    for col in abundance_df.columns:
        abundance_df[col] = pd.to_numeric(abundance_df[col], errors="coerce")

    logger.info(
        "Parsed generic proteomics matrix: %d proteins across %d runs",
        len(abundance_df),
        len(abundance_df.columns),
    )
    return abundance_df
