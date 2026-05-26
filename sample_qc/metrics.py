"""
metrics.py — Core QC metric computations for sample data (Genomics & Proteomics).

All functions accept a pandas DataFrame and return structured dictionaries
suitable for JSON serialisation and report generation.
"""

from __future__ import annotations

import os
import re
import logging
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)

# Genomics Columns expected in a standard sample sheet
STANDARD_NUMERIC_COLS = [
    "total_reads",
    "mapped_reads",
    "gc_content",
    "coverage",
    "dup_rate",
]

# Thresholds for Genomics QC flags
GC_LOW = 0.35
GC_HIGH = 0.65
MIN_MAPPING_RATE = 0.70
MAX_DUP_RATE = 0.30
MIN_COVERAGE = 10.0
OUTLIER_ZSCORE = 3.0
OUTLIER_IQR_FACTOR = 1.5


# ---------------------------------------------------------------------------
# Proteomics Metrics Implementation
# ---------------------------------------------------------------------------

def generate_qc_metrics(job_dir: str | os.PathLike, n_ms2: int | None = None, tic_area: float | None = None) -> dict[str, Any]:
    metrics = {
        "is_proteomics": True,
        "is_fragpipe_qc": True,
    }

    # 1. Process Protein-Level Data
    prot_path = os.path.join(job_dir, "combined_protein.tsv")
    if os.path.exists(prot_path):
        prots = pd.read_csv(prot_path, sep="\t")
        metrics["nProts"] = len(prots)
        metrics["nContaProts"] = prots["Protein ID"].str.contains("cont_").sum()

    # 2. Process Peptide-Level Data
    pep_path = os.path.join(job_dir, "combined_peptide.tsv")
    if os.path.exists(pep_path):
        peps = pd.read_csv(pep_path, sep="\t")
        metrics["nPeps"] = len(peps)

        if pd.notna(tic_area) and tic_area > 0:
            # Handle possible missing intensity columns safely
            intensity_col = (
                "s_1 Intensity" if "s_1 Intensity" in peps.columns else peps.columns[0]
            )
            spec_count_col = (
                "s_1 Spectral Count"
                if "s_1 Spectral Count" in peps.columns
                else peps.columns[0]
            )

            tot_pep_i = peps[intensity_col].sum(skipna=True)
            metrics["totPepI"] = tot_pep_i
            metrics["explIons"] = round(100 * tot_pep_i / tic_area, 2)

            is_contam = peps["Protein ID"].str.contains("cont_", na=False)
            metrics["contaPepI"] = peps.loc[is_contam, intensity_col].sum(
                skipna=True
            )
            metrics["contaPepSpecCt"] = peps.loc[is_contam, spec_count_col].sum(
                skipna=True
            )
            metrics["nContaPeps"] = is_contam.sum()

    # 3. Process PSM-Level Data
    psm_path = os.path.join(job_dir, "s_1", "psm.tsv")
    if not os.path.exists(psm_path):
        # Fallback to search dynamically for any psm.tsv in subdirectories
        for root, dirs, files in os.walk(job_dir):
            if "psm.tsv" in files:
                psm_path = os.path.join(root, "psm.tsv")
                break

    if os.path.exists(psm_path):
        psms = pd.read_csv(psm_path, sep="\t")

        metrics["nPsms"] = len(psms)
        metrics["idRate"] = (
            round(100 * len(psms) / n_ms2, 2) if n_ms2 else np.nan
        )
        metrics["missCl"] = round(
            100 * (psms["Number of Missed Cleavages"] > 0).mean(), 2
        )

        # Modification metrics (string matching shifts)
        mods = psms["Assigned Modifications"].fillna("")
        metrics["oxiPsms"] = mods.str.contains(r"\(15\.9949\)", regex=True).sum()
        metrics["alkyCandPsms"] = (
            psms["Peptide"].fillna("").str.contains("C", fixed=True).sum()
        )

        carbamido_psms = mods.str.contains(r"\(57\.0215\)", regex=True)
        methth_psms = mods.str.contains(r"\(45\.9877\)", regex=True)

        metrics["carbamidoPsms"] = carbamido_psms.sum()
        metrics["meththPsms"] = methth_psms.sum()
        metrics["CarbamMeththPsms"] = (carbamido_psms & methth_psms).sum()

        metrics["carbamylKpsms"] = mods.str.contains(
            r"K\(43\.0058\)", regex=True
        ).sum()
        metrics["carbamylNtPsms"] = mods.str.contains(
            r"N-term\(43\.0058\)", regex=True
        ).sum()
        metrics["contaPsms"] = psms["Protein ID"].str.contains("cont_").sum()

        # Mass deviation calculations
        within_0c5 = psms["Delta Mass"].abs() < 0.5
        if within_0c5.any():
            dms = (
                1e6
                * psms.loc[within_0c5, "Delta Mass"]
                / psms.loc[within_0c5, "Calculated Peptide Mass"]
            )

            metrics["medMassDev"] = round(dms.median(), 2)
            metrics["q25MassDev"] = round(dms.quantile(0.25), 2)
            metrics["q75MassDev"] = round(dms.quantile(0.75), 2)
            metrics["q05MassDev"] = round(dms.quantile(0.05), 2)
            metrics["q10MassDev"] = round(dms.quantile(0.10), 2)
            metrics["q90MassDev"] = round(dms.quantile(0.90), 2)
            metrics["q95MassDev"] = round(dms.quantile(0.95), 2)

    return metrics

def run_proteomics_metrics(df: pd.DataFrame) -> dict[str, Any]:
    """
    Compute comprehensive proteomics QC metrics.

    Parameters
    ----------
    df : pd.DataFrame
        Abundance matrix where index=Proteins/Peptides, columns=Samples.

    Returns
    -------
    dict
        Proteomics QC results.
    """
    logger.info("Computing proteomics QC for %d proteins across %d samples", *df.shape)

    # 1. Active identifications per sample (non-zero and non-NaN values)
    id_counts = {}
    missing_rates = {}
    intensity_stats = {}
    
    # We log-transform (log2) quantitative values for robust QC stats
    log_df = np.log2(df.replace(0, np.nan))

    for col in df.columns:
        s = df[col]
        non_zero = s.dropna()
        non_zero = non_zero[non_zero > 0]
        id_counts[col] = int(len(non_zero))
        missing_rates[col] = float(s.isna().sum() / len(df))

        # Stats on log2 values
        log_s = log_df[col].dropna()
        if not log_s.empty:
            intensity_stats[col] = {
                "mean": float(log_s.mean()),
                "std": float(log_s.std()),
                "min": float(log_s.min()),
                "q25": float(log_s.quantile(0.25)),
                "median": float(log_s.median()),
                "q75": float(log_s.quantile(0.75)),
                "max": float(log_s.max()),
            }
        else:
            intensity_stats[col] = {
                "mean": 0.0, "std": 0.0, "min": 0.0, "q25": 0.0,
                "median": 0.0, "q75": 0.0, "max": 0.0
            }

    # 2. Pearson Correlation between samples
    # Make sure we drop all-NaN rows before correlation to avoid NaN outputs
    corr_df = df.dropna(how="all").corr(method="pearson")
    correlation_matrix = {
        "samples": list(df.columns),
        "grid": corr_df.values.tolist()
    }

    # 3. PCA Coordinates using SVD (Self-contained, no Scikit-learn required)
    pca_coords = _compute_pca(log_df)

    # 4. Outlier Flagging Rules for Proteomics
    # Thresholds:
    # - Low ID count: < 60% of the median sample's ID count
    # - High Missingness: > 45% missing values
    # - Median Intensity: log2 median is an outlier (Z-score > 2.5)
    id_vals = list(id_counts.values())
    med_id_count = np.median(id_vals) if id_vals else 1
    
    med_intensities = [stats["median"] for stats in intensity_stats.values()]
    med_intensities_mean = np.mean(med_intensities) if med_intensities else 0
    med_intensities_std = np.std(med_intensities) if med_intensities else 1
    if med_intensities_std == 0:
        med_intensities_std = 1

    per_sample_flags = []
    n_fail = 0
    for col in df.columns:
        flags = []
        if id_counts[col] < 0.6 * med_id_count:
            flags.append("LOW_PROTEIN_IDS")
        if missing_rates[col] > 0.45:
            flags.append("HIGH_MISSINGNESS")
            
        z_int = (intensity_stats[col]["median"] - med_intensities_mean) / med_intensities_std
        if abs(z_int) > 2.5:
            flags.append("ABUNDANCE_OUTLIER")

        status = "FAIL" if flags else "PASS"
        if status == "FAIL":
            n_fail += 1

        per_sample_flags.append({
            "sample_id": col,
            "status": status,
            "flags": flags,
            "metrics": {
                "proteins_identified": id_counts[col],
                "missing_rate": missing_rates[col],
                "median_log2_abundance": intensity_stats[col]["median"]
            }
        })

    return {
        "is_proteomics": True,
        "n_samples": len(df.columns),
        "n_proteins": len(df),
        "identifications": id_counts,
        "missing_rates": missing_rates,
        "intensity_stats": intensity_stats,
        "correlation": correlation_matrix,
        "pca": pca_coords,
        "per_sample_flags": per_sample_flags,
        "summary": {
            "n_pass": len(df.columns) - n_fail,
            "n_fail": n_fail,
            "pass_rate": round((len(df.columns) - n_fail) / max(len(df.columns), 1), 4),
        }
    }


def _compute_pca(log_df: pd.DataFrame) -> dict[str, list[float]]:
    """Compute first two Principal Components using raw Singular Value Decomposition."""
    # Filter proteins/rows with too many NaNs to prevent SVD breakdown
    # Keep proteins that are identified in at least 50% of the runs
    min_count = max(1, len(log_df.columns) // 2)
    filtered = log_df.dropna(thresh=min_count).copy()

    if filtered.empty or len(filtered.columns) < 2:
        logger.warning("Too few columns or samples to compute PCA.")
        return {"PC1": [0.0] * len(log_df.columns), "PC2": [0.0] * len(log_df.columns)}

    # Standard median imputation for remaining NaNs
    imputed = filtered.apply(lambda row: row.fillna(row.median()), axis=1)
    
    # Fallback for rows where median is still NaN
    imputed = imputed.fillna(imputed.stack().median())

    # SVD: Center matrix (Subtract protein row means)
    centered = imputed.sub(imputed.mean(axis=1), axis=0)

    # Center samples (columns) too
    centered = centered.sub(centered.mean(axis=0), axis=1)

    # Transform: samples as rows, proteins as columns
    X = centered.values.T  # Shape: Samples x Proteins

    try:
        # Singular Value Decomposition
        U, S, Vt = np.linalg.svd(X, full_matrices=False)
        # Compute Coordinates: U * S
        coords = U * S
        pc1 = coords[:, 0].tolist()
        pc2 = coords[:, 1].tolist()
    except Exception as e:
        logger.error("PCA SVD computation failed: %s. Using dummy coordinates.", e)
        pc1 = [0.0] * len(log_df.columns)
        pc2 = [0.0] * len(log_df.columns)

    return {"PC1": pc1, "PC2": pc2}


# ---------------------------------------------------------------------------
# Genomics Basic statistics
# ---------------------------------------------------------------------------

def compute_basic_stats(df: pd.DataFrame) -> dict[str, dict[str, float]]:
    """
    Compute per-column descriptive statistics for all numeric columns.
    """
    numeric = df.select_dtypes(include="number")
    result: dict[str, Any] = {}
    for col in numeric.columns:
        s = numeric[col].dropna()
        result[col] = {
            "count": int(s.count()),
            "mean": float(s.mean()),
            "std": float(s.std()),
            "min": float(s.min()),
            "q25": float(s.quantile(0.25)),
            "median": float(s.median()),
            "q75": float(s.quantile(0.75)),
            "max": float(s.max()),
        }
    return result


def compute_missing_rate(df: pd.DataFrame) -> dict[str, Any]:
    """
    Compute missing-value rates per column and per sample (Genomics).
    """
    n_rows, n_cols = df.shape
    per_col = {col: float(df[col].isna().sum() / n_rows) for col in df.columns}
    id_col = _get_id_column(df)
    if id_col:
        per_sample = {
            str(row[id_col]): float(row.drop(labels=[id_col]).isna().sum() / (n_cols - 1))
            for _, row in df.iterrows()
        }
    else:
        per_sample = {str(i): float(row.isna().sum() / n_cols) for i, row in df.iterrows()}

    total = float(df.isna().sum().sum() / (n_rows * n_cols))
    return {
        "per_column": per_col,
        "per_sample": per_sample,
        "total_missing_fraction": total,
    }


def compute_outliers(
    df: pd.DataFrame,
    method: str = "both",
    zscore_threshold: float = OUTLIER_ZSCORE,
    iqr_factor: float = OUTLIER_IQR_FACTOR,
) -> dict[str, Any]:
    """
    Flag outlier samples using Z-score and/or IQR methods (Genomics).
    """
    numeric = df.select_dtypes(include="number")
    id_col = _get_id_column(df)
    result: dict[str, Any] = {}

    for col in numeric.columns:
        s = numeric[col].dropna()
        col_result: dict[str, list] = {}

        if method in ("zscore", "both"):
            z = np.abs(stats.zscore(s))
            outlier_idx = s.index[z > zscore_threshold].tolist()
            col_result["zscore_outliers"] = _idx_to_ids(df, id_col, outlier_idx)

        if method in ("iqr", "both"):
            q1, q3 = s.quantile(0.25), s.quantile(0.75)
            iqr = q3 - iqr_factor * iqr
            lower, upper = q1 - iqr_factor * iqr, q3 + iqr_factor * iqr
            outlier_idx = s.index[(s < lower) | (s > upper)].tolist()
            col_result["iqr_outliers"] = _idx_to_ids(df, id_col, outlier_idx)

        result[col] = col_result
    return result


def compute_duplication_rate(df: pd.DataFrame) -> dict[str, Any]:
    """
    Detect duplicate and near-duplicate samples (Genomics).
    """
    numeric = df.select_dtypes(include="number")
    id_col = _get_id_column(df)
    groups: list[list[str]] = []
    dup_mask = numeric.duplicated(keep=False)
    dup_df = df[dup_mask]

    if not dup_df.empty:
        key_cols = list(numeric.columns)
        grouped = dup_df.groupby(key_cols, dropna=False)
        for _, grp in grouped:
            if len(grp) > 1:
                ids = grp[id_col].tolist() if id_col else grp.index.tolist()
                groups.append([str(x) for x in ids])
    return {"exact_duplicates": groups, "n_exact_duplicates": sum(len(g) for g in groups)}


def compute_gc_stats(
    df: pd.DataFrame,
    gc_col: str = "gc_content",
    low: float = GC_LOW,
    high: float = GC_HIGH,
) -> dict[str, Any]:
    """
    Summarise GC-content distribution (Genomics).
    """
    if gc_col not in df.columns:
        return {"available": False}
    s = df[gc_col].dropna()
    if s.max() > 1.0:
        s = s / 100.0
    id_col = _get_id_column(df)
    flagged_idx = s.index[(s < low) | (s > high)].tolist()
    return {
        "available": True,
        "mean": float(s.mean()),
        "std": float(s.std()),
        "min": float(s.min()),
        "max": float(s.max()),
        "low_threshold": low,
        "high_threshold": high,
        "flagged_samples": _idx_to_ids(df, id_col, flagged_idx),
        "n_flagged": len(flagged_idx),
    }


def compute_coverage_stats(
    df: pd.DataFrame,
    cov_col: str = "coverage",
    min_cov: float = MIN_COVERAGE,
) -> dict[str, Any]:
    """
    Summarise sequencing coverage depth (Genomics).
    """
    if cov_col not in df.columns:
        return {"available": False}
    s = df[cov_col].dropna()
    id_col = _get_id_column(df)
    low_cov_idx = s.index[s < min_cov].tolist()
    return {
        "available": True,
        "mean": float(s.mean()),
        "std": float(s.std()),
        "min": float(s.min()),
        "max": float(s.max()),
        "min_threshold": min_cov,
        "low_coverage_samples": _idx_to_ids(df, id_col, low_cov_idx),
        "n_low_coverage": len(low_cov_idx),
    }


def _build_per_sample_flags(
    df: pd.DataFrame,
    gc_stats: dict,
    coverage_stats: dict,
    outliers: dict,
) -> list[dict]:
    """Build per-sample PASS/FAIL flag table (Genomics)."""
    id_col = _get_id_column(df)
    records = []
    gc_flagged = set(gc_stats.get("flagged_samples", []))
    cov_flagged = set(coverage_stats.get("low_coverage_samples", []))
    outlier_flagged: set[str] = set()
    for col_data in outliers.values():
        outlier_flagged.update(col_data.get("zscore_outliers", []))
        outlier_flagged.update(col_data.get("iqr_outliers", []))

    ids = df[id_col].tolist() if id_col else [str(i) for i in df.index]
    for sid in ids:
        sid_str = str(sid)
        flags = []
        if sid_str in gc_flagged:
            flags.append("GC_OUT_OF_RANGE")
        if sid_str in cov_flagged:
            flags.append("LOW_COVERAGE")
        if sid_str in outlier_flagged:
            flags.append("OUTLIER")
        records.append({
            "sample_id": sid_str,
            "status": "FAIL" if flags else "PASS",
            "flags": flags,
        })
    return records


# ---------------------------------------------------------------------------
# Master runner
# ---------------------------------------------------------------------------

def run_all_metrics(df: pd.DataFrame) -> dict[str, Any]:
    """
    Run the full QC pipeline and return a consolidated results dictionary.
    Differentiates between Genomics and Proteomics structures automatically.
    """
    # If the dataframe does not have genomics total_reads, treat as Proteomics matrix
    if "total_reads" not in df.columns:
        return run_proteomics_metrics(df)

    logger.info("Running Genomics QC on %d samples, %d columns", *df.shape)
    basic = compute_basic_stats(df)
    missing = compute_missing_rate(df)
    outliers = compute_outliers(df)
    dups = compute_duplication_rate(df)
    gc = compute_gc_stats(df)
    cov = compute_coverage_stats(df)
    per_sample = _build_per_sample_flags(df, gc, cov, outliers)
    n_fail = sum(1 for s in per_sample if s["status"] == "FAIL")

    return {
        "is_proteomics": False,
        "n_samples": len(df),
        "n_columns": len(df.columns),
        "basic_stats": basic,
        "missing": missing,
        "outliers": outliers,
        "duplicates": dups,
        "gc_stats": gc,
        "coverage_stats": cov,
        "per_sample_flags": per_sample,
        "summary": {
            "n_pass": len(per_sample) - n_fail,
            "n_fail": n_fail,
            "pass_rate": round((len(per_sample) - n_fail) / max(len(per_sample), 1), 4),
        },
    }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _get_id_column(df: pd.DataFrame) -> str | None:
    """Return the sample-ID column name if present."""
    for candidate in ("sample_id", "SampleID", "sample", "id", "ID"):
        if candidate in df.columns:
            return candidate
    return None


def _idx_to_ids(df: pd.DataFrame, id_col: str | None, idx: list) -> list[str]:
    """Convert DataFrame row indices to sample-ID strings."""
    if id_col and id_col in df.columns:
        return [str(df.loc[i, id_col]) for i in idx if i in df.index]
    return [str(i) for i in idx]

