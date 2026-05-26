"""
metrics.py — Core proteomics QC metric computations.
"""

from __future__ import annotations

import os
import logging
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


def generate_qc_metrics(job_dir: str | os.PathLike, n_ms2: int | None = None, tic_area: float | None = None) -> dict[str, Any]:
    metrics = {
        "is_proteomics": True,
        "is_fragpipe_qc": True,
    }

    prot_path = os.path.join(job_dir, "combined_protein.tsv")
    if os.path.exists(prot_path):
        prots = pd.read_csv(prot_path, sep="\t")
        metrics["nProts"] = len(prots)
        if "Protein ID" in prots.columns:
            metrics["nContaProts"] = prots["Protein ID"].str.contains("cont_").sum()
        else:
            metrics["nContaProts"] = 0

    pep_path = os.path.join(job_dir, "combined_peptide.tsv")
    if os.path.exists(pep_path):
        peps = pd.read_csv(pep_path, sep="\t")
        metrics["nPeps"] = len(peps)

        if pd.notna(tic_area) and tic_area > 0:
            # Prefer common intensity and spectral count column names, fall back to first numeric column
            intensity_col = None
            spec_count_col = None
            for c in peps.columns:
                if "Intensity" in c and intensity_col is None:
                    intensity_col = c
                if "Spectral Count" in c and spec_count_col is None:
                    spec_count_col = c

            if intensity_col is None:
                num_cols = peps.select_dtypes(include="number").columns
                intensity_col = num_cols[0] if len(num_cols) > 0 else peps.columns[0]
            if spec_count_col is None:
                num_cols = peps.select_dtypes(include="number").columns
                spec_count_col = num_cols[0] if len(num_cols) > 0 else peps.columns[0]

            peps[intensity_col] = pd.to_numeric(peps[intensity_col], errors="coerce")
            peps[spec_count_col] = pd.to_numeric(peps[spec_count_col], errors="coerce")

            tot_pep_i = peps[intensity_col].sum(skipna=True)
            metrics["totPepI"] = tot_pep_i
            metrics["explIons"] = round(100 * float(tot_pep_i) / float(tic_area), 2)

            is_contam = peps.get("Protein ID", pd.Series([], dtype=str)).fillna("").str.contains("cont_", na=False)
            metrics["contaPepI"] = peps.loc[is_contam, intensity_col].sum(skipna=True)
            metrics["contaPepSpecCt"] = peps.loc[is_contam, spec_count_col].sum(skipna=True)
            metrics["nContaPeps"] = int(is_contam.sum())

    psm_path = os.path.join(job_dir, "s_1", "psm.tsv")
    if not os.path.exists(psm_path):
        for root, _, files in os.walk(job_dir):
            if "psm.tsv" in files:
                psm_path = os.path.join(root, "psm.tsv")
                break

    if os.path.exists(psm_path):
        psms = pd.read_csv(psm_path, sep="\t")

        metrics["nPsms"] = len(psms)
        metrics["idRate"] = round(100 * len(psms) / n_ms2, 2) if n_ms2 else np.nan
        metrics["missCl"] = round(100 * (psms.get("Number of Missed Cleavages", 0) > 0).mean(), 2)

        mods = psms.get("Assigned Modifications", pd.Series([], dtype=str)).fillna("")
        metrics["oxiPsms"] = mods.str.contains(r"\(15\.9949\)", regex=True).sum()
        metrics["alkyCandPsms"] = psms.get("Peptide", pd.Series([], dtype=str)).fillna("").str.contains("C", regex=False).sum()

        carbamido_psms = mods.str.contains(r"\(57\.0215\)", regex=True)
        methth_psms = mods.str.contains(r"\(45\.9877\)", regex=True)

        metrics["carbamidoPsms"] = carbamido_psms.sum()
        metrics["meththPsms"] = methth_psms.sum()
        metrics["CarbamMeththPsms"] = (carbamido_psms & methth_psms).sum()

        metrics["carbamylKpsms"] = mods.str.contains(r"K\(43\.0058\)", regex=True).sum()
        metrics["carbamylNtPsms"] = mods.str.contains(r"N-term\(43\.0058\)", regex=True).sum()
        metrics["contaPsms"] = psms.get("Protein ID", pd.Series([], dtype=str)).fillna("").str.contains("cont_").sum()

        within_0c5 = psms.get("Delta Mass", pd.Series([], dtype=float)).abs() < 0.5
        if within_0c5.any():
            dms = 1e6 * psms.loc[within_0c5, "Delta Mass"] / psms.loc[within_0c5, "Calculated Peptide Mass"]
            metrics["medMassDev"] = round(dms.median(), 2)
            metrics["q25MassDev"] = round(dms.quantile(0.25), 2)
            metrics["q75MassDev"] = round(dms.quantile(0.75), 2)
            metrics["q05MassDev"] = round(dms.quantile(0.05), 2)
            metrics["q10MassDev"] = round(dms.quantile(0.10), 2)
            metrics["q90MassDev"] = round(dms.quantile(0.90), 2)
            metrics["q95MassDev"] = round(dms.quantile(0.95), 2)

    return metrics


def run_proteomics_metrics(df: pd.DataFrame) -> dict[str, Any]:
    logger.info("Computing proteomics QC for %d proteins across %d runs", *df.shape)

    id_counts = {}
    missing_rates = {}
    intensity_stats = {}
    log_df = np.log2(df.replace(0, np.nan))

    for col in df.columns:
        s = df[col]
        non_zero = s.dropna()
        non_zero = non_zero[non_zero > 0]
        id_counts[col] = int(len(non_zero))
        missing_rates[col] = float(s.isna().sum() / len(df))

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
                "mean": 0.0,
                "std": 0.0,
                "min": 0.0,
                "q25": 0.0,
                "median": 0.0,
                "q75": 0.0,
                "max": 0.0,
            }

    corr_df = df.dropna(how="all").corr(method="pearson")
    correlation_matrix = {"samples": list(df.columns), "grid": corr_df.values.tolist()}
    pca_coords = _compute_pca(log_df)

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
                "median_log2_abundance": intensity_stats[col]["median"],
            },
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
        },
    }


def _compute_pca(log_df: pd.DataFrame) -> dict[str, list[float]]:
    min_count = max(1, len(log_df.columns) // 2)
    filtered = log_df.dropna(thresh=min_count).copy()

    if filtered.empty or len(filtered.columns) < 2:
        logger.warning("Too few columns or samples to compute PCA.")
        return {"PC1": [0.0] * len(log_df.columns), "PC2": [0.0] * len(log_df.columns)}

    imputed = filtered.apply(lambda row: row.fillna(row.median()), axis=1)
    imputed = imputed.fillna(imputed.stack().median())
    centered = imputed.sub(imputed.mean(axis=1), axis=0)
    centered = centered.sub(centered.mean(axis=0), axis=1)
    X = centered.values.T

    try:
        U, S, Vt = np.linalg.svd(X, full_matrices=False)
        coords = U * S
        return {"PC1": coords[:, 0].tolist(), "PC2": coords[:, 1].tolist()}
    except Exception as e:
        logger.error("PCA SVD computation failed: %s. Using zero coordinates.", e)
        return {"PC1": [0.0] * len(log_df.columns), "PC2": [0.0] * len(log_df.columns)}


def run_all_metrics(df: pd.DataFrame) -> dict[str, Any]:
    return run_proteomics_metrics(df)
