import pandas as pd
import pytest

from sample_qc.metrics import (
    compute_basic_stats,
    compute_coverage_stats,
    compute_duplication_rate,
    compute_gc_stats,
    compute_missing_rate,
    compute_outliers,
    run_all_metrics,
)


@pytest.fixture
def sample_df():
    data = {
        "sample_id": ["S1", "S2", "S3", "S4", "S5"],
        "total_reads": [1000, 1200, 900, 1100, 5000],  # S5 is outlier
        "mapped_reads": [900, 1100, 800, 1000, 4800],
        "gc_content": [0.45, 0.46, 0.44, 0.20, 0.45],  # S4 is out of range
        "coverage": [30.0, 32.0, 5.0, 31.0, 35.0],    # S3 is low coverage
        "dup_rate": [0.05, 0.06, 0.05, 0.05, 0.06],
    }
    return pd.DataFrame(data)


def test_compute_basic_stats(sample_df):
    stats = compute_basic_stats(sample_df)
    assert "total_reads" in stats
    assert stats["total_reads"]["count"] == 5
    assert stats["total_reads"]["min"] == 900
    assert stats["total_reads"]["max"] == 5000


def test_compute_missing_rate(sample_df):
    # Add a missing value to test missingness calculation
    df_missing = sample_df.copy()
    df_missing.loc[0, "coverage"] = None

    rates = compute_missing_rate(df_missing)
    assert rates["total_missing_fraction"] > 0
    assert rates["per_column"]["coverage"] == 0.2  # 1 out of 5 is missing


def test_compute_outliers(sample_df):
    # High Z-score threshold to keep test deterministic, S5 should be flagged
    outliers = compute_outliers(sample_df, method="zscore", zscore_threshold=1.5)
    assert "S5" in outliers["total_reads"]["zscore_outliers"]


def test_compute_duplication_rate(sample_df):
    # Add a near duplicate
    df_dup = sample_df.copy()
    df_dup.loc[1] = df_dup.loc[0]  # Make row index 1 identical to row index 0

    dups = compute_duplication_rate(df_dup)
    assert dups["n_exact_duplicates"] == 2
    assert "S1" in dups["exact_duplicates"][0]
    assert "S2" in dups["exact_duplicates"][0]


def test_compute_gc_stats(sample_df):
    gc = compute_gc_stats(sample_df, low=0.35, high=0.65)
    assert gc["available"] is True
    assert gc["n_flagged"] == 1
    assert "S4" in gc["flagged_samples"]


def test_compute_coverage_stats(sample_df):
    cov = compute_coverage_stats(sample_df, min_cov=10.0)
    assert cov["available"] is True
    assert cov["n_low_coverage"] == 1
    assert "S3" in cov["low_coverage_samples"]


def test_run_all_metrics(sample_df):
    res = run_all_metrics(sample_df)
    assert "summary" in res
    assert "basic_stats" in res
    assert "per_sample_flags" in res
    assert len(res["per_sample_flags"]) == 5


def test_run_proteomics_metrics():
    from pathlib import Path
    from sample_qc.parser import load_proteomics_data
    
    filepath = Path(__file__).parent.parent / "data" / "mock_fragpipe.tsv"
    df = load_proteomics_data(filepath, format_type="fragpipe")
    res = run_all_metrics(df)
    
    assert res["is_proteomics"] is True
    assert res["n_samples"] == 4
    assert res["n_proteins"] == 10
    assert "correlation" in res
    assert "pca" in res
    assert "per_sample_flags" in res
    assert len(res["per_sample_flags"]) == 4
    
    # Check SVD-based PCA output
    assert "PC1" in res["pca"]
    assert "PC2" in res["pca"]
    assert len(res["pca"]["PC1"]) == 4
    
    # Check correlation matrix is 4x4
    assert len(res["correlation"]["samples"]) == 4
    assert len(res["correlation"]["grid"]) == 4
    assert len(res["correlation"]["grid"][0]) == 4
