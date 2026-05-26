from pathlib import Path

from sample_qc.metrics import generate_qc_metrics, run_all_metrics
from sample_qc.parser import load_proteomics_data


def test_run_all_metrics_proteomics_fragpipe():
    filepath = Path(__file__).parent.parent / "data" / "mock_fragpipe.tsv"
    df = load_proteomics_data(filepath, format_type="fragpipe")
    res = run_all_metrics(df)

    assert res["is_proteomics"] is True
    assert res["n_samples"] == 4
    assert res["n_proteins"] == 10
    assert "correlation" in res
    assert "pca" in res
    assert len(res["per_sample_flags"]) == 4


def test_generate_qc_metrics():
    job_dir = Path(__file__).parent.parent / "data" / "FragPipe_result1"
    res = generate_qc_metrics(job_dir, n_ms2=250000, tic_area=1e9)

    assert res["is_proteomics"] is True
    assert res["is_fragpipe_qc"] is True
    assert res["nProts"] > 0
    assert res["nPeps"] > 0
    assert res["nPsms"] > 0
    assert "medMassDev" in res
    assert "oxiPsms" in res
