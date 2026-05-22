import json
from pathlib import Path

import pytest

from sample_qc.metrics import run_all_metrics
from sample_qc.report import generate_report


@pytest.fixture
def qc_results():
    import pandas as pd
    data = {
        "sample_id": ["S1", "S2"],
        "total_reads": [1000, 1200],
        "mapped_reads": [900, 1100],
        "gc_content": [0.45, 0.46],
        "coverage": [30.0, 32.0],
        "dup_rate": [0.05, 0.06],
    }
    df = pd.DataFrame(data)
    return run_all_metrics(df)


def test_generate_report_json(qc_results, tmp_path):
    out_path = tmp_path / "report.json"
    generate_report(qc_results, output_path=out_path, fmt="json")
    
    assert out_path.exists()
    with open(out_path, "r") as f:
        loaded = json.load(f)
    assert loaded["n_samples"] == 2


def test_generate_report_html(qc_results, tmp_path):
    out_path = tmp_path / "report.html"
    generate_report(qc_results, output_path=out_path, fmt="html")
    
    assert out_path.exists()
    content = out_path.read_text()
    assert "Sample Quality Control Report" in content
    assert "DATA_PLACEHOLDER" not in content  # placeholder should be replaced


def test_generate_report_text(qc_results, tmp_path):
    # Test file-writing mode
    out_path = tmp_path / "report.txt"
    generate_report(qc_results, output_path=out_path, fmt="text")
    assert out_path.exists()
    content = out_path.read_text()
    assert "SAMPLE QUALITY CONTROL REPORT SUMMARY" in content

    # Test raw stdout mode
    raw_text = generate_report(qc_results, fmt="text")
    assert "Passed QC Filter" in raw_text


def test_generate_report_proteomics(tmp_path):
    from pathlib import Path
    from sample_qc.parser import load_proteomics_data
    from sample_qc.metrics import run_all_metrics
    
    filepath = Path(__file__).parent.parent / "data" / "mock_fragpipe.tsv"
    df = load_proteomics_data(filepath, format_type="fragpipe")
    res = run_all_metrics(df)
    
    # 1. HTML
    html_path = tmp_path / "report_prot.html"
    generate_report(res, output_path=html_path, fmt="html")
    assert html_path.exists()
    content = html_path.read_text()
    assert "Proteomics QC Dashboard" in content
    
    # 2. Text
    text_path = tmp_path / "report_prot.txt"
    generate_report(res, output_path=text_path, fmt="text")
    assert text_path.exists()
    txt_content = text_path.read_text()
    assert "PROTEOMICS SAMPLE QUALITY CONTROL REPORT SUMMARY" in txt_content
    assert "Run/Sample ID" in txt_content
