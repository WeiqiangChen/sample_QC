import json
from pathlib import Path

import pytest

from sample_qc.metrics import run_all_metrics
from sample_qc.parser import load_proteomics_data
from sample_qc.report import generate_report


@pytest.fixture
def qc_results():
    filepath = Path(__file__).parent.parent / "data" / "mock_fragpipe.tsv"
    df = load_proteomics_data(filepath, format_type="fragpipe")
    return run_all_metrics(df)


def test_generate_report_json(qc_results, tmp_path):
    out_path = tmp_path / "report.json"
    generate_report(qc_results, output_path=out_path, fmt="json")

    assert out_path.exists()
    with open(out_path, "r", encoding="utf-8") as f:
        loaded = json.load(f)
    assert loaded["n_samples"] == 4


def test_generate_report_html(qc_results, tmp_path):
    out_path = tmp_path / "report.html"
    generate_report(qc_results, output_path=out_path, fmt="html")

    assert out_path.exists()
    content = out_path.read_text(encoding="utf-8")
    assert "Proteomics QC Dashboard" in content
    assert "DATA_PLACEHOLDER" not in content


def test_generate_report_text(qc_results, tmp_path):
    out_path = tmp_path / "report.txt"
    generate_report(qc_results, output_path=out_path, fmt="text")
    assert out_path.exists()
    content = out_path.read_text(encoding="utf-8")
    assert "Proteomics QC Report" in content

    raw_text = generate_report(qc_results, fmt="text")
    assert "Passed QC" in raw_text
