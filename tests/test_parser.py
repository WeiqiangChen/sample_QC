from pathlib import Path

from sample_qc.parser import load_proteomics_data


def test_load_proteomics_data_fragpipe():
    filepath = Path(__file__).parent.parent / "data" / "mock_fragpipe.tsv"
    df = load_proteomics_data(filepath, format_type="fragpipe")

    assert len(df) == 10
    assert len(df.columns) == 4
    assert df.index.name == "Protein ID"
    assert "Sample1" in df.columns
    assert df.loc["P00761", "Sample1"] == 1250000.0


def test_load_proteomics_data_spectronaut():
    filepath = Path(__file__).parent.parent / "data" / "mock_spectronaut.tsv"
    df = load_proteomics_data(filepath, format_type="spectronaut")

    assert len(df) == 5
    assert len(df.columns) == 4
    assert df.index.name == "PG.ProteinGroups"
    assert "Run_A" in df.columns
    assert df.loc["P00761", "Run_A"] == 1250000.0


def test_load_proteomics_data_matrix(tmp_path):
    content = "Protein ID,Sample_A,Sample_B\nP1,100,200\nP2,0,300\n"
    filepath = tmp_path / "matrix.csv"
    filepath.write_text(content)

    df = load_proteomics_data(filepath)
    assert len(df) == 2
    assert list(df.columns) == ["Sample_A", "Sample_B"]
    assert df.loc["P1", "Sample_A"] == 100
