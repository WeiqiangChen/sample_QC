from pathlib import Path

import pandas as pd
import pytest

from sample_qc.parser import parse_sample_sheet, validate_schema


def test_validate_schema_valid():
    df = pd.DataFrame(
        {
            "sample_id": ["S1"],
            "total_reads": [100],
            "mapped_reads": [90],
            "gc_content": [0.4],
            "coverage": [30.0],
            "dup_rate": [0.05],
        }
    )
    # Should run fine without errors
    validate_schema(df)


def test_validate_schema_invalid():
    df = pd.DataFrame({"sample_id": ["S1"]})
    with pytest.raises(ValueError, match="missing required column"):
        validate_schema(df)


def test_parse_sample_sheet_tsv(tmp_path):
    tsv_content = (
        "sample_id\ttotal_reads\tmapped_reads\tgc_content\tcoverage\tdup_rate\n"
        "S1\t100\t90\t0.4\t30.0\t0.05\n"
    )
    file_path = tmp_path / "test_samples.tsv"
    file_path.write_text(tsv_content)

    df = parse_sample_sheet(file_path)
    assert len(df) == 1
    assert df.loc[0, "sample_id"] == "S1"
    assert df.loc[0, "total_reads"] == 100


def test_load_proteomics_data_fragpipe():
    # Test loading actual FragPipe mock data
    filepath = Path(__file__).parent.parent / "data" / "mock_fragpipe.tsv"
    df = parse_sample_sheet # dummy import check, let's use the real load function
    from sample_qc.parser import load_proteomics_data
    
    df = load_proteomics_data(filepath, format_type="fragpipe")
    assert len(df) == 10  # 10 protein rows
    assert len(df.columns) == 4  # 4 samples (Sample1, Sample2, Sample3, Sample4)
    assert df.index.name == "Protein ID"
    assert "Sample1" in df.columns
    assert df.loc["P00761", "Sample1"] == 1250000.0


def test_load_proteomics_data_spectronaut():
    # Test loading actual Spectronaut mock data
    filepath = Path(__file__).parent.parent / "data" / "mock_spectronaut.tsv"
    from sample_qc.parser import load_proteomics_data
    
    df = load_proteomics_data(filepath, format_type="spectronaut")
    assert len(df) == 5  # 5 protein rows
    assert len(df.columns) == 4  # 4 samples (Run_A, Run_B, Run_C, Run_D)
    assert df.index.name == "PG.ProteinGroups"
    assert "Run_A" in df.columns
    assert df.loc["P00761", "Run_A"] == 1250000.0


def test_load_proteomics_data_auto():
    # Test auto-detection functionality
    proj_root = Path(__file__).parent.parent
    fragpipe_path = proj_root / "data" / "mock_fragpipe.tsv"
    spectronaut_path = proj_root / "data" / "mock_spectronaut.tsv"
    genomics_path = proj_root / "data" / "demo_sample.tsv"
    
    from sample_qc.parser import load_proteomics_data
    
    df_fp = load_proteomics_data(fragpipe_path)
    assert len(df_fp) == 10
    
    df_spec = load_proteomics_data(spectronaut_path)
    assert len(df_spec) == 5
    
    df_gen = load_proteomics_data(genomics_path)
    assert len(df_gen) == 15
    assert "sample_id" in df_gen.columns
