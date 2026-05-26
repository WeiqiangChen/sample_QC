# sample_QC

A Python package for **proteomics quality control (QC)** of mass spectrometry data.
It provides fast data parsing, protein/peptide-level QC metrics, run-level flagging,
and HTML/JSON/text report generation for FragPipe, Spectronaut, and generic proteomics matrices.

---

## Features

- 🔬 **Proteomics matrix parsing** — FragPipe, Spectronaut, or generic wide abundance tables
- 📊 **Run-level QC metrics** — identified proteins, missingness, PCA, and sample correlations
- ⚠️ **Flagging of problematic runs** — low IDs, high missingness, and abundance outliers
- 🎯 **FragPipe advanced QC** — PSM, peptide, and contamination metrics for DDA searches
- 🧾 **Multi-format reports** — JSON, plain text, and HTML dashboards
- 🖥️ **CLI interface** — simple sample-qc command-line workflow

---

## Installation

`ash
pip install -e .
`

---

## Quick Start

`ash
sample-qc run data/mock_fragpipe.tsv --format all
sample-qc run data/mock_spectronaut.tsv --format html
sample-qc run data/FragPipe_result1 --format json
sample-qc version
`

---

## Input Formats

Supported input types:

- **FragPipe** combined_protein.tsv file or project directory
- **Spectronaut** exported report file (wide or long format)
- **Generic proteomics matrix** with protein/peptide identifier as the first column

Example files are available under data/.

---

## Output

- json — raw QC metrics for downstream processing
- html — interactive dashboard report
- 	ext — compact console-friendly summary

---

## Project Structure

`
sample_QC/
├── sample_qc/          # Core Python package
│   ├── __init__.py
│   ├── cli.py          # Click CLI
│   ├── metrics.py      # Proteomics QC metric computations
│   ├── parser.py       # Proteomics input parsers
│   └── report.py       # Report generation
├── tests/              # Unit tests
├── data/               # Proteomics example data
├── pyproject.toml      # Package metadata
└── README.md
`

---

## Development

`ash
pip install -e .
pytest tests/ -v
`

---

## License

MIT © 2026 WeiQiang
