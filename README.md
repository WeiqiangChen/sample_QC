# 🧬 sample_QC

A Python package for **sample quality control (QC)** of genomic / sequencing data.
It computes common QC metrics, detects outliers, and generates summary reports in
JSON, plain-text, or interactive HTML formats.

---

## Features

- 📊 **Basic statistics** — mean, std, min/max per metric column
- 🔍 **Outlier detection** — Z-score and IQR-based flagging
- 🔁 **Duplicate detection** — exact and near-duplicate samples
- ❓ **Missing-data audit** — per-sample and per-column missing rates
- 🧬 **GC-content check** — from sequence columns or pre-computed values
- 📈 **Coverage summary** — depth/coverage distribution stats
- 📝 **Multi-format reports** — JSON, plain text, interactive HTML (Plotly)
- 🖥️ **CLI interface** — easy command-line usage via `sample-qc`

---

## Installation

```bash
# Clone the repo
git clone https://github.com/WeiQiang/sample_QC.git
cd sample_QC

# Install in editable mode (recommended for development)
pip install -e ".[dev]"
```

---

## Quick Start

### CLI

```bash
# Run QC on a TSV/CSV sample sheet
sample-qc run data/demo_sample.tsv

# Generate an HTML report from saved results
sample-qc report results/qc_results.json --output report.html

# Show version
sample-qc version
```

### Python API

```python
from sample_qc import run_qc, generate_report

# Run QC
results = run_qc("data/demo_sample.tsv")

# Generate HTML report
generate_report(results, output_path="report.html", fmt="html")
```

---

## Input Format

The tool accepts **TSV or CSV** sample sheets. Required columns:

| Column         | Description                          |
|----------------|--------------------------------------|
| `sample_id`    | Unique sample identifier             |
| `total_reads`  | Total number of reads                |
| `mapped_reads` | Number of mapped reads               |
| `gc_content`   | GC content (0–100 or 0.0–1.0)        |
| `coverage`     | Mean sequencing coverage depth       |
| `dup_rate`     | Duplication rate (0.0–1.0)           |

Additional numeric columns are automatically included in statistical summaries.

See [`data/demo_sample.tsv`](data/demo_sample.tsv) for an example.

---

## Output

### JSON (`--format json`)
```json
{
  "n_samples": 20,
  "metrics": { ... },
  "outliers": { ... },
  "missing": { ... },
  "flagged_samples": [ ... ]
}
```

### HTML Report
Interactive dashboard with:
- Metric distribution histograms
- Outlier scatter plots
- Per-sample QC pass/fail table

---

## Development

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=sample_qc --cov-report=term-missing
```

---

## Project Structure

```
sample_QC/
├── sample_qc/          # Core Python package
│   ├── __init__.py
│   ├── cli.py          # Click CLI
│   ├── metrics.py      # QC metric computations
│   ├── parser.py       # Input file parsers
│   └── report.py       # Report generators
├── tests/              # Unit tests (pytest)
├── data/               # Demo data
├── docs/               # Documentation
├── pyproject.toml      # Package metadata & build config
└── README.md
```

---

## License

MIT © 2026 WeiQiang
