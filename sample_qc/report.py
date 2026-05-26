"""
report.py — Report generation for proteomics QC results.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

PROTEOMICS_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Proteomics QC Dashboard</title>
    <style>
        body { font-family: Arial, sans-serif; background: #0f172a; color: #e2e8f0; margin: 0; padding: 2rem; }
        .container { max-width: 1100px; margin: 0 auto; }
        h1, h2 { color: #f8fafc; }
        pre { background: #111827; padding: 1rem; border-radius: 12px; overflow-x: auto; }
        .card { background: #111827; border: 1px solid #334155; border-radius: 14px; padding: 1.2rem; margin-bottom: 1rem; }
        .grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 1rem; }
        .grid-item { background: #0f172a; border: 1px solid #334155; border-radius: 12px; padding: 1rem; }
        table { width: 100%; border-collapse: collapse; margin-top: 1rem; }
        th, td { padding: 0.75rem; border-bottom: 1px solid #334155; text-align: left; }
        th { background: #1e293b; }
    </style>
</head>
<body>
    <div class="container">
        <h1>Proteomics QC Dashboard</h1>
        <div class="grid">
            <div class="grid-item"><strong>Total Runs:</strong> {{TOTAL_RUNS}}</div>
            <div class="grid-item"><strong>Total Proteins:</strong> {{TOTAL_PROTEINS}}</div>
            <div class="grid-item"><strong>Passed QC:</strong> {{PASSED}}</div>
            <div class="grid-item"><strong>Failed QC:</strong> {{FAILED}}</div>
        </div>

        <div class="card">
            <h2>Run-level Quality Flags</h2>
            <table>
                <thead>
                    <tr><th>Run</th><th>Status</th><th>Proteins</th><th>Missing %</th><th>Flags</th></tr>
                </thead>
                <tbody>
                    {{TABLE_ROWS}}
                </tbody>
            </table>
        </div>

        <div class="card">
            <h2>Raw JSON Results</h2>
            <pre>{{DATA_PLACEHOLDER}}</pre>
        </div>
    </div>
</body>
</html>
"""

FRAGPIPE_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>FragPipe QC Dashboard</title>
    <style>
        body { font-family: Arial, sans-serif; background: #0f172a; color: #e2e8f0; margin: 0; padding: 2rem; }
        .container { max-width: 1100px; margin: 0 auto; }
        .card { background: #111827; border: 1px solid #334155; border-radius: 14px; padding: 1rem; margin-bottom: 1rem; }
        .grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 1rem; }
        .grid-item { background: #0f172a; border: 1px solid #334155; border-radius: 12px; padding: 1rem; }
        table { width: 100%; border-collapse: collapse; margin-top: 1rem; }
        th, td { padding: 0.75rem; border-bottom: 1px solid #334155; text-align: left; }
        th { background: #1e293b; }
    </style>
</head>
<body>
    <div class="container">
        <h1>FragPipe QC Dashboard</h1>
        <div class="grid">
            <div class="grid-item"><strong>Proteins:</strong> {{N_PROTS}}</div>
            <div class="grid-item"><strong>Peptides:</strong> {{N_PEPS}}</div>
            <div class="grid-item"><strong>PSMs:</strong> {{N_PSMS}}</div>
            <div class="grid-item"><strong>ID Rate:</strong> {{ID_RATE}}</div>
        </div>

        <div class="card">
            <h2>FragPipe QC Metrics</h2>
            <table>
                <thead>
                    <tr><th>Metric</th><th>Value</th></tr>
                </thead>
                <tbody>
                    {{METRIC_ROWS}}
                </tbody>
            </table>
        </div>

        <div class="card">
            <h2>Raw JSON Results</h2>
            <pre>{{DATA_PLACEHOLDER}}</pre>
        </div>
    </div>
</body>
</html>
"""


def _html_escape(text: str) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&#39;")
    )


def generate_json_report(results: dict[str, Any], output_path: Path) -> None:
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)


def generate_html_report(results: dict[str, Any], output_path: Path) -> None:
    data_json = json.dumps(results, indent=2, ensure_ascii=False)
    if results.get("is_fragpipe_qc", False):
        template = FRAGPIPE_HTML_TEMPLATE
        metric_rows = []
        for key in ["nProts", "nPeps", "nPsms", "idRate", "missCl", "totPepI", "explIons"]:
            if key in results:
                metric_rows.append(f"<tr><td>{_html_escape(key)}</td><td>{_html_escape(results[key])}</td></tr>")
        html = template.replace("{{N_PROTS}}", _html_escape(results.get("nProts", 0)))
        html = html.replace("{{N_PEPS}}", _html_escape(results.get("nPeps", 0)))
        html = html.replace("{{N_PSMS}}", _html_escape(results.get("nPsms", 0)))
        html = html.replace("{{ID_RATE}}", _html_escape(results.get("idRate", "N/A")))
        html = html.replace("{{METRIC_ROWS}}", "\n".join(metric_rows))
    else:
        template = PROTEOMICS_HTML_TEMPLATE
        table_rows = []
        for sample in results.get("per_sample_flags", []):
            flags = ", ".join(sample.get("flags", [])) or "-"
            metrics = sample.get("metrics", {})
            missing = f"{metrics.get('missing_rate', 0) * 100:.1f}%"
            table_rows.append(
                "<tr>"
                f"<td>{_html_escape(sample.get('sample_id', ''))}</td>"
                f"<td>{_html_escape(sample.get('status', ''))}</td>"
                f"<td>{_html_escape(metrics.get('proteins_identified', 0))}</td>"
                f"<td>{_html_escape(missing)}</td>"
                f"<td>{_html_escape(flags)}</td>"
                "</tr>"
            )
        html = template.replace("{{TOTAL_RUNS}}", _html_escape(results.get("n_samples", 0)))
        html = html.replace("{{TOTAL_PROTEINS}}", _html_escape(results.get("n_proteins", 0)))
        html = html.replace("{{PASSED}}", _html_escape(results.get("summary", {}).get("n_pass", 0)))
        html = html.replace("{{FAILED}}", _html_escape(results.get("summary", {}).get("n_fail", 0)))
        html = html.replace("{{TABLE_ROWS}}", "\n".join(table_rows))

    html = html.replace("{{DATA_PLACEHOLDER}}", _html_escape(data_json))
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)


def generate_text_report(results: dict[str, Any]) -> str:
    if results.get("is_fragpipe_qc", False):
        report = ["=" * 60, "FragPipe Proteomics QC Report", "=" * 60]
        for key in ["nProts", "nPeps", "nPsms", "idRate", "missCl", "totPepI", "explIons"]:
            if key in results:
                report.append(f"{key:<20}: {results[key]}")
        return "\n".join(report)

    if results.get("is_proteomics", False):
        report = ["=" * 60, "Proteomics QC Report", "=" * 60]
        report.append(f"Total Runs       : {results.get('n_samples', 0)}")
        report.append(f"Total Proteins   : {results.get('n_proteins', 0)}")
        report.append(f"Passed QC        : {results.get('summary', {}).get('n_pass', 0)}")
        report.append(f"Failed QC        : {results.get('summary', {}).get('n_fail', 0)}")
        report.append(f"Pass Rate        : {results.get('summary', {}).get('pass_rate', 0.0) * 100:.2f}%")
        report.append("\nRun-level details:")
        report.append(f"{'Run':<20} {'Status':<8} {'IDs':<8} {'Missing':<10} {'Flags'}")
        report.append("-" * 60)
        for sample in results.get("per_sample_flags", []):
            metrics = sample.get("metrics", {})
            missing = metrics.get("missing_rate", 0.0) * 100
            flags = ", ".join(sample.get("flags", [])) or "-"
            report.append(
                f"{sample.get('sample_id', ''):<20} {sample.get('status', ''):<8} "
                f"{metrics.get('proteins_identified', 0):<8} {missing:<9.1f}% {flags}"
            )
        return "\n".join(report)

    raise ValueError("Unsupported results type for text report.")


def generate_report(
    results: dict[str, Any],
    output_path: str | Path | None = None,
    fmt: str = "json",
) -> str | None:
    if fmt == "text" and output_path is None:
        return generate_text_report(results)
    if output_path is None:
        raise ValueError(f"An output_path must be supplied to save format '{fmt}'")

    out_p = Path(output_path)
    out_p.parent.mkdir(parents=True, exist_ok=True)

    if fmt == "json":
        generate_json_report(results, out_p)
    elif fmt == "html":
        generate_html_report(results, out_p)
    elif fmt == "text":
        text_content = generate_text_report(results)
        with open(out_p, "w", encoding="utf-8") as f:
            f.write(text_content)
    else:
        raise ValueError(f"Invalid report format '{fmt}'. Choose from: 'json', 'html', 'text'.")
    return None
