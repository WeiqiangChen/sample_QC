"""
report.py — Report generators (HTML, JSON, Text) for Genomics & Proteomics sample_qc results.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

# ===========================================================================
# Genomics HTML Dashboard Template
# ===========================================================================
GENOMICS_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>🧬 Sample Quality Control Report</title>
    <!-- Google Fonts -->
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
    <!-- Plotly.js -->
    <script src="https://cdn.plot.ly/plotly-2.24.1.min.js"></script>
    <style>
        :root {
            --bg-color: #0b0f19;
            --card-bg: rgba(17, 24, 39, 0.7);
            --card-border: rgba(255, 255, 255, 0.08);
            --text-main: #f3f4f6;
            --text-muted: #9ca3af;
            --primary: #6366f1;
            --primary-glow: rgba(99, 102, 241, 0.15);
            --success: #10b981;
            --success-glow: rgba(16, 185, 129, 0.15);
            --warning: #f59e0b;
            --warning-glow: rgba(245, 158, 11, 0.15);
            --danger: #ef4444;
            --danger-glow: rgba(239, 68, 68, 0.15);
            --font-main: 'Inter', -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            --font-mono: 'JetBrains Mono', monospace;
        }

        * { box-sizing: border-box; margin: 0; padding: 0; }
        body {
            background-color: var(--bg-color);
            color: var(--text-main);
            font-family: var(--font-main);
            line-height: 1.5;
            padding: 2rem;
            min-height: 100vh;
            background-image: radial-gradient(circle at 10% 20%, rgba(99, 102, 241, 0.05) 0%, transparent 40%),
                              radial-gradient(circle at 90% 80%, rgba(16, 185, 129, 0.04) 0%, transparent 40%);
            background-attachment: fixed;
        }
        .container { max-width: 1400px; margin: 0 auto; }
        header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 2rem;
            padding-bottom: 1.5rem;
            border-bottom: 1px solid var(--card-border);
        }
        .brand h1 {
            font-size: 2.2rem;
            font-weight: 800;
            background: linear-gradient(135deg, #a5b4fc 0%, #6366f1 50%, #4338ca 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }
        .brand p { color: var(--text-muted); font-size: 0.95rem; margin-top: 0.25rem; }
        .timestamp {
            font-family: var(--font-mono);
            font-size: 0.85rem;
            color: var(--text-muted);
            background: var(--card-border);
            padding: 0.5rem 1rem;
            border-radius: 9999px;
            border: 1px solid var(--card-border);
        }
        .summary-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 1.5rem; margin-bottom: 2rem; }
        .stat-card {
            background: var(--card-bg);
            border: 1px solid var(--card-border);
            border-radius: 16px;
            padding: 1.5rem;
            backdrop-filter: blur(12px);
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            position: relative;
            overflow: hidden;
        }
        .stat-card:hover { transform: translateY(-4px); border-color: var(--primary); box-shadow: 0 10px 20px -10px var(--primary-glow); }
        .stat-card.pass:hover { border-color: var(--success); box-shadow: 0 10px 20px -10px var(--success-glow); }
        .stat-card.fail:hover { border-color: var(--danger); box-shadow: 0 10px 20px -10px var(--danger-glow); }
        .stat-card::before { content: ''; position: absolute; top: 0; left: 0; width: 100%; height: 4px; background: var(--primary); }
        .stat-card.pass::before { background: var(--success); }
        .stat-card.fail::before { background: var(--danger); }
        .stat-card.warning::before { background: var(--warning); }
        .stat-label { font-size: 0.85rem; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; color: var(--text-muted); margin-bottom: 0.5rem; }
        .stat-value { font-size: 2.2rem; font-weight: 800; color: var(--text-main); font-family: var(--font-mono); }
        .layout-main { display: grid; grid-template-columns: 3fr 2fr; gap: 1.5rem; margin-bottom: 2rem; }
        @media (max-width: 1024px) { .layout-main { grid-template-columns: 1fr; } }
        .card { background: var(--card-bg); border: 1px solid var(--card-border); border-radius: 16px; padding: 1.5rem; backdrop-filter: blur(12px); margin-bottom: 1.5rem; }
        .card-title { font-size: 1.25rem; font-weight: 700; margin-bottom: 1rem; display: flex; align-items: center; gap: 0.5rem; border-bottom: 1px solid var(--card-border); padding-bottom: 0.75rem; }
        .chart-container { width: 100%; height: 400px; border-radius: 8px; overflow: hidden; }
        .table-wrapper { max-height: 500px; overflow-y: auto; border-radius: 8px; border: 1px solid var(--card-border); }
        table { width: 100%; border-collapse: collapse; text-align: left; font-size: 0.9rem; }
        th { background: #111827; padding: 0.75rem 1rem; font-weight: 600; color: var(--text-main); position: sticky; top: 0; z-index: 10; border-bottom: 2px solid var(--card-border); }
        td { padding: 0.75rem 1rem; border-bottom: 1px solid var(--card-border); color: var(--text-muted); font-family: var(--font-mono); }
        tr:hover td { background: rgba(255, 255, 255, 0.02); color: var(--text-main); }
        .search-container { margin-bottom: 1rem; }
        .search-input { width: 100%; padding: 0.75rem 1rem; border-radius: 8px; border: 1px solid var(--card-border); background: rgba(255, 255, 255, 0.05); color: var(--text-main); font-size: 0.9rem; transition: all 0.3s; }
        .search-input:focus { outline: none; border-color: var(--primary); box-shadow: 0 0 0 3px var(--primary-glow); }
        .badge { display: inline-flex; align-items: center; padding: 0.25rem 0.6rem; border-radius: 9999px; font-size: 0.75rem; font-weight: 600; text-transform: uppercase; }
        .badge-pass { background: var(--success-glow); color: var(--success); border: 1px solid rgba(16, 185, 129, 0.3); }
        .badge-fail { background: var(--danger-glow); color: var(--danger); border: 1px solid rgba(239, 68, 68, 0.3); animation: pulse-fail 2s infinite; }
        @keyframes pulse-fail {
            0% { box-shadow: 0 0 0 0 rgba(239, 68, 68, 0.4); }
            70% { box-shadow: 0 0 0 6px rgba(239, 68, 68, 0); }
            100% { box-shadow: 0 0 0 0 rgba(239, 68, 68, 0); }
        }
        .flag-tag { background: rgba(245, 158, 11, 0.1); color: var(--warning); border: 1px solid rgba(245, 158, 11, 0.2); padding: 0.1rem 0.4rem; border-radius: 4px; font-size: 0.7rem; margin-right: 0.3rem; display: inline-block; }
    </style>
</head>
<body>
    <div class="container">
        <header>
            <div class="brand">
                <h1>🧬 sample_QC Report</h1>
                <p>Interactive Genomic Quality Control Analysis Summary</p>
            </div>
            <div class="timestamp" id="report-time">Generated: --</div>
        </header>

        <div class="summary-grid">
            <div class="stat-card">
                <div class="stat-label">Total Samples</div>
                <div class="stat-value" id="stat-total-samples">0</div>
            </div>
            <div class="stat-card pass">
                <div class="stat-label">Passed QC</div>
                <div class="stat-value" id="stat-passed-samples">0</div>
            </div>
            <div class="stat-card fail">
                <div class="stat-label">Failed QC</div>
                <div class="stat-value" id="stat-failed-samples">0</div>
            </div>
            <div class="stat-card warning">
                <div class="stat-label">QC Pass Rate</div>
                <div class="stat-value" id="stat-pass-rate">0.0%</div>
            </div>
        </div>

        <div class="layout-main">
            <div class="plots-column">
                <div class="card">
                    <div class="card-title">📈 Metric Distribution & Outlier Analysis</div>
                    <div class="chart-container" id="plotly-distribution"></div>
                </div>
                <div class="card">
                    <div class="card-title">🔬 Sequencing Depth vs GC Content</div>
                    <div class="chart-container" id="plotly-scatter"></div>
                </div>
            </div>

            <div class="table-column">
                <div class="card">
                    <div class="card-title">📋 Per-Sample QC Flags</div>
                    <div class="search-container">
                        <input type="text" id="sample-search" class="search-input" placeholder="🔍 Search samples by ID or status...">
                    </div>
                    <div class="table-wrapper">
                        <table id="qc-table">
                            <thead>
                                <tr>
                                    <th>Sample ID</th>
                                    <th>Status</th>
                                    <th>Flagged Violations</th>
                                </tr>
                            </thead>
                            <tbody id="qc-table-body"></tbody>
                        </table>
                    </div>
                </div>
                <div class="card">
                    <div class="card-title">⚙️ Metadata Check Summary</div>
                    <div id="metadata-details" style="font-size: 0.9rem;"></div>
                </div>
            </div>
        </div>
    </div>

    <script>
        const qcData = {{DATA_PLACEHOLDER}};
        document.getElementById('report-time').innerText = `Generated: ${new Date().toLocaleString()}`;
        document.getElementById('stat-total-samples').innerText = qcData.n_samples;
        document.getElementById('stat-passed-samples').innerText = qcData.summary.n_pass;
        document.getElementById('stat-failed-samples').innerText = qcData.summary.n_fail;
        document.getElementById('stat-pass-rate').innerText = `${(qcData.summary.pass_rate * 100).toFixed(1)}%`;

        const tbody = document.getElementById('qc-table-body');
        function renderTable(filterText = '') {
            tbody.innerHTML = '';
            const filtered = qcData.per_sample_flags.filter(s => {
                const searchLower = filterText.toLowerCase();
                return s.sample_id.toLowerCase().includes(searchLower) || 
                       s.status.toLowerCase().includes(searchLower) ||
                       s.flags.some(f => f.toLowerCase().includes(searchLower));
            });

            filtered.forEach(s => {
                const tr = document.createElement('tr');
                const tdId = document.createElement('td'); tdId.innerText = s.sample_id;
                const tdStatus = document.createElement('td');
                const badge = document.createElement('span');
                badge.className = `badge badge-${s.status.toLowerCase()}`;
                badge.innerText = s.status;
                tdStatus.appendChild(badge);

                const tdFlags = document.createElement('td');
                if (s.flags.length === 0) { tdFlags.innerText = '-'; }
                else {
                    s.flags.forEach(f => {
                        const tag = document.createElement('span'); tag.className = 'flag-tag'; tag.innerText = f;
                        tdFlags.appendChild(tag);
                    });
                }
                tr.appendChild(tdId); tr.appendChild(tdStatus); tr.appendChild(tdFlags);
                tbody.appendChild(tr);
            });
        }
        document.getElementById('sample-search').addEventListener('input', (e) => { renderTable(e.target.value); });
        renderTable();

        const metaDiv = document.getElementById('metadata-details');
        metaDiv.innerHTML = `
            <div style="margin-bottom: 0.75rem;">
                <strong>Missing Data Ratio:</strong> 
                <span style="font-family: var(--font-mono); color: ${qcData.missing.total_missing_fraction > 0.05 ? 'var(--warning)' : 'var(--success)'}">
                    ${(qcData.missing.total_missing_fraction * 100).toFixed(2)}%
                </span>
            </div>
            <div style="margin-bottom: 0.75rem;">
                <strong>Exact Duplicate Group Count:</strong> 
                <span style="font-family: var(--font-mono); color: ${qcData.duplicates.n_exact_duplicates > 0 ? 'var(--danger)' : 'var(--success)'}">
                    ${qcData.duplicates.n_exact_duplicates} duplicates flagged
                </span>
            </div>
        `;

        const numericCols = Object.keys(qcData.basic_stats);
        const layoutConfig = {
            paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
            font: { color: '#f3f4f6', family: 'Inter, sans-serif' },
            xaxis: { gridcolor: 'rgba(255,255,255,0.05)', zeroline: false },
            yaxis: { gridcolor: 'rgba(255,255,255,0.05)', zeroline: false }
        };

        if (numericCols.length > 0) {
            if (qcData.basic_stats['coverage']) {
                const covStats = qcData.basic_stats['coverage'];
                const distData = [{
                    x: ['Min Coverage', 'Mean Coverage', 'Max Coverage'],
                    y: [covStats.min, covStats.mean, covStats.max],
                    type: 'bar',
                    marker: { color: ['#ef4444', '#6366f1', '#10b981'], opacity: 0.8 }
                }];
                Plotly.newPlot('plotly-distribution', distData, {
                    ...layoutConfig, title: 'Sequencing Coverage Metrics Summary',
                    yaxis: { ...layoutConfig.yaxis, title: 'Coverage (x)' }
                }, {responsive: true, displayModeBar: false});
            }
            if (qcData.gc_stats.available && qcData.coverage_stats.available) {
                const scatterData = [{
                    x: [qcData.gc_stats.min, qcData.gc_stats.mean, qcData.gc_stats.max],
                    y: [qcData.coverage_stats.min, qcData.coverage_stats.mean, qcData.coverage_stats.max],
                    mode: 'markers+lines+text', type: 'scatter', text: ['Min', 'Mean', 'Max'],
                    textposition: 'top center', marker: { size: 12, color: '#10b981' }
                }];
                Plotly.newPlot('plotly-scatter', scatterData, {
                    ...layoutConfig, title: 'GC Content Range vs Coverage Profile',
                    xaxis: { ...layoutConfig.xaxis, title: 'GC Ratio', tickformat: ',.2%' },
                    yaxis: { ...layoutConfig.yaxis, title: 'Coverage Depth (x)' }
                }, {responsive: true, displayModeBar: false});
            }
        }
    </script>
</body>
</html>
"""

# ===========================================================================
# Proteomics HTML Dashboard Template
# ===========================================================================
PROTEOMICS_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>🧬 Proteomics Quality Control Dashboard</title>
    <!-- Google Fonts -->
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
    <!-- Plotly.js -->
    <script src="https://cdn.plot.ly/plotly-2.24.1.min.js"></script>
    <style>
        :root {
            --bg-color: #080c14;
            --card-bg: rgba(15, 23, 42, 0.7);
            --card-border: rgba(255, 255, 255, 0.06);
            --text-main: #f8fafc;
            --text-muted: #94a3b8;
            --primary: #818cf8;
            --primary-glow: rgba(129, 140, 248, 0.15);
            --success: #34d399;
            --success-glow: rgba(52, 211, 153, 0.15);
            --warning: #fbbf24;
            --warning-glow: rgba(251, 191, 36, 0.15);
            --danger: #f87171;
            --danger-glow: rgba(248, 113, 113, 0.15);
            --font-main: 'Inter', -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            --font-mono: 'JetBrains Mono', monospace;
        }

        * { box-sizing: border-box; margin: 0; padding: 0; }
        body {
            background-color: var(--bg-color);
            color: var(--text-main);
            font-family: var(--font-main);
            line-height: 1.5;
            padding: 2rem;
            min-height: 100vh;
            background-image: radial-gradient(circle at 5% 5%, rgba(129, 140, 248, 0.08) 0%, transparent 40%),
                              radial-gradient(circle at 95% 95%, rgba(52, 211, 153, 0.06) 0%, transparent 40%);
            background-attachment: fixed;
        }
        .container { max-width: 1550px; margin: 0 auto; }
        header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 2rem;
            padding-bottom: 1.5rem;
            border-bottom: 1px solid var(--card-border);
        }
        .brand h1 {
            font-size: 2.2rem;
            font-weight: 800;
            background: linear-gradient(135deg, #cbd5e1 0%, #818cf8 50%, #4f46e5 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }
        .brand p { color: var(--text-muted); font-size: 0.95rem; margin-top: 0.25rem; }
        .timestamp {
            font-family: var(--font-mono);
            font-size: 0.85rem;
            color: var(--text-muted);
            background: var(--card-border);
            padding: 0.5rem 1rem;
            border-radius: 9999px;
            border: 1px solid var(--card-border);
        }
        .summary-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 1.5rem; margin-bottom: 2rem; }
        .stat-card {
            background: var(--card-bg);
            border: 1px solid var(--card-border);
            border-radius: 16px;
            padding: 1.5rem;
            backdrop-filter: blur(12px);
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            position: relative;
            overflow: hidden;
        }
        .stat-card:hover { transform: translateY(-4px); border-color: var(--primary); box-shadow: 0 10px 20px -10px var(--primary-glow); }
        .stat-card.pass:hover { border-color: var(--success); box-shadow: 0 10px 20px -10px var(--success-glow); }
        .stat-card.fail:hover { border-color: var(--danger); box-shadow: 0 10px 20px -10px var(--danger-glow); }
        .stat-card::before { content: ''; position: absolute; top: 0; left: 0; width: 100%; height: 4px; background: var(--primary); }
        .stat-card.pass::before { background: var(--success); }
        .stat-card.fail::before { background: var(--danger); }
        .stat-card.warning::before { background: var(--warning); }
        .stat-label { font-size: 0.85rem; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; color: var(--text-muted); margin-bottom: 0.5rem; }
        .stat-value { font-size: 2.2rem; font-weight: 800; color: var(--text-main); font-family: var(--font-mono); }
        
        .layout-main { display: grid; grid-template-columns: 3fr 2fr; gap: 1.5rem; margin-bottom: 2rem; }
        @media (max-width: 1150px) { .layout-main { grid-template-columns: 1fr; } }
        
        .card { background: var(--card-bg); border: 1px solid var(--card-border); border-radius: 16px; padding: 1.5rem; backdrop-filter: blur(12px); margin-bottom: 1.5rem; }
        .card-title { font-size: 1.25rem; font-weight: 700; margin-bottom: 1rem; display: flex; align-items: center; gap: 0.5rem; border-bottom: 1px solid var(--card-border); padding-bottom: 0.75rem; }
        .chart-container { width: 100%; height: 400px; border-radius: 8px; overflow: hidden; }
        .table-wrapper { max-height: 480px; overflow-y: auto; border-radius: 8px; border: 1px solid var(--card-border); }
        table { width: 100%; border-collapse: collapse; text-align: left; font-size: 0.88rem; }
        th { background: #0f172a; padding: 0.75rem 1rem; font-weight: 600; color: var(--text-main); position: sticky; top: 0; z-index: 10; border-bottom: 2px solid var(--card-border); }
        td { padding: 0.75rem 1rem; border-bottom: 1px solid var(--card-border); color: var(--text-muted); font-family: var(--font-mono); }
        tr:hover td { background: rgba(255, 255, 255, 0.02); color: var(--text-main); }
        
        .search-container { margin-bottom: 1rem; }
        .search-input { width: 100%; padding: 0.75rem 1rem; border-radius: 8px; border: 1px solid var(--card-border); background: rgba(255, 255, 255, 0.03); color: var(--text-main); font-size: 0.9rem; transition: all 0.3s; }
        .search-input:focus { outline: none; border-color: var(--primary); box-shadow: 0 0 0 3px var(--primary-glow); }
        
        .badge { display: inline-flex; align-items: center; padding: 0.25rem 0.6rem; border-radius: 9999px; font-size: 0.75rem; font-weight: 600; text-transform: uppercase; }
        .badge-pass { background: var(--success-glow); color: var(--success); border: 1px solid rgba(52, 211, 153, 0.3); }
        .badge-fail { background: var(--danger-glow); color: var(--danger); border: 1px solid rgba(248, 113, 113, 0.3); animation: pulse-fail 2s infinite; }
        @keyframes pulse-fail {
            0% { box-shadow: 0 0 0 0 rgba(248, 113, 113, 0.4); }
            70% { box-shadow: 0 0 0 6px rgba(248, 113, 113, 0); }
            100% { box-shadow: 0 0 0 0 rgba(248, 113, 113, 0); }
        }
        .flag-tag { background: rgba(251, 191, 36, 0.1); color: var(--warning); border: 1px solid rgba(251, 191, 36, 0.2); padding: 0.1rem 0.4rem; border-radius: 4px; font-size: 0.7rem; margin-right: 0.3rem; display: inline-block; }
    </style>
</head>
<body>
    <div class="container">
        <header>
            <div class="brand">
                <h1>🧬 Proteomics QC Dashboard</h1>
                <p>Interactive DDA (FragPipe) / DIA (Spectronaut) Quality Control Summary</p>
            </div>
            <div class="timestamp" id="report-time">Generated: --</div>
        </header>

        <!-- Metrics Overview -->
        <div class="summary-grid">
            <div class="stat-card">
                <div class="stat-label">Total Proteins</div>
                <div class="stat-value" id="stat-total-proteins">0</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Total Mass Spec Runs</div>
                <div class="stat-value" id="stat-total-samples">0</div>
            </div>
            <div class="stat-card pass">
                <div class="stat-label">Passed QC</div>
                <div class="stat-value" id="stat-passed-samples">0</div>
            </div>
            <div class="stat-card warning">
                <div class="stat-label font-bold text-emerald-400">QC Pass Rate</div>
                <div class="stat-value" id="stat-pass-rate">0.0%</div>
            </div>
        </div>

        <!-- Layout Split -->
        <div class="layout-main">
            <!-- Left Side: Interactive Plots -->
            <div class="plots-column">
                <div class="card">
                    <div class="card-title">🔍 Principal Component Analysis (PCA)</div>
                    <div class="chart-container" id="plotly-pca"></div>
                </div>
                <div class="card">
                    <div class="card-title">🔥 Reproducibility: Pearson Correlation Heatmap</div>
                    <div class="chart-container" id="plotly-corr-heatmap"></div>
                </div>
                <div class="card">
                    <div class="card-title">📊 Quantified Proteins Count per Run</div>
                    <div class="chart-container" id="plotly-ids"></div>
                </div>
            </div>

            <!-- Right Side: Interactive Table of Samples -->
            <div class="table-column">
                <div class="card">
                    <div class="card-title">📋 Run Quality Flags</div>
                    <div class="search-container">
                        <input type="text" id="sample-search" class="search-input" placeholder="🔍 Search run ID, status or flags...">
                    </div>
                    <div class="table-wrapper">
                        <table id="qc-table">
                            <thead>
                                <tr>
                                    <th>Run ID</th>
                                    <th>Status</th>
                                    <th>Detected IDs</th>
                                    <th>Missing %</th>
                                    <th>Violations</th>
                                </tr>
                            </thead>
                            <tbody id="qc-table-body"></tbody>
                        </table>
                    </div>
                </div>
                <div class="card">
                    <div class="card-title">📉 Run Abundance (Log2 Intensity) Profiles</div>
                    <div class="chart-container" id="plotly-intensities"></div>
                </div>
            </div>
        </div>
    </div>

    <script>
        const qcData = {{DATA_PLACEHOLDER}};
        
        // Populate stats
        document.getElementById('report-time').innerText = `Generated: ${new Date().toLocaleString()}`;
        document.getElementById('stat-total-proteins').innerText = qcData.n_proteins.toLocaleString();
        document.getElementById('stat-total-samples').innerText = qcData.n_samples;
        document.getElementById('stat-passed-samples').innerText = qcData.summary.n_pass;
        document.getElementById('stat-pass-rate').innerText = `${(qcData.summary.pass_rate * 100).toFixed(1)}%`;

        // Render searchable table
        const tbody = document.getElementById('qc-table-body');
        function renderTable(filterText = '') {
            tbody.innerHTML = '';
            const filtered = qcData.per_sample_flags.filter(s => {
                const searchLower = filterText.toLowerCase();
                return s.sample_id.toLowerCase().includes(searchLower) || 
                       s.status.toLowerCase().includes(searchLower) ||
                       s.flags.some(f => f.toLowerCase().includes(searchLower));
            });

            filtered.forEach(s => {
                const tr = document.createElement('tr');
                const tdId = document.createElement('td'); tdId.innerText = s.sample_id;
                
                const tdStatus = document.createElement('td');
                const badge = document.createElement('span');
                badge.className = `badge badge-${s.status.toLowerCase()}`;
                badge.innerText = s.status;
                tdStatus.appendChild(badge);

                const tdIds = document.createElement('td');
                tdIds.innerText = s.metrics.proteins_identified.toLocaleString();

                const tdMiss = document.createElement('td');
                tdMiss.innerText = `${(s.metrics.missing_rate * 100).toFixed(1)}%`;

                const tdFlags = document.createElement('td');
                if (s.flags.length === 0) { tdFlags.innerText = '-'; }
                else {
                    s.flags.forEach(f => {
                        const tag = document.createElement('span'); tag.className = 'flag-tag'; tag.innerText = f;
                        tdFlags.appendChild(tag);
                    });
                }

                tr.appendChild(tdId);
                tr.appendChild(tdStatus);
                tr.appendChild(tdIds);
                tr.appendChild(tdMiss);
                tr.appendChild(tdFlags);
                tbody.appendChild(tr);
            });
        }
        document.getElementById('sample-search').addEventListener('input', (e) => { renderTable(e.target.value); });
        renderTable();

        // Common layout config
        const layoutConfig = {
            paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
            font: { color: '#e2e8f0', family: 'Inter, sans-serif' },
            xaxis: { gridcolor: 'rgba(255,255,255,0.05)', zeroline: false },
            yaxis: { gridcolor: 'rgba(255,255,255,0.05)', zeroline: false }
        };

        // 1. PCA Scatter Plot
        const pcaTrace = {
            x: qcData.pca.PC1,
            y: qcData.pca.PC2,
            mode: 'markers+text',
            type: 'scatter',
            text: qcData.per_sample_flags.map(s => s.sample_id),
            textposition: 'top center',
            marker: {
                size: 14,
                color: qcData.per_sample_flags.map(s => s.status === 'PASS' ? '#34d399' : '#f87171'),
                line: { color: '#0f172a', width: 2 }
            }
        };
        Plotly.newPlot('plotly-pca', [pcaTrace], {
            ...layoutConfig,
            title: 'Sample Clustering Profile',
            xaxis: { ...layoutConfig.xaxis, title: 'Principal Component 1' },
            yaxis: { ...layoutConfig.yaxis, title: 'Principal Component 2' }
        }, {responsive: true, displayModeBar: false});

        // 2. Correlation Heatmap
        const heatmapTrace = {
            z: qcData.correlation.grid,
            x: qcData.correlation.samples,
            y: qcData.correlation.samples,
            type: 'heatmap',
            colorscale: [
                [0.0, '#312e81'], // Indigo
                [0.5, '#4f46e5'],
                [0.8, '#818cf8'],
                [1.0, '#34d399']  // Emerald green for high correlation
            ]
        };
        Plotly.newPlot('plotly-corr-heatmap', [heatmapTrace], {
            ...layoutConfig,
            title: 'Pairwise Pearson Correlation (Reproducibility)'
        }, {responsive: true, displayModeBar: false});

        // 3. Protein Identifications Bar Chart
        const runs = qcData.per_sample_flags.map(s => s.sample_id);
        const idsTrace = {
            x: runs,
            y: qcData.per_sample_flags.map(s => s.metrics.proteins_identified),
            type: 'bar',
            marker: {
                color: qcData.per_sample_flags.map(s => s.status === 'PASS' ? '#818cf8' : '#f87171'),
                opacity: 0.8
            }
        };
        Plotly.newPlot('plotly-ids', [idsTrace], {
            ...layoutConfig,
            title: 'Total Identified Proteins per Run',
            yaxis: { ...layoutConfig.yaxis, title: 'Protein Group Count' }
        }, {responsive: true, displayModeBar: false});

        // 4. Intensities Distribution Profile (Box plot summaries)
        const boxTraces = [];
        runs.forEach(run => {
            const stats = qcData.intensity_stats[run];
            boxTraces.push({
                y: [stats.min, stats.q25, stats.median, stats.q75, stats.max],
                name: run,
                type: 'box',
                boxpoints: false,
                marker: { color: qcData.per_sample_flags.find(s => s.sample_id === run).status === 'PASS' ? '#34d399' : '#f87171' }
            });
        });
        Plotly.newPlot('plotly-intensities', boxTraces, {
            ...layoutConfig,
            title: 'Log2 Protein Abundance Distribution per Run',
            yaxis: { ...layoutConfig.yaxis, title: 'Log2 Abundance' },
            showlegend: false
        }, {responsive: true, displayModeBar: false});
"""

# ===========================================================================
# FragPipe Advanced QC HTML Dashboard Template
# ===========================================================================
FRAGPIPE_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>🧬 FragPipe Advanced Quality Control Dashboard</title>
    <!-- Google Fonts -->
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
    <!-- Plotly.js -->
    <script src="https://cdn.plot.ly/plotly-2.24.1.min.js"></script>
    <style>
        :root {
            --bg-color: #0b0f19;
            --card-bg: rgba(17, 24, 39, 0.75);
            --card-border: rgba(255, 255, 255, 0.08);
            --text-main: #f3f4f6;
            --text-muted: #9ca3af;
            --primary: #818cf8;
            --primary-glow: rgba(129, 140, 248, 0.2);
            --success: #34d399;
            --success-glow: rgba(52, 211, 153, 0.2);
            --warning: #fbbf24;
            --danger: #f87171;
            --font-main: 'Inter', -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            --font-mono: 'JetBrains Mono', monospace;
        }

        * { box-sizing: border-box; margin: 0; padding: 0; }
        body {
            background-color: var(--bg-color);
            color: var(--text-main);
            font-family: var(--font-main);
            line-height: 1.5;
            padding: 2rem;
            min-height: 100vh;
            background-image: radial-gradient(circle at 10% 20%, rgba(129, 140, 248, 0.06) 0%, transparent 40%),
                              radial-gradient(circle at 90% 80%, rgba(52, 211, 153, 0.05) 0%, transparent 40%);
            background-attachment: fixed;
        }
        .container { max-width: 1450px; margin: 0 auto; }
        header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 2rem;
            padding-bottom: 1.5rem;
            border-bottom: 1px solid var(--card-border);
        }
        .brand h1 {
            font-size: 2.2rem;
            font-weight: 800;
            background: linear-gradient(135deg, #cbd5e1 0%, #818cf8 50%, #4f46e5 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }
        .brand p { color: var(--text-muted); font-size: 0.95rem; margin-top: 0.25rem; }
        .timestamp {
            font-family: var(--font-mono);
            font-size: 0.85rem;
            color: var(--text-muted);
            background: var(--card-border);
            padding: 0.5rem 1rem;
            border-radius: 9999px;
            border: 1px solid var(--card-border);
        }
        .summary-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr)); gap: 1.5rem; margin-bottom: 2rem; }
        .stat-card {
            background: var(--card-bg);
            border: 1px solid var(--card-border);
            border-radius: 16px;
            padding: 1.5rem;
            backdrop-filter: blur(12px);
            transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
            position: relative;
            overflow: hidden;
        }
        .stat-card:hover { transform: translateY(-4px); border-color: var(--primary); box-shadow: 0 10px 20px -10px var(--primary-glow); }
        .stat-card::before { content: ''; position: absolute; top: 0; left: 0; width: 100%; height: 4px; background: var(--primary); }
        .stat-card.success::before { background: var(--success); }
        .stat-card.warning::before { background: var(--warning); }
        .stat-card.danger::before { background: var(--danger); }
        .stat-label { font-size: 0.85rem; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; color: var(--text-muted); margin-bottom: 0.5rem; }
        .stat-value { font-size: 2.2rem; font-weight: 800; color: var(--text-main); font-family: var(--font-mono); }
        .stat-sub { font-size: 0.78rem; color: var(--text-muted); margin-top: 0.3rem; }
        
        .layout-main { display: grid; grid-template-columns: 1fr 1fr; gap: 1.5rem; margin-bottom: 2rem; }
        @media (max-width: 1024px) { .layout-main { grid-template-columns: 1fr; } }
        
        .card { background: var(--card-bg); border: 1px solid var(--card-border); border-radius: 16px; padding: 1.5rem; backdrop-filter: blur(12px); margin-bottom: 1.5rem; }
        .card-title { font-size: 1.25rem; font-weight: 700; margin-bottom: 1rem; display: flex; align-items: center; gap: 0.5rem; border-bottom: 1px solid var(--card-border); padding-bottom: 0.75rem; }
        .chart-container { width: 100%; height: 380px; border-radius: 8px; overflow: hidden; }
        
        .metrics-list { display: flex; flex-direction: column; gap: 0.85rem; }
        .metric-row { display: flex; justify-content: space-between; align-items: center; padding-bottom: 0.6rem; border-bottom: 1px solid rgba(255,255,255,0.03); }
        .metric-name { font-size: 0.92rem; color: var(--text-main); font-weight: 500; }
        .metric-val { font-family: var(--font-mono); font-size: 0.95rem; color: var(--primary); font-weight: 600; }
    </style>
</head>
<body>
    <div class="container">
        <header>
            <div class="brand">
                <h1>🧬 FragPipe QC Dashboard</h1>
                <p>Advanced Quality Control metrics & statistics for DDA MS proteomics searches</p>
            </div>
            <div class="timestamp" id="report-time">Generated: --</div>
        </header>

        <!-- KPI Grid -->
        <div class="summary-grid">
            <div class="stat-card">
                <div class="stat-label">Identified Proteins</div>
                <div class="stat-value" id="stat-prots">0</div>
                <div class="stat-sub" id="sub-prots">Contaminants: 0</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Identified Peptides</div>
                <div class="stat-value" id="stat-peps">0</div>
                <div class="stat-sub" id="sub-peps">Contaminant Peps: 0</div>
            </div>
            <div class="stat-card success">
                <div class="stat-label">Total PSMs</div>
                <div class="stat-value" id="stat-psms">0</div>
                <div class="stat-sub" id="sub-psms">ID Rate: 0%</div>
            </div>
            <div class="stat-card warning">
                <div class="stat-label">Missed Cleavages</div>
                <div class="stat-value" id="stat-missed">0%</div>
                <div class="stat-sub" id="sub-missed">Expl. Intensity: 0%</div>
            </div>
        </div>

        <div class="layout-main">
            <!-- Left Side: Modification & Dev charts -->
            <div class="chart-column">
                <div class="card">
                    <div class="card-title">📈 Peptide Mass Deviations (ppm)</div>
                    <div class="chart-container" id="plotly-mass-dev"></div>
                </div>
                <div class="card">
                    <div class="card-title">🔬 Assigned Modifications Breakdown</div>
                    <div class="chart-container" id="plotly-mods"></div>
                </div>
            </div>

            <!-- Right Side: Detailed Metrics Table -->
            <div class="details-column">
                <div class="card">
                    <div class="card-title">📋 Comprehensive QC Summary</div>
                    <div class="metrics-list" id="metrics-table">
                        <!-- Populated by JavaScript -->
                    </div>
                </div>
            </div>
        </div>
    </div>

    <script>
        const qcData = {{DATA_PLACEHOLDER}};
        
        document.getElementById('report-time').innerText = `Generated: ${new Date().toLocaleString()}`;
        document.getElementById('stat-prots').innerText = (qcData.nProts || 0).toLocaleString();
        document.getElementById('sub-prots').innerText = `Contaminants: ${qcData.nContaProts || 0} (${((qcData.nContaProts || 0) / (qcData.nProts || 1) * 100).toFixed(2)}%)`;
        
        document.getElementById('stat-peps').innerText = (qcData.nPeps || 0).toLocaleString();
        document.getElementById('sub-peps').innerText = qcData.nContaPeps ? `Contaminants: ${qcData.nContaPeps.toLocaleString()}` : 'Contaminant data N/A';
        
        document.getElementById('stat-psms').innerText = (qcData.nPsms || 0).toLocaleString();
        document.getElementById('sub-psms').innerText = qcData.idRate ? `ID Rate: ${qcData.idRate}%` : 'ID Rate: N/A';
        
        document.getElementById('stat-missed').innerText = qcData.missCl ? `${qcData.missCl}%` : 'N/A';
        document.getElementById('sub-missed').innerText = qcData.explIons ? `Explained Intensity: ${qcData.explIons}%` : 'Explained Intensity: N/A';

        // Populate Comprehensive QC Summary Table
        const metricsTable = document.getElementById('metrics-table');
        const rows = [
            { name: "Total Proteins", val: qcData.nProts },
            { name: "Contaminant Proteins", val: qcData.nContaProts },
            { name: "Total Peptides", val: qcData.nPeps },
            { name: "Contaminant Peptides", val: qcData.nContaPeps },
            { name: "Total PSMs", val: qcData.nPsms },
            { name: "PSM ID Rate", val: qcData.idRate ? `${qcData.idRate}%` : 'N/A' },
            { name: "Missed Cleavage Rate", val: qcData.missCl ? `${qcData.missCl}%` : 'N/A' },
            { name: "Total Peptide Intensity", val: qcData.totPepI ? qcData.totPepI.toExponential(3) : 'N/A' },
            { name: "Contaminant Peptide Intensity", val: qcData.contaPepI ? qcData.contaPepI.toExponential(3) : 'N/A' },
            { name: "Explained Intensity Ratio", val: qcData.explIons ? `${qcData.explIons}%` : 'N/A' }
        ];
        
        rows.forEach(r => {
            if (r.val !== undefined && r.val !== null) {
                const div = document.createElement('div');
                div.className = 'metric-row';
                div.innerHTML = `<span class="metric-name">${r.name}</span><span class="metric-val">${r.val}</span>`;
                metricsTable.appendChild(div);
            }
        });

        // Common layout config
        const layoutConfig = {
            paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
            font: { color: '#f3f4f6', family: 'Inter, sans-serif' },
            xaxis: { gridcolor: 'rgba(255,255,255,0.05)', zeroline: false },
            yaxis: { gridcolor: 'rgba(255,255,255,0.05)', zeroline: false }
        };

        // 1. Mass Deviation box plot/summary
        if (qcData.medMassDev !== undefined) {
            const devData = [{
                x: ['q05', 'q10', 'q25', 'median', 'q75', 'q90', 'q95'],
                y: [qcData.q05MassDev, qcData.q10MassDev, qcData.q25MassDev, qcData.medMassDev, qcData.q75MassDev, qcData.q90MassDev, qcData.q95MassDev],
                type: 'bar',
                marker: {
                    color: '#818cf8', opacity: 0.85
                }
            }];
            Plotly.newPlot('plotly-mass-dev', devData, {
                ...layoutConfig,
                title: 'Mass Deviation distribution quantiles (ppm)',
                xaxis: { ...layoutConfig.xaxis, title: 'Quantile' },
                yaxis: { ...layoutConfig.yaxis, title: 'Mass Dev (ppm)' }
            }, {responsive: true, displayModeBar: false});
        }

        // 2. Modifications Breakdown
        const modKeys = ['Oxidation (M)', 'Carbamidomethyl', 'Methylthio', 'Carbamyl (K)', 'Carbamyl (N-term)'];
        const modVals = [qcData.oxiPsms || 0, qcData.carbamidoPsms || 0, qcData.meththPsms || 0, qcData.carbamylKpsms || 0, qcData.carbamylNtPsms || 0];
        
        const modsTrace = {
            x: modKeys,
            y: modVals,
            type: 'bar',
            marker: {
                color: '#34d399', opacity: 0.85
            }
        };
        Plotly.newPlot('plotly-mods', [modsTrace], {
            ...layoutConfig,
            title: 'Assigned Modifications counts per PSM',
            yaxis: { ...layoutConfig.yaxis, title: 'PSM Count' }
        }, {responsive: true, displayModeBar: false});
    </script>
</body>
</html>
"""


def generate_json_report(results: dict[str, Any], output_path: Path) -> None:
    """Save raw QC results to a JSON file."""
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    logger.info("Saved JSON report to %s", output_path)


def generate_html_report(results: dict[str, Any], output_path: Path) -> None:
    """Save an interactive HTML dashboard report."""
    data_json_str = json.dumps(results, indent=2, ensure_ascii=False)
    
    # Choose template based on data domain
    if results.get("is_fragpipe_qc", False):
        template = FRAGPIPE_HTML_TEMPLATE
    elif results.get("is_proteomics", False):
        template = PROTEOMICS_HTML_TEMPLATE
    else:
        template = GENOMICS_HTML_TEMPLATE
    
    html_content = template.replace("{{DATA_PLACEHOLDER}}", data_json_str)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_content)
    logger.info("Saved HTML report to %s", output_path)


def generate_text_report(results: dict[str, Any]) -> str:
    """Return a plain text table summary suitable for console stdout."""
    if results.get("is_fragpipe_qc", False):
        report = []
        report.append("=" * 65)
        report.append("🧬 FRAGPIPE PROTEOMICS QUALITY CONTROL REPORT SUMMARY")
        report.append("=" * 65)
        report.append(f"Total Identified Proteins : {results.get('nProts', 'N/A')}")
        report.append(f"Contaminant Proteins       : {results.get('nContaProts', 'N/A')}")
        report.append(f"Total Identified Peptides : {results.get('nPeps', 'N/A')}")
        report.append(f"Total PSMs                : {results.get('nPsms', 'N/A')}")
        report.append(f"PSM ID Rate               : {results.get('idRate', 'N/A')}%")
        report.append(f"Missed Cleavage Rate      : {results.get('missCl', 'N/A')}%")
        if 'totPepI' in results:
            report.append(f"Total Peptide Intensity   : {results.get('totPepI', 'N/A'):.2e}")
        if 'explIons' in results:
            report.append(f"Explained Ions (%)        : {results.get('explIons', 'N/A')}%")
        report.append("-" * 65)
        
        report.append("\n📋 ASSIGNED MODIFICATIONS BREAKDOWN (PSMs):")
        report.append(f"  • Oxidation (Met) [oxiPsms]       : {results.get('oxiPsms', 0)}")
        report.append(f"  • Carbamidomethyl [carbamidoPsms]  : {results.get('carbamidoPsms', 0)}")
        report.append(f"  • Methylthio [meththPsms]          : {results.get('meththPsms', 0)}")
        report.append(f"  • Carbamyl (K) [carbamylKpsms]     : {results.get('carbamylKpsms', 0)}")
        report.append(f"  • Carbamyl (N-term) [carbamylNt]   : {results.get('carbamylNtPsms', 0)}")
        report.append("-" * 65)

        if 'medMassDev' in results:
            report.append("\n📈 PEPTIDE MASS DEVIATIONS (PPM):")
            report.append(f"  • Median Mass Dev                 : {results.get('medMassDev', 'N/A')} ppm")
            report.append(f"  • 25th - 75th Quantiles           : {results.get('q25MassDev', 'N/A')} to {results.get('q75MassDev', 'N/A')} ppm")
            report.append(f"  • 5th - 95th Quantiles             : {results.get('q05MassDev', 'N/A')} to {results.get('q95MassDev', 'N/A')} ppm")
            report.append("=" * 65)
        return "\n".join(report)

    summary = results["summary"]
    pass_rate = summary["pass_rate"] * 100
    is_proteomics = results.get("is_proteomics", False)

    report = []
    report.append("=" * 65)
    if is_proteomics:
        report.append("🧬 PROTEOMICS SAMPLE QUALITY CONTROL REPORT SUMMARY")
    else:
        report.append("🧬 GENOMICS SAMPLE QUALITY CONTROL REPORT SUMMARY")
    report.append("=" * 65)
    report.append(f"Total Samples / Runs   : {results['n_samples']}")
    if is_proteomics:
        report.append(f"Total Proteins Quant   : {results['n_proteins']}")
    report.append(f"Passed QC Filter       : {summary['n_pass']}")
    report.append(f"Failed QC Filter       : {summary['n_fail']}")
    report.append(f"QC Pass Rate           : {pass_rate:.2f}%")
    report.append("-" * 65)

    if is_proteomics:
        report.append("\n📋 INDIVIDUAL MASS SPECTROMETRY RUN DETAILS:")
        report.append(
            f"{'Run/Sample ID':<18} {'Status':<8} {'Quant IDs':<12} {'Missing %':<12} {'Median (Log2)':<15}"
        )
        report.append("-" * 65)
        for s in results["per_sample_flags"]:
            metrics = s["metrics"]
            median_val = metrics["median_log2_abundance"]
            report.append(
                f"{s['sample_id']:<18} "
                f"{s['status']:<8} "
                f"{metrics['proteins_identified']:<12} "
                f"{metrics['missing_rate']*100:<11.1f}% "
                f"{median_val:<15.4f}"
            )
    else:
        # Basic stats table for genomics
        report.append("\n📈 KEY METRICS DESCRIPTIVE STATISTICS:")
        report.append(
            f"{'Metric':<16} {'Count':<8} {'Mean':<10} {'Median':<10} {'StdDev':<10}"
        )
        report.append("-" * 65)
        for col, stats_dict in results["basic_stats"].items():
            report.append(
                f"{col:<16} "
                f"{stats_dict['count']:<8} "
                f"{stats_dict['mean']:<10.4f} "
                f"{stats_dict['median']:<10.4f} "
                f"{stats_dict['std']:<10.4f}"
            )
    report.append("-" * 65)

    # Flagged sample warnings
    failed_samples = [
        s for s in results["per_sample_flags"] if s["status"] == "FAIL"
    ]
    if failed_samples:
        report.append(f"\n⚠️  FLAGGED QC VIOLATIONS ({len(failed_samples)}):")
        for s in failed_samples:
            flags_str = ", ".join(s["flags"])
            report.append(f"  • Sample {s['sample_id']}: [{flags_str}]")
    else:
        report.append("\n🎉 All samples passed Quality Control standards.")

    report.append("=" * 65)
    return "\n".join(report)


def generate_report(
    results: dict[str, Any],
    output_path: str | Path | None = None,
    fmt: str = "json",
) -> str | None:
    """
    Main entry point for generating Sample QC reports.
    """
    if fmt == "text" and output_path is None:
        return generate_text_report(results)

    if output_path is None:
        raise ValueError(
            f"An output_path must be supplied to save format '{fmt}'"
        )

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
        logger.info("Saved Text report to %s", out_p)
    else:
        raise ValueError(
            f"Invalid report format '{fmt}'. Choose from: 'json', 'html', 'text'."
        )
    return None

