#!/usr/bin/env python3
"""Generate a standalone HTML summary report for qChIP-MS.

Embeds the LFQ volcano plot (base64 encoded) alongside a formatted table
of top enriched chromatin-binding protein hits and summary metrics.
"""

from __future__ import annotations

import base64
import sys
from pathlib import Path
import pandas as pd


def encode_image_base64(image_path: Path) -> str:
    """Encode image to base64 data URI for self-contained HTML embedding."""
    if not image_path.exists():
        return ""
    try:
        data = base64.b64encode(image_path.read_bytes()).decode("utf-8")
        return f"data:image/png;base64,{data}"
    except Exception as exc:
        print(f"Warning: Could not encode image {image_path}: {exc}")
        return ""


def main() -> None:
    if len(sys.argv) < 4:
        print(
            "Usage: summary_report.py <enriched_tsv> <annotated_bed> [volcano_png] <html_report>"
        )
        sys.exit(1)

    enriched_tsv = Path(sys.argv[1])
    annotated_bed = Path(sys.argv[2])

    if len(sys.argv) == 5:
        volcano_png = Path(sys.argv[3])
        html_report = Path(sys.argv[4])
    else:
        volcano_png = Path("results/enrichment/volcano_plot.png")
        html_report = Path(sys.argv[3])

    # 1. Parse Enriched Proteins TSV
    total_proteins = 0
    enriched_count = 0
    top_hits_table_html = "<p>No enriched proteins data available.</p>"

    try:
        if enriched_tsv.exists():
            df = pd.read_csv(enriched_tsv, sep="\t", usecols = ["Gene_Name", "Protein_ID", "Log2FC", "PValue", "FDR", "IsEnriched"])
            total_proteins = len(df)
            if "IsEnriched" in df.columns:
                enriched_df = df[df["IsEnriched"]].copy()
            else:
                enriched_df = df.copy()
            enriched_count = len(enriched_df)

            # Build HTML table for top significant hits
            if not enriched_df.empty:
                sort_col = "Log2FC" if "Log2FC" in enriched_df.columns else enriched_df.columns[0]
                top_df = enriched_df.sort_values(by=sort_col, ascending=False).head(15)

                cols = [
                    c for c in ["Gene_Name", "Protein_ID", "Log2FC", "PValue", "FDR"]
                    if c in top_df.columns
                ]
                if not cols:
                    cols = list(top_df.columns[:5])

                rows = []
                for _, r in top_df[cols].iterrows():
                    tds = []
                    for c in cols:
                        val = r[c]
                        if isinstance(val, float):
                            tds.append(f"<td>{val:.4f}</td>")
                        elif c == "Gene_Name":
                            tds.append(f"<td><strong>{val}</strong></td>")
                        else:
                            tds.append(f"<td>{val}</td>")
                    rows.append(f"<tr>{''.join(tds)}</tr>")

                headers = "".join(f"<th>{c}</th>" for c in cols)
                top_hits_table_html = f"""
                <div class="table-card">
                    <table>
                        <thead>
                            <tr>{headers}</tr>
                        </thead>
                        <tbody>
                            {''.join(rows)}
                        </tbody>
                    </table>
                </div>
                """
    except Exception as e:
        print(f"Error reading TSV: {e}")

    # 2. Parse Annotated BED File (Deduplicate using set)
    annotated_peaks = 0
    try:
        if annotated_bed.exists():
            with annotated_bed.open("r", encoding="utf-8", buffering=1024 * 1024) as f:
                unique_peaks = set()
                for line in f:
                    line_str = line.strip()
                    if line_str and not line_str.startswith("#"):
                        parts = line_str.split("\t", 3)
                        if len(parts) >= 3:
                            unique_peaks.add((parts[0], parts[1], parts[2]))
                annotated_peaks = len(unique_peaks)
    except Exception as e:
        print(f"Error reading BED: {e}")

    # 3. Base64 Encode Volcano Plot Image
    volcano_uri = encode_image_base64(volcano_png)
    if volcano_uri:
        plot_html = f'<div class="plot-card"><img src="{volcano_uri}" alt="Volcano Plot" class="volcano-img"/></div>'
    else:
        plot_html = '<div class="plot-card placeholder"><p>Volcano plot image not found.</p></div>'

    # 4. Generate Modern Responsive HTML Page
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>qChIP-MS Analysis Summary Report</title>
    <style>
        :root {{
            --bg: #0f172a;
            --surface: #1e293b;
            --border: #334155;
            --text: #f8fafc;
            --muted: #94a3b8;
            --accent: #38bdf8;
            --success: #34d399;
        }}
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            background: var(--bg);
            color: var(--text);
            padding: 40px 20px;
            line-height: 1.6;
        }}
        .container {{
            max-width: 1100px;
            margin: 0 auto;
        }}
        header {{
            margin-bottom: 28px;
            border-bottom: 1px solid var(--border);
            padding-bottom: 18px;
        }}
        h1 {{
            font-size: 1.9rem;
            color: var(--text);
            margin-bottom: 6px;
        }}
        .subtitle {{
            color: var(--muted);
            font-size: 0.95rem;
        }}
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
            gap: 18px;
            margin-bottom: 32px;
        }}
        .metric-card {{
            background: var(--surface);
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 20px;
            text-align: center;
        }}
        .metric-value {{
            font-size: 2.2rem;
            font-weight: 700;
            color: var(--accent);
            margin-bottom: 2px;
        }}
        .metric-value.success {{
            color: var(--success);
        }}
        .metric-label {{
            color: var(--muted);
            font-size: 0.85rem;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        .grid-2col {{
            display: grid;
            grid-template-columns: 1.1fr 1fr;
            gap: 24px;
            align-items: start;
            margin-bottom: 32px;
        }}
        @media (max-width: 820px) {{
            .grid-2col {{ grid-template-columns: 1fr; }}
        }}
        .section-title {{
            font-size: 1.25rem;
            margin-bottom: 12px;
            color: var(--text);
        }}
        .plot-card {{
            background: var(--surface);
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 14px;
            text-align: center;
        }}
        .volcano-img {{
            max-width: 100%;
            height: auto;
            border-radius: 6px;
            background: #ffffff;
        }}
        .placeholder {{
            padding: 50px 20px;
            color: var(--muted);
        }}
        .table-card {{
            background: var(--surface);
            border: 1px solid var(--border);
            border-radius: 10px;
            overflow-x: auto;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 0.88rem;
            text-align: left;
        }}
        th {{
            background: #243248;
            color: var(--muted);
            padding: 10px 14px;
            font-weight: 600;
            border-bottom: 1px solid var(--border);
        }}
        td {{
            padding: 10px 14px;
            border-bottom: 1px solid var(--border);
            color: var(--text);
        }}
        tr:last-child td {{
            border-bottom: none;
        }}
        tr:hover td {{
            background: rgba(255, 255, 255, 0.02);
        }}
        footer {{
            margin-top: 40px;
            text-align: center;
            color: var(--muted);
            font-size: 0.85rem;
            border-top: 1px solid var(--border);
            padding-top: 16px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>qChIP-MS Analysis Summary Report</h1>
            <p class="subtitle">Quantitative Chromatin Immunoprecipitation Mass Spectrometry Pipeline | BDB-Genomics</p>
        </header>

        <div class="metrics-grid">
            <div class="metric-card">
                <div class="metric-value">{total_proteins}</div>
                <div class="metric-label">Total Proteins Analyzed</div>
            </div>
            <div class="metric-card">
                <div class="metric-value success">{enriched_count}</div>
                <div class="metric-label">Significantly Enriched</div>
            </div>
            <div class="metric-card">
                <div class="metric-value">{annotated_peaks}</div>
                <div class="metric-label">Annotated Locus Peaks</div>
            </div>
        </div>

        <div class="grid-2col">
            <div>
                <h2 class="section-title">Enrichment Volcano Plot</h2>
                {plot_html}
            </div>
            <div>
                <h2 class="section-title">Top Enriched Proteins</h2>
                {top_hits_table_html}
            </div>
        </div>

        <footer>
            Generated automatically by <strong>BDB-Genomics / bdb-qchipms</strong> | Snakemake Pipeline
        </footer>
    </div>
</body>
</html>"""

    html_report.parent.mkdir(parents=True, exist_ok=True)
    with html_report.open("w", encoding="utf-8") as f:
        f.write(html_content)

    print(f"Report generated successfully at {html_report}")


if __name__ == "__main__":
    main()
