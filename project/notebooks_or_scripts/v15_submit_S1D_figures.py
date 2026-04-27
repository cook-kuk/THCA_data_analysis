#!/usr/bin/env python3
"""v15 S1-D: re-render the 5 v15 figures from their plotly HTML
exports into static PDF + PNG (NeurIPS submission requirement).

Approach: each plotly HTML file embeds the figure JSON as a
JavaScript Plotly.newPlot() call. We extract that JSON and recreate
the Figure object via plotly.io.from_json, then write_image() it.

Outputs:
  reports/v15_neurips/submit/figures/figure[1-5].pdf
  reports/v15_neurips/submit/figures/figure[1-5].png
  reports/v15_neurips/submit/figures/figure_export_audit.md
"""
from __future__ import annotations
import json
import re
from pathlib import Path
import sys

ROOT = Path("/opt/thyroid-dash/project")
SRC  = ROOT / "reports" / "html" / "figs_interactive" / "v15"
OUT  = ROOT / "reports" / "v15_neurips" / "submit" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# Map figure numbers (paper section) to source HTML
FIGURES = [
    ("figure1_theorem2", "theorem2_decomposition.html",       "Fig 1 — §5.1 Theorem 2 decomposition"),
    ("figure2_stress",   "stress_test_phase_diagram.html",    "Fig 2 — §5.2 stress test phase"),
    ("figure3_tta",      "tta_benchmark.html",                "Fig 3 — §5.3 TTA benchmark"),
    ("figure4_scaling",  "scaling_plot.html",                 "Fig 4 — §5.4 foundation-model scaling"),
    ("figure5_xdomain",  "cross_domain_heatmap.html",         "Fig 5 — §5.5 cross-domain"),
]

import plotly.io as pio
import plotly.graph_objects as go

audit_lines = ["# v15 figure export audit\n",
               "_Generated S1-D, 2026-04-27._\n"]

ok_count = 0
for outname, src_name, caption in FIGURES:
    src_path = SRC / src_name
    if not src_path.exists():
        audit_lines.append(f"- ❌ **{outname}** — source missing: `{src_path}`")
        continue
    html = src_path.read_text()
    # Plotly's newPlot call: Plotly.newPlot(  "div-id",  [traces], {layout}, {config}  )
    # We extract the [traces] and {layout} via regex.
    m = re.search(r'Plotly\.newPlot\(\s*"[^"]+",\s*(\[.*?\])\s*,\s*(\{.*?\})\s*,\s*\{"responsive":\s*true\}\s*\)',
                  html, re.DOTALL)
    if not m:
        # Fall back: try without config block
        m = re.search(r'Plotly\.newPlot\(\s*"[^"]+",\s*(\[.*?\])\s*,\s*(\{.*?\})\s*\)',
                      html, re.DOTALL)
    if not m:
        audit_lines.append(f"- ❌ **{outname}** — could not parse Plotly.newPlot from `{src_name}`")
        continue
    try:
        traces_json = m.group(1)
        layout_json = m.group(2)
        traces = json.loads(traces_json)
        layout = json.loads(layout_json)
        # Strip Plotly's bdata/dtype hex-encoded arrays — they break json
        # roundtrip via plotly.io.from_json. Convert via from_json instead.
        fig_json = json.dumps({"data": traces, "layout": layout})
        fig = pio.from_json(fig_json)
    except Exception as e:
        audit_lines.append(f"- ❌ **{outname}** — JSON parse error: {type(e).__name__}: {e}")
        continue
    # Write PDF and PNG
    pdf_path = OUT / f"{outname}.pdf"
    png_path = OUT / f"{outname}.png"
    try:
        fig.write_image(str(pdf_path), engine="kaleido", width=900, height=600)
        fig.write_image(str(png_path), engine="kaleido", width=1800, height=1200, scale=1)
        sz = pdf_path.stat().st_size
        audit_lines.append(f"- ✅ **{outname}** ({caption}) — PDF {sz} bytes, PNG {png_path.stat().st_size} bytes")
        ok_count += 1
    except Exception as e:
        audit_lines.append(f"- ❌ **{outname}** — write_image failed: {type(e).__name__}: {e}")

audit_lines.insert(2, f"**Result: {ok_count}/{len(FIGURES)} figures exported.**\n")
(OUT / "figure_export_audit.md").write_text("\n".join(audit_lines))
print(f"S1-D exported {ok_count}/{len(FIGURES)} figures")
sys.exit(0 if ok_count == len(FIGURES) else 1)
