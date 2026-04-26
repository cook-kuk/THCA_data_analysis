#!/usr/bin/env python3
"""Re-render reports/html/figs_interactive/v6/multi_scale_dial.html from the
current multi_resolution_dial.tsv. Bulk row red, computed pseudobulk rows
gold/orange, pending GPU rows gray.
"""
from __future__ import annotations

from pathlib import Path
import math

import pandas as pd
import plotly.graph_objects as go

PROJECT = Path("/opt/thyroid-dash/project")
TSV = PROJECT / "results" / "v6_scrna" / "dial_scrna" / "multi_resolution_dial.tsv"
OUT = PROJECT / "reports" / "html" / "figs_interactive" / "v6" / "multi_scale_dial.html"

DISPLAY = {
    "bulk_v5p1_THCA_LogReg_l2":          "Bulk<br>(v5.1)",
    "per_celltype_malignant_pseudobulk": "Malignant<br>pseudobulk",
    "per_substate_pseudobulk":           "3-substate<br>pseudobulk",
    "scgpt_cell_embedding":              "scGPT<br>embedding",
    "vega_pathway_embedding":            "VEGA<br>pathway",
}

# style buckets
COLOR_BULK    = "#E74C3C"   # red
COLOR_COMPUTED = "#F5A623"  # gold/orange
COLOR_PENDING = "#444"      # gray


def main():
    df = pd.read_csv(TSV, sep="\t")
    order = list(DISPLAY.keys())
    df["_o"] = df["layer"].apply(lambda x: order.index(x) if x in order else 99)
    df = df.sort_values("_o").reset_index(drop=True)

    xs = [DISPLAY.get(l, l) for l in df["layer"]]
    ys = []
    texts = []
    colors = []
    for _, row in df.iterrows():
        v = row["dial"]
        try:
            v_f = float(v)
            has = not math.isnan(v_f)
        except (TypeError, ValueError):
            v_f = 0.0
            has = False
        ys.append(v_f if has else 0.0)
        if has:
            texts.append(f"{v_f:.3f}")
        else:
            texts.append("pending")
        if row["layer"] == "bulk_v5p1_THCA_LogReg_l2":
            colors.append(COLOR_BULK)
        elif has:
            colors.append(COLOR_COMPUTED)
        else:
            colors.append(COLOR_PENDING)

    fig = go.Figure(
        data=[go.Bar(
            x=xs, y=ys, text=texts, textposition="outside",
            hovertemplate="%{x}<br>DIAL: %{text}<extra></extra>",
            marker=dict(color=colors, line=dict(color="#F5A623", width=1)),
        )]
    )
    fig.update_layout(
        title="Multi-scale DIAL: bulk vs single-cell foundation models",
        yaxis_title="DIAL",
        xaxis_title="",
        height=420,
        margin=dict(l=50, r=20, t=60, b=80),
        template="plotly_white",
    )
    fig.update_yaxes(range=[0, max(0.6, max(ys) * 1.25 if ys else 0.6)])

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(OUT), include_plotlyjs="cdn", full_html=True)
    print(f"[write] {OUT}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
