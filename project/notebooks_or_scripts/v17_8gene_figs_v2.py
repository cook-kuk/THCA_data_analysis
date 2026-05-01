#!/usr/bin/env python3
"""Re-make the 4 hardest-to-read figures with bigger fonts + cleaner content."""
from __future__ import annotations
from pathlib import Path
import numpy as np, pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

FIG = Path("/opt/thyroid-dash/project/reports/html/figs_interactive/v17")
DARK = dict(template="plotly_dark", paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)", font=dict(color="#EAEAEA", size=14))
GENES = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]

# ============================================================
# A · Venn (cleaner, larger fonts, less overlap)
# ============================================================
fig = go.Figure()
fig.add_shape(type="circle", x0=0, y0=0, x1=5, y1=5,
              line=dict(color="#F5A623", width=4),
              fillcolor="rgba(245,166,35,0.18)")
fig.add_shape(type="circle", x0=3.5, y0=0, x1=8.5, y1=5,
              line=dict(color="#7ccfcd", width=4),
              fillcolor="rgba(124,207,205,0.18)")
fig.add_annotation(x=1.5, y=5.5, text="<b>8-gene panel</b><br><i>(canonical biology)</i>",
                   showarrow=False, font=dict(color="#F5A623", size=18))
fig.add_annotation(x=7.0, y=5.5, text="<b>Unsupervised top 20</b><br><i>(data-driven)</i>",
                   showarrow=False, font=dict(color="#7ccfcd", size=18))
fig.add_annotation(x=1.6, y=2.5, text="<b>panel only<br>5 genes</b>",
                   showarrow=False, font=dict(color="#F5A623", size=15))
fig.add_annotation(x=1.6, y=1.4,
                   text="SLC5A5 (NIS)<br>TG · TSHR<br>PAX8 · NKX2-1",
                   showarrow=False, font=dict(color="#EAEAEA", size=14, family="JetBrains Mono"))
fig.add_annotation(x=4.25, y=2.7, text="<b>overlap<br>3 genes</b>",
                   showarrow=False, font=dict(color="#fff", size=15))
fig.add_annotation(x=4.25, y=1.6, text="TPO<br>DIO1<br>FOXE1",
                   showarrow=False, font=dict(color="#F5A623", size=14, family="JetBrains Mono"))
fig.add_annotation(x=6.9, y=2.5, text="<b>unsup only<br>17 genes</b>",
                   showarrow=False, font=dict(color="#7ccfcd", size=15))
fig.add_annotation(x=6.9, y=1.0,
                   text="<b>MAPK feedback</b>: DUSP5/6/4 · ETV4/5 · FOSL1<br>" +
                        "<b>Inflammation</b>: HLA-DRA · MMP9 · FOXP3<br>" +
                        "<b>Differentiation alt.</b>: SLC5A8 · DIO2 · IYD · SLC26A4 · THRA<br>" +
                        "<b>Senescence/RTK</b>: CDKN2A · MET · KLK10",
                   showarrow=False, font=dict(color="#EAEAEA", size=12), align="left")
fig.update_layout(title=dict(text="Panel vs Unsupervised top — overlap 3 genes",
                             font=dict(size=16)),
                  xaxis=dict(visible=False, range=[-0.5, 9]),
                  yaxis=dict(visible=False, range=[-0.5, 6.5], scaleanchor="x"),
                  height=560, margin=dict(l=20, r=20, t=70, b=20), **DARK)
fig.write_html(FIG/"v17_8gene_venn.html", include_plotlyjs="cdn", full_html=True)
print("venn re-rendered")

# ============================================================
# B · Network (larger nodes, larger labels, cleaner edges)
# ============================================================
nodes = {
    "PAX8":   (1.0, 3.0, "TF",        "#7ccfcd"),
    "NKX2-1": (1.0, 2.0, "TF",        "#7ccfcd"),
    "FOXE1":  (1.0, 1.0, "TF",        "#7ccfcd"),
    "TSHR":   (3.2, 3.4, "Receptor",  "#F5A623"),
    "SLC5A5": (4.8, 2.3, "Uptake ★",  "#E74C3C"),
    "TPO":    (6.4, 3.0, "Synthesis", "#F5A623"),
    "TG":     (6.4, 1.6, "Synthesis", "#F5A623"),
    "DIO1":   (8.0, 2.3, "Activate",  "#9b59b6"),
}
edges = [
    ("TSHR","PAX8"), ("TSHR","NKX2-1"),
    ("PAX8","SLC5A5"), ("PAX8","TPO"), ("PAX8","TG"),
    ("NKX2-1","SLC5A5"), ("NKX2-1","TPO"), ("NKX2-1","TG"),
    ("FOXE1","TPO"), ("FOXE1","TG"),
    ("SLC5A5","TPO"),
    ("TPO","TG"),
    ("TG","DIO1"),
]
edge_x, edge_y = [], []
for a,b in edges:
    edge_x += [nodes[a][0], nodes[b][0], None]
    edge_y += [nodes[a][1], nodes[b][1], None]
fig = go.Figure()
fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode="lines",
                         line=dict(width=2, color="rgba(255,255,255,0.30)"),
                         hoverinfo="skip", showlegend=False))
node_x, node_y, node_text, node_color = [], [], [], []
for n,(x,y,role,col) in nodes.items():
    node_x.append(x); node_y.append(y); node_text.append(n); node_color.append(col)
fig.add_trace(go.Scatter(x=node_x, y=node_y, mode="markers+text",
                         text=node_text, textposition="middle center",
                         textfont=dict(color="#0a0d12", size=14, family="JetBrains Mono"),
                         marker=dict(size=78, color=node_color,
                                     line=dict(color="#fff", width=3)),
                         hoverinfo="text", showlegend=False))
fig.add_annotation(x=1.0, y=4.0, text="<b>TF trio</b><br>(thyroid lineage masters)",
                   showarrow=False, font=dict(color="#7ccfcd", size=14))
fig.add_annotation(x=3.2, y=4.4, text="<b>TSH signal entry</b>",
                   showarrow=False, font=dict(color="#F5A623", size=14))
fig.add_annotation(x=4.8, y=3.4, text="<b>RAI uptake ★</b><br>(임상 핵심)",
                   showarrow=False, font=dict(color="#E74C3C", size=14))
fig.add_annotation(x=6.4, y=4.1, text="<b>Hormone synthesis</b>",
                   showarrow=False, font=dict(color="#F5A623", size=14))
fig.add_annotation(x=8.0, y=3.4, text="<b>T4→T3 활성화</b>",
                   showarrow=False, font=dict(color="#9b59b6", size=14))
fig.update_layout(title=dict(text="8-gene functional network — TSH→TF→NIS→TPO/TG→DIO1",
                             font=dict(size=16)),
                  xaxis=dict(visible=False, range=[0, 9]),
                  yaxis=dict(visible=False, range=[0.2, 5], scaleanchor="x"),
                  height=500, margin=dict(l=20, r=20, t=60, b=20), **DARK)
fig.write_html(FIG/"v17_8gene_network.html", include_plotlyjs="cdn", full_html=True)
print("network re-rendered")

# ============================================================
# C · MAPK→NIS mechanism (cleaner, no overlapping text)
# ============================================================
fig = go.Figure()
boxes = [
    (0.2, 1.7, 2.6, 3.0, "Normal thyrocyte<br>NIS↑ TPO↑ TG↑<br><b>RAI uptake 정상</b>", "#2ECC71"),
    (3.6, 1.7, 6.0, 3.0, "BRAF V600E PTC<br>MAPK constitutive ON<br><b>CREB-PAX8 회로 차단</b>", "#F5A623"),
    (7.0, 1.7, 9.4, 3.0, "Dedifferentiated<br>NIS↓ TPO↓ TG↓<br><b>RAI uptake LOST</b>", "#E74C3C"),
]
for x0,y0,x1,y1,lab,col in boxes:
    fig.add_shape(type="rect", x0=x0, y0=y0, x1=x1, y1=y1,
                  line=dict(color=col, width=3), opacity=0.20)
    fig.add_annotation(x=(x0+x1)/2, y=(y0+y1)/2, text=lab,
                       showarrow=False, font=dict(color=col, size=14), align="center")
# big arrows between boxes
for ax, x in [(2.6, 3.6), (6.0, 7.0)]:
    fig.add_annotation(x=x, y=2.35, ax=ax, ay=2.35,
                       xref="x", yref="y", axref="x", ayref="y",
                       showarrow=True, arrowhead=3, arrowsize=1.6, arrowwidth=3,
                       arrowcolor="#aaa")
# Therapeutic intervention box
fig.add_shape(type="rect", x0=3.6, y0=0.3, x1=6.0, y1=1.3,
              line=dict(color="#7ccfcd", width=3, dash="dash"), opacity=0.20)
fig.add_annotation(x=4.8, y=0.8,
                   text="<b>MAPK inhibition</b><br>selumetinib · dabrafenib + trametinib",
                   showarrow=False, font=dict(color="#7ccfcd", size=13), align="center")
fig.add_annotation(x=4.8, y=1.55, text="↑ block", showarrow=False,
                   font=dict(color="#7ccfcd", size=13))
# Citations row (single line, larger)
fig.add_annotation(x=4.7, y=-0.35,
                   text="<b>References</b>: Knauf 2003 (BRAF→NIS↓) · Chakravarty 2011 JCI (BRS52) · "
                        "Ho 2013 NEJM (selumetinib RAI 회복) · Iravani 2019 (dab+tram)",
                   showarrow=False, font=dict(color="#888", size=12), align="center")
fig.update_layout(title=dict(text="BRAF V600E → MAPK → NIS suppression → RAI uptake loss",
                             font=dict(size=16)),
                  xaxis=dict(visible=False, range=[0, 9.6]),
                  yaxis=dict(visible=False, range=[-0.7, 3.5], scaleanchor="x"),
                  height=460, margin=dict(l=20,r=20,t=60,b=40), **DARK)
fig.write_html(FIG/"v17_mapk_nis_mechanism.html", include_plotlyjs="cdn", full_html=True)
print("mechanism re-rendered")

# ============================================================
# D · Heatmap (smaller n, sorted properly, bigger gene labels)
# ============================================================
try:
    expr = np.load("/opt/thyroid-dash/project/data_processed/v5_cross_cancer/THCA/X_combined.npz")["X"]
    with open("/opt/thyroid-dash/project/data_processed/v5_cross_cancer/THCA/shared_genes.txt") as fh:
        genes = [g.strip() for g in fh]
    idx = {g:i for i,g in enumerate(genes)}
    cols = [g for g in GENES if g in idx]
    if cols:
        mat = expr[:, [idx[g] for g in cols]]
        mu = mat.mean(0); sd = mat.std(0)+1e-9
        z = (mat - mu) / sd
        # sort by mean panel z-score (low=DM1-like, high=DM2-like)
        order = np.argsort(z.mean(axis=1))
        z = z[order]
        # take alternating subset for clearer display (every 2nd sample)
        z = z[::2][:120]  # ~120 samples spread evenly
        fig = go.Figure(go.Heatmap(z=z.T, y=cols, x=list(range(z.shape[0])),
                                   colorscale="RdBu_r", zmid=0, zmin=-2.5, zmax=2.5,
                                   colorbar=dict(title=dict(text="z-score", font=dict(size=14)),
                                                 tickfont=dict(size=12))))
        # Add side annotations
        fig.add_annotation(x=0, y=-1.5, xref="x", yref="y", showarrow=False,
                           text="<b>← DM1 (low panel expr, dedifferentiated)</b>",
                           font=dict(color="#E74C3C", size=13))
        fig.add_annotation(x=z.shape[0]-1, y=-1.5, xref="x", yref="y", showarrow=False,
                           text="<b>DM2 (high panel expr, differentiated) →</b>",
                           font=dict(color="#2ECC71", size=13), xanchor="right")
        fig.update_layout(title=dict(text=f"TCGA-THCA × 8-gene panel z-score (n={z.shape[0]} samples, sorted by mean panel expression)",
                                     font=dict(size=15)),
                          yaxis=dict(autorange="reversed", tickfont=dict(size=14, family="JetBrains Mono")),
                          xaxis=dict(showticklabels=False, title=dict(text="samples (sorted)", font=dict(size=13))),
                          height=460, margin=dict(l=120, r=20, t=80, b=60), **DARK)
        fig.write_html(FIG/"v17_8gene_heatmap.html", include_plotlyjs="cdn", full_html=True)
        print(f"heatmap re-rendered (n={z.shape[0]})")
except Exception as e:
    print(f"heatmap fail: {e}")

print("\nDONE")
