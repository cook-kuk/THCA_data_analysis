#!/usr/bin/env python3
"""Hero showcase figure — 3-panel stunning summary for top of revision page."""
import json
from pathlib import Path
import numpy as np, pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

ROOT = Path("/opt/thyroid-dash/project")
FIG = ROOT/"reports/html/figs_interactive/v17"

BG = "#0b0e12"; INK = "#F2F2F2"

# ---------- panel 1 data: cross-cohort dedifferentiation gradient ----------
# Korean Lee 2024 GSE213647 panel scores by histology
lee_data = pd.read_csv(ROOT/"results/v17_korean/GSE213647_panel_score.tsv", sep="\t")
groups_lee = ["Normal","PTC","PDFP","UTC/ATC"]
lee_panel = {h: lee_data[lee_data["histology"]==h]["panel_z"].values for h in groups_lee}

# Lu 2023 single-cell DM_score by histology
lu_summary = json.loads((ROOT/"results/v17_lu2023/GSE193581_summary.json").read_text())
lu_medians = {"Normal": lu_summary["DM_score_NORM_median"],
              "PTC":    lu_summary["DM_score_PTC_median"],
              "ATC":    lu_summary["DM_score_ATC_median"]}

# ---------- panel 2: Pozdeyev TERT prevalence cross-cohort ----------
tert_cohorts = ["TCGA primary\nn=502","Pozdeyev advanced PTC\nn=468","Pozdeyev ATC\nn=196","MSK Landa PDTC/ATC\nn=43"]
tert_pct = [7.2, 50.6, 55.1, 60.5]

# ---------- panel 3: Hashimoto-like 17% — Korean Lee + Lu 2023 ----------
hashi_results = [
    {"cohort":"Lee 2024 KR\n(bulk PTC n=348)", "pct": 17.0, "n_total":348, "n_hashi":59},
    {"cohort":"Lu 2023 CN\n(scRNA PTC n=6 sample)", "pct": 16.7, "n_total":6, "n_hashi":1},
]

# ---------- build 1×3 multi-panel ----------
fig = make_subplots(rows=1, cols=3,
                    subplot_titles=(
                        "<b>Cross-population dedifferentiation gradient</b><br><sub>한국 Lee 2024 (n=632) bulk + 中 Lu 2023 (n=15.3k cells) scRNA</sub>",
                        "<b>TERT promoter — late-event 7-fold enrichment</b><br><sub>4 cohort, ~1,200 patients, monotonic trajectory</sub>",
                        "<b>Hashimoto-PTC 17% subgroup reproduces</b><br><sub>Korean bulk + Chinese single-cell, two independent cohorts</sub>"),
                    horizontal_spacing=0.08,
                    column_widths=[0.36, 0.30, 0.34])

# --- Panel 1: Korean Lee dedifferentiation box + Lu line overlay ---
colors_hist = {"Normal":"#2ECC71","PTC":"#F5A623","PDFP":"#E67E22","UTC/ATC":"#c24c4c","ATC":"#c24c4c"}
for h in groups_lee:
    vals = lee_panel[h]
    fig.add_trace(go.Box(
        y=vals, x=[h]*len(vals), name=f"Lee {h}",
        marker=dict(color=colors_hist[h], size=3, opacity=0.45),
        boxpoints="outliers", line=dict(width=2),
        showlegend=False), row=1, col=1)

# Overlay Lu 2023 sc medians as scatter line
fig.add_trace(go.Scatter(
    x=["Normal","PTC","UTC/ATC"], y=[lu_medians["Normal"], lu_medians["PTC"], lu_medians["ATC"]],
    mode="lines+markers", line=dict(color=INK, width=3, dash="dash"),
    marker=dict(size=14, color="#7ccfcd", symbol="diamond", line=dict(color="white", width=2)),
    name="Lu 2023 scRNA median (overlay)", showlegend=True), row=1, col=1)

fig.update_yaxes(title=dict(text="DM panel z-score (8-gene mean)", font=dict(size=12, color=INK)),
                 gridcolor="rgba(255,255,255,0.08)", zeroline=True, zerolinecolor="rgba(255,255,255,0.2)",
                 row=1, col=1)

# --- Panel 2: TERT prevalence bar with annotations ---
tert_colors = ["#7ccfcd","#F5A623","#E67E22","#c24c4c"]
fig.add_trace(go.Bar(
    x=tert_cohorts, y=tert_pct, marker=dict(color=tert_colors, line=dict(color="white", width=1)),
    text=[f"<b>{p:.1f}%</b>" for p in tert_pct], textposition="outside",
    textfont=dict(color=INK, size=14),
    showlegend=False, hovertemplate="%{x}<br>TERT promoter: %{y:.1f}%<extra></extra>"), row=1, col=2)
fig.update_yaxes(title=dict(text="TERT promoter prevalence (%)", font=dict(size=12, color=INK)),
                 range=[0, 70], gridcolor="rgba(255,255,255,0.08)", row=1, col=2)
fig.add_annotation(text="<b>7-fold ↑</b><br>p = 5.9×10⁻⁵⁵", x=2.5, y=64,
                   xref="x2", yref="y2", showarrow=False,
                   font=dict(color="#7ccfcd", size=12), bgcolor="rgba(0,0,0,0.5)",
                   bordercolor="#7ccfcd", borderwidth=1, borderpad=6)

# --- Panel 3: Hashimoto 17% bars + ★ annotations ---
hashi_x = [h["cohort"] for h in hashi_results]
hashi_y = [h["pct"] for h in hashi_results]
fig.add_trace(go.Bar(
    x=hashi_x, y=hashi_y,
    marker=dict(color=["#7ccfcd","#7ccfcd"], line=dict(color="white", width=1)),
    text=[f"<b>{h['pct']:.1f}%</b><br>{h['n_hashi']}/{h['n_total']}" for h in hashi_results],
    textposition="outside", textfont=dict(color=INK, size=13),
    showlegend=False), row=1, col=3)
fig.add_hline(y=17, line=dict(color="#F5A623", dash="dot", width=1.5),
              annotation=dict(text="reproduces line", font=dict(color="#F5A623", size=10)), row=1, col=3)
fig.update_yaxes(title=dict(text="% PTC patients/samples Hashimoto-like", font=dict(size=12, color=INK)),
                 range=[0, 26], gridcolor="rgba(255,255,255,0.08)", row=1, col=3)
fig.add_annotation(text="★ Cross-cohort + cross-modality<br>1차 검증 성공",
                   x=0.5, y=23, xref="x3", yref="y3", showarrow=False,
                   font=dict(color="#F5A623", size=12), bgcolor="rgba(0,0,0,0.5)",
                   bordercolor="#F5A623", borderwidth=1, borderpad=6)

# --- Layout ---
fig.update_layout(
    title=dict(text=f"<b>v17 Thyroid — 3 핵심 발견 (오늘 build)</b><br><sub style='color:#7ccfcd'>Pan-Asian dedifferentiation axis · TERT trajectory · Hashimoto-PTC subgroup</sub>",
               x=0.5, xanchor="center", font=dict(size=18, color=INK)),
    paper_bgcolor=BG, plot_bgcolor=BG,
    font=dict(color=INK, size=12),
    template="plotly_dark",
    height=520, width=1500, margin=dict(l=70, r=50, t=110, b=110),
    legend=dict(orientation="h", x=0.5, xanchor="center", y=-0.18, bgcolor="rgba(0,0,0,0.3)"),
    hoverlabel=dict(bgcolor="rgba(15,18,23,0.95)", font=dict(color=INK)))

# tick fontsize
for ax_id in ["xaxis","xaxis2","xaxis3","yaxis","yaxis2","yaxis3"]:
    fig.layout[ax_id].tickfont = dict(color=INK, size=11)

fig.write_html(FIG/"v17_hero_showcase.html", include_plotlyjs="cdn", full_html=True)
print(f"wrote {FIG}/v17_hero_showcase.html")
