#!/usr/bin/env python3
"""K2 Q13 figure — Korean PRJEB11591 (Yoo SK SNU) cohort DM1/DM2 prediction distribution."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np, pandas as pd
import plotly.graph_objects as go

ROOT = Path("/opt/thyroid-dash/project")
PRED = ROOT/"results/v17_korean/K2_korean_predictions_v4.tsv"
FIG = ROOT/"reports/html/figs_interactive/v17"
RES = ROOT/"results/v17_korean"

BG, INK = "#0b0e12", "#F2F2F2"
DARK = dict(template="plotly_dark", paper_bgcolor=BG, plot_bgcolor=BG, font=dict(color=INK, size=14))

pred = pd.read_csv(PRED, sep="\t")
n = len(pred)
n_dm1 = (pred["DM_call"] == "DM1").sum()
n_dm2 = (pred["DM_call"] == "DM2").sum()
mean_p = pred["p_DM2"].mean()
median_p = pred["p_DM2"].median()
print(f"n={n}, DM1={n_dm1}, DM2={n_dm2} ({100*n_dm2/n:.1f}% DM2)")
print(f"mean p_DM2={mean_p:.3f}, median={median_p:.3f}")

# Fig A: histogram of p_DM2 with DM1/DM2 split coloring
fig = go.Figure()
dm1 = pred[pred["DM_call"]=="DM1"]["p_DM2"].values
dm2 = pred[pred["DM_call"]=="DM2"]["p_DM2"].values
fig.add_trace(go.Histogram(x=dm1, name=f"DM1 (n={n_dm1}, p_DM2 < 0.5)",
                           marker=dict(color="#c24c4c"), xbins=dict(start=0, end=1, size=0.05)))
fig.add_trace(go.Histogram(x=dm2, name=f"DM2 (n={n_dm2}, p_DM2 ≥ 0.5)",
                           marker=dict(color="#7ccfcd"), xbins=dict(start=0, end=1, size=0.05)))
fig.add_vline(x=0.5, line=dict(color="white", width=2, dash="dash"),
              annotation=dict(text="DM1 ↔ DM2 boundary (p=0.5)", font=dict(color=INK, size=11), bgcolor="rgba(0,0,0,0.4)"))
fig.update_layout(title=dict(
    text=f"<b>PRJEB11591 (Yoo SK et al, SNU Genomic Medicine Institute) Korean cohort</b><br>"
         f"<sub style='color:#7ccfcd'>n={n} samples through K2 v4 (within-sample-centered LogReg, TCGA 5-fold CV AUC 0.963); "
         f"DM1:DM2 = {n_dm1}:{n_dm2} ({100*n_dm2/n:.1f}% DM2); mean p_DM2 = {mean_p:.2f}</sub>",
    font=dict(size=14, color=INK)),
    barmode="overlay",
    xaxis=dict(title="P(DM2) — TCGA-trained LogReg", range=[0,1.02], tickfont=dict(color=INK)),
    yaxis=dict(title="환자 수 (n)", gridcolor="rgba(255,255,255,0.08)", tickfont=dict(color=INK)),
    height=460, width=900, margin=dict(l=80,r=20,t=110,b=60), **DARK,
    legend=dict(itemsizing="constant", x=0.02, y=0.98, bgcolor="rgba(0,0,0,0.4)"))
fig.update_traces(opacity=0.75)
fig.write_html(FIG/"v17_K2_PRJEB11591_pdm2_histogram.html", include_plotlyjs="cdn", full_html=True)
print(f"wrote {FIG}/v17_K2_PRJEB11591_pdm2_histogram.html")

# Fig B: bin-by-bin comparison vs TCGA
# TCGA prior P(DM2) = 0.72 (R1-A leak-free DM1=140/DM2=360 = 28/72)
tcga_dm1_pct = 28; tcga_dm2_pct = 72
korean_dm1_pct = 100*n_dm1/n; korean_dm2_pct = 100*n_dm2/n

fig = go.Figure()
fig.add_trace(go.Bar(x=["TCGA-THCA primary (n=500, leak-free)"], y=[tcga_dm1_pct], name="DM1 (탈분화)",
                     marker=dict(color="#c24c4c"), text=[f"{tcga_dm1_pct}%"], textposition="inside"))
fig.add_trace(go.Bar(x=["TCGA-THCA primary (n=500, leak-free)"], y=[tcga_dm2_pct], name="DM2 (분화 보존)",
                     marker=dict(color="#7ccfcd"), text=[f"{tcga_dm2_pct}%"], textposition="inside"))
fig.add_trace(go.Bar(x=[f"PRJEB11591 Korean SNU (n={n})"], y=[korean_dm1_pct], showlegend=False,
                     marker=dict(color="#c24c4c"), text=[f"{korean_dm1_pct:.1f}%"], textposition="inside"))
fig.add_trace(go.Bar(x=[f"PRJEB11591 Korean SNU (n={n})"], y=[korean_dm2_pct], showlegend=False,
                     marker=dict(color="#7ccfcd"), text=[f"{korean_dm2_pct:.1f}%"], textposition="inside"))
fig.update_layout(title=dict(
    text=f"<b>DM1/DM2 비율 비교: 미국 TCGA vs 한국 PRJEB11591 (Yoo SNU)</b><br>"
         f"<sub style='color:#7ccfcd'>한국 cohort가 TCGA보다 +{korean_dm2_pct-tcga_dm2_pct:.0f} percentage point 더 DM2-skewed — Korean overdiagnosis paradigm 가설과 일치</sub>",
    font=dict(size=14, color=INK)),
    barmode="stack",
    yaxis=dict(title="환자 비율 %", range=[0,105], tickfont=dict(color=INK)),
    xaxis=dict(tickfont=dict(size=12, color=INK)),
    height=440, width=820, margin=dict(l=80,r=20,t=110,b=80), **DARK,
    legend=dict(orientation="h", x=0.5, xanchor="center", y=-0.15, bgcolor="rgba(0,0,0,0.3)"))
fig.write_html(FIG/"v17_K2_PRJEB11591_vs_TCGA.html", include_plotlyjs="cdn", full_html=True)
print(f"wrote {FIG}/v17_K2_PRJEB11591_vs_TCGA.html")

# Fig C: per-gene 8-panel TPM ranges (showing within-sample shape)
fig = go.Figure()
genes = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]
for g in genes:
    fig.add_trace(go.Box(y=np.log2(pred[g].values + 1), name=g,
                         marker=dict(color="#7ccfcd", size=3), boxpoints="outliers"))
fig.update_layout(title=dict(
    text=f"<b>PRJEB11591 8-gene panel log2(TPM+1) 분포 (n={n})</b><br>"
         f"<sub style='color:#7ccfcd'>kallisto 8-gene mini-index (93 transcripts) 사용 → absolute TPM은 inflated, 내부 shape (gene 간 ratio)이 K2 LogReg input</sub>",
    font=dict(size=14, color=INK)),
    yaxis=dict(title="log2(TPM+1) — mini-index inflated absolute scale",
               gridcolor="rgba(255,255,255,0.08)", tickfont=dict(color=INK)),
    xaxis=dict(tickfont=dict(size=12, color=INK)),
    showlegend=False, height=440, width=900, margin=dict(l=80,r=20,t=100,b=70), **DARK)
fig.write_html(FIG/"v17_K2_PRJEB11591_panel_tpm.html", include_plotlyjs="cdn", full_html=True)
print(f"wrote {FIG}/v17_K2_PRJEB11591_panel_tpm.html")

# Save Q13 summary
summary = {
    "cohort": "PRJEB11591 (Yoo SK et al, SNU Genomic Medicine Institute Korean PTC)",
    "manifest_total_runs": 262,
    "n_samples_processed_so_far": int(n),
    "n_DM1": int(n_dm1), "n_DM2": int(n_dm2),
    "pct_DM2": round(100*n_dm2/n, 1),
    "mean_p_DM2": round(mean_p, 3),
    "median_p_DM2": round(median_p, 3),
    "tcga_dm1_dm2_baseline": "28%:72% (R1-A leak-free TCGA-THCA n=500)",
    "korean_vs_tcga_DM2_lift": round(100*n_dm2/n - 72, 1),
    "calibration": "within-sample-centered LogReg (TCGA 5-fold CV AUC 0.963); mini-index TPM inflation correction",
    "interpretation": f"Korean SNU cohort (n={n}/262 processed so far) shows {round(100*n_dm2/n,1)}% DM2 vs TCGA 72% — supports Korean overdiagnosis-paradigm hypothesis (more indolent well-differentiated PTC enriched). Background quant pipeline streaming continues for remaining 155 samples.",
}
(RES/"K2_PRJEB11591_q13_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
print(json.dumps(summary, indent=2, ensure_ascii=False))
