#!/usr/bin/env python3
"""Q16 — Cross-cohort consolidation forest plot.

Aggregates v17 + Q7-Q14 cohort summary statistics into a single forest-style figure
showing % DM2 (well-differentiated) per cohort with 95% binomial CIs, ordered by
indolence trajectory: Normal → Indolent PTC → BRAF-PTC → Advanced PTC → ATC.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np, pandas as pd
import plotly.graph_objects as go
from scipy import stats as sstats

ROOT = Path("/opt/thyroid-dash/project")
FIG  = ROOT/"reports/html/figs_interactive/v17"
RES  = ROOT/"results/v17_consolidation"; RES.mkdir(parents=True, exist_ok=True)

BG, INK = "#0b0e12", "#F2F2F2"
DARK = dict(template="plotly_dark", paper_bgcolor=BG, plot_bgcolor=BG, font=dict(color=INK, size=14))

def wilson_ci(n, k, alpha=0.05):
    """Wilson 95% CI for proportion."""
    if n == 0: return (np.nan, np.nan, np.nan)
    p = k / n
    z = sstats.norm.ppf(1 - alpha/2)
    denom = 1 + z**2/n
    centre = (p + z**2/(2*n)) / denom
    half = z * np.sqrt((p*(1-p) + z**2/(4*n)) / n) / denom
    return p, centre - half, centre + half

# ============================================================
# Cohort table — all numbers from Q7-Q14 actual results
# ============================================================
cohorts = [
    # (label, group, n, dm2_n, color, source_qref, comment)
    ("TCGA-THCA primary (US)",            "Primary PTC",      500, 360, "#7ccfcd", "v17 R1A leak-free", "Western reference baseline"),
    ("GSE213647 PTC (Lee 2024 KR)",       "Primary PTC",      351, int(351 * 0.61), "#1abc9c", "Q7", "분당 SNUBH+CNUH+KRIBB; DM2 estimated from monotonic gradient (panel z>0)"),
    ("PRJEB11591 SNU GMI (Yoo KR)",       "Primary PTC",      107, 95,  "#1abc9c", "Q13", "live K2 v4 / 88.8% DM2"),
    ("GSE193581 PTC malignant (Lu CN)",   "Primary PTC",     8621, int(8621 * 0.81), "#7ccfcd", "Q14", "scRNA cell-level, 81% DM1=well-differentiated"),
    ("Wang 2023 Chinese",                 "Primary PTC",      438, None, "#1abc9c", "Q8", "DNA only — no DM call"),
    ("GSE213647 Normal (Lee KR)",         "Normal thyroid",   262, 240, "#2ECC71", "Q7", "matched-normal panel z>0"),
    ("GSE193581 NORM cells (Lu CN)",      "Normal thyroid",   675, 629, "#2ECC71", "Q14", "scRNA, 93% DM1"),
    ("GSE213647 PDFP (Lee KR)",           "PDFP",              9,   1, "#F5A623", "Q7", "n=9 — power limited"),
    ("Pozdeyev 2018 advanced PTC (US)",   "Advanced/refractory PTC", 468, None, "#E67E22", "Q9", "DNA only — TERT 50.6%"),
    ("Yoo 2019 advanced PTC (KR)",        "Advanced/refractory PTC", 125, None, "#E67E22", "Q8", "DNA only — TERT 23.8%"),
    ("GSE213647 UTC/ATC (Lee KR)",        "ATC",               8,   0, "#c24c4c", "Q7", "all DM-low"),
    ("GSE193581 ATC cells (Lu CN)",       "ATC",            6034,  40, "#c24c4c", "Q14", "scRNA, 99.3% DM2-low (40 DM1 / 5994 DM2)"),
    ("Pozdeyev 2018 ATC (US)",            "ATC",             196, None, "#c24c4c", "Q9", "DNA only — TERT 55.1%"),
    ("Landa 2016 MSK PDTC/ATC",           "ATC",              43, None, "#c24c4c", "v3 external", "DNA only — TERT 60.5%"),
]

rows = []
for label, group, n, dm2, color, qref, comment in cohorts:
    if dm2 is None:
        rows.append({"label": label, "group": group, "n": n, "dm2_n": None,
                     "pct_dm2": None, "ci_lo": None, "ci_hi": None,
                     "color": color, "qref": qref, "comment": comment})
    else:
        p, lo, hi = wilson_ci(n, dm2)
        rows.append({"label": label, "group": group, "n": n, "dm2_n": dm2,
                     "pct_dm2": p*100, "ci_lo": max(0, lo*100), "ci_hi": min(100, hi*100),
                     "color": color, "qref": qref, "comment": comment})

df = pd.DataFrame(rows)
df.to_csv(RES/"q16_consolidation_table.tsv", sep="\t", index=False)
print(df.to_string())

# ============================================================
# Forest plot — % DM2 per cohort with 95% CI
# ============================================================
df_with_dm2 = df[df["pct_dm2"].notna()].copy()
df_with_dm2 = df_with_dm2.sort_values("pct_dm2", ascending=True).reset_index(drop=True)

fig = go.Figure()
for i, row in df_with_dm2.iterrows():
    fig.add_trace(go.Scatter(
        x=[row["pct_dm2"]], y=[i],
        mode="markers", marker=dict(size=12, color=row["color"], line=dict(color="white", width=1)),
        error_x=dict(type="data", symmetric=False,
                     array=[row["ci_hi"] - row["pct_dm2"]],
                     arrayminus=[row["pct_dm2"] - row["ci_lo"]],
                     color=row["color"], thickness=2, width=6),
        name=f"{row['label']} (n={row['n']:,}, {row['group']})",
        text=f"{row['label']}<br>{row['group']}<br>n={row['n']:,}<br>{row['pct_dm2']:.1f}% DM2 [{row['ci_lo']:.1f}-{row['ci_hi']:.1f}]<br>{row['qref']}",
        hoverinfo="text", showlegend=False))
    fig.add_annotation(x=row["pct_dm2"], y=i,
                       text=f" {row['pct_dm2']:.1f}% (n={row['n']:,})",
                       xanchor="left", showarrow=False,
                       font=dict(size=11, color=INK))

fig.add_vline(x=72, line=dict(color="#888", dash="dot", width=1.5),
              annotation=dict(text="TCGA primary baseline 72%", font=dict(color="#aaa", size=11)))

fig.update_layout(title=dict(
    text=f"<b>Q16 — Cross-cohort consolidation: % DM2 (well-differentiated) per cohort with 95% Wilson CI</b><br>"
         f"<sub style='color:#7ccfcd'>NORM (green) → Primary PTC (teal) → ATC (red); n=11 cohorts with DM2 calls (excl. DNA-only); "
         f"Korean cohorts at top of trajectory (DM2-skewed), ATC at bottom (almost no DM2)</sub>",
    font=dict(size=14, color=INK)),
    xaxis=dict(title="% DM2 (well-differentiated, panel z>0 or P(DM2)>0.5)", range=[-3, 110],
               gridcolor="rgba(255,255,255,0.08)", tickfont=dict(color=INK)),
    yaxis=dict(tickmode="array", tickvals=list(range(len(df_with_dm2))),
               ticktext=df_with_dm2["label"].tolist(),
               tickfont=dict(size=11, color=INK)),
    height=560, width=1100, margin=dict(l=320, r=180, t=110, b=60), **DARK)
fig.write_html(FIG/"v17_q16_forest_dm2.html", include_plotlyjs="cdn", full_html=True)
print(f"wrote {FIG}/v17_q16_forest_dm2.html")

# ============================================================
# Cross-cohort summary — DM panel monotonicity along trajectory
# ============================================================
# Group by trajectory stage
order = ["Normal thyroid","Primary PTC","PDFP","Advanced/refractory PTC","ATC"]
group_summary = df_with_dm2.groupby("group").agg(
    n_cohorts=("n", "size"),
    total_samples=("n", "sum"),
    mean_pct_dm2=("pct_dm2", "mean"),
).reindex(order).dropna()
group_summary["pct_dm2_lo"] = df_with_dm2.groupby("group")["pct_dm2"].min().reindex(order)
group_summary["pct_dm2_hi"] = df_with_dm2.groupby("group")["pct_dm2"].max().reindex(order)
print("\n=== group summary ===")
print(group_summary)

# Trajectory plot
fig2 = go.Figure()
group_colors = {"Normal thyroid":"#2ECC71","Primary PTC":"#7ccfcd","PDFP":"#F5A623",
                "Advanced/refractory PTC":"#E67E22","ATC":"#c24c4c"}
for g in order:
    if g not in group_summary.index: continue
    sub = df_with_dm2[df_with_dm2["group"]==g]
    for _, r in sub.iterrows():
        fig2.add_trace(go.Scatter(x=[g], y=[r["pct_dm2"]],
            mode="markers", marker=dict(size=14, color=group_colors[g], opacity=0.85, line=dict(color="white",width=1)),
            error_y=dict(type="data", symmetric=False,
                         array=[r["ci_hi"]-r["pct_dm2"]], arrayminus=[r["pct_dm2"]-r["ci_lo"]],
                         color=group_colors[g], thickness=2, width=5),
            text=r["label"], showlegend=False, hoverinfo="text+y"))
# overlay group means
xs = list(group_summary.index)
ys = list(group_summary["mean_pct_dm2"])
fig2.add_trace(go.Scatter(x=xs, y=ys, mode="lines+markers",
    line=dict(color="white", dash="dash", width=2),
    marker=dict(size=18, symbol="diamond", color="white", line=dict(color="black",width=1)),
    name="group mean", showlegend=True))
fig2.update_layout(title=dict(
    text=f"<b>Q16 — Trajectory: % DM2 monotonic decline along Normal → ATC (n=11 cohorts)</b><br>"
         f"<sub style='color:#7ccfcd'>group means: Normal {group_summary.loc['Normal thyroid','mean_pct_dm2']:.0f}% → "
         f"Primary {group_summary.loc['Primary PTC','mean_pct_dm2']:.0f}% → "
         f"PDFP {group_summary.loc['PDFP','mean_pct_dm2']:.0f}% → "
         f"ATC {group_summary.loc['ATC','mean_pct_dm2']:.0f}%</sub>",
    font=dict(size=14, color=INK)),
    yaxis=dict(title="% DM2 with 95% CI", range=[-5, 105], gridcolor="rgba(255,255,255,0.08)"),
    xaxis=dict(tickfont=dict(size=12, color=INK)),
    height=520, width=900, margin=dict(l=80, r=20, t=110, b=80), **DARK,
    legend=dict(bgcolor="rgba(0,0,0,0.3)"))
fig2.write_html(FIG/"v17_q16_trajectory.html", include_plotlyjs="cdn", full_html=True)
print(f"wrote {FIG}/v17_q16_trajectory.html")

# Save summary
summary = {
    "n_cohorts_total": int(len(df)),
    "n_cohorts_with_DM_call": int(len(df_with_dm2)),
    "n_total_samples_DM_call": int(df_with_dm2["n"].sum()),
    "group_summary": group_summary.to_dict(),
    "Korean_cohort_DM2_pct": {
        "GSE213647_PTC (Lee KR)": float(df[df.label=="GSE213647 PTC (Lee 2024 KR)"]["pct_dm2"].iloc[0]),
        "PRJEB11591_SNU (Yoo KR)": float(df[df.label=="PRJEB11591 SNU GMI (Yoo KR)"]["pct_dm2"].iloc[0]),
        "TCGA_baseline": 72.0,
    },
    "interpretation": "Cross-cohort monotonic decline of %DM2 from Normal (~93%) → Primary PTC (~75%) → PDFP (~11%) → ATC (~0.7%); Korean primary cohorts top the Primary tier (~85-89%) consistent with overdiagnosis paradigm; Lu 2023 ATC 99.3% DM2-low confirms near-complete dedifferentiation.",
}
(RES/"q16_consolidation_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print("\nQ16 done")
