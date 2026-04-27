#!/usr/bin/env python3
"""v17 ACTUAL ship — Figure 8 panels + interactive plots for dashboard.

Generates:
- figures/figure8_forest_HR.pdf/png   (4 Cox HRs with 95% CI)
- figures/figure8_bootstrap_p.pdf/png (1000-iter p distribution)
- figures/figure8_LOO_p.pdf/png       (36 LOO p values)
- reports/html/figs_interactive/v17/{forest,bootstrap,loo}.html (Plotly)
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.io as pio

ROOT = Path("/opt/thyroid-dash/project")
RES  = ROOT / "results" / "v17_actual"
OUT_FIG = ROOT / "submission" / "npj" / "figures"
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_INT = ROOT / "reports" / "html" / "figs_interactive" / "v17"
OUT_INT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.size": 10, "font.family": "serif", "savefig.dpi": 300,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.06,
})

# ---------- A1-A: forest plot -----------------------------------------
cox = pd.read_csv(RES / "A1A_cox_unadjusted.tsv", sep="\t")
cox = cox.sort_values("label", key=lambda s: s.map({
    "TERT only": 0, "TERT + age": 1, "TERT + stage": 2, "TERT + stage + age + sex": 3
}))

fig, ax = plt.subplots(figsize=(7, 3.6))
ypos = np.arange(len(cox))
xerr_lo = cox["HR"].values - cox["HR_lo"].values
xerr_hi = cox["HR_hi"].values - cox["HR"].values
ax.errorbar(cox["HR"].values, ypos, xerr=[xerr_lo, xerr_hi], fmt="o",
            color="#1a4080", capsize=4, lw=1.4, mec="black", mew=0.6, ms=8)
ax.axvline(1, color="red", ls="--", lw=0.8, alpha=0.6, label="HR = 1 (no effect)")
ax.set_yticks(ypos)
ax.set_yticklabels(cox["label"], fontsize=10)
ax.set_xscale("log")
ax.set_xlabel("Hazard Ratio (95 % CI), TERT⁺ vs WT")
ax.set_xlim(0.4, 22)
ax.grid(alpha=0.3, axis="x", which="both")
# annotate p
for i, row in enumerate(cox.itertuples()):
    ax.text(row.HR_hi * 1.05, i,
            f"HR={row.HR:.2f} (CI {row.HR_lo:.2f}–{row.HR_hi:.2f}, p={row.p:.3g})",
            fontsize=8.5, va="center")
ax.set_title("Figure 8B — Forest plot of TERT⁺ Hazard Ratio (lifelines penalizer = 0.01)\n"
             "Univariate signal collapses to non-significant after stage + age adjustment",
             fontsize=10)
ax.legend(loc="upper left", frameon=False)
fig.savefig(OUT_FIG / "figure8_forest_HR.pdf")
fig.savefig(OUT_FIG / "figure8_forest_HR.png")
plt.close(fig)
print(f"wrote figure8_forest_HR.{{pdf,png}}")

# ---------- A1-B: bootstrap p distribution ----------------------------
boot = pd.read_csv(RES / "A1B_bootstrap_logrank.tsv", sep="\t")
fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
axes[0].hist(np.log10(boot["p"].clip(lower=1e-20).values), bins=40,
             color="#1a4080", alpha=0.75, edgecolor="white")
for q, c, lbl in [(0.025,"red","2.5%"), (0.5,"green","median"), (0.975,"red","97.5%")]:
    v = np.log10(boot["p"].clip(lower=1e-20).quantile(q))
    axes[0].axvline(v, color=c, ls="--", lw=1.2, label=f"{lbl}: 10^{{{v:.1f}}}")
axes[0].axvline(np.log10(0.05), color="orange", ls=":", lw=1.2, label="α = 0.05")
axes[0].set_xlabel("log₁₀(logrank p)")
axes[0].set_ylabel("Count (1 000 bootstrap iters)")
axes[0].set_title("Figure 8E.i — Bootstrap p distribution\n(median p = 2.6×10⁻⁵; 93 % iter p < 0.05)")
axes[0].legend(fontsize=8, frameon=False)
axes[0].grid(alpha=0.3)

axes[1].hist(boot["rate_ratio"].clip(upper=30).dropna().values, bins=30,
             color="#0d6765", alpha=0.75, edgecolor="white")
for q, c, lbl in [(0.025,"red","2.5%"), (0.5,"green","median"), (0.975,"red","97.5%")]:
    v = boot["rate_ratio"].quantile(q)
    axes[1].axvline(v, color=c, ls="--", lw=1.2, label=f"{lbl}: {v:.2f}×")
axes[1].axvline(1, color="orange", ls=":", lw=1.2, label="rate ratio = 1")
axes[1].set_xlabel("TERT⁺ / triple-negative event-rate ratio")
axes[1].set_ylabel("Count")
axes[1].set_title("Figure 8E.ii — Effect-size CI\n(median 4.7×, 95 % CI 1.2×–17.3×)")
axes[1].legend(fontsize=8, frameon=False)
axes[1].grid(alpha=0.3)
fig.savefig(OUT_FIG / "figure8_bootstrap_p.pdf")
fig.savefig(OUT_FIG / "figure8_bootstrap_p.png")
plt.close(fig)
print(f"wrote figure8_bootstrap_p.{{pdf,png}}")

# ---------- A1-C: LOO p values ----------------------------------------
loo = pd.read_csv(RES / "A1C_loo_logrank.tsv", sep="\t")
loo = loo.sort_values("logrank_p").reset_index(drop=True)
A1C = json.loads((RES / "A1C_summary.json").read_text())
fig, ax = plt.subplots(figsize=(7.5, 3.8))
colors = ["#c14" if e==1 else "#1a4080" for e in loo["is_event"]]
ax.bar(np.arange(len(loo)), -np.log10(loo["logrank_p"].clip(lower=1e-20)),
       color=colors, edgecolor="white", linewidth=0.4)
ax.axhline(-np.log10(0.05), color="orange", ls=":", label="α = 0.05")
ax.axhline(-np.log10(0.001), color="red", ls="--", lw=0.9, label="α = 0.001")
ax.axhline(-np.log10(A1C["p_full"]), color="green", ls="-", lw=0.7, alpha=0.6,
           label=f"full-cohort p = {A1C['p_full']:.1g}")
# Legend for bar colours
from matplotlib.patches import Patch
extra = [Patch(facecolor="#c14", label="event patient (n=6)"),
         Patch(facecolor="#1a4080", label="non-event patient (n=30)")]
ax.legend(handles=ax.get_legend_handles_labels()[0] + extra,
          loc="lower right", fontsize=8, frameon=False)
ax.set_xlabel("36 leave-one-out iterations (sorted by p)")
ax.set_ylabel("-log₁₀(logrank p)")
ax.set_title(f"Figure 8E.iii — Leave-one-out sensitivity ({len(loo)} TERT⁺ patients)\n"
             f"Worst-case p = {A1C['p_max_loo']:.4g} (all 36 LOO p < 0.001)")
ax.grid(alpha=0.3, axis="y")
fig.savefig(OUT_FIG / "figure8_LOO_p.pdf")
fig.savefig(OUT_FIG / "figure8_LOO_p.png")
plt.close(fig)
print(f"wrote figure8_LOO_p.{{pdf,png}}")

# ---------- Plotly interactive versions for dashboard -----------------
DARK = dict(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
            font=dict(color="#EAEAEA", family="Inter, sans-serif"),
            margin=dict(l=60, r=40, t=70, b=60))
AXIS = dict(gridcolor="rgba(234,234,234,0.08)", zerolinecolor="rgba(234,234,234,0.15)")

# Forest plot
fig = go.Figure()
fig.add_trace(go.Scatter(
    x=cox["HR"], y=cox["label"], mode="markers",
    error_x=dict(type="data", symmetric=False,
                 array=cox["HR_hi"] - cox["HR"], arrayminus=cox["HR"] - cox["HR_lo"]),
    marker=dict(color="#F5A623", size=12,
                line=dict(color="rgba(255,255,255,0.5)", width=1)),
    hovertemplate="%{y}<br>HR=%{x:.2f}<extra></extra>",
    name="HR (95% CI)",
))
fig.add_vline(x=1, line=dict(color="red", dash="dash", width=1),
              annotation_text="HR = 1", annotation_font_color="rgba(234,234,234,0.55)")
fig.update_layout(
    title=dict(text="TERT⁺ Cox HR — univariate signal collapses after stage + age adjustment",
               font=dict(size=15, color="#EAEAEA")),
    xaxis=dict(title="Hazard Ratio (log scale)", type="log", range=[np.log10(0.4), np.log10(22)], **AXIS),
    yaxis=dict(title="", **AXIS),
    height=360, **DARK,
)
pio.write_html(fig, str(OUT_INT / "v17_forest_HR.html"), include_plotlyjs="cdn", full_html=True)
print(f"wrote v17_forest_HR.html")

# Bootstrap p
fig = go.Figure()
fig.add_trace(go.Histogram(
    x=np.log10(boot["p"].clip(lower=1e-20).values), nbinsx=40,
    marker=dict(color="#F5A623", line=dict(color="rgba(255,255,255,0.3)", width=0.5)),
))
fig.add_vline(x=np.log10(0.05), line=dict(color="orange", dash="dot"),
              annotation_text="α=0.05", annotation_font_color="rgba(234,234,234,0.7)")
fig.update_layout(
    title=dict(text="1 000-iter bootstrap of 4-group logrank p (median 2.6×10⁻⁵, 93 % iter p<0.05)",
               font=dict(size=14, color="#EAEAEA")),
    xaxis=dict(title="log₁₀(p)", **AXIS),
    yaxis=dict(title="Count", **AXIS),
    height=380, **DARK,
)
pio.write_html(fig, str(OUT_INT / "v17_bootstrap_p.html"), include_plotlyjs="cdn", full_html=True)
print(f"wrote v17_bootstrap_p.html")

# LOO
fig = go.Figure()
fig.add_trace(go.Bar(
    x=np.arange(len(loo)),
    y=-np.log10(loo["logrank_p"].clip(lower=1e-20)),
    marker=dict(color=["#c14" if e==1 else "#0d6765" for e in loo["is_event"]],
                line=dict(color="rgba(255,255,255,0.3)", width=0.4)),
    hovertemplate="iter %{x}<br>-log₁₀(p)=%{y:.2f}<extra></extra>",
))
fig.add_hline(y=-np.log10(0.05), line=dict(color="orange", dash="dot"),
              annotation_text="α=0.05", annotation_font_color="rgba(234,234,234,0.7)")
fig.add_hline(y=-np.log10(0.001), line=dict(color="red", dash="dash"),
              annotation_text="α=0.001", annotation_font_color="rgba(234,234,234,0.7)")
fig.update_layout(
    title=dict(text=f"Leave-one-out sensitivity — all 36 LOO iterations p<0.001 (worst case {A1C['p_max_loo']:.2g})",
               font=dict(size=14, color="#EAEAEA")),
    xaxis=dict(title="36 LOO iterations (sorted by p)", **AXIS),
    yaxis=dict(title="-log₁₀(logrank p)", **AXIS),
    height=380, **DARK,
)
pio.write_html(fig, str(OUT_INT / "v17_loo_p.html"), include_plotlyjs="cdn", full_html=True)
print(f"wrote v17_loo_p.html")
