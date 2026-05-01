"""B-prompt: 4-panel figure for the 8-cell mutation × outcome matrix.

Panels:
  A: KM curves for 4 collapsed groups (BRAF_only, RAS_only, TERT+, Triple_neg)
  B: forest plot of HR (8-cell vs OTHER_TERT-) with bootstrap CI
  C: heatmap of N + event rate per cell
  D: TERT+ subgroup composition pie / stacked bar
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/audit_2026_04_29/4way"

forest = pd.read_csv(OUT / "forest_HR_8cell.tsv", sep="\t")
crosstab = pd.read_csv(OUT / "8cell_crosstab.tsv", sep="\t")
sub = pd.read_csv(OUT / "tert_subgroup_breakdown.tsv", sep="\t")

SM_PATH = ROOT / "project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv"
sm = pd.read_csv(SM_PATH, sep="\t", low_memory=False)
sm = sm[(sm["dataset"] == "TCGA-THCA") & (sm["normal_vs_tumor"] == "tumor")].copy()
raw = sm["tert_promoter_integrated"].fillna("").astype(str).str.lower().str.strip()
sm["tert_int"] = (raw == "mutated").astype(int)
sm["os_event"] = pd.to_numeric(sm["os_event"], errors="coerce")
sm["os_days"] = pd.to_numeric(sm["os_days"], errors="coerce")
drv = sm["driver_anchor"].fillna("unknown").astype(str).str.upper()
drv = drv.where(drv.isin(["BRAF", "RAS", "NTRK"]), other="OTHER")
sm["driver_simple"] = drv
sm["cell"] = sm.apply(
    lambda r: f"{r['driver_simple']}_TERT+" if r["tert_int"] == 1 else f"{r['driver_simple']}_TERT-", axis=1
)

def four_group(r):
    if r["tert_int"] == 1:
        return "TERT+ (any)"
    if r["driver_simple"] == "BRAF":
        return "BRAF only"
    if r["driver_simple"] == "RAS":
        return "RAS only"
    return "Triple_neg"

sm["four_group"] = sm.apply(four_group, axis=1)
surv = sm.dropna(subset=["os_event", "os_days"])

fig = plt.figure(figsize=(14, 11))
gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.30)

# Panel A: KM 4-group
ax = fig.add_subplot(gs[0, 0])
colors = {"BRAF only": "#d62728", "RAS only": "#2ca02c", "TERT+ (any)": "#1f77b4", "Triple_neg": "#7f7f7f"}
for grp, sub_df in surv.groupby("four_group"):
    kmf = KaplanMeierFitter()
    kmf.fit(sub_df["os_days"], sub_df["os_event"], label=f"{grp} (n={len(sub_df)}, e={int(sub_df['os_event'].sum())})")
    kmf.plot_survival_function(ax=ax, ci_show=False, color=colors.get(grp, "k"))
ax.set_title("A. Overall survival, 4 collapsed groups (TCGA-THCA)", fontsize=11, fontweight="bold")
ax.set_xlabel("Days"); ax.set_ylabel("Survival probability")
ax.set_ylim(0.85, 1.005); ax.legend(loc="lower left", fontsize=8)

# Panel B: Forest plot 8-cell
ax = fig.add_subplot(gs[0, 1])
forest_plot = forest.dropna(subset=["hr_boot_median"]).copy()
forest_plot = forest_plot[forest_plot["n_target"] >= 3]  # drop NTRK_TERT+ n=1
forest_plot = forest_plot.sort_values("hr_boot_median")
y = np.arange(len(forest_plot))
ax.errorbar(
    forest_plot["hr_boot_median"], y,
    xerr=[forest_plot["hr_boot_median"] - forest_plot["hr_boot_ci_lo"],
          forest_plot["hr_boot_ci_hi"] - forest_plot["hr_boot_median"]],
    fmt="s", capsize=4, color="black", ecolor="gray"
)
ax.axvline(1, ls="--", color="red", alpha=0.5)
ax.set_yticks(y); ax.set_yticklabels(forest_plot["cell"].tolist())
ax.set_xscale("log"); ax.set_xlabel("Hazard Ratio (vs OTHER_TERT-, log scale)")
ax.set_title("B. Bootstrap HR (1000×), 8-cell vs Triple_neg (OTHER_TERT-)", fontsize=11, fontweight="bold")
for i, row in forest_plot.reset_index(drop=True).iterrows():
    ax.text(row["hr_boot_ci_hi"] * 1.1, i, f"n={row['n_target']}, e={row['events_target']}", va="center", fontsize=8)

# Panel C: Cell heatmap
ax = fig.add_subplot(gs[1, 0])
ct = crosstab.copy()
ct["driver"] = ct["cell"].str.split("_").str[0]
ct["tert"] = ct["cell"].str.contains(r"\+").map({True: "TERT+", False: "TERT-"})
mat = ct.pivot(index="driver", columns="tert", values="event_rate_pct").reindex(["BRAF", "RAS", "NTRK", "OTHER"])
mat = mat.reindex(columns=["TERT-", "TERT+"])
n_mat = ct.pivot(index="driver", columns="tert", values="n").reindex(["BRAF", "RAS", "NTRK", "OTHER"]).reindex(columns=["TERT-", "TERT+"])
im = ax.imshow(mat.values, cmap="Reds", vmin=0, vmax=30, aspect="auto")
ax.set_xticks(range(2)); ax.set_xticklabels(["TERT-", "TERT+"])
ax.set_yticks(range(4)); ax.set_yticklabels(["BRAF", "RAS", "NTRK", "OTHER"])
for i in range(4):
    for j in range(2):
        n = n_mat.iloc[i, j]
        rate = mat.iloc[i, j]
        if pd.notna(n):
            txt = f"N={int(n)}\n{rate:.1f}%"
            ax.text(j, i, txt, ha="center", va="center", color="black" if rate < 15 else "white", fontsize=9)
plt.colorbar(im, ax=ax, label="OS event rate (%)")
ax.set_title("C. 8-cell N and event rate (driver × TERT)", fontsize=11, fontweight="bold")

# Panel D: TERT+ subgroup composition
ax = fig.add_subplot(gs[1, 1])
sub_sorted = sub.sort_values("n", ascending=True)
colors_d = {"BRAF": "#d62728", "RAS": "#2ca02c", "NTRK": "#9467bd", "OTHER": "#7f7f7f"}
bars = ax.barh(sub_sorted["driver"], sub_sorted["n"], color=[colors_d[d] for d in sub_sorted["driver"]])
for bar, frac, ev in zip(bars, sub_sorted["frac_of_TERT+"], sub_sorted["events"]):
    ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height()/2,
            f"{frac:.1f}% (e={ev})", va="center", fontsize=9)
ax.set_xlabel("Number of TERT+ patients (total = 36)")
ax.set_title("D. TERT+ population is dominated by BRAF+\n(69% are BRAF+, not 'triple-neg')", fontsize=11, fontweight="bold")
ax.set_xlim(0, 32)

fig.suptitle(
    "TERT × BRAF × RAS 8-cell re-validation, TCGA-THCA n=504\n"
    "Resolution of the 'TERT-only triple-neg-ish is worst' paradox raised by Prof. Yu (2026-04-29 AM)",
    fontsize=12, fontweight="bold", y=0.995
)
fig.savefig(OUT / "figure_4way_revalidation.pdf", bbox_inches="tight", dpi=200)
fig.savefig(OUT / "figure_4way_revalidation.png", bbox_inches="tight", dpi=150)
print(f"Saved: {OUT/'figure_4way_revalidation.pdf'}")
print(f"Saved: {OUT/'figure_4way_revalidation.png'}")
