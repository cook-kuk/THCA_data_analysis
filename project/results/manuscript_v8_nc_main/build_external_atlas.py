"""Build the External Validation Atlas — one composite figure across all 15+ cohorts.

Layout: 4 rows × 4 cols = 16 panels (15 cohorts + 1 summary score-card).
Each cell summarizes one cohort with: name, n, validation layer (RNA/HM450/sc/protein/drug),
direction statistic (Cohen's d / Spearman ρ / OR / HR), and significance.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mp
import numpy as np, pandas as pd, os

ROOT = "/home/seungho/personal/THCA_data_analysis/project/results"
OUT  = ROOT + "/manuscript_v8_nc_main"
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.titlesize":9.5,"axes.labelsize":8.5,
                     "axes.spines.top":False,"axes.spines.right":False,"legend.fontsize":7.5})

# ---- per-cohort scorecard ----
# Layer colors
layer_color = {
    "Bulk RNA":   "#4c72b0",
    "Methylation":"#8172b3",
    "Single-cell":"#55a868",
    "Protein":    "#dd8452",
    "Drug":       "#937860",
    "Spatial":    "#c44e52",
    "Pan-cancer": "#c44e52",
    "Discovery":  "#999",
}
cohorts = [
    # (label, n, layer, stat_label, stat_text, significance, color_override)
    ("TCGA-THCA primary",   504, "Bulk RNA", "Discovery",         "DM1 prevalence 28.4 %",   "anchor cohort",            "Discovery"),
    ("Lee  GSE213647",      632, "Bulk RNA", "Cohen's d (DM1−DM2)","d = +0.95",              "p < 10⁻¹⁵",                None),
    ("K2  PRJEB11591",      260, "Bulk RNA", "% DM1-like",        "37.8 %  (calibration caveat)","K-S preserved",          None),
    ("MSK-IMPACT advanced", 117, "Bulk RNA", "Cox HR (DM1 vs DM2)","HR = 2.67 [1.17, 6.10]", "p = 0.018",                None),
    ("GSE286332  PTC+HT",    18, "Bulk RNA", "MAPK × Panel ρ",    "ρ = −0.040  (HT route)",  "NS (predicted)",           None),
    ("GSE76039  PDTC+ATC",   37, "Bulk RNA", "Convergence",       "5 / 8 panel overlap",     "Landa silenced list",      None),
    ("GSE126698  Korean PTC",28, "Bulk RNA", "MAPK × Panel ρ",    "ρ = −0.18",               "trend",                    None),
    ("HM450  TCGA",         503, "Methylation","TPO d (DM1 vs DM2)","d = 2.30",              "p = 1.9 × 10⁻¹⁸",          None),
    ("Pu 2021  GSE184362",    6, "Single-cell","Per-patient r",   "r = 0.798 – 0.886",       "p < 10⁻¹⁰  Bonferroni",    None),
    ("Lu 2023  GSE193581",   23, "Single-cell","Spearman r",      "r = 0.86 median",         "p < 10⁻¹⁰",                None),
    ("GSE241184  Phase 1",   14, "Single-cell","Author-independence","r = 0.79",            "concordant",               None),
    ("GSE232237  External",  11, "Single-cell","Pseudobulk d",    "d = +1.65",               "p = 4 × 10⁻⁹",             None),
    ("Mun 2025 proteome",   336, "Protein",  "Thyroid-diff d",     "d = −1.91",              "7 / 7 sign-consistent",    None),
    ("GSE151179 post-RAI",   52, "Bulk RNA", "Pre vs post-RAI d", "d = −1.01",               "p = 1 × 10⁻⁴",             None),
    ("DepMap + PRISM",      "~766","Drug",    "AZD-0364 MEK d",   "d = −0.594",              "FDR = 5.9 × 10⁻⁷",         None),
    ("TCGA pan-cancer",   11069, "Pan-cancer","LGG / LUAD HR",    "HR 44.7 / 19.0",          "p = 1.3 × 10⁻⁵ / 8.6 × 10⁻⁴", None),
]

fig = plt.figure(figsize=(20, 13), dpi=150, facecolor="white")
gs  = fig.add_gridspec(4, 4, hspace=0.55, wspace=0.20, left=0.03, right=0.985, top=0.93, bottom=0.05)

for i, (label, n, layer, stat_label, stat_text, sig, color_override) in enumerate(cohorts):
    r, c = divmod(i, 4)
    ax = fig.add_subplot(gs[r, c])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_axis_off()
    col = layer_color.get(color_override or layer, "#999")
    # Card background
    ax.add_patch(mp.FancyBboxPatch((0.02, 0.04), 0.96, 0.92, boxstyle="round,pad=0.012", fc="#10141c", ec=col, lw=2))
    # Layer chip
    ax.add_patch(mp.FancyBboxPatch((0.05, 0.81), 0.32, 0.13, boxstyle="round,pad=0.01", fc=col, ec="none"))
    ax.text(0.21, 0.875, layer, fontsize=8.5, color="white", ha="center", va="center", fontweight="bold")
    # n chip
    ax.add_patch(mp.FancyBboxPatch((0.62, 0.81), 0.33, 0.13, boxstyle="round,pad=0.01", fc="#252b39", ec="none"))
    n_clean = str(n).replace('~','').replace('≈','').replace(',','').strip()
    if n_clean.isdigit():
        n_int = int(n_clean)
        n_str = f"n ≈ {n_int:,}" if str(n).startswith("~") else f"n = {n_int:,}"
    else:
        n_str = str(n)
    ax.text(0.785, 0.875, n_str, fontsize=8.5, color="white", ha="center", va="center", fontweight="bold")
    # Title (cohort label)
    ax.text(0.05, 0.69, label, fontsize=10.5, color="white", fontweight="bold")
    # Stat label
    ax.text(0.05, 0.55, stat_label, fontsize=8.5, color="#9aa4b2", style="italic")
    # Big stat
    ax.text(0.05, 0.36, stat_text, fontsize=13, color="#ffd166", fontweight="bold")
    # Significance line
    ax.text(0.05, 0.15, sig, fontsize=9, color="#7bd88f" if "p =" in sig or "FDR" in sig or "p <" in sig else "#bcc4d0")

fig.suptitle("External Validation Atlas · 15+ cohorts confirming the DM1/DM2 axis across bulk, methylation, single-cell, protein, drug, and pan-cancer layers",
             x=0.025, y=0.975, fontsize=14, ha="left", fontweight="bold", color="#222")
fig.text(0.025, 0.95, "Direction-consistent across all 14/15 external cohorts (GSE286332 is the predicted HT-route attenuation per the two-axis model)",
         fontsize=10, ha="left", color="#666", style="italic")
out = f"{OUT}/External_Validation_Atlas.png"
fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
fig.savefig(out.replace(".png",".pdf"), bbox_inches="tight")
plt.close(fig)
print(f"saved {out}")

# ---- also export a scorecard TSV for the HTML table ----
df = pd.DataFrame([{
    "cohort": c[0], "n": str(c[1]), "layer": c[2], "metric": c[3], "value": c[4], "significance": c[5]
} for c in cohorts])
df["main_fig_location"] = ["Fig 1 anchor","Fig 4G + Fig 5E","Fig 5E + ED14","Fig 2G + Fig 5A","Fig 3D",
                            "Fig 3D + Fig 3G","Fig 3D","Fig 3A + ED1","Fig 4B","Fig 4A + Fig 4F",
                            "Fig 4C","Fig 4D","ED11","Fig 6D","ED10","ED12"]
df.to_csv(f"{OUT}/cohort_scorecard.tsv", sep="\t", index=False)
print(f"scorecard saved {OUT}/cohort_scorecard.tsv")
print(df.to_string(index=False))
