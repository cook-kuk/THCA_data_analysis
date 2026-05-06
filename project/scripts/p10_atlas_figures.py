#!/usr/bin/env python3
"""Paper 10 — DM1 Synthetic-Lethal Vulnerability Atlas figure generator.

Reads existing TSVs from project/results/paper9_sl_first_pass/ and emits ~17
atlas-grade figures into project/results/paper10_atlas/figures/.
Adds: volcano, class facets, heatmap, forest, KM curves, class radar,
lineage-TF anchor panel, triangulation scoreboard.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy import stats

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
P9 = ROOT / "project/results/paper9_sl_first_pass"
OUT = ROOT / "project/results/paper10_atlas"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# Class palette — matched to hub teal/accent
CLASS_COLORS = {
    "NAD_salvage":             "#7B1F2A",
    "NAD_de_novo_SR_partner":  "#A04451",
    "JAK_STAT":                "#34547A",
    "JAK_STAT_upstream":       "#5E80B0",
    "Epigenetic":              "#B8893C",
    "SFK":                     "#3C6B4F",
    "Ion_channel":             "#3F7A8A",
    "DDR":                     "#962E2E",
    "Metabolism":              "#8B5E00",
    "Metabolism_glycolysis":   "#A5722B",
    "Metabolism_master":       "#5C3D00",
    "TROP2_ADC_non_SL":        "#5E5E5E",
    "Lineage_TF_anchor":       "#0F1A2E",
}
LINEAGE_TFS = ["FOXE1", "NKX2-1", "PAX8", "HHEX"]
RAI_8 = ["TPO", "DIO1", "TSHR", "PAX8", "TG", "FOXE1", "NKX2-1", "SLC5A5"]

# ====================================================================
# Load all P9 outputs
# ====================================================================
candidates = pd.read_csv(P9 / "candidate_targets.tsv", sep="\t")
ccle_state = pd.read_csv(P9 / "ccle_dm1_state.tsv", sep="\t")
ccle_corr = pd.read_csv(P9 / "ccle_dm1_target_corr.tsv", sep="\t")
dge = pd.read_csv(P9 / "dge_candidate_ranking.tsv", sep="\t")
cox_full = pd.read_csv(P9 / "tcga_thca_target_cox.tsv", sep="\t")
cox_dm1 = pd.read_csv(P9 / "tcga_thca_dm1_target_cox.tsv", sep="\t")
summary = json.loads((P9 / "summary.json").read_text())
print(f"Loaded P9 outputs: {len(candidates)} candidates · "
      f"{len(ccle_state)} CCLE lines · {len(dge)} DGE rows · "
      f"{len(cox_full)} cox(full) · {len(cox_dm1)} cox(dm1)")


def cmap(c, default="#888"):
    return CLASS_COLORS.get(c, default)


# ====================================================================
# F02 — CCLE n=13 DM1 score distribution
# ====================================================================
fig, ax = plt.subplots(figsize=(11, 4.8))
df = ccle_state.sort_values("DM1_like_score").copy()
ax.barh(df["line"], df["DM1_like_score"],
        color=["#7B1F2A" if v > 0 else "#34547A" for v in df["DM1_like_score"]])
ax.axvline(0, color="black", linewidth=0.7)
ax.axvline(df["DM1_like_score"].median(), color="#B8893C",
           linestyle="--", linewidth=1, label="median split")
ax.set_xlabel("DM1_like score (within-sample z, mean of 8 RAI genes, sign-flipped)")
ax.set_title("F02 — CCLE thyroid n=13: DM1_like state distribution\n"
             "(red bars = DM1-high; blue bars = DM1-low)")
ax.legend(loc="lower right")
plt.tight_layout()
plt.savefig(FIG / "F02_ccle_dm1_distribution.png", dpi=140)
plt.close()

# ====================================================================
# F03 — CCLE corr per-class facet
# ====================================================================
classes_in_ccle = ccle_corr.dropna(subset=["pearson_r"])["class"].unique()
n_cls = len(classes_in_ccle)
ncols = 3
nrows = int(np.ceil(n_cls / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(13, 3.3 * nrows), squeeze=False)
for i, cls in enumerate(sorted(classes_in_ccle)):
    ax = axes[i // ncols][i % ncols]
    sub = ccle_corr[ccle_corr["class"] == cls].dropna(subset=["pearson_r"])
    ax.barh(sub["gene"], sub["pearson_r"], color=cmap(cls))
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_xlim(-1, 1)
    ax.set_title(cls, fontsize=11)
    ax.set_xlabel("Pearson r vs DM1_like (CCLE n=13)")
for j in range(i + 1, nrows * ncols):
    axes[j // ncols][j % ncols].axis("off")
plt.suptitle("F03 — CCLE thyroid Pearson r vs DM1_like, per target class",
             fontsize=13, y=1.0)
plt.tight_layout()
plt.savefig(FIG / "F03_ccle_corr_class_facet.png", dpi=140)
plt.close()

# ====================================================================
# F04 — DGE volcano of all 36 candidates
# ====================================================================
fig, ax = plt.subplots(figsize=(10, 8))
dge_v = dge.dropna(subset=["log2FC_DM2vDM1", "q_BH"]).copy()
dge_v["dm1_minus_dm2"] = -dge_v["log2FC_DM2vDM1"]
dge_v["mlogq"] = -np.log10(dge_v["q_BH"].replace(0, 1e-300))
for cls, color in CLASS_COLORS.items():
    sub = dge_v[dge_v["class"] == cls]
    if len(sub):
        ax.scatter(sub["dm1_minus_dm2"], sub["mlogq"], c=color, s=70,
                   label=cls, edgecolor="white", linewidth=0.8)
        for _, r in sub.iterrows():
            ax.annotate(r["gene"], (r["dm1_minus_dm2"], r["mlogq"]),
                        fontsize=8, xytext=(4, 3),
                        textcoords="offset points", color=color)
ax.axhline(-np.log10(0.05), color="grey", linestyle="--", linewidth=0.8,
           label="q=0.05")
ax.axvline(0, color="black", linewidth=0.6)
ax.set_xlabel("log2(DM1 / DM2)  ←  DM2 enriched   |   DM1 enriched  →")
ax.set_ylabel("−log10(q_BH)")
ax.set_title("F04 — DM1 vs DM2 DGE volcano (TCGA-THCA bulk, all 36 candidates)\n"
             "(genes labeled; classes colored)")
ax.legend(loc="lower left", fontsize=8, ncol=2, framealpha=0.9)
plt.tight_layout()
plt.savefig(FIG / "F04_dge_volcano.png", dpi=140)
plt.close()

# ====================================================================
# F05 — DGE class summary boxplot
# ====================================================================
fig, ax = plt.subplots(figsize=(11, 5.5))
dge_v_sorted = dge_v.copy()
dge_v_sorted["class"] = dge_v_sorted["class"].astype(str)
class_order = (dge_v_sorted.groupby("class")["dm1_minus_dm2"].mean()
               .sort_values().index.tolist())
data = [dge_v_sorted[dge_v_sorted["class"] == c]["dm1_minus_dm2"].values
        for c in class_order]
bp = ax.boxplot(data, labels=class_order, patch_artist=True, widths=0.6,
                medianprops={"color": "black", "linewidth": 1.5})
for patch, c in zip(bp["boxes"], class_order):
    patch.set_facecolor(cmap(c))
    patch.set_alpha(0.85)
ax.axhline(0, color="black", linewidth=0.6)
ax.set_ylabel("log2(DM1/DM2) per gene in class")
ax.set_xticklabels(class_order, rotation=30, ha="right", fontsize=10)
ax.set_title("F05 — Per-class DGE distribution (DM1 vs DM2, TCGA-THCA bulk)")
plt.tight_layout()
plt.savefig(FIG / "F05_dge_class_box.png", dpi=140)
plt.close()

# ====================================================================
# F07 — TCGA-THCA Cox HR forest with class colors (whole cohort)
# ====================================================================
fig, ax = plt.subplots(figsize=(10, 11))
df = cox_full.dropna(subset=["hr"]).copy()
df["class"] = df["class"].fillna("Unknown")
df = df.sort_values("hr")
y = np.arange(len(df))
log_hr = np.log(df["hr"])
log_lo = np.log(df["hr_lo"])
log_hi = np.log(df["hr_hi"])
colors = df["class"].map(CLASS_COLORS).fillna("#888").values
ax.errorbar(log_hr, y, xerr=[log_hr - log_lo, log_hi - log_hr],
            fmt="o", ecolor="grey", capsize=3, elinewidth=0.8,
            markersize=8, mfc="white", mec="black", zorder=3)
for yi, c in zip(y, colors):
    ax.scatter([log_hr.iloc[yi]], [yi], c=c, s=60, zorder=4)
ax.axvline(0, color="black", linewidth=0.6)
ax.set_yticks(y)
ax.set_yticklabels([f"{g} ({cls})"
                    for g, cls in zip(df["gene"], df["class"])], fontsize=9)
ax.set_xlabel("log HR (per +1 expression unit, OS, lifelines Cox)")
ax.set_title("F07 — TCGA-THCA univariate Cox HR per Paper 10 candidate\n"
             "(whole cohort, n=560, 18 OS events)")
# Significance markers
for yi, p in zip(y, df["p"]):
    if pd.notna(p) and p < 0.05:
        ax.text(log_hi.iloc[yi] + 0.1, yi, " ✱  p<0.05",
                va="center", fontsize=9, color="#962E2E")
plt.tight_layout()
plt.savefig(FIG / "F07_cox_forest_full.png", dpi=140)
plt.close()

# ====================================================================
# F08 — Whole vs DM1-high HR comparison
# ====================================================================
merged = cox_full[["gene", "class", "hr", "p"]].rename(
    columns={"hr": "hr_whole", "p": "p_whole"}
).merge(
    cox_dm1[["gene", "hr", "p"]].rename(
        columns={"hr": "hr_dm1", "p": "p_dm1"}),
    on="gene", how="inner"
)
fig, ax = plt.subplots(figsize=(8, 8))
log_w = np.log(merged["hr_whole"])
log_d = np.log(merged["hr_dm1"])
colors = merged["class"].map(CLASS_COLORS).fillna("#888").values
ax.scatter(log_w, log_d, c=colors, s=80, edgecolor="white", linewidth=0.8)
lim = max(abs(log_w).max(), abs(log_d).max()) * 1.05
ax.plot([-lim, lim], [-lim, lim], color="grey", linestyle="--",
        linewidth=0.8, label="y = x")
ax.axhline(0, color="black", linewidth=0.4)
ax.axvline(0, color="black", linewidth=0.4)
for _, r in merged.iterrows():
    ax.annotate(r["gene"], (np.log(r["hr_whole"]), np.log(r["hr_dm1"])),
                fontsize=8, xytext=(4, 3), textcoords="offset points")
ax.set_xlabel("log HR (whole TCGA-THCA, n=560)")
ax.set_ylabel("log HR (DM1-high subset, n=280)")
ax.set_title("F08 — Cox HR consistency: whole cohort vs DM1-high subset\n"
             "(diagonal = invariant to subsetting; deviation = state-conditional)")
ax.legend(loc="lower right")
plt.tight_layout()
plt.savefig(FIG / "F08_cox_whole_vs_dm1high.png", dpi=140)
plt.close()

# ====================================================================
# F09 — Class × layer triangulation heatmap
# ====================================================================
classes_all = sorted(set(c for c in CLASS_COLORS) - {"NAD_de_novo_SR_partner"})
layers = ["DGE log2FC (DM1↑)", "CCLE r (DM1↑)", "Cox -log10(p)"]
mat = np.zeros((len(classes_all), len(layers)))
for i, cls in enumerate(classes_all):
    sub_dge = dge[dge["class"] == cls]
    sub_ccle = ccle_corr[ccle_corr["class"] == cls].dropna(subset=["pearson_r"])
    sub_cox = cox_full[cox_full["class"] == cls].dropna(subset=["p"])
    mat[i, 0] = -sub_dge["log2FC_DM2vDM1"].mean() if len(sub_dge) else np.nan
    mat[i, 1] = sub_ccle["pearson_r"].mean() if len(sub_ccle) else np.nan
    mat[i, 2] = (-np.log10(sub_cox["p"])).mean() if len(sub_cox) else np.nan
fig, ax = plt.subplots(figsize=(7.5, 7))
# Standardize each column for visual comparability
mat_z = np.copy(mat)
for c in range(mat_z.shape[1]):
    col = mat_z[:, c]
    finite = col[np.isfinite(col)]
    if len(finite) > 1:
        mu, sd = finite.mean(), finite.std()
        if sd > 0:
            mat_z[:, c] = (col - mu) / sd
im = ax.imshow(mat_z, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
ax.set_xticks(range(len(layers)))
ax.set_xticklabels(layers, rotation=20, ha="right")
ax.set_yticks(range(len(classes_all)))
ax.set_yticklabels(classes_all)
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        v = mat[i, j]
        if np.isfinite(v):
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center",
                    fontsize=8.5, color="black")
plt.colorbar(im, ax=ax, label="z-score across classes (per layer)")
ax.set_title("F09 — Triangulation matrix: target class × evidence layer\n"
             "(text = raw value; color = within-layer z-score)")
plt.tight_layout()
plt.savefig(FIG / "F09_triangulation_matrix.png", dpi=140)
plt.close()

# ====================================================================
# F10 — KM curves for top 5 candidates (TCGA-THCA OS)
# ====================================================================
print("\nBuilding KM curves …")
scored = pd.read_csv(
    ROOT / "project_external_st/results/extra/s_tcga_thca_scored.tsv", sep="\t"
).dropna(subset=["os_event", "os_days", "DM1_like"]).copy()
scored["os_event"] = scored["os_event"].astype(int)
scored["os_days"] = scored["os_days"].astype(float)

pancan_path = ROOT / "project/data/raw/TCGA_pancan/pancan_geneExp.gz"
top5_genes = cox_full.dropna(subset=["p"]).sort_values("p").head(5)["gene"].tolist()
print(f"Top 5 by p: {top5_genes}")

# Stream-load just these 5 genes
with gzip.open(pancan_path, "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    sample_cols = header[1:]
    rows = []
    want = set(top5_genes)
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if parts[0] in want:
            rows.append(parts)
            if len(rows) == len(want):
                break
top5_df = pd.DataFrame(rows, columns=header).set_index(header[0])
top5_df = top5_df.replace(["", "NA", "na", "NaN"], np.nan).astype(float)

pat_to_cols = {}
for c in top5_df.columns:
    pat = "-".join(c.split("-")[:3])
    pat_to_cols.setdefault(pat, []).append(c)

def pick_tumor(cols):
    tumor = [c for c in cols if c.split("-")[3].startswith("01")]
    return tumor[0] if tumor else (cols[0] if cols else None)

scored["sample_col"] = scored["patient"].map(
    lambda p: pick_tumor(pat_to_cols.get(p, []))
)
matched = scored.dropna(subset=["sample_col"]).copy()

try:
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test
    HAS_LL = True
except ImportError:
    HAS_LL = False

if HAS_LL and len(top5_genes):
    fig, axes = plt.subplots(1, len(top5_genes),
                              figsize=(4.2 * len(top5_genes), 4.5),
                              sharey=True)
    if len(top5_genes) == 1:
        axes = [axes]
    for ax, g in zip(axes, top5_genes):
        if g not in top5_df.index:
            ax.set_title(f"{g} (n/a)")
            continue
        g_expr = top5_df.loc[g, matched["sample_col"].values].values.astype(float)
        df = pd.DataFrame({
            "expr": g_expr,
            "os_event": matched["os_event"].astype(int).values,
            "os_days": matched["os_days"].astype(float).values,
        }).dropna()
        med = df["expr"].median()
        df["high"] = df["expr"] > med
        kmf = KaplanMeierFitter()
        for label, mask, color in [
            ("high (above median)", df["high"], "#7B1F2A"),
            ("low (below median)", ~df["high"], "#34547A"),
        ]:
            sub = df[mask]
            kmf.fit(sub["os_days"], event_observed=sub["os_event"], label=label)
            kmf.plot_survival_function(ax=ax, color=color, ci_show=False)
        lr = logrank_test(df.loc[df["high"], "os_days"],
                          df.loc[~df["high"], "os_days"],
                          event_observed_A=df.loc[df["high"], "os_event"],
                          event_observed_B=df.loc[~df["high"], "os_event"])
        ax.set_title(f"{g}\nlogrank p = {lr.p_value:.3f}", fontsize=11)
        ax.set_xlabel("OS days")
        ax.set_ylim(0.6, 1.02)
    fig.suptitle("F10 — TCGA-THCA OS Kaplan-Meier for top-5 Cox candidates "
                 "(median-split)", y=1.0)
    plt.tight_layout()
    plt.savefig(FIG / "F10_top5_KM.png", dpi=140)
    plt.close()
    print("Wrote F10_top5_KM.png")

# ====================================================================
# F11 — Triangulation scoreboard (figure-as-table)
# ====================================================================
fig, ax = plt.subplots(figsize=(10, 7.5))
ax.axis("off")
score_rows = [
    ["Rank", "Class / Target", "DGE", "CCLE", "Cox", "Drug", "Verdict"],
    ["1", "KCNN4 (Ion channel)", "↑↑↑ q=4e-32",
     "+0.29 (n=13)", "ns", "senicapoc / TRAM-34",
     "lead candidate"],
    ["2", "OSMR + IL6R (JAK/STAT upstream)", "↑↑ q<1e-13",
     "−0.32 OSMR (n=13)", "ns", "anti-IL6R / dual-JAK",
     "receptor-blockade SL"],
    ["3", "SFK: LYN, SRC, FYN", "↑↑ all q<<0.05",
     "+0.28 LYN (n=13)", "ns trend", "dasatinib-class",
     "consistent class signal"],
    ["4", "DDR: ATR + CHEK2", "flat",
     "n/a", "ATR HR=6.88 p=0.024",
     "ceralasertib + olaparib", "state-conditional"],
    ["5", "MYC + GLS / LDHA", "↑ MYC q=2e-9",
     "n/a", "MYC HR=1.50 trend", "CB-839 / IACS-010759",
     "metabolic SL pair"],
    ["6", "NAMPT / NAPRT", "ns",
     "+0.36 NAMPT (n=13)", "ns",
     "FK866 / KPT-9274", "needs DepMap"],
    ["7", "TROP2 (TACSTD2)", "↑↑↑ q=8e-24",
     "+0.29 (n=13)", "ns",
     "sacituzumab govitecan", "non-SL ADC carved-out"],
]
table = ax.table(cellText=score_rows, loc="center", cellLoc="left")
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 1.7)
# Color top row
for j in range(7):
    table[(0, j)].set_facecolor("#0F1A2E")
    table[(0, j)].set_text_props(color="white", weight="bold")
# Color rank cells
rank_colors = ["#3C6B4F", "#3C6B4F", "#3C6B4F", "#B8893C",
               "#B8893C", "#A5722B", "#5E5E5E"]
for i in range(1, 8):
    table[(i, 0)].set_facecolor(rank_colors[i - 1])
    table[(i, 0)].set_text_props(color="white", weight="bold")
ax.set_title("F11 — Paper 10 ATLAS triangulation scoreboard\n"
             "(verdicts derived from DGE + CCLE + Cox first-pass triangulation)",
             pad=18)
plt.tight_layout()
plt.savefig(FIG / "F11_triangulation_scoreboard.png", dpi=140)
plt.close()

# ====================================================================
# F12 — Class radar / spider plot (3-axis triangulation)
# ====================================================================
focus_classes = [
    "Ion_channel", "JAK_STAT_upstream", "SFK", "DDR",
    "Metabolism_master", "NAD_salvage", "TROP2_ADC_non_SL",
]
axes_labels = ["DGE\n(DM1↑ effect)", "CCLE\n(DM1 corr)", "Cox\n(-log10 p)"]
n_axes = len(axes_labels)
angles = np.linspace(0, 2 * np.pi, n_axes, endpoint=False).tolist()
angles += angles[:1]

fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={"polar": True})

# Per-class min-max normalization across focus_classes
metrics = []
for cls in focus_classes:
    s_dge = dge[dge["class"] == cls]
    s_ccle = ccle_corr[ccle_corr["class"] == cls].dropna(subset=["pearson_r"])
    s_cox = cox_full[cox_full["class"] == cls].dropna(subset=["p"])
    metrics.append([
        -s_dge["log2FC_DM2vDM1"].mean() if len(s_dge) else 0,
        s_ccle["pearson_r"].mean() if len(s_ccle) else 0,
        (-np.log10(s_cox["p"])).mean() if len(s_cox) else 0,
    ])
metrics = np.array(metrics)
# Min-max normalize each axis
metrics_norm = np.zeros_like(metrics)
for c in range(metrics.shape[1]):
    col = metrics[:, c]
    lo, hi = col.min(), col.max()
    metrics_norm[:, c] = (col - lo) / (hi - lo + 1e-12)

for cls, vals in zip(focus_classes, metrics_norm):
    v = vals.tolist() + [vals[0]]
    ax.plot(angles, v, color=cmap(cls), linewidth=2, label=cls)
    ax.fill(angles, v, color=cmap(cls), alpha=0.12)
ax.set_xticks(angles[:-1])
ax.set_xticklabels(axes_labels, fontsize=10)
ax.set_yticks([0.25, 0.5, 0.75, 1.0])
ax.set_yticklabels([])
ax.set_title("F12 — Class triangulation radar\n"
             "(per-axis min-max normalized across focus classes)", pad=20)
ax.legend(loc="upper right", bbox_to_anchor=(1.45, 1.05), fontsize=9)
plt.tight_layout()
plt.savefig(FIG / "F12_class_radar.png", dpi=140)
plt.close()

# ====================================================================
# F13 — Lineage TF anchor 4-panel (sanity floor)
# ====================================================================
fig, axes = plt.subplots(1, 4, figsize=(15, 4.2))
for ax, tf in zip(axes, LINEAGE_TFS):
    sub_ccle = ccle_corr[ccle_corr["gene"] == tf]
    sub_dge = dge[dge["gene"] == tf]
    r = sub_ccle["pearson_r"].iloc[0] if len(sub_ccle) else np.nan
    p = sub_ccle["p"].iloc[0] if len(sub_ccle) else np.nan
    log2fc = -sub_dge["log2FC_DM2vDM1"].iloc[0] if len(sub_dge) else np.nan
    q = sub_dge["q_BH"].iloc[0] if len(sub_dge) else np.nan
    ax.bar([0, 1], [r, log2fc],
           color=["#34547A", "#7B1F2A"])
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["CCLE r", "DGE log2(DM1/DM2)"])
    ax.set_title(f"{tf}\nCCLE r={r:.2f} p={p:.1e}\n"
                 f"DGE log2FC={log2fc:+.2f} q={q:.1e}",
                 fontsize=10)
fig.suptitle("F13 — Lineage TF anchor sanity floor\n"
             "(all 4 lineage TFs DM2-enriched: confirms DM1 = "
             "lineage-collapsed state)", y=1.02)
plt.tight_layout()
plt.savefig(FIG / "F13_lineage_tf_anchor.png", dpi=140)
plt.close()

# ====================================================================
# F14 — Paper 9 vs Paper 10 positioning schematic
# ====================================================================
fig, ax = plt.subplots(figsize=(11, 5.5))
ax.axis("off")
ax.set_xlim(0, 10)
ax.set_ylim(0, 6)
# Paper 9
p9 = plt.Rectangle((0.4, 1.2), 4.2, 3.6, edgecolor="#3F7A8A",
                    facecolor="#E5EEF0", linewidth=2)
ax.add_patch(p9)
ax.text(2.5, 4.3, "PAPER 9", ha="center", fontsize=18, fontweight="bold",
        color="#3F7A8A", family="serif")
ax.text(2.5, 3.7, "Strategic plan + first-pass",
        ha="center", fontsize=12, style="italic", color="#3F7A8A")
ax.text(2.5, 3.0, "• 14-section text-heavy plan\n"
                  "• 8 target classes designed\n"
                  "• Local-data triangulation (B2/B3/B4)\n"
                  "• Honest deferral registry",
        ha="center", fontsize=10, color="#0F1A2E")

# Paper 10
p10 = plt.Rectangle((5.4, 1.2), 4.2, 3.6, edgecolor="#7B1F2A",
                     facecolor="#FAEEED", linewidth=2)
ax.add_patch(p10)
ax.text(7.5, 4.3, "PAPER 10", ha="center", fontsize=18, fontweight="bold",
        color="#7B1F2A", family="serif")
ax.text(7.5, 3.7, "Visualization-first ATLAS",
        ha="center", fontsize=12, style="italic", color="#7B1F2A")
ax.text(7.5, 3.0, "• 17 figures · 9 sortable tables\n"
                  "• Volcano · forest · radar · KM · heatmap\n"
                  "• Per-class triangulation scoreboard\n"
                  "• Drug-class landscape · provenance map",
        ha="center", fontsize=10, color="#0F1A2E")

# Arrow
ax.annotate("", xy=(5.4, 3.0), xytext=(4.6, 3.0),
            arrowprops=dict(arrowstyle="->", color="#0F1A2E", lw=2))
ax.text(5.0, 3.25, "feeds", ha="center", fontsize=10,
        style="italic", color="#0F1A2E")

ax.text(5.0, 0.6, "Both anchored on Paper 1 DM1/RAI-lineage axis · "
                   "Same data substrate · Different communication mode",
        ha="center", fontsize=11, style="italic", color="#3A4658")
ax.text(5.0, 5.6, "Paper 9 ↔ Paper 10 positioning",
        ha="center", fontsize=14, fontweight="bold", color="#0F1A2E")
plt.savefig(FIG / "F14_paper9_vs_paper10.png", dpi=140,
            bbox_inches="tight")
plt.close()

# ====================================================================
# F15 — Data provenance diagram
# ====================================================================
fig, ax = plt.subplots(figsize=(11, 6))
ax.axis("off")
ax.set_xlim(0, 10)
ax.set_ylim(0, 7)
ax.text(5, 6.5, "F15 — Paper 10 ATLAS data provenance",
        ha="center", fontsize=14, fontweight="bold")

sources = [
    (1.5, 5, "CCLE thyroid n=13\n(local thyroid panel)", "#3C6B4F"),
    (5.0, 5, "DM1-vs-DM2 DGE\n51,711 genes (TCGA-THCA bulk)", "#34547A"),
    (8.5, 5, "TCGA-THCA pancan\nlog2(norm+1) · 11,069 samples", "#7B1F2A"),
]
for x, y, txt, color in sources:
    box = plt.Rectangle((x - 1.2, y - 0.5), 2.4, 1.0,
                         edgecolor=color, facecolor="white", linewidth=2)
    ax.add_patch(box)
    ax.text(x, y, txt, ha="center", va="center", fontsize=9, color="#0F1A2E")

ax.text(5, 3.3, "Paper 10 ATLAS\n17 figures · 9 tables",
        ha="center", fontsize=14, fontweight="bold", color="#7B1F2A",
        bbox=dict(boxstyle="round,pad=0.4", edgecolor="#7B1F2A",
                  facecolor="#FAEEED", linewidth=2))

for x, _, _, color in sources:
    ax.annotate("", xy=(5, 3.7), xytext=(x, 4.5),
                arrowprops=dict(arrowstyle="->", color=color, lw=1.5))

deferred = [
    (1.5, 1.3, "DepMap CRISPR\n(genome-wide)", "#888"),
    (5.0, 1.3, "PRISM IC50\n(/opt/thyroid-dash host)", "#888"),
    (8.5, 1.3, "Wet-lab Tier 1\n(siRNA / IC50)", "#888"),
]
ax.text(5, 2.3, "DEFERRED — gated on Paper 1 in print + marathon close",
        ha="center", fontsize=10, style="italic", color="#888")
for x, y, txt, color in deferred:
    box = plt.Rectangle((x - 1.2, y - 0.4), 2.4, 0.8,
                         edgecolor=color, facecolor="#F5F5F5",
                         linewidth=1, linestyle="--")
    ax.add_patch(box)
    ax.text(x, y, txt, ha="center", va="center", fontsize=9, color="#666")

plt.savefig(FIG / "F15_data_provenance.png", dpi=140, bbox_inches="tight")
plt.close()

# ====================================================================
# F16 — Drug-class landscape
# ====================================================================
fig, ax = plt.subplots(figsize=(12, 6))
ax.axis("off")
landscape = [
    ["Class", "Top target(s)", "Drug class", "Lead compounds", "Stage"],
    ["NAD salvage", "NAMPT (SR: NAPRT)", "NAMPT inhibitor",
     "FK866 / APO866 · KPT-9274 · OT-82", "Phase I-II"],
    ["JAK/STAT upstream", "OSMR · IL6R",
     "Anti-receptor + JAK/TYK2 inhibitor",
     "tocilizumab · ruxolitinib · deucravacitinib", "Approved (other indications)"],
    ["JAK/STAT", "JAK1/2 · TYK2 · STAT3", "Selective JAK / STAT3 inhibitor",
     "ruxolitinib · TTI-101 · napabucasin*", "Approved + clinical"],
    ["Epigenetic", "DNMT1/3A/3B · HDAC1/2 · KDM1A · EZH2",
     "DNMT inhibitor + HDACi + LSD1i",
     "decitabine · guadecitabine · vorinostat · GSK-2879552 · tazemetostat",
     "Approved + clinical"],
    ["SFK", "LYN · SRC · FYN · YES1", "SRC-family inhibitor",
     "dasatinib · bosutinib · saracatinib", "Approved (CML/Ph+ALL)"],
    ["Ion channel", "KCNN4 (KCa3.1)", "K+ channel blocker",
     "senicapoc (clinical) · TRAM-34 (tool)", "Phase II (other indications)"],
    ["DDR", "ATR · CHEK2 · PARP1/2", "ATRi + PARPi (combination)",
     "ceralasertib · olaparib · talazoparib · niraparib",
     "Approved (BRCA) + clinical"],
    ["Metabolism", "GLS · LDHA · IDH2 · MYC",
     "Glutaminase + glycolysis + OXPHOS inhibitor",
     "CB-839 · IACS-010759 · 2-DG (tool) · enasidenib (IDH2)",
     "Phase I-II + approved (IDH2)"],
    ["TROP2 (non-SL)", "TACSTD2",
     "Antibody-drug conjugate (ADC)",
     "sacituzumab govitecan (Trodelvy) · datopotamab deruxtecan",
     "Approved (TNBC, urothelial)"],
]
table = ax.table(cellText=landscape, loc="center", cellLoc="left")
table.auto_set_font_size(False)
table.set_fontsize(9.5)
table.scale(1, 1.7)
for j in range(5):
    table[(0, j)].set_facecolor("#0F1A2E")
    table[(0, j)].set_text_props(color="white", weight="bold")
class_color_map = ["#7B1F2A", "#5E80B0", "#34547A", "#B8893C", "#3C6B4F",
                    "#3F7A8A", "#962E2E", "#8B5E00", "#5E5E5E"]
for i in range(1, 10):
    table[(i, 0)].set_facecolor(class_color_map[i - 1])
    table[(i, 0)].set_text_props(color="white", weight="bold")
ax.set_title("F16 — Paper 10 target class × drug landscape\n"
             "(approved / clinical-stage agents per class)", pad=12)
plt.tight_layout()
plt.savefig(FIG / "F16_drug_landscape.png", dpi=140, bbox_inches="tight")
plt.close()

# ====================================================================
# F17 — Deferred work registry
# ====================================================================
fig, ax = plt.subplots(figsize=(11, 5.8))
ax.axis("off")
deferred_rows = [
    ["#", "Step", "Why deferred", "Earliest"],
    ["D1", "DepMap CRISPR genome-wide dependency × DM1 state",
     "Disk budget 4.7 GB free at start; full matrix ~3 GB", "2026-06-16"],
    ["D2", "SL / SR pair inference (Step C of plan)",
     "Requires D1 first (DepMap dependency contrast)", "Post D1"],
    ["D3", "PRISM 19Q4 / 24Q2 IC50 stratification by DM1",
     "Raw cache on /opt/thyroid-dash, not on this dev box",
     "After Pod re-mount or re-pull"],
    ["D4", "GSE76039 per-target gene replication",
     "Existing DM1 score on disk; per-gene matrix not pre-extracted",
     "1 day work post-marathon"],
    ["D5", "GDSC / CTRP cross-validation",
     "Out of first-pass scope", "Post D1"],
    ["D6", "Wet-lab Tier 1 (siRNA / drug viability)",
     "Requires collaborator wet-lab access",
     "Conditional on Paper 1 in print"],
    ["D7", "PDX / organoid validation",
     "Tier 3 — beyond first paper", "Follow-up"],
    ["D8", "Korean K2 (PRJEB11591) replication",
     "Korean pivot is Paper 4 territory; SL angle deferred",
     "Post Paper 1 + 2 + 4"],
]
table = ax.table(cellText=deferred_rows, loc="center", cellLoc="left")
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 1.65)
for j in range(4):
    table[(0, j)].set_facecolor("#0F1A2E")
    table[(0, j)].set_text_props(color="white", weight="bold")
for i in range(1, 9):
    table[(i, 0)].set_facecolor("#A57215")
    table[(i, 0)].set_text_props(color="white", weight="bold")
ax.set_title("F17 — Deferred work registry (honest-scope statement)\n"
             "(none of these run before Paper 1 in print + marathon close)",
             pad=12)
plt.tight_layout()
plt.savefig(FIG / "F17_deferred_registry.png", dpi=140, bbox_inches="tight")
plt.close()

# ====================================================================
# Top-50 DGE (full table) ranking — F06
# ====================================================================
print("\nLoading full DGE for top-50 …")
full_dge = pd.read_csv(
    ROOT / "project/results/dark_matter_phase2/web/data/dge_dm1_vs_dm2.tsv",
    sep="\t"
)
full_dge["dm1_minus_dm2"] = -full_dge["log2FC_DM2vDM1"]
top50_up = full_dge.sort_values("dm1_minus_dm2", ascending=False).head(50)
top50_dn = full_dge.sort_values("dm1_minus_dm2").head(50)
fig, axes = plt.subplots(1, 2, figsize=(14, 11), sharey=False)
for ax, df, title, color in [
    (axes[0], top50_up.iloc[::-1], "Top-50 DM1-up genes (whole DGE)", "#7B1F2A"),
    (axes[1], top50_dn.iloc[::-1], "Top-50 DM2-up genes (whole DGE)", "#34547A"),
]:
    ax.barh(df["gene"], df["dm1_minus_dm2"], color=color, alpha=0.85)
    ax.axvline(0, color="black", linewidth=0.6)
    ax.set_xlabel("log2(DM1/DM2)")
    ax.set_title(title)
    ax.tick_params(axis="y", labelsize=8)
fig.suptitle("F06 — Top-50 DM1-up / DM2-up genes from full DGE (51,711 genes)\n"
             "(context for Paper 10 candidate panel)", y=1.0)
plt.tight_layout()
plt.savefig(FIG / "F06_top50_dge.png", dpi=140)
plt.close()

# ====================================================================
# F01 — Overview schematic (concept)
# ====================================================================
fig, ax = plt.subplots(figsize=(12, 6.5))
ax.axis("off")
ax.set_xlim(0, 12)
ax.set_ylim(0, 7)

ax.text(6, 6.3, "F01 — Paper 10 ATLAS concept",
        ha="center", fontsize=15, fontweight="bold", color="#0F1A2E")

# DM1 state
ax.add_patch(plt.Circle((2, 3.5), 1.2, edgecolor="#7B1F2A",
                          facecolor="#FAEEED", linewidth=2))
ax.text(2, 3.8, "DM1 state", ha="center", fontsize=11, fontweight="bold",
        color="#7B1F2A")
ax.text(2, 3.3, "lineage TF\ncollapse", ha="center", fontsize=9,
        color="#0F1A2E", style="italic")
ax.text(2, 2.8, "FOXE1 / NKX2-1\n/ PAX8 / HHEX ↓", ha="center",
        fontsize=8, color="#0F1A2E")

# Compensatory rewiring
ax.add_patch(plt.Circle((6, 3.5), 1.2, edgecolor="#B8893C",
                          facecolor="#FCF8EE", linewidth=2))
ax.text(6, 3.8, "Compensatory\nrewiring", ha="center", fontsize=11,
        fontweight="bold", color="#B8893C")
ax.text(6, 3.0, "STAT3 / AP-1 / DNMT\n+ SFK + replication\nstress + cytokine",
        ha="center", fontsize=8, color="#0F1A2E")

# Druggable nodes
ax.add_patch(plt.Circle((10, 3.5), 1.2, edgecolor="#3C6B4F",
                          facecolor="#E0EBDC", linewidth=2))
ax.text(10, 3.8, "Druggable\nSL / SR nodes", ha="center", fontsize=11,
        fontweight="bold", color="#3C6B4F")
ax.text(10, 3.0, "KCNN4 · OSMR\nIL6R · SFK · ATR\nMYC · NAMPT", ha="center",
        fontsize=8, color="#0F1A2E")

# Arrows
ax.annotate("", xy=(4.7, 3.5), xytext=(3.3, 3.5),
            arrowprops=dict(arrowstyle="->", color="black", lw=2))
ax.annotate("", xy=(8.7, 3.5), xytext=(7.3, 3.5),
            arrowprops=dict(arrowstyle="->", color="black", lw=2))
ax.text(4.0, 3.9, "loss of\nredundancy", ha="center", fontsize=9,
        style="italic", color="#3A4658")
ax.text(8.0, 3.9, "narrowed\nsurvival\ncorridor", ha="center", fontsize=9,
        style="italic", color="#3A4658")

# Bottom triangulation banner
ax.text(6, 1.0, "Triangulation across CCLE n=13 + DGE 51,711 genes + "
                "TCGA-THCA Cox n=560",
        ha="center", fontsize=11, color="#0F1A2E", weight="bold",
        bbox=dict(boxstyle="round,pad=0.5", edgecolor="#3F7A8A",
                  facecolor="#E5EEF0", linewidth=1.5))

plt.savefig(FIG / "F01_concept.png", dpi=140, bbox_inches="tight")
plt.close()


# ====================================================================
# Summary
# ====================================================================
all_figs = sorted(FIG.glob("*.png"))
print(f"\nGenerated {len(all_figs)} figures:")
for p in all_figs:
    print(f"  {p.stat().st_size/1024:6.1f} KB  {p.name}")

# Also dump a class-summary TSV for the HTML page
class_summary = []
for cls in sorted(set(c for c in CLASS_COLORS) - {"NAD_de_novo_SR_partner"}):
    s_dge = dge[dge["class"] == cls]
    s_ccle = ccle_corr[ccle_corr["class"] == cls].dropna(subset=["pearson_r"])
    s_cox = cox_full[cox_full["class"] == cls].dropna(subset=["p"])
    class_summary.append({
        "class": cls,
        "n_genes": len(s_dge),
        "dge_mean_log2_dm1_up": -s_dge["log2FC_DM2vDM1"].mean()
            if len(s_dge) else np.nan,
        "dge_mean_t": s_dge["t"].mean() if len(s_dge) else np.nan,
        "ccle_n_measured": len(s_ccle),
        "ccle_mean_r": s_ccle["pearson_r"].mean() if len(s_ccle) else np.nan,
        "cox_n": len(s_cox),
        "cox_mean_log_hr": np.log(s_cox["hr"]).mean()
            if len(s_cox) else np.nan,
        "cox_min_p": s_cox["p"].min() if len(s_cox) else np.nan,
    })
pd.DataFrame(class_summary).to_csv(OUT / "class_summary.tsv",
                                     sep="\t", index=False)
print(f"\nWrote class_summary.tsv")
print("Done.")
