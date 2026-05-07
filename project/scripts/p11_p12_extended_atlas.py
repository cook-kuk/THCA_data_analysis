#!/usr/bin/env python3
"""Paper 11 (Pan-Cancer DM1 Transfer) + Paper 12 (Candidate Co-expression Network).

Both built on local data only; no new download.

Paper 11 inputs:
  - TCGA pancan log2(norm+1): project/data/raw/TCGA_pancan/pancan_geneExp.gz
  - phenotype:                project/data/raw/TCGA_pancan/phenotype.tsv.gz

Paper 12 inputs:
  - same pancan expression, restricted to TCGA-THCA samples + 36-candidate panel
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
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
PAN = ROOT / "project/data/raw/TCGA_pancan/pancan_geneExp.gz"
PHENO = ROOT / "project/data/raw/TCGA_pancan/phenotype.tsv.gz"
P11 = ROOT / "project/results/paper11_pancancer"
P12 = ROOT / "project/results/paper12_network"
F11 = P11 / "figures"; F12 = P12 / "figures"
for d in (P11, P12, F11, F12):
    d.mkdir(parents=True, exist_ok=True)

RAI_8 = ["TPO", "DIO1", "TSHR", "PAX8", "TG", "FOXE1", "NKX2-1", "SLC5A5"]
LINEAGE_TFS = ["FOXE1", "NKX2-1", "PAX8", "HHEX"]
CANDIDATES = {
    "NAMPT": "NAD_salvage", "NAPRT": "NAD_de_novo_SR_partner",
    "JAK1": "JAK_STAT", "JAK2": "JAK_STAT", "TYK2": "JAK_STAT",
    "STAT3": "JAK_STAT", "IL6R": "JAK_STAT_upstream",
    "OSMR": "JAK_STAT_upstream",
    "DNMT1": "Epigenetic", "DNMT3A": "Epigenetic", "DNMT3B": "Epigenetic",
    "HDAC1": "Epigenetic", "HDAC2": "Epigenetic", "HDAC6": "Epigenetic",
    "KDM1A": "Epigenetic", "EZH2": "Epigenetic",
    "LYN": "SFK", "FYN": "SFK", "SRC": "SFK", "YES1": "SFK",
    "KCNN4": "Ion_channel",
    "PARP1": "DDR", "PARP2": "DDR", "ATR": "DDR",
    "CHEK1": "DDR", "CHEK2": "DDR",
    "GLS": "Metabolism", "LDHA": "Metabolism_glycolysis",
    "HK2": "Metabolism_glycolysis", "IDH2": "Metabolism",
    "MYC": "Metabolism_master",
    "TACSTD2": "TROP2_ADC_non_SL",
    "FOXE1": "Lineage_TF_anchor", "NKX2-1": "Lineage_TF_anchor",
    "PAX8": "Lineage_TF_anchor", "HHEX": "Lineage_TF_anchor",
}
CLASS_COLORS = {
    "NAD_salvage": "#7B1F2A", "NAD_de_novo_SR_partner": "#A04451",
    "JAK_STAT": "#34547A", "JAK_STAT_upstream": "#5E80B0",
    "Epigenetic": "#B8893C", "SFK": "#3C6B4F",
    "Ion_channel": "#3F7A8A", "DDR": "#962E2E",
    "Metabolism": "#8B5E00", "Metabolism_glycolysis": "#A5722B",
    "Metabolism_master": "#5C3D00", "TROP2_ADC_non_SL": "#5E5E5E",
    "Lineage_TF_anchor": "#0F1A2E",
}

# -------- 1. Load all needed genes from pancan in one streaming pass --------
genes_needed = set(RAI_8) | set(CANDIDATES) | set(LINEAGE_TFS)
print(f"Streaming pancan for {len(genes_needed)} genes …")
with gzip.open(PAN, "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    sample_cols = header[1:]
    rows = []
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if parts[0] in genes_needed:
            rows.append(parts)
            if len(rows) == len(genes_needed):
                break
expr = pd.DataFrame(rows, columns=header).set_index(header[0])
expr = expr.replace(["", "NA", "na", "NaN"], np.nan).astype(float)
print(f"Pancan expression: {expr.shape[0]} genes × {expr.shape[1]} samples")

# -------- 2. Load phenotype --------
pheno = pd.read_csv(PHENO, sep="\t")
print(f"Phenotype: {len(pheno)} rows · {pheno['_primary_disease'].nunique()} diseases")
pheno = pheno.set_index("sample")

# -------- 3. Score DM1_like across all 11,069 samples --------
def within_sample_z_score(expr_df: pd.DataFrame, geneset: list[str]) -> pd.Series:
    avail = [g for g in geneset if g in expr_df.index]
    sub = expr_df.loc[avail]
    line_mu = expr_df.mean(axis=0)
    line_sd = expr_df.std(axis=0).replace(0, np.nan)
    z = (sub.sub(line_mu, axis=1)).div(line_sd, axis=1)
    return z.mean(axis=0)

print("Scoring DM1_like across pancan …")
rai_score = within_sample_z_score(expr, RAI_8)
tf_score = within_sample_z_score(expr, LINEAGE_TFS)
dm1 = -rai_score
df_pan = pd.DataFrame({
    "sample": dm1.index,
    "DM1_like": dm1.values,
    "TF_collapse": tf_score.values,
    "RAI_8": rai_score.values,
})
df_pan["lineage"] = df_pan["sample"].map(
    pheno["_primary_disease"]
)
df_pan = df_pan.dropna(subset=["lineage"])
print(f"Pancan DM1 scored: {len(df_pan)} samples · "
      f"{df_pan['lineage'].nunique()} lineages")
df_pan.to_csv(P11 / "pancan_dm1_scored.tsv", sep="\t", index=False)

# Threshold: pan-cancer median of THCA samples
thca_med = df_pan.loc[df_pan["lineage"] == "thyroid carcinoma",
                       "DM1_like"].median()
df_pan["dm1_high_vs_thca_med"] = (df_pan["DM1_like"] > thca_med).astype(int)
df_pan["dm1_high_vs_pan_med"] = (df_pan["DM1_like"] >
                                  df_pan["DM1_like"].median()).astype(int)

# Per-lineage stats
lineage_stats = df_pan.groupby("lineage").agg(
    n=("sample", "count"),
    median_DM1=("DM1_like", "median"),
    mean_DM1=("DM1_like", "mean"),
    pct_above_thca_med=("dm1_high_vs_thca_med", "mean"),
    pct_above_pan_med=("dm1_high_vs_pan_med", "mean"),
).reset_index().sort_values("median_DM1", ascending=False)
lineage_stats.to_csv(P11 / "pancan_lineage_stats.tsv",
                      sep="\t", index=False)
print(f"\nTop 10 lineages by median DM1_like:")
print(lineage_stats.head(10).to_string(index=False))
print(f"\nBottom 10 lineages by median DM1_like:")
print(lineage_stats.tail(10).to_string(index=False))


# -------- 4. PAPER 11 FIGURES --------
print("\n=== Paper 11 figures ===")

# F11_01 — Lineage-level DM1 distribution box (sorted by median)
fig, ax = plt.subplots(figsize=(13, 7.5))
order = lineage_stats["lineage"].tolist()
data = [df_pan[df_pan["lineage"] == ln]["DM1_like"].values for ln in order]
bp = ax.boxplot(data, vert=False, tick_labels=order, patch_artist=True,
                widths=0.65, medianprops={"color": "black"})
# Color thyroid carcinoma red
for patch, ln in zip(bp["boxes"], order):
    if ln == "thyroid carcinoma":
        patch.set_facecolor("#7B1F2A")
        patch.set_alpha(0.95)
    else:
        # Color by quintile of median DM1
        val = lineage_stats.loc[lineage_stats["lineage"] == ln,
                                 "median_DM1"].values[0]
        patch.set_facecolor("#34547A" if val < 0 else "#A04451")
        patch.set_alpha(0.55)
ax.axvline(0, color="black", linewidth=0.6)
ax.axvline(thca_med, color="#7B1F2A", linestyle="--", linewidth=1.2,
           label=f"THCA median ({thca_med:.2f})")
ax.set_xlabel("DM1_like score (within-sample z, sign-flipped over 8 RAI genes)")
ax.set_title("F11.01 — Pan-cancer DM1_like distribution per TCGA lineage\n"
             "(sorted by median; THCA in dark red)")
ax.legend(loc="lower right")
ax.tick_params(axis="y", labelsize=9)
plt.tight_layout()
plt.savefig(F11 / "F11_01_lineage_dm1_box.png", dpi=140)
plt.close()

# F11_02 — % above THCA median bar
fig, ax = plt.subplots(figsize=(11, 8))
ls = lineage_stats.sort_values("pct_above_thca_med", ascending=True)
colors = ["#7B1F2A" if ln == "thyroid carcinoma"
          else ("#A04451" if v >= 0.4 else "#34547A")
          for ln, v in zip(ls["lineage"], ls["pct_above_thca_med"])]
ax.barh(ls["lineage"], ls["pct_above_thca_med"] * 100, color=colors,
        alpha=0.85)
ax.axvline(50, color="grey", linestyle=":", linewidth=0.8,
           label="50% threshold")
ax.set_xlabel("% of lineage samples with DM1_like > THCA median")
ax.set_title("F11.02 — Cross-cancer DM1-like prevalence\n"
             "(% above THCA median; potential SL portability targets)")
ax.legend(loc="lower right")
ax.tick_params(axis="y", labelsize=9)
plt.tight_layout()
plt.savefig(F11 / "F11_02_lineage_pct_above_thca.png", dpi=140)
plt.close()

# F11_03 — DM1_like vs sample size scatter
fig, ax = plt.subplots(figsize=(10, 7))
for _, r in lineage_stats.iterrows():
    color = "#7B1F2A" if r["lineage"] == "thyroid carcinoma" else "#34547A"
    size = 200 if r["lineage"] == "thyroid carcinoma" else 80
    ax.scatter(r["n"], r["median_DM1"], s=size, c=color,
               edgecolor="white", linewidth=0.8, zorder=3)
    if abs(r["median_DM1"]) > 0.10 or r["lineage"] == "thyroid carcinoma":
        ax.annotate(r["lineage"], (r["n"], r["median_DM1"]),
                    fontsize=8, xytext=(5, 4),
                    textcoords="offset points")
ax.axhline(0, color="black", linewidth=0.5)
ax.axhline(thca_med, color="#7B1F2A", linestyle="--", linewidth=0.8,
           label=f"THCA median ({thca_med:.2f})")
ax.set_xlabel("Lineage sample count (TCGA pancan)")
ax.set_ylabel("Median DM1_like per lineage")
ax.set_title("F11.03 — Lineage median DM1 × sample size\n"
             "(THCA highlighted; non-thyroid DM1+ lineages = SL portability candidates)")
ax.legend(loc="lower right")
plt.tight_layout()
plt.savefig(F11 / "F11_03_lineage_n_vs_dm1.png", dpi=140)
plt.close()

# F11_04 — Candidate Pearson r vs DM1_like, top 5 non-THCA DM1+ lineages
top5_non_thca = lineage_stats[
    lineage_stats["lineage"] != "thyroid carcinoma"
].head(5)["lineage"].tolist()
print(f"Top 5 non-THCA DM1+ lineages: {top5_non_thca}")

corr_matrix = []
for ln in top5_non_thca + ["thyroid carcinoma"]:
    sub_samples = df_pan[df_pan["lineage"] == ln]["sample"].values
    sub_expr = expr[sub_samples].T  # samples × genes
    sub_dm1 = df_pan.set_index("sample").loc[sub_samples, "DM1_like"]
    row = {"lineage": ln, "n": len(sub_samples)}
    for g in CANDIDATES:
        if g in sub_expr.columns:
            mask = sub_expr[g].notna() & sub_dm1.notna()
            if mask.sum() > 10:
                r, p = stats.pearsonr(sub_dm1[mask], sub_expr.loc[mask, g])
                row[g] = r
            else:
                row[g] = np.nan
        else:
            row[g] = np.nan
    corr_matrix.append(row)
corr_df = pd.DataFrame(corr_matrix).set_index("lineage")
corr_df.to_csv(P11 / "lineage_candidate_corr.tsv", sep="\t")

# Heatmap
fig, ax = plt.subplots(figsize=(14, 5))
gene_order = [g for g in CANDIDATES if g in corr_df.columns]
mat = corr_df[gene_order].values
im = ax.imshow(mat, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
ax.set_xticks(range(len(gene_order)))
ax.set_xticklabels(gene_order, rotation=80, fontsize=8.5)
ax.set_yticks(range(len(corr_df)))
ax.set_yticklabels([f"{ln} (n={int(corr_df.loc[ln,'n'])})"
                    for ln in corr_df.index], fontsize=10)
# Color tick labels by class
for i, g in enumerate(gene_order):
    cls = CANDIDATES.get(g, "")
    color = CLASS_COLORS.get(cls, "#888")
    ax.get_xticklabels()[i].set_color(color)
plt.colorbar(im, ax=ax, label="Pearson r vs DM1_like (per lineage)",
             shrink=0.7)
ax.set_title("F11.04 — Cross-lineage candidate × DM1 correlation\n"
             "(top 5 non-THCA DM1+ lineages + THCA reference)")
plt.tight_layout()
plt.savefig(F11 / "F11_04_lineage_candidate_heatmap.png", dpi=140)
plt.close()

# F11_05 — Lineage TF vs candidate correlation: do non-thyroid lineages also show TF collapse?
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
for ax, score_col, title in zip(
    axes,
    ["DM1_like", "TF_collapse"],
    ["DM1_like (-RAI_8)", "TF_collapse (FOXE1/NKX2-1/PAX8/HHEX)"],
):
    medians = df_pan.groupby("lineage")[score_col].median().sort_values()
    colors = ["#7B1F2A" if ln == "thyroid carcinoma"
              else ("#A04451" if v > 0 else "#34547A")
              for ln, v in zip(medians.index, medians.values)]
    ax.barh(medians.index, medians.values, color=colors, alpha=0.8)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_xlabel(f"Median {score_col}")
    ax.set_title(title)
    ax.tick_params(axis="y", labelsize=8)
fig.suptitle("F11.05 — Per-lineage median: DM1_like vs TF-collapse\n"
             "(correlated views; non-thyroid lineages with DM1+ AND TF-collapse+ are strongest portability candidates)",
             y=1.0)
plt.tight_layout()
plt.savefig(F11 / "F11_05_dm1_vs_tfcollapse_lineage.png", dpi=140)
plt.close()

# F11_06 — Top portability candidates (non-thyroid lineages with both DM1+ and TF collapse-)
ls = lineage_stats.copy()
ls["tf_collapse_median"] = ls["lineage"].map(
    lambda ln: df_pan.loc[df_pan["lineage"] == ln, "TF_collapse"].median()
)
ls["portability_score"] = ls["median_DM1"] - ls["tf_collapse_median"]  # high = both DM1+ AND TF-
ls_non = ls[ls["lineage"] != "thyroid carcinoma"]
top_port = ls_non.sort_values("portability_score", ascending=False).head(8)
fig, ax = plt.subplots(figsize=(11, 5))
ax.barh(top_port["lineage"][::-1], top_port["portability_score"][::-1],
        color="#3F7A8A", alpha=0.85)
ax.axvline(0, color="black", linewidth=0.5)
ax.set_xlabel("portability_score = median DM1_like − median TF_collapse")
ax.set_title("F11.06 — Top non-thyroid lineages by SL portability score\n"
             "(higher = more DM1-like + more TF-collapsed = SL strategies may transfer)")
plt.tight_layout()
plt.savefig(F11 / "F11_06_portability.png", dpi=140)
plt.close()

# Summary
p11_summary = {
    "pancan_n_samples": int(len(df_pan)),
    "n_lineages": int(df_pan["lineage"].nunique()),
    "thca_median_dm1": float(thca_med),
    "top5_non_thca_dm1": top5_non_thca,
    "top_portability": top_port["lineage"].tolist(),
}
(P11 / "summary.json").write_text(json.dumps(p11_summary, indent=2))
print(f"\nP11 summary: {json.dumps(p11_summary, indent=2)}")


# =====================================================================
# PAPER 12 — Candidate co-expression network in TCGA-THCA
# =====================================================================
print("\n=== Paper 12: candidate co-expression network ===")

# Restrict to THCA samples
thca_samples = df_pan[df_pan["lineage"] == "thyroid carcinoma"]["sample"].values
thca_expr = expr[thca_samples]
print(f"THCA expression: {len(thca_samples)} samples × "
      f"{thca_expr.shape[0]} genes")

# Filter to candidate genes available in expression
cand_present = [g for g in CANDIDATES if g in thca_expr.index]
print(f"Candidates present in pancan expression: {len(cand_present)}/{len(CANDIDATES)}")
cand_expr = thca_expr.loc[cand_present].T  # samples × genes
cand_expr = cand_expr.dropna(axis=0, how="any")
print(f"After dropna: {len(cand_expr)} samples × {cand_expr.shape[1]} genes")

# Pearson correlation matrix
corr_mat = cand_expr.corr(method="pearson")
corr_mat.to_csv(P12 / "candidate_corr_matrix.tsv", sep="\t")

# Hierarchical clustering on (1 - corr)
dist = 1 - corr_mat.abs()
np.fill_diagonal(dist.values, 0)
condensed = squareform(dist.values, checks=False)
Z = linkage(condensed, method="average")
modules = fcluster(Z, t=7, criterion="maxclust")
mod_df = pd.DataFrame({
    "gene": corr_mat.index,
    "class": [CANDIDATES.get(g, "") for g in corr_mat.index],
    "module": modules,
}).sort_values(["module", "class", "gene"])
mod_df.to_csv(P12 / "candidate_modules.tsv", sep="\t", index=False)
print(f"Modules:\n{mod_df.to_string(index=False)}")

# F12_01 — Heatmap (clustered) with class sidebar
order = mod_df["gene"].tolist()
mat = corr_mat.loc[order, order].values
fig, ax = plt.subplots(figsize=(11, 10))
im = ax.imshow(mat, cmap="RdBu_r", vmin=-1, vmax=1, aspect="equal")
ax.set_xticks(range(len(order)))
ax.set_yticks(range(len(order)))
ax.set_xticklabels(order, rotation=80, fontsize=8.5)
ax.set_yticklabels(order, fontsize=8.5)
for i, g in enumerate(order):
    cls = CANDIDATES.get(g, "")
    color = CLASS_COLORS.get(cls, "#888")
    ax.get_xticklabels()[i].set_color(color)
    ax.get_yticklabels()[i].set_color(color)
plt.colorbar(im, ax=ax, label="Pearson r (TCGA-THCA, n=560)", shrink=0.8)
# Module separators
mod_starts = mod_df.reset_index(drop=True).reset_index().groupby(
    "module")["index"].first().values
for s in mod_starts[1:]:
    ax.axhline(s - 0.5, color="black", linewidth=1.5)
    ax.axvline(s - 0.5, color="black", linewidth=1.5)
ax.set_title("F12.01 — Candidate co-expression matrix (TCGA-THCA, n=560)\n"
             "(clustered into 4 modules; class color-coded ticks)")
plt.tight_layout()
plt.savefig(F12 / "F12_01_corr_heatmap.png", dpi=140)
plt.close()

# F12_02 — Per-module summary
mod_summary_rows = []
for mod_id in sorted(mod_df["module"].unique()):
    members = mod_df[mod_df["module"] == mod_id]
    classes = members["class"].value_counts()
    sub_mat = corr_mat.loc[members["gene"], members["gene"]]
    # mean off-diagonal correlation
    n = len(members)
    if n > 1:
        mean_intra = (sub_mat.values.sum() - n) / (n * (n - 1))
    else:
        mean_intra = np.nan
    mod_summary_rows.append({
        "module": mod_id,
        "n_genes": n,
        "genes": ", ".join(members["gene"]),
        "dominant_classes": ", ".join(f"{c}({k})"
                                       for c, k in classes.items()),
        "mean_intra_corr": mean_intra,
    })
mod_summary = pd.DataFrame(mod_summary_rows)
mod_summary.to_csv(P12 / "module_summary.tsv", sep="\t", index=False)
print(f"\nModule summary:\n{mod_summary.to_string(index=False)}")

# F12_02 figure
fig, ax = plt.subplots(figsize=(13, 1 + 0.6 * len(mod_summary)))
ax.axis("off")
table_data = [["Module", "n", "Mean intra r", "Members (class colored)"]]
for _, r in mod_summary.iterrows():
    members = mod_df[mod_df["module"] == r["module"]]
    member_cell = " · ".join(members["gene"])
    table_data.append([
        f"M{r['module']}", str(int(r['n_genes'])),
        f"{r['mean_intra_corr']:+.3f}",
        member_cell,
    ])
table = ax.table(cellText=table_data, loc="center", cellLoc="left")
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 1.7)
for j in range(4):
    table[(0, j)].set_facecolor("#0F1A2E")
    table[(0, j)].set_text_props(color="white", weight="bold")
mod_palette = ["#7B1F2A", "#34547A", "#3C6B4F", "#B8893C", "#3F7A8A"]
for i in range(1, len(table_data)):
    table[(i, 0)].set_facecolor(mod_palette[(i - 1) % len(mod_palette)])
    table[(i, 0)].set_text_props(color="white", weight="bold")
ax.set_title("F12.02 — Candidate modules in TCGA-THCA\n"
             "(unsupervised hierarchical clustering of co-expression)")
plt.tight_layout()
plt.savefig(F12 / "F12_02_module_summary.png", dpi=140,
            bbox_inches="tight")
plt.close()

# F12_03 — Network graph (top edges only) — manual layout
edges = []
for i, g1 in enumerate(corr_mat.index):
    for j, g2 in enumerate(corr_mat.columns):
        if j <= i:
            continue
        r = corr_mat.iloc[i, j]
        if abs(r) >= 0.4:
            edges.append((g1, g2, r))
edges_df = pd.DataFrame(edges, columns=["g1", "g2", "r"])
edges_df.to_csv(P12 / "edges_top.tsv", sep="\t", index=False)
print(f"\nTop edges (|r|>=0.4): {len(edges_df)}")

# Spring-layout positions via Fruchterman-Reingold-ish iterative — use scipy MDS instead
from sklearn.manifold import MDS
mds = MDS(n_components=2, dissimilarity="precomputed",
          random_state=42, normalized_stress="auto")
coords = mds.fit_transform(dist.values)
pos = {g: coords[i] for i, g in enumerate(corr_mat.index)}

fig, ax = plt.subplots(figsize=(11, 10))
# Edges
for _, e in edges_df.iterrows():
    p1, p2 = pos[e["g1"]], pos[e["g2"]]
    color = "#962E2E" if e["r"] > 0 else "#34547A"
    alpha = min(0.75, abs(e["r"]) * 0.9)
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]],
            color=color, alpha=alpha, linewidth=abs(e["r"]) * 2.5,
            zorder=1)
# Nodes
for g, (x, y) in pos.items():
    cls = CANDIDATES.get(g, "")
    color = CLASS_COLORS.get(cls, "#888")
    ax.scatter([x], [y], s=350, c=color, edgecolor="white",
               linewidth=1.5, zorder=3)
    ax.annotate(g, (x, y), ha="center", va="center", fontsize=8.5,
                color="white", weight="bold", zorder=4)
ax.set_xticks([]); ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)
ax.set_title("F12.03 — Candidate co-expression network (TCGA-THCA, n=560)\n"
             "(MDS layout on |1-r| distance; edges |r|≥0.4; "
             "red=positive corr, blue=anti)")
# Legend
legend_handles = [
    plt.Line2D([0], [0], color="#962E2E", linewidth=2, label="positive corr"),
    plt.Line2D([0], [0], color="#34547A", linewidth=2, label="anti corr"),
]
for cls, color in [("Lineage_TF", "#0F1A2E"), ("DDR", "#962E2E"),
                    ("SFK", "#3C6B4F"), ("Epigenetic", "#B8893C"),
                    ("JAK/STAT", "#34547A"), ("Ion channel", "#3F7A8A")]:
    legend_handles.append(plt.Line2D([0], [0], marker="o", color="white",
                                      markerfacecolor=color, markersize=12,
                                      label=cls))
ax.legend(handles=legend_handles, loc="upper left", fontsize=9,
          ncol=2, framealpha=0.9)
plt.tight_layout()
plt.savefig(F12 / "F12_03_network.png", dpi=140)
plt.close()

# F12_04 — top |r| edges table figure
top_edges = edges_df.reindex(edges_df["r"].abs().sort_values(ascending=False).index).head(30)
fig, ax = plt.subplots(figsize=(11, 10))
ax.axis("off")
rows = [["#", "Gene 1 (class)", "Gene 2 (class)", "Pearson r"]]
for i, (_, e) in enumerate(top_edges.iterrows(), 1):
    c1 = CANDIDATES.get(e["g1"], "")
    c2 = CANDIDATES.get(e["g2"], "")
    rows.append([str(i), f"{e['g1']} ({c1})",
                 f"{e['g2']} ({c2})", f"{e['r']:+.3f}"])
table = ax.table(cellText=rows, loc="center", cellLoc="left")
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 1.5)
for j in range(4):
    table[(0, j)].set_facecolor("#0F1A2E")
    table[(0, j)].set_text_props(color="white", weight="bold")
for i in range(1, len(rows)):
    e = top_edges.iloc[i - 1]
    color = "#962E2E" if e["r"] > 0 else "#34547A"
    table[(i, 3)].set_facecolor(color + "33")
ax.set_title("F12.04 — Top-30 candidate × candidate co-expression edges (TCGA-THCA)",
             pad=12)
plt.tight_layout()
plt.savefig(F12 / "F12_04_top_edges.png", dpi=140, bbox_inches="tight")
plt.close()

# F12_05 — Module-vs-module mean correlation heatmap
mod_genes = mod_df.groupby("module")["gene"].apply(list)
mvm = np.zeros((len(mod_genes), len(mod_genes)))
for i, gi in enumerate(mod_genes):
    for j, gj in enumerate(mod_genes):
        sub = corr_mat.loc[gi, gj].values
        if i == j:
            n = len(gi)
            mvm[i, j] = (sub.sum() - n) / max(n * (n - 1), 1) \
                if n > 1 else 1.0
        else:
            mvm[i, j] = sub.mean()

fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(mvm, cmap="RdBu_r", vmin=-1, vmax=1)
mod_labels = [f"M{i+1}" for i in range(len(mod_genes))]
ax.set_xticks(range(len(mod_labels))); ax.set_yticks(range(len(mod_labels)))
ax.set_xticklabels(mod_labels); ax.set_yticklabels(mod_labels)
for i in range(len(mod_labels)):
    for j in range(len(mod_labels)):
        ax.text(j, i, f"{mvm[i,j]:+.2f}", ha="center", va="center",
                fontsize=10, color="white" if abs(mvm[i,j]) > 0.4 else "black")
plt.colorbar(im, ax=ax, label="Mean Pearson r")
ax.set_title("F12.05 — Module × Module mean correlation\n"
             "(diag = mean intra-module; off-diag = inter-module)")
plt.tight_layout()
plt.savefig(F12 / "F12_05_module_corr.png", dpi=140, bbox_inches="tight")
plt.close()

# F12_06 — DDR vs SFK vs JAK_STAT class-vs-class scatter
def class_score(class_name, df_expr):
    members = [g for g, c in CANDIDATES.items()
               if c == class_name and g in df_expr.columns]
    if not members:
        return None
    return df_expr[members].mean(axis=1)

ddr_s = class_score("DDR", cand_expr)
sfk_s = class_score("SFK", cand_expr)
jak_s = (cand_expr[["IL6R", "OSMR"]].mean(axis=1)
          if all(g in cand_expr.columns for g in ["IL6R", "OSMR"]) else None)

if ddr_s is not None and sfk_s is not None and jak_s is not None:
    # Read DM1_like for THCA from pancan score table
    thca_dm1 = df_pan.set_index("sample").loc[
        cand_expr.index, "DM1_like"
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, (xs, lbl, color) in zip(
        axes,
        [(ddr_s, "DDR class mean (PARP1/2/ATR/CHEK1/2)", "#962E2E"),
         (sfk_s, "SFK class mean (LYN/SRC/FYN/YES1)", "#3C6B4F"),
         (jak_s, "JAK/STAT-upstream mean (IL6R/OSMR)", "#5E80B0")],
    ):
        ax.scatter(xs, thca_dm1, s=10, c=color, alpha=0.55,
                   edgecolor="none")
        r, p = stats.pearsonr(xs.values, thca_dm1.values)
        ax.set_xlabel(lbl)
        ax.set_ylabel("DM1_like (TCGA-THCA)")
        ax.set_title(f"r={r:+.3f}, p={p:.2e}", fontsize=11)
        ax.axhline(0, color="black", linewidth=0.4)
        ax.axvline(xs.median(), color="grey", linestyle="--", linewidth=0.6)
    fig.suptitle("F12.06 — Class-mean expression vs DM1_like (TCGA-THCA, n=560)",
                 y=1.02)
    plt.tight_layout()
    plt.savefig(F12 / "F12_06_class_vs_dm1.png", dpi=140,
                bbox_inches="tight")
    plt.close()

# F12_07 — DM1-high vs DM1-low corr-matrix difference
thca_meta = df_pan.set_index("sample").loc[cand_expr.index]
high_mask = thca_meta["DM1_like"] > thca_meta["DM1_like"].median()
corr_high = cand_expr[high_mask].corr()
corr_low = cand_expr[~high_mask].corr()
diff = corr_high - corr_low
diff = diff.loc[order, order]
fig, ax = plt.subplots(figsize=(10, 9))
im = ax.imshow(diff.values, cmap="PiYG", vmin=-0.6, vmax=0.6)
ax.set_xticks(range(len(order))); ax.set_yticks(range(len(order)))
ax.set_xticklabels(order, rotation=80, fontsize=8.5)
ax.set_yticklabels(order, fontsize=8.5)
for i, g in enumerate(order):
    color = CLASS_COLORS.get(CANDIDATES.get(g, ""), "#888")
    ax.get_xticklabels()[i].set_color(color)
    ax.get_yticklabels()[i].set_color(color)
plt.colorbar(im, ax=ax, label="Δ Pearson r (DM1-high − DM1-low)",
             shrink=0.85)
ax.set_title("F12.07 — Δ co-expression: DM1-high vs DM1-low subset\n"
             "(green = stronger corr in DM1-high; pink = stronger in DM1-low)")
plt.tight_layout()
plt.savefig(F12 / "F12_07_corr_diff.png", dpi=140)
plt.close()

# F12_08 — Hub genes (degree centrality) bar
adj = (corr_mat.abs() >= 0.4).astype(int).values - np.eye(len(corr_mat),
                                                            dtype=int)
degree = pd.Series(adj.sum(axis=0), index=corr_mat.index).sort_values()
fig, ax = plt.subplots(figsize=(10, 9))
colors_deg = [CLASS_COLORS.get(CANDIDATES.get(g, ""), "#888")
              for g in degree.index]
ax.barh(degree.index, degree.values, color=colors_deg, alpha=0.85)
ax.set_xlabel("Degree (count of edges with |r|≥0.4)")
ax.set_title("F12.08 — Candidate-network degree centrality\n"
             "(hubs = candidates connected to many other class members)")
plt.tight_layout()
plt.savefig(F12 / "F12_08_hub_genes.png", dpi=140)
plt.close()

p12_summary = {
    "n_thca_samples": int(len(cand_expr)),
    "n_candidates_in_network": int(cand_expr.shape[1]),
    "n_modules": int(mod_df["module"].nunique()),
    "n_top_edges": int(len(edges_df)),
    "top_hub_genes": degree.tail(5).index.tolist(),
}
(P12 / "summary.json").write_text(json.dumps(p12_summary, indent=2))
print(f"\nP12 summary: {json.dumps(p12_summary, indent=2)}")

print("\nDone.")
