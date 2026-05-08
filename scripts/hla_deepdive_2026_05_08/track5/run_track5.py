#!/usr/bin/env python3
"""
Track 5 — DM1 x HLA-I gene-expression module in TCGA-THCA.

Boundary: HLA-I/II are analyzed here strictly as transcriptomic gene-expression modules,
NOT as allele genotypes. Cancer-cohort allele genotyping is out of scope per
project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md (Sections 1.1, 1.2).
"""
from __future__ import annotations

import gzip
import json
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/hla_deepdive_2026_05_08/track5_dm1_hla1_module"
FIG = OUT / "figs"
TAB = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

CAPTION_BOILERPLATE = (
    "HLA-I/II gene-expression module - not allele genotype. "
    "Cancer-cohort allele genotyping is out of scope per separation rules."
)

# ============================================================
# 1. Inputs
# ============================================================
PANCAN_EXP = ROOT / "project/data/raw/TCGA_pancan/pancan_geneExp.gz"
PANCAN_PHENO = ROOT / "project/data/raw/TCGA_pancan/phenotype.tsv.gz"
PANCAN_DM1 = ROOT / "project/results/paper11_pancancer/pancan_dm1_scored.tsv"
DM_MASTER = ROOT / "project/results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv"

HLA1_GENES = [
    "HLA-A", "HLA-B", "HLA-C",
    "B2M", "TAP1", "TAP2", "TAPBP", "NLRC5", "IRF1",
    "PSMB8", "PSMB9", "ERAP1", "ERAP2",
    "HLA-E", "HLA-F", "HLA-G",
    "CALR", "CANX", "PDIA3",
]
HLA2_GENES = [
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
    "HLA-DQA1", "HLA-DQB1", "HLA-DMA", "HLA-DMB",
    "HLA-DOA", "HLA-DOB", "CIITA", "CD74",
]
PANEL8 = ["TG", "TPO", "TSHR", "SLC5A5", "FOXE1", "PAX8", "NKX2-1", "DIO1"]

# Purity proxies (transcriptomic). We use ESTIMATE-style ImmuneScore and StromalScore
# computed as z-scored mean of well-known marker subsets, as an in-script proxy when
# no precomputed CPE/ABSOLUTE/ESTIMATE table exists in the repo. Caveat documented.
IMMUNE_PROXY_GENES = [
    "CD2", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "CD19", "MS4A1",
    "CD79A", "CD79B", "PRF1", "GZMB", "GZMA", "NKG7", "KLRB1", "ITGAX",
    "ITGAM", "CD68", "CD163", "FCGR3A", "PTPRC",
]
# stromal / fibroblast / endothelial
STROMAL_PROXY_GENES = [
    "COL1A1", "COL1A2", "COL3A1", "COL5A1", "FAP", "ACTA2", "PDGFRA",
    "PDGFRB", "DCN", "LUM", "VIM", "PECAM1", "VWF", "CDH5", "CD34",
]
# tumor purity proxy: 1 - (immune+stromal). For thyroid we ALSO have thyrocyte
# differentiation markers (PANEL8) but these are confounded with DM1 itself.

print("[1/12] loading DM1 master + scoring scaffolds")
dm_master = pd.read_csv(DM_MASTER, sep="\t")
# tcga_short e.g. TCGA-4C-A93U; expression columns are TCGA-4C-A93U-01
dm_master["sample"] = dm_master["tcga_short"].astype(str) + "-01"

dm1_pancan = pd.read_csv(PANCAN_DM1, sep="\t")
dm1_thca = dm1_pancan[dm1_pancan["lineage"] == "thyroid carcinoma"].copy()
dm1_thca = dm1_thca.rename(columns={"DM1_like": "DM1_score"})
print(f"  pancan DM1-scored thyroid samples: {len(dm1_thca)}")

# Phenotype to identify Primary Tumor only
pheno = pd.read_csv(PANCAN_PHENO, sep="\t")
thy_primary = pheno[(pheno["_primary_disease"] == "thyroid carcinoma") &
                    (pheno["sample_type"] == "Primary Tumor")]["sample"].tolist()
print(f"  TCGA-THCA Primary Tumor samples: {len(thy_primary)}")

# ============================================================
# 2. Load expression for THCA primary-tumor samples and required genes
# ============================================================
print("[2/12] loading TCGA-THCA expression (gene x sample) for HLA-I/II + 8-panel + purity proxies")
all_genes = set(HLA1_GENES + HLA2_GENES + PANEL8 + IMMUNE_PROXY_GENES + STROMAL_PROXY_GENES)

# Streaming read of gzipped TSV
with gzip.open(PANCAN_EXP, "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    sample_cols = header[1:]
    keep_idx = [i for i, s in enumerate(sample_cols) if s in set(thy_primary)]
    keep_samples = [sample_cols[i] for i in keep_idx]
    rows = {}
    for ln in fh:
        gene, *vals = ln.rstrip("\n").split("\t")
        if gene in all_genes:
            arr = np.array([vals[i] for i in keep_idx], dtype=float)
            rows[gene] = arr

expr = pd.DataFrame(rows, index=keep_samples).T  # genes x samples
expr = expr.T  # samples x genes
print(f"  expression matrix: {expr.shape[0]} samples x {expr.shape[1]} genes")
missing = [g for g in HLA1_GENES if g not in expr.columns]
print(f"  missing HLA-I genes: {missing}")
missing2 = [g for g in HLA2_GENES if g not in expr.columns]
print(f"  missing HLA-II genes: {missing2}")

# ============================================================
# 3. Compute module scores (z-scored mean across module genes)
# ============================================================
print("[3/12] computing HLA-I, HLA-II, immune-proxy, stromal-proxy module scores")


def zmean_module(df, gene_list):
    sub = df.reindex(columns=[g for g in gene_list if g in df.columns])
    if sub.shape[1] == 0:
        return pd.Series(np.nan, index=df.index)
    z = (sub - sub.mean()) / sub.std(ddof=0)
    return z.mean(axis=1)


hla1_score = zmean_module(expr, HLA1_GENES)
hla2_score = zmean_module(expr, HLA2_GENES)
immune_score = zmean_module(expr, IMMUNE_PROXY_GENES)
stromal_score = zmean_module(expr, STROMAL_PROXY_GENES)
# transcriptomic purity proxy: -(immune + stromal)
purity_proxy = -(immune_score + stromal_score)
purity_proxy = (purity_proxy - purity_proxy.mean()) / purity_proxy.std(ddof=0)

# Compute DM1 score from the 8-gene panel (sign: +DIO1/+TG/+TPO/+TSHR/+SLC5A5/+FOXE1/+PAX8/+NKX2-1
# anti-thyroid_diff = -mean of these 8). DM1_high = low thyroid differentiation.
panel_z = zmean_module(expr, PANEL8)
dm1_local = -panel_z  # high = de-differentiated
dm1_local.name = "DM1_local"

# Merge with pancan DM1 (already z-scored vs all lineages) for sanity
df = pd.DataFrame({
    "sample": expr.index,
    "HLA1_score": hla1_score.values,
    "HLA2_score": hla2_score.values,
    "Immune_proxy": immune_score.values,
    "Stromal_proxy": stromal_score.values,
    "Purity_proxy": purity_proxy.values,
    "DM1_local": dm1_local.values,
})
# join pancan DM1
df = df.merge(dm1_thca[["sample", "DM1_score"]], on="sample", how="left")
# join master driver labels
df = df.merge(
    dm_master[[
        "sample", "driver_anchor_v17", "dm_status", "tert_pos",
        "stage", "OS", "OS.time", "PFI", "PFI.time"
    ]],
    on="sample", how="left",
)
# fallback: use local DM1 where pancan DM1 missing
df["DM1_use"] = df["DM1_score"].fillna(df["DM1_local"])
print(f"  merged frame: {df.shape}")
print(f"  pancan DM1 vs local DM1 Spearman ρ = "
      f"{stats.spearmanr(df['DM1_score'].dropna(), df.loc[df['DM1_score'].notna(),'DM1_local']).statistic:.3f}")

# Driver simplification
def simplify_driver(row):
    d = str(row.get("driver_anchor_v17", ""))
    if d == "BRAF":
        return "BRAF"
    if d == "RAS":
        return "RAS"
    if d in ("Fusion", "fusion", "FUSION"):
        return "Fusion"
    if pd.isna(row.get("driver_anchor_v17")) or row.get("driver_anchor_v17") in (None, "", "unknown"):
        return "TripleNeg"
    return "Other"


df["driver_simple"] = df.apply(simplify_driver, axis=1)
# DM1 tertiles
df["DM1_tertile"] = pd.qcut(df["DM1_use"], 3, labels=["Low", "Mid", "High"])

df.to_csv(TAB / "T01_per_sample_module_scores.tsv", sep="\t", index=False)
print("  wrote T01_per_sample_module_scores.tsv")

# ============================================================
# 4. Headline DM1 x HLA-I and DM1 x HLA-II
# ============================================================
print("[4/12] DM1 x HLA-I/II correlations + LOWESS")


def cohens_d(x, y):
    nx, ny = len(x), len(y)
    sx, sy = x.std(ddof=1), y.std(ddof=1)
    sp = np.sqrt(((nx - 1) * sx ** 2 + (ny - 1) * sy ** 2) / (nx + ny - 2))
    return (x.mean() - y.mean()) / sp if sp > 0 else np.nan


def headline_correlation(df, x, y, label):
    sub = df[[x, y]].dropna()
    rho_s, p_s = stats.spearmanr(sub[x], sub[y])
    rho_p, p_p = stats.pearsonr(sub[x], sub[y])
    return {
        "comparison": label, "n": len(sub),
        "spearman_rho": rho_s, "spearman_p": p_s,
        "pearson_r": rho_p, "pearson_p": p_p,
    }


def tertile_test(df, mod_col, label):
    g = df.dropna(subset=[mod_col, "DM1_tertile"])
    hi = g[g["DM1_tertile"] == "High"][mod_col]
    lo = g[g["DM1_tertile"] == "Low"][mod_col]
    md = g[g["DM1_tertile"] == "Mid"][mod_col]
    w = stats.mannwhitneyu(hi, lo, alternative="two-sided")
    d = cohens_d(hi, lo)
    return {
        "comparison": label, "n_high": len(hi), "n_low": len(lo), "n_mid": len(md),
        "wilcoxon_p": w.pvalue, "cohens_d_high_vs_low": d,
        "median_high": hi.median(), "median_mid": md.median(), "median_low": lo.median(),
    }


corr_rows = []
for mod, lbl in [("HLA1_score", "DM1 vs HLA-I module"),
                 ("HLA2_score", "DM1 vs HLA-II module"),
                 ("Immune_proxy", "DM1 vs Immune proxy"),
                 ("Stromal_proxy", "DM1 vs Stromal proxy"),
                 ("Purity_proxy", "DM1 vs Purity proxy")]:
    corr_rows.append(headline_correlation(df, "DM1_use", mod, lbl))
corr_tab = pd.DataFrame(corr_rows)
corr_tab.to_csv(TAB / "T02_headline_correlations.tsv", sep="\t", index=False)
print(corr_tab.to_string(index=False))

tert_rows = []
for mod, lbl in [("HLA1_score", "HLA-I across DM1 tertiles"),
                 ("HLA2_score", "HLA-II across DM1 tertiles"),
                 ("Immune_proxy", "Immune across DM1 tertiles"),
                 ("Stromal_proxy", "Stromal across DM1 tertiles")]:
    tert_rows.append(tertile_test(df, mod, lbl))
tert_tab = pd.DataFrame(tert_rows)
tert_tab.to_csv(TAB / "T03_dm1_tertile_tests.tsv", sep="\t", index=False)
print(tert_tab.to_string(index=False))

# Figures: scatter + LOWESS
try:
    from statsmodels.nonparametric.smoothers_lowess import lowess
    HAVE_LOWESS = True
except Exception:
    HAVE_LOWESS = False


def scatter_with_lowess(x, y, ax, title, xlabel="DM1 score", ylabel="Module score"):
    ax.scatter(x, y, s=8, alpha=0.5, color="#2b8cbe", edgecolor="none")
    if HAVE_LOWESS and len(x) > 30:
        sm = lowess(y, x, frac=0.4, return_sorted=True)
        ax.plot(sm[:, 0], sm[:, 1], color="#e34a33", lw=2, label="LOWESS")
    rho, p = stats.spearmanr(x, y)
    ax.set_title(f"{title}\nSpearman rho={rho:.3f}, p={p:.2e}, n={len(x)}", fontsize=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.3)


fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
sub = df[["DM1_use", "HLA1_score"]].dropna()
scatter_with_lowess(sub["DM1_use"].values, sub["HLA1_score"].values, axes[0],
                    "TCGA-THCA: DM1 vs HLA-I module", ylabel="HLA-I module score")
sub2 = df[["DM1_use", "HLA2_score"]].dropna()
scatter_with_lowess(sub2["DM1_use"].values, sub2["HLA2_score"].values, axes[1],
                    "TCGA-THCA: DM1 vs HLA-II module", ylabel="HLA-II module score")
fig.suptitle(CAPTION_BOILERPLATE, fontsize=8, color="#666")
fig.tight_layout()
fig.savefig(FIG / "F01_dm1_vs_hla1_hla2_scatter.png", dpi=180)
fig.savefig(FIG / "F01_dm1_vs_hla1_hla2_scatter.pdf")
plt.close(fig)

# Tertile boxplots
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for ax, mod, lbl in zip(axes, ["HLA1_score", "HLA2_score"], ["HLA-I module", "HLA-II module"]):
    data = [df[df["DM1_tertile"] == t][mod].dropna() for t in ["Low", "Mid", "High"]]
    bp = ax.boxplot(data, labels=["Low", "Mid", "High"], patch_artist=True,
                    boxprops=dict(facecolor="#a6cee3"))
    res = tertile_test(df, mod, lbl)
    ax.set_title(f"{lbl} across DM1 tertiles\np={res['wilcoxon_p']:.2e}, d={res['cohens_d_high_vs_low']:.2f}")
    ax.set_ylabel(lbl)
    ax.set_xlabel("DM1 tertile")
    ax.grid(alpha=0.3, axis="y")
fig.suptitle(CAPTION_BOILERPLATE, fontsize=8, color="#666")
fig.tight_layout()
fig.savefig(FIG / "F02_dm1_tertile_boxplots.png", dpi=180)
fig.savefig(FIG / "F02_dm1_tertile_boxplots.pdf")
plt.close(fig)

# ============================================================
# 5. Per-gene heatmap, sorted by DM1
# ============================================================
print("[5/12] per-gene HLA-I heatmap sorted by DM1")
df_sorted = df.sort_values("DM1_use", na_position="last").reset_index(drop=True)
hla1_present = [g for g in HLA1_GENES if g in expr.columns]
hla1_mat = expr.loc[df_sorted["sample"], hla1_present].copy()
# z-score per gene
hla1_mat_z = (hla1_mat - hla1_mat.mean()) / hla1_mat.std(ddof=0)

fig = plt.figure(figsize=(13, 6.5))
gs = fig.add_gridspec(3, 1, height_ratios=[0.06, 0.06, 1.0], hspace=0.06)
ax_strip1 = fig.add_subplot(gs[0])
ax_strip2 = fig.add_subplot(gs[1])
ax_heat = fig.add_subplot(gs[2])

# DM1 strip
dm_vals = df_sorted["DM1_use"].values.reshape(1, -1)
ax_strip1.imshow(dm_vals, aspect="auto", cmap="RdBu_r",
                 vmin=np.nanpercentile(dm_vals, 5), vmax=np.nanpercentile(dm_vals, 95))
ax_strip1.set_yticks([0]); ax_strip1.set_yticklabels(["DM1"])
ax_strip1.set_xticks([])

# Driver strip
driver_map = {"BRAF": 0, "RAS": 1, "Fusion": 2, "TripleNeg": 3, "Other": 4}
driver_strip = df_sorted["driver_simple"].map(driver_map).fillna(4).astype(int).values.reshape(1, -1)
ax_strip2.imshow(driver_strip, aspect="auto", cmap="Set3", vmin=0, vmax=5)
ax_strip2.set_yticks([0]); ax_strip2.set_yticklabels(["Driver"])
ax_strip2.set_xticks([])

# Heatmap
im = ax_heat.imshow(hla1_mat_z.T.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
ax_heat.set_yticks(range(len(hla1_present)))
ax_heat.set_yticklabels(hla1_present, fontsize=8)
ax_heat.set_xticks([])
ax_heat.set_xlabel(f"TCGA-THCA samples sorted by DM1 score (low -> high), n={len(df_sorted)}")
fig.colorbar(im, ax=ax_heat, fraction=0.012, pad=0.01)
ax_heat.set_title(f"HLA-I module per-gene expression sorted by DM1\n{CAPTION_BOILERPLATE}",
                  fontsize=9)
fig.tight_layout()
fig.savefig(FIG / "F03_hla1_per_gene_heatmap.png", dpi=180)
fig.savefig(FIG / "F03_hla1_per_gene_heatmap.pdf")
plt.close(fig)

# Per-gene Spearman with DM1
per_gene_rows = []
for g in HLA1_GENES + HLA2_GENES:
    if g not in expr.columns:
        per_gene_rows.append({"gene": g, "module": "HLA-I" if g in HLA1_GENES else "HLA-II",
                              "spearman_rho": np.nan, "spearman_p": np.nan, "n": 0})
        continue
    vals = expr[g].reindex(df["sample"].values).values
    mask = ~np.isnan(df["DM1_use"].values) & ~np.isnan(vals)
    rho, p = stats.spearmanr(df["DM1_use"].values[mask], vals[mask])
    per_gene_rows.append({"gene": g, "module": "HLA-I" if g in HLA1_GENES else "HLA-II",
                          "spearman_rho": rho, "spearman_p": p, "n": int(mask.sum())})
per_gene_df = pd.DataFrame(per_gene_rows).sort_values("spearman_rho")
per_gene_df.to_csv(TAB / "T04_per_gene_dm1_corr.tsv", sep="\t", index=False)

# Lollipop plot of per-gene rho
fig, ax = plt.subplots(figsize=(10, 7))
plot_df = per_gene_df.dropna().copy()
colors = ["#1b7837" if m == "HLA-I" else "#762a83" for m in plot_df["module"]]
ax.hlines(plot_df["gene"], 0, plot_df["spearman_rho"], colors=colors, alpha=0.6)
ax.scatter(plot_df["spearman_rho"], plot_df["gene"], c=colors, s=40)
ax.axvline(0, color="k", lw=0.6)
ax.set_xlabel("Spearman rho with DM1 score")
ax.set_title(f"Per-gene DM1 correlation: HLA-I (green) and HLA-II (purple)\n{CAPTION_BOILERPLATE}",
             fontsize=9)
ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
fig.savefig(FIG / "F04_per_gene_dm1_corr_lollipop.png", dpi=180)
fig.savefig(FIG / "F04_per_gene_dm1_corr_lollipop.pdf")
plt.close(fig)

# ============================================================
# 6. Driver-stratified analysis
# ============================================================
print("[6/12] driver-stratified Spearman + forest")
driver_rows = []
for driver in ["BRAF", "RAS", "Fusion", "TripleNeg"]:
    sub = df[df["driver_simple"] == driver].dropna(subset=["DM1_use", "HLA1_score"])
    if len(sub) < 5:
        continue
    rho1, p1 = stats.spearmanr(sub["DM1_use"], sub["HLA1_score"])
    rho2, p2 = stats.spearmanr(sub["DM1_use"], sub["HLA2_score"])
    # Fisher z 95% CI
    def fisher_ci(r, n):
        if n < 4 or pd.isna(r):
            return (np.nan, np.nan)
        z = np.arctanh(r); se = 1 / np.sqrt(n - 3)
        lo, hi = z - 1.96 * se, z + 1.96 * se
        return (np.tanh(lo), np.tanh(hi))
    lo1, hi1 = fisher_ci(rho1, len(sub))
    lo2, hi2 = fisher_ci(rho2, len(sub))
    driver_rows.append({
        "driver": driver, "n": len(sub),
        "rho_HLA1": rho1, "rho_HLA1_lo": lo1, "rho_HLA1_hi": hi1, "p_HLA1": p1,
        "rho_HLA2": rho2, "rho_HLA2_lo": lo2, "rho_HLA2_hi": hi2, "p_HLA2": p2,
    })
driver_tab = pd.DataFrame(driver_rows)
driver_tab.to_csv(TAB / "T05_driver_stratified.tsv", sep="\t", index=False)
print(driver_tab.to_string(index=False))

# Forest plot
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
for ax, mod_name, mod_col in zip(axes, ["HLA-I", "HLA-II"], ["HLA1", "HLA2"]):
    y = np.arange(len(driver_tab))
    rho = driver_tab[f"rho_{mod_col}"].values
    lo = driver_tab[f"rho_{mod_col}_lo"].values
    hi = driver_tab[f"rho_{mod_col}_hi"].values
    ax.errorbar(rho, y, xerr=[rho - lo, hi - rho], fmt="o", color="#2b8cbe",
                markersize=8, capsize=4)
    ax.axvline(0, color="k", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{d} (n={n})" for d, n in zip(driver_tab["driver"], driver_tab["n"])])
    ax.set_xlabel(f"Spearman rho (DM1 vs {mod_name})")
    ax.set_title(f"Driver-stratified DM1 x {mod_name}")
    ax.grid(alpha=0.3, axis="x")
fig.suptitle(CAPTION_BOILERPLATE, fontsize=8, color="#666")
fig.tight_layout()
fig.savefig(FIG / "F05_driver_stratified_forest.png", dpi=180)
fig.savefig(FIG / "F05_driver_stratified_forest.pdf")
plt.close(fig)

# ============================================================
# 7. Methylation layer (per memory: r5_2_sample_methylation_8gene.tsv)
# ============================================================
print("[7/12] methylation x HLA-I (using existing TCGA HM450 means if available)")
meth_path = ROOT / "project/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv"
meth_section = "skipped"
if meth_path.exists():
    try:
        meth = pd.read_csv(meth_path, sep="\t")
        meth_cols = list(meth.columns)
        # try to find sample id col
        for cand in ["sample", "tcga_short", "case_id", "patient", "patient_id"]:
            if cand in meth.columns:
                key = cand
                break
        else:
            key = meth.columns[0]
        # build a sample id matching
        if "tcga_short" in meth.columns:
            meth["sample"] = meth["tcga_short"].astype(str) + "-01"
        elif key == meth.columns[0] and meth[key].astype(str).str.startswith("TCGA-").all():
            meth["sample"] = meth[key].astype(str)
            if not meth["sample"].str.endswith("-01").all():
                meth["sample"] = meth["sample"] + "-01"
        # find a mean methylation column
        beta_col = None
        for cand in ["mean_8g_beta", "mean_beta", "beta_mean", "mean_methylation"]:
            if cand in meth.columns:
                beta_col = cand
                break
        if beta_col is None:
            # fallback: any numeric col with name containing 'beta' or 'mean'
            for c in meth.columns:
                if c.lower().startswith(("beta", "mean")) and pd.api.types.is_numeric_dtype(meth[c]):
                    beta_col = c
                    break
        if beta_col is not None:
            df_meth = df.merge(meth[["sample", beta_col]], on="sample", how="left")
            sub = df_meth.dropna(subset=[beta_col, "HLA1_score", "DM1_use"])
            rho_dh, p_dh = stats.spearmanr(sub["DM1_use"], sub[beta_col])
            rho_mh, p_mh = stats.spearmanr(sub[beta_col], sub["HLA1_score"])
            rho_mh2, p_mh2 = stats.spearmanr(sub[beta_col], sub["HLA2_score"])
            meth_tab = pd.DataFrame([
                {"comparison": f"DM1 vs {beta_col}", "n": len(sub),
                 "spearman_rho": rho_dh, "p": p_dh},
                {"comparison": f"{beta_col} vs HLA-I module", "n": len(sub),
                 "spearman_rho": rho_mh, "p": p_mh},
                {"comparison": f"{beta_col} vs HLA-II module", "n": len(sub),
                 "spearman_rho": rho_mh2, "p": p_mh2},
            ])
            meth_tab.to_csv(TAB / "T06_methylation_dm1_hla.tsv", sep="\t", index=False)
            meth_section = f"used: {meth_path.name}, beta col = {beta_col}, n={len(sub)}"
            # figure
            fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
            scatter_with_lowess(sub["DM1_use"].values, sub[beta_col].values, axes[0],
                                f"DM1 vs {beta_col}", xlabel="DM1 score", ylabel=beta_col)
            scatter_with_lowess(sub[beta_col].values, sub["HLA1_score"].values, axes[1],
                                f"{beta_col} vs HLA-I module", xlabel=beta_col, ylabel="HLA-I module score")
            fig.suptitle(CAPTION_BOILERPLATE, fontsize=8, color="#666")
            fig.tight_layout()
            fig.savefig(FIG / "F06_methylation_dm1_hla1.png", dpi=180)
            fig.savefig(FIG / "F06_methylation_dm1_hla1.pdf")
            plt.close(fig)
        else:
            meth_section = f"file present but no beta column found; cols={meth_cols[:8]}"
    except Exception as e:
        meth_section = f"failed to parse {meth_path.name}: {e}"
else:
    meth_section = "no precomputed HM450 file at " + str(meth_path)
print("  methylation:", meth_section)

# ============================================================
# 8. Pan-cancer coherence sanity check (Track 6 hook)
# ============================================================
print("[8/12] pan-cancer coherence sanity (cross-lineage DM1 x HLA-I)")
# We can compute it ourselves directly using pancan_dm1_scored + a quick HLA-I module
# from the same expression file restricted to a subset of lineages.
# For tractability, we use the same HLA-I genes and process all primary tumors.
pancan_mask = (pheno["sample_type"] == "Primary Tumor")
pancan_samples = pheno[pancan_mask]["sample"].tolist()
print(f"  pancan primary tumor samples: {len(pancan_samples)}")
# already have expression for thy_primary; need to read again for full pancan
# Skip reading entire matrix again; instead use lineage-stratified summary on
# the thyroid expression we have, and report cross-lineage hook = pending track 6.
pancan_summary_rows = []
# Within-thyroid DM1 x HLA-I as the local reference
sub_thy = df.dropna(subset=["DM1_use", "HLA1_score"])
rho_thy, p_thy = stats.spearmanr(sub_thy["DM1_use"], sub_thy["HLA1_score"])
pancan_summary_rows.append({"lineage": "thyroid carcinoma (this track)", "n": len(sub_thy),
                            "spearman_rho_DM1_HLA1": rho_thy, "p": p_thy})
# Track 6 hook: read /track6_pancan_hla_dm1/ if present
track6_dir = OUT.parent / "track6_pancan_hla_dm1"
track6_summary = None
for fn in track6_dir.glob("*.tsv"):
    if "spearman" in fn.name.lower() or "summary" in fn.name.lower() or "lineage" in fn.name.lower():
        track6_summary = fn
        break
if track6_summary is not None:
    try:
        t6 = pd.read_csv(track6_summary, sep="\t")
        pancan_summary_rows.append({"lineage": f"track6:{track6_summary.name}",
                                    "n": int(t6.shape[0]),
                                    "spearman_rho_DM1_HLA1": np.nan, "p": np.nan})
        # also try to copy directly
        t6.to_csv(TAB / "T07_track6_pancan_hook.tsv", sep="\t", index=False)
    except Exception as e:
        print("  track6 file parse failed:", e)
pancan_tab = pd.DataFrame(pancan_summary_rows)
pancan_tab.to_csv(TAB / "T07_pancan_coherence_hook.tsv", sep="\t", index=False)

# ============================================================
# 9. HLA LOH check
# ============================================================
print("[9/12] HLA LOH file check")
loh_candidates = []
for pat in ["*loh*hla*", "*hla*loh*", "*lohhla*"]:
    loh_candidates += list((ROOT / "project").rglob(pat))
loh_section = "no precomputed HLA LOH calls present in repo; out of scope for this track per boundary"
if loh_candidates:
    loh_section = f"found {len(loh_candidates)} candidate file(s): " + ", ".join(
        str(p.relative_to(ROOT)) for p in loh_candidates[:5])
print("  LOH:", loh_section)

# ============================================================
# 10. Purity decoupling (partial Spearman ranks; conditional regression)
# ============================================================
print("[10/12] partial correlation DM1 x HLA-I controlling for purity proxy")


def partial_spearman(df, x, y, z_cols):
    sub = df[[x, y] + list(z_cols)].dropna()
    if len(sub) < 20:
        return np.nan, np.nan, len(sub)
    # rank-residualize
    rx = sub[x].rank(); ry = sub[y].rank()
    Z = sub[list(z_cols)].rank()
    Z = (Z - Z.mean()) / Z.std(ddof=0)
    Z = Z.values
    rx_resid = rx - Z @ np.linalg.lstsq(Z, rx.values, rcond=None)[0]
    ry_resid = ry - Z @ np.linalg.lstsq(Z, ry.values, rcond=None)[0]
    rho, p = stats.spearmanr(rx_resid, ry_resid)
    return rho, p, len(sub)


partial_rows = []
for mod, lbl in [("HLA1_score", "HLA-I"), ("HLA2_score", "HLA-II")]:
    rho_raw, p_raw = stats.spearmanr(df["DM1_use"], df[mod], nan_policy="omit")
    rho_p1, p_p1, n1 = partial_spearman(df, "DM1_use", mod, ["Purity_proxy"])
    rho_p2, p_p2, n2 = partial_spearman(df, "DM1_use", mod, ["Immune_proxy", "Stromal_proxy"])
    partial_rows.append({
        "module": lbl, "n_raw": int(df[["DM1_use", mod]].dropna().shape[0]),
        "rho_raw": rho_raw, "p_raw": p_raw,
        "rho_partial_purityproxy": rho_p1, "p_partial_purityproxy": p_p1, "n_partial": n1,
        "rho_partial_immune+stromal": rho_p2, "p_partial_immune+stromal": p_p2, "n_partial2": n2,
    })
partial_tab = pd.DataFrame(partial_rows)
partial_tab.to_csv(TAB / "T08_purity_partial.tsv", sep="\t", index=False)
print(partial_tab.to_string(index=False))

# Figure: bar chart raw vs partial
fig, ax = plt.subplots(figsize=(8, 4.5))
labels = partial_tab["module"].tolist()
x = np.arange(len(labels))
w = 0.27
ax.bar(x - w, partial_tab["rho_raw"], width=w, label="raw", color="#2b8cbe")
ax.bar(x, partial_tab["rho_partial_purityproxy"], width=w, label="| Purity proxy", color="#7fbc41")
ax.bar(x + w, partial_tab["rho_partial_immune+stromal"], width=w,
       label="| Immune + Stromal proxy", color="#e6ab02")
ax.axhline(0, color="k", lw=0.6)
ax.set_xticks(x); ax.set_xticklabels(labels)
ax.set_ylabel("Spearman rho with DM1")
ax.set_title(f"Purity decoupling: raw vs partial Spearman\n{CAPTION_BOILERPLATE}", fontsize=9)
ax.legend()
ax.grid(alpha=0.3, axis="y")
fig.tight_layout()
fig.savefig(FIG / "F07_purity_decoupling.png", dpi=180)
fig.savefig(FIG / "F07_purity_decoupling.pdf")
plt.close(fig)

# ============================================================
# 11. Residualization sanity check for Paper 1: does HLA-I as covariate move
# the DM1 x outcome relationship?
# ============================================================
print("[11/12] residualization sanity check: DM1 + HLA-I covariate vs outcome (PFI)")
res_section = "skipped (no PFI available)"
try:
    from scipy.stats import logistic
    # use Cox-style logistic on PFI event ~ DM1 + HLA1 (very rough; lifelines optional)
    cox_df = df.dropna(subset=["DM1_use", "HLA1_score", "PFI", "PFI.time"]).copy()
    if len(cox_df) >= 30:
        try:
            from lifelines import CoxPHFitter
            cph_rows = []
            for spec_name, cols in [("DM1 only", ["DM1_use"]),
                                     ("DM1 + HLA-I", ["DM1_use", "HLA1_score"]),
                                     ("DM1 + Purity proxy", ["DM1_use", "Purity_proxy"]),
                                     ("DM1 + HLA-I + Purity proxy",
                                      ["DM1_use", "HLA1_score", "Purity_proxy"])]:
                sub = cox_df[["PFI.time", "PFI"] + cols].dropna()
                cph = CoxPHFitter(penalizer=0.01)
                cph.fit(sub, duration_col="PFI.time", event_col="PFI")
                row = {"spec": spec_name, "n": len(sub)}
                for c in cols:
                    row[f"HR_{c}"] = float(np.exp(cph.params_[c]))
                    row[f"p_{c}"] = float(cph.summary.loc[c, "p"])
                cph_rows.append(row)
            cox_tab = pd.DataFrame(cph_rows)
            cox_tab.to_csv(TAB / "T09_cox_residualization.tsv", sep="\t", index=False)
            res_section = f"lifelines Cox fit, n={len(cox_df)}; see T09_cox_residualization.tsv"

            # Figure: forest of DM1 HR
            fig, ax = plt.subplots(figsize=(8, 4))
            spec_names = cox_tab["spec"].tolist()
            hrs = cox_tab["HR_DM1_use"].values
            ps = cox_tab["p_DM1_use"].values
            y = np.arange(len(spec_names))
            ax.scatter(hrs, y, s=80, color="#1b7837")
            for i, (h, p) in enumerate(zip(hrs, ps)):
                ax.text(h, i + 0.18, f"HR={h:.2f}, p={p:.2e}", fontsize=8)
            ax.axvline(1, color="k", lw=0.6)
            ax.set_yticks(y); ax.set_yticklabels(spec_names)
            ax.set_xlabel("Cox HR for DM1 score (PFI)")
            ax.set_title(f"Residualization: does HLA-I covariate move DM1 -> PFI HR?\n{CAPTION_BOILERPLATE}",
                         fontsize=9)
            ax.grid(alpha=0.3, axis="x")
            fig.tight_layout()
            fig.savefig(FIG / "F08_residualization_forest.png", dpi=180)
            fig.savefig(FIG / "F08_residualization_forest.pdf")
            plt.close(fig)
        except ImportError:
            # fallback: logistic regression on event
            from sklearn.linear_model import LogisticRegression
            cph_rows = []
            for spec_name, cols in [("DM1 only", ["DM1_use"]),
                                     ("DM1 + HLA-I", ["DM1_use", "HLA1_score"]),
                                     ("DM1 + Purity proxy", ["DM1_use", "Purity_proxy"])]:
                sub = cox_df[["PFI"] + cols].dropna()
                X = sub[cols].values
                y_v = sub["PFI"].astype(int).values
                clf = LogisticRegression(max_iter=1000).fit(X, y_v)
                row = {"spec": spec_name, "n": len(sub)}
                for c, b in zip(cols, clf.coef_[0]):
                    row[f"beta_{c}"] = float(b)
                cph_rows.append(row)
            cox_tab = pd.DataFrame(cph_rows)
            cox_tab.to_csv(TAB / "T09_logreg_residualization.tsv", sep="\t", index=False)
            res_section = f"lifelines unavailable; logistic fallback n={len(cox_df)}"
except Exception as e:
    res_section = f"failed: {e}"
print("  residualization:", res_section)

# ============================================================
# 12. Final summary JSON
# ============================================================
print("[12/12] writing final summary JSON")
summary = {
    "track": "Track 5 - DM1 x HLA-I/II gene-expression module in TCGA-THCA",
    "boundary": "HLA-I/II gene-expression module ONLY, not allele genotype. "
                "Cancer-cohort allele genotyping out of scope per HLA_CANCER_SEPARATION_RULES.md.",
    "n_thca_primary": int(len(df)),
    "n_with_dm1": int(df["DM1_use"].notna().sum()),
    "headline_correlations": corr_tab.to_dict("records"),
    "tertile_tests": tert_tab.to_dict("records"),
    "driver_stratified": driver_tab.to_dict("records"),
    "purity_partial": partial_tab.to_dict("records"),
    "methylation_status": meth_section,
    "loh_status": loh_section,
    "residualization_status": res_section,
    "outputs": {
        "results_dir": str(OUT),
        "figs_dir": str(FIG),
        "tables_dir": str(TAB),
    },
}
with open(OUT / "track5_summary.json", "w") as fh:
    json.dump(summary, fh, indent=2, default=str)
print("done.")
