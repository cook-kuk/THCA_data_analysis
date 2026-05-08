#!/usr/bin/env python3
"""
Track 21 - TCGA-THCA HLA-I locus allelic-imbalance read-ratio proxy x DM1.

BOUNDARY (CRITICAL — see project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md):
    This script analyses HLA-A / HLA-B / HLA-C as transcriptomic *gene-expression*
    signals only. We compute per-locus normalized read counts (TPM/log-CPM-style),
    relative ratios between loci, and module-level coefficients of variation. We
    DO NOT type HLA alleles from cancer RNA-seq. The "allelic imbalance proxy" is
    a *locus*-level imbalance (HLA-A vs HLA-B vs HLA-C), not an allele-level call.
    Allele-resolution AI from RNA would require allele-typed BAMs (LOHHLA / a phasing
    pipeline), which is forbidden in Paper 1 territory.

Caption boilerplate (must appear on every figure):
    "HLA-I locus allelic imbalance read-ratio proxy - transcriptomic signal,
     not allele genotype. Cancer-cohort allele genotyping is out of scope per
     separation rules."

Goal: refine Track 5 (DM1 ↔ HLA-I module +rho) by asking whether the three
class-I loci move together or whether DM1-high samples show preferential
silencing/dominance of one locus, plus B2M re-confirmation.
"""
from __future__ import annotations

import gzip
import json
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/hla_deepdive_2026_05_08/track21_hla_ai_proxy"
FIG = OUT / "figs"
TAB = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

CAPTION = (
    "HLA-I locus allelic imbalance read-ratio proxy - transcriptomic signal, not "
    "allele genotype. Cancer-cohort allele genotyping is out of scope per separation rules."
)

# ============================================================
# 1. Inputs (mirror Track 5's data plumbing)
# ============================================================
PANCAN_EXP = ROOT / "project/data/raw/TCGA_pancan/pancan_geneExp.gz"
PANCAN_PHENO = ROOT / "project/data/raw/TCGA_pancan/phenotype.tsv.gz"
PANCAN_DM1 = ROOT / "project/results/paper11_pancancer/pancan_dm1_scored.tsv"
DM_MASTER = ROOT / "project/results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv"

# Class-I loci (the three 'A:B:C' polymorphic loci) + B2M (light chain, classic LOH target)
HLA_I_LOCI = ["HLA-A", "HLA-B", "HLA-C"]
B2M_GENE = "B2M"
# Wider HLA-I module (matched to Track 5) for CV / heterogeneity calculation
HLA_I_MODULE = [
    "HLA-A", "HLA-B", "HLA-C",
    "B2M", "TAP1", "TAP2", "TAPBP", "NLRC5", "IRF1",
    "PSMB8", "PSMB9", "ERAP1", "ERAP2",
    "HLA-E", "HLA-F", "HLA-G",
    "CALR", "CANX", "PDIA3",
]
# Polymorphic core (the three classical loci) — the only true allelic-imbalance candidates
HLA_I_POLY = ["HLA-A", "HLA-B", "HLA-C"]
# 8-gene differentiation panel (DM1 driver)
PANEL8 = ["TG", "TPO", "TSHR", "SLC5A5", "FOXE1", "PAX8", "NKX2-1", "DIO1"]

# ============================================================
# 2. Load DM1 scaffolds, identify thyroid samples
# ============================================================
print("[1/12] loading DM1 master + phenotype")
dm_master = pd.read_csv(DM_MASTER, sep="\t")
dm_master["sample"] = dm_master["tcga_short"].astype(str) + "-01"
# dm_master has duplicate rows for some patients; collapse to first per sample
dm_master = dm_master.drop_duplicates(subset=["sample"], keep="first")

dm1_pancan = pd.read_csv(PANCAN_DM1, sep="\t")
dm1_thca = dm1_pancan[dm1_pancan["lineage"] == "thyroid carcinoma"].copy()
dm1_thca = dm1_thca.rename(columns={"DM1_like": "DM1_score"})

pheno = pd.read_csv(PANCAN_PHENO, sep="\t")
thy = pheno[pheno["_primary_disease"] == "thyroid carcinoma"][["sample", "sample_type"]].copy()
print(f"  thyroid samples by type: {thy['sample_type'].value_counts().to_dict()}")

primary_set = set(thy[thy["sample_type"] == "Primary Tumor"]["sample"])
normal_set = set(thy[thy["sample_type"] == "Solid Tissue Normal"]["sample"])
keep_set = primary_set | normal_set
print(f"  total tumor+normal kept: {len(keep_set)} (primary={len(primary_set)}, normal={len(normal_set)})")

# ============================================================
# 3. Stream expression for HLA-I module + 8-panel + B2M
# ============================================================
print("[2/12] streaming TCGA pancan expression matrix (gene x sample) for HLA-I module + PANEL8")
all_genes = set(HLA_I_MODULE + PANEL8 + [B2M_GENE])

with gzip.open(PANCAN_EXP, "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    sample_cols = header[1:]
    keep_idx = [i for i, s in enumerate(sample_cols) if s in keep_set]
    keep_samples = [sample_cols[i] for i in keep_idx]
    rows = {}
    for ln in fh:
        gene, *vals = ln.rstrip("\n").split("\t")
        if gene in all_genes:
            arr = np.array([vals[i] for i in keep_idx], dtype=float)
            rows[gene] = arr

expr = pd.DataFrame(rows, index=keep_samples)  # samples x genes
print(f"  expression matrix: {expr.shape[0]} samples x {expr.shape[1]} genes")
missing = [g for g in HLA_I_MODULE if g not in expr.columns]
print(f"  missing HLA-I module genes: {missing}")

# Build sample annotation
samp = pd.DataFrame(
    {
        "sample": expr.index,
        "sample_type": [
            "Primary Tumor" if s in primary_set else "Solid Tissue Normal" for s in expr.index
        ],
    }
)

# ============================================================
# 4. Compute DM1 score (8-gene anti-thyroid-diff)
# ============================================================
print("[3/12] computing DM1 score (anti-thyroid-diff)")


def zmean_module(df, gene_list, ref_index=None):
    sub = df.reindex(columns=[g for g in gene_list if g in df.columns])
    if sub.shape[1] == 0:
        return pd.Series(np.nan, index=df.index)
    if ref_index is None:
        mu = sub.mean()
        sd = sub.std(ddof=0)
    else:
        mu = sub.loc[ref_index].mean()
        sd = sub.loc[ref_index].std(ddof=0)
    z = (sub - mu) / sd
    return z.mean(axis=1)


# DM1 score reference = primary tumors (matches Track 5 / Paper 1 convention)
panel_z = zmean_module(expr, PANEL8, ref_index=samp.loc[samp["sample_type"] == "Primary Tumor", "sample"])
dm1_local = -panel_z  # high = de-differentiated
samp["DM1_local"] = dm1_local.values
samp = samp.merge(dm1_thca[["sample", "DM1_score"]], on="sample", how="left")
samp["DM1_use"] = samp["DM1_score"].fillna(samp["DM1_local"])

# Driver labels (primary only)
samp = samp.merge(
    dm_master[["sample", "driver_anchor_v17", "dm_status", "tert_pos"]],
    on="sample",
    how="left",
)


def simplify_driver(d):
    d = str(d)
    if d == "BRAF":
        return "BRAF"
    if d == "RAS":
        return "RAS"
    if d in ("Fusion", "fusion", "FUSION"):
        return "Fusion"
    if d in ("nan", "None", "", "unknown"):
        return "TripleNeg"
    return "Other"


samp["driver_simple"] = samp["driver_anchor_v17"].apply(simplify_driver)

# ============================================================
# 5. Per-locus HLA-A / -B / -C expression (Deliverable 1)
# ============================================================
print("[4/12] per-locus HLA-A / HLA-B / HLA-C expression table")

per_locus = expr[[g for g in HLA_I_POLY if g in expr.columns] + [B2M_GENE]].copy()
per_locus.columns = [c.replace("-", "_") for c in per_locus.columns]
per_locus = per_locus.reset_index().rename(columns={"index": "sample"})
out_perlocus = samp.merge(per_locus, on="sample", how="inner")
out_perlocus.to_csv(TAB / "T01_per_locus_expression.tsv", sep="\t", index=False)
print(f"  wrote T01_per_locus_expression.tsv (n={len(out_perlocus)})")

# ============================================================
# 6. Per-locus DM1 correlation (Deliverable 2)
# ============================================================
print("[5/12] per-locus expression x DM1 (Spearman) — primary tumors only")

prim = out_perlocus[out_perlocus["sample_type"] == "Primary Tumor"].copy()
norm = out_perlocus[out_perlocus["sample_type"] == "Solid Tissue Normal"].copy()


def spearman_with_dm1(df, gene_col):
    sub = df[[gene_col, "DM1_use"]].dropna()
    if len(sub) < 5:
        return dict(n=len(sub), rho=np.nan, p=np.nan)
    rho, p = stats.spearmanr(sub[gene_col], sub["DM1_use"])
    return dict(n=int(len(sub)), rho=float(rho), p=float(p))


per_locus_corr = []
for col, label in [("HLA_A", "HLA-A"), ("HLA_B", "HLA-B"), ("HLA_C", "HLA-C"), ("B2M", "B2M")]:
    r = spearman_with_dm1(prim, col)
    r.update({"locus": label, "cohort": "TCGA-THCA Primary Tumor"})
    per_locus_corr.append(r)
per_locus_corr = pd.DataFrame(per_locus_corr)[["locus", "cohort", "n", "rho", "p"]]
per_locus_corr.to_csv(TAB / "T02_per_locus_dm1_spearman.tsv", sep="\t", index=False)
print(per_locus_corr.to_string(index=False))


def fisher_ci(r, n):
    if n < 4 or pd.isna(r):
        return (np.nan, np.nan)
    z = np.arctanh(r)
    se = 1 / np.sqrt(n - 3)
    lo, hi = z - 1.96 * se, z + 1.96 * se
    return (float(np.tanh(lo)), float(np.tanh(hi)))


per_locus_corr[["rho_lo", "rho_hi"]] = per_locus_corr.apply(
    lambda row: pd.Series(fisher_ci(row["rho"], row["n"])), axis=1
)
per_locus_corr.to_csv(TAB / "T02_per_locus_dm1_spearman.tsv", sep="\t", index=False)

# Concordance: are A, B, C rho's similar?
abc = per_locus_corr[per_locus_corr["locus"].isin(["HLA-A", "HLA-B", "HLA-C"])]
abc_range = float(abc["rho"].max() - abc["rho"].min())
print(f"  rho range across A/B/C = {abc_range:.3f}")

# Figure 1: per-locus DM1 forest
fig, ax = plt.subplots(figsize=(7, 3.6))
y = np.arange(len(per_locus_corr))
ax.errorbar(
    per_locus_corr["rho"], y,
    xerr=[per_locus_corr["rho"] - per_locus_corr["rho_lo"], per_locus_corr["rho_hi"] - per_locus_corr["rho"]],
    fmt="o", color="#2b8cbe", markersize=8, capsize=4,
)
ax.axvline(0, color="k", lw=0.7)
ax.set_yticks(y)
ax.set_yticklabels([f"{l} (n={n})" for l, n in zip(per_locus_corr["locus"], per_locus_corr["n"])])
ax.set_xlabel("Spearman rho with DM1 score")
ax.set_title(f"Per-locus HLA-I expression x DM1 (TCGA-THCA primary)\n{CAPTION}", fontsize=8)
ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
fig.savefig(FIG / "F01_per_locus_dm1_forest.png", dpi=180)
fig.savefig(FIG / "F01_per_locus_dm1_forest.pdf")
plt.close(fig)

# Per-locus scatter panel
try:
    from statsmodels.nonparametric.smoothers_lowess import lowess
    HAVE_LOWESS = True
except Exception:
    HAVE_LOWESS = False


def scatter_lowess(ax, x, y, title, xlabel, ylabel):
    ax.scatter(x, y, s=8, alpha=0.45, color="#2b8cbe", edgecolor="none")
    if HAVE_LOWESS and len(x) > 30:
        sm = lowess(y, x, frac=0.4, return_sorted=True)
        ax.plot(sm[:, 0], sm[:, 1], color="#e34a33", lw=2)
    rho, p = stats.spearmanr(x, y)
    ax.set_title(f"{title}\nrho={rho:.3f}, p={p:.2e}, n={len(x)}", fontsize=9)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.3)


fig, axes = plt.subplots(1, 4, figsize=(16, 4.0))
for ax, gene_col, title in zip(
    axes,
    ["HLA_A", "HLA_B", "HLA_C", "B2M"],
    ["HLA-A", "HLA-B", "HLA-C", "B2M"],
):
    sub = prim[[gene_col, "DM1_use"]].dropna()
    scatter_lowess(ax, sub["DM1_use"].values, sub[gene_col].values,
                   f"DM1 vs {title}", "DM1 score", f"{title} log2(norm_count+1)")
fig.suptitle(CAPTION, fontsize=8, color="#666")
fig.tight_layout()
fig.savefig(FIG / "F02_per_locus_dm1_scatter.png", dpi=180)
fig.savefig(FIG / "F02_per_locus_dm1_scatter.pdf")
plt.close(fig)

# ============================================================
# 7. Per-locus ratio HLA-A : HLA-B : HLA-C (Deliverable 3)
# ============================================================
print("[6/12] per-sample HLA-A:B:C ratio + tumor vs normal distribution")

# Use raw expression (log2(norm_count+1) per Xena pancan); treat as approximate
# log-expression. We compare *linear* ratios after de-logging (2^x - 1, floor 0).
def linearise(x):
    return np.maximum(0, 2.0 ** x - 1.0)


for col in ["HLA_A", "HLA_B", "HLA_C", "B2M"]:
    out_perlocus[col + "_lin"] = linearise(out_perlocus[col].values)

eps = 1.0
out_perlocus["A_to_B"] = (out_perlocus["HLA_A_lin"] + eps) / (out_perlocus["HLA_B_lin"] + eps)
out_perlocus["A_to_C"] = (out_perlocus["HLA_A_lin"] + eps) / (out_perlocus["HLA_C_lin"] + eps)
out_perlocus["B_to_C"] = (out_perlocus["HLA_B_lin"] + eps) / (out_perlocus["HLA_C_lin"] + eps)
out_perlocus["log_A_to_B"] = np.log2(out_perlocus["A_to_B"])
out_perlocus["log_A_to_C"] = np.log2(out_perlocus["A_to_C"])
out_perlocus["log_B_to_C"] = np.log2(out_perlocus["B_to_C"])

# Total HLA-I polymorphic locus expression
out_perlocus["HLA_I_total_lin"] = out_perlocus[["HLA_A_lin", "HLA_B_lin", "HLA_C_lin"]].sum(axis=1)
for c in ["HLA_A_lin", "HLA_B_lin", "HLA_C_lin"]:
    out_perlocus[c.replace("_lin", "_frac")] = out_perlocus[c] / out_perlocus["HLA_I_total_lin"].replace(0, np.nan)

# Ratio summary table — tumor vs normal
ratio_rows = []
for col, label in [("log_A_to_B", "log2(A/B)"), ("log_A_to_C", "log2(A/C)"), ("log_B_to_C", "log2(B/C)")]:
    t = out_perlocus.loc[out_perlocus["sample_type"] == "Primary Tumor", col].dropna()
    n = out_perlocus.loc[out_perlocus["sample_type"] == "Solid Tissue Normal", col].dropna()
    u = stats.mannwhitneyu(t, n, alternative="two-sided")
    ratio_rows.append({
        "ratio": label,
        "n_tumor": len(t), "n_normal": len(n),
        "median_tumor": float(t.median()), "median_normal": float(n.median()),
        "iqr_tumor": float(t.quantile(0.75) - t.quantile(0.25)),
        "iqr_normal": float(n.quantile(0.75) - n.quantile(0.25)),
        "abs_med_tumor": float(t.abs().median()),
        "abs_med_normal": float(n.abs().median()),
        "MW_p_tumor_vs_normal": float(u.pvalue),
    })
ratio_tab = pd.DataFrame(ratio_rows)
ratio_tab.to_csv(TAB / "T03_locus_ratio_tumor_vs_normal.tsv", sep="\t", index=False)
print(ratio_tab.to_string(index=False))

# Figure: boxplot tumor vs normal for each ratio
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
for ax, col, label in zip(
    axes,
    ["log_A_to_B", "log_A_to_C", "log_B_to_C"],
    ["log2(A/B)", "log2(A/C)", "log2(B/C)"],
):
    t = out_perlocus.loc[out_perlocus["sample_type"] == "Primary Tumor", col].dropna()
    n = out_perlocus.loc[out_perlocus["sample_type"] == "Solid Tissue Normal", col].dropna()
    bp = ax.boxplot(
        [n, t], labels=[f"Normal\n(n={len(n)})", f"Tumor\n(n={len(t)})"],
        patch_artist=True, showfliers=False,
    )
    for patch, color in zip(bp["boxes"], ["#7fbc41", "#e6ab02"]):
        patch.set_facecolor(color)
    ax.axhline(0, color="k", lw=0.7, ls="--")
    p = stats.mannwhitneyu(t, n, alternative="two-sided").pvalue
    ax.set_title(f"{label}\nMW p={p:.2e}", fontsize=9)
    ax.set_ylabel(label)
    ax.grid(alpha=0.3, axis="y")
fig.suptitle(CAPTION, fontsize=8, color="#666")
fig.tight_layout()
fig.savefig(FIG / "F03_locus_ratio_tumor_vs_normal_box.png", dpi=180)
fig.savefig(FIG / "F03_locus_ratio_tumor_vs_normal_box.pdf")
plt.close(fig)

# ============================================================
# 8. Per-sample dominance label + DM1 high vs low (Deliverables 4, 5)
# ============================================================
print("[7/12] per-sample HLA-I locus dominance + DM1-tertile crosstab")

# Compute per-locus z within primary tumors (so dominance is relative to lineage)
prim_idx = out_perlocus["sample_type"] == "Primary Tumor"
zA = (out_perlocus.loc[prim_idx, "HLA_A"] - out_perlocus.loc[prim_idx, "HLA_A"].mean()) / out_perlocus.loc[prim_idx, "HLA_A"].std(ddof=0)
zB = (out_perlocus.loc[prim_idx, "HLA_B"] - out_perlocus.loc[prim_idx, "HLA_B"].mean()) / out_perlocus.loc[prim_idx, "HLA_B"].std(ddof=0)
zC = (out_perlocus.loc[prim_idx, "HLA_C"] - out_perlocus.loc[prim_idx, "HLA_C"].mean()) / out_perlocus.loc[prim_idx, "HLA_C"].std(ddof=0)

dom_df = pd.DataFrame({"sample": out_perlocus.loc[prim_idx, "sample"].values, "zA": zA.values, "zB": zB.values, "zC": zC.values})


def label_dominance(row, thresh=1.0):
    others = {
        "A": np.median([row["zB"], row["zC"]]),
        "B": np.median([row["zA"], row["zC"]]),
        "C": np.median([row["zA"], row["zB"]]),
    }
    if row["zA"] - others["A"] > thresh:
        return "A_dominant"
    if row["zB"] - others["B"] > thresh:
        return "B_dominant"
    if row["zC"] - others["C"] > thresh:
        return "C_dominant"
    # also check silenced (z below other two by > thresh) — relevant for AI proxy
    if others["A"] - row["zA"] > thresh:
        return "A_silenced"
    if others["B"] - row["zB"] > thresh:
        return "B_silenced"
    if others["C"] - row["zC"] > thresh:
        return "C_silenced"
    return "balanced"


dom_df["locus_label"] = dom_df.apply(label_dominance, axis=1)
print(f"  locus-label distribution (primary tumors): {dom_df['locus_label'].value_counts().to_dict()}")

prim_full = prim.merge(dom_df[["sample", "locus_label"]], on="sample", how="left")

# DM1 tertiles (primary only)
prim_full["DM1_tertile"] = pd.qcut(prim_full["DM1_use"], 3, labels=["Low", "Mid", "High"])

# Crosstab + chi-square
ct = pd.crosstab(prim_full["DM1_tertile"], prim_full["locus_label"])
ct.to_csv(TAB / "T04_dm1_tertile_x_locus_dominance.tsv", sep="\t")
chi2 = stats.chi2_contingency(ct.values)
print(f"  chi2={chi2.statistic:.3f}, dof={chi2.dof}, p={chi2.pvalue:.3g}")
ct_norm = ct.div(ct.sum(axis=1), axis=0)
ct_norm.to_csv(TAB / "T05_dm1_tertile_x_locus_dominance_proportions.tsv", sep="\t")

with open(TAB / "T05b_chi2_dm1_x_locus_dominance.txt", "w") as fh:
    fh.write(f"chi2 statistic: {chi2.statistic:.4f}\n")
    fh.write(f"dof: {chi2.dof}\n")
    fh.write(f"p-value: {chi2.pvalue:.4g}\n")
    fh.write(f"crosstab:\n{ct.to_string()}\n")

# Figure: stacked bar of locus-label proportion by DM1 tertile
fig, ax = plt.subplots(figsize=(8, 4.2))
labels_order = [c for c in ["A_dominant", "B_dominant", "C_dominant",
                            "A_silenced", "B_silenced", "C_silenced", "balanced"] if c in ct_norm.columns]
bottoms = np.zeros(len(ct_norm))
palette = {
    "A_dominant": "#1b9e77", "B_dominant": "#d95f02", "C_dominant": "#7570b3",
    "A_silenced": "#a6dba0", "B_silenced": "#fdae61", "C_silenced": "#bcbddc",
    "balanced": "#bdbdbd",
}
for col in labels_order:
    vals = ct_norm[col].values
    ax.bar(ct_norm.index.astype(str), vals, bottom=bottoms, color=palette[col], edgecolor="white",
           label=col)
    bottoms += vals
ax.set_ylabel("Fraction of samples")
ax.set_xlabel("DM1 tertile")
ax.set_title(f"HLA-I locus-dominance label x DM1 tertile (TCGA-THCA primary, n={len(prim_full)})\n"
             f"chi-square p={chi2.pvalue:.3g}\n{CAPTION}", fontsize=9)
ax.legend(bbox_to_anchor=(1.01, 1.0), loc="upper left", fontsize=8)
ax.grid(alpha=0.3, axis="y")
fig.tight_layout()
fig.savefig(FIG / "F04_dm1_tertile_locus_dominance_bar.png", dpi=180)
fig.savefig(FIG / "F04_dm1_tertile_locus_dominance_bar.pdf")
plt.close(fig)

# ============================================================
# 9. B2M re-confirmation (Deliverable 6)
# ============================================================
print("[8/12] B2M x DM1 re-confirmation")
b2m_sub = prim[["B2M", "DM1_use"]].dropna()
b2m_rho, b2m_p = stats.spearmanr(b2m_sub["B2M"], b2m_sub["DM1_use"])
b2m_pearson_r, b2m_pearson_p = stats.pearsonr(b2m_sub["B2M"], b2m_sub["DM1_use"])
print(f"  B2M vs DM1: spearman rho={b2m_rho:.3f} p={b2m_p:.2e} (Track 5 baseline rho≈+0.42)")

with open(TAB / "T06_b2m_dm1_reconfirm.tsv", "w") as fh:
    fh.write("metric\tvalue\n")
    fh.write(f"n\t{len(b2m_sub)}\n")
    fh.write(f"spearman_rho\t{b2m_rho:.4f}\n")
    fh.write(f"spearman_p\t{b2m_p:.4e}\n")
    fh.write(f"pearson_r\t{b2m_pearson_r:.4f}\n")
    fh.write(f"pearson_p\t{b2m_pearson_p:.4e}\n")
    fh.write("track5_reference_rho\t+0.42\n")

# ============================================================
# 10. HLA-I module heterogeneity (CV) (Deliverable 7)
# ============================================================
print("[9/12] HLA-I module per-sample CV (heterogeneity proxy)")
hla1_present = [g for g in HLA_I_MODULE if g in expr.columns]
sub_module = expr[hla1_present].copy()
# Linearise then CV (sd / mean) per sample on linear-space normalized counts
sub_module_lin = 2.0 ** sub_module - 1.0
sub_module_lin[sub_module_lin < 0] = 0
mu = sub_module_lin.mean(axis=1)
sd = sub_module_lin.std(axis=1, ddof=0)
cv = sd / mu.replace(0, np.nan)
cv = cv.rename("HLA_I_module_CV")
prim_full = prim_full.merge(cv.rename("HLA_I_module_CV").reset_index().rename(columns={"index": "sample"}), on="sample", how="left")

cv_dm1 = prim_full[["DM1_use", "HLA_I_module_CV"]].dropna()
cv_rho, cv_p = stats.spearmanr(cv_dm1["DM1_use"], cv_dm1["HLA_I_module_CV"])
print(f"  CV(HLA-I module) vs DM1: rho={cv_rho:.3f}, p={cv_p:.2e}, n={len(cv_dm1)}")

# Also tumor vs normal CV
norm_cv = cv.reindex(samp.loc[samp["sample_type"] == "Solid Tissue Normal", "sample"]).dropna()
prim_cv = cv.reindex(samp.loc[samp["sample_type"] == "Primary Tumor", "sample"]).dropna()
cv_mw = stats.mannwhitneyu(prim_cv, norm_cv, alternative="two-sided")
print(f"  CV tumor (median={prim_cv.median():.3f}) vs normal (median={norm_cv.median():.3f}); MW p={cv_mw.pvalue:.2e}")

with open(TAB / "T07_module_cv.tsv", "w") as fh:
    fh.write("metric\tvalue\n")
    fh.write(f"n_primary\t{len(prim_cv)}\n")
    fh.write(f"median_CV_primary\t{prim_cv.median():.4f}\n")
    fh.write(f"n_normal\t{len(norm_cv)}\n")
    fh.write(f"median_CV_normal\t{norm_cv.median():.4f}\n")
    fh.write(f"MW_p_tumor_vs_normal\t{cv_mw.pvalue:.4e}\n")
    fh.write(f"spearman_rho_DM1_vs_CV\t{cv_rho:.4f}\n")
    fh.write(f"spearman_p_DM1_vs_CV\t{cv_p:.4e}\n")

# Figure: CV scatter and tumor-vs-normal box
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
axes[0].boxplot([norm_cv, prim_cv], labels=[f"Normal\n(n={len(norm_cv)})", f"Tumor\n(n={len(prim_cv)})"],
                patch_artist=True, showfliers=False)
axes[0].set_ylabel("HLA-I module CV (per sample)")
axes[0].set_title(f"Module CV tumor vs normal\nMW p={cv_mw.pvalue:.2e}", fontsize=9)
axes[0].grid(alpha=0.3, axis="y")
scatter_lowess(axes[1], cv_dm1["DM1_use"].values, cv_dm1["HLA_I_module_CV"].values,
               "DM1 vs HLA-I module CV", "DM1 score", "Module CV")
fig.suptitle(CAPTION, fontsize=8, color="#666")
fig.tight_layout()
fig.savefig(FIG / "F05_module_CV_dm1_and_normal.png", dpi=180)
fig.savefig(FIG / "F05_module_CV_dm1_and_normal.pdf")
plt.close(fig)

# ============================================================
# 11. Allelic imbalance attempt (deliverable 8) — wall-stop pivot
# ============================================================
print("[10/12] allelic imbalance attempt (BAM/SNP search) - boundary check")
ai_section = {
    "attempted": False,
    "reason": (
        "(1) Allele-resolution allelic imbalance from RNA-seq requires either (a) "
        "allele-typed haplotype-phased BAMs over HLA-A/B/C exons (LOHHLA / a phasing "
        "pipeline) or (b) per-allele typing followed by SNP read-ratio. Both routes "
        "begin with HLA allele typing on cancer-cohort BAMs, which is FORBIDDEN in "
        "Paper 1 territory per HLA_CANCER_SEPARATION_RULES.md Section 1.2. (2) No "
        "TCGA-THCA BAM files exist locally (`find ... -name '*.bam'` returned 0 hits "
        "in project/data/external/tcga_thca/, no LOHHLA / OptiType outputs either). "
        "(3) Pivot per task brief: stop at the wall and report locus-level allelic "
        "imbalance proxy (HLA-A vs -B vs -C) which is gene-expression only — already "
        "captured in deliverables 1-7."
    ),
    "wall_pivot_executed": True,
    "what_was_done_instead": (
        "Per-locus expression z-scores, A:B:C ratios, locus-dominance label, "
        "DM1-tertile chi-square crosstab, B2M re-confirmation, HLA-I module CV "
        "(coefficient of variation) as locus-level imbalance proxy."
    ),
}
print("  AI-attempt:", ai_section["reason"][:140] + "...")

# ============================================================
# 12. Thorsson cross-reference (Deliverable 9) — best effort
# ============================================================
print("[11/12] Thorsson immune-subtype cross-reference (best effort)")
thorsson_section = {"available": False, "rho_found": None}
thorsson_path = ROOT / "project/results/v17p35/tables/FIX2_thorsson_subtype.tsv"
if thorsson_path.exists():
    try:
        th = pd.read_csv(thorsson_path, sep="\t")
        if "immune_subtype" in th.columns and th["immune_subtype"].astype(str).str.contains("MISSING").all():
            thorsson_section["available"] = False
            thorsson_section["note"] = (
                "FIX2_thorsson_subtype.tsv is a placeholder ('MISSING' for all rows); "
                "Thorsson 2018 per-sample HLA expression / neoantigen burden was not "
                "fetched locally in this sprint. Cross-validation deferred."
            )
        else:
            thorsson_section["available"] = True
            thorsson_section["note"] = "Thorsson table present; merged column inspection."
    except Exception as e:
        thorsson_section["note"] = f"failed to parse: {e}"
else:
    thorsson_section["note"] = f"no Thorsson file at {thorsson_path}"
print("  Thorsson:", thorsson_section)

# ============================================================
# Driver-stratified per-locus DM1 (additional figure beyond minimum)
# ============================================================
print("[12/12] driver-stratified per-locus DM1 forest (sanity)")
driver_rows = []
for driver in ["BRAF", "RAS", "Fusion", "TripleNeg"]:
    sub = prim[prim["driver_simple"] == driver]
    if len(sub) < 5:
        continue
    for col, label in [("HLA_A", "HLA-A"), ("HLA_B", "HLA-B"), ("HLA_C", "HLA-C"), ("B2M", "B2M")]:
        s2 = sub[[col, "DM1_use"]].dropna()
        if len(s2) < 4:
            continue
        rho, p = stats.spearmanr(s2[col], s2["DM1_use"])
        lo, hi = fisher_ci(rho, len(s2))
        driver_rows.append({
            "driver": driver, "locus": label, "n": int(len(s2)),
            "rho": float(rho), "p": float(p), "rho_lo": lo, "rho_hi": hi,
        })
driver_tab = pd.DataFrame(driver_rows)
driver_tab.to_csv(TAB / "T08_driver_stratified_per_locus.tsv", sep="\t", index=False)

if len(driver_tab):
    fig, ax = plt.subplots(figsize=(9, 5.0))
    drivers = list(driver_tab["driver"].unique())
    loci = ["HLA-A", "HLA-B", "HLA-C", "B2M"]
    palette = {"HLA-A": "#1b9e77", "HLA-B": "#d95f02", "HLA-C": "#7570b3", "B2M": "#e6ab02"}
    yp = 0
    yticks, ylabs = [], []
    for d in drivers:
        for L in loci:
            row = driver_tab[(driver_tab["driver"] == d) & (driver_tab["locus"] == L)]
            if len(row) == 0:
                continue
            r = row.iloc[0]
            ax.errorbar(r["rho"], yp, xerr=[[r["rho"] - r["rho_lo"]], [r["rho_hi"] - r["rho"]]],
                        fmt="o", color=palette[L], markersize=7, capsize=3)
            yticks.append(yp); ylabs.append(f"{d}/{L} (n={r['n']})")
            yp += 1
        yp += 0.5
    ax.axvline(0, color="k", lw=0.7)
    ax.set_yticks(yticks); ax.set_yticklabels(ylabs, fontsize=7)
    ax.set_xlabel("Spearman rho with DM1 score")
    ax.set_title(f"Driver-stratified per-locus DM1 correlation\n{CAPTION}", fontsize=9)
    ax.grid(alpha=0.3, axis="x")
    fig.tight_layout()
    fig.savefig(FIG / "F06_driver_stratified_per_locus_forest.png", dpi=180)
    fig.savefig(FIG / "F06_driver_stratified_per_locus_forest.pdf")
    plt.close(fig)

# ============================================================
# Summary JSON
# ============================================================
summary = {
    "track": "Track 21 - HLA-I locus allelic imbalance read-ratio proxy x DM1 (TCGA-THCA)",
    "boundary": (
        "HLA-I per-locus expression / inter-locus ratios / module CV ONLY, "
        "not allele genotype. Allele-typing pipelines were intentionally NOT run on "
        "TCGA-THCA per HLA_CANCER_SEPARATION_RULES.md Section 1.2."
    ),
    "n_primary": int(len(prim)),
    "n_normal": int(len(norm)),
    "per_locus_dm1": per_locus_corr.to_dict("records"),
    "rho_range_across_ABC": abc_range,
    "locus_ratio_tumor_vs_normal": ratio_tab.to_dict("records"),
    "dominance_label_distribution": dom_df["locus_label"].value_counts().to_dict(),
    "dm1_tertile_x_locus_dominance_chi2": {
        "chi2": float(chi2.statistic), "dof": int(chi2.dof), "p": float(chi2.pvalue),
    },
    "b2m_dm1_reconfirm": {
        "n": int(len(b2m_sub)), "spearman_rho": float(b2m_rho), "spearman_p": float(b2m_p),
        "track5_reference_rho": 0.42,
    },
    "module_cv": {
        "median_primary": float(prim_cv.median()),
        "median_normal": float(norm_cv.median()),
        "MW_p_tumor_vs_normal": float(cv_mw.pvalue),
        "spearman_rho_DM1_vs_CV": float(cv_rho),
        "spearman_p_DM1_vs_CV": float(cv_p),
    },
    "allelic_imbalance_attempt": ai_section,
    "thorsson_status": thorsson_section,
    "outputs": {
        "results_dir": str(OUT),
        "figs_dir": str(FIG),
        "tables_dir": str(TAB),
    },
}
with open(OUT / "track21_summary.json", "w") as fh:
    json.dump(summary, fh, indent=2, default=str)
print("done.")
