#!/usr/bin/env python3
"""
Track 27 — HLA-G + HLA-E + HLA-F non-classical class-I axis x DM1.

Boundary: HLA-A/B/C/E/F/G are analyzed as transcriptomic gene-expression modules,
NOT allele genotypes. Cancer-cohort allele genotyping is out of scope per
project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md (Sections 1.1, 1.2).

Reuses Track 5 (TCGA-THCA per-sample DM1 + module scores) and Track 6 (pancan
per-sample HLA module + DM1 z-scored), and Phase C ICI cohort processed expr.
"""
from __future__ import annotations

import gzip
import json
import os
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/hla_deepdive_2026_05_08/track27_hla_g_nonclassical"
FIG = OUT / "figs"
TAB = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

CAPTION_BOILERPLATE = (
    "HLA gene-expression module - not allele genotype. "
    "Cancer-cohort allele genotyping is out of scope per separation rules."
)

PANCAN_EXP = ROOT / "project/data/raw/TCGA_pancan/pancan_geneExp.gz"
PANCAN_PHENO = ROOT / "project/data/raw/TCGA_pancan/phenotype.tsv.gz"
PANCAN_SURV = ROOT / "project/data/raw/TCGA_pancan/survival.tsv"
PANCAN_DM1 = ROOT / "project/results/paper11_pancancer/pancan_dm1_scored.tsv"
DM_MASTER = ROOT / "project/results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv"
TRACK5_TAB = ROOT / "project/results/hla_deepdive_2026_05_08/track5_dm1_hla1_module/tables/T01_per_sample_module_scores.tsv"
TRACK6_TAB = ROOT / "project/results/hla_deepdive_2026_05_08/track6_pancan_hla_dm1/tables/T01b_pancan_full_module_genes_zscored.tsv"
METH_8G = ROOT / "project/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv"
ICI_DIR = ROOT / "project/results/paper11_pancancer/phase_C_ICI/processed"

CLASSICAL = ["HLA-A", "HLA-B", "HLA-C"]
NONCLASSICAL = ["HLA-E", "HLA-F", "HLA-G"]
HLA_G_LIGS = ["LILRB1", "LILRB2", "KIR2DL4"]   # ILT2, ILT4, KIR2DL4
HLA_E_LIGS = ["KLRC1", "KLRC2", "KLRD1"]       # NKG2A, NKG2C, CD94
ALL_GENES = list(set(CLASSICAL + NONCLASSICAL + HLA_G_LIGS + HLA_E_LIGS))


def cohens_d(x, y):
    x = np.asarray(x); y = np.asarray(y)
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return np.nan
    sx, sy = x.std(ddof=1), y.std(ddof=1)
    sp = np.sqrt(((nx - 1) * sx ** 2 + (ny - 1) * sy ** 2) / (nx + ny - 2))
    return (x.mean() - y.mean()) / sp if sp > 0 else np.nan


def fisher_ci(r, n):
    if n is None or n < 4 or pd.isna(r) or abs(r) >= 0.999999:
        return (np.nan, np.nan)
    z = np.arctanh(r); se = 1 / np.sqrt(n - 3)
    return (np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se))


# ============================================================
# 1. TCGA-THCA expression for HLA + ligands (primary tumor only)
# ============================================================
print("[1/13] loading TCGA-THCA expression for non-classical HLA + ligands")
pheno = pd.read_csv(PANCAN_PHENO, sep="\t")
thy_primary = pheno[(pheno["_primary_disease"] == "thyroid carcinoma") &
                    (pheno["sample_type"] == "Primary Tumor")]["sample"].tolist()
print(f"  TCGA-THCA primary samples: {len(thy_primary)}")

with gzip.open(PANCAN_EXP, "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    sample_cols = header[1:]
    keep_idx = [i for i, s in enumerate(sample_cols) if s in set(thy_primary)]
    keep_samples = [sample_cols[i] for i in keep_idx]
    rows = {}
    for ln in fh:
        gene, *vals = ln.rstrip("\n").split("\t")
        if gene in ALL_GENES:
            arr = np.array([vals[i] for i in keep_idx], dtype=float)
            rows[gene] = arr

expr = pd.DataFrame(rows, index=keep_samples)
print(f"  THCA expr matrix: {expr.shape}")
missing = [g for g in ALL_GENES if g not in expr.columns]
print(f"  missing genes: {missing}")

# Track 5 per-sample table for DM1, driver, etc
t5 = pd.read_csv(TRACK5_TAB, sep="\t")
df = expr.copy()
df.index.name = "sample"
df = df.reset_index().merge(t5[["sample", "DM1_use", "DM1_tertile", "driver_simple",
                                 "OS", "OS.time", "PFI", "PFI.time", "stage",
                                 "Purity_proxy", "Immune_proxy"]],
                             on="sample", how="left")
print(f"  THCA merged: {df.shape}")

# Compute classical / non-classical means (z-scored per gene)
for g in CLASSICAL + NONCLASSICAL:
    if g in df.columns:
        z = (df[g] - df[g].mean()) / df[g].std(ddof=0)
        df[f"z_{g}"] = z

df["classical_mean_z"] = df[[f"z_{g}" for g in CLASSICAL if f"z_{g}" in df.columns]].mean(axis=1)
df["nonclassical_mean_z"] = df[[f"z_{g}" for g in NONCLASSICAL if f"z_{g}" in df.columns]].mean(axis=1)
df["decoupling_HLAG"] = df["z_HLA-G"] - df["classical_mean_z"]
df["decoupling_HLAE"] = df["z_HLA-E"] - df["classical_mean_z"]
df["decoupling_HLAF"] = df["z_HLA-F"] - df["classical_mean_z"]

df.to_csv(TAB / "T01_thca_per_sample.tsv", sep="\t", index=False)
print("  wrote T01_thca_per_sample.tsv")

# ============================================================
# 2. Pairwise scatter classical vs non-classical + per-locus distribution
# ============================================================
print("[2/13] classical vs non-classical pairwise scatter + distributions")

pair_rows = []
fig, axes = plt.subplots(3, 3, figsize=(12, 12))
for i, c in enumerate(CLASSICAL):
    for j, nc in enumerate(NONCLASSICAL):
        ax = axes[i, j]
        sub = df[[c, nc]].dropna()
        rho, p = stats.spearmanr(sub[c], sub[nc])
        ax.scatter(sub[c], sub[nc], s=8, alpha=0.4, color="#2b8cbe", edgecolor="none")
        ax.set_xlabel(c); ax.set_ylabel(nc)
        ax.set_title(f"{c} vs {nc}: rho={rho:.3f} p={p:.1e} n={len(sub)}", fontsize=9)
        ax.grid(alpha=0.3)
        pair_rows.append({"x": c, "y": nc, "n": len(sub),
                          "spearman_rho": rho, "spearman_p": p})
fig.suptitle(f"TCGA-THCA: classical (HLA-A/B/C) vs non-classical (HLA-E/F/G)\n{CAPTION_BOILERPLATE}",
             fontsize=10)
fig.tight_layout()
fig.savefig(FIG / "F01_classical_vs_nonclassical_pairs.png", dpi=170)
fig.savefig(FIG / "F01_classical_vs_nonclassical_pairs.pdf")
plt.close(fig)
pd.DataFrame(pair_rows).to_csv(TAB / "T02_classical_nonclassical_pair_corr.tsv",
                                sep="\t", index=False)

# Per-locus expression distribution
fig, ax = plt.subplots(figsize=(9, 4.5))
data = [df[g].dropna().values for g in CLASSICAL + NONCLASSICAL if g in df.columns]
labels = [g for g in CLASSICAL + NONCLASSICAL if g in df.columns]
bp = ax.boxplot(data, labels=labels, patch_artist=True, showfliers=False)
for patch, lbl in zip(bp["boxes"], labels):
    patch.set_facecolor("#a6cee3" if lbl in CLASSICAL else "#fdbf6f")
ax.set_ylabel("log2(norm_count+1)")
ax.set_title(f"TCGA-THCA per-locus HLA-I expression\nClassical (blue) vs non-classical (orange)\n{CAPTION_BOILERPLATE}",
             fontsize=10)
ax.grid(alpha=0.3, axis="y")
fig.tight_layout()
fig.savefig(FIG / "F02_per_locus_distribution.png", dpi=170)
fig.savefig(FIG / "F02_per_locus_distribution.pdf")
plt.close(fig)

# ============================================================
# 3. DM1 x each non-classical HLA + side-by-side classical
# ============================================================
print("[3/13] DM1 x classical vs non-classical Spearman")
dm_rows = []
for g in CLASSICAL + NONCLASSICAL:
    if g not in df.columns:
        continue
    sub = df[[g, "DM1_use"]].dropna()
    rho, p = stats.spearmanr(sub[g], sub["DM1_use"])
    lo, hi = fisher_ci(rho, len(sub))
    dm_rows.append({"gene": g, "class": "classical" if g in CLASSICAL else "nonclassical",
                    "n": len(sub), "spearman_rho": rho, "rho_lo": lo, "rho_hi": hi,
                    "spearman_p": p})
dm_corr = pd.DataFrame(dm_rows)
dm_corr.to_csv(TAB / "T03_dm1_classical_nonclassical_corr.tsv", sep="\t", index=False)
print(dm_corr.to_string(index=False))

# Forest figure
fig, ax = plt.subplots(figsize=(8, 4.5))
y = np.arange(len(dm_corr))
colors = ["#1f78b4" if c == "classical" else "#ff7f00" for c in dm_corr["class"]]
ax.errorbar(dm_corr["spearman_rho"], y,
            xerr=[dm_corr["spearman_rho"] - dm_corr["rho_lo"],
                  dm_corr["rho_hi"] - dm_corr["spearman_rho"]],
            fmt="o", ecolor="gray", capsize=4, color="black", markersize=0)
ax.scatter(dm_corr["spearman_rho"], y, c=colors, s=80, edgecolor="black", zorder=3)
ax.axvline(0, color="k", lw=0.7)
ax.set_yticks(y)
ax.set_yticklabels([f"{g} (n={n})" for g, n in zip(dm_corr["gene"], dm_corr["n"])])
ax.set_xlabel("Spearman rho with DM1 score (TCGA-THCA)")
ax.set_title(f"DM1 x classical (blue) vs non-classical (orange) HLA-I\n{CAPTION_BOILERPLATE}",
             fontsize=10)
ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
fig.savefig(FIG / "F03_dm1_classical_nonclassical_forest.png", dpi=170)
fig.savefig(FIG / "F03_dm1_classical_nonclassical_forest.pdf")
plt.close(fig)

# ============================================================
# 4. Decoupling score
# ============================================================
print("[4/13] decoupling score: z(HLA-G) - z(mean HLA-A/B/C)")
decoup_rows = []
for ncl, dcol in [("HLA-G", "decoupling_HLAG"),
                  ("HLA-E", "decoupling_HLAE"),
                  ("HLA-F", "decoupling_HLAF")]:
    sub = df[[dcol, "DM1_use"]].dropna()
    rho, p = stats.spearmanr(sub[dcol], sub["DM1_use"])
    decoup_rows.append({"nonclassical_gene": ncl, "decoupling_col": dcol,
                        "n": len(sub),
                        "decoupling_mean": float(df[dcol].mean()),
                        "decoupling_sd": float(df[dcol].std()),
                        "spearman_rho_DM1_decoupling": rho,
                        "spearman_p": p})
decoup_tab = pd.DataFrame(decoup_rows)
decoup_tab.to_csv(TAB / "T04_decoupling_dm1.tsv", sep="\t", index=False)
print(decoup_tab.to_string(index=False))

fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
for ax, (ncl, dcol) in zip(axes, [("HLA-G", "decoupling_HLAG"),
                                   ("HLA-E", "decoupling_HLAE"),
                                   ("HLA-F", "decoupling_HLAF")]):
    sub = df[[dcol, "DM1_use"]].dropna()
    ax.scatter(sub["DM1_use"], sub[dcol], s=10, alpha=0.45, color="#e34a33",
               edgecolor="none")
    rho, p = stats.spearmanr(sub["DM1_use"], sub[dcol])
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("DM1 score")
    ax.set_ylabel(f"decoupling = z({ncl}) - z(mean A/B/C)")
    ax.set_title(f"DM1 x {ncl}-decoupling\nrho={rho:.3f} p={p:.2e} n={len(sub)}",
                 fontsize=10)
    ax.grid(alpha=0.3)
fig.suptitle(f"Classical/non-classical decoupling vs DM1 (TCGA-THCA)\n{CAPTION_BOILERPLATE}",
             fontsize=10)
fig.tight_layout()
fig.savefig(FIG / "F04_decoupling_vs_dm1.png", dpi=170)
fig.savefig(FIG / "F04_decoupling_vs_dm1.pdf")
plt.close(fig)

# Decoupling distribution (histogram)
fig, ax = plt.subplots(figsize=(8, 4.5))
for col, lbl, color in [("decoupling_HLAG", "HLA-G", "#e34a33"),
                        ("decoupling_HLAE", "HLA-E", "#1b7837"),
                        ("decoupling_HLAF", "HLA-F", "#762a83")]:
    ax.hist(df[col].dropna(), bins=40, alpha=0.4, label=lbl, color=color)
ax.axvline(0, color="k", lw=0.7)
ax.set_xlabel("Decoupling = z(non-classical) - z(mean classical)")
ax.set_ylabel("# samples")
ax.set_title(f"TCGA-THCA decoupling-score distribution per non-classical locus\n{CAPTION_BOILERPLATE}",
             fontsize=10)
ax.legend()
ax.grid(alpha=0.3, axis="y")
fig.tight_layout()
fig.savefig(FIG / "F05_decoupling_distribution.png", dpi=170)
fig.savefig(FIG / "F05_decoupling_distribution.pdf")
plt.close(fig)

# ============================================================
# 5. HLA-G ligand axis (LILRB1, LILRB2, KIR2DL4)
# ============================================================
print("[5/13] HLA-G ligand co-expression")
ligG_rows = []
for lig in HLA_G_LIGS:
    if lig not in df.columns:
        continue
    sub = df[["HLA-G", lig, "DM1_use"]].dropna()
    rho_g, p_g = stats.spearmanr(sub["HLA-G"], sub[lig])
    rho_d, p_d = stats.spearmanr(sub["DM1_use"], sub[lig])
    ligG_rows.append({"ligand": lig, "n": len(sub),
                      "rho_HLAG_lig": rho_g, "p_HLAG_lig": p_g,
                      "rho_DM1_lig": rho_d, "p_DM1_lig": p_d})
ligG_tab = pd.DataFrame(ligG_rows)
ligG_tab.to_csv(TAB / "T05_hlaG_ligand_corr.tsv", sep="\t", index=False)
print(ligG_tab.to_string(index=False))

# ============================================================
# 6. HLA-E ligand axis (KLRC1=NKG2A, KLRC2=NKG2C, KLRD1=CD94)
# ============================================================
print("[6/13] HLA-E ligand co-expression")
ligE_rows = []
for lig in HLA_E_LIGS:
    if lig not in df.columns:
        continue
    sub = df[["HLA-E", lig, "DM1_use"]].dropna()
    rho_e, p_e = stats.spearmanr(sub["HLA-E"], sub[lig])
    rho_d, p_d = stats.spearmanr(sub["DM1_use"], sub[lig])
    # NKG2A:NKG2C balance
    ligE_rows.append({"ligand": lig, "n": len(sub),
                      "rho_HLAE_lig": rho_e, "p_HLAE_lig": p_e,
                      "rho_DM1_lig": rho_d, "p_DM1_lig": p_d})
ligE_tab = pd.DataFrame(ligE_rows)
ligE_tab.to_csv(TAB / "T06_hlaE_ligand_corr.tsv", sep="\t", index=False)
print(ligE_tab.to_string(index=False))

# Combined ligand-axis figure
fig, axes = plt.subplots(2, 3, figsize=(13, 8))
for i, lig in enumerate(HLA_G_LIGS):
    ax = axes[0, i]
    if lig in df.columns:
        sub = df[["HLA-G", lig]].dropna()
        ax.scatter(sub["HLA-G"], sub[lig], s=10, alpha=0.45, color="#e34a33",
                   edgecolor="none")
        rho, p = stats.spearmanr(sub["HLA-G"], sub[lig])
        ax.set_title(f"HLA-G x {lig}\nrho={rho:.3f} p={p:.1e}", fontsize=10)
    ax.set_xlabel("HLA-G expression"); ax.set_ylabel(lig)
    ax.grid(alpha=0.3)
for i, lig in enumerate(HLA_E_LIGS):
    ax = axes[1, i]
    if lig in df.columns:
        sub = df[["HLA-E", lig]].dropna()
        ax.scatter(sub["HLA-E"], sub[lig], s=10, alpha=0.45, color="#1b7837",
                   edgecolor="none")
        rho, p = stats.spearmanr(sub["HLA-E"], sub[lig])
        ax.set_title(f"HLA-E x {lig}\nrho={rho:.3f} p={p:.1e}", fontsize=10)
    ax.set_xlabel("HLA-E expression"); ax.set_ylabel(lig)
    ax.grid(alpha=0.3)
fig.suptitle(f"Non-classical HLA-I ligand co-expression (TCGA-THCA)\n{CAPTION_BOILERPLATE}",
             fontsize=10)
fig.tight_layout()
fig.savefig(FIG / "F06_ligand_axis_grid.png", dpi=170)
fig.savefig(FIG / "F06_ligand_axis_grid.pdf")
plt.close(fig)

# ============================================================
# 7. DM1-high x HLA-G-high quadrant; survival HR
# ============================================================
print("[7/13] DM1-high x HLA-G-high quadrant + survival")
df["DM1_high"] = (df["DM1_use"] > df["DM1_use"].median()).astype("Int64")
df["HLAG_high"] = (df["HLA-G"] > df["HLA-G"].median()).astype("Int64")
df["quadrant"] = (
    df.apply(lambda r: ("DM1hi_HLAGhi" if r["DM1_high"] == 1 and r["HLAG_high"] == 1 else
                        ("DM1hi_HLAGlo" if r["DM1_high"] == 1 and r["HLAG_high"] == 0 else
                         ("DM1lo_HLAGhi" if r["DM1_high"] == 0 and r["HLAG_high"] == 1 else
                          ("DM1lo_HLAGlo" if r["DM1_high"] == 0 and r["HLAG_high"] == 0 else
                           pd.NA)))),
              axis=1)
)
quad_counts = df["quadrant"].value_counts(dropna=False)
print(quad_counts)

quad_rows = []
for q in ["DM1lo_HLAGlo", "DM1lo_HLAGhi", "DM1hi_HLAGlo", "DM1hi_HLAGhi"]:
    sub = df[df["quadrant"] == q]
    quad_rows.append({"quadrant": q, "n": len(sub),
                      "median_HLA-G": sub["HLA-G"].median() if len(sub) else np.nan,
                      "median_DM1": sub["DM1_use"].median() if len(sub) else np.nan,
                      "median_classical_z": sub["classical_mean_z"].median() if len(sub) else np.nan,
                      "median_nonclassical_z": sub["nonclassical_mean_z"].median() if len(sub) else np.nan})
quad_tab = pd.DataFrame(quad_rows)
quad_tab.to_csv(TAB / "T07_quadrant_summary_THCA.tsv", sep="\t", index=False)

# Survival: PFI HR per-quadrant (THCA usually low events; report descriptively)
try:
    from lifelines import CoxPHFitter
    cox_df = df.dropna(subset=["DM1_use", "HLA-G", "PFI", "PFI.time"]).copy()
    cox_df["DM1_use"] = (cox_df["DM1_use"] - cox_df["DM1_use"].mean()) / cox_df["DM1_use"].std(ddof=0)
    cox_df["HLA-G"] = (cox_df["HLA-G"] - cox_df["HLA-G"].mean()) / cox_df["HLA-G"].std(ddof=0)
    cox_df["DM1xHLAG"] = cox_df["DM1_use"] * cox_df["HLA-G"]
    surv_rows = []
    for spec, cols in [("DM1 only", ["DM1_use"]),
                        ("HLA-G only", ["HLA-G"]),
                        ("DM1 + HLA-G", ["DM1_use", "HLA-G"]),
                        ("DM1 + HLA-G + interaction", ["DM1_use", "HLA-G", "DM1xHLAG"])]:
        sub = cox_df[["PFI.time", "PFI"] + cols].dropna()
        if len(sub) < 30 or sub["PFI"].sum() < 5:
            surv_rows.append({"spec": spec, "n": len(sub), "events": int(sub["PFI"].sum()),
                              "note": "too few events"})
            continue
        cph = CoxPHFitter(penalizer=0.05)
        cph.fit(sub, duration_col="PFI.time", event_col="PFI")
        row = {"spec": spec, "n": len(sub), "events": int(sub["PFI"].sum())}
        for c in cols:
            row[f"HR_{c}"] = float(np.exp(cph.params_[c]))
            row[f"p_{c}"] = float(cph.summary.loc[c, "p"])
        surv_rows.append(row)
    surv_tab = pd.DataFrame(surv_rows)
    surv_tab.to_csv(TAB / "T08_survival_THCA_dm1_hlaG_PFI.tsv", sep="\t", index=False)
    print(surv_tab.to_string(index=False))

    # Quadrant log-rank-style descriptive curves (KaplanMeier)
    from lifelines import KaplanMeierFitter
    fig, ax = plt.subplots(figsize=(7, 5))
    kmf = KaplanMeierFitter()
    for q, color in [("DM1lo_HLAGlo", "#1b9e77"),
                      ("DM1lo_HLAGhi", "#7570b3"),
                      ("DM1hi_HLAGlo", "#d95f02"),
                      ("DM1hi_HLAGhi", "#e7298a")]:
        sub = df[(df["quadrant"] == q) &
                 df["PFI"].notna() & df["PFI.time"].notna()]
        if len(sub) >= 5:
            kmf.fit(sub["PFI.time"], event_observed=sub["PFI"], label=f"{q} (n={len(sub)})")
            kmf.plot_survival_function(ax=ax, ci_show=False, color=color)
    ax.set_xlabel("PFI time (days)")
    ax.set_ylabel("PFI-free survival")
    ax.set_title(f"TCGA-THCA: DM1 x HLA-G quadrants vs PFI\n{CAPTION_BOILERPLATE}",
                 fontsize=10)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG / "F07_quadrant_KM_THCA.png", dpi=170)
    fig.savefig(FIG / "F07_quadrant_KM_THCA.pdf")
    plt.close(fig)
    surv_section = "lifelines Cox + KM done; T08_survival_THCA_dm1_hlaG_PFI.tsv"
except Exception as e:
    surv_section = f"survival skipped: {e}"
print("  survival:", surv_section)

# ============================================================
# 8. Pan-cancer HLA-G distribution across 33 lineages
# ============================================================
print("[8/13] pan-cancer HLA-G distribution + per-lineage DM1 x HLA-G")
# Read pancan expression streaming for HLA-G/E/F + ligands across all primary tumors
pancan_primary = pheno[pheno["sample_type"] == "Primary Tumor"]["sample"].tolist()
print(f"  pancan primary samples: {len(pancan_primary)}")
with gzip.open(PANCAN_EXP, "rt") as fh:
    header2 = fh.readline().rstrip("\n").split("\t")
    sample_cols2 = header2[1:]
    keep_idx2 = [i for i, s in enumerate(sample_cols2) if s in set(pancan_primary)]
    keep_samples2 = [sample_cols2[i] for i in keep_idx2]
    rows2 = {}
    for ln in fh:
        gene, *vals = ln.rstrip("\n").split("\t")
        if gene in (NONCLASSICAL + CLASSICAL):
            arr = np.array([vals[i] for i in keep_idx2], dtype=float)
            rows2[gene] = arr

pancan_expr = pd.DataFrame(rows2, index=keep_samples2).reset_index().rename(columns={"index": "sample"})
# join lineage
lineage_map = pheno[["sample", "_primary_disease"]].rename(columns={"_primary_disease": "lineage"})
pancan_expr = pancan_expr.merge(lineage_map, on="sample", how="left")

# pancan DM1 (lineage-specific)
pancan_dm1 = pd.read_csv(PANCAN_DM1, sep="\t").rename(columns={"DM1_like": "DM1"})
pancan_expr = pancan_expr.merge(pancan_dm1[["sample", "DM1", "lineage"]],
                                 on=["sample", "lineage"], how="left")

# Per-lineage HLA-G distribution
lineage_summary = pancan_expr.groupby("lineage").agg(
    n=("HLA-G", "count"),
    HLAG_median=("HLA-G", "median"),
    HLAG_mean=("HLA-G", "mean"),
    HLAG_p90=("HLA-G", lambda x: np.nanpercentile(x, 90)),
    HLAE_median=("HLA-E", "median"),
    HLAF_median=("HLA-F", "median"),
).reset_index().sort_values("HLAG_median", ascending=False)
lineage_summary.to_csv(TAB / "T09_pancan_per_lineage_nonclassical.tsv", sep="\t", index=False)
print("  top 10 HLA-G lineages:")
print(lineage_summary.head(10).to_string(index=False))

# Per-lineage DM1 x HLA-G correlation
lin_corr_rows = []
for lin, g in pancan_expr.groupby("lineage"):
    sub = g[["DM1", "HLA-G", "HLA-E", "HLA-F"]].dropna()
    if len(sub) < 20:
        continue
    rg, pg = stats.spearmanr(sub["DM1"], sub["HLA-G"])
    re, pe = stats.spearmanr(sub["DM1"], sub["HLA-E"])
    rf, pf = stats.spearmanr(sub["DM1"], sub["HLA-F"])
    lin_corr_rows.append({"lineage": lin, "n": len(sub),
                          "rho_DM1_HLAG": rg, "p_DM1_HLAG": pg,
                          "rho_DM1_HLAE": re, "p_DM1_HLAE": pe,
                          "rho_DM1_HLAF": rf, "p_DM1_HLAF": pf})
lin_corr = pd.DataFrame(lin_corr_rows).sort_values("rho_DM1_HLAG", ascending=False)
lin_corr.to_csv(TAB / "T10_pancan_per_lineage_dm1_nonclassical_corr.tsv",
                sep="\t", index=False)

# Figure: pancan ranked HLA-G median
fig, ax = plt.subplots(figsize=(11, 7))
ls = lineage_summary.sort_values("HLAG_median", ascending=True)
y = np.arange(len(ls))
ax.barh(y, ls["HLAG_median"], color="#e34a33")
ax.set_yticks(y); ax.set_yticklabels(ls["lineage"], fontsize=7)
ax.set_xlabel("HLA-G median expression (log2 norm count)")
ax.set_title(f"Pan-cancer HLA-G median expression by lineage (n=33 lineages)\n{CAPTION_BOILERPLATE}",
             fontsize=10)
ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
fig.savefig(FIG / "F08_pancan_HLAG_by_lineage.png", dpi=170)
fig.savefig(FIG / "F08_pancan_HLAG_by_lineage.pdf")
plt.close(fig)

# Figure: pancan per-lineage rho(DM1, HLA-G) with thyroid highlighted
fig, ax = plt.subplots(figsize=(11, 7))
lc = lin_corr.sort_values("rho_DM1_HLAG", ascending=True)
y = np.arange(len(lc))
colors = ["#e34a33" if lin == "thyroid carcinoma" else "#2b8cbe" for lin in lc["lineage"]]
ax.barh(y, lc["rho_DM1_HLAG"], color=colors)
ax.set_yticks(y); ax.set_yticklabels(lc["lineage"], fontsize=7)
ax.set_xlabel("Spearman rho(DM1, HLA-G) per lineage")
ax.axvline(0, color="k", lw=0.7)
ax.set_title(f"Pan-cancer DM1 x HLA-G correlation by lineage (THCA in red)\n{CAPTION_BOILERPLATE}",
             fontsize=10)
ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
fig.savefig(FIG / "F09_pancan_DM1_HLAG_by_lineage.png", dpi=170)
fig.savefig(FIG / "F09_pancan_DM1_HLAG_by_lineage.pdf")
plt.close(fig)

# ============================================================
# 9. ICI cohorts: HLA-G vs response
# ============================================================
print("[9/13] ICI cohorts: HLA-G x response (4 phase-C cohorts)")
ici_rows = []
ici_per_sample_rows = []
ici_dirs = ["IMvigor210", "GSE176307", "riaz_GSE91061", "MGH_GSE115821"]
for cohort in ici_dirs:
    expr_p = ICI_DIR / cohort / "expr_matrix.tsv"
    meta_p = ICI_DIR / cohort / "metadata.tsv"
    if not (expr_p.exists() and meta_p.exists()):
        continue
    e = pd.read_csv(expr_p, sep="\t").set_index("gene_symbol")
    m = pd.read_csv(meta_p, sep="\t")
    common = [s for s in e.columns if s in m["sample_id"].values]
    e = e[common]
    # extract HLA-G/E/F + classical + ligands
    needed = [g for g in ALL_GENES if g in e.index]
    em = e.loc[needed].T.copy()
    em["sample_id"] = em.index
    sub = m.merge(em, on="sample_id", how="inner")
    # binary response var
    if "response_binary_CRPR_vs_SD_PD" in sub.columns:
        sub["resp"] = pd.to_numeric(sub["response_binary_CRPR_vs_SD_PD"], errors="coerce")
    else:
        sub["resp"] = np.nan
    # filter pre-treatment if available
    if "timepoint" in sub.columns:
        sub_pre = sub[sub["timepoint"].astype(str).str.lower().isin(["pre", "baseline", "screen"])]
        if len(sub_pre) >= 10:
            sub = sub_pre
    sub["cohort"] = cohort
    # decoupling
    cl_cols = [c for c in CLASSICAL if c in sub.columns]
    if cl_cols and "HLA-G" in sub.columns:
        z_cl = sub[cl_cols].apply(lambda x: (x - x.mean()) / x.std(ddof=0), axis=0)
        z_g = (sub["HLA-G"] - sub["HLA-G"].mean()) / sub["HLA-G"].std(ddof=0)
        sub["decoupling_HLAG"] = z_g - z_cl.mean(axis=1)
    ici_per_sample_rows.append(sub)
    # response stats
    s2 = sub.dropna(subset=["resp"])
    n_R = int((s2["resp"] == 1).sum()); n_NR = int((s2["resp"] == 0).sum())
    if n_R >= 3 and n_NR >= 3 and "HLA-G" in s2.columns:
        R = s2[s2["resp"] == 1]["HLA-G"]
        NR = s2[s2["resp"] == 0]["HLA-G"]
        u, p = stats.mannwhitneyu(R, NR, alternative="two-sided")
        d = cohens_d(R.values, NR.values)
        # decoupling
        if "decoupling_HLAG" in s2.columns:
            R_d = s2[s2["resp"] == 1]["decoupling_HLAG"]
            NR_d = s2[s2["resp"] == 0]["decoupling_HLAG"]
            u_d, p_d = stats.mannwhitneyu(R_d, NR_d, alternative="two-sided")
            d_decoup = cohens_d(R_d.values, NR_d.values)
        else:
            p_d, d_decoup = (np.nan, np.nan)
        ici_rows.append({"cohort": cohort, "n_R": n_R, "n_NR": n_NR,
                          "median_HLAG_R": float(R.median()),
                          "median_HLAG_NR": float(NR.median()),
                          "MW_p_HLAG": float(p), "cohens_d_HLAG_R_vs_NR": d,
                          "MW_p_decoupling": p_d, "cohens_d_decoupling": d_decoup})
ici_tab = pd.DataFrame(ici_rows)
ici_tab.to_csv(TAB / "T11_ici_HLAG_response.tsv", sep="\t", index=False)
if len(ici_per_sample_rows):
    pd.concat(ici_per_sample_rows, ignore_index=True, sort=False) \
        .to_csv(TAB / "T11b_ici_per_sample_nonclassical.tsv", sep="\t", index=False)
print(ici_tab.to_string(index=False))

# Forest figure ICI
if len(ici_tab):
    fig, ax = plt.subplots(figsize=(8, 4.5))
    y = np.arange(len(ici_tab))
    d = ici_tab["cohens_d_HLAG_R_vs_NR"].astype(float).values
    ax.scatter(d, y, s=80, color="#e34a33", edgecolor="black", zorder=3)
    ax.axvline(0, color="k", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{c} (R={r}, NR={nr})"
                       for c, r, nr in zip(ici_tab["cohort"], ici_tab["n_R"], ici_tab["n_NR"])])
    ax.set_xlabel("Cohen's d (HLA-G in Responder vs Non-responder)")
    ax.set_title(f"HLA-G expression x ICI response (4 cohorts)\n{CAPTION_BOILERPLATE}",
                 fontsize=10)
    ax.grid(alpha=0.3, axis="x")
    fig.tight_layout()
    fig.savefig(FIG / "F10_ici_HLAG_response_forest.png", dpi=170)
    fig.savefig(FIG / "F10_ici_HLAG_response_forest.pdf")
    plt.close(fig)

# ============================================================
# 10. Methylation x HLA-G (8-gene mean methylation as global hypomethylation proxy)
# ============================================================
print("[10/13] methylation x HLA-G (8-gene mean as proxy for global hypomethylation)")
meth_section = "skipped"
try:
    meth = pd.read_csv(METH_8G, sep="\t")
    meth["sample"] = meth["sample_short"].astype(str) + "-01"
    df_meth = df.merge(meth[["sample", "mean_8g_beta"]], on="sample", how="left")
    sub = df_meth.dropna(subset=["mean_8g_beta", "HLA-G", "DM1_use"])
    rows_m = []
    for g in NONCLASSICAL + CLASSICAL:
        if g not in sub.columns:
            continue
        rho, p = stats.spearmanr(sub["mean_8g_beta"], sub[g])
        rows_m.append({"gene": g, "n": len(sub),
                       "rho_meth_gene": rho, "p": p,
                       "class": "classical" if g in CLASSICAL else "nonclassical"})
    rho_dm_meth, p_dm_meth = stats.spearmanr(sub["DM1_use"], sub["mean_8g_beta"])
    rho_dm_g, p_dm_g = stats.spearmanr(sub["DM1_use"], sub["HLA-G"])
    rho_g_meth, p_g_meth = stats.spearmanr(sub["HLA-G"], sub["mean_8g_beta"])
    meth_tab = pd.DataFrame(rows_m)
    meth_tab.to_csv(TAB / "T12_meth_per_gene.tsv", sep="\t", index=False)

    summary_rows = [
        {"comparison": "DM1 vs mean_8g_beta", "n": len(sub), "rho": rho_dm_meth, "p": p_dm_meth},
        {"comparison": "DM1 vs HLA-G", "n": len(sub), "rho": rho_dm_g, "p": p_dm_g},
        {"comparison": "mean_8g_beta vs HLA-G", "n": len(sub), "rho": rho_g_meth, "p": p_g_meth},
    ]
    pd.DataFrame(summary_rows).to_csv(TAB / "T13_meth_dm1_hlaG_summary.tsv",
                                       sep="\t", index=False)

    # Partial: HLA-G ~ DM1 controlling for mean_8g_beta (rank-based)
    rx = sub["HLA-G"].rank(); ry = sub["DM1_use"].rank(); rz = sub["mean_8g_beta"].rank()
    Z = (rz - rz.mean()) / rz.std(ddof=0)
    rx_resid = rx - Z * np.linalg.lstsq(Z.values.reshape(-1, 1), rx.values, rcond=None)[0][0]
    ry_resid = ry - Z * np.linalg.lstsq(Z.values.reshape(-1, 1), ry.values, rcond=None)[0][0]
    rho_part, p_part = stats.spearmanr(rx_resid, ry_resid)
    summary_rows.append({"comparison": "HLA-G vs DM1 partial | mean_8g_beta",
                         "n": len(sub), "rho": rho_part, "p": p_part})
    pd.DataFrame(summary_rows).to_csv(TAB / "T13_meth_dm1_hlaG_summary.tsv",
                                       sep="\t", index=False)

    # Figure
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    for ax, (x, y, lbl) in zip(axes, [
        ("DM1_use", "mean_8g_beta", "DM1 vs 8-gene mean methylation"),
        ("mean_8g_beta", "HLA-G", "mean_8g_beta vs HLA-G"),
        ("DM1_use", "HLA-G", "DM1 vs HLA-G")]):
        s = sub[[x, y]].dropna()
        ax.scatter(s[x], s[y], s=10, alpha=0.45, color="#2b8cbe", edgecolor="none")
        rho, p = stats.spearmanr(s[x], s[y])
        ax.set_xlabel(x); ax.set_ylabel(y)
        ax.set_title(f"{lbl}\nrho={rho:.3f} p={p:.2e} n={len(s)}", fontsize=10)
        ax.grid(alpha=0.3)
    fig.suptitle(f"TCGA-THCA: methylation x DM1 x HLA-G (8-gene mean methylation as global hypomethylation proxy)\n{CAPTION_BOILERPLATE}",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "F11_methylation_dm1_hlaG.png", dpi=170)
    fig.savefig(FIG / "F11_methylation_dm1_hlaG.pdf")
    plt.close(fig)
    meth_section = (
        f"used 8-gene mean methylation as proxy (HLA-G-specific HM450 not in repo); "
        f"n={len(sub)}; HLA-G vs mean_8g_beta rho={rho_g_meth:.3f} p={p_g_meth:.2e}"
    )
except Exception as e:
    meth_section = f"failed: {e}"
print("  methylation:", meth_section)

# ============================================================
# 11. Driver-stratified HLA-G
# ============================================================
print("[11/13] driver-stratified HLA-G")
drv_rows = []
for drv in ["BRAF", "RAS", "Fusion", "TripleNeg"]:
    sub = df[df["driver_simple"] == drv].dropna(subset=["HLA-G", "DM1_use"])
    if len(sub) < 5:
        continue
    rho, p = stats.spearmanr(sub["DM1_use"], sub["HLA-G"])
    drv_rows.append({"driver": drv, "n": len(sub),
                      "median_HLAG": float(sub["HLA-G"].median()),
                      "mean_HLAG": float(sub["HLA-G"].mean()),
                      "median_decoupling_HLAG": float(sub["decoupling_HLAG"].median()),
                      "rho_DM1_HLAG": rho, "p_DM1_HLAG": p})
drv_tab = pd.DataFrame(drv_rows)
drv_tab.to_csv(TAB / "T14_driver_stratified_HLAG.tsv", sep="\t", index=False)
print(drv_tab.to_string(index=False))

# pairwise tests TripleNeg vs BRAF
test_rows = []
ref = df[df["driver_simple"] == "TripleNeg"].dropna(subset=["HLA-G"])
for drv in ["BRAF", "RAS", "Fusion"]:
    cur = df[df["driver_simple"] == drv].dropna(subset=["HLA-G"])
    if len(cur) >= 3 and len(ref) >= 3:
        u, p = stats.mannwhitneyu(ref["HLA-G"], cur["HLA-G"], alternative="two-sided")
        d = cohens_d(ref["HLA-G"].values, cur["HLA-G"].values)
        test_rows.append({"contrast": f"TripleNeg vs {drv}",
                          "n_ref": len(ref), "n_cur": len(cur),
                          "MW_p": p, "cohens_d_TripleNeg_minus_other": d})
pd.DataFrame(test_rows).to_csv(TAB / "T15_driver_pairwise_HLAG.tsv", sep="\t", index=False)

# Figure
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
order = ["BRAF", "RAS", "Fusion", "TripleNeg"]
data_g = [df[df["driver_simple"] == d]["HLA-G"].dropna().values for d in order]
data_d = [df[df["driver_simple"] == d]["decoupling_HLAG"].dropna().values for d in order]
bp = axes[0].boxplot(data_g, labels=order, patch_artist=True,
                      boxprops=dict(facecolor="#fdbf6f"), showfliers=False)
axes[0].set_ylabel("HLA-G expression"); axes[0].set_title("HLA-G by driver class")
axes[0].grid(alpha=0.3, axis="y")
bp2 = axes[1].boxplot(data_d, labels=order, patch_artist=True,
                       boxprops=dict(facecolor="#a6cee3"), showfliers=False)
axes[1].set_ylabel("decoupling_HLAG = z(HLA-G) - z(mean A/B/C)")
axes[1].set_title("HLA-G decoupling by driver class")
axes[1].grid(alpha=0.3, axis="y")
fig.suptitle(f"TCGA-THCA driver-stratified HLA-G\n{CAPTION_BOILERPLATE}",
             fontsize=10)
fig.tight_layout()
fig.savefig(FIG / "F12_driver_stratified_HLAG.png", dpi=170)
fig.savefig(FIG / "F12_driver_stratified_HLAG.pdf")
plt.close(fig)

# ============================================================
# 12. Final summary
# ============================================================
print("[12/13] writing track27_summary.json")
summary = {
    "track": "Track 27 - HLA-G + HLA-E + HLA-F non-classical class-I axis x DM1",
    "boundary": (
        "HLA gene-expression module ONLY, not allele genotype. "
        "Cancer-cohort allele genotyping out of scope per HLA_CANCER_SEPARATION_RULES.md."
    ),
    "n_thca_primary": int(len(df)),
    "dm1_classical_nonclassical_corr": dm_corr.to_dict("records"),
    "decoupling_dm1": decoup_tab.to_dict("records"),
    "hlaG_ligand": ligG_tab.to_dict("records"),
    "hlaE_ligand": ligE_tab.to_dict("records"),
    "quadrant_summary_THCA": quad_tab.to_dict("records"),
    "survival_status": surv_section,
    "pancan_top10_HLAG_lineages": lineage_summary.head(10).to_dict("records"),
    "pancan_per_lineage_dm1_HLAG_top5": lin_corr.head(5).to_dict("records"),
    "pancan_per_lineage_dm1_HLAG_thyroid_row": lin_corr[
        lin_corr["lineage"] == "thyroid carcinoma"
    ].to_dict("records"),
    "ici_HLAG_response": ici_tab.to_dict("records"),
    "methylation_status": meth_section,
    "driver_stratified_HLAG": drv_tab.to_dict("records"),
    "outputs": {"figs": str(FIG), "tables": str(TAB)},
}
with open(OUT / "track27_summary.json", "w") as fh:
    json.dump(summary, fh, indent=2, default=str)
print("[13/13] done.")
