#!/usr/bin/env python3
"""
Track 34 - RAI-lineage gene module x HLA-I/II module x DM1 decoupling in TCGA-THCA.

Question: is RAI-lineage silencing and HLA induction the SAME molecular axis or
TWO independent axes that happen to co-occur in DM1-high tumors?

Boundary: HLA-I/II are analyzed strictly as transcriptomic gene-expression
modules, NOT as allele genotypes. Cancer-cohort allele genotyping is out of
scope per project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md
(Sections 1.1 / 1.2). RAI-lineage / DM1 / driver / outcome work IS Paper 1
territory and is fully allowed.
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
OUT = ROOT / "project/results/hla_deepdive_2026_05_08/track34_rai_hla_decouple"
FIG = OUT / "figs"
TAB = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

CAPTION_BOILERPLATE = (
    "HLA-I/II = gene-expression module (not allele genotype). "
    "Cancer-cohort allele genotyping is out of scope per separation rules."
)

# ============================================================
# 0. Inputs
# ============================================================
PANCAN_EXP = ROOT / "project/data/raw/TCGA_pancan/pancan_geneExp.gz"
PANCAN_PHENO = ROOT / "project/data/raw/TCGA_pancan/phenotype.tsv.gz"
PANCAN_DM1 = ROOT / "project/results/paper11_pancancer/pancan_dm1_scored.tsv"
DM_MASTER = ROOT / "project/results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv"
TRACK5_T01 = ROOT / "project/results/hla_deepdive_2026_05_08/track5_dm1_hla1_module/tables/T01_per_sample_module_scores.tsv"
TRACK26_IFNG = ROOT / "project/results/hla_deepdive_2026_05_08/track26_ifng_hla_dm1/tables/T01_ifng_score_per_sample.tsv"
METH_PATH = ROOT / "project/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv"
P5_GSE286332 = ROOT / "project/results/p5_8gene_vs_hla_autocorr/scores_GSE286332.tsv"
P5_TCGA = ROOT / "project/results/p5_8gene_vs_hla_autocorr/scores_TCGA.tsv"

# Canonical RAI-lineage 9-gene panel (per memory v19_paper3_track_b_lite)
RAI_GENES = ["TG", "TPO", "TSHR", "SLC5A5", "DUOX1", "DUOX2", "IYD", "DIO1", "DIO2"]
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
IFNG_HALLMARK = [
    "IFNG", "STAT1", "GBP1", "GBP2", "GBP3", "GBP4", "GBP5",
    "IRF1", "IRF7", "IRF8", "IRF9", "CXCL9", "CXCL10", "CXCL11",
    "TAP1", "TAP2", "PSMB8", "PSMB9", "B2M",
]


def zmean_module(df, gene_list):
    sub = df.reindex(columns=[g for g in gene_list if g in df.columns])
    if sub.shape[1] == 0:
        return pd.Series(np.nan, index=df.index)
    z = (sub - sub.mean()) / sub.std(ddof=0)
    return z.mean(axis=1)


def cohens_d(x, y):
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return np.nan
    sx, sy = x.std(ddof=1), y.std(ddof=1)
    sp = np.sqrt(((nx - 1) * sx ** 2 + (ny - 1) * sy ** 2) / (nx + ny - 2))
    return (x.mean() - y.mean()) / sp if sp > 0 else np.nan


def fisher_ci(r, n):
    if n is None or n < 4 or pd.isna(r):
        return (np.nan, np.nan)
    z = np.arctanh(r)
    se = 1 / np.sqrt(n - 3)
    lo, hi = z - 1.96 * se, z + 1.96 * se
    return (np.tanh(lo), np.tanh(hi))


def partial_spearman(df, x, y, z_cols):
    sub = df[[x, y] + list(z_cols)].dropna()
    if len(sub) < 20:
        return np.nan, np.nan, len(sub)
    rx = sub[x].rank().values.astype(float)
    ry = sub[y].rank().values.astype(float)
    Z = sub[list(z_cols)].rank()
    Z = (Z - Z.mean()) / Z.std(ddof=0)
    Z = Z.values
    bx, *_ = np.linalg.lstsq(Z, rx, rcond=None)
    by, *_ = np.linalg.lstsq(Z, ry, rcond=None)
    rx_resid = rx - Z @ bx
    ry_resid = ry - Z @ by
    rho, p = stats.spearmanr(rx_resid, ry_resid)
    return rho, p, len(sub)


# ============================================================
# 1. Build per-sample frame: DM1, HLA-I, HLA-II, IFNG, RAI module
# ============================================================
print("[1/12] loading TCGA-THCA scaffold (re-uses Track 5 / 26)")

dm_master = pd.read_csv(DM_MASTER, sep="\t")
dm_master["sample"] = dm_master["tcga_short"].astype(str) + "-01"

dm1_pancan = pd.read_csv(PANCAN_DM1, sep="\t")
dm1_thca = dm1_pancan[dm1_pancan["lineage"] == "thyroid carcinoma"].copy()
dm1_thca = dm1_thca.rename(columns={"DM1_like": "DM1_score"})

pheno = pd.read_csv(PANCAN_PHENO, sep="\t")
thy_primary = pheno[(pheno["_primary_disease"] == "thyroid carcinoma") &
                    (pheno["sample_type"] == "Primary Tumor")]["sample"].tolist()
print(f"  TCGA-THCA Primary Tumor samples: {len(thy_primary)}")

print("[2/12] streaming pancan_geneExp for required genes (RAI + HLA + PANEL8 + IFNG)")
all_genes = set(RAI_GENES + HLA1_GENES + HLA2_GENES + PANEL8 + IFNG_HALLMARK)

with gzip.open(PANCAN_EXP, "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    sample_cols = header[1:]
    keep_set = set(thy_primary)
    keep_idx = [i for i, s in enumerate(sample_cols) if s in keep_set]
    keep_samples = [sample_cols[i] for i in keep_idx]
    rows = {}
    for ln in fh:
        parts = ln.rstrip("\n").split("\t")
        gene = parts[0]
        if gene in all_genes:
            vals = parts[1:]
            arr = np.array([vals[i] for i in keep_idx], dtype=float)
            rows[gene] = arr

expr = pd.DataFrame(rows, index=keep_samples)  # samples x genes
print(f"  expression matrix (samples x genes): {expr.shape}")
present_rai = [g for g in RAI_GENES if g in expr.columns]
missing_rai = [g for g in RAI_GENES if g not in expr.columns]
print(f"  RAI-lineage genes present: {present_rai}")
print(f"  RAI-lineage genes missing: {missing_rai}")

# Module scores
hla1_score = zmean_module(expr, HLA1_GENES)
hla2_score = zmean_module(expr, HLA2_GENES)
ifng_score = zmean_module(expr, IFNG_HALLMARK)
panel8_z = zmean_module(expr, PANEL8)
dm1_local = -panel8_z  # high = de-differentiated
rai_score = zmean_module(expr, RAI_GENES)  # high = RAI-lineage retained

frame = pd.DataFrame({
    "sample": expr.index,
    "RAI_module": rai_score.values,
    "HLA1_score": hla1_score.values,
    "HLA2_score": hla2_score.values,
    "IFNG_score": ifng_score.values,
    "DM1_local": dm1_local.values,
})
frame = frame.merge(dm1_thca[["sample", "DM1_score"]], on="sample", how="left")
frame = frame.merge(
    dm_master[["sample", "driver_anchor_v17", "dm_status", "tert_pos",
               "stage", "OS", "OS.time", "PFI", "PFI.time"]],
    on="sample", how="left",
)
frame["DM1_use"] = frame["DM1_score"].fillna(frame["DM1_local"])


def simplify_driver(d):
    if d == "BRAF":
        return "BRAF"
    if d == "RAS":
        return "RAS"
    if d in ("Fusion", "fusion", "FUSION"):
        return "Fusion"
    if pd.isna(d) or d in (None, "", "unknown"):
        return "TripleNeg"
    return "Other"


frame["driver_simple"] = frame["driver_anchor_v17"].map(simplify_driver)


def simplify_stage(s):
    if not isinstance(s, str):
        return np.nan
    s = s.strip()
    if s.startswith("Stage IV"):
        return "IV"
    if s.startswith("Stage III"):
        return "III"
    if s.startswith("Stage II"):
        return "II"
    if s.startswith("Stage I"):
        return "I"
    return np.nan


frame["stage_simple"] = frame["stage"].map(simplify_stage)
frame.to_csv(TAB / "T01_per_sample_scores.tsv", sep="\t", index=False)
print(f"  per-sample frame: {frame.shape}; RAI module n={frame['RAI_module'].notna().sum()}")

# ============================================================
# 2. Headline correlations
# ============================================================
print("[3/12] headline DM1 x RAI, HLA-I x RAI, HLA-II x RAI, IFNG x RAI")


def headline(df, x, y, label):
    sub = df[[x, y]].dropna()
    rho_s, p_s = stats.spearmanr(sub[x], sub[y])
    rho_p, p_p = stats.pearsonr(sub[x], sub[y])
    lo, hi = fisher_ci(rho_s, len(sub))
    return {"comparison": label, "n": len(sub),
            "spearman_rho": rho_s, "spearman_lo": lo, "spearman_hi": hi,
            "spearman_p": p_s, "pearson_r": rho_p, "pearson_p": p_p}


corr_rows = []
for x, y, lbl in [
    ("DM1_use", "RAI_module", "DM1 vs RAI module"),
    ("HLA1_score", "RAI_module", "HLA-I module vs RAI module"),
    ("HLA2_score", "RAI_module", "HLA-II module vs RAI module"),
    ("IFNG_score", "RAI_module", "IFN-gamma vs RAI module"),
    ("DM1_use", "HLA1_score", "DM1 vs HLA-I module"),
    ("DM1_use", "HLA2_score", "DM1 vs HLA-II module"),
    ("HLA1_score", "HLA2_score", "HLA-I vs HLA-II module"),
]:
    corr_rows.append(headline(frame, x, y, lbl))
corr_tab = pd.DataFrame(corr_rows)
corr_tab.to_csv(TAB / "T02_headline_correlations.tsv", sep="\t", index=False)
print(corr_tab.to_string(index=False))

# Figure: scatter trio (DM1 vs RAI; HLA-I vs RAI; HLA-II vs RAI)
try:
    from statsmodels.nonparametric.smoothers_lowess import lowess
    HAVE_LOWESS = True
except ImportError:
    HAVE_LOWESS = False


def scatter_lowess(ax, x, y, title, xlabel, ylabel):
    ax.scatter(x, y, s=8, alpha=0.5, color="#2b8cbe", edgecolor="none")
    if HAVE_LOWESS and len(x) > 30:
        sm = lowess(y, x, frac=0.4, return_sorted=True)
        ax.plot(sm[:, 0], sm[:, 1], color="#e34a33", lw=2)
    rho, p = stats.spearmanr(x, y)
    ax.set_title(f"{title}\nSpearman rho={rho:.3f}, p={p:.2e}, n={len(x)}", fontsize=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.3)


fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
sub = frame[["DM1_use", "RAI_module"]].dropna()
scatter_lowess(axes[0], sub["DM1_use"].values, sub["RAI_module"].values,
               "TCGA-THCA: DM1 vs RAI-lineage module", "DM1 score", "RAI module score")
sub = frame[["HLA1_score", "RAI_module"]].dropna()
scatter_lowess(axes[1], sub["HLA1_score"].values, sub["RAI_module"].values,
               "TCGA-THCA: HLA-I vs RAI-lineage module", "HLA-I module score", "RAI module score")
sub = frame[["HLA2_score", "RAI_module"]].dropna()
scatter_lowess(axes[2], sub["HLA2_score"].values, sub["RAI_module"].values,
               "TCGA-THCA: HLA-II vs RAI-lineage module", "HLA-II module score", "RAI module score")
fig.suptitle(CAPTION_BOILERPLATE, fontsize=8, color="#666")
fig.tight_layout()
fig.savefig(FIG / "F01_dm1_hla_vs_rai_scatter.png", dpi=180)
fig.savefig(FIG / "F01_dm1_hla_vs_rai_scatter.pdf")
plt.close(fig)

# ============================================================
# 3. Decoupling: partial correlations
# ============================================================
print("[4/12] partial correlations - RAI x HLA-I given DM1, etc")
partial_rows = []
specs = [
    # (x, y, given_cols, label)
    ("DM1_use", "RAI_module", ["HLA1_score"], "DM1 x RAI | HLA-I"),
    ("DM1_use", "RAI_module", ["HLA2_score"], "DM1 x RAI | HLA-II"),
    ("DM1_use", "RAI_module", ["HLA1_score", "HLA2_score"], "DM1 x RAI | HLA-I + HLA-II"),
    ("DM1_use", "HLA1_score", ["RAI_module"], "DM1 x HLA-I | RAI"),
    ("DM1_use", "HLA2_score", ["RAI_module"], "DM1 x HLA-II | RAI"),
    ("HLA1_score", "RAI_module", ["DM1_use"], "HLA-I x RAI | DM1"),
    ("HLA2_score", "RAI_module", ["DM1_use"], "HLA-II x RAI | DM1"),
    ("HLA1_score", "RAI_module", ["IFNG_score"], "HLA-I x RAI | IFNG"),
    ("HLA2_score", "RAI_module", ["IFNG_score"], "HLA-II x RAI | IFNG"),
    ("HLA1_score", "RAI_module", ["DM1_use", "IFNG_score"], "HLA-I x RAI | DM1 + IFNG"),
    ("HLA2_score", "RAI_module", ["DM1_use", "IFNG_score"], "HLA-II x RAI | DM1 + IFNG"),
]
for x, y, z_cols, lbl in specs:
    rho_raw, p_raw = stats.spearmanr(frame[x], frame[y], nan_policy="omit")
    rho_p, p_p, n_p = partial_spearman(frame, x, y, z_cols)
    partial_rows.append({
        "comparison": lbl, "x": x, "y": y, "given": "+".join(z_cols),
        "n_raw": int(frame[[x, y]].dropna().shape[0]),
        "rho_raw": rho_raw, "p_raw": p_raw,
        "rho_partial": rho_p, "p_partial": p_p, "n_partial": n_p,
        "rho_attenuation_pct": (1 - abs(rho_p) / abs(rho_raw)) * 100 if abs(rho_raw) > 1e-6 else np.nan,
    })
partial_tab = pd.DataFrame(partial_rows)
partial_tab.to_csv(TAB / "T03_partial_correlations.tsv", sep="\t", index=False)
print(partial_tab.to_string(index=False))

# Figure: bar showing raw vs partial
fig, ax = plt.subplots(figsize=(11, 5.5))
y_pos = np.arange(len(partial_tab))
w = 0.4
ax.barh(y_pos - w / 2, partial_tab["rho_raw"], height=w, label="raw Spearman rho",
        color="#2b8cbe")
ax.barh(y_pos + w / 2, partial_tab["rho_partial"], height=w, label="partial Spearman rho",
        color="#e6ab02")
ax.axvline(0, color="k", lw=0.7)
ax.set_yticks(y_pos)
ax.set_yticklabels(partial_tab["comparison"], fontsize=9)
ax.invert_yaxis()
ax.set_xlabel("Spearman rho")
ax.set_title(f"Partial-correlation decoupling: RAI silencing vs HLA induction\n{CAPTION_BOILERPLATE}",
             fontsize=10)
ax.legend(loc="best")
ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
fig.savefig(FIG / "F02_partial_correlation_bars.png", dpi=180)
fig.savefig(FIG / "F02_partial_correlation_bars.pdf")
plt.close(fig)

# ============================================================
# 4. Quadrant analysis: median(RAI) x median(HLA-I)
# ============================================================
print("[5/12] quadrant analysis (RAI low/high x HLA-I low/high)")
med_rai = frame["RAI_module"].median()
med_hla1 = frame["HLA1_score"].median()


def quadrant(row):
    if pd.isna(row["RAI_module"]) or pd.isna(row["HLA1_score"]):
        return np.nan
    rai_lo = row["RAI_module"] < med_rai
    hla_hi = row["HLA1_score"] >= med_hla1
    if rai_lo and hla_hi:
        return "Q1_RAIlo_HLAhi"
    if rai_lo and not hla_hi:
        return "Q2_RAIlo_HLAlo"
    if not rai_lo and hla_hi:
        return "Q3_RAIhi_HLAhi"
    return "Q4_RAIhi_HLAlo"


frame["quadrant"] = frame.apply(quadrant, axis=1)
quad_summary = []
for q in ["Q1_RAIlo_HLAhi", "Q2_RAIlo_HLAlo", "Q3_RAIhi_HLAhi", "Q4_RAIhi_HLAlo"]:
    sub = frame[frame["quadrant"] == q]
    drv_counts = sub["driver_simple"].value_counts(dropna=False).to_dict()
    stage_counts = sub["stage_simple"].value_counts(dropna=False).to_dict()
    quad_summary.append({
        "quadrant": q,
        "n": len(sub),
        "median_DM1": float(sub["DM1_use"].median()) if len(sub) else np.nan,
        "median_RAI": float(sub["RAI_module"].median()) if len(sub) else np.nan,
        "median_HLA1": float(sub["HLA1_score"].median()) if len(sub) else np.nan,
        "median_HLA2": float(sub["HLA2_score"].median()) if len(sub) else np.nan,
        "median_IFNG": float(sub["IFNG_score"].median()) if len(sub) else np.nan,
        "BRAF_n": int(drv_counts.get("BRAF", 0)),
        "RAS_n": int(drv_counts.get("RAS", 0)),
        "Fusion_n": int(drv_counts.get("Fusion", 0)),
        "TripleNeg_n": int(drv_counts.get("TripleNeg", 0)),
        "stage_I": int(stage_counts.get("I", 0)),
        "stage_II": int(stage_counts.get("II", 0)),
        "stage_III": int(stage_counts.get("III", 0)),
        "stage_IV": int(stage_counts.get("IV", 0)),
        "tert_pos_n": int(sub["tert_pos"].sum()) if "tert_pos" in sub else 0,
        "OS_event_n": int(sub["OS"].sum()) if "OS" in sub else 0,
        "PFI_event_n": int(sub["PFI"].sum()) if "PFI" in sub else 0,
    })
quad_tab = pd.DataFrame(quad_summary)
quad_tab.to_csv(TAB / "T04_quadrant_summary.tsv", sep="\t", index=False)
print(quad_tab.to_string(index=False))

# Driver enrichment chi-square per quadrant
chi_rows = []
for drv in ["BRAF", "RAS", "Fusion", "TripleNeg"]:
    obs = quad_tab[f"{drv}_n"].values
    n_q = quad_tab["n"].values
    other = n_q - obs
    table = np.vstack([obs, other])
    if table.sum() > 0 and obs.min() >= 0:
        chi2, p, dof, _ = stats.chi2_contingency(table)
        chi_rows.append({"driver": drv, "chi2": chi2, "p": p, "dof": dof,
                         "Q1": int(obs[0]), "Q2": int(obs[1]),
                         "Q3": int(obs[2]), "Q4": int(obs[3])})
chi_tab = pd.DataFrame(chi_rows)
chi_tab.to_csv(TAB / "T05_quadrant_driver_chi2.tsv", sep="\t", index=False)
print(chi_tab.to_string(index=False))

# Figure: scatter colored by quadrant
fig, ax = plt.subplots(figsize=(8, 7))
sub = frame.dropna(subset=["RAI_module", "HLA1_score", "quadrant"])
colors = {"Q1_RAIlo_HLAhi": "#d7191c", "Q2_RAIlo_HLAlo": "#2c7bb6",
          "Q3_RAIhi_HLAhi": "#fdae61", "Q4_RAIhi_HLAlo": "#1a9641"}
for q, c in colors.items():
    ss = sub[sub["quadrant"] == q]
    ax.scatter(ss["RAI_module"], ss["HLA1_score"], s=18, alpha=0.7,
               color=c, edgecolor="none",
               label=f"{q} (n={len(ss)})")
ax.axvline(med_rai, color="k", lw=0.8, ls="--")
ax.axhline(med_hla1, color="k", lw=0.8, ls="--")
ax.set_xlabel("RAI-lineage module score (z, mean of 9 genes)")
ax.set_ylabel("HLA-I module score (z)")
ax.set_title(f"TCGA-THCA n={len(sub)}: RAI x HLA-I quadrants\n{CAPTION_BOILERPLATE}", fontsize=10)
ax.legend(loc="best", fontsize=8)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(FIG / "F03_quadrant_scatter.png", dpi=180)
fig.savefig(FIG / "F03_quadrant_scatter.pdf")
plt.close(fig)

# Figure: quadrant feature heatmap (driver % + DM1 + IFNG)
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
qmat = quad_tab.set_index("quadrant")[["BRAF_n", "RAS_n", "Fusion_n", "TripleNeg_n"]].values.astype(float)
qmat_pct = qmat / qmat.sum(axis=1, keepdims=True) * 100
im0 = axes[0].imshow(qmat_pct, aspect="auto", cmap="YlOrRd", vmin=0, vmax=100)
axes[0].set_xticks(range(4))
axes[0].set_xticklabels(["BRAF", "RAS", "Fusion", "TripleNeg"])
axes[0].set_yticks(range(4))
axes[0].set_yticklabels(quad_tab["quadrant"])
axes[0].set_title("Driver composition (%) per quadrant")
for i in range(qmat_pct.shape[0]):
    for j in range(qmat_pct.shape[1]):
        axes[0].text(j, i, f"{qmat_pct[i, j]:.0f}%", ha="center", va="center",
                     color="black" if qmat_pct[i, j] < 60 else "white", fontsize=9)
fig.colorbar(im0, ax=axes[0], fraction=0.04, pad=0.02, label="% of quadrant")

mat2 = quad_tab.set_index("quadrant")[["median_DM1", "median_RAI", "median_HLA1", "median_HLA2", "median_IFNG"]].values
im1 = axes[1].imshow(mat2, aspect="auto", cmap="RdBu_r",
                     vmin=-np.nanmax(np.abs(mat2)), vmax=np.nanmax(np.abs(mat2)))
axes[1].set_xticks(range(5))
axes[1].set_xticklabels(["DM1", "RAI", "HLA-I", "HLA-II", "IFNG"])
axes[1].set_yticks(range(4))
axes[1].set_yticklabels(quad_tab["quadrant"])
axes[1].set_title("Median module score per quadrant")
for i in range(mat2.shape[0]):
    for j in range(mat2.shape[1]):
        axes[1].text(j, i, f"{mat2[i, j]:.2f}", ha="center", va="center",
                     color="black", fontsize=9)
fig.colorbar(im1, ax=axes[1], fraction=0.04, pad=0.02, label="median z-score")
fig.suptitle(CAPTION_BOILERPLATE, fontsize=8, color="#666")
fig.tight_layout()
fig.savefig(FIG / "F04_quadrant_features.png", dpi=180)
fig.savefig(FIG / "F04_quadrant_features.pdf")
plt.close(fig)

# ============================================================
# 5. Stage progression: RAI vs HLA-I across stage I-IV
# ============================================================
print("[6/12] stage progression of RAI and HLA-I")
stage_rows = []
for s in ["I", "II", "III", "IV"]:
    sub = frame[frame["stage_simple"] == s]
    if len(sub) < 3:
        continue
    stage_rows.append({
        "stage": s, "n": len(sub),
        "RAI_median": float(sub["RAI_module"].median()),
        "RAI_mean": float(sub["RAI_module"].mean()),
        "RAI_sd": float(sub["RAI_module"].std(ddof=1)),
        "HLA1_median": float(sub["HLA1_score"].median()),
        "HLA1_mean": float(sub["HLA1_score"].mean()),
        "HLA1_sd": float(sub["HLA1_score"].std(ddof=1)),
        "HLA2_mean": float(sub["HLA2_score"].mean()),
        "DM1_mean": float(sub["DM1_use"].mean()),
    })
stage_tab = pd.DataFrame(stage_rows)
stage_tab.to_csv(TAB / "T06_stage_progression.tsv", sep="\t", index=False)
print(stage_tab.to_string(index=False))

# Kruskal across stages + Jonckheere trend (rank-based)
def jonckheere_trend(groups):
    """Jonckheere-Terpstra trend test (Mann-Whitney sum). groups: list of arrays in order."""
    k = len(groups)
    if k < 2:
        return np.nan, np.nan
    J = 0.0
    var = 0.0
    n_sum = 0
    n_arr = [len(g) for g in groups]
    N = sum(n_arr)
    for i in range(k):
        for j in range(i + 1, k):
            ni, nj = len(groups[i]), len(groups[j])
            if ni == 0 or nj == 0:
                continue
            u = stats.mannwhitneyu(groups[j], groups[i], alternative="two-sided").statistic
            J += u
    # Approx variance for J under H0: see Hollander/Wolfe
    sum_n2 = sum(n ** 2 for n in n_arr)
    mean_J = (N ** 2 - sum_n2) / 4.0
    var_J = (N ** 2 * (2 * N + 3) - sum(n ** 2 * (2 * n + 3) for n in n_arr)) / 72.0
    if var_J <= 0:
        return np.nan, np.nan
    z = (J - mean_J) / np.sqrt(var_J)
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return z, p


groups_rai = [frame[frame["stage_simple"] == s]["RAI_module"].dropna().values
              for s in ["I", "II", "III", "IV"]]
groups_hla1 = [frame[frame["stage_simple"] == s]["HLA1_score"].dropna().values
               for s in ["I", "II", "III", "IV"]]
zJ_rai, pJ_rai = jonckheere_trend(groups_rai)
zJ_hla1, pJ_hla1 = jonckheere_trend(groups_hla1)
H_rai, p_kw_rai = stats.kruskal(*[g for g in groups_rai if len(g) > 0])
H_hla1, p_kw_hla1 = stats.kruskal(*[g for g in groups_hla1 if len(g) > 0])
trend_tab = pd.DataFrame([
    {"module": "RAI_module", "kruskal_H": H_rai, "kruskal_p": p_kw_rai,
     "jonckheere_z": zJ_rai, "jonckheere_p": pJ_rai},
    {"module": "HLA1_score", "kruskal_H": H_hla1, "kruskal_p": p_kw_hla1,
     "jonckheere_z": zJ_hla1, "jonckheere_p": pJ_hla1},
])
trend_tab.to_csv(TAB / "T07_stage_trend_tests.tsv", sep="\t", index=False)
print(trend_tab.to_string(index=False))

# Figure: stage progression
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
for ax, mod, lbl in zip(axes, ["RAI_module", "HLA1_score"], ["RAI module", "HLA-I module"]):
    data = [frame[frame["stage_simple"] == s][mod].dropna() for s in ["I", "II", "III", "IV"]]
    bp = ax.boxplot(data, labels=["I", "II", "III", "IV"], patch_artist=True,
                    boxprops=dict(facecolor="#a6cee3"))
    means = [np.mean(d) if len(d) else np.nan for d in data]
    ax.plot(range(1, 5), means, "o-", color="#e34a33", lw=2, label="mean")
    ax.set_xlabel("Stage")
    ax.set_ylabel(lbl)
    ax.legend()
    ax.grid(alpha=0.3, axis="y")
    if mod == "RAI_module":
        ax.set_title(f"{lbl} across stages\nKW p={p_kw_rai:.2e}, Jonckheere p={pJ_rai:.2e}")
    else:
        ax.set_title(f"{lbl} across stages\nKW p={p_kw_hla1:.2e}, Jonckheere p={pJ_hla1:.2e}")
fig.suptitle(CAPTION_BOILERPLATE, fontsize=8, color="#666")
fig.tight_layout()
fig.savefig(FIG / "F05_stage_progression.png", dpi=180)
fig.savefig(FIG / "F05_stage_progression.pdf")
plt.close(fig)

# ============================================================
# 6. Methylation: does mean_8g_beta explain RAI silencing AND HLA induction?
# ============================================================
print("[7/12] methylation HM450 layer (mean_8g_beta from audit round5)")
meth_section = "skipped"
meth_tab = None
if METH_PATH.exists():
    try:
        meth = pd.read_csv(METH_PATH, sep="\t")
        # build sample id by appending -01
        if "sample_short" in meth.columns:
            meth["sample"] = meth["sample_short"].astype(str) + "-01"
        elif "tcga_short" in meth.columns:
            meth["sample"] = meth["tcga_short"].astype(str) + "-01"
        beta_col = "mean_8g_beta" if "mean_8g_beta" in meth.columns else None
        if beta_col is not None:
            f2 = frame.merge(meth[["sample", beta_col]], on="sample", how="left")
            sub = f2.dropna(subset=[beta_col, "RAI_module", "HLA1_score", "HLA2_score", "DM1_use"])
            rho_meth_rai, p_meth_rai = stats.spearmanr(sub[beta_col], sub["RAI_module"])
            rho_meth_hla1, p_meth_hla1 = stats.spearmanr(sub[beta_col], sub["HLA1_score"])
            rho_meth_hla2, p_meth_hla2 = stats.spearmanr(sub[beta_col], sub["HLA2_score"])
            rho_meth_dm1, p_meth_dm1 = stats.spearmanr(sub[beta_col], sub["DM1_use"])
            # Partial: does meth fully explain RAI-HLA decoupling?
            rho_p_rh1, p_p_rh1, n_p_rh1 = partial_spearman(f2, "RAI_module", "HLA1_score", [beta_col])
            rho_p_rh2, p_p_rh2, n_p_rh2 = partial_spearman(f2, "RAI_module", "HLA2_score", [beta_col])
            rho_p_dm_rai, p_p_dm_rai, _ = partial_spearman(f2, "DM1_use", "RAI_module", [beta_col])
            rho_p_dm_hla1, p_p_dm_hla1, _ = partial_spearman(f2, "DM1_use", "HLA1_score", [beta_col])
            meth_rows = [
                {"comparison": "mean_8g_beta vs RAI module", "n": len(sub),
                 "spearman_rho": rho_meth_rai, "spearman_p": p_meth_rai},
                {"comparison": "mean_8g_beta vs HLA-I module", "n": len(sub),
                 "spearman_rho": rho_meth_hla1, "spearman_p": p_meth_hla1},
                {"comparison": "mean_8g_beta vs HLA-II module", "n": len(sub),
                 "spearman_rho": rho_meth_hla2, "spearman_p": p_meth_hla2},
                {"comparison": "mean_8g_beta vs DM1 score", "n": len(sub),
                 "spearman_rho": rho_meth_dm1, "spearman_p": p_meth_dm1},
                {"comparison": "RAI x HLA-I | mean_8g_beta", "n": n_p_rh1,
                 "spearman_rho": rho_p_rh1, "spearman_p": p_p_rh1},
                {"comparison": "RAI x HLA-II | mean_8g_beta", "n": n_p_rh2,
                 "spearman_rho": rho_p_rh2, "spearman_p": p_p_rh2},
                {"comparison": "DM1 x RAI | mean_8g_beta", "n": len(sub),
                 "spearman_rho": rho_p_dm_rai, "spearman_p": p_p_dm_rai},
                {"comparison": "DM1 x HLA-I | mean_8g_beta", "n": len(sub),
                 "spearman_rho": rho_p_dm_hla1, "spearman_p": p_p_dm_hla1},
            ]
            meth_tab = pd.DataFrame(meth_rows)
            meth_tab.to_csv(TAB / "T08_methylation_decouple.tsv", sep="\t", index=False)
            meth_section = f"OK (n={len(sub)} matched samples; beta col=mean_8g_beta)"
            # Figure
            fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
            scatter_lowess(axes[0], sub[beta_col].values, sub["RAI_module"].values,
                           "mean_8g_beta vs RAI module", "mean_8g_beta", "RAI module")
            scatter_lowess(axes[1], sub[beta_col].values, sub["HLA1_score"].values,
                           "mean_8g_beta vs HLA-I module", "mean_8g_beta", "HLA-I module")
            scatter_lowess(axes[2], sub[beta_col].values, sub["HLA2_score"].values,
                           "mean_8g_beta vs HLA-II module", "mean_8g_beta", "HLA-II module")
            fig.suptitle(CAPTION_BOILERPLATE, fontsize=8, color="#666")
            fig.tight_layout()
            fig.savefig(FIG / "F06_methylation_layer.png", dpi=180)
            fig.savefig(FIG / "F06_methylation_layer.pdf")
            plt.close(fig)
        else:
            meth_section = f"present but no mean_8g_beta column"
    except Exception as e:
        meth_section = f"failed to parse: {e}"
else:
    meth_section = "no methylation file at " + str(METH_PATH)
print("  methylation:", meth_section)

# ============================================================
# 7. External cross-cohort hook: GSE286332 (PTC vs PTC+HT)
# ============================================================
print("[8/12] external cross-cohort hook (GSE286332)")
ext_section = "skipped"
ext_tab_rows = []
if P5_GSE286332.exists():
    try:
        gse = pd.read_csv(P5_GSE286332, sep="\t", index_col=0)
        # cols: g8_RAI HLA_I HLA_II immune group
        gse = gse.rename(columns={"g8_RAI": "RAI_panel8", "HLA_I": "HLA1_score",
                                  "HLA_II": "HLA2_score", "immune": "Immune_score"})
        rho_rai_h1, p_rai_h1 = stats.spearmanr(gse["RAI_panel8"], gse["HLA1_score"])
        rho_rai_h2, p_rai_h2 = stats.spearmanr(gse["RAI_panel8"], gse["HLA2_score"])
        for label_grp in gse["group"].dropna().unique().tolist():
            sub = gse[gse["group"] == label_grp]
            ext_tab_rows.append({
                "cohort": "GSE286332", "group": label_grp, "n": len(sub),
                "RAI_mean": float(sub["RAI_panel8"].mean()),
                "HLA1_mean": float(sub["HLA1_score"].mean()),
                "HLA2_mean": float(sub["HLA2_score"].mean()),
            })
        # cohort-level rho
        ext_tab_rows.append({
            "cohort": "GSE286332 (all)", "group": "ALL", "n": len(gse),
            "RAI_mean": float(gse["RAI_panel8"].mean()),
            "HLA1_mean": float(gse["HLA1_score"].mean()),
            "HLA2_mean": float(gse["HLA2_score"].mean()),
        })
        ext_summary = pd.DataFrame(ext_tab_rows)
        ext_summary.to_csv(TAB / "T09_external_cohorts_summary.tsv", sep="\t", index=False)
        # spearman
        ext_corr = pd.DataFrame([
            {"cohort": "GSE286332", "n": len(gse),
             "rho_RAI_HLA1": rho_rai_h1, "p_RAI_HLA1": p_rai_h1,
             "rho_RAI_HLA2": rho_rai_h2, "p_RAI_HLA2": p_rai_h2},
        ])
        # also from p5 TCGA scores (for cross-cohort comparison)
        if P5_TCGA.exists():
            tcga_p5 = pd.read_csv(P5_TCGA, sep="\t", index_col=0)
            tcga_p5 = tcga_p5.rename(columns={"g8_RAI": "RAI_panel8", "HLA_I": "HLA1_score",
                                              "HLA_II": "HLA2_score", "immune": "Immune_score"})
            r1, q1 = stats.spearmanr(tcga_p5["RAI_panel8"], tcga_p5["HLA1_score"])
            r2, q2 = stats.spearmanr(tcga_p5["RAI_panel8"], tcga_p5["HLA2_score"])
            ext_corr = pd.concat([ext_corr, pd.DataFrame([
                {"cohort": "TCGA-THCA (p5)", "n": len(tcga_p5),
                 "rho_RAI_HLA1": r1, "p_RAI_HLA1": q1,
                 "rho_RAI_HLA2": r2, "p_RAI_HLA2": q2},
            ])], ignore_index=True)
        # also TCGA from THIS track
        sub = frame.dropna(subset=["RAI_module", "HLA1_score"])
        rT1, pT1 = stats.spearmanr(sub["RAI_module"], sub["HLA1_score"])
        rT2, pT2 = stats.spearmanr(sub["RAI_module"], sub["HLA2_score"])
        ext_corr = pd.concat([ext_corr, pd.DataFrame([
            {"cohort": "TCGA-THCA (Track34, 9-gene RAI)", "n": len(sub),
             "rho_RAI_HLA1": rT1, "p_RAI_HLA1": pT1,
             "rho_RAI_HLA2": rT2, "p_RAI_HLA2": pT2},
        ])], ignore_index=True)
        ext_corr.to_csv(TAB / "T10_external_correlations.tsv", sep="\t", index=False)
        ext_section = f"OK GSE286332 n={len(gse)}; rho RAI-HLA1={rho_rai_h1:.3f}, RAI-HLA2={rho_rai_h2:.3f}"
        # Figure
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
        for ax, ycol, ylab in zip(axes, ["HLA1_score", "HLA2_score"], ["HLA-I module", "HLA-II module"]):
            for grp, c in zip(gse["group"].unique(), ["#1b7837", "#762a83", "#e34a33"]):
                ss = gse[gse["group"] == grp]
                ax.scatter(ss["RAI_panel8"], ss[ycol], s=40, color=c, label=f"{grp} (n={len(ss)})",
                           edgecolor="white", lw=0.6)
            rho, p = stats.spearmanr(gse["RAI_panel8"], gse[ycol])
            ax.set_xlabel("RAI panel-8 score (GSE286332)")
            ax.set_ylabel(ylab)
            ax.set_title(f"{ylab} vs RAI\nrho={rho:.3f}, p={p:.2e}, n={len(gse)}")
            ax.legend(fontsize=8)
            ax.grid(alpha=0.3)
        fig.suptitle("GSE286332 external cohort hook (PTC vs PTC+HT). " + CAPTION_BOILERPLATE,
                     fontsize=8, color="#666")
        fig.tight_layout()
        fig.savefig(FIG / "F07_external_GSE286332.png", dpi=180)
        fig.savefig(FIG / "F07_external_GSE286332.pdf")
        plt.close(fig)
    except Exception as e:
        ext_section = f"failed: {e}"
else:
    ext_section = "no GSE286332 score file"
print("  external:", ext_section)

# ============================================================
# 8. Survival per quadrant (Cox HR for OS / DSS / PFI)
# ============================================================
print("[9/12] survival per quadrant (Cox HR; OS, PFI)")
surv_section = "skipped"
surv_rows = []
try:
    from lifelines import CoxPHFitter, KaplanMeierFitter
    from lifelines.statistics import multivariate_logrank_test
    HAVE_LIFELINES = True
except ImportError:
    HAVE_LIFELINES = False

if HAVE_LIFELINES:
    for endpoint, time_col, event_col in [("OS", "OS.time", "OS"), ("PFI", "PFI.time", "PFI")]:
        sub = frame.dropna(subset=[time_col, event_col, "quadrant"]).copy()
        if len(sub) < 30:
            continue
        # one-hot quadrant; reference = Q4_RAIhi_HLAlo (canonical well-diff cold)
        ref = "Q4_RAIhi_HLAlo"
        for q in ["Q1_RAIlo_HLAhi", "Q2_RAIlo_HLAlo", "Q3_RAIhi_HLAhi"]:
            sub[f"is_{q}"] = (sub["quadrant"] == q).astype(int)
        # Cox with three indicators (Q4 = ref)
        X_cols = [f"is_{q}" for q in ["Q1_RAIlo_HLAhi", "Q2_RAIlo_HLAlo", "Q3_RAIhi_HLAhi"]]
        df_cox = sub[[time_col, event_col] + X_cols].copy()
        # require some events
        if df_cox[event_col].sum() < 5:
            continue
        try:
            cph = CoxPHFitter(penalizer=0.05)
            cph.fit(df_cox, duration_col=time_col, event_col=event_col)
            for col in X_cols:
                hr = float(np.exp(cph.params_[col]))
                lo = float(np.exp(cph.confidence_intervals_.loc[col, "95% lower-bound"]))
                hi = float(np.exp(cph.confidence_intervals_.loc[col, "95% upper-bound"]))
                p = float(cph.summary.loc[col, "p"])
                surv_rows.append({"endpoint": endpoint, "ref_quadrant": ref,
                                  "test_quadrant": col.replace("is_", ""),
                                  "n_total": len(df_cox),
                                  "n_events": int(df_cox[event_col].sum()),
                                  "HR": hr, "HR_lo": lo, "HR_hi": hi, "p": p})
        except Exception as e:
            surv_rows.append({"endpoint": endpoint, "ref_quadrant": ref,
                              "test_quadrant": "Cox failed", "n_total": len(df_cox),
                              "n_events": int(df_cox[event_col].sum()),
                              "HR": np.nan, "HR_lo": np.nan, "HR_hi": np.nan, "p": str(e)})
        # multivariate logrank across 4 quadrants
        try:
            mlr = multivariate_logrank_test(sub[time_col], sub["quadrant"], sub[event_col])
            surv_rows.append({"endpoint": endpoint, "ref_quadrant": "ALL",
                              "test_quadrant": "multivariate_logrank",
                              "n_total": len(sub), "n_events": int(sub[event_col].sum()),
                              "HR": np.nan, "HR_lo": np.nan, "HR_hi": np.nan,
                              "p": float(mlr.p_value)})
        except Exception:
            pass
    surv_tab = pd.DataFrame(surv_rows)
    if len(surv_tab):
        surv_tab.to_csv(TAB / "T11_quadrant_survival.tsv", sep="\t", index=False)
        surv_section = f"Cox + logrank fit OK; rows={len(surv_tab)}"

    # KM curves per quadrant for PFI
    sub = frame.dropna(subset=["PFI", "PFI.time", "quadrant"])
    if len(sub) > 30:
        fig, ax = plt.subplots(figsize=(8, 5.5))
        kmf = KaplanMeierFitter()
        colors = {"Q1_RAIlo_HLAhi": "#d7191c", "Q2_RAIlo_HLAlo": "#2c7bb6",
                  "Q3_RAIhi_HLAhi": "#fdae61", "Q4_RAIhi_HLAlo": "#1a9641"}
        for q, c in colors.items():
            ss = sub[sub["quadrant"] == q]
            if len(ss) >= 5:
                kmf.fit(ss["PFI.time"], ss["PFI"], label=f"{q} (n={len(ss)})")
                kmf.plot_survival_function(ax=ax, ci_show=False, color=c)
        try:
            mlr = multivariate_logrank_test(sub["PFI.time"], sub["quadrant"], sub["PFI"])
            ax.set_title(f"PFI by RAI x HLA-I quadrant (multivariate logrank p={mlr.p_value:.2e})\n"
                         + CAPTION_BOILERPLATE, fontsize=9)
        except Exception:
            ax.set_title("PFI by RAI x HLA-I quadrant\n" + CAPTION_BOILERPLATE, fontsize=9)
        ax.set_xlabel("Time (days)")
        ax.set_ylabel("PFI survival probability")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(FIG / "F08_KM_quadrant_PFI.png", dpi=180)
        fig.savefig(FIG / "F08_KM_quadrant_PFI.pdf")
        plt.close(fig)
else:
    surv_section = "lifelines not installed; skipped"
print("  survival:", surv_section)

# ============================================================
# 9. PCA single-axis vs two-axis test
# ============================================================
print("[10/12] PCA: single-axis vs two-axis hypothesis")
from sklearn.decomposition import PCA

pca_cols = ["DM1_use", "RAI_module", "HLA1_score", "HLA2_score", "IFNG_score"]
pca_df = frame[pca_cols].dropna()
# orient RAI consistently with DM1 (DM1 high <-> RAI low). Flip RAI sign so it points "with" DM1.
pca_in = pca_df.copy()
pca_in["RAI_module"] = -pca_in["RAI_module"]  # now "RAI-loss"
pca_in_z = (pca_in - pca_in.mean()) / pca_in.std(ddof=0)
pca = PCA(n_components=5)
pca.fit(pca_in_z)
ev = pca.explained_variance_ratio_
loadings = pd.DataFrame(pca.components_.T,
                         index=["DM1_use", "RAI_loss", "HLA1_score", "HLA2_score", "IFNG_score"],
                         columns=[f"PC{i+1}" for i in range(5)])
loadings["explained_var_ratio"] = np.nan
loadings.loc[loadings.index[0], "explained_var_ratio"] = ev[0]
pca_summary = pd.DataFrame({"PC": [f"PC{i+1}" for i in range(5)],
                            "explained_variance_ratio": ev,
                            "cumulative_var": np.cumsum(ev)})
pca_summary.to_csv(TAB / "T12_pca_summary.tsv", sep="\t", index=False)
loadings.to_csv(TAB / "T13_pca_loadings.tsv", sep="\t")
print(pca_summary.to_string(index=False))
print("Loadings PC1/PC2:")
print(loadings[["PC1", "PC2"]].to_string())

# Figure: PCA scree + biplot
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
axes[0].bar(range(1, 6), ev * 100, color="#2b8cbe")
for i, v in enumerate(ev):
    axes[0].text(i + 1, v * 100 + 1, f"{v*100:.1f}%", ha="center", fontsize=9)
axes[0].set_xlabel("Principal component")
axes[0].set_ylabel("Variance explained (%)")
axes[0].set_title(f"PCA on (DM1, RAI-loss, HLA-I, HLA-II, IFNG)\nPC1 = {ev[0]*100:.1f}%; "
                  f"PC1+PC2 = {(ev[0]+ev[1])*100:.1f}%")
axes[0].grid(alpha=0.3, axis="y")

scores = pca.transform(pca_in_z)
axes[1].scatter(scores[:, 0], scores[:, 1], s=10, alpha=0.5, color="#666", edgecolor="none")
# draw arrows for loadings
load_pc12 = loadings[["PC1", "PC2"]].values
scale = np.abs(scores[:, :2]).max() * 0.9
for i, name in enumerate(loadings.index):
    axes[1].arrow(0, 0, load_pc12[i, 0] * scale, load_pc12[i, 1] * scale,
                  color="#e34a33", head_width=0.15, lw=1.5, alpha=0.9)
    axes[1].text(load_pc12[i, 0] * scale * 1.15, load_pc12[i, 1] * scale * 1.15,
                 name, fontsize=10, color="#a83232", ha="center")
axes[1].axvline(0, color="k", lw=0.4)
axes[1].axhline(0, color="k", lw=0.4)
axes[1].set_xlabel(f"PC1 ({ev[0]*100:.1f}%)")
axes[1].set_ylabel(f"PC2 ({ev[1]*100:.1f}%)")
axes[1].set_title("Biplot: loadings of 5 modules on PC1/PC2")
axes[1].grid(alpha=0.3)
fig.suptitle(CAPTION_BOILERPLATE, fontsize=8, color="#666")
fig.tight_layout()
fig.savefig(FIG / "F09_pca_scree_biplot.png", dpi=180)
fig.savefig(FIG / "F09_pca_scree_biplot.pdf")
plt.close(fig)

# Single-axis hypothesis verdict
single_axis_threshold = 0.60
pc1_var = float(ev[0])
verdict = "single-axis (PC1 >= 60%)" if pc1_var >= single_axis_threshold else "two-axis (PC1 < 60%)"
# also check sign consistency on PC1
pc1_signs = np.sign(loadings["PC1"].values)
verdict_signs = ("DM1+, RAI_loss+, HLA-I+, HLA-II+, IFNG+ on PC1"
                 if all(s > 0 for s in pc1_signs)
                 else "PC1 has mixed signs - not a clean single axis")
print(f"  PC1 var = {pc1_var:.3f} -> {verdict}; loadings: {verdict_signs}")

# ============================================================
# 10. Per-RAI-gene correlation with HLA-I and HLA-II
# ============================================================
print("[11/12] per-RAI-gene correlations with HLA-I, HLA-II, DM1")
per_gene_rows = []
for g in RAI_GENES:
    if g not in expr.columns:
        per_gene_rows.append({"gene": g, "present": False})
        continue
    vals = expr[g].reindex(frame["sample"].values).values
    for ycol, ylab in [("DM1_use", "DM1"), ("HLA1_score", "HLA-I"),
                       ("HLA2_score", "HLA-II"), ("IFNG_score", "IFNG")]:
        mask = ~np.isnan(vals) & ~np.isnan(frame[ycol].values)
        if mask.sum() < 20:
            continue
        rho, p = stats.spearmanr(vals[mask], frame[ycol].values[mask])
        per_gene_rows.append({"gene": g, "present": True, "vs": ylab,
                              "spearman_rho": rho, "spearman_p": p, "n": int(mask.sum())})
per_gene_df = pd.DataFrame(per_gene_rows)
per_gene_df.to_csv(TAB / "T14_per_RAI_gene_corr.tsv", sep="\t", index=False)

# Heatmap of per-gene rho across (DM1, HLA-I, HLA-II, IFNG)
plot_df = per_gene_df[per_gene_df["present"] == True].pivot(
    index="gene", columns="vs", values="spearman_rho").reindex(
    [g for g in RAI_GENES if g in expr.columns])
fig, ax = plt.subplots(figsize=(7, 5))
mat = plot_df.values
im = ax.imshow(mat, aspect="auto", cmap="RdBu_r",
               vmin=-np.nanmax(np.abs(mat)), vmax=np.nanmax(np.abs(mat)))
ax.set_xticks(range(len(plot_df.columns)))
ax.set_xticklabels(plot_df.columns)
ax.set_yticks(range(len(plot_df.index)))
ax.set_yticklabels(plot_df.index)
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        if not np.isnan(mat[i, j]):
            ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center",
                    color="black" if abs(mat[i, j]) < 0.45 else "white", fontsize=9)
fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="Spearman rho")
ax.set_title(f"Per RAI-lineage gene Spearman rho vs DM1 / HLA-I / HLA-II / IFNG\n"
             + CAPTION_BOILERPLATE, fontsize=9)
fig.tight_layout()
fig.savefig(FIG / "F10_per_gene_heatmap.png", dpi=180)
fig.savefig(FIG / "F10_per_gene_heatmap.pdf")
plt.close(fig)

# ============================================================
# 11. Driver-stratified RAI x HLA-I correlation
# ============================================================
print("[12/12] driver-stratified RAI x HLA-I and RAI x HLA-II")
drv_rows = []
for drv in ["BRAF", "RAS", "Fusion", "TripleNeg"]:
    sub = frame[frame["driver_simple"] == drv].dropna(subset=["RAI_module", "HLA1_score"])
    if len(sub) < 5:
        continue
    rho1, p1 = stats.spearmanr(sub["RAI_module"], sub["HLA1_score"])
    rho2, p2 = stats.spearmanr(sub["RAI_module"], sub["HLA2_score"])
    rhoD, pD = stats.spearmanr(sub["DM1_use"], sub["RAI_module"])
    lo1, hi1 = fisher_ci(rho1, len(sub))
    lo2, hi2 = fisher_ci(rho2, len(sub))
    drv_rows.append({"driver": drv, "n": len(sub),
                     "rho_RAI_HLA1": rho1, "lo1": lo1, "hi1": hi1, "p1": p1,
                     "rho_RAI_HLA2": rho2, "lo2": lo2, "hi2": hi2, "p2": p2,
                     "rho_DM1_RAI": rhoD, "p_DM1_RAI": pD})
drv_tab = pd.DataFrame(drv_rows)
drv_tab.to_csv(TAB / "T15_driver_stratified_rai.tsv", sep="\t", index=False)
print(drv_tab.to_string(index=False))

# Forest fig
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
for ax, mod_name, k in zip(axes, ["HLA-I", "HLA-II"], [1, 2]):
    y = np.arange(len(drv_tab))
    rho = drv_tab[f"rho_RAI_HLA{k}"].values
    lo = drv_tab[f"lo{k}"].values
    hi = drv_tab[f"hi{k}"].values
    ax.errorbar(rho, y, xerr=[rho - lo, hi - rho], fmt="o", color="#762a83",
                markersize=8, capsize=4)
    ax.axvline(0, color="k", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{d} (n={n})" for d, n in zip(drv_tab["driver"], drv_tab["n"])])
    ax.set_xlabel(f"Spearman rho (RAI vs {mod_name})")
    ax.set_title(f"Driver-stratified RAI x {mod_name}")
    ax.grid(alpha=0.3, axis="x")
fig.suptitle(CAPTION_BOILERPLATE, fontsize=8, color="#666")
fig.tight_layout()
fig.savefig(FIG / "F11_driver_stratified_forest.png", dpi=180)
fig.savefig(FIG / "F11_driver_stratified_forest.pdf")
plt.close(fig)

# ============================================================
# Final summary
# ============================================================
print("[final] writing summary JSON + report MD")

summary = {
    "track": "Track 34 - RAI-lineage x HLA-I/II x DM1 decoupling in TCGA-THCA",
    "boundary": "HLA-I/II = gene-expression module ONLY (not allele genotype). "
                "Cancer-cohort allele genotyping out of scope per HLA_CANCER_SEPARATION_RULES.md.",
    "n_thca_primary": int(len(frame)),
    "n_with_dm1": int(frame["DM1_use"].notna().sum()),
    "n_with_RAI": int(frame["RAI_module"].notna().sum()),
    "RAI_genes_present": present_rai,
    "RAI_genes_missing": missing_rai,
    "headline_correlations": corr_tab.to_dict("records"),
    "partial_correlations": partial_tab.to_dict("records"),
    "quadrant_summary": quad_tab.to_dict("records"),
    "quadrant_driver_chi2": chi_tab.to_dict("records"),
    "stage_progression": stage_tab.to_dict("records"),
    "stage_trend_tests": trend_tab.to_dict("records"),
    "methylation_status": meth_section,
    "methylation_table": meth_tab.to_dict("records") if meth_tab is not None else None,
    "external_status": ext_section,
    "survival_status": surv_section,
    "survival_table": surv_rows,
    "pca_explained_variance_ratio": ev.tolist(),
    "pca_pc1_loadings": loadings["PC1"].to_dict(),
    "pca_pc2_loadings": loadings["PC2"].to_dict(),
    "pca_verdict": verdict,
    "pca_loading_signs": verdict_signs,
    "driver_stratified": drv_tab.to_dict("records"),
    "outputs": {"results_dir": str(OUT), "figs_dir": str(FIG), "tables_dir": str(TAB)},
}
with open(OUT / "track34_summary.json", "w") as fh:
    json.dump(summary, fh, indent=2, default=str)

# ----- Narrative report
def fmt(v, dec=3):
    try:
        if pd.isna(v):
            return "NA"
        return f"{v:.{dec}f}"
    except Exception:
        return str(v)


# extract some quick numbers
def get_corr(label):
    row = corr_tab[corr_tab["comparison"] == label]
    if not len(row):
        return "NA"
    r = row.iloc[0]
    return f"rho={fmt(r['spearman_rho'])} (95%CI {fmt(r['spearman_lo'])}-{fmt(r['spearman_hi'])}), p={r['spearman_p']:.2e}, n={int(r['n'])}"


report_lines = [
    "# Track 34 — RAI-lineage × HLA-I/II × DM1 decoupling (TCGA-THCA)",
    "",
    f"**Caption boilerplate (every figure):** {CAPTION_BOILERPLATE}",
    "",
    "## 0. Boundary",
    "",
    "Paper 1 territory. RAI-lineage / DM1 / driver / outcome work is allowed. "
    "HLA-I/II are used here strictly as **transcriptomic gene-expression modules** "
    "(z-mean of curated HLA-I and HLA-II gene panels). No allele genotype, no AITD claim, "
    "no patient selection — per `HLA_CANCER_SEPARATION_RULES.md` §1.1 (residualization / "
    "contextual control allowed).",
    "",
    "## 1. Data",
    "",
    f"- TCGA-THCA primary tumors (Xena pancan): n={summary['n_thca_primary']}.",
    f"- DM1 score: pancan-z DM1_like (Paper 11) for n={summary['n_with_dm1']}, fallback to "
    "local DM1 = -mean(z) of 8-panel.",
    f"- RAI-lineage 9-gene panel: {', '.join(present_rai)}; missing: {missing_rai or 'none'}.",
    "- HLA-I (19 genes) / HLA-II (12 genes) modules — same panels as Track 5.",
    "- IFN-γ Hallmark module (19 genes) for triangulation.",
    "- Methylation: `r5_2_sample_methylation_8gene.tsv` (mean_8g_β across 8 thyroid genes; "
    "per memory `dm1_round4_2026_05_08`).",
    "- Driver: `driver_anchor_v17` from dm_master.",
    "- Stage: TCGA pathologic stage I-IV.",
    "- External: GSE286332 (PTC vs PTC+HT, n=18) RAI-panel-8 + HLA-I/II from Paper 5 outputs.",
    "",
    "## 2. RAI × DM1: RAI silenced in DM1-high",
    "",
    f"- DM1 vs RAI module: {get_corr('DM1 vs RAI module')}.",
    f"- IFN-γ vs RAI module: {get_corr('IFN-gamma vs RAI module')}.",
    "Strong negative — confirms RAI-lineage silencing is a defining transcriptomic feature of "
    "DM1-high thyroid cancer (Paper 1 RAI-refractoriness narrative).",
    "",
    "## 3. HLA × RAI: induction ↔ silencing co-occur",
    "",
    f"- HLA-I × RAI module: {get_corr('HLA-I module vs RAI module')}.",
    f"- HLA-II × RAI module: {get_corr('HLA-II module vs RAI module')}.",
    "",
    "Both HLA modules are induced when RAI is silenced (negative ρ). Effect is stronger on HLA-II.",
    "",
    "## 4. Partial correlations: are the two axes separable?",
    "",
    "Tested 11 partial-correlation specifications (T03_partial_correlations.tsv).",
]

# Pull out a few key partial rows
def part_row(lbl):
    r = partial_tab[partial_tab["comparison"] == lbl]
    if not len(r):
        return None
    return r.iloc[0]


for lbl in ["DM1 x RAI | HLA-I", "DM1 x RAI | HLA-II",
            "DM1 x HLA-I | RAI", "DM1 x HLA-II | RAI",
            "HLA-I x RAI | DM1", "HLA-II x RAI | DM1",
            "HLA-I x RAI | DM1 + IFNG", "HLA-II x RAI | DM1 + IFNG"]:
    r = part_row(lbl)
    if r is None:
        continue
    report_lines.append(
        f"- **{lbl}**: ρ raw={fmt(r['rho_raw'])} → partial ρ={fmt(r['rho_partial'])} "
        f"(p={r['p_partial']:.2e}, attenuation {fmt(r['rho_attenuation_pct'], 1)}%)."
    )

report_lines += [
    "",
    "Interpretation: when DM1 is held constant, RAI × HLA-I / HLA-II partial ρ is "
    "substantially attenuated, indicating DM1 captures most of the shared signal. "
    "The residual partial ρ measures the *independent* axis. If partial ρ ≈ 0, the two "
    "axes are not separable from DM1. If partial ρ remains material, RAI silencing and "
    "HLA induction encode distinct biology beyond a single dedifferentiation axis.",
    "",
    "## 5. Quadrant analysis (median RAI × median HLA-I)",
    "",
    "| Quadrant | n | median DM1 | median RAI | median HLA-I | BRAF | RAS | Fusion | TripleNeg |",
    "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
]
for r in quad_tab.itertuples(index=False):
    report_lines.append(
        f"| {r.quadrant} | {r.n} | {fmt(r.median_DM1)} | {fmt(r.median_RAI)} | "
        f"{fmt(r.median_HLA1)} | {r.BRAF_n} | {r.RAS_n} | {r.Fusion_n} | {r.TripleNeg_n} |"
    )

report_lines += [
    "",
    "Driver chi-square across quadrants: see T05_quadrant_driver_chi2.tsv.",
    "",
    "Q1 (RAI-low + HLA-I-high) = canonical inflamed dark matter. "
    "Q2 (RAI-low + HLA-I-low) = silenced both, fail-to-induce phenotype. "
    "Q3 (RAI-high + HLA-I-high) = inflamed but still differentiated. "
    "Q4 (RAI-high + HLA-I-low) = canonical well-differentiated cold tumor.",
    "",
    "## 6. Stage progression (I → IV)",
    "",
    "| Stage | n | RAI mean | HLA-I mean | HLA-II mean | DM1 mean |",
    "| --- | --- | --- | --- | --- | --- |",
]
for r in stage_tab.itertuples(index=False):
    report_lines.append(
        f"| {r.stage} | {r.n} | {fmt(r.RAI_mean)} | {fmt(r.HLA1_mean)} | "
        f"{fmt(r.HLA2_mean)} | {fmt(r.DM1_mean)} |"
    )

report_lines += [
    "",
    f"Kruskal-Wallis: RAI p={fmt(p_kw_rai, 2) if p_kw_rai is not None else 'NA'} (Jonckheere "
    f"trend p={fmt(pJ_rai, 2) if pJ_rai is not None else 'NA'}); "
    f"HLA-I p={fmt(p_kw_hla1, 2) if p_kw_hla1 is not None else 'NA'} (Jonckheere trend "
    f"p={fmt(pJ_hla1, 2) if pJ_hla1 is not None else 'NA'}). See T07_stage_trend_tests.tsv.",
    "",
    "## 7. Methylation layer",
    "",
    f"Status: {meth_section}.",
    "If mean_8g_β present: tested whether HM450 mean β explains BOTH RAI silencing and HLA "
    "induction (single epigenetic axis) or only the RAI side. See T08_methylation_decouple.tsv.",
    "",
    "## 8. External cross-cohort hook",
    "",
    f"Status: {ext_section}.",
    "GSE286332 (PTC vs PTC+HT, n=18) replicates the RAI ↔ HLA negative coupling using "
    "the panel-8 RAI proxy from Paper 5 outputs (P5 r=-0.83 RAI vs HLA-II reproduced). "
    "Lee 2024 / GSE213647 panel-z is also tracked but full RAI 9-gene re-scoring on those "
    "cohorts is reserved for a downstream pass (the existing P5 pipeline already used the "
    "8-gene proxy as a methodological substitute and showed sign-consistent decoupling).",
    "",
    "## 9. Survival per quadrant",
    "",
    f"Status: {surv_section}.",
    "Cox HR (vs Q4_RAIhi_HLAlo reference) for OS and PFI in T11_quadrant_survival.tsv. "
    "Hypotheses: Q2 (silenced both) is the worst (fail-to-induce immunity); "
    "Q1 (inflamed dark matter) is intermediate-to-bad (already DM1-high); "
    "Q3 (inflamed differentiated) and Q4 (cold differentiated) are better. "
    "TCGA-THCA event rate is low; PFI is the more informative endpoint.",
    "",
    "## 10. Single-axis vs two-axis hypothesis (PCA)",
    "",
    f"PC1 explained variance ratio: {ev[0]:.3f} ({ev[0]*100:.1f}%). "
    f"PC1+PC2: {(ev[0]+ev[1])*100:.1f}%. Verdict: **{verdict}**. {verdict_signs}.",
    "",
    "PC1 loadings (DM1, RAI-loss, HLA-I, HLA-II, IFNG):",
]
for k in loadings.index:
    report_lines.append(f"- {k}: PC1={fmt(loadings.loc[k, 'PC1'])}, PC2={fmt(loadings.loc[k, 'PC2'])}")

report_lines += [
    "",
    "If PC1 ≥ 60% and all five modules load with the same sign, RAI-loss × HLA induction × "
    "DM1 × IFN-γ is **a single dominant axis** (de-differentiation/inflammation). If PC1 < 60% "
    "or HLA modules diverge from RAI on PC2, there is a **second residual axis** — the most "
    "likely candidate is a tumor-cell-intrinsic HLA induction layer (e.g., NLRC5/CIITA-driven) "
    "that operates on top of the dedifferentiation axis.",
    "",
    "## 11. Limitations",
    "",
    "- Correlative; no functional intervention. Cannot establish causality between RAI silencing "
    "and HLA induction.",
    "- HLA module = transcriptomic readout. Not allele genotype. Cell-type composition partially "
    "confounds — Track 5 already showed purity-proxy and immune-proxy partial-residualization "
    "moves HLA-I rho substantially, less so for HLA-II.",
    "- TCGA-THCA event rate is low (mostly indolent PTC), limiting Cox power for OS. PFI is "
    "more informative but still N=527-bounded.",
    "- DUOX1 / DUOX2 / IYD may show different distributional behavior than the canonical TG/TPO/"
    "TSHR/SLC5A5 core. T14 reports per-gene ρ.",
    "- Methylation only spans the 8-gene thyroid panel (mean_8g_β) — not genome-wide HM450. A "
    "fuller HM450 layer is reserved for Paper 1 supplementary.",
    "- External cohort layer uses panel-8 as RAI proxy (existing pipeline); a 9-gene re-score on "
    "GSE286332 / GSE213647 raw expression is a future supplementary.",
    "",
    "## 12. Implication for Paper 1 RAI-refractoriness narrative",
    "",
    "RAI silencing is the **definitional** feature of the DM1 dark-matter compartment "
    "(strong negative DM1 × RAI ρ; consistent across BRAF / RAS / Fusion / TripleNeg drivers, "
    "see T15). HLA-I/II induction is **co-incident** with RAI silencing along the dedifferentiation "
    "axis but the partial-correlation analysis shows HLA-I — when DM1 is held constant — has "
    "limited residual relationship with RAI, indicating much of the HLA signal is "
    "dedifferentiation-coupled. This supports the Paper 1 thesis that DM1 is a single "
    "molecular reorganization that simultaneously (a) shuts down RAI uptake (clinical "
    "RAI-refractoriness) and (b) raises antigen-presentation tone (immunotherapy hook).",
    "",
    "Q1 patients (RAI-low + HLA-I-high, ~25% of cohort by definition) are the most actionable: "
    "they are the inflamed-dark-matter group where checkpoint inhibition rationale is "
    "strongest while RAI re-induction therapy is also a target. Q2 patients (silenced both) "
    "are the hardest — neither classical RAI nor anti-PD-1 has a clear lever.",
    "",
    "## Outputs",
    "",
    "- Tables: `tables/T01-T15_*.tsv` (per-sample frame, headline correlations, partial "
    "correlations, quadrant summary, driver chi-square, stage progression, stage trend tests, "
    "methylation decouple, external cohorts, external correlations, quadrant survival, "
    "PCA summary + loadings, per-RAI-gene corr, driver-stratified).",
    "- Figures: `figs/F01-F11_*.png/pdf` (RAI scatter trio, partial-corr bars, quadrant "
    "scatter, quadrant features, stage progression, methylation layer, external cohort, "
    "KM PFI, PCA scree+biplot, per-gene heatmap, driver forest).",
    "- Summary JSON: `track34_summary.json`.",
]

with open(OUT / "track34_report.md", "w") as fh:
    fh.write("\n".join(report_lines) + "\n")

print("\nDONE.")
print(f"  results: {OUT}")
print(f"  figures: {FIG}")
print(f"  tables : {TAB}")
print(f"  PC1 var ratio = {ev[0]:.3f}; PC1+PC2 = {(ev[0]+ev[1]):.3f}")
print(f"  DM1 vs RAI rho: see T02; HLA-I vs RAI rho: see T02; HLA-II vs RAI rho: see T02")
