#!/usr/bin/env python3
"""
Track 26 - IFN-gamma pathway x HLA-I/II module x DM1 mediation in TCGA-THCA.

Boundary: HLA-I/II are gene-expression modules ONLY (transcript-level), not allele
genotypes. Per project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md, no allele
genotyping is performed on cancer cohorts in this track.

Question: is HLA module elevation along DM1 mediated by IFN-gamma signaling, or do
DM1 -> HLA paths exist independent of IFN-gamma? Answered with Baron-Kenny + Sobel,
direct/indirect decomposition (bootstrap), and conditional partial correlation.

Mediation framing here is observational, not interventional. Wording stays
'consistent with mediation', not 'proves mediation'.
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
import pingouin as pg
import statsmodels.api as sm
from scipy import stats

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/hla_deepdive_2026_05_08/track26_ifng_hla_dm1"
FIG = OUT / "figs"
TAB = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

CAPTION_BOILERPLATE = (
    "HLA-I/II gene-expression module - not allele genotype. "
    "Cancer-cohort allele genotyping is out of scope per separation rules."
)

# ============================================================
# 0. Inputs
# ============================================================
PANCAN_EXP = ROOT / "project/data/raw/TCGA_pancan/pancan_geneExp.gz"
PANCAN_PHENO = ROOT / "project/data/raw/TCGA_pancan/phenotype.tsv.gz"
TRACK5_T01 = (ROOT / "project/results/hla_deepdive_2026_05_08/"
              "track5_dm1_hla1_module/tables/T01_per_sample_module_scores.tsv")
METH_FILE = (ROOT / "project/results/audit_2026_04_30/round5/"
             "r5_2_sample_methylation_8gene.tsv")

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

# Ayers 2017 18-gene tumor inflammation signature (TIS), preferred IFN-gamma set
# for tumor RNA-seq (NEJM/JCI 2017; gene order = canonical).
AYERS_TIS = [
    "IFNG", "STAT1", "CCR5", "CXCL9", "CXCL10", "CXCL11",
    "IDO1", "PRF1", "GZMA", "HLA-DRA", "CD8A", "CD27",
    "HLA-E", "NKG7", "PSMB10", "CMKLR1", "TIGIT", "LAG3",
]

# HALLMARK_INTERFERON_GAMMA_RESPONSE (MSigDB 200-gene set, current release).
# Verified against MSigDB v2024.1.Hs core members (commonly cited subset);
# we include the main canonical members. Missing genes silently dropped.
HALLMARK_IFNG = [
    "ADAR", "APOL6", "ARID5B", "ARL4A", "AUTS2", "B2M", "BANK1", "BATF2", "BPGM",
    "BST2", "BTG1", "C1R", "C1S", "CASP1", "CASP3", "CASP4", "CASP7", "CASP8",
    "CCL2", "CCL5", "CCL7", "CD274", "CD38", "CD40", "CD69", "CD74", "CD86",
    "CDKN1A", "CFB", "CFH", "CIITA", "CMKLR1", "CMPK2", "CMTR1", "CSF2RB",
    "CXCL10", "CXCL11", "CXCL9", "DDX58", "DDX60", "DHX58", "EIF2AK2", "EIF4E3",
    "EPSTI1", "FAS", "FCGR1A", "FGL2", "FPR1", "GBP4", "GBP6", "GCH1", "GPR18",
    "GZMA", "HELZ2", "HERC6", "HIF1A", "HLA-A", "HLA-B", "HLA-DMA", "HLA-DQA1",
    "HLA-DRB1", "HLA-G", "ICAM1", "IDO1", "IFI27", "IFI30", "IFI35", "IFI44",
    "IFI44L", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFITM2", "IFITM3", "IFNAR2",
    "IL10RA", "IL15", "IL15RA", "IL18BP", "IL2RB", "IL4R", "IL6", "IL7", "IRF1",
    "IRF2", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "ISG15", "ISG20", "ISOC1",
    "ITGB7", "JAK2", "KLRK1", "LAP3", "LATS2", "LCP2", "LGALS3BP", "LY6E",
    "LYSMD2", "MARCHF1", "METTL7B", "MT2A", "MTHFD2", "MVP", "MX1", "MX2",
    "MYD88", "NAMPT", "NCOA3", "NFKB1", "NFKBIA", "NLRC5", "NMI", "NOD1",
    "NUP93", "OAS2", "OAS3", "OASL", "OGFR", "P2RY14", "PARP12", "PARP14",
    "PDE4B", "PELI1", "PFKP", "PIM1", "PLA2G4A", "PLSCR1", "PML", "PNP",
    "PNPT1", "PSMA2", "PSMA3", "PSMB10", "PSMB2", "PSMB8", "PSMB9", "PSME1",
    "PSME2", "PTGS2", "PTPN1", "PTPN2", "PTPN6", "RAPGEF6", "RBCK1", "RIPK1",
    "RIPK2", "RNF31", "RSAD2", "RTP4", "SAMD9L", "SAMHD1", "SECTM1", "SELP",
    "SERPING1", "SLAMF7", "SLC25A28", "SOCS1", "SOCS3", "SOD2", "SP110",
    "SPPL2A", "SRI", "SSPN", "ST3GAL5", "ST8SIA4", "STAT1", "STAT2", "STAT3",
    "STAT4", "TAP1", "TAPBP", "TDRD7", "TNFAIP2", "TNFAIP3", "TNFAIP6",
    "TNFSF10", "TOR1B", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26",
    "TXNIP", "UBE2L6", "UPP1", "USP18", "VAMP5", "VAMP8", "VCAM1", "WARS1",
    "XAF1", "XCL1", "ZBP1", "ZNFX1",
]

# Pretty per-locus splits to study which HLA-I genes are most/least IFNG-mediated.
HLA1_PER_GENE_FOCUS = [
    "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9",
    "NLRC5", "IRF1", "ERAP1", "ERAP2", "HLA-E", "HLA-G",
    "CALR", "CANX", "PDIA3",  # chaperones, expected less IFNG-driven
]

# ============================================================
# 1. Load existing per-sample module scores (Track 5)
# ============================================================
print("[1/12] loading Track 5 per-sample module scores")
t5 = pd.read_csv(TRACK5_T01, sep="\t")
print(f"  Track 5 frame: {t5.shape}")

# ============================================================
# 2. Pull IFN-gamma genes from TCGA pancan expression for THCA primary
# ============================================================
print("[2/12] streaming IFN-gamma + extra HLA gene expression for THCA primary samples")
pheno = pd.read_csv(PANCAN_PHENO, sep="\t")
thy_primary = set(pheno[(pheno["_primary_disease"] == "thyroid carcinoma") &
                        (pheno["sample_type"] == "Primary Tumor")]["sample"])
print(f"  THCA primary samples in pheno: {len(thy_primary)}")

needed_genes = set(HALLMARK_IFNG) | set(AYERS_TIS) | set(HLA1_GENES) | set(HLA2_GENES)

with gzip.open(PANCAN_EXP, "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    sample_cols = header[1:]
    keep_idx = [i for i, s in enumerate(sample_cols) if s in thy_primary]
    keep_samples = [sample_cols[i] for i in keep_idx]
    rows = {}
    for ln in fh:
        gene, *vals = ln.rstrip("\n").split("\t")
        if gene in needed_genes:
            arr = np.array([vals[i] for i in keep_idx], dtype=float)
            rows[gene] = arr

expr = pd.DataFrame(rows, index=keep_samples).T  # genes x samples -> wait
expr = expr.T  # samples x genes
print(f"  expr shape (samples x genes): {expr.shape}")
print(f"  Hallmark IFNG genes present: "
      f"{len([g for g in HALLMARK_IFNG if g in expr.columns])}/{len(HALLMARK_IFNG)}")
print(f"  Ayers TIS genes present: "
      f"{len([g for g in AYERS_TIS if g in expr.columns])}/{len(AYERS_TIS)}")


def zmean_module(df, gene_list):
    sub = df.reindex(columns=[g for g in gene_list if g in df.columns])
    if sub.shape[1] == 0:
        return pd.Series(np.nan, index=df.index)
    z = (sub - sub.mean()) / sub.std(ddof=0)
    return z.mean(axis=1)


ifng_hallmark = zmean_module(expr, HALLMARK_IFNG)
ifng_ayers = zmean_module(expr, AYERS_TIS)
hla1_recompute = zmean_module(expr, HLA1_GENES)
hla2_recompute = zmean_module(expr, HLA2_GENES)

ifng_df = pd.DataFrame({
    "sample": expr.index,
    "IFNG_hallmark": ifng_hallmark.values,
    "IFNG_ayers": ifng_ayers.values,
    "HLA1_recompute": hla1_recompute.values,
    "HLA2_recompute": hla2_recompute.values,
})
ifng_df.to_csv(TAB / "T01_ifng_score_per_sample.tsv", sep="\t", index=False)
print(f"  IFN-gamma score Spearman (Hallmark vs Ayers): "
      f"{stats.spearmanr(ifng_df['IFNG_hallmark'], ifng_df['IFNG_ayers']).statistic:.3f}")

# Merge into Track 5 frame
df = t5.merge(ifng_df, on="sample", how="inner")
# Sanity: Track 5 HLA1_score vs recomputed HLA1 should be ~1
sanity_rho = stats.spearmanr(df["HLA1_score"], df["HLA1_recompute"]).statistic
print(f"  Track 5 vs recomputed HLA1 Spearman = {sanity_rho:.3f} (sanity)")
print(f"  merged Track 26 frame: {df.shape}")

# ============================================================
# 3. DM1 x IFN-gamma headline correlations + scatter
# ============================================================
print("[3/12] DM1 x IFN-gamma correlations")


def headline_correlation(d, x, y, label):
    sub = d[[x, y]].dropna()
    rho_s, p_s = stats.spearmanr(sub[x], sub[y])
    rho_p, p_p = stats.pearsonr(sub[x], sub[y])
    return {"comparison": label, "n": len(sub),
            "spearman_rho": rho_s, "spearman_p": p_s,
            "pearson_r": rho_p, "pearson_p": p_p}


corr_rows = []
for pair, lbl in [
    (("DM1_use", "IFNG_hallmark"), "DM1 vs IFN-gamma Hallmark"),
    (("DM1_use", "IFNG_ayers"), "DM1 vs IFN-gamma Ayers TIS"),
    (("IFNG_hallmark", "HLA1_score"), "IFN-gamma Hallmark vs HLA-I module"),
    (("IFNG_hallmark", "HLA2_score"), "IFN-gamma Hallmark vs HLA-II module"),
    (("IFNG_ayers", "HLA1_score"), "IFN-gamma Ayers vs HLA-I module"),
    (("IFNG_ayers", "HLA2_score"), "IFN-gamma Ayers vs HLA-II module"),
    (("DM1_use", "HLA1_score"), "DM1 vs HLA-I module (reference)"),
    (("DM1_use", "HLA2_score"), "DM1 vs HLA-II module (reference)"),
]:
    corr_rows.append(headline_correlation(df, pair[0], pair[1], lbl))

corr_tab = pd.DataFrame(corr_rows)
corr_tab.to_csv(TAB / "T02_headline_correlations.tsv", sep="\t", index=False)
print(corr_tab.to_string(index=False))

# Figure 1: scatter DM1 vs IFNG (both signatures)
try:
    from statsmodels.nonparametric.smoothers_lowess import lowess
    HAVE_LOWESS = True
except Exception:
    HAVE_LOWESS = False


def scatter_with_lowess(x, y, ax, title, xlabel, ylabel, color="#2b8cbe"):
    ax.scatter(x, y, s=8, alpha=0.5, color=color, edgecolor="none")
    if HAVE_LOWESS and len(x) > 30:
        sm_lo = lowess(y, x, frac=0.4, return_sorted=True)
        ax.plot(sm_lo[:, 0], sm_lo[:, 1], color="#e34a33", lw=2)
    rho, p = stats.spearmanr(x, y)
    ax.set_title(f"{title}\nSpearman rho={rho:.3f}, p={p:.2e}, n={len(x)}", fontsize=9)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.3)


fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
sub = df[["DM1_use", "IFNG_hallmark"]].dropna()
scatter_with_lowess(sub["DM1_use"].values, sub["IFNG_hallmark"].values, axes[0],
                    "DM1 vs IFN-gamma Hallmark", "DM1 score", "IFN-gamma Hallmark z-mean")
sub = df[["DM1_use", "IFNG_ayers"]].dropna()
scatter_with_lowess(sub["DM1_use"].values, sub["IFNG_ayers"].values, axes[1],
                    "DM1 vs IFN-gamma Ayers TIS", "DM1 score", "Ayers TIS z-mean",
                    color="#33a02c")
fig.suptitle(CAPTION_BOILERPLATE, fontsize=8, color="#666")
fig.tight_layout()
fig.savefig(FIG / "F01_dm1_vs_ifng_scatter.png", dpi=180)
fig.savefig(FIG / "F01_dm1_vs_ifng_scatter.pdf")
plt.close(fig)

# Figure 2: scatter IFNG vs HLA-I/II (Hallmark)
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
sub = df[["IFNG_hallmark", "HLA1_score"]].dropna()
scatter_with_lowess(sub["IFNG_hallmark"].values, sub["HLA1_score"].values, axes[0],
                    "IFN-gamma Hallmark vs HLA-I module", "IFN-gamma Hallmark", "HLA-I module")
sub = df[["IFNG_hallmark", "HLA2_score"]].dropna()
scatter_with_lowess(sub["IFNG_hallmark"].values, sub["HLA2_score"].values, axes[1],
                    "IFN-gamma Hallmark vs HLA-II module", "IFN-gamma Hallmark", "HLA-II module",
                    color="#762a83")
fig.suptitle(CAPTION_BOILERPLATE, fontsize=8, color="#666")
fig.tight_layout()
fig.savefig(FIG / "F02_ifng_vs_hla_scatter.png", dpi=180)
fig.savefig(FIG / "F02_ifng_vs_hla_scatter.pdf")
plt.close(fig)


# ============================================================
# 4. Mediation analysis: DM1 -> IFNG -> HLA-I/II (Baron-Kenny + Sobel + bootstrap)
# ============================================================
print("[4/12] Baron-Kenny + Sobel mediation: DM1 -> IFNG -> HLA modules")


def run_mediation(dframe, x, m, y, label, n_boot=2000, seed=42):
    """Run pingouin mediation_analysis (bias-corrected non-parametric bootstrap).

    Returns a one-row dict with total/direct/indirect, % mediated, Sobel-style z,
    and bootstrap CI for indirect effect.
    """
    sub = dframe[[x, m, y]].dropna().copy()
    n = len(sub)
    res = pg.mediation_analysis(data=sub, x=x, m=m, y=y, n_boot=n_boot, seed=seed)
    # pingouin 0.6 path labels: 'm ~ X', 'Y ~ m', 'Total', 'Direct', 'Indirect'.
    # CI columns are CI2.5 / CI97.5.
    a = res[res["path"] == f"{m} ~ X"].iloc[0]
    b = res[res["path"] == f"Y ~ {m}"].iloc[0]
    total = res[res["path"] == "Total"].iloc[0]
    direct = res[res["path"] == "Direct"].iloc[0]
    indirect = res[res["path"] == "Indirect"].iloc[0]
    sobel_z = a["coef"] * b["coef"] / np.sqrt(
        (b["coef"] ** 2) * (a["se"] ** 2) + (a["coef"] ** 2) * (b["se"] ** 2)
    )
    sobel_p = 2 * (1 - stats.norm.cdf(abs(sobel_z)))
    pct_mediated = indirect["coef"] / total["coef"] * 100 if total["coef"] != 0 else np.nan
    return {
        "label": label, "x": x, "m": m, "y": y, "n": n,
        "a_X_to_M": a["coef"], "a_se": a["se"], "a_p": a["pval"],
        "b_M_to_Y_given_X": b["coef"], "b_se": b["se"], "b_p": b["pval"],
        "c_total": total["coef"], "c_total_p": total["pval"],
        "c_prime_direct": direct["coef"], "c_prime_p": direct["pval"],
        "ab_indirect": indirect["coef"], "ab_indirect_p": indirect["pval"],
        "ab_CI95_lo": indirect["CI2.5"], "ab_CI95_hi": indirect["CI97.5"],
        "sobel_z": sobel_z, "sobel_p": sobel_p,
        "pct_mediated": pct_mediated,
    }


med_rows = []
for x_var in ["DM1_use"]:
    for m_var, m_lbl in [("IFNG_hallmark", "IFNG Hallmark"), ("IFNG_ayers", "Ayers TIS")]:
        for y_var, y_lbl in [("HLA1_score", "HLA-I"), ("HLA2_score", "HLA-II")]:
            lbl = f"DM1 -> {m_lbl} -> {y_lbl}"
            try:
                med_rows.append(run_mediation(df, x_var, m_var, y_var, lbl))
            except Exception as e:
                print(f"  mediation failed for {lbl}: {e}")

med_tab = pd.DataFrame(med_rows)
med_tab.to_csv(TAB / "T03_mediation_baron_kenny_sobel.tsv", sep="\t", index=False)
print(med_tab[[
    "label", "n", "c_total", "c_prime_direct", "ab_indirect",
    "pct_mediated", "sobel_z", "sobel_p"
]].to_string(index=False))

# Figure 3: bar of % mediated for each combination
fig, ax = plt.subplots(figsize=(9, 4.5))
y_pos = np.arange(len(med_tab))
ax.barh(y_pos, med_tab["pct_mediated"].values, color="#2b8cbe", alpha=0.8)
for i, (pct, sob) in enumerate(zip(med_tab["pct_mediated"], med_tab["sobel_p"])):
    ax.text(pct + 0.5, i, f"{pct:.1f}% (Sobel p={sob:.1e})", va="center", fontsize=8)
ax.axvline(0, color="k", lw=0.6)
ax.set_yticks(y_pos)
ax.set_yticklabels(med_tab["label"].tolist())
ax.set_xlabel("% mediated through IFN-gamma (indirect / total)")
ax.set_title(f"DM1 -> IFN-gamma -> HLA mediation in TCGA-THCA\n{CAPTION_BOILERPLATE}",
             fontsize=9)
ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
fig.savefig(FIG / "F03_mediation_pct_mediated.png", dpi=180)
fig.savefig(FIG / "F03_mediation_pct_mediated.pdf")
plt.close(fig)

# Figure 4: bootstrap CI plot for indirect effects
fig, ax = plt.subplots(figsize=(9, 4.5))
y_pos = np.arange(len(med_tab))
xerr_lo = (med_tab["ab_indirect"] - med_tab["ab_CI95_lo"]).values
xerr_hi = (med_tab["ab_CI95_hi"] - med_tab["ab_indirect"]).values
ax.errorbar(med_tab["ab_indirect"].values, y_pos,
            xerr=[xerr_lo, xerr_hi], fmt="o", color="#1b7837",
            markersize=8, capsize=4)
ax.axvline(0, color="k", lw=0.7)
ax.set_yticks(y_pos)
ax.set_yticklabels(med_tab["label"].tolist())
ax.set_xlabel("Indirect effect (a*b) with bootstrap 95% CI")
ax.set_title(f"Bootstrap 95% CI for indirect effect (DM1 -> IFNG -> HLA)\n{CAPTION_BOILERPLATE}",
             fontsize=9)
ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
fig.savefig(FIG / "F04_indirect_effect_bootstrap_CI.png", dpi=180)
fig.savefig(FIG / "F04_indirect_effect_bootstrap_CI.pdf")
plt.close(fig)


# ============================================================
# 5. Conditional independence: partial Spearman DM1 vs HLA | IFN-gamma
# ============================================================
print("[5/12] partial Spearman / conditional independence DM1 x HLA-I | IFN-gamma")


def partial_spearman(dframe, x, y, z_cols):
    sub = dframe[[x, y] + list(z_cols)].dropna()
    if len(sub) < 20:
        return np.nan, np.nan, len(sub)
    rx = sub[x].rank().values
    ry = sub[y].rank().values
    Z = sub[list(z_cols)].rank().values
    Z = (Z - Z.mean(axis=0)) / Z.std(axis=0, ddof=0)
    rx_resid = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
    ry_resid = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    rho, p = stats.spearmanr(rx_resid, ry_resid)
    return rho, p, len(sub)


cond_rows = []
for x, y, ylab in [
    ("DM1_use", "HLA1_score", "HLA-I"),
    ("DM1_use", "HLA2_score", "HLA-II"),
]:
    raw_rho, raw_p = stats.spearmanr(df[x], df[y], nan_policy="omit")
    rho_h, p_h, n_h = partial_spearman(df, x, y, ["IFNG_hallmark"])
    rho_a, p_a, n_a = partial_spearman(df, x, y, ["IFNG_ayers"])
    rho_p, p_p, n_p = partial_spearman(df, x, y, ["Purity_proxy"])
    rho_hp, p_hp, n_hp = partial_spearman(df, x, y, ["IFNG_hallmark", "Purity_proxy"])
    cond_rows.append({
        "module": ylab, "n_raw": int(df[[x, y]].dropna().shape[0]),
        "rho_raw": raw_rho, "p_raw": raw_p,
        "rho_partial_IFNG_hallmark": rho_h, "p_partial_IFNG_hallmark": p_h, "n_h": n_h,
        "rho_partial_IFNG_ayers": rho_a, "p_partial_IFNG_ayers": p_a, "n_a": n_a,
        "rho_partial_Purity": rho_p, "p_partial_Purity": p_p, "n_p": n_p,
        "rho_partial_IFNG+Purity": rho_hp, "p_partial_IFNG+Purity": p_hp, "n_hp": n_hp,
    })
cond_tab = pd.DataFrame(cond_rows)
cond_tab.to_csv(TAB / "T04_conditional_independence.tsv", sep="\t", index=False)
print(cond_tab.to_string(index=False))

fig, ax = plt.subplots(figsize=(9, 4.5))
labels = cond_tab["module"].tolist()
x_pos = np.arange(len(labels))
w = 0.16
ax.bar(x_pos - 2 * w, cond_tab["rho_raw"], width=w, label="raw", color="#2b8cbe")
ax.bar(x_pos - w, cond_tab["rho_partial_IFNG_hallmark"], width=w,
       label="| IFNG Hallmark", color="#1b7837")
ax.bar(x_pos, cond_tab["rho_partial_IFNG_ayers"], width=w,
       label="| IFNG Ayers", color="#33a02c")
ax.bar(x_pos + w, cond_tab["rho_partial_Purity"], width=w,
       label="| Purity proxy", color="#e6ab02")
ax.bar(x_pos + 2 * w, cond_tab["rho_partial_IFNG+Purity"], width=w,
       label="| IFNG + Purity", color="#762a83")
ax.axhline(0, color="k", lw=0.6)
ax.set_xticks(x_pos)
ax.set_xticklabels(labels)
ax.set_ylabel("Spearman rho with DM1")
ax.set_title(f"Conditional independence: DM1 x HLA module | IFN-gamma\n{CAPTION_BOILERPLATE}",
             fontsize=9)
ax.legend(fontsize=8, loc="best")
ax.grid(alpha=0.3, axis="y")
fig.tight_layout()
fig.savefig(FIG / "F05_conditional_independence.png", dpi=180)
fig.savefig(FIG / "F05_conditional_independence.pdf")
plt.close(fig)


# ============================================================
# 6. Driver-stratified mediation
# ============================================================
print("[6/12] driver-stratified mediation (BRAF / RAS / Fusion / TripleNeg)")
driver_rows = []
for driver in ["BRAF", "RAS", "Fusion", "TripleNeg"]:
    sub = df[df["driver_simple"] == driver]
    n_sub = len(sub.dropna(subset=["DM1_use", "IFNG_hallmark", "HLA1_score"]))
    if n_sub < 20:
        print(f"  {driver}: n={n_sub} too small, skipping")
        continue
    for y_var, y_lbl in [("HLA1_score", "HLA-I"), ("HLA2_score", "HLA-II")]:
        try:
            r = run_mediation(sub, "DM1_use", "IFNG_hallmark", y_var,
                              f"{driver} | DM1 -> IFNG_hallmark -> {y_lbl}",
                              n_boot=1000)
            r["driver"] = driver
            r["target"] = y_lbl
            driver_rows.append(r)
        except Exception as e:
            print(f"  {driver} {y_lbl} failed: {e}")

drv_med_tab = pd.DataFrame(driver_rows)
drv_med_tab.to_csv(TAB / "T05_driver_stratified_mediation.tsv", sep="\t", index=False)
print(drv_med_tab[[
    "driver", "target", "n", "c_total", "c_prime_direct", "ab_indirect",
    "pct_mediated", "sobel_p"
]].to_string(index=False))

# Figure 6: driver x target % mediated heatmap
if len(drv_med_tab) > 0:
    pivot_pct = drv_med_tab.pivot(index="driver", columns="target", values="pct_mediated")
    fig, ax = plt.subplots(figsize=(7, 4.5))
    im = ax.imshow(pivot_pct.values, aspect="auto", cmap="RdBu_r",
                   vmin=-100, vmax=100)
    ax.set_xticks(range(pivot_pct.shape[1]))
    ax.set_xticklabels(pivot_pct.columns)
    ax.set_yticks(range(pivot_pct.shape[0]))
    ax.set_yticklabels(pivot_pct.index)
    for i in range(pivot_pct.shape[0]):
        for j in range(pivot_pct.shape[1]):
            v = pivot_pct.values[i, j]
            ax.text(j, i, f"{v:.0f}%" if not np.isnan(v) else "n/a",
                    ha="center", va="center", fontsize=10,
                    color="white" if abs(v) > 50 else "black")
    plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="% mediated")
    ax.set_title(f"Driver-stratified mediation: % of DM1->HLA via IFN-gamma\n{CAPTION_BOILERPLATE}",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "F06_driver_stratified_mediation_heatmap.png", dpi=180)
    fig.savefig(FIG / "F06_driver_stratified_mediation_heatmap.pdf")
    plt.close(fig)

# ============================================================
# 7. Per-gene split: which HLA-I genes are most/least IFN-gamma-mediated
# ============================================================
print("[7/12] per-gene IFN-gamma mediation across HLA-I/II genes")

# We need per-gene expression of HLA-I and HLA-II focus genes; expr already has them
per_gene_med_rows = []
for g in HLA1_GENES + HLA2_GENES:
    if g not in expr.columns:
        continue
    g_vals = expr[g].reindex(df["sample"].values).values
    sub = df.copy()
    sub[f"_{g}"] = g_vals
    sub_clean = sub[["DM1_use", "IFNG_hallmark", f"_{g}"]].dropna()
    if len(sub_clean) < 50:
        continue
    try:
        r = run_mediation(sub_clean.rename(columns={f"_{g}": "Y"}),
                          "DM1_use", "IFNG_hallmark", "Y", f"{g}", n_boot=500)
        r["gene"] = g
        r["module"] = "HLA-I" if g in HLA1_GENES else "HLA-II"
        per_gene_med_rows.append(r)
    except Exception as e:
        pass

pg_med_tab = pd.DataFrame(per_gene_med_rows)
pg_med_tab = pg_med_tab.sort_values("pct_mediated", ascending=False)
pg_med_tab.to_csv(TAB / "T06_per_gene_mediation.tsv", sep="\t", index=False)
print(pg_med_tab[["gene", "module", "n", "c_total", "ab_indirect",
                   "pct_mediated", "sobel_p"]].to_string(index=False))

# Figure 7: per-gene % mediated lollipop
fig, ax = plt.subplots(figsize=(10, 8))
plot_df = pg_med_tab.dropna(subset=["pct_mediated"]).copy()
# Cap extreme values for readability (mediation pct can blow up if c_total is tiny)
plot_df["pct_capped"] = plot_df["pct_mediated"].clip(-150, 200)
colors = ["#1b7837" if m == "HLA-I" else "#762a83" for m in plot_df["module"]]
y_pos = np.arange(len(plot_df))
ax.hlines(y_pos, 0, plot_df["pct_capped"].values, colors=colors, alpha=0.6)
ax.scatter(plot_df["pct_capped"].values, y_pos, c=colors, s=50)
for i, (g, p) in enumerate(zip(plot_df["gene"], plot_df["sobel_p"])):
    ax.text(plot_df["pct_capped"].values[i] + 3, i,
            f"p={p:.1e}", fontsize=7, va="center")
ax.axvline(0, color="k", lw=0.6)
ax.axvline(100, color="grey", lw=0.5, linestyle="--")
ax.set_yticks(y_pos)
ax.set_yticklabels(plot_df["gene"].tolist(), fontsize=8)
ax.set_xlabel("% of DM1 -> gene effect mediated by IFN-gamma Hallmark (capped at -150/+200)")
ax.set_title(f"Per-gene IFN-gamma mediation: HLA-I (green) and HLA-II (purple)\n{CAPTION_BOILERPLATE}",
             fontsize=9)
ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
fig.savefig(FIG / "F07_per_gene_mediation_lollipop.png", dpi=180)
fig.savefig(FIG / "F07_per_gene_mediation_lollipop.pdf")
plt.close(fig)


# ============================================================
# 8. Methylation 3-way: DM1 -> methylation -> HLA-I and DM1 -> IFNG -> HLA-I
# ============================================================
print("[8/12] methylation 3-way mediation")
meth_section = "skipped"
df_meth = None
if METH_FILE.exists():
    meth = pd.read_csv(METH_FILE, sep="\t")
    meth["sample"] = meth["sample_short"].astype(str) + "-01"
    df_meth = df.merge(meth[["sample", "mean_8g_beta"]], on="sample", how="inner")
    print(f"  joined methylation: {len(df_meth)} samples")
    # Single-mediator: DM1 -> mean_8g_beta -> HLA-I
    # And dual-mediator: DM1 -> [mean_8g_beta, IFNG_hallmark] -> HLA-I
    rows = []
    for y_var, y_lbl in [("HLA1_score", "HLA-I"), ("HLA2_score", "HLA-II")]:
        # Single mediator: methylation
        try:
            r = run_mediation(df_meth, "DM1_use", "mean_8g_beta", y_var,
                              f"DM1 -> mean_8g_beta -> {y_lbl}", n_boot=2000)
            rows.append(r)
        except Exception as e:
            print(f"  meth single failed {y_lbl}: {e}")
        # Single mediator: IFNG (in same intersected sample set, for fair comparison)
        try:
            r = run_mediation(df_meth, "DM1_use", "IFNG_hallmark", y_var,
                              f"DM1 -> IFNG_hallmark -> {y_lbl} (meth-intersect)",
                              n_boot=2000)
            rows.append(r)
        except Exception as e:
            print(f"  meth IFNG-single failed {y_lbl}: {e}")
    # Dual mediator: pingouin supports list of mediators
    dual_rows = []
    for y_var, y_lbl in [("HLA1_score", "HLA-I"), ("HLA2_score", "HLA-II")]:
        sub = df_meth[["DM1_use", "mean_8g_beta", "IFNG_hallmark", y_var]].dropna()
        try:
            res = pg.mediation_analysis(
                data=sub, x="DM1_use", m=["mean_8g_beta", "IFNG_hallmark"],
                y=y_var, n_boot=2000, seed=42,
            )
            # res rows include 'Indirect mean_8g_beta' and 'Indirect IFNG_hallmark'
            res = res.assign(target=y_lbl, n=len(sub))
            dual_rows.append(res)
        except Exception as e:
            print(f"  dual mediator failed {y_lbl}: {e}")
    if dual_rows:
        dual = pd.concat(dual_rows, ignore_index=True)
        dual.to_csv(TAB / "T07_dual_mediator_meth_ifng.tsv", sep="\t", index=False)
    if rows:
        meth_med_tab = pd.DataFrame(rows)
        meth_med_tab.to_csv(TAB / "T08_methylation_single_mediator.tsv", sep="\t", index=False)
        meth_section = f"used n={len(df_meth)}"
        print(meth_med_tab[["label", "n", "c_total", "ab_indirect",
                             "pct_mediated", "sobel_p"]].to_string(index=False))

    # Partial: DM1 x HLA-I | (mean_8g_beta + IFNG)
    pir_rows = []
    for y_var, y_lbl in [("HLA1_score", "HLA-I"), ("HLA2_score", "HLA-II")]:
        raw_rho, raw_p = stats.spearmanr(df_meth["DM1_use"], df_meth[y_var],
                                         nan_policy="omit")
        rho_meth, p_meth, n1 = partial_spearman(df_meth, "DM1_use", y_var,
                                                 ["mean_8g_beta"])
        rho_ifng, p_ifng, n2 = partial_spearman(df_meth, "DM1_use", y_var,
                                                 ["IFNG_hallmark"])
        rho_both, p_both, n3 = partial_spearman(df_meth, "DM1_use", y_var,
                                                 ["mean_8g_beta", "IFNG_hallmark"])
        pir_rows.append({
            "module": y_lbl, "n": len(df_meth.dropna(subset=["DM1_use", y_var])),
            "rho_raw": raw_rho, "p_raw": raw_p,
            "rho_partial_meth": rho_meth, "p_partial_meth": p_meth,
            "rho_partial_IFNG": rho_ifng, "p_partial_IFNG": p_ifng,
            "rho_partial_both": rho_both, "p_partial_both": p_both,
        })
    pir_tab = pd.DataFrame(pir_rows)
    pir_tab.to_csv(TAB / "T09_partial_meth_ifng.tsv", sep="\t", index=False)
    print(pir_tab.to_string(index=False))

    # Figure 8: dual-mediator decomposition stacked bars
    if dual_rows:
        fig, ax = plt.subplots(figsize=(9, 5))
        # Pull out indirect rows per target
        cats = []
        for tgt in ["HLA-I", "HLA-II"]:
            d = dual[dual["target"] == tgt]
            ind_meth = d[d["path"].str.contains("Indirect mean_8g_beta")]["coef"]
            ind_ifng = d[d["path"].str.contains("Indirect IFNG_hallmark")]["coef"]
            direct = d[d["path"] == "Direct"]["coef"]
            total = d[d["path"] == "Total"]["coef"]
            cats.append({
                "target": tgt,
                "indirect_meth": float(ind_meth.values[0]) if len(ind_meth) else 0,
                "indirect_ifng": float(ind_ifng.values[0]) if len(ind_ifng) else 0,
                "direct": float(direct.values[0]) if len(direct) else 0,
                "total": float(total.values[0]) if len(total) else 0,
            })
        cat_df = pd.DataFrame(cats)
        x_pos = np.arange(len(cat_df))
        bars = ["indirect_meth", "indirect_ifng", "direct"]
        colors_b = ["#762a83", "#1b7837", "#666666"]
        labels_b = ["Indirect via methylation", "Indirect via IFN-gamma", "Direct (residual)"]
        bottom = np.zeros(len(cat_df))
        for bn, c, lb in zip(bars, colors_b, labels_b):
            vals = cat_df[bn].values
            ax.bar(x_pos, vals, bottom=bottom, color=c, label=lb)
            bottom += vals
        # total marker
        ax.scatter(x_pos, cat_df["total"].values, color="black", s=80,
                   marker="D", label="Total c", zorder=5)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(cat_df["target"].tolist())
        ax.set_ylabel("Effect size (regression coefficient)")
        ax.set_title(f"Dual-mediator decomposition: DM1 -> [methylation, IFN-gamma] -> HLA module\n{CAPTION_BOILERPLATE}",
                     fontsize=9)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3, axis="y")
        fig.tight_layout()
        fig.savefig(FIG / "F08_dual_mediator_decomposition.png", dpi=180)
        fig.savefig(FIG / "F08_dual_mediator_decomposition.pdf")
        plt.close(fig)
else:
    print("  methylation file not present; skipping 3-way mediation")

# ============================================================
# 9. Direct vs indirect summary chart (HLA-I and HLA-II focus)
# ============================================================
print("[9/12] direct vs indirect summary chart")
fig, ax = plt.subplots(figsize=(9, 5))
focus = med_tab[med_tab["m"] == "IFNG_hallmark"].copy()
x_pos = np.arange(len(focus))
direct_vals = focus["c_prime_direct"].values
indirect_vals = focus["ab_indirect"].values
total_vals = focus["c_total"].values
w = 0.27
ax.bar(x_pos - w, total_vals, width=w, label="Total c", color="#2b8cbe")
ax.bar(x_pos, direct_vals, width=w, label="Direct c'", color="#666666")
ax.bar(x_pos + w, indirect_vals, width=w, label="Indirect ab (via IFN-gamma)",
       color="#1b7837")
for i, (t, d_v, ind) in enumerate(zip(total_vals, direct_vals, indirect_vals)):
    ax.text(i - w, t + 0.005, f"{t:.2f}", ha="center", fontsize=7)
    ax.text(i, d_v + 0.005, f"{d_v:.2f}", ha="center", fontsize=7)
    ax.text(i + w, ind + 0.005, f"{ind:.2f}", ha="center", fontsize=7)
ax.axhline(0, color="k", lw=0.6)
ax.set_xticks(x_pos)
ax.set_xticklabels([f.replace("DM1 -> IFNG Hallmark -> ", "") for f in focus["label"]])
ax.set_ylabel("Effect size (regression coefficient)")
ax.set_title(f"DM1 -> HLA module: total vs direct vs indirect (Hallmark IFN-gamma)\n{CAPTION_BOILERPLATE}",
             fontsize=9)
ax.legend()
ax.grid(alpha=0.3, axis="y")
fig.tight_layout()
fig.savefig(FIG / "F09_direct_indirect_summary.png", dpi=180)
fig.savefig(FIG / "F09_direct_indirect_summary.pdf")
plt.close(fig)


# ============================================================
# 10. Final summary JSON
# ============================================================
print("[10/12] writing summary JSON")
summary = {
    "track": "Track 26 - IFN-gamma x HLA-I/II x DM1 mediation in TCGA-THCA",
    "boundary": ("HLA-I/II gene-expression module ONLY (no allele genotyping). "
                 "Cancer-cohort allele genotyping out of scope per "
                 "HLA_CANCER_SEPARATION_RULES.md (Sections 1.1, 1.2). "
                 "Mediation framing observational, not interventional."),
    "n_samples": int(len(df)),
    "ifng_signatures": {
        "Hallmark_present": int(sum(g in expr.columns for g in HALLMARK_IFNG)),
        "Hallmark_total": len(HALLMARK_IFNG),
        "Ayers_present": int(sum(g in expr.columns for g in AYERS_TIS)),
        "Ayers_total": len(AYERS_TIS),
        "Hallmark_vs_Ayers_spearman": float(stats.spearmanr(
            ifng_df["IFNG_hallmark"], ifng_df["IFNG_ayers"]).statistic),
    },
    "headline_correlations": corr_tab.to_dict("records"),
    "mediation_baron_kenny_sobel": med_tab.to_dict("records"),
    "conditional_independence": cond_tab.to_dict("records"),
    "driver_stratified_mediation": drv_med_tab.to_dict("records"),
    "per_gene_mediation_top": pg_med_tab.head(10).to_dict("records"),
    "per_gene_mediation_bot": pg_med_tab.tail(10).to_dict("records"),
    "methylation_status": meth_section,
    "outputs": {
        "results_dir": str(OUT),
        "figs_dir": str(FIG),
        "tables_dir": str(TAB),
    },
}
with open(OUT / "track26_summary.json", "w") as fh:
    json.dump(summary, fh, indent=2, default=str)

print("[11/12] writing per-sample combined frame")
df.to_csv(TAB / "T10_per_sample_combined.tsv", sep="\t", index=False)
print("[12/12] done.")
