#!/usr/bin/env python3
"""
Track 28 - DM1-high x IFN-gamma-low "fail-to-induce" THCA subgroup.

Pan-cancer Track 6' signal (DM1-high x HLA-I-low BRCA OS HR=2.20) failed to
replicate in TCGA-THCA because DM1 x HLA-I is positively correlated under
inflamed-induction biology. This track tests whether a residual subgroup of
THCA samples that are DM1-high but FAIL to induce the canonical inflamed
pattern (low IFN-gamma) exists, and if so, whether they look biologically
distinct (driver, methylation, RAI-dediff axis) or clinically distinct (worse
OS / DSS / PFI).

Boundary: HLA-I/II are gene-expression modules ONLY. Per
project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md, no allele
genotyping is performed on cancer cohorts in this track. Caption boilerplate
is appended to every figure.
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
import seaborn as sns
import statsmodels.api as sm
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/hla_deepdive_2026_05_08/track28_immune_cold_thca"
FIG = OUT / "figs"
TAB = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

CAPTION = (
    "HLA-I/II = gene-expression module (transcript level) only; "
    "not allele genotype. Cancer-cohort allele genotyping out of scope."
)

# -----------------------------------------------------------------------------
# 0. Inputs
# -----------------------------------------------------------------------------
TRACK26_T10 = (ROOT / "project/results/hla_deepdive_2026_05_08/"
               "track26_ifng_hla_dm1/tables/T10_per_sample_combined.tsv")
TERT_MASTER = (ROOT / "project/results/v17_tert_recovery/v2/"
               "sample_master_v17_tert_v2.tsv")
PANCAN_SURV = ROOT / "project/data/raw/TCGA_pancan/survival.tsv"
METH_FILE = (ROOT / "project/results/audit_2026_04_30/round5/"
             "r5_2_sample_methylation_8gene.tsv")
CLIN_EXT = ROOT / "project/results/tables/tcga_thca_clinical_extended.tsv"
PANCAN_EXP = ROOT / "project/data/raw/TCGA_pancan/pancan_geneExp.gz"
PANCAN_PHENO = ROOT / "project/data/raw/TCGA_pancan/phenotype.tsv.gz"
RAI_DELTA = (ROOT / "project/results/paper3_ici_track_b_lite/"
             "gse151179_rai_after_vs_before.tsv")

# Marker-gene proxies for cell-type composition (used because no per-sample
# CIBERSORT for full TCGA-THCA in this repo). Each list is z-meaned.
CELLTYPE_MARKERS = {
    "CD8_T_cell": ["CD8A", "CD8B", "GZMB", "GZMK", "PRF1"],
    "CD4_T_helper": ["CD4", "IL7R", "CCR7"],
    "T_reg": ["FOXP3", "IL2RA", "IKZF2"],
    "B_cell": ["CD19", "CD20", "MS4A1", "CD79A", "CD79B"],
    "Macrophage_M1": ["CD68", "NOS2", "IL1B", "TNF"],
    "Macrophage_M2": ["CD163", "MRC1", "MS4A4A"],
    "NK_cell": ["NCAM1", "NCR1", "KLRD1", "NKG7"],
    "Dendritic": ["ITGAX", "CLEC9A", "BATF3"],
    "Stromal_fibroblast": ["FAP", "PDPN", "COL1A1"],
    "Endothelial": ["PECAM1", "VWF", "CDH5"],
    "Neutrophil": ["FCGR3B", "CXCR2", "S100A8", "S100A9"],
    "Cytolytic": ["GZMA", "PRF1"],
}

# RAI-dediff signature direction from GSE151179 post-RAI vs pre-RAI
# (memory v19_paper3_rai_dediff_axis_2026_05_06):
# thyroid_diff DOWN (d=-1.01), HLA-II UP, myeloid_suppressive UP. The
# composite score below is sign-aligned to "post-RAI dediff state".
# We re-build it from the same module families.
THYROID_DIFF_GENES = ["TG", "TPO", "TSHR", "SLC5A5", "DIO1", "DIO2",
                      "FOXE1", "PAX8", "NKX2-1", "DUOX1", "DUOX2", "IYD"]
MYELOID_SUPP_GENES = ["CD163", "MRC1", "ARG1", "IL10", "VSIG4", "MS4A4A",
                      "CD68", "CSF1R"]
HLA2_GENES = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1",
              "HLA-DQB1", "HLA-DMA", "HLA-DMB", "CIITA", "CD74"]

# Hallmark gene sets for GSEA-lite (precomputed canonical short lists).
# We do not pull from MSigDB here; we use compact, well-annotated proxy sets.
HALLMARK_PROXY = {
    "INTERFERON_GAMMA_RESPONSE": [
        "STAT1", "IRF1", "IRF7", "IRF8", "GBP1", "GBP4", "GBP5",
        "CXCL9", "CXCL10", "CXCL11", "MX1", "MX2", "OAS1", "OAS2", "OAS3",
        "ISG15", "IFI6", "IFI27", "IFIT1", "IFIT3", "TAP1", "TAP2",
        "PSMB8", "PSMB9", "B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E",
        "NLRC5", "CIITA",
    ],
    "INTERFERON_ALPHA_RESPONSE": [
        "ISG15", "IFI6", "IFI27", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2",
        "OAS1", "OAS2", "OAS3", "RSAD2", "IFI44", "IFI44L", "DDX60",
        "USP18", "STAT1",
    ],
    "TNFA_SIGNALING_VIA_NFKB": [
        "TNF", "TNFAIP3", "ICAM1", "VCAM1", "NFKB1", "NFKB2", "RELA",
        "RELB", "BIRC3", "CXCL1", "CXCL2", "CXCL8", "IL6", "PTGS2",
        "ATF3", "EGR1", "FOS", "JUNB", "JUN",
    ],
    "INFLAMMATORY_RESPONSE": [
        "IL6", "IL1B", "CXCL8", "CCL2", "CCL5", "CXCL9", "CXCL10",
        "TLR2", "TLR4", "MYD88", "NFKB1", "PTGS2", "NLRP3",
    ],
    "EPITHELIAL_MESENCHYMAL_TRANSITION": [
        "VIM", "FN1", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1",
        "CDH2", "MMP2", "MMP9", "COL1A1", "COL1A2", "COL3A1", "FAP",
    ],
    "OXIDATIVE_PHOSPHORYLATION": [
        "ATP5A1", "ATP5B", "ATP5C1", "COX4I1", "COX5A", "COX5B",
        "NDUFA1", "NDUFA2", "NDUFA4", "NDUFB1", "NDUFB2", "SDHA", "SDHB",
        "UQCRB", "UQCRC1", "UQCRC2",
    ],
    "FATTY_ACID_METABOLISM": [
        "ACAA1", "ACAA2", "ACADM", "ACADL", "ACADVL", "CPT1A", "CPT2",
        "HADHA", "HADHB", "ACSL1", "ACSL3", "ACSL4", "FASN", "ACACA",
    ],
    "MYC_TARGETS_V1": [
        "MYC", "ODC1", "EIF4E", "PA2G4", "NPM1", "EIF2S1", "EIF3K",
        "EIF4A1", "RPL3", "RPL4", "RPS5", "RPS9", "NOP56", "PPRC1",
    ],
    "G2M_CHECKPOINT": [
        "CDK1", "CDK2", "CCNB1", "CCNB2", "CCNA2", "AURKA", "AURKB",
        "MKI67", "TOP2A", "PLK1", "BIRC5", "BUB1", "MCM2", "MCM4",
    ],
    "E2F_TARGETS": [
        "E2F1", "E2F2", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7",
        "PCNA", "CDC6", "CDC20", "CCNE1", "CCNE2",
    ],
    "P53_PATHWAY": [
        "TP53", "CDKN1A", "MDM2", "BAX", "BBC3", "PMAIP1", "ZMAT3",
        "GADD45A", "FAS", "DDB2",
    ],
    "ANGIOGENESIS": [
        "VEGFA", "VEGFB", "VEGFC", "FLT1", "KDR", "FLT4", "PECAM1",
        "VWF", "ANGPT1", "ANGPT2", "TEK", "PDGFA", "PDGFB",
    ],
    "ALLOGRAFT_REJECTION": [
        "CD3D", "CD3E", "CD8A", "CD8B", "GZMB", "PRF1", "IFNG", "IL2",
        "HLA-A", "HLA-B", "HLA-DRA", "CIITA", "TAP1", "TAP2",
    ],
}

# -----------------------------------------------------------------------------
# 1. Load the master Track 26 frame (per-sample modules, DM1, IFN-g, HLA, etc.)
# -----------------------------------------------------------------------------
print("[1/12] loading Track 26 per-sample combined frame")
df = pd.read_csv(TRACK26_T10, sep="\t")
print(f"  Track 26 frame: {df.shape}")
# unify TCGA short id (drop -01 suffix)
df["sample_short"] = df["sample"].str[:12]

# -----------------------------------------------------------------------------
# 2. Quadrant definition (DM1 median x IFN-gamma median)
# -----------------------------------------------------------------------------
print("[2/12] defining 4-way quadrant")
dm1_med = df["DM1_use"].median()
ifng_med = df["IFNG_hallmark"].median()
df["DM1_high"] = df["DM1_use"] > dm1_med
df["IFNG_high"] = df["IFNG_hallmark"] > ifng_med

def assign_quad(r):
    if r["DM1_high"] and r["IFNG_high"]:
        return "Q_inflamed_dark"        # DM1-high, IFNg-high (canonical)
    if r["DM1_high"] and not r["IFNG_high"]:
        return "Q_fail_to_induce"       # DM1-high, IFNg-low (cold dark matter)
    if (not r["DM1_high"]) and r["IFNG_high"]:
        return "Q_inflamed_bright"      # DM1-low, IFNg-high
    return "Q_quiet_bright"             # DM1-low, IFNg-low

df["quadrant"] = df.apply(assign_quad, axis=1)
quad_counts = df["quadrant"].value_counts().rename_axis("quadrant").reset_index(name="n")
quad_counts["pct"] = (quad_counts["n"] / len(df) * 100).round(2)
quad_counts.to_csv(TAB / "T01_quadrant_counts.tsv", sep="\t", index=False)
print(quad_counts.to_string(index=False))

# Also a tertile-based 3x3 view for figure 1
df["DM1_t"] = pd.qcut(df["DM1_use"], 3, labels=["DM1_low", "DM1_mid", "DM1_hi"])
df["IFNG_t"] = pd.qcut(df["IFNG_hallmark"], 3, labels=["IFNG_low", "IFNG_mid", "IFNG_hi"])
tertile_table = pd.crosstab(df["DM1_t"], df["IFNG_t"])
tertile_table.to_csv(TAB / "T02_dm1_x_ifng_tertile_counts.tsv", sep="\t")

# -----------------------------------------------------------------------------
# 3. Driver enrichment per quadrant (Fisher exact for each driver vs all-other)
# -----------------------------------------------------------------------------
print("[3/12] driver enrichment per quadrant")
df["driver"] = df["driver_simple"].fillna("Unknown")
fish_rows = []
quads = ["Q_fail_to_induce", "Q_inflamed_dark", "Q_inflamed_bright", "Q_quiet_bright"]
drivers = sorted(df["driver"].unique())
for q in quads:
    for drv in drivers:
        a = ((df["quadrant"] == q) & (df["driver"] == drv)).sum()
        b = ((df["quadrant"] == q) & (df["driver"] != drv)).sum()
        c = ((df["quadrant"] != q) & (df["driver"] == drv)).sum()
        d = ((df["quadrant"] != q) & (df["driver"] != drv)).sum()
        if a + b == 0 or c + d == 0:
            continue
        odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
        fish_rows.append({"quadrant": q, "driver": drv, "n_in_quad_with_driver": a,
                          "n_in_quad_other_driver": b,
                          "n_outside_quad_with_driver": c,
                          "n_outside_quad_other_driver": d,
                          "OR": odds, "p": p})
fish = pd.DataFrame(fish_rows)
fish["fdr"] = multipletests(fish["p"], method="fdr_bh")[1]
fish.to_csv(TAB / "T03_driver_enrichment_per_quadrant.tsv", sep="\t", index=False)
print("  driver enrichment top hits:")
print(fish.sort_values("p").head(8).to_string(index=False))

# -----------------------------------------------------------------------------
# 4. Clinical features per quadrant
# -----------------------------------------------------------------------------
print("[4/12] clinical features per quadrant")
clin = pd.read_csv(CLIN_EXT, sep="\t")
clin["sample_short"] = clin["sample_id"].str[:12]
clin_keep = clin[["sample_short", "stage", "age_at_diagnosis", "gender",
                  "tumor_size_mm"]].drop_duplicates(subset="sample_short")
df = df.merge(clin_keep, on="sample_short", how="left", suffixes=("", "_ext"))

# Pull stage + AJCC / hist subtype + sex / age from TERT master too
tm = pd.read_csv(TERT_MASTER, sep="\t")
tm["sample_short"] = tm["sample_id"].str[:12]
tm_keep = tm[["sample_short", "histology_subtype", "ajcc_stage_group",
              "ata_risk_proxy", "tert_promoter_integrated",
              "molecular_subtype", "tds_group", "dediff_flag",
              "rai_score_v17", "tds16_score_v17"]].drop_duplicates(
                  subset="sample_short")
df = df.merge(tm_keep, on="sample_short", how="left")

# tabulate categorical clinical features per quadrant
def cat_table(d, col):
    out = pd.crosstab(d["quadrant"], d[col].fillna("NA"))
    out["total"] = out.sum(axis=1)
    return out

clin_cats = {}
for col in ["stage", "gender", "histology_subtype",
            "molecular_subtype", "tds_group", "dediff_flag", "ata_risk_proxy"]:
    if col in df.columns:
        clin_cats[col] = cat_table(df, col)
        clin_cats[col].to_csv(TAB / f"T04_clin_cat_{col}.tsv", sep="\t")

# continuous: age_at_diagnosis, tumor_size, RAI/TDS proxies
cont_rows = []
for col in ["age_at_diagnosis", "tumor_size_mm", "rai_score_v17",
            "tds16_score_v17", "DM1_use", "IFNG_hallmark", "IFNG_ayers",
            "HLA1_score", "HLA2_score"]:
    if col not in df.columns:
        continue
    grp = df.groupby("quadrant")[col].agg(["count", "mean", "median", "std"])
    grp["feature"] = col
    cont_rows.append(grp.reset_index())
cont = pd.concat(cont_rows, ignore_index=True) if cont_rows else pd.DataFrame()
cont.to_csv(TAB / "T05_continuous_per_quadrant.tsv", sep="\t", index=False)

# Fail-to-induce vs inflamed-dark continuous tests (Mann-Whitney)
ft = df[df["quadrant"] == "Q_fail_to_induce"]
ind = df[df["quadrant"] == "Q_inflamed_dark"]
test_rows = []
for col in ["age_at_diagnosis", "tumor_size_mm", "rai_score_v17",
            "tds16_score_v17", "DM1_use", "IFNG_hallmark",
            "HLA1_score", "HLA2_score"]:
    if col not in df.columns:
        continue
    a = ft[col].dropna()
    b = ind[col].dropna()
    if len(a) < 5 or len(b) < 5:
        continue
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    pooled = np.sqrt(((len(a) - 1) * a.std() ** 2 + (len(b) - 1) * b.std() ** 2)
                     / (len(a) + len(b) - 2)) if (len(a) + len(b) - 2) > 0 else np.nan
    d = (a.mean() - b.mean()) / pooled if pooled and pooled > 0 else np.nan
    test_rows.append({"feature": col, "n_fail": len(a), "n_inflamed_dark": len(b),
                      "mean_fail": a.mean(), "mean_inflamed_dark": b.mean(),
                      "cohens_d": d, "MW_U": u, "MW_p": p})
ttab = pd.DataFrame(test_rows)
if not ttab.empty:
    ttab["fdr"] = multipletests(ttab["MW_p"], method="fdr_bh")[1]
ttab.to_csv(TAB / "T06_fail_vs_inflamed_dark_continuous.tsv", sep="\t", index=False)
print(ttab.to_string(index=False))

# -----------------------------------------------------------------------------
# 5. TERT promoter mutation per quadrant
# -----------------------------------------------------------------------------
print("[5/12] TERT promoter mutation per quadrant")
df["tert_mut"] = (df["tert_promoter_integrated"].fillna("wildtype") == "mutated").astype(int)
tert_tab = df.groupby("quadrant")["tert_mut"].agg(["sum", "count"])
tert_tab["pct"] = (tert_tab["sum"] / tert_tab["count"] * 100).round(2)
tert_tab.to_csv(TAB / "T07_tert_per_quadrant.tsv", sep="\t")

# Fisher fail-to-induce vs all others
def fisher_quad(d, q, col):
    a = ((d["quadrant"] == q) & (d[col] == 1)).sum()
    b = ((d["quadrant"] == q) & (d[col] == 0)).sum()
    c = ((d["quadrant"] != q) & (d[col] == 1)).sum()
    e = ((d["quadrant"] != q) & (d[col] == 0)).sum()
    odds, p = stats.fisher_exact([[a, b], [c, e]])
    return {"quadrant": q, "col": col, "a": a, "b": b, "c": c, "d": e,
            "OR": odds, "p": p}

tert_fish = pd.DataFrame([fisher_quad(df, q, "tert_mut") for q in quads])
tert_fish.to_csv(TAB / "T08_tert_fisher_per_quadrant.tsv", sep="\t", index=False)
print(tert_fish.to_string(index=False))

# -----------------------------------------------------------------------------
# 6. Methylation per quadrant
# -----------------------------------------------------------------------------
print("[6/12] methylation per quadrant")
meth = pd.read_csv(METH_FILE, sep="\t")
meth_cols = [c for c in meth.columns if c not in ("sample_short",)]
df_meth = df.merge(meth, on="sample_short", how="left")
meth_summary = df_meth.groupby("quadrant")[meth_cols].agg(["mean", "median", "count"])
meth_summary.to_csv(TAB / "T09_methylation_per_quadrant.tsv", sep="\t")

# -----------------------------------------------------------------------------
# 7. Cell-type composition + cytolytic + neoantigen-proxy + extra HLA modules
# -----------------------------------------------------------------------------
print("[7/12] streaming TCGA-THCA expression for cell-type marker modules + GSEA")

needed_genes = set()
for g in CELLTYPE_MARKERS.values():
    needed_genes.update(g)
for g in HALLMARK_PROXY.values():
    needed_genes.update(g)
needed_genes.update(THYROID_DIFF_GENES)
needed_genes.update(MYELOID_SUPP_GENES)
needed_genes.update(HLA2_GENES)
# neoantigen burden proxy: APM components + NLRC5 / IRF1 (Rooney 2015 cytolytic)
needed_genes.update(["GZMA", "PRF1", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9",
                     "NLRC5", "IRF1"])

pheno = pd.read_csv(PANCAN_PHENO, sep="\t")
thy_primary = set(pheno[(pheno["_primary_disease"] == "thyroid carcinoma") &
                        (pheno["sample_type"] == "Primary Tumor")]["sample"])

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
expr = pd.DataFrame(rows, index=keep_samples)
print(f"  expr shape: {expr.shape}")


def zmean(d, gl):
    sub = d.reindex(columns=[g for g in gl if g in d.columns])
    if sub.shape[1] == 0:
        return pd.Series(np.nan, index=d.index)
    z = (sub - sub.mean()) / sub.std(ddof=0)
    return z.mean(axis=1)


# Cell-type proxies
ct_scores = pd.DataFrame({k: zmean(expr, v) for k, v in CELLTYPE_MARKERS.items()})
ct_scores["sample"] = ct_scores.index
df = df.merge(ct_scores, on="sample", how="left")

# Hallmark scores
hm_scores = pd.DataFrame({k: zmean(expr, v) for k, v in HALLMARK_PROXY.items()})
hm_scores["sample"] = hm_scores.index
df = df.merge(hm_scores, on="sample", how="left")

# RAI-dediff axis composite (sign-aligned to post-RAI dediff state)
thyroid_diff_z = zmean(expr, THYROID_DIFF_GENES)
myeloid_supp_z = zmean(expr, MYELOID_SUPP_GENES)
hla2_z = zmean(expr, HLA2_GENES)
rai_dediff = (-thyroid_diff_z + 0.5 * myeloid_supp_z + 0.5 * hla2_z)
rai_score_df = pd.DataFrame({
    "sample": rai_dediff.index,
    "thyroid_diff_score": thyroid_diff_z.values,
    "myeloid_supp_score": myeloid_supp_z.values,
    "hla2_recomp_score": hla2_z.values,
    "RAI_dediff_score": rai_dediff.values,
})
df = df.merge(rai_score_df, on="sample", how="left")

# Cytolytic activity (Rooney 2015 - GZMA + PRF1 log mean)
df["cyt_activity"] = df["Cytolytic"]

# Neoantigen-presentation proxy (APM components z-mean)
apm_genes = ["B2M", "TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5", "IRF1"]
df["APM_score"] = zmean(expr, apm_genes).reindex(df["sample"]).values

# Per-quadrant summary of immune + RAI + APM + cell composition
imm_cols = (list(CELLTYPE_MARKERS.keys()) + list(HALLMARK_PROXY.keys())
            + ["RAI_dediff_score", "thyroid_diff_score", "myeloid_supp_score",
               "hla2_recomp_score", "cyt_activity", "APM_score"])
imm_summary = df.groupby("quadrant")[imm_cols].mean().round(3)
imm_summary.to_csv(TAB / "T10_immune_celltype_per_quadrant.tsv", sep="\t")

# Compare fail-to-induce vs inflamed-dark for these features
ft = df[df["quadrant"] == "Q_fail_to_induce"]
ind = df[df["quadrant"] == "Q_inflamed_dark"]
imm_test_rows = []
for col in imm_cols:
    a = ft[col].dropna()
    b = ind[col].dropna()
    if len(a) < 5 or len(b) < 5:
        continue
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    pooled = np.sqrt(((len(a) - 1) * a.std() ** 2 + (len(b) - 1) * b.std() ** 2)
                     / (len(a) + len(b) - 2))
    d = (a.mean() - b.mean()) / pooled if pooled > 0 else np.nan
    imm_test_rows.append({"feature": col, "n_fail": len(a), "n_inflamed_dark": len(b),
                          "mean_fail": a.mean(), "mean_inflamed_dark": b.mean(),
                          "cohens_d_fail_minus_inflamed": d, "MW_p": p})
imm_t = pd.DataFrame(imm_test_rows)
imm_t["fdr"] = multipletests(imm_t["MW_p"], method="fdr_bh")[1]
imm_t = imm_t.sort_values("MW_p")
imm_t.to_csv(TAB / "T11_immune_fail_vs_inflamed_dark.tsv", sep="\t", index=False)
print("  fail-to-induce vs inflamed-dark top features:")
print(imm_t.head(10).to_string(index=False))

# -----------------------------------------------------------------------------
# 8. Survival per quadrant (OS/DSS/PFI/DFI), with Cox HR
# -----------------------------------------------------------------------------
print("[8/12] survival per quadrant")
surv = pd.read_csv(PANCAN_SURV, sep="\t")
surv_thca = surv[surv["cancer type abbreviation"] == "THCA"][[
    "sample", "OS", "OS.time", "DSS", "DSS.time", "DFI", "DFI.time",
    "PFI", "PFI.time"]].copy()
# df already has OS/OS.time/PFI/PFI.time from Track 5 frame; drop them so
# merge brings the canonical pancan survival columns and adds DSS/DFI.
for col in ["OS", "OS.time", "PFI", "PFI.time"]:
    if col in df.columns:
        df = df.drop(columns=[col])
df_s = df.merge(surv_thca, on="sample", how="left")

# Cox HR per quadrant (vs Q_inflamed_dark reference) for each endpoint
cox_rows = []
for endp, time_col, ev_col in [
        ("OS", "OS.time", "OS"), ("DSS", "DSS.time", "DSS"),
        ("PFI", "PFI.time", "PFI"), ("DFI", "DFI.time", "DFI")]:
    sub = df_s[["sample", "quadrant", time_col, ev_col,
                "age_at_diagnosis", "stage", "driver"]].dropna(
                    subset=[time_col, ev_col])
    sub = sub[sub[time_col] > 0]
    if sub[ev_col].sum() < 3:
        continue
    sub = sub.copy()
    sub["q"] = sub["quadrant"].astype("category")
    sub["q"] = sub["q"].cat.set_categories(
        ["Q_inflamed_dark", "Q_fail_to_induce",
         "Q_inflamed_bright", "Q_quiet_bright"], ordered=False)
    # Univariate Cox
    cph = CoxPHFitter()
    Xq = pd.get_dummies(sub["q"], drop_first=True).astype(float)
    Xq[time_col] = sub[time_col].values
    Xq[ev_col] = sub[ev_col].values
    try:
        cph.fit(Xq, duration_col=time_col, event_col=ev_col)
        for name, row in cph.summary.iterrows():
            cox_rows.append({"endpoint": endp, "model": "univariate",
                             "term": name, "HR": row["exp(coef)"],
                             "HR_low": row["exp(coef) lower 95%"],
                             "HR_high": row["exp(coef) upper 95%"],
                             "p": row["p"], "n": len(sub),
                             "events": int(sub[ev_col].sum())})
    except Exception as exc:
        print(f"  univariate Cox failed for {endp}: {exc}")

    # Multivariate adjusted for stage (high vs low), age, driver(BRAF y/n)
    sub["stage_high"] = sub["stage"].astype(str).str.contains(
        "III|IV", regex=True, na=False).astype(int)
    sub["age"] = pd.to_numeric(sub["age_at_diagnosis"], errors="coerce")
    sub["BRAF_yn"] = (sub["driver"] == "BRAF").astype(int)
    Xm = pd.concat([pd.get_dummies(sub["q"], drop_first=True),
                    sub[["stage_high", "age", "BRAF_yn"]]], axis=1).astype(float)
    Xm[time_col] = sub[time_col].values
    Xm[ev_col] = sub[ev_col].values
    Xm = Xm.dropna()
    if len(Xm) < 30 or Xm[ev_col].sum() < 3:
        continue
    cph_m = CoxPHFitter(penalizer=0.01)
    try:
        cph_m.fit(Xm, duration_col=time_col, event_col=ev_col)
        for name, row in cph_m.summary.iterrows():
            cox_rows.append({"endpoint": endp, "model": "multivariate_adj",
                             "term": name, "HR": row["exp(coef)"],
                             "HR_low": row["exp(coef) lower 95%"],
                             "HR_high": row["exp(coef) upper 95%"],
                             "p": row["p"], "n": len(Xm),
                             "events": int(Xm[ev_col].sum())})
    except Exception as exc:
        print(f"  multivariate Cox failed for {endp}: {exc}")

cox_tab = pd.DataFrame(cox_rows)
cox_tab.to_csv(TAB / "T12_cox_per_quadrant.tsv", sep="\t", index=False)
print(cox_tab.head(20).to_string(index=False))

# Multivariate logrank
lr_rows = []
for endp, time_col, ev_col in [
        ("OS", "OS.time", "OS"), ("DSS", "DSS.time", "DSS"),
        ("PFI", "PFI.time", "PFI"), ("DFI", "DFI.time", "DFI")]:
    sub = df_s.dropna(subset=[time_col, ev_col])
    sub = sub[sub[time_col] > 0]
    if sub[ev_col].sum() < 3:
        continue
    res = multivariate_logrank_test(
        sub[time_col], sub["quadrant"], sub[ev_col])
    lr_rows.append({"endpoint": endp, "test_stat": res.test_statistic,
                    "p": res.p_value, "n": len(sub),
                    "events": int(sub[ev_col].sum())})
lr_tab = pd.DataFrame(lr_rows)
lr_tab.to_csv(TAB / "T13_logrank_quadrant.tsv", sep="\t", index=False)
print(lr_tab.to_string(index=False))

# -----------------------------------------------------------------------------
# 9. Gene-level DEG: fail-to-induce vs inflamed-dark
# -----------------------------------------------------------------------------
print("[9/12] DEG fail-to-induce vs inflamed-dark — streaming whole expression")

ft_samples = set(df[df["quadrant"] == "Q_fail_to_induce"]["sample"])
ind_samples = set(df[df["quadrant"] == "Q_inflamed_dark"]["sample"])
print(f"  n fail: {len(ft_samples)}, n inflamed_dark: {len(ind_samples)}")

# stream expression file again, retain all genes for THCA primary samples
with gzip.open(PANCAN_EXP, "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    sample_cols = header[1:]
    keep_idx_ft = [i for i, s in enumerate(sample_cols) if s in ft_samples]
    keep_idx_ind = [i for i, s in enumerate(sample_cols) if s in ind_samples]
    deg_rows = []
    for ln in fh:
        parts = ln.rstrip("\n").split("\t")
        gene = parts[0]
        try:
            ft_v = np.array([parts[i + 1] for i in keep_idx_ft], dtype=float)
            in_v = np.array([parts[i + 1] for i in keep_idx_ind], dtype=float)
        except ValueError:
            continue
        ft_v = ft_v[~np.isnan(ft_v)]
        in_v = in_v[~np.isnan(in_v)]
        if len(ft_v) < 5 or len(in_v) < 5:
            continue
        if ft_v.std() == 0 and in_v.std() == 0:
            continue
        # Welch's t-test on log-normalised pancan values (already log2-norm)
        t, p = stats.ttest_ind(ft_v, in_v, equal_var=False)
        pooled_sd = np.sqrt(
            ((len(ft_v) - 1) * ft_v.std(ddof=1) ** 2 +
             (len(in_v) - 1) * in_v.std(ddof=1) ** 2) /
            (len(ft_v) + len(in_v) - 2))
        d = (ft_v.mean() - in_v.mean()) / pooled_sd if pooled_sd > 0 else np.nan
        deg_rows.append({"gene": gene, "n_fail": len(ft_v), "n_ind": len(in_v),
                         "mean_fail": ft_v.mean(), "mean_inflamed_dark": in_v.mean(),
                         "delta_mean": ft_v.mean() - in_v.mean(),
                         "cohens_d": d, "t": t, "p": p})
deg = pd.DataFrame(deg_rows)
deg["fdr"] = multipletests(deg["p"].fillna(1.0), method="fdr_bh")[1]
deg = deg.sort_values("p")
deg.to_csv(TAB / "T14_deg_fail_vs_inflamed_dark.tsv.gz",
           sep="\t", index=False, compression="gzip")
print(f"  total genes tested: {len(deg)}")
print("  top 10 by p:")
print(deg.head(10).to_string(index=False))

# Save top 50 up + down for the report
top_up = deg[deg["delta_mean"] > 0].head(50)
top_dn = deg[deg["delta_mean"] < 0].head(50)
top_up.to_csv(TAB / "T15_deg_top50_up_in_fail.tsv", sep="\t", index=False)
top_dn.to_csv(TAB / "T16_deg_top50_dn_in_fail.tsv", sep="\t", index=False)

# Hallmark-proxy GSEA (one-sided rank-sum on log2FC ranks)
deg_clean = deg.dropna(subset=["t"]).copy()
deg_clean = deg_clean.set_index("gene")
gsea_rows = []
all_t = deg_clean["t"].values
for hm_name, gl in HALLMARK_PROXY.items():
    in_set = deg_clean.index.isin(gl)
    if in_set.sum() < 3:
        continue
    in_t = deg_clean.loc[in_set, "t"].values
    out_t = deg_clean.loc[~in_set, "t"].values
    # Mann-Whitney on signed t-stat: positive means hallmark genes shifted UP
    # in fail-to-induce relative to inflamed-dark
    u, p = stats.mannwhitneyu(in_t, out_t, alternative="two-sided")
    sign = "up_in_fail" if in_t.mean() > out_t.mean() else "down_in_fail"
    gsea_rows.append({"hallmark": hm_name, "n_genes_in_set": int(in_set.sum()),
                      "mean_t_in_set": in_t.mean(),
                      "mean_t_outside": out_t.mean(),
                      "direction": sign, "MW_U": u, "p": p})
gsea = pd.DataFrame(gsea_rows)
gsea["fdr"] = multipletests(gsea["p"], method="fdr_bh")[1]
gsea = gsea.sort_values("p")
gsea.to_csv(TAB / "T17_hallmark_gsea_fail_vs_inflamed_dark.tsv",
            sep="\t", index=False)
print("  Hallmark GSEA-lite (fail-to-induce vs inflamed-dark):")
print(gsea.to_string(index=False))

# -----------------------------------------------------------------------------
# 10. RAI-refractoriness signature compare across quadrants (extra view)
# -----------------------------------------------------------------------------
print("[10/12] RAI-dediff signature across quadrants")
rai_summary = df.groupby("quadrant")[
    ["RAI_dediff_score", "thyroid_diff_score", "myeloid_supp_score",
     "hla2_recomp_score"]].agg(["mean", "median", "count"])
rai_summary.to_csv(TAB / "T18_rai_dediff_per_quadrant.tsv", sep="\t")

# fail-to-induce vs inflamed-dark on RAI score
a = df[df["quadrant"] == "Q_fail_to_induce"]["RAI_dediff_score"].dropna()
b = df[df["quadrant"] == "Q_inflamed_dark"]["RAI_dediff_score"].dropna()
u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
pooled = np.sqrt(((len(a) - 1) * a.std() ** 2 + (len(b) - 1) * b.std() ** 2)
                 / (len(a) + len(b) - 2))
rai_d = (a.mean() - b.mean()) / pooled if pooled > 0 else np.nan
rai_test = pd.DataFrame([{
    "feature": "RAI_dediff_score", "n_fail": len(a), "n_inflamed_dark": len(b),
    "mean_fail": a.mean(), "mean_inflamed_dark": b.mean(),
    "cohens_d_fail_minus_inflamed": rai_d, "MW_p": p,
}])
rai_test.to_csv(TAB / "T19_rai_dediff_fail_vs_inflamed_dark.tsv",
                sep="\t", index=False)
print(rai_test.to_string(index=False))

# -----------------------------------------------------------------------------
# 11. Per-sample export
# -----------------------------------------------------------------------------
df_export_cols = [
    "sample", "sample_short", "DM1_use", "IFNG_hallmark", "IFNG_ayers",
    "HLA1_score", "HLA2_score", "DM1_high", "IFNG_high", "quadrant",
    "driver", "tert_promoter_integrated", "tert_mut",
    "stage", "age_at_diagnosis", "gender", "histology_subtype",
    "molecular_subtype", "tds_group", "dediff_flag",
    "rai_score_v17", "tds16_score_v17",
    "RAI_dediff_score", "thyroid_diff_score", "myeloid_supp_score",
    "hla2_recomp_score", "cyt_activity", "APM_score",
] + list(CELLTYPE_MARKERS.keys()) + list(HALLMARK_PROXY.keys())
df_export_cols = [c for c in df_export_cols if c in df.columns]
df[df_export_cols].to_csv(TAB / "T20_per_sample_full.tsv",
                          sep="\t", index=False)

# =============================================================================
# FIGURES
# =============================================================================
print("[11/12] figures")
sns.set_style("white")
QCOL = {
    "Q_inflamed_dark": "#e34a33",   # red
    "Q_fail_to_induce": "#2c7fb8",  # blue (highlight subgroup)
    "Q_inflamed_bright": "#fdae61",
    "Q_quiet_bright": "#bdbdbd",
}
QORDER = ["Q_quiet_bright", "Q_inflamed_bright", "Q_inflamed_dark",
          "Q_fail_to_induce"]


def add_caption(fig):
    fig.text(0.5, 0.005, CAPTION, ha="center", fontsize=7, color="#666")


# Figure 1: 4-quadrant scatter DM1 vs IFN-gamma
fig, ax = plt.subplots(figsize=(7, 6))
for q in QORDER:
    sub = df[df["quadrant"] == q]
    ax.scatter(sub["DM1_use"], sub["IFNG_hallmark"], s=14,
               c=QCOL[q], alpha=0.7, edgecolor="none",
               label=f"{q} (n={len(sub)})")
ax.axvline(dm1_med, color="#888", lw=0.7, ls="--")
ax.axhline(ifng_med, color="#888", lw=0.7, ls="--")
ax.set_xlabel("DM1 score (Track 5)")
ax.set_ylabel("IFN-gamma Hallmark z-mean")
ax.set_title("TCGA-THCA: DM1 x IFN-gamma quadrants")
ax.legend(fontsize=8, loc="lower right")
add_caption(fig)
fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(FIG / "F01_quadrant_scatter.png", dpi=180)
fig.savefig(FIG / "F01_quadrant_scatter.pdf")
plt.close(fig)

# Figure 2: Driver bar chart per quadrant
driver_share = (pd.crosstab(df["quadrant"], df["driver"], normalize="index") * 100)
driver_share = driver_share.reindex(QORDER)
fig, ax = plt.subplots(figsize=(8, 5))
driver_share.plot(kind="bar", stacked=True, ax=ax, colormap="tab10")
ax.set_ylabel("% of quadrant")
ax.set_xlabel("")
ax.set_title("Driver composition per quadrant (%, TCGA-THCA primary)")
ax.legend(title="driver", fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
plt.xticks(rotation=20, ha="right")
add_caption(fig)
fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(FIG / "F02_driver_per_quadrant.png", dpi=180)
fig.savefig(FIG / "F02_driver_per_quadrant.pdf")
plt.close(fig)

# Figure 3: TERT promoter mutation rate per quadrant
fig, ax = plt.subplots(figsize=(6, 4))
tert_rate = df.groupby("quadrant")["tert_mut"].mean().reindex(QORDER) * 100
tert_n = df.groupby("quadrant")["tert_mut"].sum().reindex(QORDER)
bars = ax.bar(range(len(QORDER)), tert_rate.values,
              color=[QCOL[q] for q in QORDER])
for i, (rate, n) in enumerate(zip(tert_rate.values, tert_n.values)):
    ax.text(i, rate + 0.5, f"{rate:.1f}%\n(n={int(n)})",
            ha="center", fontsize=8)
ax.set_xticks(range(len(QORDER)))
ax.set_xticklabels(QORDER, rotation=20, ha="right")
ax.set_ylabel("TERT promoter mutation rate (%)")
ax.set_title("TERT promoter mutation per quadrant")
add_caption(fig)
fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(FIG / "F03_tert_per_quadrant.png", dpi=180)
fig.savefig(FIG / "F03_tert_per_quadrant.pdf")
plt.close(fig)

# Figure 4: Boxplots HLA1, HLA2, cyt_activity, APM_score per quadrant
fig, axes = plt.subplots(1, 4, figsize=(15, 4.5), sharey=False)
for ax, col in zip(axes, ["HLA1_score", "HLA2_score", "cyt_activity", "APM_score"]):
    parts = [df[df["quadrant"] == q][col].dropna().values for q in QORDER]
    bp = ax.boxplot(parts, patch_artist=True,
                    boxprops=dict(linewidth=0.5),
                    medianprops=dict(color="k", lw=1.0))
    for patch, q in zip(bp["boxes"], QORDER):
        patch.set_facecolor(QCOL[q])
        patch.set_alpha(0.7)
    ax.set_xticklabels([q.replace("Q_", "") for q in QORDER],
                       rotation=20, ha="right", fontsize=8)
    ax.set_title(col, fontsize=9)
fig.suptitle("HLA / cytolytic / APM modules per quadrant", fontsize=11)
add_caption(fig)
fig.tight_layout(rect=[0, 0.03, 1, 0.96])
fig.savefig(FIG / "F04_hla_cyt_apm_per_quadrant.png", dpi=180)
fig.savefig(FIG / "F04_hla_cyt_apm_per_quadrant.pdf")
plt.close(fig)

# Figure 5: Methylation mean_8g per quadrant
if "mean_8g_beta" in df_meth.columns:
    fig, ax = plt.subplots(figsize=(6, 4))
    parts = [df_meth[df_meth["quadrant"] == q]["mean_8g_beta"].dropna().values
             for q in QORDER]
    bp = ax.boxplot(parts, patch_artist=True,
                    medianprops=dict(color="k"))
    for patch, q in zip(bp["boxes"], QORDER):
        patch.set_facecolor(QCOL[q])
        patch.set_alpha(0.7)
    ax.set_xticklabels([q.replace("Q_", "") for q in QORDER],
                       rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("mean_8g_beta (HM450)")
    ax.set_title("8-gene RAI-lineage methylation per quadrant")
    add_caption(fig)
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.savefig(FIG / "F05_meth_per_quadrant.png", dpi=180)
    fig.savefig(FIG / "F05_meth_per_quadrant.pdf")
    plt.close(fig)

# Figure 6: KM curves OS / DSS / PFI / DFI by quadrant
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
endpoints = [("OS", "OS.time", "OS"), ("DSS", "DSS.time", "DSS"),
             ("PFI", "PFI.time", "PFI"), ("DFI", "DFI.time", "DFI")]
for ax, (lbl, time_col, ev_col) in zip(axes.ravel(), endpoints):
    sub = df_s.dropna(subset=[time_col, ev_col])
    sub = sub[sub[time_col] > 0]
    for q in QORDER:
        ss = sub[sub["quadrant"] == q]
        if len(ss) < 5 or ss[ev_col].sum() < 1:
            continue
        kmf = KaplanMeierFitter()
        kmf.fit(ss[time_col] / 365.25, ss[ev_col],
                label=f"{q.replace('Q_', '')} n={len(ss)} ev={int(ss[ev_col].sum())}")
        kmf.plot_survival_function(ax=ax, ci_show=False, color=QCOL[q])
    ax.set_title(f"{lbl} by quadrant (TCGA-THCA)")
    ax.set_xlabel("years")
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=7)
add_caption(fig)
fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(FIG / "F06_km_quadrant.png", dpi=180)
fig.savefig(FIG / "F06_km_quadrant.pdf")
plt.close(fig)

# Figure 7: Cox HR forest (univariate) - non-reference quadrants vs Q_inflamed_dark
fig, ax = plt.subplots(figsize=(8, 6))
forest = cox_tab[cox_tab["model"] == "univariate"].copy()
forest = forest.sort_values(["endpoint", "term"])
y = np.arange(len(forest))
ax.errorbar(np.log(forest["HR"].values), y,
            xerr=[np.log(forest["HR"].values) - np.log(forest["HR_low"].values),
                  np.log(forest["HR_high"].values) - np.log(forest["HR"].values)],
            fmt="o", color="#444")
ax.axvline(0, color="#888", lw=0.7, ls="--")
labels = [f"{r['endpoint']} | {r['term']} | p={r['p']:.2g}"
          for _, r in forest.iterrows()]
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=7)
ax.set_xlabel("ln HR (vs Q_inflamed_dark)")
ax.set_title("Univariate Cox HR — quadrant terms vs inflamed_dark (THCA)")
add_caption(fig)
fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(FIG / "F07_cox_forest.png", dpi=180)
fig.savefig(FIG / "F07_cox_forest.pdf")
plt.close(fig)

# Figure 8: Cell-type composition heatmap (mean per quadrant, z-scored)
ct_mean = df.groupby("quadrant")[list(CELLTYPE_MARKERS.keys())].mean()
ct_mean = ct_mean.reindex(QORDER)
fig, ax = plt.subplots(figsize=(10, 4.5))
sns.heatmap(ct_mean, cmap="RdBu_r", center=0, ax=ax,
            cbar_kws={"label": "z-mean"}, annot=True, fmt=".2f",
            annot_kws={"size": 7})
ax.set_title("Cell-type marker module mean per quadrant")
add_caption(fig)
fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(FIG / "F08_celltype_heatmap.png", dpi=180)
fig.savefig(FIG / "F08_celltype_heatmap.pdf")
plt.close(fig)

# Figure 9: Hallmark GSEA bar (signed -log10 p)
fig, ax = plt.subplots(figsize=(8, 6))
gs = gsea.copy()
gs["sign_log_p"] = -np.log10(gs["p"]) * np.where(gs["direction"] == "up_in_fail",
                                                  1, -1)
gs = gs.sort_values("sign_log_p")
colors = ["#2c7fb8" if v > 0 else "#e34a33" for v in gs["sign_log_p"]]
ax.barh(gs["hallmark"], gs["sign_log_p"], color=colors)
ax.axvline(0, color="#888", lw=0.7)
ax.set_xlabel("signed -log10(p)  +  = up in fail-to-induce")
ax.set_title("Hallmark-proxy GSEA: fail-to-induce vs inflamed-dark")
add_caption(fig)
fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(FIG / "F09_hallmark_gsea.png", dpi=180)
fig.savefig(FIG / "F09_hallmark_gsea.pdf")
plt.close(fig)

# Figure 10: RAI-dediff score boxplot across quadrants
fig, ax = plt.subplots(figsize=(6, 4.5))
parts = [df[df["quadrant"] == q]["RAI_dediff_score"].dropna().values for q in QORDER]
bp = ax.boxplot(parts, patch_artist=True, medianprops=dict(color="k"))
for patch, q in zip(bp["boxes"], QORDER):
    patch.set_facecolor(QCOL[q])
    patch.set_alpha(0.7)
ax.set_xticklabels([q.replace("Q_", "") for q in QORDER],
                   rotation=20, ha="right", fontsize=8)
ax.set_ylabel("RAI-dediff composite (-thyroid_diff + 0.5*myeloid + 0.5*HLA-II)")
ax.set_title("RAI-dediff axis (post-RAI proxy) per quadrant")
add_caption(fig)
fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(FIG / "F10_rai_dediff_per_quadrant.png", dpi=180)
fig.savefig(FIG / "F10_rai_dediff_per_quadrant.pdf")
plt.close(fig)

# -----------------------------------------------------------------------------
# 12. Summary JSON
# -----------------------------------------------------------------------------
print("[12/12] writing summary JSON")
n_fail = int((df["quadrant"] == "Q_fail_to_induce").sum())
n_inflamed = int((df["quadrant"] == "Q_inflamed_dark").sum())
n_total = int(len(df))

# Driver enrichment for fail-to-induce
fish_fail = fish[fish["quadrant"] == "Q_fail_to_induce"].copy().sort_values("p")

# Top GSEA up/down
gsea_up_top = gsea[gsea["direction"] == "up_in_fail"].head(5).to_dict("records")
gsea_dn_top = gsea[gsea["direction"] == "down_in_fail"].head(5).to_dict("records")

cox_fail = cox_tab[
    (cox_tab["term"].astype(str).str.contains("fail_to_induce"))
].to_dict("records")

summary = {
    "boundary": CAPTION,
    "n_total_THCA": n_total,
    "n_fail_to_induce": n_fail,
    "n_inflamed_dark": n_inflamed,
    "DM1_median": dm1_med,
    "IFNG_median": ifng_med,
    "fail_to_induce_top_drivers": fish_fail.head(5).to_dict("records"),
    "TERT_per_quadrant": tert_tab.reset_index().to_dict("records"),
    "logrank_quadrant": lr_tab.to_dict("records"),
    "cox_fail_vs_inflamed_dark": cox_fail,
    "gsea_up_in_fail_top": gsea_up_top,
    "gsea_dn_in_fail_top": gsea_dn_top,
    "rai_dediff_fail_minus_inflamed_d": float(rai_d),
    "rai_dediff_p": float(p),
}

with open(OUT / "track28_summary.json", "w") as fh:
    json.dump(summary, fh, indent=2, default=str)

print("done.")
print(json.dumps(summary, indent=2, default=str)[:2500])
