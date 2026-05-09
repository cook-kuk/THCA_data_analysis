#!/usr/bin/env python3
"""
H19 — TMB / neoantigen / immune-evasion mechanisms in DM1 vs DM2 BRAF-cPTC.

Answers the WHY of DM1's immune-hot signal:
  (1) TMB / non-silent / mutation-class burden
  (2) HLA-I + HLA-II + antigen-presentation machinery (B2M / TAP1 / TAP2)
  (3) Immune-escape signature (B2M / HLA-I downregulation, JAK1/JAK2)
  (4) Neoantigen load (Thorsson 2018) + presentation-efficiency proxy
  (5) Korean K2 arcasHLA allele DM1-vs-DM2 enrichment (HT-priming bridge)
  (6) Immune-gene methylation / expression (CTLA4, PDCD1, FOXP3, IFNG)

Stratum: BRAF-driven, classical PTC histology (BRAF-cPTC), TCGA-THCA.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------------
# paths
# ---------------------------------------------------------------------------
ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/p2_braf_nature_sprint_2026_05_09/h19_tmb_neoantigen"
OUT.mkdir(parents=True, exist_ok=True)

CLIN = Path(
    "/data/thca/repo_results/p2_braf_nature_sprint_2026_05_09/h6_survival/h6_merged_clinical.tsv"
)
TRACK29 = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track29_neoantigen_hla/tables/T01_per_sample_merged.tsv"
)
TCGA_Z = Path("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv")
THORSSON_NEO = Path(
    "/home/seungho/personal/THCA_data_analysis/project/data/external/thorsson_2018/snv_neoantigens.tsv"
)
THORSSON_INDEL = Path(
    "/home/seungho/personal/THCA_data_analysis/project/data/external/thorsson_2018/indel_neoantigens.tsv"
)
K2_HLA = Path("/data/thca/_repo_offload/arcasHLA/K2_arcasHLA_FINAL.tsv")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def cohen_d(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    s = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    if s == 0:
        return np.nan
    return float((a.mean() - b.mean()) / s)


def wilcox(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    try:
        return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
    except Exception:
        return np.nan


def compare(df: pd.DataFrame, col: str, group: str = "dm") -> dict:
    """DM1 vs DM2 effect size + Wilcoxon for one column."""
    a = df.loc[df[group] == "DM1", col].dropna().values
    b = df.loc[df[group] == "DM2", col].dropna().values
    return {
        "feature": col,
        "n_DM1": len(a),
        "n_DM2": len(b),
        "mean_DM1": float(a.mean()) if len(a) else np.nan,
        "mean_DM2": float(b.mean()) if len(b) else np.nan,
        "cohen_d_DM1_vs_DM2": cohen_d(a, b),
        "wilcoxon_p": wilcox(a, b),
    }


# ---------------------------------------------------------------------------
# load + build BRAF-cPTC stratum
# ---------------------------------------------------------------------------
print("[H19] loading clinical + track29")
clin = pd.read_csv(CLIN, sep="\t")
clin["sample15"] = clin["sample_id"].str.slice(0, 15)
# 3-class label: DM1 / DM2 / not_DM. For mechanism comparison we use DM1 vs true-DM2.
# 'not_DM' samples are unclassified middle band — kept on disk but excluded from primary contrast.

t29 = pd.read_csv(TRACK29, sep="\t")
t29["sample15"] = t29["sample"].astype(str).str.slice(0, 15)

# inner-join on sample15
mer = clin.merge(
    t29.drop(columns=[c for c in t29.columns if c in clin.columns and c not in {"sample15"}], errors="ignore"),
    on="sample15", how="left",
)

# BRAF-cPTC stratum (full table for sample-level export)
braf_cptc_full = mer[(mer["molecular_subtype"] == "BRAF_like") & (mer["histology_subtype"] == "cPTC")].copy()
# Primary mechanism contrast: DM1 vs true-DM2 (drop 'not_DM' middle band)
braf_cptc = braf_cptc_full[braf_cptc_full["dm"].isin(["DM1", "DM2"])].copy()
print(f"  BRAF-cPTC total={len(braf_cptc_full)}; DM1-vs-DM2 contrast n={len(braf_cptc)} "
      f"(DM1={(braf_cptc.dm=='DM1').sum()}, DM2={(braf_cptc.dm=='DM2').sum()}, "
      f"not_DM excluded={(braf_cptc_full.dm=='not_DM').sum()})")

# track29 has Thorsson neoantigen / HLA TPM / CD274 etc already merged; report which columns are populated
key_cols_present = [c for c in
                    ["thorsson_nonsil_perMb", "thorsson_silent_perMb", "n_nonsyn_snv",
                     "n_immunogenic_mut", "n_binding_pMHC", "n_binding_expr_pMHC",
                     "n_indel", "n_immunogenic_indel", "n_neoantigen_indel",
                     "HLA_A_tpm", "HLA_B_tpm", "HLA_C_tpm", "B2M_tpm",
                     "HLA1_score", "HLA2_score", "leukocyte_frac", "cytolytic",
                     "CD274", "PDCD1", "CTLA4", "IDO1", "IFNG_hallmark",
                     "classI_nonsilent_count", "classI_lof_count", "nonsilent_mut", "total_mut", "tmb_a5"]
                    if c in braf_cptc.columns]


# ---------------------------------------------------------------------------
# (1) TMB
# ---------------------------------------------------------------------------
print("[H19] (1) TMB / mutation burden")
tmb_features = [
    "thorsson_nonsil_perMb", "thorsson_silent_perMb",
    "nonsilent_mut", "total_mut",
    "n_nonsyn_snv", "n_indel",
    "classI_nonsilent_count", "classI_lof_count",
]
tmb_features = [c for c in tmb_features if c in braf_cptc.columns]
tmb_rows = [compare(braf_cptc, c) for c in tmb_features]
tmb_df = pd.DataFrame(tmb_rows)
tmb_df.insert(0, "category", "TMB_mutation_burden")
tmb_df.to_csv(OUT / "h19_tmb_results.tsv", sep="\t", index=False)
print(tmb_df.to_string(index=False))


# ---------------------------------------------------------------------------
# (2) HLA expression preservation (need extra HLA-II / TAP1/TAP2 / JAK from z-score file)
# ---------------------------------------------------------------------------
print("[H19] (2) HLA expression preservation — pulling z-score for HLA-II + TAP + JAK")
hla_genes = [
    "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2",
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
    "JAK1", "JAK2", "STAT1", "IRF1",
    "CTLA4", "PDCD1", "FOXP3", "IFNG", "CD274", "IDO1",
    "GZMA", "PRF1", "CD8A", "CD8B", "CXCL9", "CXCL10", "CXCL13",
]
# stream-read z-score: rows=gene, cols=sample
import csv
zhead = pd.read_csv(TCGA_Z, sep="\t", nrows=0).columns.tolist()
z_samples = zhead[1:]
keep_idx = []
keep_genes = []
with open(TCGA_Z) as fh:
    r = csv.reader(fh, delimiter="\t")
    next(r)
    for row in r:
        if row[0] in hla_genes:
            keep_idx.append(row)
            keep_genes.append(row[0])
        if len(keep_genes) == len(hla_genes):
            break
zsel = pd.DataFrame(keep_idx, columns=zhead).set_index(zhead[0]).astype(float)
print(f"  loaded z-scores for {len(zsel)} genes × {zsel.shape[1]} samples")

# transpose to sample × gene
zsamp = zsel.T.reset_index().rename(columns={"index": "sample_z"})
zsamp["sample15"] = zsamp["sample_z"].astype(str).str.slice(0, 15)
# average duplicates
zsamp_avg = zsamp.drop(columns=["sample_z"]).groupby("sample15").mean().reset_index()

braf_cptc = braf_cptc.merge(zsamp_avg, on="sample15", how="left", suffixes=("", "_z"))

# composite HLA-I / HLA-II / antigen-presentation z scores
class_I = ["HLA-A", "HLA-B", "HLA-C"]
class_II = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1"]
ag_machinery = ["B2M", "TAP1", "TAP2"]
braf_cptc["HLA_I_z_mean"] = braf_cptc[[g for g in class_I if g in braf_cptc.columns]].mean(axis=1)
braf_cptc["HLA_II_z_mean"] = braf_cptc[[g for g in class_II if g in braf_cptc.columns]].mean(axis=1)
braf_cptc["AgPresentation_z_mean"] = braf_cptc[[g for g in ag_machinery if g in braf_cptc.columns]].mean(axis=1)

hla_features = (
    ["HLA_I_z_mean", "HLA_II_z_mean", "AgPresentation_z_mean"]
    + class_I + class_II + ag_machinery
    + ["JAK1", "JAK2", "STAT1", "IRF1"]
    + ["IFNG", "CTLA4", "PDCD1", "FOXP3", "CD274", "IDO1",
       "GZMA", "PRF1", "CD8A", "CD8B", "CXCL9", "CXCL10", "CXCL13"]
)
hla_features = [c for c in hla_features if c in braf_cptc.columns]
hla_rows = [compare(braf_cptc, c) for c in hla_features]
hla_df = pd.DataFrame(hla_rows)
hla_df.insert(0, "category", "HLA_AgPresentation")
hla_df.to_csv(OUT / "h19_hla_preservation.tsv", sep="\t", index=False)
print(hla_df.to_string(index=False))


# ---------------------------------------------------------------------------
# (3) Immune-escape signature
# ---------------------------------------------------------------------------
print("[H19] (3) Immune escape — counting low-HLA-I / low-B2M cases")
# define escape as HLA_I_z_mean < -0.5 OR B2M < -0.5 (rough Thorsson-style)
escape_thr = -0.5
braf_cptc["escape_HLA_I_low"] = (braf_cptc["HLA_I_z_mean"] < escape_thr).astype(int)
braf_cptc["escape_B2M_low"] = (braf_cptc["B2M"] < escape_thr).astype(int) if "B2M" in braf_cptc.columns else np.nan
braf_cptc["escape_score"] = braf_cptc[["escape_HLA_I_low", "escape_B2M_low"]].sum(axis=1, min_count=1)

# escape rates per group
esc_rows = []
for grp in ["DM1", "DM2"]:
    sub = braf_cptc[braf_cptc["dm"] == grp]
    n = len(sub)
    rate_hla = sub["escape_HLA_I_low"].mean() if n else np.nan
    rate_b2m = sub["escape_B2M_low"].mean() if n else np.nan
    rate_either = ((sub["escape_HLA_I_low"] == 1) | (sub["escape_B2M_low"] == 1)).mean() if n else np.nan
    esc_rows.append({"group": grp, "n": n,
                     "HLA_I_low_rate": rate_hla,
                     "B2M_low_rate": rate_b2m,
                     "either_low_rate": rate_either})
esc_df = pd.DataFrame(esc_rows)
# Fisher for either-low
from scipy.stats import fisher_exact
def fisher_p(grp_col: str, esc_col: str) -> float:
    tab = pd.crosstab(braf_cptc[grp_col], braf_cptc[esc_col])
    if tab.shape != (2, 2):
        return np.nan
    return float(fisher_exact(tab.values, alternative="two-sided")[1])

esc_df.attrs = {}
esc_p = {}
braf_cptc["either_low"] = ((braf_cptc["escape_HLA_I_low"] == 1) | (braf_cptc["escape_B2M_low"] == 1)).astype(int)
esc_p["fisher_p_HLA_I_low"] = fisher_p("dm", "escape_HLA_I_low")
esc_p["fisher_p_B2M_low"] = fisher_p("dm", "escape_B2M_low")
esc_p["fisher_p_either_low"] = fisher_p("dm", "either_low")
print("  escape rates:", esc_df.to_dict("records"))
print("  Fisher p:", esc_p)
esc_df.to_csv(OUT / "h19_immune_escape_score.tsv", sep="\t", index=False)
with open(OUT / "h19_immune_escape_fisher.json", "w") as f:
    json.dump(esc_p, f, indent=2)


# ---------------------------------------------------------------------------
# (4) Neoantigen load + presentation efficiency
# ---------------------------------------------------------------------------
print("[H19] (4) Neoantigen / presentation efficiency")
neo_features = [
    "n_nonsyn_snv", "n_immunogenic_mut", "n_binding_pMHC", "n_binding_expr_pMHC",
    "n_neoantigen_indel", "total_neoantigen_pMHC",
    "n_binding_pMHC_log1p", "n_binding_expr_pMHC_log1p", "total_neoantigen_log1p",
]
neo_features = [c for c in neo_features if c in braf_cptc.columns]
neo_rows = [compare(braf_cptc, c) for c in neo_features]

# compute presentation efficiency = n_immunogenic_mut / max(n_nonsyn_snv,1) per sample
if {"n_nonsyn_snv", "n_immunogenic_mut"}.issubset(braf_cptc.columns):
    braf_cptc["presentation_eff"] = braf_cptc["n_immunogenic_mut"] / braf_cptc["n_nonsyn_snv"].replace(0, np.nan)
    neo_rows.append(compare(braf_cptc, "presentation_eff"))

# composite immunogenicity proxy = nonsil_perMb × HLA_I_z_mean (presented neoantigen rate)
if {"thorsson_nonsil_perMb", "HLA_I_z_mean"}.issubset(braf_cptc.columns):
    braf_cptc["TMB_x_HLAI"] = braf_cptc["thorsson_nonsil_perMb"] * braf_cptc["HLA_I_z_mean"]
    neo_rows.append(compare(braf_cptc, "TMB_x_HLAI"))

neo_df = pd.DataFrame(neo_rows)
neo_df.insert(0, "category", "Neoantigen_load")
neo_combined = pd.concat([tmb_df, hla_df, neo_df], ignore_index=True)
neo_combined.to_csv(OUT / "h19_all_features_combined.tsv", sep="\t", index=False)
print(neo_df.to_string(index=False))


# ---------------------------------------------------------------------------
# (5) Korean K2 arcasHLA — DM1-vs-DM2 allele frequency comparison
# ---------------------------------------------------------------------------
print("[H19] (5) Korean K2 arcasHLA DM1 vs DM2 allele frequencies")
k2 = pd.read_csv(K2_HLA, sep="\t")
print(f"  K2 n={len(k2)}, DM1={(k2.DM_call=='DM1').sum()}, DM2={(k2.DM_call=='DM2').sum()}")

# 4-digit allele columns
allele_cols = [c for c in k2.columns if c.endswith("_4d")]
# build long table
long = []
for col in allele_cols:
    locus = col.replace("_a1_4d", "").replace("_a2_4d", "").replace("_4d", "")
    locus = locus.split("_")[0]
    sub = k2[[col, "DM_call"]].dropna()
    sub = sub.rename(columns={col: "allele"})
    sub["locus"] = locus
    long.append(sub)
long_df = pd.concat(long, ignore_index=True)
long_df = long_df[long_df["DM_call"].isin(["DM1", "DM2"])]

# carrier frequency per allele per group: each sample contributes once per locus
# we use simple allele-count rate (per allele copy)
n_dm1 = (k2["DM_call"] == "DM1").sum()
n_dm2 = (k2["DM_call"] == "DM2").sum()
total_alleles_per_sample = 2  # diploid
allele_rows = []
for (locus, allele), grp in long_df.groupby(["locus", "allele"]):
    n_dm1_carry = (grp["DM_call"] == "DM1").sum()
    n_dm2_carry = (grp["DM_call"] == "DM2").sum()
    # frequency = allele_copies / (n_samples * 2)
    freq_dm1 = n_dm1_carry / (n_dm1 * 2) if n_dm1 else np.nan
    freq_dm2 = n_dm2_carry / (n_dm2 * 2) if n_dm2 else np.nan
    if (n_dm1_carry + n_dm2_carry) < 5:
        continue  # skip rare alleles
    # 2x2 chi-square / Fisher: carry vs not carry
    a = n_dm1_carry; b = n_dm1 * 2 - n_dm1_carry
    c = n_dm2_carry; d = n_dm2 * 2 - n_dm2_carry
    try:
        p = float(fisher_exact([[a, b], [c, d]], alternative="two-sided")[1])
    except Exception:
        p = np.nan
    allele_rows.append({
        "locus": locus, "allele": allele,
        "n_DM1_copies": int(n_dm1_carry), "n_DM2_copies": int(n_dm2_carry),
        "freq_DM1": freq_dm1, "freq_DM2": freq_dm2,
        "log2_FC_DM1_over_DM2": np.log2((freq_dm1 + 1e-3) / (freq_dm2 + 1e-3)),
        "fisher_p": p,
    })
allele_df = pd.DataFrame(allele_rows).sort_values("fisher_p")
# BH-FDR
from statsmodels.stats.multitest import multipletests
if len(allele_df):
    allele_df["fdr_bh"] = multipletests(allele_df["fisher_p"].fillna(1.0).values, method="fdr_bh")[1]
allele_df.to_csv(OUT / "h19_hla_alleles_dm1.tsv", sep="\t", index=False)
print(allele_df.head(15).to_string(index=False))


# ---------------------------------------------------------------------------
# headline json
# ---------------------------------------------------------------------------
def get_d_p(df: pd.DataFrame, feat: str):
    row = df[df["feature"] == feat]
    if not len(row):
        return None, None
    return float(row.iloc[0]["cohen_d_DM1_vs_DM2"]), float(row.iloc[0]["wilcoxon_p"])

headline = {
    "stratum": "BRAF_like × cPTC (TCGA-THCA)",
    "n_DM1": int((braf_cptc["dm"] == "DM1").sum()),
    "n_DM2": int((braf_cptc["dm"] == "DM2").sum()),
    "tmb_nonsil_perMb_d_p": get_d_p(tmb_df, "thorsson_nonsil_perMb"),
    "n_nonsyn_snv_d_p": get_d_p(tmb_df, "n_nonsyn_snv"),
    "classI_nonsilent_d_p": get_d_p(tmb_df, "classI_nonsilent_count"),
    "HLA_I_z_d_p": get_d_p(hla_df, "HLA_I_z_mean"),
    "HLA_II_z_d_p": get_d_p(hla_df, "HLA_II_z_mean"),
    "AgPresentation_z_d_p": get_d_p(hla_df, "AgPresentation_z_mean"),
    "B2M_d_p": get_d_p(hla_df, "B2M"),
    "TAP1_d_p": get_d_p(hla_df, "TAP1"),
    "TAP2_d_p": get_d_p(hla_df, "TAP2"),
    "JAK1_d_p": get_d_p(hla_df, "JAK1"),
    "JAK2_d_p": get_d_p(hla_df, "JAK2"),
    "STAT1_d_p": get_d_p(hla_df, "STAT1"),
    "IFNG_d_p": get_d_p(hla_df, "IFNG"),
    "FOXP3_d_p": get_d_p(hla_df, "FOXP3"),
    "CTLA4_d_p": get_d_p(hla_df, "CTLA4"),
    "PDCD1_d_p": get_d_p(hla_df, "PDCD1"),
    "n_binding_pMHC_d_p": get_d_p(neo_df, "n_binding_pMHC"),
    "presentation_eff_d_p": get_d_p(neo_df, "presentation_eff") if "presentation_eff" in neo_df["feature"].values else None,
    "TMB_x_HLAI_d_p": get_d_p(neo_df, "TMB_x_HLAI") if "TMB_x_HLAI" in neo_df["feature"].values else None,
    "escape_rates": esc_df.to_dict("records"),
    "escape_fisher_p": esc_p,
    "K2_arcasHLA_n_DM1": int(n_dm1),
    "K2_arcasHLA_n_DM2": int(n_dm2),
    "K2_top_DM1_enriched_alleles_at_p_lt_0p05": (
        allele_df[(allele_df.fisher_p < 0.05) & (allele_df.log2_FC_DM1_over_DM2 > 0)].head(5).to_dict("records")
    ),
}
with open(OUT / "h19_headline.json", "w") as f:
    json.dump(headline, f, indent=2, default=str)

# sample-level merged for record (full BRAF-cPTC including not_DM)
braf_cptc_full.to_csv(OUT / "h19_braf_cptc_merged.tsv", sep="\t", index=False)

print("\n[H19] DONE — outputs in", OUT)
print(json.dumps(headline, indent=2, default=str)[:3000])
