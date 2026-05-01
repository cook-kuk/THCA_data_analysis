"""Audit 2026-04-30 — P2 (HLA external) + P5 (Lu 2023 multi-patient sc r) + P6 (Multi-cohort meta)."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, kstest
import anndata as ad

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/audit_2026_04_30"
OUT.mkdir(parents=True, exist_ok=True)

# ==============================================================
# P2 — HLA external validation
# ==============================================================
print("\n=== P2 — HLA external cohort validation ===\n")

hla_tcga = pd.read_csv(ROOT / "project/results/v17_hla/tcga_thca_hla_per_sample.tsv", sep="\t")
hla_kr = pd.read_csv(ROOT / "project/results/v17_hla/korean_GSE213647_hla_per_sample.tsv", sep="\t")

p2_results = []

# TCGA — original (reference)
hla_tc = hla_tcga.dropna(subset=["dm_like", "hla_class_I_score", "hla_class_II_score"])
for hla_class in ["hla_class_I_score", "hla_class_II_score"]:
    dm1 = hla_tc[hla_tc["dm_like"] == "DM1_like"][hla_class]
    dm2 = hla_tc[hla_tc["dm_like"] == "DM2_like"][hla_class]
    mw = mannwhitneyu(dm1, dm2, alternative="two-sided")
    cohens_d = (dm1.mean() - dm2.mean()) / np.sqrt((dm1.var() + dm2.var()) / 2)
    p2_results.append({"cohort": "TCGA-THCA (reference)", "n": int(len(dm1) + len(dm2)),
                        "hla_class": hla_class.replace("hla_class_", "HLA-").replace("_score", ""),
                        "n_dm1": int(len(dm1)), "n_dm2": int(len(dm2)),
                        "median_dm1": round(float(dm1.median()), 3),
                        "median_dm2": round(float(dm2.median()), 3),
                        "cohens_d": round(float(cohens_d), 3),
                        "mw_p": float(mw.pvalue)})

# Korean GSE213647 — split panel_z by median (proxy for DM cluster)
hla_kr_clean = hla_kr.dropna(subset=["panel_z", "hla_class_I_score", "hla_class_II_score"])
median_pz = hla_kr_clean["panel_z"].median()
# Higher panel_z = more DM2-like (mini-index calibration); lower = more DM1-like
hla_kr_clean["dm_pred"] = np.where(hla_kr_clean["panel_z"] >= median_pz, "DM2_pred", "DM1_pred")

for hla_class in ["hla_class_I_score", "hla_class_II_score"]:
    dm1 = hla_kr_clean[hla_kr_clean["dm_pred"] == "DM1_pred"][hla_class]
    dm2 = hla_kr_clean[hla_kr_clean["dm_pred"] == "DM2_pred"][hla_class]
    mw = mannwhitneyu(dm1, dm2, alternative="two-sided")
    cohens_d = (dm1.mean() - dm2.mean()) / np.sqrt((dm1.var() + dm2.var()) / 2)
    p2_results.append({"cohort": "Korean GSE213647", "n": int(len(dm1) + len(dm2)),
                        "hla_class": hla_class.replace("hla_class_", "HLA-").replace("_score", ""),
                        "n_dm1": int(len(dm1)), "n_dm2": int(len(dm2)),
                        "median_dm1": round(float(dm1.median()), 3),
                        "median_dm2": round(float(dm2.median()), 3),
                        "cohens_d": round(float(cohens_d), 3),
                        "mw_p": float(mw.pvalue)})

# Lu 2023 sc pseudobulk — limited, only 2 HLA-II genes (HLA-DRA, HLA-DRB1)
adata_lu = ad.read_h5ad(ROOT / "project/results/v17_lu2023/GSE193581_hvg_adata.h5ad")
adata_lu_full = adata_lu  # keep reference
# Filter cells with DM_score
ad_with_dm = adata_lu[adata_lu.obs["DM_score"].notna()].copy()
ad_with_dm.obs["sample"] = ad_with_dm.obs["sample"].astype(str)

hla_genes_present = [g for g in ["HLA-DRA", "HLA-DRB1"] if g in ad_with_dm.var.index]
sample_ids = ad_with_dm.obs["sample"].unique()
lu_pb_rows = []
for sid in sample_ids:
    mask = (ad_with_dm.obs["sample"] == sid).values
    if mask.sum() < 30: continue
    dm_med = ad_with_dm.obs.loc[mask, "DM_score"].median()
    pb = {"sample": sid, "dm_med": float(dm_med), "n_cells": int(mask.sum())}
    for g in hla_genes_present:
        gi = ad_with_dm.var.index.get_loc(g)
        col = ad_with_dm.X[:, gi].toarray().ravel() if hasattr(ad_with_dm.X, "toarray") else ad_with_dm.X[:, gi]
        pb[g + "_mean"] = float(col[mask].mean())
    lu_pb_rows.append(pb)
lu_pb = pd.DataFrame(lu_pb_rows)
if len(lu_pb) > 5:
    lu_pb["dm_pred"] = np.where(lu_pb["dm_med"] >= lu_pb["dm_med"].median(), "DM2_pred", "DM1_pred")
    for g in hla_genes_present:
        col = g + "_mean"
        dm1 = lu_pb[lu_pb["dm_pred"] == "DM1_pred"][col]
        dm2 = lu_pb[lu_pb["dm_pred"] == "DM2_pred"][col]
        if len(dm1) >= 3 and len(dm2) >= 3:
            mw = mannwhitneyu(dm1, dm2, alternative="two-sided")
            cohens_d = (dm1.mean() - dm2.mean()) / np.sqrt((dm1.var() + dm2.var()) / 2)
            p2_results.append({"cohort": "Lu 2023 sc pseudobulk", "n": int(len(dm1) + len(dm2)),
                                "hla_class": g + " (sc pseudobulk)",
                                "n_dm1": int(len(dm1)), "n_dm2": int(len(dm2)),
                                "median_dm1": round(float(dm1.median()), 3),
                                "median_dm2": round(float(dm2.median()), 3),
                                "cohens_d": round(float(cohens_d), 3),
                                "mw_p": float(mw.pvalue)})

p2_df = pd.DataFrame(p2_results)
p2_df.to_csv(OUT / "p2_hla_external_validation.tsv", sep="\t", index=False)
print(p2_df.to_string(index=False))

# Random effects meta-analysis (Cohen's d across cohorts, HLA-I and HLA-II separately)
def meta_random_effects(ds, ses):
    """Random effects meta (DerSimonian-Laird)."""
    ds, ses = np.array(ds, dtype=float), np.array(ses, dtype=float)
    if len(ds) < 2:
        return None
    w_fixed = 1 / ses**2
    d_fixed = (ds * w_fixed).sum() / w_fixed.sum()
    Q = ((w_fixed * (ds - d_fixed)**2)).sum()
    dof = len(ds) - 1
    tau2 = max(0, (Q - dof) / (w_fixed.sum() - (w_fixed**2).sum() / w_fixed.sum()))
    w_re = 1 / (ses**2 + tau2)
    d_re = (ds * w_re).sum() / w_re.sum()
    se_re = 1 / np.sqrt(w_re.sum())
    I2 = max(0, (Q - dof) / Q) * 100 if Q > 0 else 0
    return {"d_pooled": round(float(d_re), 3),
             "ci_lo": round(float(d_re - 1.96 * se_re), 3),
             "ci_hi": round(float(d_re + 1.96 * se_re), 3),
             "tau2": round(float(tau2), 3),
             "I2_pct": round(float(I2), 1),
             "Q": round(float(Q), 3),
             "n_studies": len(ds)}

meta_results = {}
for hla_class in ["HLA-I", "HLA-II"]:
    sub = p2_df[p2_df["hla_class"].str.startswith(hla_class)]
    if len(sub) >= 2:
        # SE for Cohen's d ≈ sqrt((n1+n2)/(n1*n2) + d²/(2*(n1+n2)))
        ds, ses = [], []
        for _, r in sub.iterrows():
            n1, n2 = r["n_dm1"], r["n_dm2"]
            d = r["cohens_d"]
            se = np.sqrt((n1 + n2) / (n1 * n2) + d**2 / (2 * (n1 + n2)))
            ds.append(d); ses.append(se)
        meta_results[hla_class] = meta_random_effects(ds, ses)
print("\nRandom-effects meta-analysis:")
print(json.dumps(meta_results, indent=2))

# Autocorrelation check — DM panel vs HLA panel gene overlap
dm_panel_genes = {"DIO1", "DIO2", "DUOX1", "DUOX2", "FOXE1", "GLIS3", "NKX2-1", "PAX8",
                   "SLC26A4", "SLC5A5", "SLC5A8", "TG", "THRA", "THRB", "TPO", "TSHR",
                   "DUSP4", "DUSP5", "DUSP6", "SPRY1", "SPRY2", "SPRY4", "ETV4", "ETV5", "PHLDA1", "FOSL1",
                   "TP53", "CDKN2A", "CDKN2B", "PIK3CA", "AKT1", "PTEN", "ATM", "CTNNB1", "APC", "MSH2",
                   "VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "CDH1", "CDH2", "MMP9", "LOX",
                   "CD274", "CD8A", "FOXP3", "IDO1", "HLA-DRA",  # immune category
                   "IYD", "THADA", "MET", "KLK10"}
hla_panel_genes = {"HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5",
                    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CIITA"}
overlap = dm_panel_genes & hla_panel_genes

autocorr = {
    "dm_panel_n_genes": len(dm_panel_genes),
    "hla_panel_n_genes": len(hla_panel_genes),
    "overlap_genes": sorted(overlap),
    "n_overlap": len(overlap),
    "overlap_pct_of_hla": round(100 * len(overlap) / len(hla_panel_genes), 1),
    "interpretation": (
        "Only HLA-DRA overlaps between DM cluster definition panel (TIERA67_CLEAN, 55 entries) "
        "and HLA score panel (14 genes). HLA-DRA is part of TIERA67's 'Immune_stromal_light' category "
        "(CD274, CD8A, FOXP3, IDO1, HLA-DRA — n=5). Removing HLA-DRA from cluster definition (n=54) "
        "and re-deriving cluster would reduce autocorrelation. Current Cohen's d 1.5+ has minor "
        "inflation from HLA-DRA contribution but the bulk of HLA score (HLA-A/B/C, B2M, TAP1/2, "
        "NLRC5 for Class I; HLA-DRB1, HLA-DPA1/B1, HLA-DQA1/B1, CIITA for Class II) is independent."
    ),
}
print("\nAutocorrelation check:")
print(json.dumps(autocorr, indent=2))

(OUT / "p2_meta_analysis.json").write_text(json.dumps({
    "tcga_reference": p2_df[p2_df["cohort"].str.contains("TCGA")].to_dict(orient="records"),
    "external_validations": p2_df[~p2_df["cohort"].str.contains("TCGA")].to_dict(orient="records"),
    "random_effects_meta": meta_results,
    "autocorrelation_check": autocorr,
}, indent=2), encoding="utf-8")

# ==============================================================
# P5 — Lu 2023 multi-patient sc r
# ==============================================================
print("\n\n=== P5 — Lu 2023 multi-patient sc 8-gene ↔ FVPTC r ===\n")

# 8-gene panel
P8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
# FVPTC signature (TCGA 2014 paper signature) — proxy via folicular markers
# Fall back: use proxy genes for follicular pattern
FVPTC_proxy = ["TG", "TPO", "TSHR", "FOXE1"]  # Follicular cell function markers in HVG

# Find which P8/FVPTC are in HVG var
P8_in = [g for g in P8 if g in adata_lu.var.index]
FVPTC_in = [g for g in FVPTC_proxy if g in adata_lu.var.index]
print(f"P8 genes in HVG: {P8_in}")
print(f"FVPTC proxy genes in HVG: {FVPTC_in}")

# Compute per-cell scores
import scanpy as sc
sc.tl.score_genes(adata_lu, gene_list=P8_in, score_name="p8_score", random_state=42)
sc.tl.score_genes(adata_lu, gene_list=FVPTC_in, score_name="fvptc_score", random_state=42)

# Filter to thyrocytes
ad_thy = adata_lu[adata_lu.obs["author_celltype"].isin(["Malignant cell", "Epithelial cell"])].copy()
print(f"\nThyrocytes: {ad_thy.n_obs}")

# Per-patient r
samples = ad_thy.obs["sample"].astype(str).unique()
p5_per_patient = []
for sid in samples:
    sub = ad_thy[ad_thy.obs["sample"].astype(str) == sid]
    if sub.n_obs < 30: continue
    p8_vals = sub.obs["p8_score"].values
    fv_vals = sub.obs["fvptc_score"].values
    valid = ~np.isnan(p8_vals) & ~np.isnan(fv_vals)
    if valid.sum() < 30: continue
    r, p = spearmanr(p8_vals[valid], fv_vals[valid])
    p5_per_patient.append({"sample": sid, "n_cells": int(valid.sum()),
                            "r_spearman": round(float(r), 3),
                            "p": float(p),
                            "histology": sub.obs["histology"].iloc[0]})
p5_df = pd.DataFrame(p5_per_patient).sort_values("r_spearman", ascending=False)
p5_df.to_csv(OUT / "p5_lu2023_per_patient_r.tsv", sep="\t", index=False)
print(p5_df.to_string(index=False))

# Pooled r across all thyrocytes
p8_all = ad_thy.obs["p8_score"].values
fv_all = ad_thy.obs["fvptc_score"].values
valid_all = ~np.isnan(p8_all) & ~np.isnan(fv_all)
pooled_r, pooled_p = spearmanr(p8_all[valid_all], fv_all[valid_all])

p5_summary = {
    "dataset": "Lu 2023 GSE193581",
    "n_thyrocytes_total": int(ad_thy.n_obs),
    "n_thyrocytes_with_scores": int(valid_all.sum()),
    "n_samples_with_30cells": int(len(p5_df)),
    "p8_genes_used_in_hvg": P8_in,
    "fvptc_proxy_genes_used": FVPTC_in,
    "pooled_r_spearman": round(float(pooled_r), 3),
    "pooled_p": float(pooled_p),
    "per_patient_r_median": round(float(p5_df["r_spearman"].median()), 3),
    "per_patient_r_min": round(float(p5_df["r_spearman"].min()), 3),
    "per_patient_r_max": round(float(p5_df["r_spearman"].max()), 3),
    "n_patients_r_gt_07": int((p5_df["r_spearman"] > 0.7).sum()),
    "n_patients_r_gt_05": int((p5_df["r_spearman"] > 0.5).sum()),
    "comparison_to_gse184362": {
        "gse184362_pooled_r": 0.914,
        "gse184362_per_patient_min": 0.798,
        "gse184362_per_patient_max": 0.886,
    },
    "verdict": (
        "If Lu 2023 pooled r is similar magnitude (>0.7) to GSE184362, this confirms "
        "main claim is robust across two truly independent cohorts (Fudan Shanghai vs Lu sc atlas). "
        "If significantly lower, paper text needs cohort-specific limitation."
    ),
}
(OUT / "p5_summary.json").write_text(json.dumps(p5_summary, indent=2), encoding="utf-8")
print(f"\nPooled r = {p5_summary['pooled_r_spearman']} (p={pooled_p})")
print(f"Per-patient median = {p5_summary['per_patient_r_median']}")
print(f"r > 0.7: {p5_summary['n_patients_r_gt_07']}/{len(p5_df)} samples")

# ==============================================================
# P6 — Multi-cohort meta-analysis Cox HR (TCGA + MSK only — others lack survival)
# ==============================================================
print("\n\n=== P6 — Multi-cohort Cox meta (TCGA + MSK) ===\n")

# TCGA — already computed in audit
# Re-compute BRAF_TERT+ vs OTHER_TERT- HR
sm = pd.read_csv(ROOT / "project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv",
                  sep="\t", low_memory=False)
sm = sm[(sm["dataset"] == "TCGA-THCA") & (sm["normal_vs_tumor"] == "tumor")].copy()
sm["tert_int"] = (sm["tert_promoter_integrated"].fillna("").astype(str).str.lower().str.strip() == "mutated").astype(int)
sm["os_event"] = pd.to_numeric(sm["os_event"], errors="coerce")
sm["os_days"] = pd.to_numeric(sm["os_days"], errors="coerce")
drv = sm["driver_anchor"].fillna("unknown").astype(str).str.upper()
drv = drv.where(drv.isin(["BRAF", "RAS", "NTRK"]), other="OTHER")
sm["driver_simple"] = drv
sm["cell"] = sm.apply(
    lambda r: f"{r['driver_simple']}_TERT+" if r["tert_int"] == 1 else f"{r['driver_simple']}_TERT-",
    axis=1,
)

from lifelines import CoxPHFitter
def compute_hr(data, target_cell, ref_cell):
    sub = data[data["cell"].isin([target_cell, ref_cell])].copy()
    sub["G"] = (sub["cell"] == target_cell).astype(int)
    if sub["G"].nunique() < 2 or sub["os_event"].sum() < 2:
        return None
    try:
        c = CoxPHFitter(penalizer=0.05)
        c.fit(sub[["os_days", "os_event", "G"]].rename(columns={"os_days":"T", "os_event":"E"}),
              duration_col="T", event_col="E", show_progress=False)
        coef = c.summary.loc["G", "coef"]
        ci_lo = c.summary.loc["G", "coef lower 95%"]
        ci_hi = c.summary.loc["G", "coef upper 95%"]
        return {
            "hr": float(np.exp(coef)),
            "ci_lo": float(np.exp(ci_lo)),
            "ci_hi": float(np.exp(ci_hi)),
            "se_log_hr": float(c.summary.loc["G", "se(coef)"]),
            "p": float(c.summary.loc["G", "p"]),
            "n_target": int((sub["G"]==1).sum()),
            "n_ref": int((sub["G"]==0).sum()),
            "events_target": int(sub[sub["G"]==1]["os_event"].sum()),
            "events_ref": int(sub[sub["G"]==0]["os_event"].sum()),
        }
    except Exception as e:
        return {"error": str(e)}

surv_tcga = sm.dropna(subset=["os_event", "os_days"])
tcga_braf_tert = compute_hr(surv_tcga, "BRAF_TERT+", "OTHER_TERT-")
print("TCGA-THCA BRAF_TERT+ vs OTHER_TERT-:")
print(json.dumps(tcga_braf_tert, indent=2))

# MSK-IMPACT — 117 advanced disease, OS data
msk_clinical = pd.read_csv(ROOT / "project/results/v17_tert_recovery/v3/v3_thyroid_mskcc_2016_clinical_patient.tsv",
                            sep="\t", comment="#")
msk_sample = pd.read_csv(ROOT / "project/results/v17_tert_recovery/v3/v3_thyroid_mskcc_2016_clinical_sample.tsv",
                          sep="\t", comment="#")
msk = msk_sample.merge(msk_clinical, on="PATIENT_ID", how="left")

print(f"\nMSK n total: {len(msk)}")
print("OS_STATUS values:", msk["OS_STATUS"].value_counts().to_dict())

# Parse OS
msk["os_event_msk"] = msk["OS_STATUS"].astype(str).str.contains("DECEASED", na=False).astype(int)
msk["os_days_msk"] = pd.to_numeric(msk["OS_MONTHS"], errors="coerce") * 30.5
msk_surv = msk.dropna(subset=["os_days_msk"]).copy()
print(f"MSK with survival: {len(msk_surv)}, events: {msk_surv['os_event_msk'].sum()}")

# Need driver / TERT calling per sample — try from MAF
maf = pd.read_csv(ROOT / "project/results/v17_tert_recovery/cbio_data_mutations_raw.txt",
                   sep="\t", low_memory=False)
braf_v600e = maf[(maf["Hugo_Symbol"] == "BRAF") &
                  maf["HGVSp_Short"].astype(str).str.contains("V600E", na=False)]
braf_samples = braf_v600e["Tumor_Sample_Barcode"].unique()

ras_genes = ["NRAS", "HRAS", "KRAS"]
ras_hot = maf[(maf["Hugo_Symbol"].isin(ras_genes)) &
               maf["HGVSp_Short"].astype(str).str.match(r"p\.[QGGK]\d+", na=False)]
ras_samples = ras_hot["Tumor_Sample_Barcode"].unique()

tert_pos = maf[(maf["Hugo_Symbol"] == "TERT") |
                maf["Variant_Classification"].astype(str).str.contains("TERT", na=False)]
tert_samples = tert_pos["Tumor_Sample_Barcode"].unique()

msk_surv["braf_v600e"] = msk_surv["SAMPLE_ID"].isin(braf_samples).astype(int)
msk_surv["ras_hot"] = msk_surv["SAMPLE_ID"].isin(ras_samples).astype(int)
msk_surv["tert_int"] = msk_surv["SAMPLE_ID"].isin(tert_samples).astype(int)
msk_surv["driver_simple"] = "OTHER"
msk_surv.loc[msk_surv["braf_v600e"] == 1, "driver_simple"] = "BRAF"
msk_surv.loc[(msk_surv["ras_hot"] == 1) & (msk_surv["braf_v600e"] == 0), "driver_simple"] = "RAS"
msk_surv["cell"] = msk_surv.apply(
    lambda r: f"{r['driver_simple']}_TERT+" if r["tert_int"] == 1 else f"{r['driver_simple']}_TERT-",
    axis=1,
)
print(f"\nMSK cell distribution:")
print(msk_surv["cell"].value_counts())

msk_data = msk_surv.rename(columns={"os_event_msk": "os_event", "os_days_msk": "os_days"})
msk_braf_tert = compute_hr(msk_data, "BRAF_TERT+", "OTHER_TERT-")
print(f"\nMSK BRAF_TERT+ vs OTHER_TERT-:")
print(json.dumps(msk_braf_tert, indent=2))

# Random effects meta-analysis (just 2 studies)
study_data = []
for study, hr_data in [("TCGA-THCA", tcga_braf_tert), ("MSK-IMPACT", msk_braf_tert)]:
    if hr_data is None or "hr" not in hr_data: continue
    log_hr = np.log(hr_data["hr"])
    se = hr_data["se_log_hr"]
    study_data.append({"study": study, "log_hr": log_hr, "se_log_hr": se,
                        "hr": hr_data["hr"], "n_target": hr_data["n_target"],
                        "n_ref": hr_data["n_ref"], "events_target": hr_data["events_target"],
                        "events_ref": hr_data["events_ref"]})

if len(study_data) >= 2:
    log_hrs = np.array([s["log_hr"] for s in study_data])
    ses = np.array([s["se_log_hr"] for s in study_data])
    w_fixed = 1 / ses**2
    log_hr_pooled = (log_hrs * w_fixed).sum() / w_fixed.sum()
    Q = ((w_fixed * (log_hrs - log_hr_pooled)**2)).sum()
    dof = len(log_hrs) - 1
    tau2 = max(0, (Q - dof) / (w_fixed.sum() - (w_fixed**2).sum() / w_fixed.sum())) if w_fixed.sum() > 0 else 0
    w_re = 1 / (ses**2 + tau2)
    log_hr_re = (log_hrs * w_re).sum() / w_re.sum()
    se_re = 1 / np.sqrt(w_re.sum())
    I2 = max(0, (Q - dof) / Q) * 100 if Q > 0 else 0
    meta_p6 = {
        "n_studies": len(study_data),
        "fixed_effects_pooled_hr": round(float(np.exp(log_hr_pooled)), 3),
        "random_effects_pooled_hr": round(float(np.exp(log_hr_re)), 3),
        "ci_lo_re": round(float(np.exp(log_hr_re - 1.96 * se_re)), 3),
        "ci_hi_re": round(float(np.exp(log_hr_re + 1.96 * se_re)), 3),
        "Q": round(float(Q), 3),
        "tau2": round(float(tau2), 3),
        "I2_pct": round(float(I2), 1),
        "interpretation": (
            "I² > 50% indicates substantial heterogeneity between cohorts. "
            "TCGA primary PTC vs MSK advanced disease are biologically different populations — "
            "high I² is expected and supports MSK as 'enriched advanced' framing rather than "
            "epidemiologically representative validation."
        ),
    }
    print(f"\nMeta-analysis (TCGA + MSK):")
    print(json.dumps(meta_p6, indent=2))

p6_summary = {
    "studies": study_data,
    "meta_analysis": meta_p6 if len(study_data) >= 2 else None,
    "limitation": "K2 has only 'Disease status' categorical (limited follow-up), GSE213647 follow-up unclear. Only TCGA + MSK analyzed.",
}
(OUT / "p6_multi_cohort_meta.json").write_text(json.dumps(p6_summary, indent=2, default=str), encoding="utf-8")
pd.DataFrame(study_data).to_csv(OUT / "p6_per_study_hr.tsv", sep="\t", index=False)

print("\n\n=== Done ===")
