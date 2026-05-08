#!/usr/bin/env python3
"""
Rapid high-impact topic pilot screen.

Scope:
- Includes the current CNV residual-class idea.
- Adds 10 additional candidate paper directions.
- Uses only local public-data-derived outputs already present in the project.
- Writes a Korean decision report plus machine-readable TSVs.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/high_impact_topic_pilots_2026_05_08"
REPORT = ROOT / "project/reports/2026_05_08_high_impact_11_topic_pilots.md"
OUT.mkdir(parents=True, exist_ok=True)


def cohen_d(a, b) -> float:
    a = pd.Series(a).dropna().astype(float)
    b = pd.Series(b).dropna().astype(float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    pooled = ((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2)
    if pooled <= 0 or np.isnan(pooled):
        return np.nan
    return float((a.mean() - b.mean()) / math.sqrt(pooled))


def mw(a, b) -> float:
    a = pd.Series(a).dropna().astype(float)
    b = pd.Series(b).dropna().astype(float)
    if len(a) < 1 or len(b) < 1:
        return np.nan
    return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)


def fisher_flag(df: pd.DataFrame, group_col: str, group_a, group_b, flag_col: str) -> dict:
    sub = df[df[group_col].isin([group_a, group_b])].copy()
    sub = sub.dropna(subset=[flag_col])
    a = sub[sub[group_col] == group_a][flag_col].astype(bool)
    b = sub[sub[group_col] == group_b][flag_col].astype(bool)
    table = [[int(a.sum()), int((~a).sum())], [int(b.sum()), int((~b).sum())]]
    if min(len(a), len(b)) == 0:
        odds, p = np.nan, np.nan
    else:
        odds, p = stats.fisher_exact(table)
    return {
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "rate_a": float(a.mean()) if len(a) else np.nan,
        "rate_b": float(b.mean()) if len(b) else np.nan,
        "odds_ratio": float(odds) if np.isfinite(odds) else np.inf,
        "p": float(p) if not np.isnan(p) else np.nan,
    }


def cont_metric(df: pd.DataFrame, group_col: str, group_a, group_b, value_col: str) -> dict:
    sub = df[df[group_col].isin([group_a, group_b])].copy()
    a = sub[sub[group_col] == group_a][value_col]
    b = sub[sub[group_col] == group_b][value_col]
    return {
        "n_a": int(a.dropna().shape[0]),
        "n_b": int(b.dropna().shape[0]),
        "mean_a": float(a.mean()) if len(a.dropna()) else np.nan,
        "mean_b": float(b.mean()) if len(b.dropna()) else np.nan,
        "diff_a_minus_b": float(a.mean() - b.mean()) if len(a.dropna()) and len(b.dropna()) else np.nan,
        "cohen_d_a_minus_b": cohen_d(a, b),
        "p": mw(a, b),
    }


def add_metric(rows: list[dict], topic_id: str, topic: str, metric: str, result: dict, note: str = "") -> None:
    flat = {
        "topic_id": topic_id,
        "topic": topic,
        "metric": metric,
        "note": note,
    }
    flat.update(result)
    rows.append(flat)


def p_text(p: float) -> str:
    if p is None or np.isnan(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.4f}"


def grade(effect_strength: str, validation: str, risk: str) -> tuple[str, int]:
    table = {
        ("strong", "ready", "low"): ("A", 95),
        ("strong", "ready", "medium"): ("A-", 90),
        ("strong", "partial", "medium"): ("B+", 84),
        ("strong", "partial", "high"): ("B", 78),
        ("moderate", "ready", "medium"): ("B+", 82),
        ("moderate", "partial", "medium"): ("B", 75),
        ("moderate", "partial", "high"): ("B-", 68),
        ("weak", "partial", "high"): ("C", 55),
        ("weak", "none", "high"): ("NO-GO", 35),
        ("negative", "ready", "medium"): ("negative-use-only", 60),
    }
    return table.get((effect_strength, validation, risk), ("BORDER", 62))


def main() -> None:
    metrics: list[dict] = []
    topics: list[dict] = []

    driver = pd.read_csv(ROOT / "project/results/aggressive_sprint_2026_05_06/tcga_driver_class_mechanism_matrix.tsv", sep="\t")
    class6 = pd.read_csv(ROOT / "project/results/aggressive_sprint_2026_05_06/class6_cnv_signature_samples.tsv", sep="\t")
    class6_arm = pd.read_csv(ROOT / "project/results/aggressive_sprint_2026_05_06/class6_dm1_dm2_arm_cnv.tsv", sep="\t")
    methyl_gene = pd.read_csv(ROOT / "project/results/audit_2026_04_30/round5/r5_2_per_gene_methylation_DM.tsv", sep="\t")
    ht = pd.read_csv(ROOT / "project/results/d4p2_tcga_hashimoto_signature/tcga_signature_scores.tsv", sep="\t")
    cox = pd.read_csv(ROOT / "project/results/03_pathology_poc/spark_tcga_cox_modules.tsv", sep="\t")
    cox_adj = pd.read_csv(ROOT / "project/results/03_pathology_poc/spark_tcga_multivariate_cox_DSS.tsv", sep="\t")
    cci_driver = pd.read_csv(ROOT / "project/results/03_pathology_poc/spark_tcga_driver_class_CCI.tsv", sep="\t")
    nichenet = pd.read_csv(ROOT / "project/results/03_pathology_poc/spark_nichenet_lite_pooled.tsv", sep="\t")
    lr_stage = pd.read_csv(ROOT / "project/results/03_pathology_poc/spark_cci_LR_activity_per_stage.tsv", sep="\t")
    spatial64 = pd.read_csv(ROOT / "project/results/spatial_full_2026_05_06/spatial_F64_28sample_multiaxis_morans.tsv", sep="\t")
    spatial47 = pd.read_csv(ROOT / "project/results/spatial_full_2026_05_06/spatial_F47_per_axis_moran_summary_v7plus.tsv", sep="\t")
    spatial_ht = pd.read_csv(ROOT / "project/results/spatial_full_2026_05_06/spatial_C_HT_TLS_spatial.tsv", sep="\t")
    spatial_v12 = pd.read_csv(ROOT / "project/results/spatial_full_2026_05_06/spatial_v12_stage_gradient.tsv", sep="\t")
    geomx = pd.read_csv(ROOT / "project/results/spatial_full_2026_05_06/GSE301163_per_ROI_scores.tsv", sep="\t")
    geomx_stats = pd.read_csv(ROOT / "project/results/spatial_full_2026_05_06/GSE301163_crossplatform_stats.tsv", sep="\t")
    closure = pd.read_csv(ROOT / "project/results/03_pathology_poc/closure_battery_metrics.tsv", sep="\t")

    # T00 — current topic: CNV residual class.
    topic = "T00 CNV-defined true-driver-negative residual class"
    add_metric(metrics, "T00", topic, "Class6-DM2 vs Class6-DM1 CNV-signature rate", fisher_flag(class6, "v17_dark_cluster", "DM2", "DM1", "cnv_signature_any"), "7-arm quick signature")
    for arm in class6_arm.head(7).itertuples(index=False):
        add_metric(metrics, "T00", topic, f"{arm.arm}_{arm.call} DM2 vs DM1", {"rate_a": arm.DM2_pct / 100, "rate_b": arm.DM1_pct / 100, "p": arm.p, "bh_q": arm.bh_q}, "top Class6 arm CNV")
    g, score = grade("strong", "partial", "medium")
    topics.append({
        "topic_id": "T00", "topic": topic, "rank_hint": 1, "grade": g, "score": score,
        "headline": "True driver-negative thyroid cancers contain a recurrent arm-CNV residual class concentrated in DM2/FVPTC-like tumors.",
        "pilot_result": "CNV signature marks Class6-DM2 far more than Class6-DM1; top arms include 7p/7q gain, 2p/2q loss, 16p/16q gain, 12q gain.",
        "next_gate": "Validate from independent GDC/cBioPortal GISTIC source; adjust for FVPTC/Hurthle histology and purity.",
    })

    # T01 — GeoMx miRNA-biogenesis progression.
    topic = "T01 miRNA-biogenesis-defective progression route (DICER1/DGCR8 GeoMx)"
    panck = geomx[geomx["label"] == "PanCK+"].copy()
    for val in ["DM1_axis", "RAI_8", "Eight_gene_DM1", "HLA_II", "Cell_cycle", "Hypoxia", "EMT"]:
        add_metric(metrics, "T01", topic, f"PDTC PanCK+ vs microPTC PanCK+ {val}", cont_metric(panck, "histology", "PDTC", "microPTC", val), "GSE301163 GeoMx")
    for row in geomx_stats.itertuples(index=False):
        if "PDTC PanCK+ vs microPTC" in row.comparison or "Eight_gene_DM1" in row.comparison:
            add_metric(metrics, "T01", topic, row.comparison, {"n_a": row.n_panck, "n_b": row.n_vim, "mean_a": row.panck_mean, "mean_b": row.vim_mean, "p": row.MW_p}, "precomputed GeoMx stat")
    g, score = grade("strong", "partial", "medium")
    topics.append({
        "topic_id": "T01", "topic": topic, "rank_hint": 2, "grade": g, "score": score,
        "headline": "DICER1/DGCR8-associated lesions trace a compartment-resolved route from microPTC to PDTC with lineage loss and immune/stromal divergence.",
        "pilot_result": "GeoMx n=78 ROI supports a PDTC-vs-microPTC lineage drop and HLA/DM-axis shifts, but genotype/histology are partially confounded.",
        "next_gate": "Patient-stratified model; separate DICER1 from DGCR8 and histology; add CNV/DM overlay if possible.",
    })

    # T02 — TROP2 spatial niche.
    topic = "T02 tumor-specific TROP2 spatial niche"
    spatial64["cancer_ptc_lptc"] = spatial64["condition"].isin(["PTC", "LPTC"])
    spatial64["negative_control"] = spatial64["condition"].isin(["PT", "CONTROL", "HT", "GD"])
    sub = spatial64[spatial64["cancer_ptc_lptc"] | spatial64["negative_control"]].copy()
    sub["group"] = np.where(sub["cancer_ptc_lptc"], "PTC_LPTC", "NEG_CTRL")
    add_metric(metrics, "T02", topic, "TROP2 Moran PTC/LPTC vs negative controls", cont_metric(sub, "group", "PTC_LPTC", "NEG_CTRL", "TROP2"), "28-slide Visium corpus")
    add_metric(metrics, "T02", topic, "TROP2 niche threshold rate >0.25", fisher_flag(sub.assign(trop2_high=sub["TROP2"] > 0.25), "group", "PTC_LPTC", "NEG_CTRL", "trop2_high"), "Moran's I threshold")
    g, score = grade("strong", "ready", "medium")
    topics.append({
        "topic_id": "T02", "topic": topic, "rank_hint": 3, "grade": g, "score": score,
        "headline": "Bulk-negative TROP2 is rescued as a tumor-specific spatial niche rather than a whole-tumor expression marker.",
        "pilot_result": "PTC/LPTC Visium slides show high TROP2 spatial autocorrelation versus normal/autoimmune controls.",
        "next_gate": "Keep as spatial-niche biology, not ADC vulnerability; add GeoMx/IF/IHC if pursuing high-tier.",
    })

    # T03 — HT TLS spatial niche.
    topic = "T03 Hashimoto-overlap PTC TLS / B-cell / IGHV spatial niche"
    for gene_set in ["HLA_II", "B_cell", "TLS", "IGHV_AICDA"]:
        vals = spatial_ht[spatial_ht["gene_set"] == gene_set]["morans_I"].astype(float)
        add_metric(metrics, "T03", topic, f"GSE230424 {gene_set} Moran summary", {"n": int(vals.shape[0]), "mean": float(vals.mean()), "min": float(vals.min()), "max": float(vals.max()), "rate_gt_0_24": float((vals > 0.24).mean())}, "4 PTC+HT slides")
    ptht = spatial47[(spatial47["condition"] == "PTC_HT") & (spatial47["axis"].isin(["HLA_II", "B_cell", "TLS", "IGHV_AICDA"]))]
    other = spatial47[(spatial47["condition"].isin(["PTC", "LPTC", "PT", "ATC", "CONTROL", "HT", "GD"])) & (spatial47["axis"].isin(["HLA_II", "B_cell", "TLS", "IGHV_AICDA"]))]
    add_metric(metrics, "T03", topic, "PTC_HT immune-Moran mean vs others", {"n_a": int(ptht["count"].sum()), "n_b": int(other["count"].sum()), "mean_a": float(ptht["mean"].mean()), "mean_b": float(other["mean"].mean()), "diff_a_minus_b": float(ptht["mean"].mean() - other["mean"].mean())}, "axis-level summary")
    g, score = grade("strong", "partial", "high")
    topics.append({
        "topic_id": "T03", "topic": topic, "rank_hint": 4, "grade": g, "score": score,
        "headline": "Hashimoto-overlap PTC contains spatially organized TLS/B-cell/IGHV niches visible directly in tissue.",
        "pilot_result": "4/4 GSE230424 slides show organized HLA-II/B-cell/TLS/IGHV Moran patterns; n=4 is the main limitation.",
        "next_gate": "Independent PTC+HT cohort or local IF/IHC; maintain no causal HT-to-cancer language.",
    })

    # T04 — CAF/ECM ligand axis.
    topic = "T04 CAF/ECM ligand remodeling as a suppressor context for thyroid-lineage/RAI programs"
    top_lig = nichenet.sort_values("med_r", ascending=False).head(8)
    for row in top_lig.itertuples(index=False):
        add_metric(metrics, "T04", topic, f"spatial ligand {row.ligand}", {"n": row.n_slides, "median_r": row.med_r, "q25": row.q25, "q75": row.q75}, "NicheNet-lite spatial ligand screen")
    neg = cci_driver[cci_driver["bucket"] == "NEG"]
    for module in ["CAF", "M1_M2", "TLS", "CD8"]:
        x = neg[neg["module"] == module]
        if not x.empty:
            add_metric(metrics, "T04", topic, f"TCGA NEG {module} vs RAI", {"n": int(x.iloc[0]["n"]), "spearman_r": float(x.iloc[0]["r_vs_RAI"])}, "BRAF/RAS-negative bulk replication")
    g, score = grade("moderate", "ready", "medium")
    topics.append({
        "topic_id": "T04", "topic": topic, "rank_hint": 5, "grade": g, "score": score,
        "headline": "A conserved CAF/ECM ligand program tracks spatial thyroid-lineage loss and bulk RAI suppression in driver-negative tumors.",
        "pilot_result": "DCN/COL1A2/SFRP2/COL1A1 are top spatial ligands; TCGA driver-negative CAF/M1_M2 modules correlate negatively with RAI.",
        "next_gate": "Disentangle stromal abundance from tumor-cell intrinsic loss; validate with GeoMx VIM vs PanCK compartments.",
    })

    # T05 — immune checkpoint spatial routing.
    topic = "T05 spatial immune-checkpoint ligand-receptor routes in aggressive thyroid cancer"
    top_lr = lr_stage.sort_values(["DM1_top_minus_bot_mean", "mean"], ascending=False).head(12)
    for row in top_lr.itertuples(index=False):
        add_metric(metrics, "T05", topic, f"{row.stage} {row.lr_pair}", {"n": row.n_slides, "mean_activity": row.mean, "dm1_top_minus_bot": row.DM1_top_minus_bot_mean}, "spatial ligand-receptor activity")
    g, score = grade("moderate", "partial", "high")
    topics.append({
        "topic_id": "T05", "topic": topic, "rank_hint": 9, "grade": g, "score": score,
        "headline": "ATC/LPTC progression may route through spatially localized LGALS9-HAVCR2, CXCL10-CXCR3, and CCL5-CCR5 immune-checkpoint circuits.",
        "pilot_result": "ATC dominates the top ligand-receptor activity rows, but this is spatial-signature only and sample-limited.",
        "next_gate": "Cross-check against the 2026 JCI Insight atlas and avoid ICI-response claims without treated cohorts.",
    })

    # T06 — NIS/SLC5A5 non-methylation split.
    topic = "T06 NIS/SLC5A5 non-methylation silencing exception"
    for row in methyl_gene.itertuples(index=False):
        add_metric(metrics, "T06", topic, f"{row.gene} promoter methylation DM1-DM2", {"n_a": row.n_DM1, "n_b": row.n_DM2, "diff_a_minus_b": row.Δβ, "cohen_d_a_minus_b": row.cohens_d, "p": row.mw_p}, "HM450 per-gene methylation")
    slc = methyl_gene[methyl_gene["gene"] == "SLC5A5"].iloc[0]
    tpo = methyl_gene[methyl_gene["gene"] == "TPO"].iloc[0]
    add_metric(metrics, "T06", topic, "TPO methylation effect minus SLC5A5 exception", {"tpo_d": tpo.cohens_d, "slc5a5_d": slc.cohens_d, "tpo_p": tpo.mw_p, "slc5a5_p": slc.mw_p, "d_gap": tpo.cohens_d - slc.cohens_d}, "mechanism split")
    g, score = grade("moderate", "partial", "medium")
    topics.append({
        "topic_id": "T06", "topic": topic, "rank_hint": 6, "grade": g, "score": score,
        "headline": "RAI-lineage silencing splits into methylation-driven hormone-biosynthesis loss and non-methylation NIS/SLC5A5 loss.",
        "pilot_result": "TPO/TG/TSHR/PAX8 show strong methylation effects; SLC5A5/NIS is the outlier with nonsignificant promoter methylation.",
        "next_gate": "Add histone/miRNA/TF motif evidence; this is mechanism-follow-up, not standalone clinical biomarker yet.",
    })

    # T07 — immune/TLS protective paradox.
    topic = "T07 immune/TLS protective paradox in thyroid cancer"
    for module in ["M1_M2", "CD8", "TLS", "CAF", "RAI"]:
        x = cox[(cox["endpoint"] == "DSS") & (cox["module"] == module)]
        if not x.empty:
            add_metric(metrics, "T07", topic, f"DSS Cox {module}", {"n": int(x.iloc[0]["n"]), "HR": float(x.iloc[0]["HR"]), "p": float(x.iloc[0]["p"])}, "univariate TCGA module Cox")
    for module in ["M1_M2", "CD8", "TLS"]:
        x = cox_adj[(cox_adj["module"] == module) & (cox_adj["adjust"] == "age+stage")]
        if not x.empty:
            add_metric(metrics, "T07", topic, f"DSS Cox age+stage {module}", {"n": int(x.iloc[0]["n"]), "HR": float(x.iloc[0]["HR"]), "p": float(x.iloc[0]["p"]), "HR_lo": float(x.iloc[0]["HR_lo"]), "HR_hi": float(x.iloc[0]["HR_hi"])}, "adjusted TCGA module Cox")
    g, score = grade("moderate", "ready", "medium")
    topics.append({
        "topic_id": "T07", "topic": topic, "rank_hint": 7, "grade": g, "score": score,
        "headline": "Immune-rich/TLS-like thyroid cancers may represent a protective indolence axis rather than a generic aggressive inflammation axis.",
        "pilot_result": "DSS HRs for M1_M2/CD8 are protective in TCGA; TLS trends protective but needs morphology-level validation.",
        "next_gate": "Validate in independent survival cohort and separate Hashimoto-like immune context from tumor immune escape.",
    })

    # T08 — DICER1/EIF1AX rare driver class.
    topic = "T08 DICER1/EIF1AX rare-driver FVPTC-like residual class"
    driver["is_DM2"] = driver["v17_dark_cluster"].eq("DM2")
    driver["is_class4"] = driver["driver_class"].eq("Class4_DICER1_EIF1AX")
    d4 = driver[driver["driver_class"].isin(["Class4_DICER1_EIF1AX", "Class6_True_driver_neg"])].copy()
    add_metric(metrics, "T08", topic, "Class4 vs Class6 DM2 rate", fisher_flag(d4, "driver_class", "Class4_DICER1_EIF1AX", "Class6_True_driver_neg", "is_DM2"), "TCGA driver class")
    for val in ["rai_score_v17", "tds16_score_v17", "mean_8g_beta", "age"]:
        add_metric(metrics, "T08", topic, f"Class4 vs Class6 {val}", cont_metric(d4, "driver_class", "Class4_DICER1_EIF1AX", "Class6_True_driver_neg", val), "TCGA driver class")
    g, score = grade("moderate", "partial", "high")
    topics.append({
        "topic_id": "T08", "topic": topic, "rank_hint": 8, "grade": g, "score": score,
        "headline": "DICER1/EIF1AX-mutant PTC forms a rare FVPTC-like differentiated residual class that may bridge genomics and miRNA-biogenesis biology.",
        "pilot_result": "TCGA n=7 is small but heavily DM2-like with high RAI/TDS and low methylation; K2 rare-driver evidence exists but needs harmonization.",
        "next_gate": "Pool TCGA + K2 + GeoMx DICER1 lesion data; do not pitch as standalone until n improves.",
    })

    # T09 — TERT-only microclass.
    topic = "T09 TERT-only triple-negative high-risk microclass"
    d5 = driver[driver["driver_class"].isin(["Class5_TERT_only", "Class6_True_driver_neg"])].copy()
    d5["PFI_flag"] = d5["PFI"].fillna(0).astype(bool)
    add_metric(metrics, "T09", topic, "Class5 TERT-only vs Class6 PFI event rate", fisher_flag(d5, "driver_class", "Class5_TERT_only", "Class6_True_driver_neg", "PFI_flag"), "TCGA sparse pilot")
    for val in ["rai_score_v17", "mean_8g_beta", "age"]:
        add_metric(metrics, "T09", topic, f"Class5 vs Class6 {val}", cont_metric(d5, "driver_class", "Class5_TERT_only", "Class6_True_driver_neg", val), "TCGA sparse pilot")
    g, score = grade("weak", "partial", "high")
    topics.append({
        "topic_id": "T09", "topic": topic, "rank_hint": 11, "grade": g, "score": score,
        "headline": "A tiny TERT-only residual class may carry disproportionate recurrence risk.",
        "pilot_result": "PFI event rate appears high but TCGA n=5 is too small for a paper.",
        "next_gate": "External TERT-rich cohort or MSK/GENIE-level expansion required; hold as section only.",
    })

    # T10 — histology-invisible molecular axis.
    topic = "T10 morphology-invisible RAI/DM molecular axis"
    target_rows = closure[closure["section"].isin(["A", "B"])].copy()
    best_abs = target_rows["pooled_spearman_r"].abs().max()
    best_auc = target_rows[["pooled_auroc_q75", "pooled_auroc_q50"]].max().max()
    add_metric(metrics, "T10", topic, "Closure battery best absolute Spearman", {"n_models": int(target_rows.shape[0]), "max_abs_spearman": float(best_abs), "max_auroc": float(best_auc)}, "H&E/Visium closure battery")
    for row in target_rows.head(8).itertuples(index=False):
        add_metric(metrics, "T10", topic, f"{row.section} {row.experiment} {row.model}", {"pooled_n": row.pooled_n, "spearman_r": row.pooled_spearman_r, "auroc_q75": row.pooled_auroc_q75, "auroc_q50": row.pooled_auroc_q50}, "negative result")
    g, score = grade("negative", "ready", "medium")
    topics.append({
        "topic_id": "T10", "topic": topic, "rank_hint": 10, "grade": g, "score": score,
        "headline": "Routine morphology does not recover the RAI/DM transcriptomic axis, arguing for molecular readouts rather than H&E surrogates.",
        "pilot_result": "Closure battery stays near random (max |Spearman| < 0.1; AUROC near 0.55).",
        "next_gate": "Use as negative-control/methods paper or supplement; not a wet-lab biology headline.",
    })

    # T11 — ATC spatial coherence collapse/super-organization.
    topic = "T11 ATC spatial coherence collapse with single-sample super-organization"
    for axis in ["Thyroid_differen", "Cell_cycle", "T_cell", "Macrophage_TAM", "Hypoxia", "EMT_program", "Cellular_stress"]:
        add_metric(metrics, "T11", topic, f"ATC vs LPTC stage-gradient {axis}", cont_metric(spatial_v12, "condition", "ATC", "LPTC", axis), "16-slide GSE250521 spatial v12")
    atc = spatial_v12[spatial_v12["condition"] == "ATC"].copy()
    if not atc.empty:
        high_axes = ["Cell_cycle", "T_cell", "Macrophage_TAM", "Hypoxia", "EMT_program", "Cellular_stress"]
        atc["super_org_mean"] = atc[high_axes].mean(axis=1)
        best = atc.sort_values("super_org_mean", ascending=False).iloc[0]
        add_metric(metrics, "T11", topic, "top ATC super-organization sample", {"sample": best["condition"], "super_org_mean": float(best["super_org_mean"]), "n_atc": int(atc.shape[0])}, "case-like outlier screen")
    g, score = grade("moderate", "partial", "high")
    topics.append({
        "topic_id": "T11", "topic": topic, "rank_hint": 9, "grade": g, "score": score,
        "headline": "ATC progression collapses thyroid differentiation spatial coherence while creating stress/immune/hypoxia super-organized regions.",
        "pilot_result": "Spatial v12 shows differentiation collapse in ATC and an outlier-like highly organized ATC sample, but n=4 ATC and n=1 outlier caveats dominate.",
        "next_gate": "Keep supplementary unless an independent ATC spatial cohort or GeoMx PDTC/ATC validation supports it.",
    })

    metrics_df = pd.DataFrame(metrics)
    topics_df = pd.DataFrame(topics).sort_values(["score", "rank_hint"], ascending=[False, True])
    metrics_df.to_csv(OUT / "pilot_metrics_long.tsv", sep="\t", index=False)
    topics_df.to_csv(OUT / "topic_scorecard.tsv", sep="\t", index=False)

    # Markdown report.
    lines: list[str] = []
    lines.append("# 11개 고IF 후보 주제 pilot screen — 2026-05-08\n")
    lines.append("범위: 현재 CNV residual class 포함 + 추가 10개 후보. 기존 원고/figure는 수정하지 않았고, 로컬 public-data-derived 결과만 사용.\n")
    lines.append("## 최종 순위\n")
    lines.append("| Rank | Topic | Grade | Score | 한 줄 판정 | Next gate |")
    lines.append("|---:|---|---:|---:|---|---|")
    for i, row in enumerate(topics_df.itertuples(index=False), start=1):
        lines.append(f"| {i} | {row.topic_id} · {row.topic} | {row.grade} | {row.score} | {row.headline} | {row.next_gate} |")
    lines.append("\n## Topic별 pilot 요약\n")
    for row in topics_df.itertuples(index=False):
        lines.append(f"### {row.topic_id} — {row.topic}\n")
        lines.append(f"**Headline:** {row.headline}\n")
        lines.append(f"**Pilot result:** {row.pilot_result}\n")
        lines.append(f"**Next gate:** {row.next_gate}\n")
        sub = metrics_df[metrics_df["topic_id"] == row.topic_id].head(12)
        if not sub.empty:
            lines.append("\n| Metric | n/rate/effect | p | Note |")
            lines.append("|---|---|---:|---|")
            for m in sub.itertuples(index=False):
                bits = []
                for key in ["n", "n_a", "n_b", "rate_a", "rate_b", "mean_a", "mean_b", "diff_a_minus_b", "cohen_d_a_minus_b", "spearman_r", "median_r", "HR", "max_abs_spearman", "max_auroc"]:
                    if hasattr(m, key):
                        val = getattr(m, key)
                        if pd.notna(val):
                            if isinstance(val, float):
                                bits.append(f"{key}={val:.3g}")
                            else:
                                bits.append(f"{key}={val}")
                p = getattr(m, "p", np.nan) if hasattr(m, "p") else np.nan
                lines.append(f"| {m.metric} | {'; '.join(bits) if bits else '-'} | {p_text(p)} | {m.note} |")
        lines.append("")
    lines.append("## 산출 파일\n")
    lines.append(f"- `{OUT / 'topic_scorecard.tsv'}`\n")
    lines.append(f"- `{OUT / 'pilot_metrics_long.tsv'}`\n")
    lines.append(f"- `{REPORT}`\n")
    REPORT.write_text("\n".join(lines), encoding="utf-8")

    print(f"[write] {OUT / 'topic_scorecard.tsv'}")
    print(f"[write] {OUT / 'pilot_metrics_long.tsv'}")
    print(f"[write] {REPORT}")


if __name__ == "__main__":
    main()
