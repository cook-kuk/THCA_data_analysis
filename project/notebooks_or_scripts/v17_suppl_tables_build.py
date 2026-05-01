#!/usr/bin/env python3
"""Build paper-ready supplementary TSV files (S1-S8).

Voice-protected sections skipped. Fact-only data assembly.
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd

PROJ = Path("/opt/thyroid-dash/project")
OUT = PROJ / "submission_data"
OUT.mkdir(parents=True, exist_ok=True)

# ============================================================
# S1 — Cohort assembly
# ============================================================
print("=== S1 — Cohort assembly ===")
s1 = pd.DataFrame([
    {"cohort_id": "TCGA-THCA", "source": "TCGA Cancer Network 2014",
     "n_total": 500, "n_DM1": 140, "n_DM2": 360,
     "modality": "Illumina RNA-seq (HiSeq + GA)",
     "population": "EUR-dominant (mixed)",
     "HLA_avail": "HLA-LA imputed 4-digit",
     "use_case": "Discovery, classifier training, Pillars 1-5"},
    {"cohort_id": "K2 (PRJEB11591)", "source": "Yoo SK 2016 Mol Cell Biol; SNU-GMI",
     "n_total": 260, "n_DM1": 14, "n_DM2": 246,
     "modality": "Illumina RNA-seq",
     "population": "Korean",
     "HLA_avail": "arcasHLA 4-digit (235 valid, 5/2 burst)",
     "use_case": "Korean PTC validation, Pillar 1"},
    {"cohort_id": "Lee 2024 (GSE213647)", "source": "Lee Y et al. 2024",
     "n_total": 632, "n_DM1": "proxy via 8-gene", "n_DM2": "proxy",
     "modality": "Illumina RNA-seq (Macrogen)",
     "population": "Korean",
     "HLA_avail": "arcasHLA 4-digit (630 valid)",
     "use_case": "Korean PTC replication, Hashimoto-like generalization, Pillars 1+5"},
    {"cohort_id": "GSE286332", "source": "Lim DW 2025; Dongguk Univ; PMID 41113708",
     "n_total": 18, "n_DM1": 0, "n_DM2": 18,
     "modality": "Illumina NovaSeq X (Macrogen)",
     "population": "Korean",
     "HLA_avail": "arcasHLA 4-digit (18/18, 5/2 burst)",
     "use_case": "PTC vs PTC+HT discovery, Pillars 2+5"},
    {"cohort_id": "Chu 2018 (PMC 6161647)", "source": "Chu X et al. 2018 J Med Genet 55(10):685-692",
     "n_total": 2958, "n_DM1": "NA (GD vs ctrl)", "n_DM2": "NA",
     "modality": "SNP2HLA Pan-Asian reference panel",
     "population": "Han Chinese",
     "HLA_avail": "Published summary statistics (4-digit)",
     "use_case": "Pan-Asian forest meta replication, Pillar 1"},
])
s1.to_csv(OUT / "S1_cohort_assembly.tsv", sep="\t", index=False)
print(f"  ✓ S1: {len(s1)} cohorts")

# ============================================================
# S2 — TIERA67 7-category 67-gene candidate pool
# ============================================================
print("\n=== S2 — TIERA67 candidate pool ===")
s2_rows = []
categories = {
    "TDS_core": ["DIO1", "DIO2", "DUOX1", "DUOX2", "FOXE1", "GLIS3", "NKX2-1", "PAX8",
                  "SLC26A4", "SLC5A5", "SLC5A8", "TG", "THRA", "THRB", "TPO", "TSHR"],
    "MAPK_output_ERK": ["DUSP4", "DUSP5", "DUSP6", "SPRY1", "SPRY2", "SPRY4",
                         "ETV4", "ETV5", "PHLDA1", "FOSL1"],
    "Driver_anchor": ["BRAF", "NRAS", "HRAS", "KRAS", "RET", "NTRK1", "NTRK3", "ALK",
                       "PAX8", "PPARG", "TERT", "EIF1AX"],
    "Aggressive_marker": ["TP53", "CDKN2A", "CDKN2B", "PIK3CA", "AKT1", "PTEN",
                           "ATM", "CTNNB1", "APC", "MSH2"],
    "Dediff_invasion": ["VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1",
                         "CDH1", "CDH2", "MMP9", "LOX"],
    "Immune_stromal_light": ["CD274", "CD8A", "FOXP3", "IDO1", "HLA-DRA"],
    "Thyroid_lineage_extra": ["IYD", "THADA", "MET", "KLK10"],
}
GENE_8 = ['SLC5A5', 'TPO', 'TG', 'TSHR', 'PAX8', 'NKX2-1', 'FOXE1', 'DIO1']
for cat, genes in categories.items():
    for g in genes:
        s2_rows.append({"category": cat, "gene": g, "in_8gene_panel": g in GENE_8})
s2 = pd.DataFrame(s2_rows)
s2.to_csv(OUT / "S2_TIERA67_candidate_pool.tsv", sep="\t", index=False)
print(f"  ✓ S2: {len(s2)} genes; 8-gene = {s2['in_8gene_panel'].sum()}")

# ============================================================
# S3 — GSE286332 top 300 DEGs
# ============================================================
print("\n=== S3 — GSE286332 top 300 DEGs ===")
deg = pd.read_csv(PROJ / "results/p3_gse286332/deg_ptcht_vs_ptc.tsv", sep="\t")
deg_clean = deg.dropna(subset=["padj", "log2FoldChange"]).copy()
deg_clean = deg_clean.sort_values("padj").head(300)
deg_clean = deg_clean[["gene", "baseMean", "log2FoldChange", "stat", "pvalue", "padj"]]
deg_clean["category_predicted"] = deg_clean["gene"].apply(lambda g:
    "Ig V/J/C" if g.startswith(("IGHV", "IGKV", "IGLV", "IGHJ", "IGKJ", "IGLJ", "IGHC")) else
    "HLA-II" if g.startswith("HLA-D") else
    "B-cell" if g in {"BLK", "MS4A1", "CD19", "CD79A", "CD79B", "BANK1"} else
    "T-cell" if g in {"CD3D", "CD3E", "CD8A", "CD4", "EOMES"} else
    "Other immune" if "CXCL" in g or "CCL" in g or "STAT" in g else
    "Other")
deg_clean.to_csv(OUT / "S3_gse286332_top300_DEGs.tsv", sep="\t", index=False)
print(f"  ✓ S3: top 300 DEGs (full 29,672 in main DEG file)")

# ============================================================
# S4 — Pan-Asian HLA per-allele forest (with Korean sub-cohort)
# ============================================================
print("\n=== S4 — Pan-Asian HLA forest ===")
forest = pd.read_csv(PROJ / "results/p2_pillar1_forest/forest_meta_results.tsv", sep="\t")
pooled = pd.read_csv(PROJ / "results/p2_pillar1_forest/random_effects_pooled.tsv", sep="\t")
sub = pd.read_csv(PROJ / "results/p2_pillar1_forest/korean_PTC_pool_per_subcohort.tsv", sep="\t")
het = pd.read_csv(PROJ / "results/p2_pillar1_forest/korean_subcohort_heterogeneity.tsv", sep="\t")

# Merge per-allele
s4 = forest.merge(pooled[["allele", "pooled_OR", "pooled_ci_lo", "pooled_ci_hi", "I2_pct", "tau2"]],
                    on="allele", how="left")
s4 = s4.merge(het[["allele", "K2", "Lee", "GSE286", "Cochran_Q", "I2_pct"]].rename(
    columns={"I2_pct": "Korean_I2_pct"}), on="allele", how="left")
s4.to_csv(OUT / "S4_panasian_HLA_forest_full.tsv", sep="\t", index=False)
print(f"  ✓ S4: {len(s4)} alleles × full forest + Korean sub-cohort + Cochran's Q")

# ============================================================
# S5 — Mediation analysis full (HLA-II, g8_RAI, immune, HLA-I)
# ============================================================
print("\n=== S5 — Mediation analysis ===")
med_full = json.loads((PROJ / "results/d3p5_pdm1_gradient/mediation_results.json").read_text())
s5_rows = []
for med, vals in med_full.items():
    s5_rows.append({"mediator": med, **vals})
s5 = pd.DataFrame(s5_rows)
s5.to_csv(OUT / "S5_mediation_analysis.tsv", sep="\t", index=False)
print(f"  ✓ S5: {len(s5)} mediators (HLA-II 140%, g8_RAI 63%, immune 88%, HLA-I 87%)")

# Also export OLS decomposition
decomp = json.loads((PROJ / "results/d3p5_pdm1_gradient/decomposition.json").read_text())
ols_rows = []
for k, v in decomp["coefficients"].items():
    ols_rows.append({"predictor": k, **v})
ols_df = pd.DataFrame(ols_rows)
ols_df.to_csv(OUT / "S5b_OLS_decomposition.tsv", sep="\t", index=False)
print(f"  ✓ S5b: OLS R²={decomp['R2']}")

# ============================================================
# S6 — TCGA Hashimoto-like × DM cluster cross-tabulation
# ============================================================
print("\n=== S6 — TCGA Hashimoto-like × DM ===")
d4p2 = json.loads((PROJ / "results/d4p2_tcga_hashimoto_signature/D4P2_summary.json").read_text())
s6_rows = []
for method, ct in d4p2["crosstabs"].items():
    s6_rows.append({
        "threshold_method": method,
        "n_hashi_pos": d4p2["hashi_calls"][method.replace("hashi_", "").replace("resid_otsu", "resid_otsu")] if method.replace("hashi_", "") in d4p2["hashi_calls"] else None,
        "DM1_hashi_pct": ct["DM1_hashi_pct"],
        "DM2_hashi_pct": ct["DM2_hashi_pct"],
        "OR": ct["OR"],
        "Fisher_p": ct["fisher_p"],
    })
s6 = pd.DataFrame(s6_rows)
s6.to_csv(OUT / "S6_tcga_hashimoto_DM_crosstab.tsv", sep="\t", index=False)
print(f"  ✓ S6: {len(s6)} thresholds")

# ============================================================
# S7 — Cross-cohort generalization (TCGA + Korean replication + sub-B)
# ============================================================
print("\n=== S7 — Cross-cohort replication ===")
d8b = json.loads((PROJ / "results/d8b_korean_replication/D8B_summary.json").read_text())
d8c = json.loads((PROJ / "results/d8c_dm1_subB_x_K2_NBNR/D8C_summary.json").read_text())
s7 = pd.DataFrame([
    {"cohort": "TCGA-THCA", "n": 500,
     "Hashimoto_GMM_pct": 18.0, "Hashimoto_Otsu_pct": 19.6,
     "subB_GMM_pct": "NA (DM1-only re-cluster)", "subB_Otsu_pct": "NA"},
    {"cohort": "Korean GSE213647 (Lee 2024)", "n": 632,
     "Hashimoto_GMM_pct": d8b["hashi_pct"]["GMM"], "Hashimoto_Otsu_pct": d8b["hashi_pct"]["Otsu"],
     "subB_GMM_pct": d8c["Korean_subB_rate_GMM_pct"], "subB_Otsu_pct": d8c["Korean_subB_rate_Otsu_pct"]},
    {"cohort": "GSE286332 (PTC+HT only)", "n": 9,
     "Hashimoto_GMM_pct": 100.0, "Hashimoto_Otsu_pct": 100.0,
     "subB_GMM_pct": "NA", "subB_Otsu_pct": "NA"},
])
s7.to_csv(OUT / "S7_cross_cohort_generalization.tsv", sep="\t", index=False)
print(f"  ✓ S7: {len(s7)} cohorts")

# ============================================================
# S8 — DM1 sub-cluster A vs B characterization
# ============================================================
print("\n=== S8 — DM1 sub-A vs sub-B ===")
d6p7 = json.loads((PROJ / "results/d6p7_dm1_subcluster/D6P7_summary.json").read_text())
sub_score = pd.DataFrame(d6p7["score_profile"])
sub_score.to_csv(OUT / "S8a_dm1_subcluster_score_profile.tsv", sep="\t", index=False)
sub_clin = pd.DataFrame(d6p7.get("clinical", []))
sub_clin.to_csv(OUT / "S8b_dm1_subcluster_clinical.tsv", sep="\t", index=False)

# Mutation × sub-cluster
sub_label = pd.read_csv(PROJ / "results/d6p7_dm1_subcluster/dm1_subcluster_labels.tsv", sep="\t", index_col=0)
muts = pd.read_csv(PROJ / "results/tables/tcga_thca_mutation_groups.tsv", sep="\t", index_col=0)
hashi = pd.read_csv(PROJ / "results/d4p2_tcga_hashimoto_signature/tcga_signature_scores.tsv", sep="\t", index_col=0)
join = sub_label.join(muts[["has_braf_v600e", "has_ras_mut"]]).join(hashi[["hashi_otsu"]])
mut_xtab = pd.DataFrame({
    "sub_cluster": ["sub_A", "sub_B"],
    "n": [(join["sub_cluster"]=="sub_A").sum(), (join["sub_cluster"]=="sub_B").sum()],
    "BRAF_pos": [int(((join["sub_cluster"]=="sub_A") & (join["has_braf_v600e"]==1)).sum()),
                  int(((join["sub_cluster"]=="sub_B") & (join["has_braf_v600e"]==1)).sum())],
    "RAS_pos": [int(((join["sub_cluster"]=="sub_A") & (join["has_ras_mut"]==1)).sum()),
                 int(((join["sub_cluster"]=="sub_B") & (join["has_ras_mut"]==1)).sum())],
    "Hashimoto_pos_Otsu": [int(((join["sub_cluster"]=="sub_A") & (join["hashi_otsu"]==1)).sum()),
                            int(((join["sub_cluster"]=="sub_B") & (join["hashi_otsu"]==1)).sum())],
})
mut_xtab["mut_neg_pct"] = (1 - (mut_xtab["BRAF_pos"] + mut_xtab["RAS_pos"]) / mut_xtab["n"]) * 100
mut_xtab["Hashimoto_pct"] = mut_xtab["Hashimoto_pos_Otsu"] / mut_xtab["n"] * 100
mut_xtab.to_csv(OUT / "S8c_dm1_subcluster_mutation_hashimoto.tsv", sep="\t", index=False)
print(f"  ✓ S8a/b/c: score profile + clinical + mutation×Hashi")

# ============================================================
# Summary
# ============================================================
print("\n=== Summary ===")
files = sorted(OUT.glob("S*.tsv"))
for f in files:
    print(f"  {f.name}: {f.stat().st_size:>8d} bytes")
print(f"\n✓ Total {len(files)} suppl TSVs at {OUT}")
