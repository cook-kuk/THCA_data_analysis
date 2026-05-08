#!/usr/bin/env python3
"""Build curated figure inventory TSV for spatial-full 100 figures.
Output: project/reports/spatial_figure_curation_2026_05_07.tsv
"""
from pathlib import Path
import csv

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT/"project/reports/spatial_figure_curation_2026_05_07.tsv"
OUT.parent.mkdir(exist_ok=True, parents=True)

# version source map by figure id
VER = {}
for n in range(1, 7): VER[n]="v1"
for n in range(7,19): VER[n]="v2"
for n in range(19,27): VER[n]="v3"
for n in range(27,33): VER[n]="v4"
for n in range(33,41): VER[n]="v5"
for n in range(41,49): VER[n]="v6"
for n in range(49,57): VER[n]="v7"
for n in range(57,65): VER[n]="v8"
for n in range(65,73): VER[n]="v9"
for n in range(73,81): VER[n]="v10"
for n in range(81,89): VER[n]="v11"
for n in range(89,101): VER[n]="v12"
for n in range(101,107): VER[n]="v13_geomx"

# cohort default per figure id (overridden where mixed)
def default_cohort(n):
    if n in (3,8,18,20,24,29,37,42,53,58,62,65): return "GSE230424"
    if n in (4,15,36): return "GSE248205"
    if n in (6,11,16,22,26,38,40,41,47,49,55,56,63,64,66,67,68,71,99,100): return "ALL_28_3cohort"
    if n in (13,30,61): return "GSE250521+TCGA_h&e"
    if n in (101,102,103,104,105): return "GSE301163_GeoMx"
    if n == 106: return "Visium_28+GeoMx_78_crossplatform"
    return "GSE250521"

# Hand-curated metadata for each figure
# (axis, paper_role, tier, claim_strength, risk_level, include_reason, exclude_reason, sentence_draft, needs_recalc, notes)
M = {}
def set_(n, axis, role, tier, claim, risk, inc, exc, sent, recalc="NO", notes=""):
    M[n] = (axis, role, tier, claim, risk, inc, exc, sent, recalc, notes)

# --- v1 (S_F1-S_F6) — foundation
set_(1, "DM1_axis", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "Early DM1 niche density visualization", "Superseded by S_F19 stage representative DM1 maps",
     "[AUTHOR VOICE NEEDED: DM1 niche 2D density overview]")
set_(2, "Driver_orthogonal", "Paper2_Pillar1_support", "SUPPLEMENTARY_FIGURE", "MODERATE_SUPPORT", "LOW",
     "Driver-orthogonal axis exists at spot level", "",
     "[AUTHOR VOICE NEEDED: Driver-orthogonal axis box]")
set_(3, "HT_TLS", "Paper2_Pillar2_main", "WEB_ONLY_ARCHIVE", "MODERATE_SUPPORT", "LOW",
     "Early HT TLS panel; superseded by S_F8 16-slide layout", "",
     "[AUTHOR VOICE NEEDED: HT TLS preview]")
set_(4, "Autoimmune_baseline", "Paper2_neg_control", "WEB_ONLY_ARCHIVE", "MODERATE_SUPPORT", "LOW",
     "Early GSE248205 baseline; superseded by S_F36", "",
     "[AUTHOR VOICE NEEDED: GSE248205 negative-control baseline]")
set_(5, "TROP2_niche", "Paper2_Pillar1_support", "SUPPLEMENTARY_FIGURE", "MODERATE_SUPPORT", "LOW",
     "Per-sample TROP2 spot distribution (full cohort)", "",
     "[AUTHOR VOICE NEEDED: TROP2 spot distribution per sample]")
set_(6, "Cross_cohort_TROP2", "Paper2_Pillar1_support", "WEB_ONLY_ARCHIVE", "MODERATE_SUPPORT", "LOW",
     "Earlier 3-cohort integration; superseded by S_F11", "",
     "[AUTHOR VOICE NEEDED: cross-cohort TROP2 overview]")

# --- v2 (S_F7-S_F18) — TROP2 niche identification core + HT TLS core
set_(7, "TROP2_niche", "Paper2_Pillar1_main", "MAIN_FIGURE", "STRONG_SUPPORT", "LOW",
     "Visual evidence: TROP2 spatially organized in 8/8 PTC+LPTC slides", "",
     "[AUTHOR VOICE NEEDED: TROP2 16-slide spatial maps demonstrate tumor-specific niche organization in GSE250521]")
set_(8, "HT_TLS", "Paper2_Pillar2_main", "MAIN_FIGURE", "STRONG_SUPPORT", "MEDIUM",
     "HT-overlap PTC TLS niche directly visible 4/4 slides (n=4 risk)", "",
     "[AUTHOR VOICE NEEDED: HT TLS 4-slide spatial maps support Pillar II HT immune niche]",
     notes="n=4 PTC+HT — flag in reviewer risk")
set_(9, "TROP2_niche", "Paper2_Pillar1_support", "SUPPLEMENTARY_FIGURE", "STRONG_SUPPORT", "LOW",
     "Per-sample ranked Moran's I quantification of S_F7", "",
     "[AUTHOR VOICE NEEDED: ranked Moran's I per sample (post-fix)]",
     recalc="DONE_2026_05_07", notes="REMADE 2026-05-07 with hex 6-nbr; values now from spatial_morans_v7plus_2026_05_07.tsv")
set_(10, "TROP2_x_DM1", "Paper2_Pillar1_support", "SUPPLEMENTARY_FIGURE", "MODERATE_SUPPORT", "LOW",
     "Spatial co-localization TROP2 × DM1 axis", "",
     "[AUTHOR VOICE NEEDED: TROP2-DM1 spatial co-localization]")
set_(11, "Cross_cohort_TROP2", "Paper2_Pillar1_main", "MAIN_FIGURE", "STRONG_SUPPORT", "LOW",
     "Quantitative: TROP2 niche tumor-specific across 28 samples / 3 cohorts (post-fix)", "",
     "[AUTHOR VOICE NEEDED: cross-cohort 28×12 axis Moran's I heatmap (post-fix hex 6-nbr)]",
     recalc="DONE_2026_05_07", notes="REMADE 2026-05-07; PTC 0.33-0.58, LPTC 0.25-0.57 confirmed; PTC_HT TROP2 0.02-0.19 (lower) flagged")
set_(12, "DM1_axis", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "Representative DM1 maps; superseded by S_F19 16-slide", "",
     "[AUTHOR VOICE NEEDED: DM1 representative maps]")
set_(13, "Closure_battery", "Paper2_Pillar2_support", "SUPPLEMENTARY_FIGURE", "STRONG_SUPPORT", "LOW",
     "Closure battery summary (NO-GO record)", "",
     "[AUTHOR VOICE NEEDED: H&E → DM1 closure battery LOSO summary]")
set_(14, "TROP2_niche", "Paper2_Pillar1_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "Per-stage violin; redundant with S_F11 heatmap", "",
     "[AUTHOR VOICE NEEDED: TROP2 violin per stage]")
set_(15, "Autoimmune_baseline", "Paper2_neg_control", "WEB_ONLY_ARCHIVE", "MODERATE_SUPPORT", "LOW",
     "GSE248205 baseline maps; superseded by S_F36", "",
     "[AUTHOR VOICE NEEDED: autoimmune-only baseline maps]")
set_(16, "Marker_x_condition", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "MODERATE_SUPPORT", "LOW",
     "Fixed (v7) 2-panel; redundant with S_F11 + S_F18", "Redundant after main + supp pair",
     "[AUTHOR VOICE NEEDED: marker × condition Moran summary]",
     notes="originally broken; fixed in v7; demoted to web")
set_(17, "TROP2_niche", "Paper2_Pillar1_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "TROP2 high-spot fraction; redundant with S_F11", "",
     "[AUTHOR VOICE NEEDED: TROP2 high-fraction]")
set_(18, "HT_TLS", "Paper2_Pillar2_main", "MAIN_FIGURE", "STRONG_SUPPORT", "MEDIUM",
     "HT TLS Moran I heatmap; quantitative companion to S_F8 (post-fix)", "",
     "[AUTHOR VOICE NEEDED: HT TLS Moran's I 4 slides × 4 immune axes (post-fix hex 6-nbr); HLA_II 0.48-0.81, IGHV 0.40-0.85, TLS 0.28-0.61]",
     recalc="DONE_2026_05_07", notes="REMADE 2026-05-07; n=4 risk shared with S_F8")

# --- v3 (S_F19-S_F26)
set_(19, "DM1_axis", "Paper2_Pillar2_main", "MAIN_FIGURE", "STRONG_SUPPORT", "LOW",
     "Spatial demonstration of DM1 region heterogeneity (signature-level)", "",
     "[AUTHOR VOICE NEEDED: DM1 16-slide spatial maps demonstrate regional heterogeneity in GSE250521]")
set_(20, "DM1_axis_HT", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "MODERATE_SUPPORT", "LOW",
     "DM1 axis on PTC+HT cohort", "Redundant with S_F19 main + S_F49 supp",
     "[AUTHOR VOICE NEEDED: DM1 maps on GSE230424]")
set_(21, "Microenv_TROP2", "Paper2_Pillar1_support", "WEB_ONLY_ARCHIVE", "MODERATE_SUPPORT", "LOW",
     "TROP2 microenvironment context", "Redundant with S_F10 colocalization",
     "[AUTHOR VOICE NEEDED: TROP2 microenv]")
set_(22, "Moran_per_condition", "Paper2_Pillar1_support", "SUPPLEMENTARY_FIGURE", "STRONG_SUPPORT", "LOW",
     "Moran's I distribution per condition (post-fix)", "",
     "[AUTHOR VOICE NEEDED: TROP2 Moran's I per condition box (post-fix hex 6-nbr)]",
     recalc="DONE_2026_05_07", notes="REMADE 2026-05-07")
set_(23, "TROP2_distribution", "Paper2_Pillar1_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "TROP2 distribution histogram", "",
     "[AUTHOR VOICE NEEDED: TROP2 hist]")
set_(24, "HT_immune_coloc", "Paper2_Pillar2_support", "SUPPLEMENTARY_FIGURE", "MODERATE_SUPPORT", "LOW",
     "HT immune marker co-localization", "",
     "[AUTHOR VOICE NEEDED: HT immune co-localization]")
set_(25, "TROP2_clusters", "Paper2_Pillar1_support", "SUPPLEMENTARY_FIGURE", "TECHNICAL_QC", "LOW",
     "Niche cluster sizes — Visium hex 6-nbr fixed (v7)", "",
     "[AUTHOR VOICE NEEDED: TROP2 niche cluster size distribution (hex 6-neighbor adjacency)]",
     notes="post-fix; original v3 version DEPRECATED")
set_(26, "Mega_summary", "Paper2_Pillar1_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "Mega summary (v3 era); superseded by S_F56/S_F64/S_F100", "",
     "[AUTHOR VOICE NEEDED: mega summary v3]")

# --- v4 (S_F27-S_F32)
set_(27, "RGB_overlay", "Paper2_Pillar1_support", "WEB_ONLY_ARCHIVE", "MODERATE_SUPPORT", "LOW",
     "TROP2 + DM1 + Epithelial RGB overlay", "Visually appealing but redundant with S_F10",
     "[AUTHOR VOICE NEEDED: RGB tri-axis overlay]")
set_(28, "Driver_genes_spatial", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "MODERATE_SUPPORT", "LOW",
     "Driver gene spatial expression maps", "Driver-orthogonal already captured in S_F2",
     "[AUTHOR VOICE NEEDED: driver gene maps]")
set_(29, "TROP2_HT", "Paper2_Pillar1_support", "SUPPLEMENTARY_FIGURE", "STRONG_SUPPORT", "LOW",
     "TROP2 niche on GSE230424 PTC+HT 4 slides", "",
     "[AUTHOR VOICE NEEDED: TROP2 maps on GSE230424]")
set_(30, "Closure_battery", "Paper2_Pillar2_main", "MAIN_FIGURE", "STRONG_SUPPORT", "LOW",
     "H&E → DM1 closure NO-GO LOSO ResNet50 ρ ≈ 0.06 (fixed v7)", "",
     "[AUTHOR VOICE NEEDED: closure battery LOSO per-slide breakdown shows H&E → DM1 closure NO-GO]",
     notes="schema fix in v7; main NO-GO evidence")
set_(31, "Niche_spacing", "Paper2_Pillar1_support", "SUPPLEMENTARY_FIGURE", "TECHNICAL_QC", "LOW",
     "Niche cluster spacing — REMADE 2026-05-07 with hex 6-nbr union-find", "",
     "[AUTHOR VOICE NEEDED: TROP2 niche cluster pairwise distance distribution (post-fix Visium hex union-find)]",
     recalc="DONE_2026_05_07", notes="REMADE 2026-05-07; Bug A resolved (visium_hex_clusters union-find); promoted from web back to supp")
set_(32, "ATC_collapse", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "MEDIUM",
     "ATC niche degradation visual (n=4 ATC risk)", "Subsumed by S_F91 (v12) which is the supp citation",
     "[AUTHOR VOICE NEEDED: ATC niche degradation overview]",
     notes="n=4 ATC; reviewer risk")

# --- v5 (S_F33-S_F40)
set_(33, "RAI8_individual", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "MODERATE_SUPPORT", "LOW",
     "RAI8 individual gene maps (dedifferentiation context)", "Differentiation captured in S_F91 supp",
     "[AUTHOR VOICE NEEDED: RAI8 per-gene maps]")
set_(34, "Epithelial_axis", "Paper2_Pillar1_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "Epithelial axis 16-slide; supplements TROP2-Epi colocalization", "",
     "[AUTHOR VOICE NEEDED: epithelial maps]")
set_(35, "Proliferation_axis", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "Proliferation 16-slide", "",
     "[AUTHOR VOICE NEEDED: proliferation maps]")
set_(36, "Autoimmune_full", "Paper2_neg_control", "MAIN_FIGURE", "STRONG_SUPPORT", "LOW",
     "Negative control: autoimmune-only cohort lacks TROP2/TLS niche signal", "",
     "[AUTHOR VOICE NEEDED: GSE248205 full maps autoimmune-only negative control]")
set_(37, "HT_combined_immune", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "MODERATE_SUPPORT", "LOW",
     "HT combined immune density", "Captured in S_F8 main + S_F18 supp",
     "[AUTHOR VOICE NEEDED: HT combined immune density]")
set_(38, "Per_sample_QC", "Paper2_methods", "SUPPLEMENTARY_FIGURE", "TECHNICAL_QC", "LOW",
     "Per-sample technical QC (depth, n_spots, fraction expressed)", "",
     "[AUTHOR VOICE NEEDED: per-sample QC supplementary]")
set_(39, "Microenv_trajectory", "Paper2_Pillar1_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "Microenvironment trajectory; redundant", "",
     "[AUTHOR VOICE NEEDED: microenv trajectory]")
set_(40, "Moran_axis_corr", "Paper2_methods", "SUPPLEMENTARY_FIGURE", "TECHNICAL_QC", "LOW",
     "Cross-axis Spearman correlation of Moran's I (post-fix)", "",
     "[AUTHOR VOICE NEEDED: cross-axis Moran's I Spearman correlation across 28 samples (post-fix)]",
     recalc="DONE_2026_05_07", notes="REMADE 2026-05-07")

# --- v6 (S_F41-S_F48)
set_(41, "TROP2_28sample", "Paper2_Pillar1_main", "SUPPLEMENTARY_FIGURE", "STRONG_SUPPORT", "LOW",
     "All 28 samples TROP2 mega view", "",
     "[AUTHOR VOICE NEEDED: TROP2 28-sample mega panel]")
set_(42, "HT_individual", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "MODERATE_SUPPORT", "MEDIUM",
     "HT individual marker maps (CXCL13, CD20, IgG)", "Subsumed by S_F65 (v9) chemokine receptor supp",
     "[AUTHOR VOICE NEEDED: HT individual immune markers]",
     notes="depends on n=4 PTC+HT cohort")
set_(43, "Hypoxia_v6", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "Hypoxia 16-slide v6; superseded by v12 S_F89", "",
     "[AUTHOR VOICE NEEDED: hypoxia spatial maps]")
set_(44, "CAF_ECM", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "CAF/ECM 16-slide", "",
     "[AUTHOR VOICE NEEDED: CAF ECM maps]")
set_(45, "EMT_v6", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "EMT 16-slide v6; superseded by v12 S_F90", "",
     "[AUTHOR VOICE NEEDED: EMT v6 maps]")
set_(46, "Niche_threshold", "Paper2_methods", "SUPPLEMENTARY_FIGURE", "TECHNICAL_QC", "LOW",
     "Niche threshold sensitivity sweep (Q3/Q4/percentile)", "",
     "[AUTHOR VOICE NEEDED: threshold sensitivity supplementary]")
set_(47, "Per_axis_moran_summary", "Paper2_Pillar1_support", "SUPPLEMENTARY_FIGURE", "STRONG_SUPPORT", "LOW",
     "Per-axis Moran summary across 8 conditions (post-fix)", "",
     "[AUTHOR VOICE NEEDED: per-axis Moran's I summary by condition (post-fix hex 6-nbr)]",
     recalc="DONE_2026_05_07", notes="REMADE 2026-05-07; mean ± SD per condition × 12 axes")
set_(48, "Depth_confound", "Paper2_methods", "SUPPLEMENTARY_FIGURE", "TECHNICAL_QC", "LOW",
     "Depth confounding control (Moran I vs sequencing depth)", "",
     "[AUTHOR VOICE NEEDED: depth confounding QC]")

# --- v7 (S_F49-S_F56) — fixes + 4 KEY collage
set_(49, "DM1_28sample", "Paper2_Pillar2_main", "SUPPLEMENTARY_FIGURE", "STRONG_SUPPORT", "LOW",
     "All 28 samples DM1 mega view (companion to S_F19)", "",
     "[AUTHOR VOICE NEEDED: DM1 28-sample mega panel]")
set_(50, "Niche_robustness", "Paper2_methods", "SUPPLEMENTARY_FIGURE", "TECHNICAL_QC", "LOW",
     "Niche robustness post-fix (uses Visium hex 6-nbr)", "",
     "[AUTHOR VOICE NEEDED: niche robustness]")
set_(51, "TROP2_DM1_scatter", "Paper2_Pillar1_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "TROP2 vs DM1 Moran scatter; redundant with S_F71", "",
     "[AUTHOR VOICE NEEDED: TROP2-DM1 scatter]")
set_(52, "Niche_progression", "Paper2_Pillar1_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "Niche progression by stage; superseded by S_F97", "",
     "[AUTHOR VOICE NEEDED: niche progression]")
set_(53, "HT_full_axis_morans", "Paper2_Pillar2_support", "SUPPLEMENTARY_FIGURE", "STRONG_SUPPORT", "MEDIUM",
     "HT cohort all-axis Moran I (TLS, IGHV, AICDA, HLA-II)", "",
     "[AUTHOR VOICE NEEDED: HT full-axis Moran summary]",
     notes="n=4 risk")
set_(54, "TROP2_jaccard", "Paper2_Pillar1_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "TROP2 niche Jaccard between samples; superseded by S_F66", "",
     "[AUTHOR VOICE NEEDED: TROP2 Jaccard]")
set_(55, "Niche_vs_n", "Paper2_methods", "WEB_ONLY_ARCHIVE", "TECHNICAL_QC", "LOW",
     "Niche detection sensitivity vs cohort size", "QC; redundant with S_F46 + S_F48",
     "[AUTHOR VOICE NEEDED: niche vs n_samples QC]")
set_(56, "Four_key_collage", "Paper2_overview", "MAIN_FIGURE", "STRONG_SUPPORT", "LOW",
     "Graphical abstract — answers Q1+Q2+Pillar I+Pillar II in one figure", "",
     "[AUTHOR VOICE NEEDED: four-key-finding collage as graphical abstract]")

# --- v8 (S_F57-S_F64)
set_(57, "Per_sample_violin", "Paper2_Pillar1_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "Per-sample violin; redundant with S_F22/S_F11", "",
     "[AUTHOR VOICE NEEDED: per-sample violin]")
set_(58, "GSE230424_RGB", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "GSE230424 RGB overlay; redundant with S_F27 + S_F8", "",
     "[AUTHOR VOICE NEEDED: GSE230424 RGB]")
set_(59, "TROP2_CDF", "Paper2_Pillar1_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "TROP2 CDF / density; redundant with S_F23", "",
     "[AUTHOR VOICE NEEDED: TROP2 CDF]")
set_(60, "Niche_cluster_IDs", "Paper2_Pillar1_support", "WEB_ONLY_ARCHIVE", "TECHNICAL_QC", "LOW",
     "Niche cluster ID color visualization (post-fix)", "",
     "[AUTHOR VOICE NEEDED: cluster ID visualization]")
set_(61, "Closure_target_model", "Paper2_Pillar2_main", "SUPPLEMENTARY_FIGURE", "STRONG_SUPPORT", "LOW",
     "Closure target × model heatmap (extends S_F30)", "",
     "[AUTHOR VOICE NEEDED: closure target × model heatmap]")
set_(62, "GSE230424_joint_axes", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "MODERATE_SUPPORT", "MEDIUM",
     "GSE230424 P3 joint 4-axis scatter", "Single slide; case-like; web only",
     "[AUTHOR VOICE NEEDED: GSE230424 P3 joint axis scatter]",
     notes="single slide P3")
set_(63, "Sample_table", "Paper2_methods", "SUPPLEMENTARY_FIGURE", "TECHNICAL_QC", "LOW",
     "28-sample comprehensive metadata table figure", "",
     "[AUTHOR VOICE NEEDED: 28-sample summary table]")
set_(64, "Master_moran_matrix", "Paper2_overview", "SUPPLEMENTARY_FIGURE", "STRONG_SUPPORT", "LOW",
     "Master 28-sample × 11-axis Moran I heatmap", "",
     "[AUTHOR VOICE NEEDED: master 28-sample multiaxis Moran's I matrix]")

# --- v9 (S_F65-S_F72)
set_(65, "TLS_individual_genes", "Paper2_Pillar2_support", "SUPPLEMENTARY_FIGURE", "MODERATE_SUPPORT", "MEDIUM",
     "TLS individual gene maps on PTC+HT P3 (CXCL13/CCL19/CCL21/CXCR5/CCR7/SELL)", "",
     "[AUTHOR VOICE NEEDED: TLS individual chemokine receptor maps]",
     notes="single slide; n=1 case-like; flag")
set_(66, "Sample_similarity", "Paper2_overview", "WEB_ONLY_ARCHIVE", "MODERATE_SUPPORT", "LOW",
     "28-sample similarity matrix in axis-Moran space", "Subsumed by S_F67 PCA + S_F64 master matrix",
     "[AUTHOR VOICE NEEDED: 28-sample similarity matrix]")
set_(67, "Sample_PCA", "Paper2_overview", "SUPPLEMENTARY_FIGURE", "MODERATE_SUPPORT", "LOW",
     "Sample PCA in axis-Moran space", "",
     "[AUTHOR VOICE NEEDED: sample PCA in axis-Moran space]")
set_(68, "Per_sample_trajectory", "Paper2_overview", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "Per-sample axis trajectory line plot", "",
     "[AUTHOR VOICE NEEDED: per-sample axis trajectory]")
set_(69, "Stage_3axis_spatial", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "Stage representative DM1/Epi/Prol 3-axis spatial", "",
     "[AUTHOR VOICE NEEDED: stage 3-axis spatial]")
set_(70, "Niche_compactness", "Paper2_methods", "WEB_ONLY_ARCHIVE", "TECHNICAL_QC", "LOW",
     "TROP2 niche compactness (size vs bbox)", "",
     "[AUTHOR VOICE NEEDED: niche compactness]")
set_(71, "Within_slide_rho", "Paper2_overview", "SUPPLEMENTARY_FIGURE", "STRONG_SUPPORT", "LOW",
     "Within-slide cross-axis Spearman heatmap", "",
     "[AUTHOR VOICE NEEDED: within-slide cross-axis Spearman summary]")
set_(72, "Ultimate_summary", "Paper2_overview", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "Ultimate summary infographic; superseded by S_F100", "",
     "[AUTHOR VOICE NEEDED: ultimate summary v9]")

# --- v10 (S_F73-S_F80) non-TROP2 deep biology
set_(73, "Thyroid_TF", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "MEDIUM",
     "Thyroid TF (PAX8/NKX2-1/FOXE1) maps", "",
     "[AUTHOR VOICE NEEDED: thyroid TF maps]",
     notes="signature-level, related to differentiation collapse")
set_(74, "Cell_cycle", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "Cell-cycle marker maps", "",
     "[AUTHOR VOICE NEEDED: cell cycle maps]")
set_(75, "Tcell_axis", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "T-cell marker maps", "",
     "[AUTHOR VOICE NEEDED: T-cell maps]")
set_(76, "Macrophage_TAM", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "Macrophage/TAM maps", "",
     "[AUTHOR VOICE NEEDED: macrophage TAM maps]")
set_(77, "Stromal_endo", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "Stromal/endothelial maps", "",
     "[AUTHOR VOICE NEEDED: stromal endothelial maps]")
set_(78, "Immune_escape", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "MEDIUM",
     "Immune escape (PD-L1/IDO1/B2M/HLA-A) — claim risk", "",
     "[AUTHOR VOICE NEEDED: immune escape maps — keep as exploratory only]",
     notes="HIGH risk for biomarker drift; gate as supplementary at most")
set_(79, "Glycolysis_OxPhos", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "LOW",
     "Glycolysis vs OxPhos maps", "",
     "[AUTHOR VOICE NEEDED: metabolism maps]")
set_(80, "Deep_geneset_morans", "Paper2_overview", "WEB_ONLY_ARCHIVE", "MODERATE_SUPPORT", "LOW",
     "16-sample × 8-deep-geneset Moran heatmap", "Subsumed by S_F64 master + S_F97 stage gradient",
     "[AUTHOR VOICE NEEDED: deep gene-set Moran summary]")

# --- v11 (S_F81-S_F88) cancer pathways
for fid, pw in zip(range(81,89),
        ["WNT","NOTCH","Hippo_YAP","RAS_MAPK","PI3K_AKT","MYC","TGF_beta","NFkB"]):
    set_(fid, f"Pathway_{pw}", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "MEDIUM",
         f"{pw} pathway spatial maps; exploratory only", "",
         f"[AUTHOR VOICE NEEDED: {pw} pathway spatial]",
         notes="signature-level; do not draw mechanistic claim")

# --- v12 (S_F89-S_F100)
set_(89, "Hypoxia_v12", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "MEDIUM",
     "Hypoxia signature spatial maps (v12 update)", "",
     "[AUTHOR VOICE NEEDED: hypoxia maps v12]")
set_(90, "EMT_v12", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "MEDIUM",
     "EMT program spatial maps (v12 update)", "",
     "[AUTHOR VOICE NEEDED: EMT maps v12]")
set_(91, "Thyroid_diff_collapse", "Paper2_Pillar2_support", "SUPPLEMENTARY_FIGURE", "EXPLORATORY", "HIGH",
     "v12 finding: differentiation spatial coherence collapse PT→ATC (n=4 ATC risk)", "",
     "[AUTHOR VOICE NEEDED: thyroid differentiation spatial coherence — exploratory; n=4 ATC]",
     notes="exploratory; do not center main claim on n=4 ATC")
set_(92, "Stemness_lineage", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "HIGH",
     "Stemness/lineage plasticity (case-like)", "",
     "[AUTHOR VOICE NEEDED: stemness maps]",
     notes="HIGH claim risk")
set_(93, "Senescence_SASP", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "MEDIUM",
     "Senescence/SASP maps", "",
     "[AUTHOR VOICE NEEDED: senescence SASP maps]")
set_(94, "Apoptosis_DDR", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "MEDIUM",
     "Apoptosis + DDR maps", "",
     "[AUTHOR VOICE NEEDED: apoptosis DDR maps]")
set_(95, "Angiogenesis", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "MEDIUM",
     "Angiogenesis maps", "",
     "[AUTHOR VOICE NEEDED: angiogenesis maps]")
set_(96, "Cellular_stress", "Paper2_Pillar2_support", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "MEDIUM",
     "Cellular stress quartet (ferroptosis/autophagy/UPR/lipid)", "",
     "[AUTHOR VOICE NEEDED: cellular stress quartet]")
set_(97, "Stage_gradient_24axes", "Paper2_overview", "SUPPLEMENTARY_FIGURE", "EXPLORATORY", "MEDIUM",
     "v12 finding: 24-axis stage gradient lineplot", "",
     "[AUTHOR VOICE NEEDED: 24-axis stage Moran gradient — exploratory]",
     notes="signature-level; flag as exploratory in caption")
set_(98, "ATC2_outlier_deep_dive", "Paper2_Pillar2_support", "SUPPLEMENTARY_FIGURE", "CASE_LIKE_OUTLIER", "HIGH",
     "ATC-2 single-sample 8-axis deep dive (n=1 outlier)", "",
     "[AUTHOR VOICE NEEDED: ATC-2 case-like outlier — strictly exploratory]",
     notes="n=1; HIGH claim risk; never center main claim")
set_(99, "Cross_axis_jaccard", "Paper2_overview", "SUPPLEMENTARY_FIGURE", "MODERATE_SUPPORT", "LOW",
     "Cross-axis spatial co-occurrence Jaccard heatmap", "",
     "[AUTHOR VOICE NEEDED: cross-axis spatial co-occurrence Jaccard]")
set_(100, "Centennial_summary", "Paper2_overview", "WEB_ONLY_ARCHIVE", "MODERATE_SUPPORT", "LOW",
     "Centennial findings infographic (web only)", "",
     "[AUTHOR VOICE NEEDED: centennial infographic]")

# --- v13 GeoMx DSP cross-platform external validation (S_F101-S_F106) ---
set_(101, "GeoMx_ROI_heatmap", "Paper2_external_validation", "SUPPLEMENTARY_FIGURE", "MODERATE_SUPPORT", "MEDIUM",
     "78 ROIs × 14 axes z-scored heatmap (cross-platform external validation)", "",
     "[AUTHOR VOICE NEEDED: GSE301163 GeoMx 78 ROI × 14 axes z-score heatmap (cross-platform external)]",
     notes="ROI-level pooled, NOT spot-level Visium; n=2 patients caveat")
set_(102, "GeoMx_TROP2_compartment", "Paper2_external_validation", "WEB_ONLY_ARCHIVE", "EXPLORATORY", "MEDIUM",
     "TROP2 PanCK+ vs VIM+ compartment by histology (honest negative)", "ROI-mean ns; web archive (does NOT refute spot-level Visium niche)",
     "[AUTHOR VOICE NEEDED: GSE301163 TROP2 compartment box (ROI-level ns, separate measurement layer)]",
     notes="ROI-mean cannot detect spatial niche; expected ns; flag honestly")
set_(103, "GeoMx_DM1_histology", "Paper2_external_validation", "SUPPLEMENTARY_FIGURE", "STRONG_SUPPORT", "MEDIUM",
     "DM1 axis by histology — PDTC vs microPTC PanCK+ p=0.001 cross-platform consistency", "",
     "[AUTHOR VOICE NEEDED: GSE301163 DM1 axis collapse PDTC vs microPTC (cross-platform consistent with Visium v12)]",
     notes="STRONG cross-platform validation of v12 differentiation collapse finding")
set_(104, "GeoMx_8gene_DM1", "Paper2_external_validation", "SUPPLEMENTARY_FIGURE", "MODERATE_SUPPORT", "MEDIUM",
     "Paper 1 8-gene mini-index by histology + genotype (cross-platform reactivity)", "",
     "[AUTHOR VOICE NEEDED: GSE301163 8-gene mini-index by histology and genotype]",
     notes="microPTC vs Normal PanCK+ p=0.017; cross-paper sanity, do NOT mix Paper 1/2 framing")
set_(105, "GeoMx_differentiation_3axis", "Paper2_external_validation", "SUPPLEMENTARY_FIGURE", "STRONG_SUPPORT", "MEDIUM",
     "3-axis (Thyroid_TF / RAI_8 / DM1_axis) by histology — TF preserved, function collapsed (Landa 2016 consistent)", "",
     "[AUTHOR VOICE NEEDED: GSE301163 3-axis differentiation collapse PDTC vs microPTC PanCK+]",
     notes="STRONG; Landa 2016 (Disc 3.1 cite) literature consistent")
set_(106, "Crossplatform_consistency", "Paper2_overview", "MAIN_FIGURE", "STRONG_SUPPORT", "MEDIUM",
     "★ Cross-platform summary: Visium 28 spot-level Moran ↔ GeoMx 78 ROI mean", "",
     "[AUTHOR VOICE NEEDED: 4-panel cross-platform consistency — Visium Moran's I trajectory vs GeoMx ROI-mean trajectory]",
     notes="MAIN candidate (M9?); v13 reach venue gate; n=2 GeoMx patients caveat")

# Build TSV
with open(OUT, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["figure_id","file_name","version_source","biological_axis","cohort",
                "paper_role","recommended_tier","claim_strength","risk_level",
                "include_reason","exclude_or_archive_reason",
                "manuscript_usage_sentence_draft","needs_recalculation","notes"])
    asset_dir = ROOT/"project/papers_hub_2026_05_04/assets/spatial_full"
    files = sorted([p.name for p in asset_dir.glob("S_F*.png")])
    name_by_id = {}
    for fn in files:
        try: nid = int(fn.split("_")[1].lstrip("F"))
        except: continue
        name_by_id[nid] = fn
    for nid in range(1, 107):
        fn = name_by_id.get(nid, "")
        if nid not in M:
            continue
        ax, role, tier, claim, risk, inc, exc, sent, recalc, notes = M[nid]
        w.writerow([f"S_F{nid}", fn, VER[nid], ax, default_cohort(nid),
                    role, tier, claim, risk, inc, exc, sent, recalc, notes])

print(f"WROTE {OUT}")

# Summary counts
counts = {"MAIN_FIGURE":0,"SUPPLEMENTARY_FIGURE":0,"WEB_ONLY_ARCHIVE":0,"DROP_OR_REMAKE":0}
for nid in M:
    counts[M[nid][2]] += 1
print("Tier counts:", counts)
print("Total in M:", len(M))

# (placeholder — see next overlay edit for proper 101-106 entries)
