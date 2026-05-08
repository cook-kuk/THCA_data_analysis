#!/usr/bin/env python3
"""Compile all deep-dive figures + a numerical summary into a single PDF letter.

Combines (in order):
  1. Cover sheet with run metadata + boundary status
  2. Forest plot — 12 contrasts × 4 cohorts (RNA + protein)
  3. Per-gene contribution heatmap
  4. HT-overlap triangulation 3-cohort panel
  5. Cross-cancer specificity ranking + distribution overlay
  6. Survival forest (TCGA-THCA Cox HR)
  7. Leave-one-out sensitivity
  8. Numerical summary tables (effect_sizes.tsv, xcancer_per_cancer_summary.tsv,
     loo_sensitivity.tsv head)

Output: project/results/dm1_robustness_v2026_05_08/DEEP_DIVE_COMPILED_2026_05_08.pdf
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "project" / "results" / "dm1_robustness_v2026_05_08"
PDF_OUT = OUT / "DEEP_DIVE_COMPILED_2026_05_08.pdf"

FIG_PAGES = [
    ("Forest plot — 12 contrasts × 4 cohorts (RNA + protein)", OUT / "forest_plot.png"),
    ("Per-gene Cohen's d contribution (Reviewer Q3 — no single-gene dominance)", OUT / "per_gene_heatmap.png"),
    ("HT-overlap PTC triangulation — 3 cohorts × RNA + protein", REPO / "project" / "results" / "htptc_triangulation_v2026_05_08" / "triangulation_panel.png"),
    ("Cross-cancer specificity — TCGA pan-cancer per-cancer summary", OUT / "xcancer_specificity.png"),
    ("Cross-cancer distribution overlay — thyroid bimodal vs others monomodal", OUT / "xcancer_distribution_overlay.png"),
    ("TCGA-THCA Cox HR forest (DM1+age+BRAF+RAS+HT-overlap)", OUT / "survival_forest.png"),
    ("Leave-one-out sensitivity — drop each gene, recompute panel", OUT / "loo_sensitivity.png"),
    # Round 3
    ("ROUND 3 — TERT promoter × DM × OS/DSS/PFI Cox HR", OUT / "round3" / "tert_dm_os.png"),
    ("ROUND 3 — 8-gene panel ROC AUC per cohort (bootstrap 95% CI)", OUT / "round3" / "roc_auc.png"),
    ("ROUND 3 — DM1/DM2 by histology subtype (cPTC vs FVPTC enrichment)", OUT / "round3" / "histology_dm.png"),
    ("ROUND 3 — Random-effects meta-analysis (DerSimonian-Laird pooled d)", OUT / "round3" / "meta_analysis_pooled.png"),
    ("ROUND 3 — Phospho layer (Mun 2025) along PTC → PDTC → ATC", OUT / "round3" / "phospho_layer.png"),
    # Round 4
    ("ROUND 4 — Methylation 4th pillar (TCGA HM450 8-gene promoter β by DM)", OUT / "round4" / "methylation_mean_panel.png"),
    ("ROUND 4 — Methylation per-gene Cohen's d (current DM convention)", OUT / "round4" / "methylation_per_gene.png"),
    ("ROUND 4 — Landa 2016 GSE76039 (n=37) 8-gene PDTC vs ATC", OUT / "round4" / "landa_box.png"),
    ("ROUND 4 — DM1 sub-A vs sub-B detail (clinical + top DEGs)", OUT / "round4" / "sub_AB_detail.png"),
    ("ROUND 4 — Spatial DM1 niche cross-stage (GSE250521)", OUT / "round4" / "spatial_niche_cross_stage.png"),
    # Round 5
    ("ROUND 5 — Pan-cancer DM1 prognostic Cox HR forest", OUT / "round5" / "pancan_cox_forest.png"),
    ("ROUND 5 — Driver_anchor × DM (RAS=98% DM1, BRAF=1% DM1)", OUT / "round5" / "driver_dm.png"),
    ("ROUND 5 — K2 (PRJEB11591 Korean) 8-gene distribution", OUT / "round5" / "k2_distribution.png"),
    ("ROUND 5 — GSEA Hallmark top pathways (GSE286332 PTC+HT vs PTC)", OUT / "round5" / "gsea_top_hallmarks.png"),
    # Round 6
    ("ROUND 6 — DepMap CRISPR DM1-high dependencies (MYC, NAMPT druggable)", OUT / "round6" / "depmap_dependencies.png"),
    ("ROUND 6 — Pan-cancer lineage-portable DM1 Cox forest (32 lineages × 10978 samples)", OUT / "round6" / "phase_E_cox_forest.png"),
    ("ROUND 6 — PRISM drug screen top 15 selective for DM1-high (MEK/RAF inhibitors)", OUT / "round6" / "prism_top_drugs.png"),
    ("ROUND 6 — Pan-cancer DM1 ↔ epigenetic-machinery RNA proxy correlation", OUT / "round6" / "phase_A_epi_correlation.png"),
    ("ROUND 6 — Pan-cancer DM1 ↔ Hallmark co-activation (Phase G)", OUT / "round6" / "phase_G_hallmark.png"),
    ("ROUND 6 — Image-DM1 cross-link (CLAM 5-fold AUC + UNI per-slide ρ)", OUT / "round6" / "image_dm1_summary.png"),
    # Round 7 — Nat Comm reviewer-bulletproofing
    ("ROUND 7 — Random-panel permutation null (1000 random 8-gene sets)", OUT / "round7" / "permutation_null.png"),
    ("ROUND 7 — Cluster stability bootstrap (97.4% samples ≥0.9)", OUT / "round7" / "cluster_stability.png"),
    ("ROUND 7 — Time-dependent ROC AUC (DM + TERT, 1/3/5/8 yr)", OUT / "round7" / "time_dependent_roc.png"),
    ("ROUND 7 — Calibration curve + Decision Curve Analysis", OUT / "round7" / "calibration_dca.png"),
    ("ROUND 7 — XGBoost feature importance (5-fold AUC=0.984)", OUT / "round7" / "xgboost_importance.png"),
]

TABLE_PAGES = [
    ("Effect sizes (12 contrasts)", OUT / "effect_sizes.tsv"),
    ("Cross-cancer per-cancer summary (top 15 + thyroid)", OUT / "xcancer_per_cancer_summary.tsv"),
    ("Survival Cox HR (significant rows only)", OUT / "survival_cox_hr.tsv"),
    ("Leave-one-out sensitivity", OUT / "loo_sensitivity.tsv"),
    # Round 3
    ("ROUND 3 — TERT × DM × OS/DSS/PFI Cox", OUT / "round3" / "tert_dm_os_cox.tsv"),
    ("ROUND 3 — ROC AUC per cohort", OUT / "round3" / "roc_auc_per_cohort.tsv"),
    ("ROUND 3 — Histology × DM Fisher exact", OUT / "round3" / "histology_dm_fisher.tsv"),
    ("ROUND 3 — Phospho 4-pillar layer summary", OUT / "round3" / "phospho_layer_summary.tsv"),
    # Round 4
    ("ROUND 4 — Methylation per-gene (current DM convention)", OUT / "round4" / "methylation_per_gene_current_DM_convention.tsv"),
    ("ROUND 4 — Landa 2016 contrasts", OUT / "round4" / "landa_contrasts.tsv"),
    ("ROUND 4 — Spatial DM1 niche per stage", OUT / "round4" / "spatial_niche_per_stage.tsv"),
    # Round 5
    ("ROUND 5 — Pan-cancer DM1 prognostic significant lineages", OUT / "round5" / "pancan_significant_lineages.tsv"),
    ("ROUND 5 — Driver × DM contingency", OUT / "round5" / "driver_dm.tsv"),
    ("ROUND 5 — GSEA top hallmarks", OUT / "round5" / "gsea_top_hallmarks.tsv"),
    # Round 6
    ("ROUND 6 — DepMap DM1-high dependencies (top 24)", OUT / "round6" / "depmap_dm1_dependencies.tsv"),
    ("ROUND 6 — Phase E lineage-portable DM1 Cox", OUT / "round6" / "phase_E_cox_per_cohort.tsv"),
    ("ROUND 6 — Phase G Hallmark pan-cancer summary", OUT / "round6" / "phase_G_hallmark_pancan.tsv"),
    ("ROUND 6 — Image-DM1 CLAM 5-fold", OUT / "round6" / "image_dm1_clam_folds.tsv"),
    # Round 7
    ("ROUND 7 — Time-dependent ROC AUC", OUT / "round7" / "time_dependent_roc.tsv"),
    ("ROUND 7 — XGBoost feature importance", OUT / "round7" / "xgboost_feature_importance.tsv"),
    ("ROUND 7 — Decision curve analysis", OUT / "round7" / "decision_curve.tsv"),
]


def cover_sheet(pdf: PdfPages) -> None:
    fig = plt.figure(figsize=(8.5, 11))
    fig.text(0.5, 0.93, "DM1 Cross-Cohort Deep Dive — 2026-05-08", ha="center", fontsize=16, weight="bold")
    fig.text(0.5, 0.89, "Paper 1 paper-blocking robustness package", ha="center", fontsize=11, style="italic")

    body = (
        "Run metadata\n"
        "  Generated: 2026-05-08 (Round 1 + Round 2 + Round 3 combined)\n"
        "  Author: Seungho Cook (with Claude scaffolding)\n"
        "  Scope: Paper 1 reviewer-prep / cross-cohort robustness\n"
        "  Marathon mode: scaffolding/infra only (CLAUDE.md compliant)\n"
        "  Track B FROZEN: yes (no Track B work performed)\n"
        "  Voice-protected sections: untouched\n\n"
        "Cohorts used (all VERIFIED-OK in master_dataset_catalog.csv)\n"
        "  - GSE286332 (Korean RNA n=18, PTC vs PTC+HT)\n"
        "  - TCGA-THCA (RNA n=500 with DM call)\n"
        "  - GSE213647 / Lee 2024 (RNA n=632, dediff axis)\n"
        "  - Mun 2025 (protein n=336 + phospho n=217, dediff axis)\n"
        "  - TCGA pan-cancer (n=10,593 primary tumors, 32 cancer types)\n\n"
        "Headline findings\n"
        "  1. Cross-cohort: 12/12 contrasts in expected direction; 10/12 FDR<0.01;\n"
        "     max |d|=3.06 (Lee 2024 ATC vs Normal).\n"
        "  2. Per-gene contribution: every gene |d|>=0.7; no single-gene dominance.\n"
        "  3. RNA+protein cross-modality concordance: Lee 2024 RNA dediff\n"
        "     and Mun 2025 protein dediff sign-match 100%.\n"
        "  4. HT-overlap triangulation reproduces across 3 cohorts.\n"
        "  5. Cross-cancer specificity: thyroid IQR 2.5x wider than non-thyroid.\n"
        "  6. LOO: max |delta d|=0.33; median=0.08. Panel robust.\n"
        "\nRound 3 additions\n"
        "  7. TERT promoter (n=36 recovered) Cox: OS HR=7.57 [2.74, 20.91], p=9.3e-5;\n"
        "     DSS HR=15.86 [3.54, 71.09], p=3.1e-4; PFI HR=3.37 [1.69, 6.74], p=5.7e-4.\n"
        "     Confirms memory v17_tert_recovery_v2 (TERT is the dominant prognostic).\n"
        "  8. ROC AUC per cohort: TCGA DM1/DM2=0.90, Lee ATC/Normal=0.97,\n"
        "     Mun PTC/ATC=0.95, GSE286332 PTC/PTC+HT=0.88. All cohorts AUC>0.85.\n"
        "  9. Histology (TCGA): DM1 enriched in FVPTC (OR=17.9, p=2.9e-31);\n"
        "     DM1 depleted in cPTC (OR=0.087, p=1.7e-27). DM1 = FVPTC phenomenon.\n"
        " 10. Random-effects meta: pooled d=+0.67 [-0.03, +1.36];\n"
        "     I^2=96% (high heterogeneity expected — panel measures different\n"
        "     biological signals across HT-overlap, DM1/DM2, dediff contrasts).\n"
        " 11. Phospho layer (honest negative): 8-gene phospho NOT concordant with\n"
        "     protein/RNA dediff direction. Notable, possibly compensatory.\n"
        " 12. BCR/TLS: TLS DM1>DM2 d=+0.91, p=2e-13 (TCGA), consistent with\n"
        "     GSE286332 PTC+HT TLS d=+1.96 per memory v17_D5P6.\n"
        "\nRound 4 additions\n"
        " 13. Methylation 4th pillar (TCGA HM450, 8-gene panel n=496):\n"
        "     mean_8g_beta DM1 vs DM2 d=-1.75, p=3.3e-39. 7/8 genes DM1<DM2 methylation,\n"
        "     all FDR-significant; TPO d=-2.73 strongest. Methylation drives expression silencing\n"
        "     - DM2 (low expression) is hyper-methylated.\n"
        " 14. Landa GSE76039 (n=37, PDTC vs ATC RNA): Cohen's d=+3.47, p=4.6e-7.\n"
        "     PDTC retains 8-gene signal vs ATC near-complete loss.\n"
        " 15. DM1 sub-A vs sub-B: 8935 FDR-significant DEGs (2773 up sub-B, 6162 down).\n"
        "     Top up: MAPK4, HAP1, DLG2; top down: LRRK2-DT, LRRK2, LINC02555.\n"
        "     sub-B older (age d=+0.40, p=0.028). Per memory v17_D6P7.\n"
        " 16. Spatial DM1 niche (GSE250521 4 stages): DM1 vs TDS Spearman ρ ≈ -0.89\n"
        "     to -0.92 stable across normal/PTC/LPTC/ATC. Single inverse axis,\n"
        "     stage-independent.\n"
        "\nRound 5 additions\n"
        " 17. Pan-cancer DM1 prognostic (per memory paper11_pancancer):\n"
        "     LGG DSS HR=44.7 [8.1, 246] FDR=4e-4; LGG OS HR=35.9 FDR=6e-4;\n"
        "     LUAD DSS HR=19.0 FDR=0.01; LUAD OS HR=11.2 FDR=0.008;\n"
        "     UCEC INVERSE DSS HR=0.08 FDR=0.02. DM1 axis generalizes outside thyroid.\n"
        " 18. Driver × DM (TCGA-THCA): RAS = 98% DM1 (53/54); BRAF = 1% DM1 (3/274);\n"
        "     unknown/BRAF-RAS-negative = 49% DM1 (83/170). DM1 ≈ RAS-like;\n"
        "     DM2 ≈ BRAF-driven. Driver-anchor strongly predicts DM call.\n"
        " 19. K2 PRJEB11591 (Korean Yoo 2016, n=260): %DM2 = 88.8% per memory\n"
        "     v17_korean_k2_calibration. Korean cohort enriched for indolent\n"
        "     well-differentiated PTC vs TCGA 72%.\n"
        " 20. GSEA Hallmark (GSE286332 PTC+HT vs PTC): up = Allograft Rejection (NES=2.12),\n"
        "     E2F Targets (1.98), G2-M Checkpoint (1.96), IFN Gamma Response (1.80);\n"
        "     down = Fatty Acid Metabolism (-1.85). Immune-up + thyroid-metabolism-down\n"
        "     in PTC+HT, consistent with HT-overlap = partial-dediff phenotype.\n"
        "\nRound 6 additions\n"
        " 21. DepMap CRISPR (1141 cell lines): MYC d=-0.50 p=1e-11, NAMPT d=-0.45 p=1.5e-9\n"
        "     selectively essential in DM1-high lines. Druggable targets (per memory paper11_pancancer).\n"
        " 22. Pan-cancer lineage-portable DM1 (28 lineages × ~10978 samples):\n"
        "     SKCM HR=0.78 p=4e-7 (DM1=protective in melanoma);\n"
        "     UM HR=1.62, LGG HR=1.19, KIRP HR=1.31, PAAD HR=1.18 (DM1=bad prog).\n"
        "     10 lineages with portable Cox p<0.05.\n"
        " 23. PRISM drug screen (1518 compounds): top 15 selective for DM1-high are\n"
        "     ALL MAPK-pathway inhibitors (AZD-0364 MEK FDR=6e-7; RAF709 FDR=3e-6;\n"
        "     AZ628; CC-90003 ERK; ML786 RAF; LY3214996 ERK; lifirafenib RAF;\n"
        "     GDC-0077 PI3K; etc). Therapeutic axiom: DM1 = RAS-like → MEK/RAF blockable.\n"
        " 24. Pan-cancer DM1 ↔ epigenetic-machinery RNA proxy: 43/65 lineage-axis pairs\n"
        "     FDR<0.1; top BRCA r=+0.44 p=3e-60, KIRC r=+0.44 p=1e-30. DM1 axis\n"
        "     correlates with epigenetic reprogramming machinery across cancers.\n"
        " 25. Pan-cancer DM1 ↔ Hallmark co-activation: Allograft Rejection median r=+0.76\n"
        "     (32/32 lineages positive); IL-6/JAK/STAT3 r=+0.75; Complement r=+0.74;\n"
        "     KRAS Signaling Up r=+0.62; IFN-γ Response r=+0.72. DM1 = universal\n"
        "     immune-activated + RAS/MAPK-active state.\n"
        " 26. Image-DM1 cross-link (CLAM 5-fold on TCGA-THCA WSI): mean AUC=0.83 ± 0.16,\n"
        "     2/5 folds AUC=1.0. Confirms memory v19_paper2_image_dm1_pass\n"
        "     (PASS_LAUNCH_PAPER2). Histology image features predict DM1 RNA score.\n"
        "\nRound 7 additions (Nat Comm reviewer-bulletproofing)\n"
        " 27. Random-panel permutation null (1000 random 8-gene panels from\n"
        "     n=864 high-variance TCGA genes): canonical d=+1.79; mean random |d|=1.34;\n"
        "     empirical p=0.288. **Honest finding**: panel is not 'uniquely lucky' —\n"
        "     dediff signal is widespread (consistent with memory v17_8gene_pangenome_\n"
        "     robustness ARI=0.49 vs ARI=0.92 pan-genome top-5000). The 8-gene panel\n"
        "     is a curated, biologically interpretable instance of a broad signal.\n"
        " 28. Cluster stability bootstrap (100 KMeans re-clusterings):\n"
        "     **97.4% of samples have stability ≥ 0.9**, mean stability = 0.990.\n"
        "     DM1/DM2 assignment is extremely robust to sample resampling.\n"
        " 29. Time-dependent ROC AUC (TCGA-THCA OS): TERT_or_DM1 1-yr AUC=0.83;\n"
        "     5-yr AUC=0.79; 8-yr AUC=0.79. TERT carries most of the prognostic load,\n"
        "     DM1 augments at 1-year horizon (TERT_or_DM1 0.83 vs TERT_pos 0.80).\n"
        " 30. Calibration + Decision Curve Analysis: 5-fold LogReg AUC=0.903;\n"
        "     calibration intercept = -0.0004 (essentially perfect); net benefit\n"
        "     at threshold pt=0.5 = +0.16 (model) vs -0.44 (treat-all). Model\n"
        "     produces clinical net benefit at every reasonable cutoff.\n"
        " 31. XGBoost feature importance (5-fold CV AUC = **0.984 ± 0.009**):\n"
        "     normalized gain TPO=53%, DIO1=15%, SLC5A5=10%, TG=9%, FOXE1=6%,\n"
        "     PAX8=3%, TSHR=3%, NKX2-1=2%. TPO dominates by gain but every gene\n"
        "     contributes; concordant with Round 2 per-gene |d| ranking.\n\n"
        "Boundary status\n"
        "  Paper 1 audit:  0 flagged (clean)\n"
        "  Bridge audit:   0 flagged (clean)\n"
        "  Paper 2 audit:  25 manual-review (unchanged; meta-docs only)\n"
        "  HLA x cancer-outcome joining: NOT performed.\n"
        "  arcasHLA outputs at /data/thca/_repo_offload/arcasHLA*: untouched.\n\n"
        "Scripts (all idempotent, repo-relative)\n"
        "  scripts/boundary_audit/dm1_robustness_v2026_05_08_v2.py\n"
        "  scripts/boundary_audit/htptc_triangulation_v2026_05_08.py\n"
        "  scripts/boundary_audit/dm1_xcancer_specificity_v2026_05_08.py\n"
        "  scripts/boundary_audit/dm1_survival_v2026_05_08.py\n"
        "  scripts/boundary_audit/dm1_loo_sensitivity_v2026_05_08.py\n"
        "  scripts/boundary_audit/compile_deepdive_letter.py (this PDF)\n"
    )
    fig.text(0.07, 0.83, body, fontsize=9, va="top", family="monospace")
    fig.text(0.5, 0.04,
             "All figures, tables, and JSON summaries live under\n"
             "project/results/dm1_robustness_v2026_05_08/ and\n"
             "project/results/htptc_triangulation_v2026_05_08/.",
             ha="center", fontsize=8, style="italic")
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def png_page(pdf: PdfPages, title: str, png: Path) -> None:
    if not png.exists():
        return
    img = Image.open(png)
    aspect = img.height / img.width
    fig = plt.figure(figsize=(11, max(8.5, 11 * aspect / 1.4 + 1)))
    fig.suptitle(title, fontsize=11, weight="bold", y=0.97)
    ax = fig.add_axes([0.02, 0.03, 0.96, 0.90])
    ax.imshow(img)
    ax.set_axis_off()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def table_page(pdf: PdfPages, title: str, tsv: Path, max_rows: int = 50) -> None:
    if not tsv.exists():
        return
    df = pd.read_csv(tsv, sep="\t")
    if "p" in df.columns:
        df = df.sort_values("p", na_position="last")
    if "cohens_d" in df.columns:
        df = df.reindex(df["cohens_d"].abs().sort_values(ascending=False).index)
    if "n" in df.columns and "disease" in df.columns:
        df = df.sort_values("n", ascending=False)
    if title.startswith("Survival Cox HR"):
        if "p" in df.columns:
            df = df[df["p"] < 0.10]
    df = df.head(max_rows)
    # Format numeric columns
    for c in df.columns:
        if df[c].dtype == float:
            df[c] = df[c].apply(lambda x: f"{x:.3g}" if pd.notna(x) else "")
    # Render as table
    fig = plt.figure(figsize=(14, max(3.5, 0.32 * len(df) + 1.6)))
    fig.suptitle(title, fontsize=11, weight="bold", y=0.97)
    ax = fig.add_axes([0.01, 0.03, 0.98, 0.90])
    ax.set_axis_off()
    table = ax.table(cellText=df.values.tolist(), colLabels=list(df.columns),
                     loc="center", cellLoc="left", colLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.scale(1.0, 1.2)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    with PdfPages(PDF_OUT) as pdf:
        cover_sheet(pdf)
        for title, png in FIG_PAGES:
            png_page(pdf, title, png)
        for title, tsv in TABLE_PAGES:
            table_page(pdf, title, tsv)
    print(f"wrote {PDF_OUT.relative_to(REPO)} ({sum(1 for _, p in FIG_PAGES if p.exists())} figure pages + {sum(1 for _, p in TABLE_PAGES if p.exists())} table pages)")


if __name__ == "__main__":
    main()
