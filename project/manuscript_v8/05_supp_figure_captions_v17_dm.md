---
title: "Paper 1 — Supplementary figure captions (v17 audit phase + Dark matter phase 1/2 + Landa 2016 evidence)"
date: 2026-05-04
status: CAPTION DRAFT — filename + project memory context 기반 추정. 본인 실제 figure 직접 보고 panel content / stats 정확 표현 검증 필요.
target_format: Cell Press style (bold title + (A)/(B)/.. panel descriptions + n + test + p; 본인 voice 영역 아님)
parent_caption_file: project/manuscript_v8/05_figure_captions.md (Fig 1-8 main + S1-S9 main supp)
verification_markers: ⚠ = 본인 figure 직접 보고 검증 필요 / ✅ = filename + memory 에서 명확
---

# Paper 1 — Supplementary figure captions (v17 audit + Dark matter phases + Landa 2016)

These supplementary figures support Paper 1 Section 2 (Results) main claims with sanity-check, robustness, multisite, and per-patient evidence. Source: `project/results/audit_2026_04_29/`, `project/results/dark_matter_phase1/`, `project/results/dark_matter_phase2/`, `project/results/landa2016_evidence/` (paths inferred from naming convention — verify).

---

## Part A — v17 Audit phase (2026-04-29) sanity-check figures

Memory cross-ref: `v17_audit_session_2026_04_29.md` — 17 angles A-J full sweep, 8-gene paper not blocked.

### Figure SA1. **4-way revalidation of the 8-gene DM1/DM2 cluster across driver mutation strata.**

`figure_4way_revalidation`

(A-D) Stratum-specific overall-survival outcome and DM cluster cross-tabulation across the 4-way decomposition of TCGA-THCA primary tumors by BRAF V600E status × TERT promoter mutation status. (A) BRAF+ / TERT+. (B) BRAF+ / TERT−. (C) BRAF− / TERT+ (n = 4, small-N caveat). (D) BRAF− / TERT− (the dark matter compartment, sub-stratified by DM1/DM2). DM1 versus DM2 sub-strata are colored as in Figure 1B. ⚠ Panel-specific n and Cox HR per stratum to be filled from `4way_revalidation` source TSV.

Statistical tests: stratified Cox proportional hazards; log-rank for KM curves.

### Figure SA2. **FFPE versus fresh-frozen tissue compatibility QC for the 8-gene panel.**

`figure_ffpe_qc`

(A) Distribution of 8-gene panel score by tissue processing type (FFPE vs FF) across all available cohorts; Kolmogorov-Smirnov p = 0.44 (no detectable distributional shift). (B) ⚠ Per-gene paired comparison (FFPE matched samples vs FF matched samples) — Bland-Altman or scatter to verify per-gene preservation. (C) ⚠ DM1/DM2 call concordance between FFPE and FF for samples with both modalities (if applicable). 

Statistical tests: Kolmogorov-Smirnov for distribution; per-gene Wilcoxon signed-rank for paired analysis.

### Figure SA3. **Panel size sensitivity — 8 / 10 / 12 / 16-gene variants converge on the same DM1/DM2 axis.**

`figure_panel_coverage`

(A) TCGA-THCA 5-fold cross-validated AUC for DM1/DM2 classification across panel sizes: 8-gene (0.962), 10-gene (+DIO2, IYD; 0.972), 12-gene (+SLC26A4, SLC5A8; 0.969), 16-gene full TDS-core (0.975). ΔAUC 8 vs 16 = 0.013, NS. (B) ⚠ Pairwise concordance matrix (8/10/12/16-gene cluster calls). (C) ⚠ Per-panel ARI vs reference (Memory: P4 ARI ladder).

Statistical tests: 5-fold cross-validated AUC; adjusted Rand index for cluster concordance.

### Figure SA4. **Single-cell wrap-up — thyrocyte-intrinsic 8-gene signal across external sc cohorts.**

`figure_sc_wrapup`

(A) Per-cohort summary of 8-gene signal in thyrocyte-restricted populations (KRT8+ KRT19+ EPCAM+) from GSE184362 (Pu et al., 2021; n = 6 PTC), GSE193581 (Lu 2023; n = 23 samples), and GSE241184 (n = 1 Phase 1). (B) Per-patient Spearman correlation between tumor and adjacent normal thyrocyte 8-gene scores in GSE184362 (r = 0.798 to 0.886, all p < 10⁻¹⁰). (C) ⚠ Author-independence cross-cohort concordance check (GSE184362 vs Lu 2023). 

Statistical tests: Spearman r per patient with Bonferroni-adjusted p; cross-cohort concordance via Pearson r.

### Figure SA5. **Differentiation trajectory along the 8-gene panel score (TCGA-THCA primary tumors).**

`figure_trajectory`

(A) Pseudotime / monotonic ordering of TCGA-THCA primary tumors along the 8-gene panel score axis. (B) ⚠ Driver mutation distribution (BRAF V600E, RAS hotspots, RET fusion, NTRK fusion) along the trajectory. (C) ⚠ Differentiation-related gene expression (e.g., NIS/SLC5A5, TPO, TG, TSHR, PAX8) along the trajectory. (D) DM1/DM2 cluster boundary along the score axis.

Statistical tests: Spearman r between trajectory rank and per-gene expression.

### Figure SA6. **MSK-IMPACT bias panel — advanced-disease enrichment audit.**

`msk_bias_panel`

(A) MSK-IMPACT thyroid cohort (Landa et al., 2016; n = 117) composition: PDTC (n = 84) vs ATC (n = 33), explicitly labeled as **advanced disease cohort** (not primary). (B) ⚠ Mutation frequency comparison MSK-IMPACT versus TCGA-THCA primary tumors (BRAF, RAS, TERT, kinase fusion partners). (C) ⚠ Stage / clinical aggressiveness distribution differential. (D) ⚠ Cross-cohort DM1 prevalence — MSK vs TCGA — with bias disclosure.

Statistical tests: Fisher exact for mutation frequency comparison; chi² for clinical strata.

---

## Part B — Dark matter phase 1 — single-cell UMAP foundations

Memory cross-ref: `v17_audit_session_2026_04_29.md` (P2-A finding: GSE184362 r > 0.79 in 6/6 patients, TRUE INDEPENDENT).

### Figure SB1. **Single-cell UMAP overview of GSE184362 (Pu et al., 2021).**

`dm_phase1_umap_overview`

(A) UMAP embedding of all sequenced cells from GSE184362 (n = 6 PTC patients, Fudan University; total cell count ⚠ verify), colored by patient identity. (B) UMAP colored by canonical cell-type annotation (thyrocyte / immune / stromal / endothelial). (C) ⚠ Tumor versus adjacent normal labeling per cell. Scale bars or cluster boundary marks as appropriate.

Method note: Seurat or scanpy-based clustering; per Pu et al., 2021 published preprocessing; thyrocyte annotation via KRT8 + KRT19 + EPCAM.

### Figure SB2. **Single-cell UMAP — molecular signatures projected onto cell embedding.**

`dm_phase1_umap_signatures`

(A) UMAP from SB1A colored by 8-gene panel score per cell (continuous gradient). (B) UMAP colored by HLA-II module signature (HLA-DRA, HLA-DRB1, HLA-DPB1, HLA-DQB1) — note: Paper 1 uses HLA-II only as residualization control (Methods); deeper HLA-II analysis is Paper 2 territory. (C) UMAP colored by canonical thyroid differentiation transcripts (TG, TPO, TSHR averaged). (D) ⚠ UMAP colored by proliferation signature (MKI67, TOP2A) for completeness.

Statistical notes: per-cell module score = z-mean of normalized expression; gradient color scale 0-1 or per-module z-bounds.

### Figure SB3. **Single-cell UMAP — thyrocyte-restricted DM1 signal.**

`dm_phase1_umap_thyrocytes`

(A) UMAP subset to thyrocyte-annotated cells (KRT8+ KRT19+ EPCAM+) from SB1B. (B) Thyrocyte UMAP colored by 8-gene panel score per cell. (C) Thyrocyte UMAP colored by tumor versus adjacent normal labeling (paired tissue per patient). (D) ⚠ Per-patient thyrocyte distribution along the 8-gene score axis (violin or ridge plot).

Statistical notes: Mann-Whitney U for tumor vs normal thyrocyte scores within each patient; per-patient Spearman r reported in Figure SC1 (per_patient_r_forest).

---

## Part C — Dark matter phase 2 — multisite + per-patient + external validation

Memory cross-ref: `v17_audit_session_2026_04_29.md` (P2-A multi-patient sc PASS); Pu et al., 2021 (GSE184362 6 patients with paired tumor + normal).

### Figure SC1. **Per-patient Spearman r forest — tumor vs adjacent normal thyrocyte 8-gene scores in GSE184362.**

`dm_phase2_per_patient_r_forest`

Forest plot of per-patient Spearman r between mean 8-gene panel score in tumor thyrocytes versus adjacent normal thyrocytes for each of the 6 GSE184362 PTC patients (Pu et al., 2021). Per-patient r ranges from 0.798 to 0.886, all p < 10⁻¹⁰; all 6 patients direction-consistent. ⚠ Per-patient n_cells (tumor / normal) annotated alongside each forest entry.

Statistical tests: Spearman r per patient with Bonferroni adjustment for n = 6 multiple comparisons.

### Figure SC2. **Multisite trajectory — DM1/DM2 score reproducibility across cohorts.**

`dm_phase2_multisite_trajectory`

Trajectory of 8-gene panel scores across the differentiation axis in TCGA-THCA primary tumors (n = 504), MSK-IMPACT advanced disease (Landa et al., 2016; n = 117), Korean K2 / PRJEB11591 (n = 260), and Korean Lee / GSE213647 (n = 632). (A) ⚠ Per-cohort distribution of 8-gene scores (violin or density). (B) DM1/DM2 cluster boundary preserved across cohorts. (C) ⚠ Cross-cohort score distribution correlation matrix.

Statistical tests: Kolmogorov-Smirnov for distribution comparison across cohorts; Pearson r for cross-cohort distribution alignment.

### Figure SC3. **Multisite UMAP — joint cell embedding across single-cell cohorts.**

`dm_phase2_umap_multisite`

(A) Joint UMAP of single-cell data from GSE184362 (Pu 2021; n = 6 patients), GSE193581 (Lu 2023; n = 23 samples), and GSE241184 (Phase 1; n = 1) after batch-corrected integration (e.g., Harmony or scVI; ⚠ verify integration method). (B) Joint UMAP colored by source dataset. (C) Joint UMAP colored by 8-gene panel score per cell — demonstrates dataset-independent gradient. (D) ⚠ Thyrocyte sub-cluster preservation across datasets.

Statistical notes: ARI between per-dataset clusters and joint clusters; batch effect quantification (e.g., kBET or LISI).

### Figure SC4. **GSE184362 (Pu 2021) summary table — 6-patient external validation overview.**

`dm_phase2_gse184362_summary`

Tabular figure summarizing GSE184362 per-patient metadata: patient ID, age, sex, BRAF V600E status, RAS status, fusion status (if known from published metadata), tumor stage, n_cells per condition (tumor / adjacent normal), per-patient Spearman r (from SC1), DM1 vs DM2 prediction. ⚠ Direct extraction from Pu 2021 metadata + this paper's per-patient pseudobulk classification.

Statistical notes: descriptive only; no inferential test in this summary panel.

### Figure SC5. **External validation pooled — joint scatter of tumor vs normal thyrocyte 8-gene scores.**

`dm_phase2_p2a_external_validation`

(A) Per-patient pseudobulk 8-gene panel score in tumor versus adjacent normal thyrocytes from GSE184362 (Pu 2021), one point per patient (n = 6). (B) ⚠ Same plot for Lu 2023 GSE193581 (if applicable cohort overlap). (C) Combined external pooled scatter with per-cohort markers and pooled regression line.

Statistical tests: Pearson r and slope for tumor-vs-normal pooled comparison.

### Figure SC6. **Pooled scatter — DM1/DM2 axis across external single-cell cohorts.**

`dm_phase2_p2a_pooled_scatter`

Pooled scatter of per-patient (or per-sample) 8-gene panel score versus DM1 probability across GSE184362 (Pu 2021) and Lu 2023 GSE193581. (A) Per-patient point colored by source cohort. (B) ⚠ Marginal density distribution per axis. (C) Linear regression with 95% CI.

Statistical tests: Pearson r with 95% CI; bootstrap regression.

### Figure SC7. **K2 (PRJEB11591) versus Yoo 2016 reference — 8-gene score concordance.**

`dm_phase2_k2_vs_yoo`

(A) Per-sample 8-gene panel score in K2 / PRJEB11591 (n = 260) computed via this paper's pipeline versus Yoo 2016's published TDS-core summary (if directly comparable; ⚠ verify whether Yoo 2016 reports per-sample TDS-core scores or only cohort-level statistics). (B) Cluster call concordance (DM1 vs DM2). (C) ⚠ K2 mini-index calibration mismatch diagnostic — within-sample-centered profile vs raw TPM form.

Statistical tests: Spearman r for score concordance; Cohen's kappa for cluster call concordance; correlation diagnostics for calibration mismatch (R4-4 audit).

---

## Part D — Landa 2016 evidence (Discussion §3.1 cite save)

Memory cross-ref: `v17_landa2016_cite_save.md` — Discussion §3.1 cite correction (Krishnamoorthy 2025 Nat Comm misattribution → Landa 2016 JCI 126(3):1052-1066). 5/8 panel gene overlap with Landa 2016 ATC-silenced gene list.

### Figure SD1. **Landa 2016 GSE76039 — differentiation transcript suppression in advanced thyroid cancer.**

`landa2016_gse76039_heatmap`

(A) Per-sample heatmap of canonical thyroid differentiation transcripts (TG, TSHR, TPO, PAX8, SLC26A4, DIO1, DUOX2 — the 7-gene list reported as profoundly suppressed in ATC by Landa et al., 2016) across the GSE76039 transcriptome subset of the Landa 2016 cohort (PDTC n = 17 + ATC n = 20; ⚠ verify exact n from Landa 2016). (B) ⚠ Side-by-side comparison with this paper's 8-gene panel (5 of 8 overlap: TG, TSHR, TPO, PAX8, DIO1). (C) ⚠ Convergence note: Yoo 2016 (this paper's panel origin) and Landa 2016 (advanced-disease ATC list) reach the same 5-gene differentiation core via independent paths. Convergent biological validation, not panel-driven circular inference (per Discussion §3.1 reverse-causality framing).

Statistical tests: per-sample z-scoring per gene; heatmap clustering by Ward linkage.

Cite: Landa I, Ibrahimpasic T, Boucai L, Sinha R, Knauf JA, Shah RH, Dogan S, Ricarte-Filho JC, **Krishnamoorthy GP** (9th coauthor), Xu B, Schultz N, Berger MF, Sander C, Taylor BS, Ghossein R, Ganly I, Fagin JA. *J Clin Invest.* 2016;126(3):1052-1066. doi:10.1172/JCI85271. PMID: 26878173. PMC4767360.

---

## Caption verification queue (본인 figure 직접 read 후 confirm)

| Figure ID | Verification needed |
|---|---|
| SA1 4way_revalidation | Per-stratum n, Cox HR table |
| SA2 ffpe_qc | Panels B + C exact analysis (paired Bland-Altman? per-gene Wilcoxon?) |
| SA3 panel_coverage | Panels B + C (concordance matrix + ARI ladder) |
| SA4 sc_wrapup | Panel C author-independence specifics |
| SA5 trajectory | Panels B + C + D (driver/diff gene/DM cluster overlay) |
| SA6 msk_bias_panel | Panels B-D (mutation freq comparison + stage + DM prevalence) |
| SB1 umap_overview | Total cell count, panel C tumor/normal label |
| SB2 umap_signatures | Panel D proliferation signature presence |
| SB3 umap_thyrocytes | Panel D per-patient distribution form |
| SC1 per_patient_r_forest | Per-patient n_cells annotation |
| SC2 multisite_trajectory | Panels A + C exact form |
| SC3 umap_multisite | Integration method (Harmony / scVI / other), Panel D thyrocyte preservation |
| SC4 gse184362_summary | Direct Pu 2021 metadata extraction |
| SC5 p2a_external_validation | Lu 2023 cohort overlap status |
| SC6 p2a_pooled_scatter | Marginal density form |
| SC7 k2_vs_yoo | Yoo 2016 per-sample availability; calibration diagnostic |
| SD1 landa2016_gse76039 | Exact GSE76039 n; convergence note phrasing |

---

## Marathon mode compliance

- ✅ No analysis execution
- ✅ No figure generation
- ✅ No manuscript prose (captions = factual technical descriptions, NOT voice-protected sections per `v17_sprint_vs_marathon_violation`)
- ✅ Cite verification deferred (Landa 2016 confirmed via memory; SC4/Pu 2021 metadata verify TBD)
- ✅ Captions ready for review against actual figures; user keyboard for any narrative additions

# END SUPP CAPTIONS — 17 figures captioned (templates with verification markers). Edit per figure-by-figure verification.
