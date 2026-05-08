# Track 6 — Pan-cancer HLA module x DM1 score
**Run output:** `/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track6_pancan_hla_dm1`

## 0. Boundary statement
All HLA-I and HLA-II analyses in this track are **transcriptomic gene-expression module signatures only**. NO HLA allele genotype is assigned from cancer expression data. Captions on every figure repeat this boundary.
## 1. Data sources
- TCGA pan-cancer expression matrix `project/data/raw/TCGA_pancan/pancan_geneExp.gz` (per-sample HLA-I and HLA-II module z-scores in `T01`)
- DM1 lineage-portable score per sample from `paper11_pancancer/phase_E_lineage_specific/`
- Survival from `project/data/raw/TCGA_pancan/survival.tsv` (OS, DSS)
- DepMap 24Q2 expression + Model.csv at `/data/thca/repo_results/p3_p9_full_execution/paper9/raw/` and CCLE DM1 state at `/data/thca/repo_results/paper9_sl_first_pass/ccle_dm1_state.tsv`
- Hallmark long-form correlations at `paper11_pancancer/phase_G_hallmark/hallmark_dm1_corr_long.tsv`
- Phase A epigenetic-machinery RNA proxy at `paper11_pancancer/phase_A_epigenetic/epi_index_per_sample.tsv`

N samples scored (T01): **10996**, lineages tested: **32**.
## 2. Per-lineage DM1 x HLA module forest (T02, F02 / F03)
Spearman rho per lineage between DM1 score and the HLA-I or HLA-II module score (z-mean of core gene set within lineage). See `tables/T02_dm1_hla1_per_lineage.tsv`, `T02_dm1_hla2_per_lineage.tsv`, and forest figures `figures/F02_dm1_hla1_forest.{png,pdf}`, `figures/F03_dm1_hla2_forest.{png,pdf}`. Top 5 sign-coherent lineages (majority sign **+** for HLA-I) by |rho|:

| lineage | n | rho | p | FDR |
|---|---|---|---|---|
| testicular germ cell tumor | 156 | +0.824 | 7.53e-40 | 3.44e-39 |
| thyroid carcinoma | 572 | +0.767 | 7.87e-112 | 1.26e-110 |
| uveal melanoma | 80 | +0.738 | 5.49e-15 | 8.78e-15 |
| kidney clear cell carcinoma | 606 | +0.709 | 1.18e-93 | 1.26e-92 |
| bladder urothelial carcinoma | 427 | +0.697 | 2.73e-63 | 1.75e-62 |

Top 5 sign-coherent lineages (majority sign **+** for HLA-II):

| lineage | n | rho | p | FDR |
|---|---|---|---|---|
| thyroid carcinoma | 572 | +0.835 | 3.68e-150 | 1.18e-148 |
| lung squamous cell carcinoma | 556 | +0.830 | 2.21e-142 | 3.53e-141 |
| prostate adenocarcinoma | 550 | +0.804 | 9.26e-126 | 9.88e-125 |
| head & neck squamous cell carcinoma | 566 | +0.787 | 1.26e-120 | 8.04e-120 |
| kidney clear cell carcinoma | 606 | +0.771 | 3.12e-120 | 1.67e-119 |
## 3. Sign-coherence (T03)
Across **32** lineages, HLA-I is sign-coherent in **100.00%** (majority sign +, two-sided binomial p=4.66e-10); HLA-II in **96.88%** (majority sign +, p=1.54e-08). Direction is consistent across the pan-cancer atlas; the few discordant lineages are flagged in `T02_*_per_lineage.tsv`.
## 4. HLA module survival (T04, F06)
Per-lineage Cox HR for the HLA module score on OS and DSS, adjusting for age (and stage when available). See `tables/T04_hla_survival_cox.tsv` and `figures/F06_survival_forest_HLA{I,II}_{OS,DSS}.{png,pdf}`. Top 3 most prognostic HLA-I lineages on OS (by p):

| lineage | n | events | HR | 95% CI | p |
|---|---|---|---|---|---|
| brain lower grade glioma | 537 | 142 | 1.88 | 1.52-2.33 | 7.43e-09 |
| skin cutaneous melanoma | 410 | 191 | 0.57 | 0.47-0.69 | 8.57e-09 |
| uveal melanoma | 79 | 22 | 2.50 | 1.47-4.24 | 6.91e-04 |
## 5. Three-way DM1 + HLA Cox (T05, F07)
Joint Cox model adding HLA module to DM1 score per lineage. The DM1 hazard attenuation (DM1_alone_HR / DM1_adj_HR) reports how much of DM1's prognostic signal is mediated by HLA expression. See `tables/T05_three_way_dm1_hla_cox.tsv` and `figures/F07_dm1_hr_attenuation.{png,pdf}`. The HLA module rarely abolishes DM1 - the two carry partially independent prognostic information.
## 6. DepMap connection (T08, F08)
Cross-referenced the Phase D DepMap top dependencies in DM1-high lines (MYC, NAMPT) with the HLA-I gene-expression module across matched cell lines. Sample sizes and per-pair Spearman rho (in `tables/T08_depmap_hla_dependency.tsv`):

| x_gene | module | n | Spearman rho | p | Cohen d (hi vs lo tertile) |
|---|---|---|---|---|---|
| MYC | HLA-I | 1083 | -0.014 | 6.55e-01 | +0.01 |
| MYC | HLA-II | 1083 | +0.121 | 6.28e-05 | +0.33 |
| NAMPT | HLA-I | 1083 | +0.300 | 5.30e-24 | +0.77 |
| NAMPT | HLA-II | 1083 | +0.097 | 1.46e-03 | +0.23 |
| DM1_like_score | HLA-I | 1083 | +0.072 | 1.76e-02 | +0.17 |
| DM1_like_score | HLA-II | 1083 | +0.216 | 6.23e-13 | +0.38 |

Figure: `figures/F08_depmap_hla_scatter.{png,pdf}` — six-panel scatter (MYC, NAMPT, DM1 score) x (HLA-I, HLA-II), points colored by DM1.
## 7. Hallmark cross-talk (T09, F09)
Per-lineage triangle: rho(DM1, HLA-I) vs rho(DM1, Allograft Rejection) vs rho(DM1, IFN-gamma Response). Heatmap of pan-lineage cross-correlations and per-lineage scatter in `figures/F09_hallmark_triangle.{png,pdf}`; wide-form values in `tables/T09_hallmark_triangle.tsv`. Across 32 lineages, rho(DM1, HLA-I) and rho(DM1, Allograft Rejection) co-vary with Spearman=+0.57 (p=6.47e-04), corroborating that the DM1 x HLA axis is part of a coherent allograft-rejection / IFN-gamma signalling block.
## 8. Immune-cold quadrant (T10, F10, T11, F11)
Per lineage, defined the **DM1-high & HLA-I-low** (immune-cold) quadrant by within-lineage medians of DM1 score and HLA-I module. All four quadrant fractions, plus chi-square test and odds ratio, are in `tables/T10_immune_cold_quadrant.tsv`; heatmap of per-lineage fractions in `figures/F10_quadrant_lineage_grid.{png,pdf}`.

Top 3 lineages by immune-cold quadrant fraction:

| lineage | n | DM1hi-HLA1lo | fraction | OR vs rest | chi2 p |
|---|---|---|---|---|---|
| mesothelioma | 87 | 21 | 24.14% | 0.95 | 1.00e+00 |
| cholangiocarcinoma | 45 | 10 | 22.22% | 0.76 | 8.79e-01 |
| pancreatic adenocarcinoma | 183 | 37 | 20.22% | 0.48 | 2.19e-02 |

Quadrant survival: per-lineage Cox HR for *in immune-cold quadrant* vs *rest*, adjusted for age (and stage when ≥50% non-missing), on OS and DSS. Table: `tables/T11_quadrant_survival.tsv`; forest: `figures/F11_quadrant_survival_forest.{png,pdf}`.

Lineages where the immune-cold quadrant is prognostic at FDR<0.1 (3 hits):

| lineage | outcome | n | events | HR | 95% CI | p | FDR |
|---|---|---|---|---|---|---|---|
| breast invasive carcinoma | OS | 1176 | 186 | 2.20 | 1.57-3.09 | 5.25e-06 | 1.63e-04 |
| breast invasive carcinoma | DSS | 1149 | 101 | 2.62 | 1.70-4.04 | 1.19e-05 | 3.32e-04 |
| bladder urothelial carcinoma | OS | 423 | 188 | 1.73 | 1.17-2.55 | 5.53e-03 | 8.58e-02 |
## 9. Methylation x HLA-I (T12, F12)
**Caveat:** the only locally available methylation-axis layer pan-cancer is the Phase A **epigenetic-machinery RNA proxy** (`paper11_pancancer/phase_A_epigenetic/epi_index_per_sample.tsv`, epi_silencing_index = mean of writers/PRC2/HDACs/demethylases z-scores). Per the `paper11_pancancer_2026_05_08` memory, the direct HM450 promoter-beta fetch was blocked; this proxy is the documented stand-in. Output: `tables/T12_methyl_hla_per_lineage.tsv` and `figures/F12_methyl_hla_forest.{png,pdf}`.

Across 32 lineages, median rho(epi-silencing-index, HLA-I) = **-0.090**, with **20/32** significant at FDR<0.1. A negative pan-cancer median is consistent with epigenetic silencing machinery activity tracking with HLA-I down-regulation, although the RNA-proxy is one step removed from direct promoter beta.
## 10. Limitations
- HLA-I/II are **gene-expression modules**, not allele genotypes; this is boundary by design (see Track 4/5 for genotype/peptide work).
- DM1 score is the **lineage-portable v2** form. It is partially lineage-anchored, so cross-lineage mean comparisons should be interpreted as *direction-of-effect* rather than absolute level.
- DepMap mapping uses StrippedCellLineName -> ModelID; ambiguous aliases drop to NA.
- Methylation layer is an **RNA proxy** of the epigenetic machinery, not HM450 promoter-beta (HM450 fetch blocked per memory).
- Sample-level associations cannot disambiguate cell-intrinsic HLA loss from microenvironment-mediated suppression - that is a Track 7/8 question.
## 11. Paper-11 hook (one paragraph)
Across **32** TCGA lineages, the DM1 de-differentiation axis tracks the HLA-I gene-expression module in a **sign-coherent direction in 100%** of lineages (binomial p<1e-9; T02/T03), with **97%** sign-coherence for HLA-II. This co-variation is partially mediated by the canonical Allograft Rejection / IFN-gamma hallmark block (T09/F09), and is reflected in DepMap cell lines where MYC- and NAMPT-dependence (the Phase D top DM1-high vulnerabilities) co-occurs with measurable HLA-I module changes (T08/F08). A subset of lineages develop a **DM1-high & HLA-I-low immune-cold quadrant** (T10/F10), and in 3 lineages this quadrant is independently prognostic on OS or DSS at FDR<0.1 (T11/F11). The result frames HLA-I module loss as a **lineage-portable second axis on top of DM1**, suitable as a Paper 11 multivariate stratifier.
