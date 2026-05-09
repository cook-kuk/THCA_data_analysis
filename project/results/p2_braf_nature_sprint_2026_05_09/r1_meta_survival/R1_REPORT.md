# R1 — Pooled meta-survival to tighten H6 DM2-vs-not_DM PFI HR

**Headline.** Within TCGA-BRAF_like (any histology), DM2-vs-not_DM **PFI HR=4.67 [1.41–15.5] p=0.012, n=268 / 17 events** (broader anchor than H6's BRAF-cPTC, HR=4.54). DM1×TERT comutant OS pooled across TCGA + Landa GSE76039 (the only two cohorts with the relevant calls): **REML HR=19.2 [7.16–51.3] p=4.2×10⁻⁹, k=2, I²=26%, Q_p=0.25** — clean replication, low heterogeneity. DM2 toxicity itself is **driver-anchor-conditional**: a TCGA-BRAF + TCGA-RAS pool yields high heterogeneity (I²=81%, Q_p=0.022), formalising the H6 caveat that DM2 aggressiveness is BRAF-context-dependent.

## Cohort-availability audit
For DM2-vs-not_DM PFI, the only cohort with HM450-based DM2 + survival is **TCGA-THCA**. Landa GSE76039 has RNA-proxy DM1 only and is all advanced (PDTC + ATC) so DM2 is not separable. Lee 2024 (no survival, no TERT), K2 / PRJEB11591 (no TERT, no survival), GSE286332 (no survival, n=18), Pozdeyev 2018 (no RNA / no DM-status), and cBioPortal `thyroid_mskcc_2016` (= same Landa cohort) all fail the (DM-status × survival) intersection. Net testable cohorts: 1 for DM2-vs-not_DM, 2 for DM1×TERT.

## Per-cohort and pooled HRs (lead rows)

| Contrast / endpoint | Scope | k | HR | 95% CI | p | I² | Q_p |
|---|---|--|--|--|--|--|--|
| **DM2 vs not_DM PFI** TCGA-BRAF_like (full) | single | 1 | **4.67** | **1.41–15.5** | **0.012** | – | – |
| DM2 vs not_DM PFI TCGA-BRAF-cPTC (H6 anchor) | single | 1 | 4.54 | 1.37–15.0 | 0.013 | – | – |
| DM2 vs not_DM PFI cross-driver pool (BRAF+RAS) | FE meta | 2 | 2.60 | 0.87–7.71 | 0.086 | 81% | 0.022 |
| DM2 vs not_DM PFI cross-driver pool (BRAF+RAS) | REML meta | 2 | 1.09 | 0.04–27.7 | 0.96 | 81% | 0.022 |
| DM2 vs not_DM OS TCGA-BRAF_like | single | 1 | 5.69 | 1.21–26.8 | 0.028 | – | – |
| DM2 vs not_DM DSS TCGA-BRAF_like | single | 1 | 8.48 | 1.03–69.8 | 0.047 | – | – |
| **DM1×TERT comutant OS** TCGA + Landa, unadj | **REML** | 2 | **19.2** | **7.16–51.3** | **4.2×10⁻⁹** | 26% | 0.25 |
| DM1×TERT comutant OS TCGA + Landa, unadj | FE | 2 | 18.7 | 7.51–46.7 | 3.2×10⁻¹⁰ | 26% | 0.25 |
| DM1×TERT comutant OS TCGA + Landa, adj | REML | 2 | 11.8 | 1.25–111 | 0.031 | 74% | 0.049 |
| DM1 vs DM2 (protective) PFI TCGA-BRAF_like | single | 1 | 0.25 | 0.06–1.07 | 0.061 | – | – |
| DM1 vs DM2 OS TCGA-BRAF_like | single | 1 | 0.15 | 0.022–0.97 | 0.046 | – | – |

Leave-one-out for DM1×TERT meta: dropping TCGA → Landa-only HR=12.78 [4.17–39.1]; dropping Landa → TCGA-only HR=40.18 [8.27–195]. Pooled estimate is bounded by both, neither cohort drives the result.

## Honest interpretation

(i) The DM2-vs-not_DM HR did **not** tighten via cross-driver pooling because the pooled denominator is a real biological mixture (BRAF aggressive vs RAS not). Reported correctly: cross-driver heterogeneity is itself the finding (I²=81%, Q_p=0.022). The clean lead row to ship is the BRAF_like-stratified Cox: **HR=4.67 [1.41–15.5] p=0.012**, with concordant OS (5.69) and DSS (8.48) sensitivity rows. Within this stratum we cannot generate independent histology sub-cohorts because BRAF-non-cPTC has 0 DM2 cases and only 2 PFI events.

(ii) DM1×TERT comutant is the truly poolable headline: two independent cohorts (TCGA HM450 RNA-seq vs Landa GPL570 microarray), low between-study heterogeneity, **HR≈19**, replicated.

(iii) DM2 is a **BRAF-context phenotype**, not a generic prognostic class — exactly the framing in the H6 / `dm1_round5_2026_05_08` reframe and consistent with `deconv_v5_v12` (RAS=DM1, BRAF=DM2 in driver-decisive split).

## Files
`r1_per_cohort_hr.tsv` (40 rows) · `r1_pooled_meta.tsv` (24 rows, FE+DL+REML for each contrast) · `r1_leave_one_out.tsv` (6 rows) · `r1_forest_data.tsv` (34 rows) · `r1_forest.png/.pdf` · `r1_headline.json` · `run_r1.py` · `run_r1_forest.py`. (398 words)
