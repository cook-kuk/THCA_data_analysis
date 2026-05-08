---
title: "Paper 1 manuscript v8 — Reviewer Q&A pre-empt"
date: 2026-05-04
author: Seungho Cook
target_format: paper supplementary or revision response template
status: clean draft with author voice placeholder
purpose: 12 likely reviewer questions 사전 답변 — paper revision 시 정확 표현 재사용 가능
---

# Reviewer Q&A pre-empt (12 items)

---

## Q1. Why exclude BRAF/RAS/TERT from the candidate gene pool?

**A1.** They were *not* excluded. The 67-gene candidate pool (TIERA67) explicitly includes BRAF, NRAS, HRAS, KRAS, RET, NTRK1/3, ALK, PAX8, PPARG, TERT, and EIF1AX in the `[Driver_anchor]` sub-category (12 genes; Supplementary Table S1). The 8-gene panel was selected from the `[TDS-core]` sub-category (16 genes representing canonical thyroid differentiation transcripts) based on Yoo et al. (2016) RAI-uptake biology, not driver mutational status. This curation choice reflects implicit biological prior (differentiation axis vs MAPK driver axis), not exclusion. The 8-gene panel + Driver_anchor concatenation was tested and yielded ARI = 0.49 alone vs 0.90 with full TIERA67 — Driver_anchor genes alone yielded ARI = −0.007, confirming they cannot define the DM cluster axis.

## Q2. Did you bias the panel toward iodine metabolism genes?

**A2.** Yes, by explicit design. The 8-gene panel was selected from `[TDS-core]` because canonical RAI uptake biology (NIS/SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) is the clinical differentiation axis with the strongest biological prior (Yoo et al., 2016). This bias is appropriate for the clinical question (RAI decision-making) and is leak-free in cross-validation: 5-fold AUC = 0.962 (8-gene) vs 0.975 (16-gene full TDS-core), ΔAUC = 0.013 NS. R1-B leak control analysis confirmed AUC = 0.925 in held-out cohorts.

## Q3. How does this panel differ from the BRAF-RAS Score (BRS, 273 genes)?

**A3.** Spearman ρ between 8-gene RAI score and Yoo 2016 BRS in the K2 cohort is 0.49 — partial correlation (not redundant). Concordance is 77%; discordant cases (23%) constitute the BRAF/RAS-negative dark matter axis. The 8-gene panel captures a differentiation axis orthogonal to BRS. Within Xing 2014 BRAF/TERT-negative dark matter, the 8-gene panel sub-stratified 131 of 180 (73%) into a defined molecular axis (DM1 vs DM2).

## Q4. Why didn't BRAF V600E rank near the top in unsupervised analysis?

**A4.** BRAF V600E is a mutation, not a transcript. BRAF transcript expression is comparable in V600E mutant versus wild-type tumors (Cohen's d = −0.044, Mann-Whitney p = 0.57; n = 273 vs 182 in TCGA-THCA). Single-feature DM classification AUC for driver transcripts ranged 0.50–0.61 across BRAF, HRAS, NRAS, KRAS, and TERT — none reaching the discriminative threshold of the 8-gene panel (AUC = 0.962). The 8-gene panel reflects a transcriptional differentiation axis, not driver mutation status.

## Q5. Is the DM1 fusion finding robust to structural variant (SV) missingness?

**A5.** SV-tested status was available for 542 of 557 TCGA-THCA tumors (97.3%); the missing 15 tumors were tested for missingness × DM cluster correlation (chi² p = 0.56, missing-at-random). Three sensitivity analyses imputing missing SV status (best-case fusion-negative, worst-case fusion-positive, case-control matched) yielded DM1 vs DM2 fusion ORs of 7.18, 9.07, and 7.41, respectively — all retaining statistical significance (p < 10⁻¹⁰).

## Q6. Why are DM1 fusion-positive vs fusion-negative tumors so different in clinical phenotype?

**A6.** This heterogeneity is the central finding of our R4-2 analysis (Results 2.4b). Fusion-positive DM1 tumors are younger (37.3 vs 51.3 yr), less advanced (15.3% vs 44.4% stage III/IV), and less immune-infiltrated (CD8/IFN-γ Cohen's d −0.5 to −0.6 vs sub-B). We interpret this as evidence for two distinct DM1 sub-populations: a fusion-driven, well-differentiated, actionable-by-targeted-therapy subgroup (sub-A); and an older-onset, immune-overlap subgroup (sub-B) whose deeper mechanism is reserved for companion-study analysis.

## Q7. Is the 8-gene RNA score specific for RET-fusion-positive PTC?

**A7.** A positive DM1 score captures 27 of 33 (81.8%) of TCGA RET-fusion-positive primary tumors. The remaining 6 RET-positive tumors classified as DM2 were predominantly small (intrathyroidal, T1) tumors with retained differentiation transcript expression — likely reflecting early-stage RET fusion before clonal expansion of the dedifferentiation phenotype. Within DM1, RET fusions account for 33 of 63 (52.4%) fusion-positive cases; the remaining DM1 fusion-positive cases are NTRK1/3 (n = 10), ALK (n = 4), and BRAF fusions (n = 5). Thus DM1 is enriched for, but not specific to, RET fusions — the broader fusion repertoire (76.8% capture) supports the algorithm.

## Q8. Why did one small external Korean cohort show inverse direction in some sub-analyses?

**A8.** This was the central result of our R4-4 calibration mismatch analysis. The K2 mini-index (raw log2(TPM+1) form) does not match the TCGA-trained classifier scaling (within-sample-centered profile form). Inflation factors range 4.9–12.5× across panel genes. The within-sample-centered profile yields concordant direction with TCGA (Spearman ρ = +0.84, p = 1.4 × 10⁻⁵; Methods, Supplementary Figure S9). The K2 mini-index TPM form remains useful for downstream clinical-platform calibration but should not be applied directly with the TCGA-trained absolute-form classifier.

## Q9. What is the mechanism of differentiation gene silencing in DM1?

**A9.** We treat the methylation finding as correlative rather than mechanistic. The TCGA HM450 data establish that DM1 tumors carry markedly elevated promoter β-values across the 8-gene panel (mean β 0.385 vs DM2 0.253; TPO Cohen's d = 2.30), and that within DM1 this signature is fusion-independent at the cohort level (n = 63 fusion-positive vs n = 19 fusion-negative; d = −0.36, NS). What our data do not establish is the upstream regulator. Plausible candidates include canonical de novo methyltransferases (DNMT3A/3B), histone-methylation-coupled silencing complexes such as SETDB1, and MAPK-effector chromatin-remodeling pathways previously implicated in thyroid dedifferentiation; we did not perform DNA methyltransferase profiling, ChIP-seq for repressive marks (H3K9me3/H3K27me3), or CRISPRi knockdown of candidate regulators in this study. We therefore restrict our claim to the observation of a fusion-independent, DM1-associated promoter hypermethylation pattern that motivates — rather than confirms — hypomethylating-agent-based RAI re-induction strategies. Causal attribution of the silencing program to a specific upstream effector is left to dedicated mechanism work in cell-line and PDX systems, and we explicitly hold the line that promoter methylation in our data is suggested as a parallel mechanism, not proven as a causal one.

## Q10. Could the DM1 axis be a generic immune-infiltration artifact?

**A10.** No. We performed residualization analysis on TCGA-THCA: after residualizing the 8-gene panel score on Stromal score + a generic immune-proxy signature (CD8 + IFN-γ + checkpoint), DM1 vs DM2 Cohen's d remained 1.00 (down from raw d = 1.78). The DM1 signal is therefore partially independent of generic immune contamination, with 56% of the original effect retained after residualization. Single-cell analysis confirms that the DM1 axis is thyrocyte-intrinsic rather than stromal-driven (per-patient r 0.798–0.886 in Pu 2021 PTC + adjacent normal pairs).

## Q11. Why do some small external Korean reference samples map strongly to DM2?

**A11.** This is consistent with — not contradictory to — our Paper 1 thesis. The DM1/DM2 axis is continuous, not binary in its underlying biology, and the small external Korean reference set occupies the low-score end of the same classifier used in TCGA. In that set, P(DM1) values range from 0.005 to 0.304, which is compatible with preserved differentiation-program expression rather than failure of the classifier. Because this reference set is used here only for calibration and score-portability checks, we do not assign additional biological interpretation in the present paper.

## Q12. Why is sub-B 96% mutation-negative? Could this be a clustering artifact?

**A12.** Sub-B is the fusion-negative component of DM1 identified by unsupervised clustering within the DM1 compartment, with silhouette support (0.584) and reproducible cross-cohort DM score geometry in Korean validation data. Its mutation-negativity indicates that the DM1 axis is not reducible to canonical BRAF/RAS driver status. We therefore interpret sub-B as a biologically distinct transcriptional state within DM1, while treating its more specific mechanism as hypothesis-generating and reserving deeper characterization for companion work.

---
