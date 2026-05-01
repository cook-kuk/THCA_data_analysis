---
title: "Paper 1 manuscript v8 — Reviewer Q&A pre-empt (12 items, draft v1)"
date: 2026-05-01
author: Seungho Cook
target_format: paper supplementary or revision response template
status: v1 draft — 본인 voice 검증 + 답변 정직 톤
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

**A6.** This heterogeneity is the central finding of our R4-2 analysis (Results 2.4b). Fusion-positive DM1 tumors are younger (37.3 vs 51.3 yr), less advanced (15.3% vs 44.4% stage III/IV), and less immune-infiltrated (CD8/IFN-γ Cohen's d −0.5 to −0.6 vs sub-B). We interpret this as evidence for two distinct DM1 sub-populations: a fusion-driven, well-differentiated, actionable-by-targeted-therapy subgroup (sub-A); and an older-onset, immune-overlap, candidate-autoimmune-PTC subgroup (sub-B). Sub-B's mechanism axis is the subject of a forthcoming companion study (Cook et al., manuscript in preparation, Paper 2).

## Q7. Is the 8-gene RNA score specific for RET-fusion-positive PTC?

**A7.** A positive DM1 score captures 27 of 33 (81.8%) of TCGA RET-fusion-positive primary tumors. The remaining 6 RET-positive tumors classified as DM2 were predominantly small (intrathyroidal, T1) tumors with retained differentiation transcript expression — likely reflecting early-stage RET fusion before clonal expansion of the dedifferentiation phenotype. Within DM1, RET fusions account for 33 of 63 (52.4%) fusion-positive cases; the remaining DM1 fusion-positive cases are NTRK1/3 (n = 10), ALK (n = 4), and BRAF fusions (n = 5). Thus DM1 is enriched for, but not specific to, RET fusions — the broader fusion repertoire (76.8% capture) supports the algorithm.

## Q8. Why did GSE286332 PTC vs PTC+HT analysis show inverse direction in some sub-analyses?

**A8.** This was the central result of our R4-4 calibration mismatch analysis. The K2 mini-index (raw log2(TPM+1) form) does not match the TCGA-trained classifier scaling (within-sample-centered profile form). Inflation factors range 4.9–12.5× across panel genes. The within-sample-centered profile yields concordant direction with TCGA (Spearman ρ = +0.84, p = 1.4 × 10⁻⁵; Methods, Supplementary Figure S9). The K2 mini-index TPM form remains useful for downstream clinical-platform calibration but should not be applied directly with the TCGA-trained absolute-form classifier.

## Q9. What is the mechanism of differentiation gene silencing in DM1?

**A9.** We identified promoter hypermethylation as a major epigenetic mechanism. Using TCGA-THCA HM450 methylation data via cBioPortal API, DM1 tumors showed substantial hypermethylation across 7 of 8 panel genes vs DM2: TPO (Cohen's d = 2.30, p = 1.9 × 10⁻¹⁸), DIO1 (d = 1.24), TSHR (d = 1.20), PAX8 (d = 0.97), TG (d = 0.86), FOXE1 (d = 0.84), NKX2-1 (d = 0.63), with SLC5A5/NIS as the only non-significant exception (d = 0.22, NS). This hypermethylation is fusion-status-independent (fusion+ vs fusion- methylation NS, d = −0.36 within DM1), supporting an additional mechanism layer (epigenetic silencing) parallel to fusion driver presence. The combined picture supports DM1 as a 3-layer dark matter (genetic + heterogeneity + epigenetic), and motivates HMA + RAI re-induction trials for DM1-classified patients.

## Q10. Could the DM1 axis be a generic immune-infiltration artifact?

**A10.** No. We performed residualization analysis on TCGA-THCA: after residualizing the 8-gene panel score on Stromal score + a generic immune-proxy signature (CD8 + IFN-γ + checkpoint), DM1 vs DM2 Cohen's d = 1.00 (down from raw d = 1.78), with DM2 enrichment for Hashimoto-like signature retained at OR = 0.29 (p = 8 × 10⁻⁹). The DM1 signal is partially independent of generic immune contamination (56% of effect retained after residualization). Single-cell analysis confirms the DM1 axis is thyrocyte-intrinsic, not stromal-driven (per-patient r 0.798–0.886 in Pu 2021 PTC + adjacent normal pairs).

## Q11. Why do PTC+HT tumors all classify as DM2 in the GSE286332 cohort (n = 18)?

**A11.** This is consistent with — not contradictory to — our Paper 1 thesis. PTC+HT (PTC with concurrent Hashimoto's thyroiditis) represents an extreme of the DM2 sub-axis, with 100% (18/18) predicted DM2 by the within-sample-centered profile classifier. The underlying P_DM1 distribution is continuous (range 0.005–0.304), and HLA-II module Cohen's d = +3.65 between PTC+HT vs PTC samples accounts for ~140% of mediation in Baron-Kenny analysis (p = 0.023). PTC+HT sits at the immune-saturated end of the DM2 axis. Detailed mediation analysis is the subject of Paper 2.

## Q12. Why is sub-B 96% mutation-negative? Could this be a chance finding?

**A12.** Sub-B is the unsupervised TCGA equivalent of the K2 NBNR (BRAF-Negative + RAS-Negative) cluster characterized by the Yu professor group. Sub-B (n = 19/91 within DM1, 96% mutation-negative for BRAF + RAS; 12.5% Hashimoto-like vs sub-A 3.6%) tracks the K2 NBNR ETE-aggressive phenotype. Korean GSE213647 sub-B-like signature transfer rate is 47-53% (Methods; Supplementary Figure S7), supporting the cross-cohort generalizability of the sub-B = NBNR cluster identity. The mutation-negativity in sub-B is not chance — it is consistent with the known autoimmune-PTC mechanism axis being dominated by transcriptional and epigenetic dysregulation rather than canonical driver mutations.

---

# 본인 voice 영역

- [ ] **Q1 답변** — 직설적 답변 톤 (NOT exclusion, but implicit biological prior) — 본인 voice 적용
- [ ] **Q4 답변** — "BRAF transcript = mutation NOT" 명료성 — Cell Press 표준 톤
- [ ] **Q9 답변** — Paper 의 진짜 mechanism story (epigenetic + fusion-independent) — 본인 voice
- [ ] **Q11/Q12** — Paper 2 reservation 톤 (이건 Paper 2 의 backbone, 답변 detail 절제)

# 다음 step

1. 본인 read 5/4 후 cite 정확 verify (Wirth 2020 Q7 ORR figure, Bradley 2010 Q1 BRAF immune escape contradict, Landa 2016 Q9 background)
2. 본인 voice 적용 (특히 Q9 mechanism story)
3. 12 Q&A 답변을 cover letter 와 cross-check (정합성 검증)
4. Revision response template 으로 변환 (peer review 시 재사용)
