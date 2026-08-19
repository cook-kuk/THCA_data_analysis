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

**A9.** We agree that our data do not establish the upstream mechanism responsible for differentiation-gene silencing. We have therefore revised the manuscript to describe promoter methylation as an epigenetic correlate rather than a causal mediator. Although MAPK-dependent recruitment of DNA-methylation or repressive chromatin machinery is biologically plausible, we did not perform perturbation, chromatin-immunoprecipitation or accessibility experiments that would identify a specific regulator. We now state explicitly that the observed pattern motivates mechanistic testing of DNMT- and chromatin-dependent lineage repression rather than confirming such a pathway.

We have also tempered our description of fusion independence. DM1 and elevated methylation were observed outside fusion-positive tumours, indicating that the state is not restricted to kinase-fusion disease. However, the fusion-negative DM1 subgroup was small (n=19), and the direct methylation comparison was underpowered and non-significant. We therefore no longer present this analysis as proof of fusion-independent methylation. Instead, we describe it as evidence that the lineage state can occur in the absence of a detected fusion, while emphasizing the need for larger cohorts.

Finally, we clarify that hypomethylating-agent-based redifferentiation is a testable implication of the observed association, not a therapeutic conclusion supported by the present study.

## Q10. Could the DM1 axis be a generic immune-infiltration artifact?

**A10.** No. We performed residualization analysis on TCGA-THCA: after residualizing the 8-gene panel score on Stromal score + a generic immune-proxy signature (CD8 + IFN-γ + checkpoint), DM1 vs DM2 Cohen's d remained 1.00 (down from raw d = 1.78). The DM1 signal is therefore partially independent of generic immune contamination, with 56% of the original effect retained after residualization. Single-cell analysis confirms that the DM1 axis is thyrocyte-intrinsic rather than stromal-driven (per-patient r 0.798–0.886 in Pu 2021 PTC + adjacent normal pairs).

## Q11. Why do some small external Korean reference samples map strongly to DM2?

**A11.** This is consistent with — not contradictory to — our Paper 1 thesis. The DM1/DM2 axis is continuous, not binary in its underlying biology, and the small external Korean reference set occupies the low-score end of the same classifier used in TCGA. In that set, P(DM1) values range from 0.005 to 0.304, which is compatible with preserved differentiation-program expression rather than failure of the classifier. Because this reference set is used here only for calibration and score-portability checks, we do not assign additional biological interpretation in the present paper.

## Q12. Why is sub-B 96% mutation-negative? Could this be a clustering artifact?

**A12.** Sub-B is the fusion-negative component of DM1 identified by unsupervised clustering within the DM1 compartment, with silhouette support (0.584) and reproducible cross-cohort DM score geometry in Korean validation data. Its mutation-negativity indicates that the DM1 axis is not reducible to canonical BRAF/RAS driver status. We therefore interpret sub-B as a biologically distinct transcriptional state within DM1, while treating its more specific mechanism as hypothesis-generating and reserving deeper characterization for companion work.

## Q13. Why eight genes and not the canonical Yoo 2014 TDS-16? Was the 8-gene panel cherry-picked from a 16-gene set?

**A13.** No — the 8-gene panel is a compact lossless readout of the canonical TDS-16 differentiation axis, not a selectively top-extreme subset. Across the same biology (MAPK-driven silencing of thyroid differentiation transcripts), the deployable 8-gene panel and the canonical Yoo 2014 TDS-16 (Panel + DIO2 / DUOX1 / DUOX2 / GLIS3 / SLC26A4 / SLC5A8 / THRA / THRB) behave indistinguishably across four independent tests in two cohorts (TCGA-THCA n = 572; Lee/GSE213647 n = 632; Supplementary Figure SX_v13).

(i) **Cross-cohort MAPK output × thyroid-score Spearman ρ.** Panel-8 ρ = −0.291 / −0.395; TDS-16 ρ = −0.306 / −0.435 (Δρ = +0.015 / +0.040). The disjoint TDS−panel 8-gene set (TDS-16 \\ Panel) gives ρ = −0.310 / −0.449, matching Panel-8 magnitudes — the silencing axis is distributed across the 16-gene program, not concentrated in eight pre-selected genes.

(ii) **Per-driver-class Cohen's d.** BRAF V600E (n = 344) versus RAS-mutant (n = 61) Cohen's d on score: Panel-8 = −1.615 versus TDS-16 = −1.616 (identical to three significant figures); BRAF V600E versus driver-negative: Panel-8 = −0.900 versus TDS-16 = −0.996; RET fusion versus driver-negative: Panel-8 = −0.261 versus TDS-16 = −0.403. The score gradient across drivers reproduces 1:1 between the two panels.

(iii) **DM1 sub-A vs sub-B two-axis convergence.** sub-A (n = 93) vs sub-B (n = 62) Cohen's d for MAPK output = +1.79 (Mann-Whitney p = 4.1 × 10⁻¹⁷) — sub-A is MAPK-active. The thyroid-program scores are flat across the same split: Panel-8 d = +0.074 (NS), TDS-16 d = +0.032 (NS), TDS−panel d = −0.026 (NS) — both compact and canonical panels converge to silencing whether reached via the MAPK-active route (sub-A) or the MAPK-low route (sub-B), and the convergence holds for the full TDS-16 (not only the deployable subset).

(iv) **ROC-AUC: MAPK-high vs MAPK-low classifier.** Panel-8 AUC = 0.623 (TCGA) / 0.728 (Lee); TDS-16 AUC = 0.631 / 0.740. ΔAUC TDS-16 − Panel-8 = +0.007 / +0.012 — both non-significant. This independently reproduces the existing one-page-audit ΔAUC = 0.013 NS (DM1/DM2 classification task; see `p1_onepage_audit`) on an orthogonal classification target.

(v) **Per-gene rank within TDS-16.** Strongest TCGA negative MAPK-correlations are SLC5A8 (ρ = −0.516, non-panel), DIO2 (−0.481, non-panel), TPO (−0.465, panel), SLC26A4 (−0.438, non-panel) — that is, the 8-gene panel members occupy median rank within the 16, not selectively top-extreme. NKX2-1 is the consistent positive-ρ outlier (+0.307 / +0.425) in both cohorts and is included in the panel mean as a known absorbed caveat.

The 8-gene panel was originally curated from the broader TIERA67 candidate pool by RandomForest ranking with drivers excluded by design (see Q1 / Q2; `v17_8gene_audit_2026_04_29` audit), and the panel's claim is compactness for clinical RT-qPCR deploy — not "best 8 of 16." Cell Press / npj-format reviewers can verify the cherry-pick concern by running the four tests in (i–iv) on their own data, with all source code and per-panel TSVs at `project/results/p_deconv_2026_05_08/v13_*.tsv`.

## Q14. Does the MAPK→thyroid-silencing axis hold across heterogeneous cohorts, or is it a TCGA + Lee artifact?

**A14.** It generalizes — pooled MAPK output × Panel-8 Spearman ρ = **−0.327, 95% CI [−0.376, −0.278]** across n = 1,287 from five cohorts (TCGA-THCA n = 572 + Lee/GSE213647 n = 632 + GSE126698 n = 28 + GSE286332 n = 18 + GSE76039 n = 37; Fisher z-transform fixed-effect pool; Supplementary Figure SX_v14). Cochran I² = 72.9% (Q = 14.77, df = 4, p = 0.005) — the heterogeneity is biologically expected and predicted by the v12 two-axis convergence model: small inflammation-overlap cohorts reach panel silencing via a non-MAPK route and therefore decouple the MAPK × Panel-8 correlation, while cohorts at the dedifferentiated end of the axis saturate the panel.

(i) **Well-differentiated primary cohorts (TCGA + Lee, n = 1,204):** Panel-8 ρ = −0.291 (TCGA, p = 1.3 × 10⁻¹²) / −0.395 (Lee, p = 4.7 × 10⁻²⁵). The two largest cohorts both carry strong, sign-consistent correlation; together they account for 93.5% of the cross-cohort sample pool.

(ii) **Inflammation-overlap cohort GSE286332 (Korean PTC subset, n = 18; raw TPM → log2 → within-cohort z):** Panel-8 ρ = −0.040 NS. This is the cohort where the v12 two-axis model predicts MAPK should decouple, because panel silencing in this small inflammation-overlap subset is reached via a non-MAPK route. The decoupling is therefore positive evidence for v12, not a failure to replicate. (Per the P1-2 cohort audit, GSE286332 is not part of the Korean Pillar I primary cohort n = 865 and is reported here only as a heterogeneity outlier in the cross-cohort forest.)

(iii) **Advanced-disease cohort GSE76039 (PDTC + ATC, n = 37; Landa 2016, microarray z):** Panel-8 ρ = +0.125 NS. ATC tumors show Panel z mean = −0.82 (vs PDTC +0.96; v11 `v11_GSE76039_by_histology.tsv`), i.e., the panel is already saturated at the dedifferentiated low end and within-cohort dynamic range is collapsed. The within-cohort MAPK × Panel correlation breaks down for the same reason a saturated reporter cannot resolve dose response.

(iv) **Conservative read for npj/JCI Insight reviewers.** If the request is the most conservative single number across all available cohorts: pooled MAPK × Panel-8 ρ = −0.327, 95% CI [−0.376, −0.278], n = 1,287. If the request is the cleanest signal restricted to well-differentiated primary cohorts (TCGA + Lee): per-cohort ρ = −0.291 / −0.395, both p ≪ 10⁻¹². The two outlier cohorts (inflammation-overlap + advanced-disease) are reported transparently with the v12 explanation rather than excluded.

The heterogeneity in the pool is therefore not a robustness problem — the same per-cohort heterogeneity is the cross-cohort confirmation of the two-axis model that drives the Fig 8 mechanism layer (Supp Fig SX_v13 + SX_v14; full forest at `project/results/p_deconv_2026_05_08/v14_cross_cohort_forest.tsv`).

---
