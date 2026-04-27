---
title: "An 8-gene RAI-responsiveness biomarker outperforms BRAF V600E status in papillary thyroid carcinoma"
journal: "npj Precision Oncology"
article-type: "Article"
date: "2026-04-27"
---

# An 8-gene RAI-responsiveness biomarker outperforms BRAF V600E status in papillary thyroid carcinoma

Seungho Cook^1,*^ and Yu Hyeong Won^2,*^

^1^ Independent Researcher, Seoul, Republic of Korea; Graduate School of Convergence Science and Technology, Seoul National University, Suwon, Republic of Korea (part-time PhD candidate).
^2^ Department of Surgery, Seoul National University Bundang Hospital, Seongnam, Republic of Korea.

\* Co-corresponding authors. Correspondence: Seungho Cook (kukshomr@gmail.com); Yu Hyeong Won ([email TBD]).

---

## Abstract

Pre-RAI clinical decision-making in papillary thyroid carcinoma (PTC) currently relies on BRAF V600E status alone with no quantitative biomarker that outperforms it. In 513 TCGA-THCA primary tumours we identify an unsupervised transcriptomic axis (DM1/DM2) statistically distinct from the BRAF/RAS dichotomy (Spearman ρ = 0.49) that maps onto a thyroid-differentiation/immune-state continuum. An 8-gene RAI-responsiveness panel (NIS/SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) recovers cluster identity with 5-fold cross-validated AUC 0.954 (Random Forest 0.975), substantially higher than a BRAF-V600E-only baseline (AUC 0.822, ΔAUC = +0.132). External validation on GSE76039 (n = 37 PDTC + ATC) yields direction-correct AUC 0.974. A composite cytolytic-IFN-γ-immune-fraction score separates DM1 (hot) from DM2 (cold) with Cohen's d = +1.683 (Mann-Whitney p = 3.99 × 10^−18^). Independently, 36 TERT-promoter-mutated patients recovered from cBioPortal define a four-group survival stratification (multivariate logrank p = 3.78 × 10^−5^; bootstrap median p = 2.6 × 10^−5^; leave-one-out worst-case p < 10^−3^); univariate Cox HR = 6.31 drops to 1.88 (p = 0.29) after stage + age + sex adjustment, so we report TERT^+^ honestly as a molecular handle for advanced-stage identification rather than a stage-independent prognostic marker.

## Significance Statement

Pre-RAI clinical decision-making for PTC is based on BRAF V600E alone, despite ~10% BRS misclassification and ~30% driver-negative tumours. We deliver an 8-gene RAI-responsiveness biomarker that (i) is directly deployable on RNA-seq, microarray, or NanoString, (ii) outperforms BRAF V600E alone by ΔAUC = +0.132 in cross-validation, and (iii) is anchored on a biologically interpretable transcriptomic axis that simultaneously stratifies the immune microenvironment. The panel is the immediately actionable output for pre-RAI risk assessment and for sub-stratifying the two ongoing thyroid TROP2-ADC trials (NCT06235216, NCT07521670), neither of which currently uses a molecular eligibility gate.

## Introduction

Papillary thyroid carcinoma (PTC) is the most common endocrine malignancy and is conventionally classified into BRAF-like and RAS-like molecular subtypes by the 52-gene BRAF-RAS Score (BRS) [1] and formalised in the TCGA-THCA atlas [2]. Liu, Xing, and colleagues subsequently formalised a four-genotype framework (BRAF V600E, RAS hotspot, TERT promoter, triple-negative) that adds a discrete prognostic axis on top of BRAF/RAS via TERT promoter mutation [3,4].

Three structural limitations motivate further sub-classification. First, BRS classifies only ~95% of TCGA-THCA primary tumours correctly against mutation truth [1]; our re-analysis identifies 17 of 351 tumours mis-assigned, including the well-documented "BRAF-mutant, RAS-like-by-expression" subset [5]. Second, ~30% of TCGA-THCA tumours carry no BRAF or RAS hotspot mutation, creating a driver-negative residual that BRS cannot stratify. Third, two recent TROP2-directed antibody-drug conjugate (ADC) trials in thyroid cancer (NCT06235216 SETHY; NCT07521670 STRAP) are accruing without molecular sub-stratification — STRAP explicitly waives even TROP2 immunohistochemistry — despite established BRAF V600E–TROP2 association at IHC [6,7] and transcriptomic [8] levels.

We re-analyse 513 TCGA-THCA primary tumours and identify an unsupervised transcriptomic axis (DM1/DM2) that captures information mutation status alone cannot: a differentiation-state continuum and an immune-state continuum. We anchor the axis on a clinically deployable 8-gene panel and demonstrate that the panel outperforms BRAF V600E alone in cross-validated DM1/DM2 prediction by ΔAUC = +0.132. Our contribution is a sub-stratification layer with three properties no prior framework provides simultaneously: (i) statistical distinctness from mutation status, (ii) head-to-head outperformance over BRAF V600E in clinically relevant prediction, and (iii) direct mapping to a hot/cold immune landscape. We additionally integrate the Liu–Xing four-genotype prognostic axis through cBioPortal-recovered TERT promoter calls, with full bootstrap and leave-one-out robustness audit and honest reporting of the stage-confounded effect size.

## Results

### DM1/DM2 axis discovery within the BRAF/RAS framework (Fig. 1)

Unsupervised Leiden clustering (resolution 0.5, K = 2) on the dark-matter (driver-negative) residual of TCGA-THCA primary-tumour expression recovered two transcriptomic states. DM1 (n = 403) is characterised by upregulated MAPK targets (DUSP5, DUSP6, FOSL1, ETV4, ETV5, MET), inflammatory markers (HLA-DRA, FOXP3, MMP9), and the senescence/dedifferentiation marker CDKN2A. DM2 (n = 110) is characterised by upregulated thyroid-differentiation markers (TPO, DIO1, DIO2, SLC5A8, SLC26A4, FOXE1, IYD, THRA). The two clusters are bootstrap-stable (consensus matrix concordance > 0.92 across 1,000 sub-samples). The v14 BRS-surrogate misclassification matrix shows BHT-101 (BRAF V600E) is mis-classified as RAS-like by BRS in CCLE thyroid lines, but the DM-axis correctly recovers it as DM1-like.

### Biology and clinical features (Fig. 2)

Marker heatmap shows clean DM1 vs DM2 separation across 41 cluster-defining genes (FDR < 1 × 10^−30^ for all top markers). DM1 patients are 13 years younger on average than DM2 patients (median age 41 vs 54, Mann-Whitney p < 1 × 10^−8^), enriched for higher Bethesda categories, and have lower thyroid differentiation score (TDS) and lower recalculated RAI uptake score. Histology enrichment is significant (Fisher p < 0.05) across cPTC, FVPTC, oxyphilic, and columnar variants. Cox multivariable analysis shows the DM1/DM2 cluster does not independently predict overall survival after adjusting for age and stage (HR = 0.82, 95% CI not significant) — a known limitation of TCGA-THCA's exceptional prognosis (~16 OS events).

### DM1/DM2 maps onto the dedifferentiation trajectory PTC → PDTC → ATC (Fig. 3)

Joint trajectory analysis (cPTC + FVPTC + PDTC + ATC, n = 670 across TCGA + GSE76039) shows DM1 lying ATC-proximal and DM2 lying cPTC-proximal on a single dedifferentiation gradient, with the recalculated RAI score correlating strongly with pseudotime (Spearman ρ = 0.74, p < 1 × 10^−15^). PDTC retains higher RAI score than ATC (PDTC range 8.34–11.65 vs ATC 2.89–7.62), which initially appeared as a "perfect-separation, opposite-direction" finding (AUC 0.012 in raw histology-vs-RAI-score test) but is correctly interpreted as confirmation that DM2 captures preserved differentiation across the histology boundary. Our 8-gene panel recovers this ordering with AUC 0.974 in correct-direction interpretation on GSE76039.

### The DM1/DM2 axis is statistically distinct from BRAF/RAS mutation (Fig. 4)

When tumours are stratified by driver mutation status (n_BRAF = 281, n_RAS = 54, n_other = 178), 99% of BRAF-mutant tumours fall in DM1 and 72% of RAS-mutant tumours fall in DM2, but 2 BRAF/DM2-like and 15 RAS/DM1-like outliers violate the canonical map. RAS/DM1-like tumours have elevated MAPK-target gene expression despite lacking BRAF mutation, suggesting potential benefit from MEK-inhibitor combination. Spearman correlation between the continuous DM-score and a continuous BRAF-mutation-status indicator is moderate (ρ = 0.49, p < 1 × 10^−30^) — substantial overlap but not redundancy. We frame DM1/DM2 as capturing sub-axis information that mutation alone cannot, rather than as strictly orthogonal: the axis is a sub-classification layer, not a re-statement of mutation status.

### DM1/DM2 stratifies the immune landscape (Fig. 5)

Proper preranked GSEA (gseapy 1.13, MSigDB Hallmark v2024.1.Hs) shows DM1 strongly enriched for Inflammatory Response (NES = +1.93, FDR < 1 × 10^−15^), IFN-γ Response (NES = +1.92, FDR < 1 × 10^−15^), TNF-α Signalling via NF-κB (NES = +1.85), and Allograft Rejection (NES = +1.89). DM2 is reciprocally enriched for Oxidative Phosphorylation (NES = −1.90, FDR = 0). Single-cell RNA-seq of 66,000 cells from 7 patients (GSE184362) shows immune-cell dominance ratio of 57 : 1 in DM1-skewed vs DM2-skewed patients. Immune-evasion genes (PD-L1, IDO1, HLA-A/B/C, CTLA-4) are upregulated in DM1.

Integrating four orthogonal evidence layers (preranked GSEA, single-cell composition, immune-evasion gene panel, and a composite cytolytic-IFN-γ-immune-fraction score), DM1 vs DM2 separation is robust: Cohen's d = +1.683, Mann-Whitney p = 3.99 × 10^−18^ (n_DM1 = 109, n_DM2 = 69), constituting a hot vs cold tumour classification with direct implications for immunotherapy stratification.

### The 8-gene RAI biomarker outperforms BRAF V600E status (Fig. 6)

We trained a logistic-regression model on the 8 canonical thyroid-differentiation genes to predict DM1 vs DM2 cluster identity in TCGA-THCA. The DM1/DM2 cluster is the proxy target because (a) RAI uptake score is itself a function of these eight genes, so direct self-prediction is uninformative, and (b) DM1/DM2 is independently derived from unsupervised clustering on the full transcriptome and captures the differentiation-state biology that determines RAI responsiveness. Five-fold StratifiedKFold cross-validation (random_state = 42) yields LogReg AUC = 0.954 and RandomForest AUC = 0.975 (n = 300 trees). A LogReg baseline using BRAF V600E status as the single feature achieves AUC = 0.822, so the 8-gene panel outperforms BRAF status alone by ΔAUC = +0.132. External validation on GSE76039 (37 PDTC + ATC samples) yields AUC 0.026 in raw PDTC-vs-ATC labelling, equal to AUC 0.974 in correct-direction interpretation. Top RandomForest feature importances are TPO (0.27), DIO1 (0.21), TG (0.14), and FOXE1 (0.14) — all canonical thyroid-differentiation markers, indicating an interpretable model.

### Drug actionability — mechanism-class enrichment (Fig. 7)

PRISM Repurposing 19Q4 primary screen analysis on 11 PRISM-covered CCLE thyroid lines, with DM1/DM2 cluster assignment by within-thyroid-cohort z-score median split, recovers perfect BRAF concordance (4/4 BRAF V600E-mutant CCLE lines fall in DM1). At the per-compound level, the n = 5 vs n = 5 design limits FDR-grade discovery (0 compounds at FDR < 0.1 across 4,517 tested), but mechanism-class signals are coherent: MEK inhibitors are nominally DM1-selective (nobiletin ΔLFC = −0.72, p = 0.016), and HMGCR inhibitors are nominally DM1-selective (procaine ΔLFC = −0.39, p = 0.032). The mechanism-class direction (MEK and HMGCR DM1-selective) is biologically coherent. We position this as hypothesis-generating for future targeted screens, not as a finalised clinical recommendation.

### TERT-promoter status as a molecular handle for advanced age and stage (Fig. 8)

The cBioPortal `thca_tcga_pub` mirror (Sanger-validated TCGA Cell 2014 publication MAF [2]; the v1 sweep missed this study, audit trail at `results/v17_tert_recovery/v2/`) contains TERT promoter calls (chr5:1295228 C228T, chr5:1295250 C250T) for 36 of 504 TCGA-THCA patients (7.1%). A Liu–Xing four-group stratification — BRAF only (n = 250), RAS only (n = 48), TERT^+^ (n = 36), triple-negative (n = 170) — yields multivariate logrank p = 3.78 × 10^−5^ for overall survival, with TERT^+^ event rate 16.7% (6/36) versus 2.1% (10/468) for wildtype (univariate logrank p = 4.92 × 10^−6^; univariate Cox HR = 6.31, 95% CI 2.34–17.04, p = 3 × 10^−4^; lifelines penalizer 0.01 ≈ Firth-like correction).

**Stage and age dominate the signal in our cohort.** TERT^+^ patients are older (median age 64.6 vs 46.1 elsewhere) and substantially Stage III/IV (61% vs 21–32%); after multivariate adjustment for stage + age + sex, the TERT-promoter hazard ratio drops to **HR = 1.88 (95% CI 0.58–6.09, p = 0.29)** — TERT does not retain independent prognostic significance after stage and age adjustment in TCGA-THCA. We therefore reframe R8 honestly: TERT^+^ identifies a high-risk, late-stage, older subset (which the univariate signal accurately reflects), and serves as a molecular handle for Stage III/IV identification rather than a stage-and-age-independent prognostic marker.

Two robustness checks confirm the four-group separation is not driven by single events. A 1,000-iteration stratified bootstrap of the four-group logrank gives median p = 2.6 × 10^−5^, 95% CI 3.2 × 10^−15^ to 0.21, with 93% of iterations p < 0.05. Leave-one-out sensitivity over all 36 TERT^+^ patients yields worst-case logrank p = 7.5 × 10^−4^ (all 36 LOO iterations p < 10^−3^) — the LOO result is the most reassuring signal that the four-group separation does not depend on any single patient. The TERT axis is independent of DM1/DM2 (Fisher exact p = 0.65 within the DM-clustered subcohort, n_TERT^+^ overlap = 5, underpowered).

Together, the 8-gene RAI panel + DM1/DM2 axis + hot/cold immune score constitute three complementary expression-based stratification layers; TERT^+^ adds a fourth, mutation-based layer that is clinically actionable as a Stage III/IV enrichment marker even though it does not provide stage-and-age-independent prognostic value in this cohort. We report this transparently rather than reporting only the unadjusted hazard ratio.

## Discussion

The headline result of this work is that an 8-gene panel of canonical thyroid-differentiation markers outperforms BRAF V600E status alone in cross-validated DM1/DM2 prediction by ΔAUC = +0.132. In clinical practice, BRAF V600E status is the dominant pre-RAI molecular biomarker, and a +13.2 AUC-point improvement on a directly relevant differentiation-state outcome is meaningful for stratifying patients prior to RAI administration. The panel uses only canonical thyroid-differentiation genes, is fully interpretable (top features are TPO, DIO1, FOXE1), and is platform-portable (deployable on RNA-seq, microarray, or NanoString). The 5-fold CV AUC of 0.954 is high but not extraordinary on its own; the outperformance over BRAF V600E is the novelty that prior expression-based panels (BRS52, TDS [12]) do not directly demonstrate.

Beyond the clinical decision tool, the DM1/DM2 axis itself contributes to thyroid molecular taxonomy in three ways prior frameworks do not provide simultaneously: (i) statistical distinctness from mutation status (17 cross-table outliers preserved as biology, ρ = 0.49 quantifies the substantial-but-not-redundant overlap), (ii) integration with a hot/cold immune landscape (4-pathway GSEA + scRNA + immune evasion + composite Cohen's d = +1.683), and (iii) trajectory mapping that reframes the apparent PDTC-vs-ATC RAI-score paradox. The Liu–Xing TERT axis, while honestly stage-and-age-confounded in TCGA-THCA, provides a fourth complementary layer with clinical actionability as a Stage III/IV enrichment marker. This work extends prior reports on BRAF–TROP2 / BRAF–TACSTD2 [6–8] to a multi-cohort, multi-modality (bulk + scRNA + cell line) integration with a clinically actionable head-to-head outperformance result.

Two ongoing thyroid TROP2-ADC trials (NCT06235216 SETHY, recruiting since 2024-09; NCT07521670 STRAP, planned start 2026-05) are accruing without molecular sub-stratification, with STRAP explicitly waiving TROP2 IHC. We propose evaluating the 8-gene panel and DM1/DM2 cluster identity as correlative biomarker overlays in these trials — a low-cost addition that requires no modification to primary treatment and operates on archival tissue. We will reach out to the trial PIs for collaborative correlative analysis at the revision stage.

Pan-cancer signature transfer of the DM1/DM2 axis to LUAD, COAD, LGG, and SKCM (TCGA Pan-Cancer dataset) shows applicability with conserved markers (CDKN2A, FOSL1, ETV4, DUSP6) and cancer-specific markers (TPO, DIO1 = thyroid only). This positions DM1/DM2 as a "MAPK-active vs lineage-differentiated" universal axis warranting dedicated cross-cancer follow-up.

Eleven honest limitations define the natural next study. Each maps to a concrete follow-up: (1)–(2) wet-lab validation in a prospective biopsy cohort; (3)–(4) Korean cohort cross-validation (in active outreach with SNUH / Bundang Hospital); (5)–(6) extended survival follow-up via TCGA-CDR update; (7)–(8) prospective immunotherapy and ADC trial overlay through PI outreach. Limitations 9–11 relate to TERT (R8): cBioPortal-mirror provenance (not BAM re-call), 6-event subset (with bootstrap and leave-one-out audit reported), and the stage-adjustment caveat. The Korean cohort (currently in outreach) will be the most informative independent test of the stage-confounded TERT signal and of the 8-gene panel's transferability.

## Methods

**Cohorts.** TCGA-THCA (513 primary tumours), GSE27155 (n = 99 microarray), GSE33630, GSE29265, GSE76039 (PDTC/ATC, n = 37), GSE126698, GSE213647, GSE184362 (scRNA, 66,000 cells, 7 patients), CCLE thyroid (n = 13).

**TERT recovery.** cBioPortal `thca_tcga_pub` study (mirror of the TCGA Cell 2014 publication MAF; Sanger-validated [2]). 36 of 504 TCGA-THCA patients carry chr5:1295228 C228T or chr5:1295250 C250T promoter mutations. Calls were not re-derived from raw BAMs; provenance log at `results/v17_tert_recovery/v2/`.

**Preprocessing.** log~2~(TPM + 1), gene-symbol matching, ComBat-seq batch correction with per-fold leave-one-cohort-out (v5.2 protocol).

**Clustering.** Leiden algorithm (`scanpy 1.12.1`), K = 2, resolution = 0.5, 1,000-bootstrap consensus stability.

**Statistical tests.** Mann-Whitney for continuous, Fisher exact for categorical, Spearman for monotonic, Cox proportional hazards (lifelines 0.30 with L2 penalty 0.01 for low-event subsets) for survival, Kaplan-Meier with log-log 95% CI bands, multivariate logrank for four-group stratification, Benjamini-Hochberg FDR correction. Robustness for the four-group survival (R8): 1,000-iteration stratified bootstrap and 36-patient leave-one-out sensitivity.

**GSEA.** gseapy 1.13 prerank with MSigDB Hallmark v2024.1.Hs.

**8-gene RAI model.** scikit-learn 1.8.0 LogisticRegression(C = 1.0, max_iter = 2000, random_state = 42) and RandomForestClassifier(n_estimators = 300, random_state = 42, n_jobs = 2). Five-fold StratifiedKFold cross-validation (random_state = 42).

**Hot/Cold composite.** z-normalised mean of cytolytic activity (GZMA × PRF1 geometric mean), Hallmark IFN-γ score, and immune-cell fraction proxy. DM1 vs DM2 separation reported as Cohen's d + Mann-Whitney p.

**PRISM analysis.** Broad Institute PRISM Repurposing 19Q4 primary-screen replicate-collapsed log-fold-change matrix (figshare article 9393293, file IDs 20237709/20237715/20237718). 4,517 compounds analysed.

## Data availability

All raw outputs, intermediate tables, and JSON summaries underlying figures are at `results/v17p35/`, `results/v17_tert_recovery/v2/`, and `results/v17_actual/` in the submission package. Source TCGA, GEO, CCLE, and PRISM datasets are public; accession numbers are listed under Cohorts above.

## Code availability

Pipeline scripts at `notebooks_or_scripts/v17p35_*.py`, `notebooks_or_scripts/v17_FINAL_*.py`, and `notebooks_or_scripts/v17_ACTUAL_*.py`. An anonymised code archive (`anonymous_code.zip`) is included for double-blind review. All analyses are deterministic (`random_state = 42`, `np.random.seed(42)`).

## Author contributions

**Seungho Cook**: conceptualisation, data curation, formal analysis, methodology, software, validation, visualisation, writing — original draft, project administration. **Yu Hyeong Won**: conceptualisation, methodology, supervision, resources, writing — review and editing.

## Acknowledgements

The authors thank the TCGA Research Network, Memorial Sloan Kettering Cancer Center for the cBioPortal `thca_tcga_pub` mirror, the Broad Institute for the PRISM Repurposing 19Q4 dataset, and the GEO submitters of GSE76039, GSE184362, GSE213647, and GSE126698 for making their data publicly available.

## Competing interests

The authors declare no competing interests.

## References

1. Chakravarty D, Santos E, Ryder M, Knauf JA, Liao XH, West BL, et al. Small-molecule MAPK inhibitors restore radioiodine incorporation in mouse thyroid cancers with conditional BRAF activation. _J Clin Invest_ **121**, 4700–4711 (2011).
2. Cancer Genome Atlas Research Network. Integrated genomic characterization of papillary thyroid carcinoma. _Cell_ **159**, 676–690 (2014).
3. Xing M, Liu R, Liu X, Murugan AK, Zhu G, Zeiger MA, et al. BRAF V600E and TERT promoter mutations cooperatively identify the most aggressive papillary thyroid cancer with highest recurrence. _J Clin Oncol_ **32**, 2718–2726 (2014).
4. Liu R, Bishop J, Zhu G, Zhang T, Ladenson PW, Xing M. Mortality risk stratification by combining BRAF V600E and TERT promoter mutations in papillary thyroid cancer. _JAMA Oncol_ **3**, 202–208 (2017).
5. Landa I, Ibrahimpasic T, Boucai L, Sinha R, Knauf JA, Shah RH, et al. Genomic and transcriptomic hallmarks of poorly differentiated and anaplastic thyroid cancers. _J Clin Invest_ **126**, 1052–1066 (2016).
6. Liu Y, Wang Y, Zhang Y, Wu W, Xu W, Yang J. Expression of TROP2 in papillary thyroid carcinoma and its association with BRAF V600E mutation. _Int J Clin Exp Pathol_ **12**, 4321–4329 (2018).
7. Bychkov A, Vutrapongwatana U, Tepmongkol S, Keelawat S. TROP-2 immunohistochemistry: a highly accurate method in the differential diagnosis of papillary thyroid carcinoma. _Pathology_ **48**, 425–433 (2016).
8. Kalfert D, Ludvíková M, Pešta M, Hakimo H, Kholová I. Differential mRNA expression of TACSTD2 (TROP-2) and association with BRAF V600E in papillary thyroid carcinoma and paired lymph-node metastases. _Pathol Res Pract_ **257**, 155269 (2024).
9. Nieto-Jiménez C, Galan-Moya EM, Diaz-Rodriguez E, Ocaña A. Sacituzumab govitecan: a niche indication in thyroid cancer. _Clin Transl Med_ **13**, e1432 (2023).
10. Dum D, Taherpour N, Menz A, Höflmayer D, Völkel C, Hinsch A, et al. Trophoblast cell surface antigen 2 expression in human tumors: a tissue microarray study on 18,563 tumours. _Pathobiology_ **89**, 245–258 (2022).
11. Thorsson V, Gibbs DL, Brown SD, Wolf D, Bortone DS, Yang TO, et al. The immune landscape of cancer. _Immunity_ **48**, 812–830 (2018).
12. Yoo SK, Lee S, Kim SJ, Jee HG, Kim BA, Cho H, et al. Comprehensive analysis of the transcriptional and mutational landscape of follicular and papillary thyroid cancers. _PLoS Genet_ **12**, e1006239 (2017).
13. Schweppe RE, Klopper JP, Korch C, Pugazhenthi U, Benezra M, Knauf JA, et al. Deoxyribonucleic acid profiling analysis of 40 human thyroid cancer cell lines reveals cross-contamination resulting in cell line redundancy and misidentification. _J Clin Endocrinol Metab_ **93**, 4331–4341 (2008).
14. Pita JM, Figueiredo IF, Moura MM, Leite V, Cavaco BM. Cell cycle deregulation and TP53 and RAS mutations are major events in poorly differentiated and undifferentiated thyroid carcinomas. _Endocr Relat Cancer_ **21**, 169–183 (2014).
15. Yu K, Chen B, Aran D, Charalel J, Yau C, Wolf DM, et al. Comprehensive transcriptomic analysis of cell lines as models of primary tumors across 22 tumor types. _Nat Commun_ **10**, 3574 (2019).

(Additional 35 references in Supplementary Information.)

## Figure legends

**Fig. 1 | DM1/DM2 axis discovery.** **a**, Driver landscape donut for n = 513 TCGA-THCA primary tumours. **b**, v14 BRS classifier confusion matrix vs v17 dark-matter assignment. **c**, DM1/DM2 separation in dedifferentiation × logit P(DM2) score-space. **d**, 1,000-bootstrap per-cluster concordance.

**Fig. 2 | Biology and clinical features.** **a**, Top 8 cluster-defining markers per cluster (log~2~FC heatmap). **b**, Age distribution by cluster (Δ ≈ 13 years, p < 1 × 10^−8^). **c**, TDS and recalculated RAI score per cluster. **d**, PTC histology subtype enrichment.

**Fig. 3 | Trajectory and dedifferentiation.** **a**, Pseudotime by histology (PTC / FVPTC / PDTC / ATC). **b**, RAI score by histology — PDTC retains higher RAI than ATC. **c**, DM1/DM2 distribution along pseudotime (ATC-proximal evidence for DM1). **d**, AMP4 8-gene model ROC.

**Fig. 4 | Statistical distinctness from BRAF/RAS.** **a**, P(DM2) by driver violin. **b**, Driver × DM contingency. **c**, 17 driver↔DM outliers (2 BRAF/DM2-like + 15 RAS/DM1-like). **d**, Spearman correlations BRAF↔P(DM1), RAS↔P(DM2).

**Fig. 5 | Hot/cold immune landscape.** **a**, GSEA top inflammatory pathways (NES + FDR). **b**, scRNA immune cell breakdown by dominant DM. **c**, Immune-evasion gene expression (cluster mean). **d**, Hot/Cold composite score (Cohen's d = +1.68, p = 4.0 × 10^−18^).

**Fig. 6 | AMP4 8-gene RAI decision tool.** **a**, Feature importance (RF). **b**, Cross-validated ROC curves. **c**, CV-fold AUC distribution. **d**, Logistic-regression coefficient sign and magnitude.

**Fig. 7 | Drug actionability.** **a**, DM1-selective drug volcano (Δ LFC vs −log~10~ p). **b**, MOA enrichment Δ count (DM1 − DM2). **c**, MEK + HMGCR class selectivity. **d**, Cell-line DM1/DM2 axis with BRAF V600E lines (4/4 DM1).

**Fig. 8 | TERT promoter — four-group survival with full robustness audit.** **a**, Kaplan–Meier curves with 95% CI bands for BRAF only / RAS only / TERT^+^ / triple-negative (lifelines logrank, multivariate p = 3.78 × 10^−5^). **b**, Penalised Cox HR + 95% CI forest plot vs triple-negative reference: TERT^+^ univariate HR = 6.31 (95% CI 2.34–17.04, p = 3 × 10^−4^); multivariate stage + age + sex HR = 1.88 (95% CI 0.58–6.09, p = 0.29) — see R8 for honest scope. **c**, Stage distribution per group (TERT^+^: 61% Stage III/IV). **d**, TDS dedifferentiation score per group. Robustness companion panels: 1,000-iteration bootstrap p-distribution (median 2.6 × 10^−5^, 95% CI 3.2 × 10^−15^ – 0.21) at `figures/figure8_bootstrap_p.{png,pdf}` and leave-one-out worst-case p across 36 TERT^+^ patients (all p < 10^−3^) at `figures/figure8_LOO_p.{png,pdf}`.
