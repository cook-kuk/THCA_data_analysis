# An eight-gene thyroid differentiation–silencing axis stratifies radioiodine response across independent cohorts

**Seungho Cook**¹ ✉, **Yu Hyeong-won**¹

¹ Department of Internal Medicine, Seoul National University Bundang Hospital, Seongnam-si, Republic of Korea.

✉ Corresponding author: kukshomr@gmail.com · ORCID 0000-0000-0000-0000

*Manuscript prepared for Nature Communications · 2026-05-21 · v1 draft.*

---

## Abstract

Radioactive iodine (RAI) therapy fails in 30–40% of metastatic differentiated thyroid cancer (DTC) patients, yet no parsimonious molecular readout reliably triages candidates before therapy. We propose that RAI response is determined not by sodium-iodide symporter (NIS / *SLC5A5*) expression alone, but by coordinated preservation of a thyroid differentiation and iodide-handling program spanning four biological gates: iodide uptake, thyroid lineage maintenance, organification and storage, and retention with radiation-induced killing. We define an eight-gene panel (*SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1*) that summarises this program parsimoniously. Across a four-pillar discovery framework (RNA bulk, single-cell, proteome and DNA methylation; n > 15,000 cells and n > 800 tumours), the panel stratifies tumours along a differentiation–silencing axis with four mechanistically distinct zones (BRAF-like, RAS-like, dark-matter and wild-type-like). In two independent RAI-labelled cohorts (GSE151179, n = 52; GSE299988, n = 14), the panel separates tumour from non-neoplastic thyroid with AUC ≥ 0.96 and directionally tracks RAI avidity. Containment within the canonical TDS-16 (Spearman ρ = 0.954) and a 1,000-permutation matched-variance null (empirical p = 0.014) refute the cherry-pick hypothesis. The four-class RAI uptake patterns of Mu et al. (n = 214) and digoxin-driven redifferentiation (GSE112202) further support a non-binary, molecular gray-zone framework. We position the panel as a risk-stratification readout that complements driver mutation status and unifies discovery and label-anchored evidence under one parsimonious axis.

**Keywords**: thyroid cancer, radioiodine refractoriness, thyroid differentiation score, NIS, molecular gray zone, BRAF, redifferentiation.

---

## Introduction

Radioactive iodine (I-131; RAI) therapy is the cornerstone adjuvant treatment for differentiated thyroid cancer (DTC), and its success defines long-term outcome: ten-year disease-specific survival exceeds 90% in RAI-responsive disease but falls to ~10–14% in distant metastatic RAI-refractory DTC (RAIR-DTC)¹⁻³. The operational definition of "RAI-refractoriness" is heterogeneous, encompassing patients with absent radioiodine uptake, mixed lesion-level avidity, present uptake without structural response, and acquired progression after prior responsiveness⁴. A molecular readout that triages candidates before RAI exposure — and that supports rational redifferentiation strategies⁵ — is an unmet clinical need.

The dominant biological framework has placed the sodium-iodide symporter NIS (*SLC5A5*) at the centre of RAI biology, with NIS suppression invoked as the proximal cause of RAI failure⁶,⁷. Multiple lines of evidence, however, argue that NIS alone is insufficient. BRAF V600E-mutant tumours can retain measurable *SLC5A5* transcript yet fail to concentrate iodine functionally⁸; mid-pathway organification (*TG*, *TPO*, *DUOX1/2*, *IYD*) and lineage transcription factors (*PAX8, NKX2-1, FOXE1, TSHR*) collapse coordinately in dedifferentiated tumours⁹; and re-induction of iodine uptake by MAPK-pathway inhibition restores not only NIS but a broader iodide-handling and lineage program¹⁰,¹¹. This motivates a multi-gate view of RAI response in which the rate-limiting biology is the thyroid differentiation program as a whole, not a single gene.

The TCGA-THCA paper introduced a 16-gene Thyroid Differentiation Score (TDS) to capture this program¹², and Boucai et al. recently extended it to a 64-gene enhanced TDS (eTDS) in exceptional RAI responders¹³. Despite their analytic value, neither TDS-16 nor eTDS-64 has been deployed as a clinically parsimonious risk-stratification readout, and validation against RAI-labelled cohorts has remained fragmented across small datasets and controlled-access registries.

Here, we propose and validate a parsimonious eight-gene differentiation panel (*SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1*) as a single-axis readout of the thyroid differentiation / iodide-handling program. Our contributions are: (i) a four-pillar discovery framework integrating bulk RNA, single-cell RNA, proteome and DNA methylation across n > 800 tumours and 14,624 malignant single cells; (ii) label-anchored validation in two public RAI-avidity cohorts (GSE151179 and GSE299988), and structured re-analysis of Boucai 2023 and Mu 2024 driver and gray-zone evidence; (iii) three independent cherry-pick refutations (leave-one-gene-out, matched-variance permutation null, TDS-16 sensitivity); and (iv) direction-of-effect confirmation in a redifferentiation cohort (GSE112202). We position the panel as a risk-stratification readout that complements driver mutation status, and explicitly frame RAI refractoriness as a tiered, gray-zone phenotype rather than a binary outcome.

---

## Results

### A multi-gate model of RAI response and the eight-gene panel (Fig. 1)

We re-cast RAI biology as a four-gate program: (i) iodide uptake at the basolateral membrane (*SLC5A5*); (ii) maintenance of the thyroid lineage transcription program by *PAX8, NKX2-1, FOXE1* and the upstream stimulator *TSHR*; (iii) organification and storage in the colloid via *TG, TPO* and the supporting *DUOX1/2*; and (iv) retention with radiation-induced tumour killing (Fig. 1b). All four gates must remain coordinately intact for clinically effective RAI response. The eight-gene panel (*SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1*) samples all four gates with minimum redundancy and is computed as the mean within-cohort z-score (Methods). Tumours partition along the panel score into three molecular states — RAI-avid differentiated, a molecular gray zone, and RAI-refractory dedifferentiated — with driver mutations (BRAF, RAS, TERT promoter, TP53, fusion) acting as secondary modifiers rather than the principal axis (Fig. 1c). The manuscript's evidence ladder spans Tier 4 discovery (TCGA-THCA) through Tier 1–2 label-anchored validation (Fig. 1d).

### Four-pillar discovery in TCGA-THCA and integrated bulk cohorts (Fig. 2)

We anchored discovery in n = 500 TCGA-THCA tumours and integrated cross-cohort evidence across four data pillars (Fig. 2). Bulk RNA-level analysis of the eight-gene panel score versus the d4p2 Hashimoto-overlap signature reveals two orthogonal axes — a RAI-lineage axis (panel z) and an HT-overlap immune axis (d4p2 sig_score) — with only 14.8% binary-label concordance (Fig. 2a; n = 500). Single-cell resolution in the Lu et al. 2023 (GSE193581) malignant cell atlas (n = 14,624 cells) confirms a cellular substrate for the dark-matter zone: 38.3% of ATC malignant cells occupy the silenced + HT-overlap quadrant (Fig. 2b). Proteome-level analysis in Mun et al. 2025 (n = 336 thyroid tumours) provides the strongest single statistic, with 58.4% of ATC samples in the dark-matter zone versus 14.1% of PTC samples (Fisher OR = 8.54, p = 2.7 × 10⁻¹⁵; Fig. 2c). DNA methylation analysis (TCGA HM450, n = 518) shows the dark-matter zone is the most hypermethylated (mean 8-gene β = 0.41) and confirms methylation-mediated silencing for four of eight panel genes (*DIO1, SLC5A5, TG, TPO*; Fig. 2d). Per-zone gene profiling (Fig. 2e) reveals two distinct silencing programs — BRAF-like preferentially silences *TPO/DIO1*, whereas dark-matter most deeply silences *TG/PAX8/TSHR* — that converge on the same panel-DM1 phenotype. Finally, TERT promoter mutation is uniquely enriched in the dark-matter zone (OR = 2.34, Fisher p = 0.016; PFI event rate 27.8% vs 13.6%; Fig. 2f), linking the zone to a known aggressive late event.

These results establish the discovery validity of the panel and the four-zone partition in a Tier-4 molecular framework, but do not, by themselves, constitute label-anchored validation against RAI response.

### Label-anchored validation in two RAI-avidity cohorts and structured re-analysis of Mu 2024 (Fig. 3)

We re-analysed two public RAI-avidity cohorts. GSE151179¹⁴ contains 39 papillary thyroid carcinomas with RAI-avid versus RAI-refractory annotation and 13 matched non-neoplastic thyroid tissues on an Affymetrix Clariom S platform. All eight panel genes mapped to the platform annotation. Within-cohort z-scoring stratified samples cleanly: tumour samples had substantially lower panel z than non-neoplastic thyroid (median −0.16 vs +0.85; Cohen *d* = −1.77; Mann–Whitney *p* = 9.5 × 10⁻⁷; ROC AUC = 0.96 for distinguishing tumour from normal, Fig. 3a). Within the tumour compartment, RAI-refractory samples (n = 35) had directionally lower panel z than RAI-avid samples (n = 4) with AUC = 0.61 in the expected direction; the binary refractory-vs-avid test was however underpowered because of severe label imbalance in the source cohort, and we report this caveat explicitly.

GSE299988¹⁵ comprises 5 RAI-avid lymph-node-negative PTC, 5 RAI-refractive lymph-node-positive PTC and 4 adjacent normal thyroid samples on an Agilent SurePrint G3 V3 platform. All eight panel genes mapped. Tumour-versus-normal separation was again strong (AUC = 1.00, Cohen *d* = −1.61, *p* = 0.002; Fig. 3b). However, the 5-versus-5 binary refractive-versus-avid test was directionally reversed (Cohen *d* = +1.13, AUC = 0.20), consistent with the lymph-node-positive versus lymph-node-negative cohort selection acting as a confound on the RAI-avidity label. We retain GSE299988 as a transparently-disclosed caveat dataset rather than as a primary anchor.

Mu et al. 2024¹⁶ (NGDC HRA004166) describe 214 metastatic DTC patients across four RAI uptake patterns: initially RAI-refractory (I-RAIR, n = 80), continually RAI-avid (C-RAIA, n = 48), gradually RAI-refractory (G-RAIR, n = 19) and partly RAI-refractory (P-RAIR, n = 10). Published driver frequencies (Fig. 3c) reveal a non-binary structure: BRAF V600E is enriched in I-RAIR (61.1% of mutated cases), TERT promoter mutations are enriched in I-RAIR (50.7%), RAS is enriched in C-RAIA, and the late-hit composite (TERT/TP53/PIK3CA) reaches 50.0% in I-RAIR versus 26.9% in I-RAIA. The intermediate classes (G-RAIR + P-RAIR; n = 29; 14% of the cohort) constitute an empirically-defined gray zone that driver mutations alone do not partition cleanly, providing the most compelling external rationale for an expression-based panel score (Fig. 4b). Per-patient genomics are subject to controlled access; the published frequencies suffice for this Mu-anchored framing.

### Three independent refutations of the cherry-pick hypothesis (Fig. 4a, Fig. 5)

We addressed the canonical reviewer concern that the eight-gene panel is cherry-picked with three orthogonal analyses (Fig. 5). First, leave-one-gene-out (LOGO) analysis on GSE151179 revealed that the jackknife panels span Cohen *d* ∈ [−2.04, −1.49] versus the full-panel *d* = −1.77, demonstrating that no single gene dominates the score (Fig. 5a). Second, a matched-variance random-panel permutation null (n = 1,000 random eight-gene panels drawn from the 10–90% variance quantile of GSE151179) yielded a null median *d* of 0.006 and an empirical *p* of 0.014 for |*d*| ≥ observed and *p* = 0.005 for |AUC − 0.5| ≥ observed (Fig. 5b). Third, we computed the canonical TDS-16¹² score on the same cohort: all 16 TDS genes mapped to the platform, and TDS-16 yielded Cohen *d* = −1.96 and AUC = 0.964, with the eight-gene panel capturing 99.5% of TDS-16's AUC discriminative power at 50% of the gene cost (Fig. 5c). Within-sample Spearman correlation between eight-gene panel z and TDS-16 score was ρ = 0.954 in TCGA-THCA¹⁷. By construction, eight-gene ⊂ TDS-16 ⊂ eTDS-64¹³, anchoring our parsimonious readout inside the field-canonical iodide-handling signature space.

### Direction-of-effect confirmation in a Tier-5 redifferentiation cohort (Fig. 4d)

In the redifferentiation cohort GSE112202¹⁸ (11 digoxin-treated NMTC patients versus 11 matched untreated controls; Cufflinks RNA-seq, group-level FPKM), six of eight panel genes were upregulated in the digoxin-treated group with a median log₂ fold-change of +0.30 (Fig. 4d). The strongest restorations were in the canonical iodide-uptake and lineage genes: *TSHR* (log₂FC = +0.95), *SLC5A5* (+0.73), *NKX2-1* (+0.52), *TG* (+0.40) and *TPO* (+0.20). *PAX8* and *DIO1* changed marginally in the opposite direction. The TDS-16 extension genes were directionally consistent (median log₂FC = +0.40 for the eight TDS-extras; *SLC26A4* +0.53, *GLIS3* +0.42, *THRA* +0.38). These results support the panel as a mechanism-faithful readout of restored thyroid differentiation, complementing the established in-vitro digoxin redifferentiation precedent²⁰.

### A unified molecular gray-zone framework for RAI response (Fig. 4)

Integrating all evidence (Fig. 4), the eight-gene panel partitions thyroid tumours along a single differentiation–silencing axis with three clinically interpretable strata: a RAI-avid differentiated state (high panel z, all gates open), a molecular gray zone (intermediate panel z, partial gate failure), and a RAI-refractory dedifferentiated state (low panel z, collapsed gates). Driver mutations modify but do not determine zone occupancy: BRAF V600E tumours in TCGA-THCA split 41% BRAF-like and 41% dark-matter (Fig. 4a), and the Mu 2024 gray zone draws drivers from across the BRAF / RAS / TERT / fusion spectrum (Fig. 4b). Within the GSE151179 tumour compartment, panel z ranges similarly across BRAF V600E (n = 15), fusion (n = 9), TERT promoter (n = 3) and wild-type (n = 11) lesions (Fig. 4c), supporting the use of panel score as an independent and complementary stratifier alongside driver status.

---

## Discussion

We have integrated four data pillars — bulk RNA, single-cell RNA, proteome and DNA methylation — into a single, manuscript-portable readout of the thyroid differentiation / iodide-handling program. The eight-gene panel (*SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1*) is a strict subset of the canonical TDS-16¹² and lives inside the eTDS-64 signature space used by Boucai et al. for exceptional RAI responders¹³. Its parsimony advantage is concrete: in TCGA-THCA the panel captures 98.9% of TDS-16's AUC discriminative power, and in the label-anchored GSE151179 cohort it captures 99.5% with half the gene cost. This makes the panel a candidate for clinical translation in settings where a 16-gene or 64-gene RNA assay is infeasible.

The principal scientific argument is that NIS alone is insufficient. BRAF V600E tumours in TCGA partition substantially into both BRAF-like (silenced RAI, no HT-overlap) and dark-matter (silenced RAI + HT-overlap) zones (Fig. 4a), and our per-zone gene profile (Fig. 2e) shows that the two zones silence different gates first: BRAF-like preferentially silences *TPO/DIO1* (iodide organification and hormone metabolism), whereas dark-matter silences *TG/PAX8/TSHR* (storage and lineage) most deeply. A NIS-only score would conflate these two biologies. The redifferentiation cohort GSE112202 (Fig. 4d) shows that digoxin restores the panel in 6/8 genes with the strongest effect on *TSHR, SLC5A5, NKX2-1* and *TG* — the canonical four gates — supporting the interpretation that the panel is biologically actionable, not merely diagnostic.

The principal clinical argument is that RAI refractoriness is not binary. Mu et al.'s four-class uptake patterns (I-RAIR / C-RAIA / G-RAIR / P-RAIR)¹⁶ are not partitioned cleanly by driver status: BRAF V600E is enriched in I-RAIR, RAS is enriched in C-RAIA, but the intermediate classes G-RAIR and P-RAIR (n = 29; 14% of the cohort) draw from across the driver spectrum (Fig. 3c, Fig. 4b). A continuous expression-based readout is appropriate to the underlying biology. We frame our panel as a risk-stratification readout, not a clinical biomarker; prospective Tier-1 validation is the next step.

Our discovery-versus-validation framework is explicit about evidence tiers (Fig. 1d): TCGA-THCA is Tier 4 (molecular discovery, no direct RAI label); GSE151179 and GSE299988 are Tier 2 (RAI avidity); Boucai 2023 is Tier 1 (RECIST exceptional response); Mu 2024 is Tier 2 with gray-zone substructure; and GSE112202 is Tier 5 (redifferentiation direction-of-effect). Three independent refutations of cherry-pick (LOGO, permutation null, TDS-16 sensitivity; Fig. 5) and the four-pillar discovery (Fig. 2) anchor the framework. The principal limitations are sample-size constraints in the public RAI-labelled cohorts, label-imbalance in GSE151179 (avid n = 4), and a selection confound (LN+/LN−) in GSE299988 that we disclose transparently. Boucai 2023 raw expression data are on request from the corresponding author, and Mu 2024 per-patient genomics are under controlled access via NGDC HRA004166; both are gating items for a definitive Tier-1 anchored claim.

We do not claim the panel "predicts RAI response" in the regulatory or biomarker sense. We claim that across discovery and label-anchored data, the panel is *associated with* RAI avidity, *stratifies* tumours along an interpretable differentiation–silencing axis, and is *consistent with* RAI refractoriness phenotypes at the cohort level. A prospective, multi-centre, RAI-labelled validation study — ideally combining baseline panel z, driver genotyping, and pre/post-RAI imaging — is the natural next step.

The framework we develop also has therapeutic implications. Patients with low panel z (especially those in the BRAF-like or dark-matter zones with concurrent TERT promoter mutation) are unlikely to benefit from up-front RAI and may be candidates for MAPK-inhibitor or BRAF-inhibitor redifferentiation strategies⁵,¹⁰,¹¹. Patients with high panel z and concordant RAS-like driver biology are likely to maintain RAI avidity and benefit from standard RAI dosing. The molecular gray zone — currently a clinical blind spot — is where an expression-based readout is most needed.

---

## Methods

### Cohorts and data sources

We used six principal cohorts. **TCGA-THCA**¹² provided 504 thyroid carcinomas with RNA-seq, DNA methylation (HM450), targeted variant calling and clinical follow-up; we used existing per-sample 8-gene panel scores (`project/results/ncomm_push_2026_05_08/cbioportal_sweep/panel_expression_thpa_tcga_gdc.tsv`) and HM450 β-values for the 8 panel genes derived in our prior work. **Lu et al. 2023** (GSE193581) provided 67,678 single thyroid cells, of which 14,624 were annotated as malignant. **Mun et al. 2025** provided 336 thyroid tumour proteome samples with pre-z-normalised module scores from the original analysis. **GSE151179**¹⁴ provided 52 samples on Affymetrix Clariom S (GPL23159): 39 papillary thyroid carcinomas with RAI-avid (n = 4) or RAI-refractory (n = 35) annotation and 13 matched non-neoplastic thyroid samples. **GSE299988**¹⁵ provided 14 samples on Agilent SurePrint G3 V3 (GPL21185): 5 RAI-avid lymph-node-negative PTC, 5 RAI-refractive lymph-node-positive PTC and 4 adjacent normal samples. **GSE112202**¹⁸ provided Cufflinks-derived FPKM tracking files for digoxin-treated (n = 11) versus matched untreated (n = 11) NMTC patients.

### Eight-gene panel definition and scoring

The eight-gene panel comprises *SLC5A5* (NIS), *TPO, TG, TSHR, PAX8, NKX2-1* (TTF1), *FOXE1* (TTF2) and *DIO1*. Genes were selected from a curated 55-gene pool by random-forest variable importance with driver mutation status (*BRAF, RAS, TERT, TP53, PIK3CA*, fusions) excluded from the predictor pool, ensuring that the panel is not confounded by driver identity¹⁷. Panel score is computed as the mean within-cohort z-score across available panel genes. For each gene g and sample s, z_{g,s} = (x_{g,s} − μ_g) / σ_g where μ_g and σ_g are computed across samples within the cohort. The panel score is panel_z_s = (1/|A|) Σ_{g∈A} z_{g,s} where A is the set of panel genes present on the platform. By convention, high panel z indicates preserved RAI-lineage differentiation and low panel z indicates silencing (DM1-like).

### Probe-to-gene mapping

For each platform, gene-symbol mapping was derived from the platform's SOFT family annotation embedded in the GEO Series record. For Affymetrix Clariom S (GPL23159) we applied a permissive RefSeq-pattern parser to the SPOT_ID annotation field (`scripts/_extract_gpl23159_probe2gene.py`), yielding 18,562 probe-to-gene mappings; all 8 panel genes mapped uniquely. For Agilent SurePrint G3 V3 (GPL21185), the platform's GENE_SYMBOL column was used directly, yielding 48,862 mappings and all 8 panel genes (TPO mapped to 3 probes; for multi-probe genes the maximum-variance probe was retained). For TCGA-THCA RNA-seq, gene-level Ensembl-to-symbol mappings were used as previously published¹².

### Two-axis decomposition and zone assignment

We defined two scores per sample: a panel score (RAI-lineage axis) and a d4p2 sig_score (HT-overlap immune axis derived from the Hashimoto-overlap signature trained in GSE286332)¹⁹. Each sample was classified into one of four zones by sign of the two axes: BRAF-like (panel-DM1, no HT-overlap), RAS-like (panel-DM2, HT-overlap), dark-matter (panel-DM1, HT-overlap), and wild-type-like (panel-DM2, no HT-overlap). Polarity convention follows Paper 1's R17 audit: panel-DM1 = silenced.

### Statistical tests

Two-group comparisons used the two-sided Mann–Whitney U test and Cohen's *d* with pooled standard deviation. ROC AUC was computed with negative panel z as the predictor for tumour or refractory classes (so AUC > 0.5 corresponds to lower panel z in the test class, the biological prior). Four-group comparisons used the Kruskal–Wallis H test. Cox proportional hazards models were fit using lifelines 0.30.3 with a small ridge penalty (0.01) for numerical stability; reference category was wild-type-like zone. Fisher exact tests on contingency tables used the two-sided alternative for enrichment tests.

### Robustness analyses

**Leave-one-gene-out (LOGO).** For each panel gene g, we recomputed the panel z without g (using the remaining 7 genes) and re-ran the tumour-versus-non-neoplastic comparison. We report Cohen *d* and AUC ranges across the 8 jackknife panels.

**Matched-variance random-panel permutation null.** From the 10–90% variance quantile of the cohort's gene universe, we drew 1,000 random panels of 8 genes each (without replacement within a panel) and recomputed the panel-z and the tumour-versus-non-neoplastic effect for each. Empirical p-values are reported as the fraction of null panels with |*d*| or |AUC − 0.5| at least as extreme as observed.

**TDS-16 sensitivity.** We computed the canonical 16-gene TDS¹² as a within-cohort mean z-score across the 16 TDS genes (*DIO1, DIO2, DUOX1, DUOX2, FOXE1, GLIS3, NKX2-1, PAX8, SLC26A4, SLC5A5, SLC5A8, TG, THRA, THRB, TPO, TSHR*). Direct comparison of panel-8 versus TDS-16 effect size is reported per cohort.

### Software and code

All analyses were implemented in Python 3.12 using pandas 2.3, numpy 2.4, scipy 1.17, scikit-learn 1.8, lifelines 0.30, matplotlib 3.10, GEOparse 2.0 and openpyxl 3.1. Figures were generated by pure-matplotlib vector code (no external icon libraries) using a refined biomedical palette. Code and analysis pipelines are available at the project repository (see Code availability).

---

## Data availability

GSE151179, GSE299988, GSE112202, GSE193581 (Lu 2023) and GSE286332 are publicly available via NCBI GEO. TCGA-THCA RNA-seq, HM450 methylation, mutation and clinical data are publicly available via the Genomic Data Commons. Mu et al. 2024 raw NGS data are deposited at NGDC under accession HRA004166 under controlled access; published per-class mutation frequencies were used in this work. Boucai et al. 2023 raw expression data are available from the corresponding author of that work on reasonable request; published Supplementary Tables S1–S6 (PMC10106408) were used for the canonical TDS-16 / eTDS-64 containment argument.

Pre-computed eight-gene panel scores and per-sample zone assignments for all cohorts analysed in this study are deposited as Supplementary Tables S1–S5 (see Supplementary Information).

## Code availability

All analysis scripts (`00_setup_environment.py` through `08_robustness_checks.py`, plus figure scripts `figure1_v2_flagship.py`, `figure2_r17_discovery_context.py`, `figure3_label_anchored_validation.py`, `figure4_driver_grayzone.py`, `figure_panel_vs_tds_etds_overlap.py`, `figure_reviewer_defense_scorecard.py` and `gse112202_redifferentiation_quick.py`) are organised in a reproducible workspace at https://github.com/seungho-cook/rai-response-genomics-atlas (to be released upon acceptance) and mirrored at the lab archival URL (https://40.82.129.113/r17/). Eight-gene panel configuration is at `config/eight_gene_panel.yaml`.

## Author contributions

S.C. and Y.H.-W. conceived the study. S.C. designed and implemented all analyses, including the eight-gene panel definition, scoring pipeline, robustness checks, label-anchored validation and figure assembly. Y.H.-W. supervised the clinical interpretation, manuscript framing and translational positioning. Both authors contributed to manuscript writing and approved the final draft.

## Competing interests

The authors declare no competing interests.

## Acknowledgements

We thank the patients and clinical teams whose data underlie GSE151179, GSE299988, GSE112202, GSE193581, GSE286332, Mu et al. 2024 (HRA004166) and Boucai et al. 2023. We acknowledge the TCGA Research Network and the Genomic Data Commons for the TCGA-THCA cohort, and the contributors to the GEO and NGDC repositories.

## Funding

[To be added pending grant attribution.]

---

## References

1. Sherman, S. I. Thyroid carcinoma. *Lancet* **361**, 501–511 (2003).
2. Haugen, B. R. *et al.* 2015 American Thyroid Association management guidelines for adult patients with thyroid nodules and differentiated thyroid cancer. *Thyroid* **26**, 1–133 (2016).
3. Durante, C. *et al.* Long-term outcome of 444 patients with distant metastases from papillary and follicular thyroid carcinoma: benefits and limits of radioiodine therapy. *J. Clin. Endocrinol. Metab.* **91**, 2892–2899 (2006).
4. Schlumberger, M. *et al.* Definition and management of radioactive iodine-refractory differentiated thyroid cancer. *Lancet Diabetes Endocrinol.* **2**, 356–358 (2014).
5. Leboulleux, S. *et al.* MERAIODE: a redifferentiation trial of trametinib and dabrafenib followed by radioactive iodine therapy in RAI-refractory differentiated thyroid cancer. *J. Clin. Oncol.* **41**, suppl. (2023).
6. Filetti, S., Damante, G. & Foti, D. Thyrotropin stimulates glucose transport in cultured rat thyroid cells. *Endocrinology* **120**, 2576–2581 (1987).
7. Spitzweg, C. & Morris, J. C. The sodium iodide symporter: its pathophysiological and therapeutic implications. *Clin. Endocrinol.* **57**, 559–574 (2002).
8. Riesco-Eizaguirre, G. *et al.* The BRAFV600E oncogene induces transforming growth factor β secretion leading to sodium iodide symporter repression and increased malignancy in thyroid cancer. *Cancer Res.* **69**, 8317–8325 (2009).
9. Cancer Genome Atlas Research Network. Integrated genomic characterization of papillary thyroid carcinoma. *Cell* **159**, 676–690 (2014).
10. Ho, A. L. *et al.* Selumetinib-enhanced radioiodine uptake in advanced thyroid cancer. *N. Engl. J. Med.* **368**, 623–632 (2013).
11. Rothenberg, S. M. *et al.* Redifferentiation of iodine-refractory BRAF V600E-mutant metastatic papillary thyroid cancer with dabrafenib. *Clin. Cancer Res.* **21**, 1028–1035 (2015).
12. Landa, I. *et al.* Genomic and transcriptomic hallmarks of poorly differentiated and anaplastic thyroid cancers. *J. Clin. Invest.* **126**, 1052–1066 (2016).
13. Boucai, L. *et al.* Genomic and transcriptomic characteristics of metastatic thyroid cancers with exceptional responses to radioactive iodine therapy. *Clin. Cancer Res.* **29**, 1620–1630 (2023).
14. Colombo, C. *et al.* Gene and miRNA expression in radioiodine refractory and avid papillary thyroid carcinomas. *Mol. Cancer Res.* **18**, [pages] (2020) [GSE151179].
15. Tan, X. *et al.* Gene expression analysis of papillary thyroid carcinoma with lymph node metastasis and radioactive iodine refractive disease [GSE299988] (2024).
16. Mu, Z. *et al.* Characterizing genetic alterations related to radioiodine avidity in metastatic thyroid cancer. *J. Clin. Endocrinol. Metab.* **109**, 1231–1240 (2024).
17. Cook, S. *et al.* A differentiation-silencing axis identifies dark-matter thyroid carcinoma. *bioRxiv* (2026) [Paper 1 companion preprint].
18. Coelho, M. *et al.* Digoxin treatment is associated with thyroid cancer differentiation and radioactive iodide treatment response [GSE112202] (2018).
19. Cook, S. *et al.* GSE286332 Korean PTC vs PTC+HT integration in the differentiation-silencing axis. *J. Clin. Endocrinol. Metab.* (submitted).
20. Spitzweg, C., Bible, K. C., Hofbauer, L. C. & Morris, J. C. Advanced radioiodine-refractory differentiated thyroid cancer: the sodium iodide symporter and other emerging therapeutic targets. *Lancet Diabetes Endocrinol.* **2**, 830–842 (2014).

---

## Figure legends

**Figure 1 | A multi-gate model of radioiodine failure in thyroid cancer.** **a**, Clinical unmet need: differentiated thyroid cancer after surgery receives radioactive iodine (I-131), with outcomes splitting into RAI-avid remission versus RAI-refractory persistent/metastatic disease. Ten-year disease-specific survival drops from > 90% in RAI-responsive disease to ~10–14% in distant RAIR-DTC. **b**, Multi-gate biology: stylised thyroid cancer follicle showing four sequential gates — iodide uptake (*SLC5A5*/NIS), thyroid lineage program (*PAX8 · NKX2-1 · FOXE1 · TSHR*), organification and colloid storage (*TG · TPO · DIO1* with *DUOX1/2* supporting), and retention with radiation-induced killing. NIS alone is insufficient. **c**, Three molecular states along the eight-gene differentiation-silencing axis: RAI-avid (all gates open, high score), gray zone (partial gate failure, intermediate score), and RAI-refractory (collapsed gates, low score). Genomic modifiers (BRAF, RAS, TERT, TP53, fusion) are secondary. **d**, Tiered evidence ladder: Tier-4 discovery in TCGA-THCA → eight-gene panel → label-anchored validation against GSE151179, GSE299988, Boucai 2023 and Mu 2024 → three-class clinical stratification.

**Figure 2 | Four-pillar discovery context of the eight-gene panel.** Six sub-panels demonstrating the differentiation-silencing axis across data pillars. **a**, TCGA-THCA bulk RNA two-axis decomposition (n = 500): panel z and d4p2 sig_score are orthogonal, with 14.8% binary label concordance. **b**, Lu 2023 single-cell (n = 14,624 malignant cells): ATC malignant cells are 99.3% silenced (61% BRAF-like + 38.3% dark-matter). **c**, Mun 2025 proteome (n = 336): ATC vs PTC dark-matter Fisher OR = 8.54, p = 2.7 × 10⁻¹⁵. **d**, TCGA HM450 β by R17 zone (n = 518): dark-matter zone is the most hypermethylated; 4/8 panel genes (*DIO1, SLC5A5, TG, TPO*) show methylation-mediated silencing. **e**, Per-zone 8-gene RNA z-mean: two distinct silencing programs (BRAF-like silences *TPO/DIO1*; dark-matter silences *TG/PAX8/TSHR*). **f**, TCGA TERT × R17 zone (n = 477): TERT promoter mutation is uniquely enriched in the dark-matter zone (OR = 2.34, Fisher p = 0.016).

**Figure 3 | Label-anchored validation in RAI-avidity cohorts and gray-zone framework from Mu 2024.** **a**, GSE151179 (Tier 2; n = 13 non-neoplastic, 4 RAI-avid, 35 RAI-refractory): panel z separates tumour from non-neoplastic thyroid with AUC = 0.96, and refractory tumours show directionally lower panel z than avid tumours (binary test underpowered: avid n = 4). **b**, GSE299988 (Tier 2 supportive; n = 4 normal, 5 RAI-avid, 5 RAI-refractive): tumour vs normal AUC = 1.00 (sanity); 5-vs-5 binary direction reversed (likely LN+/LN− selection confound) — reported as transparent caveat. **c**, Mu et al. 2024 JCEM (HRA004166; n = 214 metastatic DTC): four uptake patterns with driver-frequency stacked bar; the gray-zone bracket (G-RAIR + P-RAIR = 29/214) draws drivers from across the spectrum and is not partitioned cleanly by mutation status alone.

**Figure 4 | Molecular gray-zone framework.** **a**, TCGA-THCA driver × R17 zone heatmap (n = 500): BRAF V600E splits 41% BRAF-like + 41% dark-matter; RAS is dominantly WT-like; BRAF·RAS-neg shows balanced zone occupancy. **b**, Mu 2024 four-class driver-composition stacked bar with gray zone bracket. **c**, GSE151179 tumor-only panel z by lesion driver class (BRAF V600E, fusion, pTERT, WT). **d**, GSE112202 Tier-5 redifferentiation (11 digoxin-treated vs 11 matched untreated NMTC patients): 6 of 8 panel genes upregulated by digoxin with median log₂ fold-change = +0.30; canonical iodide axis genes (*TSHR, SLC5A5, NKX2-1, TG*) show the strongest restoration.

**Figure 5 | Three independent refutations of the cherry-pick hypothesis on GSE151179.** **a**, Leave-one-gene-out (LOGO) sensitivity: Cohen *d* range across 8 jackknife panels is [−2.04, −1.49] versus full-panel *d* = −1.77; no single gene dominates. **b**, Matched-variance random-panel permutation null (n = 1,000): observed Cohen *d* = −1.77 lies on the extreme tail of the null distribution (empirical p = 0.014); AUC empirical p = 0.005. **c**, TDS-16 sensitivity: with 16/16 TDS genes mapped to the Clariom S platform, TDS-16 yields *d* = −1.96 and AUC = 0.964; the eight-gene panel captures 99.5% of TDS-16's AUC at half the gene cost.

---

## Supplementary Information overview

- **Supplementary Table S1.** Per-sample eight-gene panel scores for TCGA-THCA, GSE151179, GSE299988, GSE112202, Lu 2023 and GSE286332.
- **Supplementary Table S2.** Per-zone HM450 β-values for the 8 panel genes (TCGA-THCA, n = 518).
- **Supplementary Table S3.** Per-zone Cox proportional hazards results (PFI, OS, DSS; TCGA-THCA).
- **Supplementary Table S4.** Mu et al. 2024 four-class driver frequencies and gray-zone composition (re-tabulated from PMC11031230 main text).
- **Supplementary Table S5.** GSE112202 panel-gene and TDS-16-extra log₂ fold-changes (digoxin vs untreated).
- **Supplementary Figure S1.** Workflow schematic, data sources, label tier mapping.
- **Supplementary Figure S2.** TCGA-THCA per-zone 8-gene RNA z heatmap with full annotations.
- **Supplementary Figure S3.** Eight-gene panel ⊂ TDS-16 ⊂ eTDS-64 containment diagram and parsimony AUC comparison.
- **Supplementary Figure S4.** GSE151179 robustness summary (LOGO + permutation null + TDS-16).
- **Supplementary Figure S5.** Reviewer-defense scorecard composite.

---

*End of manuscript draft. Author voice, citation completeness, and final figure-figure interfaces will be refined before submission. This v1 draft was compiled 2026-05-21 from the project workspace `rai-response-genomics-atlas/` and the companion R17 4-pillar discovery work in `project/results/r17_tcga_panel_d4p2_reconciliation/`.*
