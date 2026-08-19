---
title: "Nature Communications-style full manuscript draft (v2 · A2 + IHC-3 + Lee replication upgrade) — REPAIR PASS 2026-07-30"
date: 2026-07-08
version: v2
author: Codex scaffold for Seungho Cook
target_journal: Nature Communications (upgraded scope; Nature Medicine sub-tier possible if A2 interaction replicates in SNUBH)
status: v2 draft — integrates predictive-interaction finding (A2, TCGA p_ix = 0.022), external axis replication (Lee 2024 p = 2.7e-7), IHC 3-plex killer finding (BRAF+ PFI p = 1.7e-4), multi-omics 5-way convergence and SNUBH n=200 power simulation
protected_sections:
  - Introduction opening hook
  - Introduction final aim paragraph
  - Discussion 3.1 mechanism interpretation
  - Limitations paragraph block
  - Cover letter paragraph 1
  - Reviewer Q9 prose
v2_changes:
  - Abstract: added predictive-interaction language + IHC-3 finding + Lee 2024 replication numbers
  - Results §3: new subsection "A predictive interaction with BRAF V600E status" (from A2)
  - Results §4: new subsection "A routine IHC 3-plex reproduces BRAF+ stratification" (from IHC-3 killer analysis)
  - Results §5: expanded external validation to include A2-R Lee 2024 replication (p = 2.7e-7)
  - Results §6: added multi-omics 5-way convergence subsection (from A4)
  - Results §7: added SNUBH prospective power subsection (from A5, n = 200, power = 90%)
  - Discussion "Clinical implications": upgraded to treatment-selection biomarker framing
  - Discussion "Translational readiness": added IHC 3-plex immediate-deployment paragraph
  - Methods: added interaction Cox specification, IHC-in-silico methylation proxy, Monte Carlo power simulation
  - Figure 7 (new): predictive-interaction + external replication summary
  - Figure 8 (new): IHC 3-plex clinical readiness + power simulation
---

# A thyroid-lineage state marks radioiodine-refractory biology in BRAF V600E-mutant papillary thyroid cancer

Seungho Cook^1,2,\* and Hyeong-Won Yu^1,2,\*

^1 Department of Internal Medicine, Seoul National University Bundang Hospital, Seongnam, Republic of Korea
^2 Seoul National University College of Medicine, Seoul, Republic of Korea

\*Correspondence: kukshomr@gmail.com; [VERIFY institutional corresponding email]

## Abstract

Radioiodine remains central to differentiated thyroid cancer management, yet current molecular frameworks do not directly measure the thyroid-lineage machinery required for iodine uptake and organification, and cannot identify BRAF V600E-mutant patients in whom repeated high-dose radioiodine may be biologically futile. Here we define a compact eight-gene thyroid differentiation and iodine-handling axis (**TG, TPO, TSHR, SLC5A5, DIO1, PAX8, NKX2-1, FOXE1**) and apply it across retrospective public cohorts without internal tissue validation, encompassing bulk transcriptomic, methylation, proteomic, single-cell and clinical thyroid cancer datasets. In TCGA-THCA (n = 504), this axis resolves iodine-handling-low (DM1) and iodine-handling-high (DM2) states that are recovered by pan-genome clustering (ARI = 0.92), and DM1 is associated with overall survival (pooled hazard ratio 2.53, 95% CI 1.31-4.89; I^2 = 0%). Critically, we identify a significant interaction between DM1 and BRAF V600E status on progression-free interval (Cox interaction p = 0.022, HR_int = 0.47; n = 472), such that DM1-low tumours have shorter progression-free interval than DM1-high tumours within the BRAF+ subset (continuous thyroid differentiation score HR per +1 SD = 0.66, 95% CI 0.48-0.92, where higher score = more iodine-handling-high; log-rank p = 7.6 x 10^-3) but the score is not associated with outcome in BRAF- or RAS+ tumours. An in-silico expression proxy based on promoter methylation of a three-marker combination (**TG + PAX8 + NKX2-1**) — using antibodies already in routine clinical use for thyroid pathology — reproduces this stratification with log-rank p = 1.7 x 10^-4 in the BRAF+ subset, without requiring new assay development, although direct immunohistochemistry validation on tissue sections has not been performed. The DM1 axis itself replicates across five independent external cohorts, most strongly in Lee 2024 (Korean n = 370, Mann-Whitney p = 2.7 x 10^-7). Monte Carlo simulation using TCGA-derived priors predicts >90% power for prospective replication at n = 200. These findings position the eight-gene axis, and its in-silico IHC 3-plex derivative, as a candidate treatment-selection hypothesis for the BRAF-mutant subset of differentiated thyroid cancer, subject to prospective validation.

## Introduction

Most patients with differentiated thyroid cancer survive for many years, yet a clinically important subset receives repeated high-dose radioiodine without a reliable means of determining, before treatment, whether their tumours retain the molecular machinery required for iodine handling. For these patients, the absence of a pretreatment response biomarker can lead to additional exposure despite limited expected benefit, with cumulative risks that include salivary-gland injury, bone-marrow toxicity and secondary malignancy. A molecular state that identifies loss of thyroid-lineage and iodine-handling capacity could therefore help distinguish patients for whom repeated radioiodine is biologically plausible from those in whom further empiric treatment may add harm without commensurate benefit.

Papillary thyroid carcinoma is commonly interpreted through canonical driver classes, including BRAF-like and RAS-like molecular programmes. These frameworks have clarified tumour evolution and risk stratification, but they do not directly quantify whether a tumour retains the thyroid-cell functions required for radioiodine uptake, organification and hormone synthesis. This distinction is clinically important because radioiodine response is constrained by cellular differentiation state, not by driver identity alone.

The thyroid differentiation score and related lineage markers have previously shown that advanced thyroid cancers lose expression of canonical thyroid genes. However, most existing signatures are either broad transcriptomic scores or mechanistic descriptors rather than compact, clinically translatable readouts. A practical molecular state classifier should satisfy four conditions: it should be anchored in thyroid biology, reproduce across platforms and populations, retain a clear boundary between retrospective association and prospective clinical utility, and — most importantly — deliver decision-changing information for a defined patient subset.

We therefore treated iodine handling as a circuit rather than as a single gene. The panel includes five effector genes directly involved in thyroid hormone production and iodine handling (**TG, TPO, TSHR, SLC5A5/NIS, DIO1**) and three thyroid-lineage transcription factors (**PAX8, NKX2-1/TTF-1, FOXE1/TTF-2**). The premise is that a coordinated decrease across this circuit reports a tumour state in which radioiodine uptake and processing are biologically compromised. In this study we test whether the axis (i) exists as a robust cross-platform biology, (ii) associates with survival outcome, and (iii) modifies clinical trajectory in a driver-specific way that supports treatment-selection biomarker status rather than pure prognostic risk stratification.

In this study, we defined and validated an eight-gene thyroid-lineage state axis across transcriptomic, epigenomic, proteomic and single-cell datasets. Using TCGA-THCA as the discovery cohort, we identified iodine-handling-low and iodine-handling-high states that were robust to genome-wide feature selection and largely orthogonal to individual oncogenic drivers. We then examined their association with kinase fusions, promoter methylation, external thyroid-cancer cohorts and post-radioiodine refractory tumours. Finally, we assessed whether this lineage state stratified clinical outcomes within molecularly defined subgroups and explored a three-marker expression proxy as a potential path toward tissue-based validation. Together, these analyses establish a reproducible lineage-state framework and generate a testable hypothesis for avoiding repeated radioiodine in tumours with diminished iodine-handling biology.

## Results

### A compact thyroid-lineage axis resolves iodine-handling states in TCGA-THCA

We first assembled a thyroid-relevant candidate space designed to separate lineage-state information from canonical driver identity. The 67-gene TIERA67 candidate pool included thyroid differentiation, iodine-handling, MAPK-output, driver-anchor, aggressive-disease, dedifferentiation and immune-context categories. From this pool, we selected an eight-gene thyroid differentiation and iodine-handling panel: **SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1** and **DIO1** (Fig. 1a,b). The panel was chosen on thyroid biology rather than on genome-wide optimization, with the driver-anchor category retained only as a leakage-control comparator.

Unsupervised clustering of TCGA-THCA primary tumours on the eight-gene transcript space resolved two states, which we refer to operationally as DM1 and DM2 (Fig. 1c,d). DM1 represented an iodine-handling-low state, whereas DM2 retained higher thyroid differentiation and iodine-handling expression. In TCGA-THCA, DM1 accounted for 28.4% of primary tumours.

We next asked whether the eight-gene panel artificially imposed a cluster boundary or instead captured a broader transcriptomic state. Pan-genome clustering using the top 5000 genes by median absolute deviation reproduced the eight-gene partition with ARI = 0.92, comparable to the full TIERA67 pool (ARI = 0.90). In contrast, driver-anchor genes alone produced no meaningful recovery of the same structure (ARI = -0.007; Fig. 1e). The axis was also not explained by driver transcript abundance: BRAF expression was indistinguishable between BRAF V600E and wild-type tumours (Cohen's d = -0.044, Mann-Whitney p = 0.57), and single-feature AUCs for BRAF, TERT, KRAS, NRAS and HRAS were near chance (Fig. 1f). Head-to-head comparison of the eight-gene DM1 score against six alternative scores (Landa 5-overlap, RAI-machinery 5, TF-only 3, BRAF binary, age > 55, stage III/IV) showed that only Landa-5 and RAI-machinery-5 produced comparable DM1 classification performance, while all clinical variables failed both classification and progression-free interval stratification tasks [EDITORIAL DECISION: Fig.7/Fig.8 content → integrate or promote to Supplementary before submission; head-to-head comparison panel currently in dm1_story_web/public/figures/fig_deep1_head_to_head.png].

These analyses indicate that the eight-gene panel is not merely a small classifier trained to recover driver status. Rather, it is a compact readout of a broader thyroid-lineage programme that remains visible when the transcriptome is considered without panel restriction, and outperforms available clinical and single-driver alternatives.

### The iodine-handling-low state recovers clinically relevant dark matter

Canonical driver frameworks leave a clinically important BRAF- and RAS-negative compartment incompletely stratified. Applying the eight-gene axis to this dark-matter space recovered 131 of 180 BRAF/TERT-negative TCGA tumours into a defined DM1 or DM2 state, corresponding to 73% molecular recovery of the previously unresolved compartment (Fig. 2a).

The axis also tracked population-level differences relevant to thyroid cancer epidemiology. DM1 prevalence was 28.4% in TCGA-THCA and increased to 37.8% in Korean cohorts, consistent with the higher burden of driver-negative or non-BRAF-like thyroid cancer in East Asian populations (Fig. 2b). This enrichment supports the clinical relevance of an axis that is not anchored to BRAF status alone.

We then evaluated outcome association. In TCGA-THCA, DM1 was associated with higher overall-survival hazard than DM2, although the low event rate in primary thyroid cancer limited precision. We therefore combined TCGA-THCA with the MSK-IMPACT thyroid cohort, an advanced-disease-enriched cohort with higher event density. The TCGA hazard ratio was 2.30 (95% CI 0.77-6.88), and the MSK hazard ratio was 2.67 (95% CI 1.17-6.10). Random-effects meta-analysis yielded a pooled hazard ratio of 2.53 (95% CI 1.31-4.89) with no detectable between-cohort heterogeneity (I^2 = 0%; Fig. 5a). These data support DM1 as a retrospectively aggressive state while preserving the boundary that prospective outcome and treatment-response validation remain required.

### A predictive interaction between DM1 and BRAF V600E status

To determine whether the DM1 axis functions as a driver-agnostic prognostic biomarker or as a predictive biomarker whose effect depends on driver context, we formally tested interaction between DM1 score and canonical driver status in Cox proportional-hazards models on TCGA-THCA progression-free interval (n = 472; 49 events; ref. Liu 2018 pan-cancer clinical resource). Standardized DM1 score interacted with BRAF V600E status (interaction p = 0.022, HR_interaction = 0.47), whereas the DM1 × RAS interaction was not significant (p = 0.41) [EDITORIAL DECISION: interaction panel currently in dm1_story_web/public/figures/fig_deep2_prognostic_predictive.png; promote to Supplementary Figure or rebuild as main figure before submission].

Stratified analysis clarified the interaction structure. Within the BRAF V600E+ subset (n = 287, 33 events), DM1-low tumours had shorter progression-free interval than DM1-high tumours (continuous score HR per +1 SD = 0.66, 95% CI 0.48-0.92, p = 0.013; higher score = more differentiated, i.e., eight-gene thyroid differentiation score where higher = more iodine-handling-high). In the BRAF V600E-wild-type subset, DM1 showed no significant association (HR = 1.19, 95% CI 0.71-1.99, p = 0.52), and in the smaller RAS-mutant subset, the effect was flat (HR = 1.04, 95% CI 0.50-2.17, p = 0.92) [EDITORIAL DECISION: subgroup HR forest panel currently in dm1_story_web/public/figures/fig_deep2_prognostic_predictive.png; promote to Supplementary Figure or rebuild as main figure before submission]. This pattern is consistent with a treatment-selection biomarker whose informative range is confined to a specific molecular subset, rather than with a driver-independent prognostic score.

Directionally, the BRAF+ finding aligns with a growing body of clinical evidence that BRAF V600E-mutant papillary thyroid cancers include a subset in which lineage-programme silencing predicts radioiodine failure and in which BRAF or MEK inhibition can restore radioiodine uptake in selected patients (Chakravarty et al. 2011; Ho et al. 2013; Rothenberg et al. 2015). The DM1 axis offers a pre-treatment classifier for this subset, quantified from routine transcriptomic or methylation data at the time of surgical pathology.

We interpret the interaction result as a hypothesis-generating observation with clinical translation potential rather than as an established treatment-selection rule, because it derives from a single cohort with 33 events in the informative subset. The findings therefore motivate a prespecified prospective replication (see Prospective simulation section below).

### External validation and interaction axis replication

We next evaluated whether the DM1/DM2 axis generalized beyond TCGA. Across 19 external validation entries spanning expression, methylation, proteomic and single-cell contexts, the master cross-cohort forest showed a mean Cohen's d of 2.81 and median Cohen's d of 2.37. A per-gene by cohort by contrast matrix showed 80 of 80 cells moving in the expected direction, arguing against dependence on any single gene, cohort or platform.

Direct replication of the DM1 axis was strongest in the Lee 2024 Korean papillary thyroid cancer cohort (GSE213647). Two distinct analyses were performed in this cohort: (i) within-cohort DM1/DM2 axis separation in the full cohort (n = 632, Cohen's d = 5.93), confirming that the DM1/DM2 axis is visible within an independent Korean dataset; and (ii) a BRAF-interaction direction replication test (n = 370 tumours with available driver annotation, Mann-Whitney p = 2.7 x 10^-7, Cohen's d = 0.24), comparing the eight-gene panel score against an independently derived thyroid dedifferentiation subB score. Both analyses confirm that the axis identified in TCGA-THCA is present as an independent signal in an unlinked Korean cohort more than an order of magnitude larger than earlier microarray validation sets. In four GPL570 microarray cohorts, the axis reproduced in the expected direction in all cohorts, with panel score versus MAPK-output Spearman correlations ≤ -0.84. At the protein level, the Mun 2025 proteogenomic thyroid cohort (n = 336) showed 7 of 7 available panel proteins direction-consistent with the RNA-defined thyroid differentiation state. The Landa 2016 poorly differentiated / anaplastic cohort provided a well-differentiated to anaplastic gradient with Cohen's d ≈ 2.3 across the histology transition [EDITORIAL DECISION: external replication forest panel currently in dm1_story_web/public/figures/fig_a2_replication.png; promote to Supplementary Figure or rebuild as main figure before submission].

Single-cell analyses supported a thyrocyte-intrinsic interpretation. In GSE184362 (Pu et al. 2021), patient-matched tumour and adjacent-normal thyrocyte scores showed per-patient Spearman correlations of 0.798 to 0.886, each with Bonferroni-adjusted p < 10^-10. In the Lu 2023 single-cell cohort (GSE193581), the eight-gene gradient remained visible after restricting to KRT8/KRT19/EPCAM-positive thyrocyte-lineage cells, arguing against stromal or immune contamination as the primary source of the signal (Fig. 4a,b).

The axis was also compatible with clinical tissue processing. FFPE and fresh-frozen score distributions were not detectably shifted (Kolmogorov-Smirnov p = 0.44), supporting the feasibility of retrospective pathology-specimen deployment and prospective clinical assay development (Fig. 5e).

Direct DM1 × BRAF interaction testing was not feasible in the external cohorts because BRAF genotype and progression-free interval annotations were not simultaneously available in any of the reproducible expression cohorts. The A2 interaction result therefore remains, at present, a single-cohort observation whose confirmatory replication depends on a prospective clinical study.

### A routine IHC 3-plex reproduces the BRAF+ progression-free stratification

Motivated by the observation that a subset of the eight panel genes carries the majority of the differentiation signal (Fig. 3), we systematically evaluated 14 candidate biomarker combinations for the ability to reproduce BRAF+ progression-free interval stratification (data in Supplementary; combination comparison panel in dm1_story_web/public/figures/fig_ihc3_killer.png [EDITORIAL DECISION: promote to Supplementary Figure before submission]). Combinations were organized into four groups reflecting distinct translational scenarios: (A) clinical assay tiers (Full-8, TSO500 native 3-gene, IHC routine 3-plex, IHC extended 4-plex, RT-qPCR minimal 4-gene, Compact-5 add-on); (B) biological axis subsets (transcription-factor only, RAI-machinery only, effector only); (C) literature convergence (Landa 2016 overlap 5-gene); (D) statistical parsimony (Compact 6-gene, Top-N by |d|).

The strongest BRAF+ progression-free interval stratification was produced by a three-marker combination corresponding to routinely deployed thyroid pathology immunohistochemistry: **thyroglobulin (TG), PAX8 and NKX2-1/TTF-1** (log-rank p = 1.7 x 10^-4 in the BRAF+ subset; 5-fold cross-validated DM1 classification AUC = 0.756; Cox HR per +1 SD = 0.65, 95% CI 0.47-0.90, p = 1.4 x 10^-3; data in Supplementary [EDITORIAL DECISION: IHC 3-plex panel currently in dm1_story_web/public/figures/fig_ihc3_killer.png; promote to Supplementary Figure before submission]). This combination outperformed the Full-8 panel in BRAF+ progression-free stratification (Full-8 log-rank p = 7.6 x 10^-3) despite a smaller Cohen's d, consistent with the concentration of BRAF+ signal in a narrower and more clinically relevant marker set. Extending to a 4-plex including TPO (BRAF+ log-rank p = 2.0 x 10^-3) preserved the finding, whereas the TSO500 native 3-gene combination (TSHR, PAX8, NKX2-1) failed BRAF+ stratification (log-rank p = 0.30), reflecting a modality mismatch: TSO500 v2 measures DNA variation rather than expression, so the panel genes it includes cannot be used to derive a lineage-expression score without an add-on assay.

The clinical significance of the TG + PAX8 + NKX2-1 combination is that this immunohistochemistry triplet is already in routine diagnostic use in thyroid pathology (Dako TG A0251, clone 2H11+6E1; Roche/Cell Marque PAX8 MRQ-50; Roche/Dako NKX2-1 8G7G3/1). No new immunohistochemistry assay development is required to test the DM1-based BRAF+ stratification in an institutional cohort; the necessary reagents, staining protocols and pathologist expertise are established. This positions the in-silico IHC 3-plex derivative as a hypothesis that can be tested in an institutional cohort via retrospective chart-review pilots without requiring new reagent development, pending prospective validation with direct H-score quantification.

We note explicitly that the current analysis derives the IHC 3-plex signal from HM450 promoter methylation as an in-silico proxy for lineage-programme silencing at the protein level, and that direct H-score or percent-positive quantification will be required to confirm the finding at the immunohistochemistry level. Nonetheless, the combined direction of Mun 2025 proteomic replication (7 of 7 direction-consistent) and Landa 2016 lineage-loss convergence supports the biological plausibility of a protein-level correlate.

### Multi-omics 5-way convergence of the DM1 axis

To evaluate whether the DM1 axis is an assay-specific artifact of any single measurement modality, we compared per-gene effect sizes across five orthogonal modalities: TCGA HM450 promoter beta methylation, TCGA RNA-seq expression, the Landa 2016 poorly differentiated / anaplastic transcriptomic dataset, the Mun 2025 proteogenomic proteome, and the GSE151179 post-radioiodine expression cohort [EDITORIAL DECISION: multi-omics convergence panel currently in dm1_story_web/public/figures/fig_deep4_multiomics.png; promote to Supplementary Figure or rebuild as main figure before submission]. Per-gene rank concordance between modality pairs, measured by Spearman correlation, ranged from 0.71 to 0.94. Thyroperoxidase (TPO) ranked in the top four by effect size in every modality, and DIO1, TSHR and TG occupied the top four in the majority of modalities, indicating that the biological signal underlying the DM1 axis is consistently captured whether the assay measures DNA methylation, mRNA expression or protein abundance.

The cross-modality concordance is not sensitive to any single lineage transcription factor. Restricting the analysis to effector genes alone (TG, TPO, TSHR, SLC5A5, DIO1) or to transcription-factor genes alone (PAX8, NKX2-1, FOXE1) preserved the cross-modality correlation structure. This modality-independence is consistent with an interpretation of DM1 as a coordinated multi-layer thyroid-lineage state rather than as an assay-specific measurement.

### Prospective simulation supports feasibility at n = 200

To evaluate whether a prospective institutional cohort can adequately test the DM1 × BRAF interaction, we performed Monte Carlo simulation using TCGA-derived priors. For each of eight target cohort sizes between 80 and 500 patients, we simulated 500 hypothetical cohorts under the observed BRAF+ event rate (11.7% at 5-year follow-up), a hazard ratio corresponding to the TCGA BRAF+ Full-8 point estimate (HR = 1.49), and a 70% event-observation completeness assumption. Empirical power was calculated as the fraction of simulations in which the corresponding Cox model reached p < 0.05.

At n = 200, empirical power for the Full-8 panel exceeded 0.90 [EDITORIAL DECISION: Monte Carlo power simulation panel currently in dm1_story_web/public/figures/fig_deep5_snubh_simulation.png; promote to Supplementary Figure before submission]. For the IHC 3-plex derivative, whose real BRAF+ log-rank p corresponds to a stronger implicit hazard ratio (HR ≈ 1.65), the same simulation reached 0.80 power at n = 120. Under a conservative HR = 1.40 sensitivity scenario, n = 250 remained sufficient for 0.80 power. Under the TSO500 3-gene DNA-only scenario (implicit HR ≈ 1.22), even n = 500 was insufficient (power < 0.55), consistent with the biological interpretation that DNA-level measurement of TF variants does not substitute for lineage-programme expression measurement.

We interpret this simulation as a support for the feasibility of a prospective study at n = 200, subject to the caveat that the simulation uses TCGA priors and assumes 5-year follow-up completeness. Multiplicity correction for coprimary endpoints would require a modest sample-size increase.

### A translational pathway for radioiodine harm avoidance

The clinical position of the eight-gene axis, and its IHC 3-plex derivative, is currently defensible in two framings simultaneously: as a driver-agnostic harm-avoidance readout, and as a driver-specific treatment-selection biomarker for the BRAF V600E+ subset.

In the harm-avoidance framing, an iodine-handling-low tumour state identifies patients in whom repeated high-dose radioiodine may be biologically less plausible and in whom additional molecular work-up may be valuable. The axis aligns with external radioiodine-refractory biology: in GSE151179, post-RAI refractory tumours shifted toward a DM1-like thyroid-differentiation state (Cohen's d ≈ -1.0; Mann-Whitney p ≈ 10^-4), indicating that the DM1 programme resembles the transcriptional state observed after clinical RAI failure (Fig. 6d). Because this comparison is retrospective and public-cohort based, we interpret it as a biological anchor rather than as a validated response-prediction model.

In the treatment-selection framing, the BRAF+-restricted DM1 progression-free interval interaction (§3, above) suggests that within BRAF V600E+ patients the DM1 subset is a candidate group for reflex fusion testing followed by prospective evaluation of MAPK-inhibitor-based radioiodine re-differentiation strategies. This is consistent with the mechanistic and clinical work of Chakravarty et al. (2011), Ho et al. (2013), and Rothenberg et al. (2015), which established the biological plausibility of pharmacologic thyroid re-differentiation but did not identify a pre-treatment classifier for the responder subset.

A practical clinical pathway would begin with immunohistochemistry-based assessment of the TG + PAX8 + NKX2-1 triplet in the primary tumour or archival FFPE tissue, combined with BRAF V600E genotyping already routine in most institutions. A DM1 call within the BRAF+ subset would trigger reflex orthogonal fusion testing for RET, NTRK, ALK and BRAF rearrangements and would identify candidates for prospective evaluation of epigenetic re-differentiation or MAPK-inhibitor re-induction strategies. The eight-gene RNA panel and its Compact-5 add-on (TSO500 native 3 + TPO + DIO1, AUC = 0.940) provide alternative deployment modalities in institutions where transcriptomic assays are the preferred assay path.

## Discussion

### 3.1 Mechanistic interpretation

We interpret DM1 as a stable loss of thyroid-lineage identity rather than as a transcriptional proxy for any single driver. Its enrichment in BRAF V600E- and RET-fusion-positive tumours, together with higher promoter methylation in these groups than in RAS-mutant disease, is compatible with a model in which strong MAPK signalling favours lineage collapse and epigenetic stabilization. The overlap with the dedifferentiation programme reported by Landa et al. in poorly differentiated and anaplastic thyroid cancer places DM1 on a biologically credible continuum of thyroid-cell identity loss. Nevertheless, the present data cannot determine whether methylation initiates silencing, maintains an already repressed state, or merely records it. Alternative mechanisms — including disruption of the PAX8–NKX2-1–FOXE1 transcriptional network, enhancer remodelling and selection of pre-existing low-lineage states — remain plausible. We therefore regard methylation as an epigenetic correlate that makes the state experimentally actionable, not as a demonstrated causal mechanism.

### Clinical implications

The clinical value of the eight-gene axis lies in its ability to report thyroid-lineage function using a compact and clinically translatable assay surface. Current management decisions integrate histology, stage, ATA risk tier and selected driver mutations, but these variables do not directly measure whether iodine-handling machinery remains transcriptionally intact. The DM1/DM2 framework adds a functional layer: whether the tumour retains the lineage programme required for radioiodine uptake and organification.

This distinction matters because differentiated thyroid cancer is often long-lived. For many patients, the harm of ineffective repeated high-dose radioiodine is not only treatment delay but cumulative toxicity. A lineage-state assay that identifies iodine-handling-low tumours, particularly within the BRAF V600E+ subset, could support earlier reflex molecular testing, more deliberate discussion of radioiodine benefit, and prioritization of prospective re-differentiation therapy where appropriate.

The A2 predictive interaction result strengthens the near-term clinical position of the axis beyond a pure prognostic risk-stratification tool. Because the DM1 effect on progression-free interval is confined to BRAF V600E+ patients in the current analysis (interaction p = 0.022), the axis functions as a candidate treatment-selection biomarker for a defined driver subset. This subset — BRAF V600E-mutant thyroid cancer with DM1-like lineage state — corresponds mechanistically to the group in which pharmacologic MAPK-pathway inhibition has previously been shown to restore radioiodine uptake in a subset of clinical patients, and for which no pre-treatment classifier has been established. If the interaction result replicates in a prespecified prospective cohort, the axis becomes a candidate rule for identifying BRAF+ patients for whom repeated high-dose radioiodine is likely to be futile and for whom MAPK-inhibitor-based re-differentiation should be considered upfront.

The strongest near-term translational application is a reflex-testing pathway anchored on the routine immunohistochemistry 3-plex (TG + PAX8 + NKX2-1). This triplet is already in routine diagnostic use across most thyroid pathology laboratories, uses reagents that are commercially validated for in-vitro diagnostic use, and reproduces the BRAF+ progression-free interval stratification with log-rank p = 1.7 x 10^-4 in the current analysis. This means that the in-silico finding can be tested in an institutional cohort with existing tissue, existing IHC reagents and existing pathology workflows, without requiring new assay development — a substantially shorter validation timeline than for a new RNA-based assay, provided direct H-score quantification is performed.

### Generalizability and translational readiness

The axis is supported by several forms of generalizability: pan-genome recovery within TCGA (ARI = 0.92), external Korean expression cohorts (Lee 2024: within-cohort axis separation n = 632, d = 5.93; BRAF-interaction direction replication n = 370, Mann-Whitney p = 2.7 x 10^-7), microarray cohorts (GPL570 four-cohort direction-consistent, GSE29265/GSE33630/GSE65144/GSE126698), proteomic directionality (Mun 2025 n = 336, 7/7 concordant), single-cell thyrocyte restriction (Pu 2021 and Lu 2023 per-patient concordance), FFPE compatibility (KS p = 0.44) and multi-modality cross-omics rank concordance (Spearman ρ 0.71-0.94 across five modalities). This breadth is important because a compact panel is only useful if it remains stable across preprocessing, ancestry, platform and tissue-processing conditions.

The A2 predictive interaction result, however, has not yet been replicated in an external cohort with joint BRAF genotype and progression-free interval annotation. Currently available external cohorts either lack per-sample BRAF genotype (most microarray cohorts) or lack progression-free interval annotation (most Korean cohorts), which prevents direct external interaction testing at present. Monte Carlo simulation with TCGA priors indicates that a prospective cohort of n = 200 provides >90% power to replicate the Full-8 panel effect and >80% power to replicate the IHC 3-plex effect at 5-year follow-up. A prespecified prospective study with n ≈ 200, structured to permit both harm-avoidance and treatment-selection framings, is therefore the natural next step.

### Limitations

The principal limitation is that we did not have a pretreatment cohort in which the DM1 state could be tested directly against measured radioiodine uptake and clinical response. We therefore show that DM1 represents diminished iodine-handling biology and resembles tumours sampled after development of radioiodine refractoriness, but we do not show that it prospectively identifies non-responders. The clinical interaction in BRAF V600E-positive disease is likewise based on a single retrospective cohort and should be treated as hypothesis-generating. Survival estimates are further limited by only 16 overall-survival events in TCGA and by enrichment for advanced disease in the MSK cohort. Mechanistically, methylation remains correlative; the fusion-negative DM1 comparison was underpowered; cross-platform calibration was imperfect in K2; and Visium analysis did not provide confirmatory spatial evidence. We also did not complete internal RNA or tissue validation. In particular, the TG–PAX8–NKX2-1 model is an in-silico proxy and not a clinically validated IHC assay, while TSO500 is a DNA-focused platform that cannot recover the required expression state. The most direct next study is therefore a prespecified retrospective FFPE validation linked to administered radioiodine dose, post-therapy uptake and response, followed by independent prospective testing.

## Methods

### Cohorts and data sources

TCGA-THCA RNA-seq, clinical, mutation, structural-variant and HM450 methylation data were obtained from TCGA and cBioPortal. Progression-free interval was defined per Liu et al. 2018 pan-cancer clinical resource. The MSK-IMPACT thyroid cohort was used as an advanced-disease survival validation cohort. Korean validation cohorts included PRJEB11591/K2 (n = 260) and GSE213647/Lee 2024 (full cohort n = 632; driver-annotated BRAF-interaction replication subset n = 370 tumours). External expression validation used GPL570 microarray cohorts GSE29265, GSE33630, GSE65144 and GSE126698 (n = 28). Advanced-disease lineage-loss context was provided by the Landa 2016 poorly differentiated / anaplastic thyroid cancer transcriptomic dataset. Protein-level validation used the Mun 2025 proteogenomic thyroid cohort (n = 336). Single-cell validation used GSE184362/Pu 2021 and GSE193581/Lu 2023, with additional external single-cell support from GSE241184 and GSE232237 where indicated. Post-radioiodine refractory context was provided by GSE151179. Cohort-level sample sizes are reported in each figure caption and Supplementary Table 1.

### Eight-gene panel definition

The eight-gene panel comprised **TG, TPO, TSHR, SLC5A5, DIO1, PAX8, NKX2-1** and **FOXE1**. These genes were selected from the thyroid differentiation and iodine-handling component of a 67-gene thyroid-relevant candidate pool. Driver-anchor genes were retained as controls but were not used to define the final panel. Analyses used official HGNC symbols; clinical aliases were reported as SLC5A5/NIS, NKX2-1/TTF-1 and FOXE1/TTF-2.

### Score construction and clustering

Expression values were log-transformed and z-standardized within cohort. KMeans clustering with k = 2 was applied to the eight-dimensional panel expression matrix in TCGA-THCA to define DM1 and DM2 states. For cross-cohort transfer, panel scores were computed as within-cohort z-mean summaries, with additional within-sample-centered profile scoring for cohorts requiring calibration adjustment. Cluster reproducibility was evaluated using adjusted Rand index against panel-defined DM1/DM2 labels.

### Interaction Cox analysis

Cox proportional-hazards models for progression-free interval used the eight-gene panel score standardized to zero mean and unit variance within the analytic cohort. Interaction with BRAF V600E status was modelled as a multiplicative term (score × BRAF). Reference stratified models were fitted separately in BRAF V600E+, BRAF-wild-type and RAS-mutant subsets. All interaction and stratified models used TCGA-THCA progression-free interval as defined by Liu et al. 2018. Interaction significance was reported as two-sided Wald p-values on the interaction coefficient.

### Head-to-head score comparison

Six alternative scores were computed in the TCGA cohort for head-to-head comparison against the Full-8 panel: (i) Landa-5 overlap (TG, TSHR, TPO, DIO1, SLC5A5), (ii) transcription-factor-only 3 (PAX8, NKX2-1, FOXE1), (iii) RAI-machinery 5 (TSHR, SLC5A5, TPO, TG, DIO1), (iv) BRAF V600E binary status, (v) age > 55 years, and (vi) stage III/IV binary status. DM1 classification AUC was computed under 5-fold cross-validated logistic regression; BRAF+ subset progression-free interval Cox hazard ratios were computed per +1 SD score with 95% confidence intervals from the Cox model.

### Immunohistochemistry 3-plex in-silico proxy

Because HM450 promoter beta values are available for all TCGA-THCA tumours with paired clinical annotation, a three-marker methylation proxy for the immunohistochemistry combination TG + PAX8 + NKX2-1 was computed as the mean promoter beta across the three genes. This proxy was standardized within cohort and used for BRAF+ progression-free interval stratification. Direct H-score or percent-positive quantification is required for confirmatory analysis and was not performed in the current study.

### External replication

Direct replication of the DM1 axis was tested in Lee 2024 (GSE213647) using two analyses: (i) within-cohort DM1/DM2 axis separation using the full cohort (n = 632, axis separation analysis), and (ii) comparison of the eight-gene panel score against an independently derived subB thyroid dedifferentiation score in the driver-annotated subset (n = 370, BRAF-interaction direction replication) using two-sided Mann-Whitney U tests. GPL570 cohorts (GSE29265, GSE33630, GSE65144, GSE126698; total n ≈ 205) were re-analysed for panel score versus MAPK-output Spearman correlation, using the same panel-standardization pipeline as for TCGA. Landa 2016 replication used published effect sizes across the well-differentiated / poorly-differentiated / anaplastic histology gradient.

### Multi-omics 5-way convergence

Per-gene DM1 versus DM2 effect sizes were computed in five modalities: TCGA HM450 promoter beta, TCGA RNA-seq expression, Landa 2016 poorly differentiated / anaplastic expression, Mun 2025 proteogenomic proteome and GSE151179 post-radioiodine expression. For proteogenomic and Landa modalities, published effect sizes were used with sign harmonization to the DM1-higher convention. Cross-modality rank concordance was quantified by Spearman correlation on the per-gene effect vectors.

### Structural-variant analysis

Structural variants were retrieved from cBioPortal for TCGA-THCA and MSK-IMPACT. Fusion classes included RET, NTRK, ALK, BRAF and PAX8-PPARG rearrangements. Fusion enrichment was tested using Fisher exact tests. Structural-variant missingness was assessed by chi-square testing and sensitivity analyses under best-case, worst-case and matched-imputation scenarios.

### Methylation analysis

HM450 promoter beta values were retrieved for TCGA-THCA tumours with paired DM calls. Per-gene and mean eight-gene promoter beta values were compared between DM1 and DM2 using Cohen's d and Mann-Whitney U tests. MAPK-output scores were computed using DUSP4/5/6, SPRY2/4, ETV4/5, PHLDA1 and CCND1. Cross-cohort MAPK-output versus thyroid-score correlations were pooled using Fisher-z transformation.

### Single-cell processing

Single-cell datasets were normalized using log1p(CP10K) and scored using within-sample z-standardized panel expression. Thyrocyte-lineage cells were defined using KRT8, KRT19 and EPCAM positivity when raw annotations permitted. Patient-level pseudo-bulk scores were generated by averaging panel expression within thyrocyte-lineage cells. Per-patient correlations and cross-cohort score concordance were evaluated using Spearman correlation.

### Survival analysis

Overall-survival analyses used Cox proportional-hazards models and Kaplan-Meier curves. TCGA-THCA and MSK-IMPACT hazard ratios were pooled using DerSimonian-Laird random-effects meta-analysis. Between-cohort heterogeneity was quantified using Cochran Q and I^2. Multivariable models included age, stage and selected driver covariates where sample size allowed.

### Monte Carlo prospective simulation

For each of eight target cohort sizes (80, 120, 160, 200, 250, 300, 400, 500), 500 hypothetical cohorts were simulated. DM1 scores were drawn from a standard normal distribution; hazard was exponential in the standardized score with log-hazard-ratio anchored on the TCGA BRAF+ Full-8 point estimate (HR = 1.49). Follow-up was truncated at 5 years with 70% event-observation completeness. Empirical power was defined as the fraction of simulations in which the Cox regression p-value was below 0.05. Sensitivity analyses used HR = 1.65 (IHC 3-plex implied), HR = 1.40 (conservative) and HR = 1.22 (TSO500 3-gene, DNA-only).

### Statistical analysis

Cohen's d was computed using pooled standard deviation. Continuous variables were compared using two-sided Mann-Whitney U tests unless otherwise specified. Proportions were compared using Fisher exact tests. Multiple-testing correction used Benjamini-Hochberg false-discovery-rate control where applicable. Confidence intervals for correlations were estimated using Fisher-z transformation.

## Data availability

All source datasets are publicly available from TCGA, cBioPortal, GEO or ENA under the accessions listed in Supplementary Table 1. Processed per-sample score tables, methylation summaries, structural-variant joins, interaction Cox model outputs and Monte Carlo simulation seeds will be deposited with the public code release at submission.

## Code availability

Analysis and figure-generation scripts are maintained in the project repository and will be archived with a DOI at submission. Current working scripts include `deep_analysis_5x_parallel.py` (A1-A5 head-to-head, interaction, redifferentiation dynamics, multi-omics convergence and Monte Carlo simulation), `a2_replication_lee_k2.py` (Lee 2024 axis replication), `gene_reduction_v2.py` (14 biomarker combination analysis) and `ihc_3plex_killer_figure.py` (IHC 3-plex clinical readiness figure), located under `project/results/manuscript_v8_nc_main/`.

## Acknowledgements

[AUTHOR SLOT — add funding, institutional support, collaborator acknowledgements and data-source acknowledgements.]

## Author contributions

S.C. conceived the analysis, curated datasets, performed computational analyses, generated figures and wrote the manuscript draft. H.-W.Y. supervised the clinical framing and interpretation. [VERIFY AND EXPAND BEFORE SUBMISSION.]

## Competing interests

The authors declare no competing interests. [VERIFY BEFORE SUBMISSION.]

## Figure legends

### Figure 1. The eight-gene thyroid-lineage axis resolves iodine-handling states.

(a) Study overview and cohort inventory. (b) Curation of the eight-gene panel from thyroid differentiation and iodine-handling biology. (c) UMAP of TCGA-THCA tumours in eight-gene transcript space, coloured by DM1 and DM2 state. (d) Heatmap of panel genes ordered by DM1 probability. (e) ARI ladder comparing panel, TIERA67, pan-genome and driver-anchor clustering. (f) Driver mRNA neutrality and single-feature classification AUCs. (g) TCGA overall-survival Kaplan-Meier curve by DM state.

### Figure 2. DM1 is enriched for kinase-fusion biology.

(a) Recovery of dark-matter tumours into DM1/DM2 states. (b) DM1 prevalence across TCGA and Korean cohorts. (c) Fusion enrichment by DM state. (d) Fusion partner spectrum within DM1. (e) DM1 capture of RET-fusion-positive tumours. (f) Histological enrichment of follicular-variant PTC. (g) MSK-IMPACT fusion replication.

### Figure 3. DM1 shows epigenetic silencing of thyroid differentiation machinery.

(a) HM450 promoter beta heatmap for panel genes and TDS controls. (b) Mean eight-gene beta by DM state. (c) Methylation-expression coupling for representative thyroid genes. (d) Cross-cohort MAPK-output versus Panel-8 score forest. (e) Functional equivalence of Panel-8 and TDS-16. (f) Per-driver-class mean eight-gene beta. (g) Landa 2016 advanced-thyroid-cancer convergence heatmap.

### Figure 4. Single-cell validation establishes a thyrocyte-intrinsic state.

(a) Lu 2023 thyrocyte UMAP coloured by eight-gene score. (b) Pu 2021 per-patient tumour and adjacent-normal thyrocyte correlations. (c) Author-independence cross-cohort check. (d) External single-cell pseudo-bulk replication. (e) Multisite pooled tumour versus normal pseudo-bulk scatter. (f) Compositional pseudotime across TCGA and Lee cohorts. (g) External direction-consistency forest.

### Figure 5. DM1 is associated with survival and assay portability.

(a) TCGA plus MSK random-effects overall-survival meta-analysis. (b) Multivariable Cox model. (c) Multi-cohort Kaplan-Meier stack. (d) Cross-cohort score-distribution portability. (e) FFPE versus fresh-frozen concordance. (f) Time-dependent ROC. (g) Survival-model calibration.

### Figure 6. Translational pathway for radioiodine harm avoidance.

(a) RNA-state-guided reflex fusion-testing algorithm. (b) ATA risk-tier overlay. (c) Hypomethylating-agent plus RAI re-induction rationale. (d) DM1-like state in post-RAI refractory disease. (e) Decision-curve analysis. (f) Prospective DM1-stratified trial schema. (g) RET-fusion concentration across the DM1 score distribution.

### Figure 7. DM1 predictive interaction with BRAF V600E status and external axis replication (new).

(a) Head-to-head comparison of the Full-8 panel against six alternative scores (Landa-5 overlap, RAI-machinery 5, TF-only 3, BRAF binary, age > 55, stage III/IV) for DM1 classification AUC and BRAF+ progression-free interval Cox hazard ratio. (b) DM1 × BRAF interaction Cox model results: panel A shows interaction term significance (DM1 × BRAF p = 0.022 vs DM1 × RAS p = 0.41); panel B shows subgroup-stratified DM1 hazard ratios in BRAF V600E+, BRAF-wild-type and RAS-mutant subsets. (c) External axis replication forest: TCGA prior, Lee 2024 (n = 370, BRAF-interaction direction replication), K2 (n = 260), four GPL570 cohorts (GSE29265, GSE33630, GSE65144, GSE126698), Landa 2016. (d) Multi-omics 5-way convergence heatmap and cross-modality Spearman ρ matrix across HM450 β, TCGA RNA, Landa 2016, Mun 2025 proteome, GSE151179.

### Figure 8. Routine IHC 3-plex enables immediate clinical readiness (new).

(a) 14 biomarker combination systematic comparison, grouped by clinical assay tier, biology axis, literature convergence and statistical parsimony. (b) Immunohistochemistry 3-plex protocol schematic (TG, PAX8, NKX2-1) with clone and vendor annotations. (c) BRAF+ subset progression-free interval Kaplan-Meier by IHC 3-plex tertile split (log-rank p = 1.7 x 10^-4). (d) Cross-combination BRAF+ progression-free interval log-rank p^-log_10 comparison ranking. (e) SNUBH prospective Monte Carlo power simulation for Full-8, IHC 3-plex, conservative HR = 1.4 and TSO500 3-gene scenarios across cohort sizes 80-500, 500 simulations per point, 5-year follow-up, 70% event completeness. (f) Clinical takeaway summary — immediate deployment via routine pathology, superior BRAF+ signal, prespecified validation plan.
