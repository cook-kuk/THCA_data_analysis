# Reviewer Defense v3 — 25 Anticipated Attacks

**Manuscript**: v6 ULTIMATE — "An 8-gene RAI-responsiveness biomarker is the deployable molecular quantification of the WHO 2022 morphologic axis in PTC"
**Last updated**: 2026-04-27 KST · post-ULTIMATE sprint
**v3 changes**: 5 new attacks (A21–A25) addressing literature integration + de-circularization

---

## A1–A20: Inherited from v2 (preserved)

[Numbered list A1–A20 covers: methodology, data sources, statistics, sample size, BRAF baseline, single-cohort risk, p-value cherry-picking, alternative endpoints, drug-actionability prematurity, scRNA scale, immune scoring methodology, TERT recovery, multivariate adjustment, decision curve protocol, NIFTP exclusion, ComBat-seq protocol, cross-platform z-score, cluster reproducibility (bootstrap), Korean BRAF rate generalisability, eight-gene rationale.]

Detailed answers preserved in `v17p35_REVIEWER_DEFENSE.md` (v2). The five additions below address the ULTIMATE-sprint findings.

---

## A21. "Why K=2 not K=3, K=4, or K=5?"

**Anticipated form**: "Pu et al. (Oncogene 2022) and others have demonstrated four molecular subtypes in PTC. Why does this manuscript collapse to a binary axis?"

**Answer**:
We explicitly evaluated K=2, 3, 4, 5 on the variant-A leak-free panel (47 genes, n=500 TCGA primary tumours). Internal-validity metrics favour K=2 monotonically: silhouette score 0.186 (K=2) > 0.164 (K=3) > 0.104 (K=4) > 0.091 (K=5); Calinski–Harabasz 82.2 (K=2), best across the range; Davies–Bouldin K=4 is worst (2.21).

We deploy K=2 because (i) it is statistically superior on every internal metric, (ii) the clinical decision is binary (RAI-responsive vs RAI-poor), and (iii) K=4 nests cleanly inside K=2 — each K=4 partition assigns to one of the two K=2 partitions plus a Pu 2022 canonical subtype name (Stromal, CNV-enriched, Immune-enriched, BRAF-enriched). We provide K=4 as a research overlay (Supplementary Fig. S20) reproducing Pu 2022 in a single self-consistent analysis. *Manuscript location*: Results, "K=2 vs K=4" subsection; Methods, "Cluster number selection".

---

## A22. "How does this differ from Han SC et al. 2023 (BL/RL-PTC)?"

**Anticipated form**: "Han SC, Park YJ et al. (ENM 2023) already described 'BL-PTC' and 'RL-PTC' transcriptomic phenotypes with CAF/immune characterisation. What is the novel contribution of this manuscript?"

**Answer**:
We agree the transcriptomic axis itself is not novel — it was described by Han 2023, the canonical BRAF-like/RAS-like distinction (Landa 2016), Pu 2021 at single-cell resolution, and Pu 2022 at four-subtype resolution. **Our contribution is translation to deployment**:

1. **An eight-gene minimal panel** (TPO, TG, NIS, TSHR, PAX8, NKX2-1, FOXE1, DIO1 — all NanoString or qPCR-measurable) that recovers the axis with cross-validation AUC 0.962, reproducible across three independent leak-free constructions (variants A, B, D).
2. **Decision curve analysis** demonstrating clinical net benefit dominating BRAF V600E baseline across threshold probability range 5%-95% — a clinical-utility argument absent from the discovery papers.
3. **Cross-modality methylation validation** in GSE97466 (independent cohort, 100% specificity for aggressive histologies → met_DM2).
4. **Explicit WHO 2022 mapping** — operationalising the 2022 morphologic reclassification (cPTC + infiltrative FVPTC = BRAF-like, IEFVPTC = RAS-like) molecularly.
5. **Korean parallel work positioning** — collaborative validation with the Han / Park YJ group is in progress (revision-round commitment).

We reframe in Discussion: "Our work is the deployable minimal version of axes previously described in pieces by Han 2023, Pu 2021, Pu 2022 — extending these by adding decision curve analysis, methylation cross-modality, leak-free meta-analysis, and explicit WHO 2022 alignment."
*Manuscript location*: Introduction (last paragraph) + Discussion (first three paragraphs).

---

## A23. "Pu et al. 2021 already described the BRAF-like-B subtype. Where is the novelty?"

**Anticipated form**: "Pu W et al. Nat Commun 2021 identified a 'BRAF-like-B' subtype with dedifferentiation predominance, CAF enrichment, and immunotherapy candidacy. The DM2 cluster appears to be the bulk-level proxy of BRAF-like-B. Is this not redundant?"

**Answer**:
We cite Pu 2021 explicitly as the discovery-level antecedent, and we agree DM2 is the bulk-level proxy of BRAF-like-B. The 88.9% direction concordance with Pu 2021's published dedifferentiation signature (8/9 marker genes; TG/TPO/DIO2 down, MYC/SOX4 up in DM2) is itself one of our manuscript's *evidence-of-validity* anchors — we use it as cross-study replication of the Pu 2021 finding at bulk scale.

**Our novelty over Pu 2021**: (i) a deployable eight-gene panel rather than a 158k-cell single-cell discovery requiring scRNA infrastructure; (ii) cross-cohort meta-analysis at n>2,000 instead of n=11 patients / 23 samples; (iii) decision-curve clinical-utility analysis; (iv) WHO 2022 morphologic mapping; (v) cross-modality methylation validation. Pu 2021 demonstrated *what* the axis is at single-cell resolution; we demonstrate *how* a clinician can read it without a single-cell platform.

*Manuscript location*: Discussion, "Korean parallel work / Pu 2021 / Pu 2022" paragraph; Fig. S21 (DM2 ↔ Pu BRAF-like-B 88.9% direction concordance).

---

## A24. "The eight-gene panel was used in defining the DM1/DM2 cluster — is the AUC circular?"

**Anticipated form**: "Reading the v5 preprint, the original DM1/DM2 cluster was constructed using a 67-gene curated panel (TIERA67) which includes the eight RAI-responsiveness markers. Predicting that cluster from those eight genes is autocorrelation, not prediction."

**Answer**:
**This is correct and we identified and corrected this issue ourselves before submission.** The v5 preprint's 0.954 cross-validation AUC was an upper bound, not a generalisable predictive estimate. The pre-submission self-audit pipeline (REALFIX R1–R4) reconstructed the DM1/DM2 axis four times with leak-free gene sets:

- **Variant A** (TIERA67 minus eight-gene, 47 genes, biology-preserving): leak-free AUC 0.962 [0.940-0.979]
- **Variant B** (MAPK + immune + EMT + cell-cycle, 30 genes, mechanistically motivated, zero overlap): AUC 0.925 [0.895-0.952]
- **Variant C** (5-feature immune composite, 5 genes): AUC 0.664 — the panel collapses on a non-differentiation axis, demonstrating biology specificity
- **Variant D** (BRS71 minus eight-gene, 26 genes, Chakravarty 2011 framework): AUC 0.957 [0.936-0.976]

Three of four leak-free variants recover AUC > 0.92 with 95% CIs not crossing the BRAF V600E baseline; ΔAUC vs BRAF in the [+0.111, +0.130] range. The eight-gene panel's predictive value is real, just not 0.954. We report the honest range transparently.

This is also our defence against the strongest competing claim — that the panel is merely an autocorrelation. **Three independent leak-free constructions agree.** *Manuscript location*: Results, "A circular-validation issue was identified and corrected"; Fig. R1–R4; supplementary tables U1A.

---

## A25. "WHO 2022 already separates IEFVPTC. Marginal molecular gain?"

**Anticipated form**: "If WHO 2022 already operationalises the BRAF-like / RAS-like axis morphologically (cPTC + infiltrative FVPTC = BRAF-like, IEFVPTC = RAS-like), what additional value does an eight-gene molecular panel provide?"

**Answer**:
WHO 2022 morphologic classification requires a histopathology consensus on encapsulation status, capsular invasion, and follicular pattern percentage. This consensus is reader-dependent (κ for IEFVPTC vs FVPTC infiltrative typically 0.6-0.7 in published inter-pathologist studies), is performed only at primary diagnosis on resected tissue, and cannot be repeated longitudinally on the same patient (e.g. on a recurrence biopsy or fine-needle aspirate).

The eight-gene panel provides:

1. **Consensus-independent reading**: NanoString or qPCR yields a numeric score with reproducible thresholds (P(DM1) > 0.5 = differentiated/RAS-like).
2. **Repeatable on FFPE/FNA**: applicable to FNA cytology specimens (where histopathology cannot confirm encapsulation status) and to recurrence biopsies, enabling longitudinal axis tracking.
3. **Quantitative output**: a continuous 0-1 score rather than a binary morphologic call, enabling clinicians to read the strength of evidence (a sample with P(DM1) = 0.85 carries different therapeutic implications than P(DM1) = 0.55).
4. **Direct link to RAI uptake**: the eight genes are the canonical thyroid-differentiation gene-set whose products are the molecular machinery of iodide uptake and organification (NIS = Na-I symporter, TPO = peroxidase, TG = thyroglobulin substrate). Higher panel score → mechanistically explained higher RAI uptake.
5. **Concordance with WHO 2022 quantified**: accuracy 84%, OR 20.4, p = 2.5×10⁻³³ (n=476 evaluable). The molecular reading is consistent with the morphologic axis where both apply, but applies in additional clinical contexts.

We do not displace WHO 2022 — we operationalise it for clinical contexts where morphology is unavailable, ambiguous, or non-repeatable. *Manuscript location*: Discussion, "WHO 2022 alignment"; Discussion limitations (TCGA does not split FVPTC encap vs infiltrative — main ceiling on observed concordance).

---

## A26 (bonus, internal): "The cluster orientation flipped between v5 and v6 — why should reviewers trust this?"

**Internal note for cover letter / response to reviewers**:

The v5 → v6 ULTIMATE re-derivation revealed that the v5 preprint's labelling convention ("DM1 = aggressive, DM2 = differentiated") was inconsistent with the actual leak-free cluster identity. Three independent lines of evidence in v6 agree on the corrected convention (DM1 = differentiated, DM2 = aggressive):

1. WHO 2022 mapping (U2A): DM2 contains 85.6% BRAF-like cPTC; DM1 contains 77.5% RAS-like FVPTC.
2. Pu 2021 dediff-signature direction (U2C): 88.9% direction concordance with DM2 (i.e. DM2 has dediff signature up, DM1 has differentiation signature up).
3. Methylation cross-modality (U1D): 100% of aggressive histologies (ATC, PDTC, FTC, Hürthle, mFTC) cluster in met_DM2.
4. GSE213647 (n=632, U1C): cPTC has higher 8-gene-DM1-score than ATC/PDTC, again consistent with DM1 = differentiated.

The biology — well-differentiated vs dedifferentiated axis — is preserved across all four constructions. The convention correction is purely a relabeling and is reported transparently in Results. We frame this as "the de-circularization pipeline also surfaced a labelling correction" rather than as a substantive change of conclusions.

---

## Submission readiness

| Attack class | n | Status |
|---|---|---|
| A1–A20 (v2 inherited) | 20 | answered, manuscript-cited |
| A21–A25 (ULTIMATE additions) | 5 | answered, manuscript-cited |
| A26 (internal cover-letter only) | 1 | for cover letter v5, not for response-to-reviewers |
| **Total** | **26** | **submission-ready** |

Cover letter v5 must include explicit acknowledgement of:
- Pre-submission self-audit revealing circular validation
- Honest re-evaluation across four leak-free constructions
- Cluster orientation correction with three independent lines of confirming evidence
- Korean parallel work outreach in progress (Park YJ / Bundang SNU)
