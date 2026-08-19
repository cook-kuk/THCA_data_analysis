---
title: Clinical review
audit_date: 2026-07-30
auditor: manuscript-audit-agent (clinical reviewer mode)
---

# 03 — Clinical review

## 1. "RAI harm avoidance" framing — is it honest?

**Study design:** Retrospective, public-cohort, in-silico analysis. No direct measurement of pretreatment radioiodine uptake or post-treatment response. No patient-level RAI dosing data. No comparison of outcomes in patients who did vs did not receive RAI stratified by DM status.

**Assessment:** The "harm avoidance" framing is clinically plausible as a motivating hypothesis but the current data do not establish it. The manuscript correctly frames GSE151179 as a "biological anchor" and "retrospective alignment" rather than prediction. However, the abstract (v2 line 40) uses the phrase "immediate translational deployment without new assay development" — this implies a clinical readiness that the evidence does not support. A clinician reading only the abstract would reasonably interpret the IHC result as closer to deployment than it actually is.

**Required correction:** The abstract must explicitly state the retrospective in-silico design and that no RAI response was directly measured. The CLAIM_LANGUAGE_MATRIX.md (Abstract rule) and DM1_PROJECT_STATE.md (Clinical motivation) both require this disclosure. The current abstract does not fulfill this requirement.

---

## 2. DM1 × BRAF interaction — "predictive biomarker" vs "hypothesis-generating"

**The evidence:** A genotype-stratified Cox PH interaction (p=0.022) on PFI in one cohort (TCGA-THCA, n=287 BRAF+, 33 events). No treatment assignment or randomization. No external cohort has simultaneously available BRAF genotype and PFI annotation that permits replication.

**Critical distinction:** A "predictive biomarker" formally means the biomarker predicts DIFFERENTIAL BENEFIT from a specific treatment in different patient subgroups (e.g., treatment A vs treatment B). The current finding is a PROGNOSTIC interaction — DM1 modifies the prognostic association between BRAF status and PFI, but without a treatment comparison, it cannot establish predictive value for any treatment decision.

**What the manuscript says:** The Introduction (v2 line 50) states the study tests whether the axis "modifies clinical trajectory in a driver-specific way that supports treatment-selection biomarker status rather than pure prognostic risk stratification." The Discussion (v2 lines 122, 142) uses "candidate treatment-selection biomarker." The Results (v2 line 78) appropriately qualifies: "consistent with a treatment-selection biomarker... rather than a driver-independent prognostic score." The v2 line 82 qualifies: "hypothesis-generating observation with clinical translation potential rather than as an established treatment-selection rule."

**Assessment:** The manuscript oscillates between appropriately qualified and overclaiming language. The framing is most problematic in:
- The abstract: "candidate treatment-selection framework" without adequate qualification
- The Discussion clinical implications section (v2 lines 122, 142-144): the language "the axis becomes a candidate rule for identifying BRAF+ patients for whom repeated high-dose radioiodine is likely to be futile" goes beyond what a single interaction term on PFI supports

**Required correction:** Add "hypothesis-generating" or "exploratory" to every use of "treatment-selection biomarker" in the abstract and Discussion. The Introduction framing of testing for "treatment-selection biomarker status" is forward-looking and permissible if the Results clearly state this was not achieved.

---

## 3. IHC 3-plex — "validated assay" vs "in-silico expression proxy"

**What the analysis does:** The IHC 3-plex result is derived from the mean HM450 promoter beta values for TG, PAX8, and NKX2-1 in TCGA-THCA. This is a DNA methylation proxy — not measured protein expression, H-score, or percent-positivity by immunohistochemistry.

**What the manuscript says:**
- v2 Methods (line 180) correctly states: "a three-marker methylation proxy for the immunohistochemistry combination TG + PAX8 + NKX2-1 was computed as the mean promoter beta across the three genes... Direct H-score or percent-positive quantification is required for confirmatory analysis and was not performed in the current study."
- v2 Results (line 104) correctly states: "the current analysis derives the IHC 3-plex signal from HM450 promoter methylation as an in-silico proxy for lineage-programme silencing at the protein level, and that direct H-score or percent-positive quantification will be required to confirm the finding at the immunohistochemistry level."
- However, the abstract (v2 line 40) states: "Reducing the panel to a three-marker immunohistochemistry combination (**TG + PAX8 + NKX2-1**) ... already in routine clinical use for thyroid pathology — reproduces this stratification with log-rank p = 1.7 x 10^-4 in the BRAF+ subset, positioning the finding for immediate translational deployment without new assay development."
- The Discussion (v2 line 144) states: "This triplet is already in routine diagnostic use across most thyroid pathology laboratories, uses reagents that are commercially validated for in-vitro diagnostic use, and reproduces the BRAF+ progression-free interval stratification."

**Assessment:** The abstract and Discussion sections imply the IHC 3-plex result is more clinically ready than it is. The Methods and Results sections are appropriately caveated. The abstract does not mention "methylation proxy" — it reads as if actual IHC testing was done. This is a major clinical inaccuracy. A clinician reading the abstract would not know that the IHC 3-plex result is from DNA methylation data.

**Required correction:**
- Abstract must say "an in-silico expression proxy based on methylation data for three markers (TG, PAX8, NKX2-1) corresponding to routinely used IHC antibodies" rather than simply "a three-marker immunohistochemistry combination."
- "immediate translational deployment without new assay development" must be softened to "motivates a retrospective IHC-based pilot study" — the deployment step still requires actual IHC quantification.
- The Discussion line noting that the triplet "reproduces the BRAF+ progression-free interval stratification" must clarify that this reproduction is via a methylation proxy, not via actual IHC staining.

---

## 4. GSE151179 — "prediction" vs "concordance" distinction

**What the data show:** GSE151179 is a post-RAI-refractory dataset. The comparison is cross-sectional between post-RAI refractory tumors and non-refractory controls. It is NOT a pretreatment dataset. The manuscript cannot claim the DM1 state predicts who will become RAI-refractory — it can only claim that post-RAI-refractory tumors share transcriptional features with DM1 tumors.

**What the manuscript says:**
- v2 line 124: "indicating that the DM1 programme resembles the transcriptional state observed after clinical RAI failure... we interpret it as a biological anchor rather than as a validated response-prediction model." — This is correct.
- v2 Results §7 (line 122): "the DM1 programme resembles the transcriptional state observed after clinical RAI failure." — Correct.
- The Limitations slot (voice-protected, unfilled) lists "prospective SNUBH/IHC validation pending" — this is inaccurate since the DM1_PROJECT_STATE.md confirms validation is terminated.

**Assessment:** The GSE151179 framing is mostly handled correctly in the body of the text. The critical risk is in the title ("predicts radioiodine-refractoriness") and abstract, where the boundary between retrospective concordance and prospective prediction can blur. The Limitations slot must be filled before submission and must accurately state that internal validation is not proceeding.

---

## 5. Methylation language — "silencing" vs "correlate"

**CLAIM_LANGUAGE_MATRIX.md** requires: "promoter methylation correlated with reduced expression," "epigenetic correlate" — prohibits "methylation causes silencing" or "mechanism proven."

**What the manuscript says:**
- The Results (04_results.md §2.5a, line 78) states: "promoter hypermethylation ... provides a rationale for hypomethylating agents (decitabine, azacitidine) as a candidate epigenetic-targeted RAI re-induction strategy."
- This implies a causal connection (HMA would reverse the silencing). The word "rationale" is on the boundary — it is motivational language, not a causal claim, but it could be read as implying that methylation IS the silencing mechanism.
- The v2 Fig. 3 legend uses "DM1 epigenetically silences thyroid differentiation machinery" — the word "silences" implies causality.
- The Q9 reviewer response (09_reviewer_qa.md lines 47-48) correctly states: "We treat the methylation finding as correlative rather than mechanistic."

**Required correction:**
- Change "DM1 epigenetically silences thyroid differentiation machinery" (Fig 3 title, v2 caption line 246 and reviewer-first plan) to "DM1 is associated with epigenetic correlates of thyroid differentiation gene silencing" or "DM1 shows promoter hypermethylation of thyroid differentiation genes."
- In 04_results.md (line 84), the phrase "provides a rationale for hypomethylating agents" is acceptable if explicitly preceded by "correlative, not proven causal."

---

## 6. Abstract overclaim vs study design

**Issues in the abstract (v2 line 40):**

| Abstract phrase | Problem | Correction |
|---|---|---|
| "positioning the finding for immediate translational deployment without new assay development" | IHC result is methylation proxy; "immediate deployment" exceeds data | "motivates retrospective IHC pilot validation using existing clinical antibodies" |
| "candidate treatment-selection framework" | No treatment comparison; single-cohort interaction | "candidate pretreatment risk stratifier subject to prospective validation" |
| "These findings position the eight-gene axis, and its clinical-grade IHC 3-plex derivative" | "clinical-grade" overstates status of methylation proxy | Remove "clinical-grade" OR specify "with potential for clinical translation via IHC" |
| No mention of retrospective design | NC checklist requires disclosure in abstract | Add "In this retrospective in-silico analysis using public cohorts..." |

---

## 7. Limitations — completeness and honesty

The Limitations section is voice-protected and has not yet been authored. The Limitations slot text (v2 line 154) lists suggested ingredients including "prospective SNUBH/IHC validation pending." Per DM1_PROJECT_STATE.md (dated 2026-07-30), wet-lab validation is definitively terminated. The Limitations must state:

- "Internal validation was not completed. TSO500 RNA-level expression analysis was not feasible on the DNA-oriented TSO500 platform. Routine RNA-seq validation was not attempted. Research IHC antibody reliability concerns precluded a dedicated IHC validation study. Although TG, PAX8, and NKX2-1 are clinically used antibodies, no retrospective IHC validation has been completed at our institution."

Any phrase implying that validation is "pending," "forthcoming," "will be added during review," or "underway" is **inaccurate** and must not appear in the submitted manuscript.

**Other required Limitations disclosures:**
- Retrospective public cohorts only
- TCGA primary PTC has very low event rate (3.2%) limiting precision of OS estimates
- MSK cohort is advanced-disease-enriched and not representative of primary PTC
- Methylation is correlational; no causal proof of silencing mechanism
- n=19 fusion-negative DM1 subgroup for within-DM1 methylation comparison is underpowered
- K2 calibration mismatch documented; within-sample-centered rescaling was required
- Spatial Visium analysis was non-confirmatory and is included only as a caveat
- A2 interaction result is single-cohort only; no external replication available with both BRAF genotype and PFI annotation
- IHC 3-plex result is an in-silico methylation proxy, not a validated pathology test
- GSE151179 is post-RAI-refractory concordance, not prospective pretreatment prediction
