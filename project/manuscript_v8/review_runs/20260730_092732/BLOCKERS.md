---
title: Submission blockers
audit_date: 2026-07-30
auditor: manuscript-audit-agent
---

# BLOCKERS — Grouped by category

---

## FATAL SCIENTIFIC (must fix before submission)

### FS-01: HR direction unverified for S-04 and OS pooled estimate

The survival analysis reference direction (which state is the test group, which is the reference) has NOT been confirmed from source code output. Claim S-04 (BRAF+ PFI HR=0.66 per +1 SD DM1 score) is particularly ambiguous: HR<1 per +1 SD of iodine-handling-low score would mean MORE lineage loss → LOWER hazard of progression, which is biologically counterintuitive. This could indicate a continuous score direction issue (higher score = more intact lineage, not more loss). Until confirmed from lifelines model output, the biological narrative cannot be finalized.

**Files:** AUDIT_LOCKED_RESULTS.md rows S-01 to S-04; v2 Results lines 72-78; SUBMISSION_CONTROL/DM1_PROJECT_STATE.md item 7.

**Action required:** Author must run Cox model summary for both OS and PFI, confirm which group is reference, and confirm HR direction is correctly labeled in every figure and text location.

---

### FS-02: Figure architecture unresolved (5-fig vs 6-fig vs 8-fig)

Three parallel figure architectures exist simultaneously:
- v2 draft: 8 main figures (Fig 1-8), with Fig 7 and Fig 8 in a separate web figure directory
- NC captions file (05_figure_captions_NC.md): 6 main figures (Fig 1-6)
- Reviewer-first plan (REVIEWER_FIRST_NC_REBUILD_2026_07_08.md): 5 main figures

No editorial decision has been documented selecting one architecture. This means the current manuscript body (8 figures), the current captions file (6 figures), and the cited Fig 7/Fig 8 panels in the draft text are internally inconsistent. Reviewers would receive a manuscript where Fig 7 and Fig 8 references in the text have no corresponding captions file.

**Files:** NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md, 05_figure_captions_NC.md, REVIEWER_FIRST_NC_REBUILD_2026_07_08.md.

**Action required:** EDITORIAL DECISION REQUIRED. Choose one architecture. Update all captions, figure files, and text references to match. Favored recommendation per reviewer-first plan: 5 main figures + Supplementary Figures.

---

### FS-03: Title overclaims ("predicts radioiodine-refractoriness")

The working title "A thyroid-lineage state predicts radioiodine-refractoriness in BRAF V600E-mutant papillary thyroid cancer" uses "predicts" for a retrospective study without direct RAI response measurement. The CLAIM_LANGUAGE_MATRIX.md title rule flags this. SUBMISSION_CONTROL/DM1_PROJECT_STATE.md item 3 explicitly lists this as a known conflict.

**Action required:** AUTHOR DECISION REQUIRED. The CLAIM_LANGUAGE_MATRIX.md offers three alternatives:
1. "marks radioiodine-refractory biology"
2. "stratifies outcomes and aligns with radioiodine-refractory PTC"
3. "defines an adverse BRAF-context in PTC"

---

### FS-04: Abstract overclaims that exceed study design

Three specific abstract phrases are overclaims (see 03_clinical_review.md and 01_claim_evidence_ledger.tsv):

1. "immediate translational deployment without new assay development" — the IHC result is a methylation proxy; actual IHC testing has not been done.
2. "candidate treatment-selection framework" — no treatment comparison exists; this is at most a hypothesis-generating prognostic interaction.
3. The abstract does not disclose the retrospective in-silico study design. NC checklist requires this disclosure.

**Action required:** Rewrite abstract to (a) disclose retrospective in-silico design, (b) replace "immediate translational deployment" with "motivates retrospective IHC pilot validation," (c) qualify "treatment-selection" as "hypothesis-generating" or "candidate."

---

### FS-05: Limitations section is entirely absent (voice-protected, unfilled)

The Limitations section in the v2 draft (line 154) is a voice-protected placeholder slot. No actual limitations text exists. NC requires a complete Limitations section. Critically, the Limitations must explicitly state that internal wet-lab validation is NOT proceeding (wet-lab termination per DM1_PROJECT_STATE.md dated 2026-07-30), and must NOT use phrases like "SNUBH/IHC validation pending."

**Action required:** AUTHOR INPUT REQUIRED. See VOICE_PROTECTED_SLOTS.md V4 for required elements.

---

### FS-06: IHC 3-plex presented as clinical-ready in abstract despite being methylation proxy

The abstract says "three-marker immunohistochemistry combination (TG + PAX8 + NKX2-1) — already in routine clinical use for thyroid pathology — reproduces this stratification." This implies actual IHC testing was done. The body correctly states "in-silico proxy" but the abstract does not.

**Action required:** Abstract must add "in-silico expression proxy based on promoter methylation" before the IHC 3-plex description.

---

## MAJOR SCIENTIFIC (must fix or explicitly limit)

### MS-01: Lee 2024 cohort n discrepancy (n=370 vs n=632)

Two Lee 2024 analyses produce different n values (370 vs 632) and are not clearly distinguished. Status file lists both separately: main validation at n=632 d=5.93, and A2-R replication at n=370 p=2.7e-7.

**Action:** Clearly label these as distinct analyses with distinct denominators and methodological differences in both text and Methods. Current v2 Methods (line 160) states n=370 without explaining the full cohort size.

---

### MS-02: GPL570 cohort count mismatch (3 vs 4 in v2 Methods)

v2 text says "four GPL570 microarray cohorts" (lines 88, 148) but v2 Methods (line 160) lists only three accessions (GSE29265, GSE33630, GSE65144). The fourth cohort accession is missing.

**Action:** Identify and add the fourth GPL570 cohort accession.

---

### MS-03: "80/80 cells" denominator undefined

"80 of 80 cells moving in the expected direction" (v2 line 86) — the denominator "cell" is undefined. With 8 genes × 10 cohorts = 80, there must be exactly 10 cohorts counted. But neither the v2 draft text nor the Methods lists these 10 cohorts explicitly.

**Action:** State explicitly which cohorts × which comparison constitute the 80 cells, or replace with "all gene-by-cohort comparisons in 8 cohorts" (which would be 64 if 8 genes × 8 cohorts).

---

### MS-04: "19 external validation entries" vs "14 forest entries"

v2 Results (line 86) says "19 external validation entries." AUDIT_LOCKED_RESULTS.md V-01 says "14" master forest entries. The manuscript must reconcile these by defining what constitutes an "entry" in each count.

**Action:** Provide a single list of all external validation entries in Supplementary Table format. Distinguish between "forest entries" (with effect size) and "validation entries" (broader category including qualitative direction checks).

---

### MS-05: "MAR confirmed" language

04_results.md (line 66) says "consistent with missing-at-random (MAR)." CLAIM_LANGUAGE_MATRIX.md prohibits "MAR confirmed." A non-significant chi-square (p=0.56) does not prove MAR.

**Action:** Replace with "no evidence that SV availability differed by DM cluster status."

---

### MS-06: Methylation "silences" causal language

v2 Figure 3 description (NC captions file, line 57) uses "DM1 epigenetically silences thyroid differentiation machinery" as the figure title. CLAIM_LANGUAGE_MATRIX.md prohibits "methylation causes silencing."

**Action:** Change to "DM1 is associated with promoter hypermethylation of thyroid differentiation genes" or "epigenetic correlate of lineage-gene silencing in DM1."

---

### MS-07: GSE151179 approximate values

The v2 text uses "d ≈ -1.0" and "p ≈ 10^-4" for GSE151179. Approximate values should not remain in the final submitted manuscript.

**Action:** Replace with exact values from the GSE151179 analysis output.

---

## MANUSCRIPT CONSISTENCY (cross-file conflicts)

### MC-01: v2 draft vs older section files

The v2 draft is the authoritative document. However, the older section files (01_abstract.md, 04_results.md, 05_figure_captions.md/NC.md, 06_discussion.md) were written for a different figure architecture and contain different claims, different figure numbers, and different citation formats. If these files are used in submission, they will introduce inconsistencies.

**Action:** The v2 draft should be the sole source for submission. The section files should be archived or explicitly synchronized with v2 before compilation.

---

### MC-02: SNUBH validation language

Several files (v2 Fig 8 caption, Limitations slot, STATUS file) reference "SNUBH" validation as prospective, future, or pending. DM1_PROJECT_STATE.md (dated 2026-07-30) states wet-lab validation is terminated. Any reference to SNUBH/IHC/RNA validation as "pending" or "forthcoming" must be corrected.

**Action:** Search and replace all instances of "SNUBH validation pending," "internal validation underway," or similar phrases.

---

### MC-03: STAR methods vs v2 Methods section

07_star_methods.md is written in Cell Press STAR Methods format (Key Resources Table, STAR sections). The v2 draft uses a journal-specific NC Methods section. These are different formats for the same content. Only one should be submitted to NC.

**Action:** The v2 Methods section is formatted for NC. The STAR Methods file should be archived or converted. Confirm which format NC accepts for the Methods.

---

### MC-04: K2 n=260 in v2 Methods vs n=235 in STAR Methods

v2 Methods (line 160): "PRJEB11591/K2 (n = 260)." STAR Methods (line 59): "K2 / PRJEB11591 (n = 235; primary Korean PTC, post-QC)." These differ: n=260 appears to be pre-QC, n=235 is post-QC.

**Action:** Standardize to post-QC n=235 throughout, or explain the distinction (total=260, analyzed=235) where the larger number is used.

---

## FIGURES AND SOURCE DATA (missing assets, numbering conflicts)

### FD-01: Fig 3C missing (FATAL for figure completeness)

The beta-expression scatter plot for TPO, DIO1, TSHR, TG (Fig 3C per all figure architectures) does not exist. This panel is described in v2 Fig 3 caption (line 246) and the reviewer-first plan identifies it as "the missing mechanism bridge panel."

**Action:** Build Fig 3C from TCGA HM450 beta and RNA-seq paired data (n=503). Script template: `fig3c_beta_expression_scatter.py` noted as needed in STATUS file.

---

### FD-02: Extended Data → Supplementary Figure rename required

15 files named ED1-ED15_*.png in `project/results/manuscript_v8_nc_main/` use "Extended Data" terminology, which is prohibited for Nature Communications.

**Action:** Rename all ED files to SFig format. Update all references in the manuscript.

---

### FD-03: No source data manifest exists

No document maps each figure panel to a specific source data file and generating script.

**Action:** Create a Supplementary Table listing: figure number, panel letter, source file path, generating script, key analysis parameters.

---

### FD-04: Fig 7 and Fig 8 exist only in web figure directory

The v2 draft references Fig 7 and Fig 8 extensively, but the assets live in `project/dm1_story_web/public/figures/` (web resolution/format), not in the NC submission figure directory.

**Action:** Copy or regenerate Fig 7 and Fig 8 panels at 180mm / 600 dpi PDF for NC submission. OR: resolve the figure architecture (FS-02) and reframe this content as Supplementary Figures.

---

## REPRODUCIBILITY (methods gaps)

### RP-01: KMeans random seed not specified

**Action:** Specify `random_state=42` (or actual seed used) in Methods and code.

### RP-02: HM450 probe aggregation not specified

**Action:** State which probes (TSS200, 1stExon, body?) and aggregation method (mean? max? specific probe set?) were used for each gene's "promoter beta."

### RP-03: ATA risk-tier operationalization absent

**Action:** State which TCGA clinical fields were used to assign ATA 2015 low/intermediate/high tiers to TCGA-THCA samples.

### RP-04: GSE126698 absent from v2 Methods cohort list

**Action:** Add GSE126698 to the Methods cohort table with n and source citation.

### RP-05: Mun 2025 proteogenomic accession absent

**Action:** Add public accession number (ProteomeXchange PXD or equivalent) for Mun 2025 dataset.

### RP-06: Two Zenodo/repo TODO placeholders remain

**Action:** Create public repository before submission; obtain Zenodo DOI; replace placeholders.

---

## JOURNAL COMPLIANCE (NC formatting requirements)

### JC-01: "Extended Data" terminology throughout

NC uses "Supplementary Figure" not "Extended Data." Must be corrected in all files.

### JC-02: No Life Sciences Reporting Summary

NC requires a completed Life Sciences Reporting Summary for biological studies. Not mentioned in any manuscript file.

**Action:** Download NC Life Sciences Reporting Summary form; complete; attach to submission.

### JC-03: Statistics and Reproducibility subsection in Methods

NC requires a dedicated "Statistics and Reproducibility" subsection within Methods. The STAR Methods format has a "Quantification and Statistical Analysis" section but NC does not use STAR format. The v2 Methods section (lines 210-213) has a brief "Statistical analysis" section; it needs to be formatted as a dedicated NC-required subsection.

### JC-04: Data availability statement must include accession numbers in body

The current v2 Data Availability statement (line 216) says accessions are "listed in Supplementary Table 1." NC prefers that key accessions appear directly in the Data Availability statement, not only in a Supplementary Table.

### JC-05: Author Contributions statement incomplete

v2 line 228: "H.-W.Y. supervised the clinical framing and interpretation. [VERIFY AND EXPAND BEFORE SUBMISSION.]"

**Action:** Author must expand and finalize all contributions per NC CRediT taxonomy.

### JC-06: Competing Interests and Acknowledgements

v2 line 232: "The authors declare no competing interests. [VERIFY BEFORE SUBMISSION.]"
v2 line 224: Acknowledgements slot unfilled.

**Action:** Author must verify and complete both.

---

## AUTHOR INPUT REQUIRED (voice sections, admin)

### AI-01: Introduction opening hook (V1) — UNFILLED

The clinical hook paragraph is absent. A manuscript cannot be submitted with this slot empty.

### AI-02: Introduction final aim paragraph (V2) — UNFILLED

The study aims paragraph is absent.

### AI-03: Discussion mechanism interpretation §3.1 (V3) — UNFILLED

The mechanistic Discussion section is absent. This is the most scientifically sensitive voice section.

### AI-04: Limitations paragraph (V4) — UNFILLED + CONTENT CONSTRAINT

Must explicitly state wet-lab validation termination. See FS-05.

### AI-05: Cover letter paragraph 1 (V5) — UNFILLED

Cannot submit without a cover letter.

### AI-06: Reviewer Q9 prose (V6) — UNFILLED

Q9 is voice-protected; author must draft.

---

## SUBMISSION ADMINISTRATION

### SA-01: Institutional corresponding email unverified

v2 line 36: "[VERIFY institutional corresponding email]"

**Action:** Confirm institutional email for the corresponding author (Hyeong-Won Yu at SNUBH). The personal email (kukshomr@gmail.com) may not be accepted by NC for correspondence.

### SA-02: Zenodo DOI not obtained

Two placeholder TODO items remain in STAR Methods for the Zenodo archive DOI.

### SA-03: Public code repository not established

The code repository URL placeholder in STAR Methods is empty.

### SA-04: Author affiliations need final verification

v2 line 31-32: Two affiliations listed for both authors. Confirm current institutional affiliation and appointment status.

### SA-05: Funding acknowledgement absent

Acknowledgements section (v2 line 224) is an author-placeholder slot. NC requires funding statement including grant numbers.
