---
title: Repair summary — nc-submission-repair-agent
repair_date: 2026-07-30
draft_file: NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md
---

# Repair summary

## Repairs completed

### REPAIR 1 — Title (FS-03): DONE
Changed H1 title from "A thyroid-lineage state **predicts** radioiodine-refractoriness..." to "A thyroid-lineage state **marks** radioiodine-refractory biology..." (Option A). "Predicts" was dropped because no direct pretreatment RAI response measurement exists. YAML header annotated with repair pass date.

### REPAIR 2 — Abstract overclaims (FS-04 and FS-06): DONE
Four changes made to the abstract:
- (a) Added study design disclosure: "across retrospective public cohorts without internal tissue validation"
- (b) IHC 3-plex sentence now reads: "An in-silico expression proxy based on promoter methylation of a three-marker combination...without requiring new assay development, although direct immunohistochemistry validation on tissue sections has not been performed"
- (c) "candidate treatment-selection framework" changed to "candidate treatment-selection hypothesis"
- (d) "clinical-grade IHC 3-plex derivative" changed to "in-silico IHC 3-plex derivative"

### REPAIR 3 — Figure architecture (FS-02): DONE (editorial notes placed)
All nine Fig. 7 and Fig. 8 citations in the Results text have been annotated with [EDITORIAL DECISION: ...] notes that:
- Name the current web-figure-only asset location
- Direct to either "promote to Supplementary Figure" or "rebuild as main figure before submission"
- Where possible (Fig. 8 panels), text is rephrased to "data in Supplementary" pending the editorial decision

Note: The underlying figures exist as web assets in `dm1_story_web/public/figures/` but are NOT in `manuscript_v8_nc_main/` and have no captions in `05_figure_captions_NC.md`. A definitive 6-figure-plan editorial decision is still required before submission. This repair pass marks every orphan citation clearly so the author cannot miss them.

### REPAIR 4 — HR direction clarification (FS-01): DONE
Two locations fixed:
- Abstract: "DM1 stratifies progression-free interval within BRAF+ tumours (HR = 0.66...)" replaced with explicit direction: "DM1-low tumours have shorter progression-free interval than DM1-high tumours within the BRAF+ subset (continuous thyroid differentiation score HR per +1 SD = 0.66, where higher score = more iodine-handling-high...)"
- Results §3: "DM1 score was significantly associated with progression-free interval (HR per +1 SD = 0.66...)" replaced with "DM1-low tumours had shorter progression-free interval than DM1-high tumours (continuous score HR per +1 SD = 0.66, 95% CI 0.48-0.92, p = 0.013; higher score = more differentiated, i.e., eight-gene thyroid differentiation score where higher = more iodine-handling-high)"

### REPAIR 5 — Extended Data → Supplementary Figure: NOT NEEDED IN v2 DRAFT
The v2 draft (NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md) contains no "Extended Data" occurrences. This repair was a no-op for the target file. However, `05_figure_captions_NC.md`, `REVIEWER_FIRST_NC_REBUILD_2026_07_08.md`, and asset files in `manuscript_v8_nc_main/` (ED1-ED15 named files) still need this correction — these are outside scope of this repair pass.

### REPAIR 6 — Lee 2024 cohort disambiguation (MS-01): DONE
Four locations updated to distinguish:
- Within-cohort DM1/DM2 axis separation: n=632, Cohen's d=5.93
- BRAF-interaction direction replication: n=370, Mann-Whitney p=2.7e-7, d=0.24
Updated in: Results §4 opening paragraph, Methods Cohorts section, Methods External Replication section, Discussion Generalizability paragraph, Fig. 7c caption.

### REPAIR 7 — GPL570 cohort count (MS-02, RP-04): DONE
Added GSE126698 (n=28) as the fourth GPL570 cohort across three locations:
- Methods Cohorts section: now lists "GSE29265, GSE33630, GSE65144 and GSE126698 (n = 28)"
- Methods External Replication: now lists all four accessions with total n ≈ 205
- Fig. 7c caption: changed "three GPL570 cohorts" to "four GPL570 cohorts (GSE29265, GSE33630, GSE65144, GSE126698)"
- Discussion Generalizability: added all four accession numbers

**Source confirmation:** GSE126698 (n=28) confirmed from reviewer_qa.md A14, 05_figure_captions_NC.md Fig 3D entry, and AUDIT_LOCKED_RESULTS.md V-08 which lists "four cohorts, total n=205."

---

## Additional fixes applied (not in original repair list)

Two body-text sentences carried the same "immediate translational deployment" overclaim as the abstract (FS-04/FS-06). Both were fixed for consistency:
- Results §5 IHC section closing sentence: "immediate translational deployment" → "hypothesis that can be tested...pending prospective validation with direct H-score quantification"
- Discussion Clinical implications: "existing pathology workflows" paragraph: added "in-silico finding" qualifier and "without requiring new assay development" phrasing with IHC confirmation caveat

---

## Repairs deferred

| Repair | Reason for deferral |
|---|---|
| FS-05: Limitations section (voice-protected) | Content is entirely author-protected; cannot be filled by this agent. Author must draft. |
| AI-01 through AI-06: All voice-protected slots | Protected sections left exactly as-is. |
| MS-03: "80/80 cells" denominator | Requires author or analyst to enumerate the exact cohort × gene matrix; not safe to guess |
| MS-04: "19 external validation entries" vs "14 forest entries" | Requires a Supplementary Table to be built; content decision needed |
| MS-05: "MAR confirmed" in 04_results.md | 04_results.md is a legacy section file, not the v2 canonical draft; out of scope for this pass |
| MS-06: Fig 3 caption "epigenetically silences" causal language | In 05_figure_captions_NC.md, not in the v2 draft; out of scope |
| MS-07: GSE151179 approximate values (d≈-1.0, p≈10^-4) | Exact values require analyst to re-run or retrieve output; values not changed |
| FD-01: Fig 3C missing | Requires figure build; out of scope for text repair |
| FD-02: ED1-ED15 file rename | Requires file system operation on nc_main/ assets; out of scope for this pass |
| FD-03: Source data manifest | Requires Supplementary Table build |
| RP-01 through RP-03: KMeans seed, HM450 probe aggregation, ATA risk-tier | Methods additions requiring author/analyst input |
| RP-05: Mun 2025 accession | Author must obtain ProteomeXchange PXD |
| RP-06: Zenodo DOI placeholders | Admin task; out of scope |
| JC-02 through JC-06: Life Sciences Reporting Summary, Stats section, Data availability accessions, Author contributions, Competing interests | Require author action |

---

## Remaining known blockers after repair pass

### Still fatal before submission:
1. **FS-02 (Figure architecture)**: Fig. 7 and Fig. 8 have been annotated with [EDITORIAL DECISION] notes but the underlying architectural decision (5-fig vs 6-fig vs 8-fig) has NOT been made. Every [EDITORIAL DECISION] note in the Results must be resolved before submission.
2. **FS-05 (Limitations section)**: Voice-protected placeholder still unfilled. No actual limitations text exists.
3. **FD-01 (Fig 3C missing)**: The methylation-expression scatter plot still does not exist.
4. **AI-01 through AI-06 (Voice-protected sections)**: Introduction hook, Introduction aim, Discussion §3.1, Limitations, Cover letter ¶1, Reviewer Q9 — all unfilled.

### Major but repairable:
5. **MS-07**: GSE151179 approximate values (d≈-1.0, p≈10^-4) remain approximate in body text.
6. **FD-02**: ED1-ED15 file names need renaming to SFig1-SFig15; the v2 draft itself is clean but other files and assets are not.
7. **RP-01/RP-02/RP-03**: KMeans seed, HM450 probe aggregation, ATA tier operationalization absent from Methods.

---

## Recommendation: re-run hostile reviewer?

**Yes — but only after the author resolves FS-02 (figure architecture) and FS-05 (Limitations).**

Rationale: The repair pass eliminated FS-03 and FS-04/FS-06 language-level overclaims and added disambiguating language for FS-01 direction. However, a hostile reviewer will immediately see:
- The [EDITORIAL DECISION] notes where Fig. 7 and Fig. 8 were cited — these are unfilled gaps that no reviewer should encounter in a submission
- The empty Limitations placeholder
- The still-approximate GSE151179 values

A re-run hostile review at the current state would produce the same FS-02 and FS-05 fatals. Once the author fills the Limitations slot and decides on figure architecture (6-figure plan is strongly recommended), the text will be ready for a productive second hostile review.
