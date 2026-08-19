---
title: Change log — repair pass 2026-07-30
repair_agent: nc-submission-repair-agent
draft_file: NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md
---

# Change log — repair pass 2026-07-30

Format: file | location | old text (first 30 chars) | new text (first 30 chars) | reason | audit_issue_id

---

## REPAIR 1 — Title (FS-03)

| file | location | old text | new text | reason | audit_issue_id |
|---|---|---|---|---|---|
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | YAML title field | "Nature Communications-style full" | "Nature Communications-style full ... — REPAIR PASS 2026-07-30" | Added repair pass marker to YAML | FS-03 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | H1 title line 29 | "A thyroid-lineage state predicts" | "A thyroid-lineage state marks ra" | "predicts" overclaims for retrospective study with no direct RAI response measurement; replaced with "marks" per Option A | FS-03 |

---

## REPAIR 2 — Abstract overclaims (FS-04 and FS-06)

| file | location | old text | new text | reason | audit_issue_id |
|---|---|---|---|---|---|
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Abstract, study design disclosure | "apply it across bulk transcriptom" | "apply it across retrospective pub" | Added "retrospective public cohorts without internal tissue validation" disclosure per NC requirement | FS-04 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Abstract, IHC claim | "Reducing the panel to a three-mar" | "An in-silico expression proxy bas" | Changed to clarify IHC 3-plex is an in-silico methylation proxy, not actual IHC; added "without requiring new assay development" replacing "immediate translational deployment" | FS-04, FS-06 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Abstract, DM1 stratification direction | "such that DM1 stratifies progress" | "such that DM1-low tumours have sh" | Clarified direction of HR=0.66 (higher score = more iodine-handling-high) and removed ambiguity | FS-01, FS-04 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Abstract, treatment-selection | "as a candidate treatment-selectio" | "as a candidate treatment-selectio" | Changed "framework" to "hypothesis" | FS-04 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Abstract, IHC proxy labeling | "its clinical-grade IHC 3-plex der" | "its in-silico IHC 3-plex derivati" | Changed "clinical-grade" to "in-silico" to accurately reflect the analysis type | FS-06 |

---

## REPAIR 3 — Figure architecture (FS-02)

| file | location | old text | new text | reason | audit_issue_id |
|---|---|---|---|---|---|
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Results §1 line ~62, Fig. 7a citation | "(Fig. 7a)" | "[EDITORIAL DECISION: Fig.7/Fig.8..." | Fig. 7a has no asset in nc_main; cited asset is web-figure-only; editorial note added | FS-02 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Results §3 line ~76, Fig. 7b panel A citation | "(Fig. 7b, panel A)" | "[EDITORIAL DECISION: interaction ..." | Fig. 7b not in NC submission assets; editorial note added | FS-02 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Results §3 line ~78, Fig. 7b panel B citation | "(Fig. 7b, panel B)" | "[EDITORIAL DECISION: subgroup HR ..." | Fig. 7b not in NC submission assets; editorial note added | FS-02 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Results §4 line ~88, Fig. 7c citation | "(Fig. 7c)" | "[EDITORIAL DECISION: external rep..." | Fig. 7c not in NC submission assets; editorial note added | FS-02 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Results §5 IHC line ~98, Fig. 8a citation | "(Fig. 8a)" | "(data in Supplementary [EDITORIAL..." | Fig. 8a not in NC submission assets; redirected to Supplementary with editorial note | FS-02 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Results §5 IHC line ~100, Fig. 8b,c citation | "(Fig. 8b,c)" | "(data in Supplementary [EDITORIAL..." | Fig. 8b,c not in NC submission assets; redirected to Supplementary with editorial note | FS-02 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Results §5 IHC, trailing "(Fig. 8c)" | "(Fig. 8c)" | removed (integrated into preceding note) | Duplicate Fig. 8c reference removed | FS-02 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Results §6 multi-omics, Fig. 7d citation | "(Fig. 7d)" | "[EDITORIAL DECISION: multi-omics ..." | Fig. 7d not in NC submission assets; editorial note added | FS-02 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Results §7 simulation, Fig. 8d citation | "(Fig. 8d)" | "[EDITORIAL DECISION: Monte Carlo ..." | Fig. 8d not in NC submission assets; editorial note added | FS-02 |

---

## REPAIR 4 — HR direction clarification (FS-01)

| file | location | old text | new text | reason | audit_issue_id |
|---|---|---|---|---|---|
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Abstract, DM1 HR direction | "such that DM1 stratifies progress" | "such that DM1-low tumours have sh" | Added explicit direction: lower score = lower differentiation = higher hazard = shorter PFI; clarified HR<1 means higher score = better outcome | FS-01 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Results §3, BRAF+ stratified analysis | "Within the BRAF V600E+ subset (n =" | "Within the BRAF V600E+ subset (n =" | Rephrased to "DM1-low tumours had shorter progression-free interval than DM1-high tumours" with parenthetical clarifying score direction (higher = more iodine-handling-high) | FS-01 |

---

## REPAIR 5 — Extended Data → Supplementary Figure

| file | location | old text | new text | reason | audit_issue_id |
|---|---|---|---|---|---|
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | — | — | — | No "Extended Data" found in v2 draft text; repair not needed in this file. Other files (05_figure_captions_NC.md, REVIEWER_FIRST_NC_REBUILD_2026_07_08.md, file names in nc_main/) contain ED terminology — those are out of scope for this repair pass but should be addressed before submission. | FD-02, JC-01 |

---

## REPAIR 6 — Lee 2024 cohort disambiguation (MS-01)

| file | location | old text | new text | reason | audit_issue_id |
|---|---|---|---|---|---|
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Results §4, Lee 2024 opening sentence | "Direct replication of the DM1 axi" | "Direct replication of the DM1 axi" | Added explicit disambiguation: (i) n=632 within-cohort axis separation, Cohen's d=5.93; (ii) n=370 BRAF-interaction direction replication, p=2.7e-7. | MS-01 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Methods, Cohorts and data sources | "GSE213647/Lee 2024 (n = 370 tumou" | "GSE213647/Lee 2024 (full cohort n" | Added n=632 for full cohort and clarified n=370 is driver-annotated subset | MS-01, MS-02 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Methods, External replication | "Direct replication of the DM1 axi" | "Direct replication of the DM1 axi" | Added two-analysis specification: (i) n=632 axis separation, (ii) n=370 BRAF-interaction direction replication | MS-01 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Discussion, Generalizability | "Lee 2024 Mann-Whitney p = 2.7 x 1" | "Lee 2024: within-cohort axis sepa" | Added both analyses with n and statistic | MS-01 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Fig. 7c caption | "Lee 2024 (n = 370, Mann-Whitney p" | "Lee 2024 (n = 370, BRAF-interacti" | Disambiguated to BRAF-interaction direction replication | MS-01 |

---

## REPAIR 7 — GPL570 cohort count (MS-02, RP-04)

| file | location | old text | new text | reason | audit_issue_id |
|---|---|---|---|---|---|
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Methods, Cohorts and data sources | "GPL570 microarray cohorts GSE29265, GSE33630 and GSE65144" | "GPL570 microarray cohorts GSE29265, GSE33630, GSE65144 and GSE126698 (n = 28)" | Added missing fourth GPL570 cohort accession GSE126698 (confirmed from reviewer QA and figure captions SX_v14); resolves three-vs-four mismatch | MS-02, RP-04 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Methods, External replication section | "GPL570 cohorts were re-analysed" | "GPL570 cohorts (GSE29265, GSE33630, GSE65144, GSE126698; total n ≈ 205)" | Listed all four cohort accessions with aggregate n | MS-02, RP-04 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Fig. 7c caption | "three GPL570 cohorts" | "four GPL570 cohorts (GSE29265, GSE33630, GSE65144, GSE126698)" | Corrected count from three to four; listed accessions | MS-02 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Discussion, Generalizability | "(GPL570 4-cohort direction-consist" | "(GPL570 four-cohort direction-cons..." | Added accession list to the generalizability summary | MS-02 |

---

## Additional fixes applied during repair pass

| file | location | old text | new text | reason | audit_issue_id |
|---|---|---|---|---|---|
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Results §5, IHC translational claim | "This positions the IHC 3-plex der" | "This positions the in-silico IHC 3-plex derivative as a hypothesis..." | Replaced overclaim about "immediate translational deployment" in body text; consistent with FS-04/FS-06 fix | FS-04, FS-06 |
| NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md | Discussion, Clinical implications | "This means that the finding can b" | "This means that the in-silico fin" | Added "in-silico" qualifier and "without requiring new assay development" with IHC confirmation caveat | FS-04, FS-06 |
