# RAI claims registry — P0/P1 summary

Date: 2026-08-06
Full registry: `06_CLAIMS_REGISTRY.csv` (101 rows: 15 P0 fatal, 8 P1 major, 13 P2 moderate, 65 P3 minor)

## Which manuscript is current

**`project/manuscript_v8/NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md`** is the current manuscript.

Evidence:
- `project/manuscript_v8/STATUS_INSILICO_FINALIZATION_2026_07_30.md` (dated 2026-07-30, the newest status note in the folder) states explicitly: "가장 최신 draft: `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md`" ("most recent draft").
- Its own front matter is stamped `REPAIR PASS 2026-07-30`.
- `project/manuscript_v8/CLAUDE.md` (the project constitution) lists it as canonical file #6 and instructs: "Do not silently substitute an older draft."
- `08_cover_letter.md` and `09_reviewer_qa.md` were both edited the same day (2026-07-30) as the NC v2 repair pass, confirming they are its companion submission documents.
- By contrast, `10_full_manuscript_compiled.md` is stamped `date: 2026-05-13`, targets `Cell Reports Medicine (1순위)` (an earlier, superseded venue plan), and was last regenerated 2026-05-13 — before the NC v2 upgrade. It is superseded but still present on disk and was audited per the task instruction ("and `10_full_manuscript_compiled.md` if separate").

Also audited as submission-critical companion files (evidence-lock rule covers "manuscript, caption, cover-letter, reviewer-response, and submission files"): `08_cover_letter.md`, `09_reviewer_qa.md`, `05_figure_captions_NC.md`.

## Ground truth used to judge claims (supplied, not re-derived)

1. TCGA-THCA **does** contain radioiodine dose (mCi) and per-course response (GDC BCR Biotab radiation file): 235 RAI-treated patients have a panel score, 167 have an evaluable response.
2. The panel shows **no** detectable association with initial recorded RAI response (Cohen's d = −0.03 to −0.00, P = 0.73–0.91 across five index-course rules).
3. The panel **is** associated with post-treatment structural events *within RAI-treated patients*: Firth OR 0.221 (0.062–0.784), P = 0.019 for new tumour event; OR 0.224 (0.061–0.826), P = 0.025 for persistent disease. Prognostic-within-treated-population only.
4. Treatment × score interaction is **not estimable** (2 events among non-RAI patients, P = 0.077) → no predictive/treatment-selection claim is permissible.
5. Driver class (BRAF/RAS status) does **not** consistently discriminate refractoriness across 3 cohorts (n=294; pooled OR 0.75, 95% CI 0.30–1.87, P = 0.53, I² = 63%); driver-negative compartment sits at 46–48%, near the cohort mean.
6. Zhang 2026 proteomic subtype = orthogonal conceptual support only, **not** validation of the 8-gene panel (not currently cited anywhere in the audited files).

## P0 fatal claims (15)

### A — "No cohort contains linked RAI outcome data" (contradicts ground truth #1)

| Claim ID | File:line | Text | Fix |
|---|---|---|---|
| NCv2-041 | `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md:154` | "The principal limitation is that we did not have a pretreatment cohort in which the DM1 state could be tested directly against measured radioiodine uptake and clinical response." | Replace with the actual TCGA-THCA RAI dose/response finding: n=235 scored / 167 evaluable; no association with initial response (d=−0.03 to −0.00, P=0.73–0.91); prognostic within RAI-treated patients (OR 0.221/0.224); treatment×score interaction not estimable (P=0.077). |
| NCv2-045 | `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md:154` | "The most direct next study is therefore a prespecified retrospective FFPE validation linked to administered radioiodine dose, post-therapy uptake and response, followed by independent prospective testing." | State that TCGA-THCA dose/response data already exist and were analyzed; reserve the FFPE/prospective study for independent replication, not first access to linked data. |
| COMP-024 | `10_full_manuscript_compiled.md:511` (source `06_discussion.md:47`) | "Sixth, the panel name explicitly references RAI-responsiveness biology, but no cohort in this study contains prospective post-thyroidectomy RAI outcome data linked to per-sample DM1 calls..." | Same fix as NCv2-041; also fix the "RAI-responsiveness" name (see group B). |

### B — "RAI-responsiveness panel" naming (contradicts ground truth #2: no association with RAI response)

The panel's own name asserts it measures/predicts RAI responsiveness. It does not appear in the current NC v2 main text (already corrected there to "thyroid differentiation and iodine-handling axis"), but survives, unfixed, in several companion documents:

| Claim ID | File:line | Text |
|---|---|---|
| COMP-003 | `10_full_manuscript_compiled.md:67` (`01_abstract.md:28`) | "We applied an 8-gene RAI-responsiveness panel ... to TCGA-THCA (n=504), MSK-IMPACT (n=117), and Korean cohorts (n=865)." |
| COMP-009 | `10_full_manuscript_compiled.md:158` (`03_introduction.md:37`) | "We applied an 8-gene RAI-responsiveness panel — independently selected from canonical thyroid differentiation biology..." |
| COMP-017 | `10_full_manuscript_compiled.md:261` | "The 8-gene RAI-responsiveness panel resolves a DM1/DM2 cluster within papillary thyroid carcinoma." (Figure 1 title) |
| **CL-004** | `08_cover_letter.md:21` | "In this study, we address that gap using an 8-gene RAI-responsiveness panel grounded in canonical thyroid differentiation biology and **validated** across TCGA-THCA, MSK-IMPACT thyroid cancer, Korean cohorts..." — **this is in the current, 2026-07-30-dated cover letter**, i.e. it is stale boilerplate that was not updated when the cover letter was otherwise revised alongside NC v2 on the same day. |
| FIGNC-001 | `05_figure_captions_NC.md:21` | "...to the final 8-gene RAI-responsiveness panel: SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1." |

Fix (all rows): rename to "thyroid differentiation and iodine-handling panel" throughout, matching the terminology NC v2's main body already uses; in CL-004 also replace "validated" with "reproduced" per `SUBMISSION_CONTROL/CLAIM_LANGUAGE_MATRIX.md`'s own rule.

### C — "Predictive interaction" / "treatment-selection biomarker" framing (contradicts ground truth #4: treatment×score interaction is not estimable)

The manuscript tests a DM1 × **BRAF-genotype** interaction on progression-free interval (p=0.022) — a real, internally supported result. But Results §3, Discussion "Clinical implications," and Figure 7 repeatedly relabel this genotype interaction as a "predictive interaction" and elevate it to "treatment-selection biomarker" status. This conflates genotype-context prognosis with radioiodine-treatment prediction, which the project's own `CLAUDE.md` warns against ("A genotype interaction on PFI is not automatically a treatment-predictive biomarker") — and which today's audit now directly falsifies: the actual DM1-score × radioiodine-treatment interaction is not estimable (P=0.077, only 2 events among non-RAI patients).

| Claim ID | File:line | Text |
|---|---|---|
| NCv2-023 | `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md:74` | Section header: "### A predictive interaction between DM1 and BRAF V600E status" |
| NCv2-024 | `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md:76` | "...whether the DM1 axis functions as a driver-agnostic prognostic biomarker or as a **predictive biomarker** whose effect depends on driver context..." |
| NCv2-026 | `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md:78` | "...This pattern is consistent with **a treatment-selection biomarker** whose informative range is confined to a specific molecular subset, rather than with a driver-independent prognostic score." |
| NCv2-027 | `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md:80` | "...lineage-programme silencing **predicts radioiodine failure**... The DM1 axis offers **a pre-treatment classifier** for this subset..." |
| NCv2-037 | `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md:142` | "...the axis **functions as a candidate treatment-selection biomarker** for a defined driver subset." (single strongest overclaim in the manuscript) |
| NCv2-010 | `NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md:40` (Abstract) | "...position the eight-gene axis...as **a candidate treatment-selection hypothesis** for the BRAF-mutant subset..." |

Fix (all rows): relabel as a "genotype-context-dependent prognostic" finding; explicitly state that a radioiodine treatment-selection role requires a treatment × score interaction, which was tested in TCGA-THCA and found not estimable, and that initial RAI response itself is not associated with the score.

### D — Overreach on the GSE151179 anchor (goes beyond ground truth #3's prognostic-only scope)

| Claim ID | File:line | Text |
|---|---|---|
| FIGNC-005 | `05_figure_captions_NC.md:125` | "The post-RAI dedifferentiation profile recapitulates the DM1 transcriptional state, providing an external clinical anchor that the DM1 axis **names the same biology that defines RAI-refractory progression**." |

This caption asserts biological identity/causation from a cross-sectional group comparison, going beyond the (correctly hedged) main-text language in NC v2 itself ("we interpret it as a biological anchor rather than as a validated response-prediction model," line 124). Fix: match the main-text hedge.

## P1 major claims (8) — brief list

- **NCv2-007** (Abstract, line 40): headline PFI-interaction sentence needs an explicit clause distinguishing driver-genotype interaction from radioiodine-treatment interaction (not estimable).
- **NCv2-015** (Introduction, line 50): study-aim sentence claims the axis was tested for "treatment-selection biomarker status" — only genotype interaction was tested.
- **NCv2-034** (Results, line 126): "treatment-selection framing" label applied to the BRAF-genotype PFI interaction.
- **NCv2-038** (Discussion, line 142): conditional "candidate rule for identifying BRAF+ patients for whom repeated high-dose radioiodine is likely to be futile" skips the missing treatment-interaction precondition.
- **NCv2-040** (Discussion, line 150): "A2 predictive interaction result" — same mislabel as group C, applied to external-replication status.
- **NCv2-042** (Limitations, line 154): non-responder-identification disclaimer needs updating once the TCGA RAI-response/structural-event data (ground truth #1–#4) are incorporated.
- **NCv2-054** (Figure 7 title, line 260): "DM1 predictive interaction with BRAF V600E status" — same mislabel in a figure title.
- **FIGNC-003** (Fig 3 caption, line 121): "supporting DM1 as an orthogonal molecular axis" overreaches from an ATA-tier distribution finding to implied RAI-decision utility.

## Notable non-issues (kept for contrast)

- NC v2 lines 124 (biological-anchor hedge for GSE151179) and 82/94 (single-cohort caveat on the BRAF interaction) already use the correct, defensible hedged language and should be used as the template for fixing the P0/P1 rows above.
- Ground truth #5 (driver class does not consistently discriminate refractoriness) is not directly contradicted anywhere in the current text — the manuscript never claims BRAF/RAS-negative status alone is intrinsically more RAI-refractory. Motivational sentences (e.g. COMP-002, COMP-006) that say the driver-negative compartment "complicates radioiodine decisions" are defensible as background but were flagged P2 with an optional strengthening (cite pooled OR 0.75, P=0.53 as evidence driver status alone is insufficient — which is precisely why a lineage-based panel is needed).
- Ground truth #6 (Zhang 2026) is not cited anywhere in the audited files, so there is no existing claim to correct; noted so it is not introduced as a "validation" claim later.

## Files delivered

- `/home/seungho/personal/THCA_data_analysis/audit/rai_integration_20260806/06_CLAIMS_REGISTRY.csv` — full 101-row registry, all RAI/refractory/predictive/prognostic/survival/etc. sentences found in the current manuscript and its companion submission documents.
- `/home/seungho/personal/THCA_data_analysis/audit/rai_integration_20260806/06_CLAIMS_SUMMARY.md` — this file.

No manuscript files were edited (per file-safety rule); this is an audit-only deliverable.
