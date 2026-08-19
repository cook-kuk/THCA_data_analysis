# DM1 manuscript project constitution

## Mission

Bring Paper 1 (DM1) to an honest, internally consistent, reproducible, submission-ready Nature Communications package. Act as a managing editor, not a fluent autocomplete system.

## Canonical project files

Before editing manuscript content, read:

1. `manuscript_v8/SUBMISSION_CONTROL/DM1_PROJECT_STATE.md`
2. `manuscript_v8/SUBMISSION_CONTROL/AUDIT_LOCKED_RESULTS.md`
3. `manuscript_v8/SUBMISSION_CONTROL/CLAIM_LANGUAGE_MATRIX.md`
4. `manuscript_v8/SUBMISSION_CONTROL/VOICE_PROTECTED_SLOTS.md`
5. `manuscript_v8/SUBMISSION_CONTROL/NATURE_COMMUNICATIONS_CHECKLIST.md`
6. `manuscript_v8/NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md`
7. `manuscript_v8/REVIEWER_FIRST_NC_REBUILD_2026_07_08.md`

If a listed file is absent, record it as a blocker. Do not silently substitute an older draft.

## Non-negotiable scientific boundaries

- This is a retrospective public-cohort in-silico study with no completed internal wet-lab validation.
- No analysis directly links pretreatment RNA to measured radioiodine uptake.
- GSE151179 is a post-RAI-refractory transcriptional alignment and biological anchor, not prospective prediction.
- Methylation findings are correlational. Never use causal verbs unless a causal experiment exists.
- The three-marker IHC result is an in-silico expression proxy using clinically available antibodies, not a validated pathology assay.
- A genotype interaction on PFI is not automatically a treatment-predictive biomarker.
- Never imply that high-dose RAI should be withheld in clinical practice based on this study alone.
- The defensible clinical frame is hypothesis-generating harm avoidance: identify a candidate lineage-low state that may support future pretreatment validation and reduce futile repeat RAI after validation.

## Evidence discipline

- Every quantitative manuscript claim must map to an analysis output, table, script, or immutable audit record.
- Never invent or estimate missing n, confidence intervals, p-values, software versions, accession numbers, DOI, PMID, IRB numbers, or author affiliations.
- A numeric value may be locked while its semantic interpretation remains unverified. Check event coding, comparator, reference group, transformation, and model formula before wording it.
- If two files disagree, preserve both, identify the conflict, and resolve from primary output or code. Do not choose the more attractive value.
- Use `[AUTHOR CONFIRMATION REQUIRED]`, `[SOURCE OUTPUT REQUIRED]`, or `[EDITORIAL DECISION REQUIRED]` rather than guessing.

## Editing protocol

- Audit before repair.
- Make a timestamped review run under `manuscript_v8/review_runs/`.
- Before editing, record git branch, HEAD, working-tree status, target draft checksum, and figure manifest.
- Never overwrite raw data, original analysis outputs, or the latest manuscript without a recoverable diff.
- Prefer targeted edits. Preserve a change log with old wording, new wording, reason, and evidence.
- Do not perform git commit, push, tag, Zenodo release, repository publication, manuscript upload, or email without explicit user authorization.

## Manuscript hierarchy

Use this order of authority:

1. Primary analysis outputs and scripts
2. Audit-locked result ledger
3. Latest v2 full draft
4. Section files and figure captions
5. Older drafts and narrative notes

Older prose cannot overrule primary analysis.

## Writing standard

- Write clear, restrained scientific English.
- State the result first, then the evidence, then the interpretation.
- Keep claims proportional to design.
- Distinguish discovery, validation, replication, concordance, and mechanistic support.
- Avoid hype words such as “transformative,” “definitive,” “unprecedented,” “clinically actionable,” and “precision treatment” unless explicitly justified.
- Do not write “data not shown.”
- Do not use “validated” for an analysis that was only repeated on a related public dataset without prespecified lock and independent endpoint measurement.
- Preserve gene symbols and clinical aliases consistently: SLC5A5 (NIS), NKX2-1 (TTF-1), FOXE1 (TTF-2).

## Protected author voice

Do not independently finalize the six sections listed in `VOICE_PROTECTED_SLOTS.md`. First obtain author language through `/voice-interview`. You may polish supplied author text without changing intent or adding unsupported motivation.

## Figures

- One coherent, sequential figure architecture must be selected. Do not leave orphan numbering such as Fig. 7 and Fig. 8 without a resolved Fig. 6.
- Reconcile the five-figure reviewer-first plan with the v2 multi-figure draft through an explicit editorial decision log.
- Main-text and supplementary figure labels, captions, panel letters, sample sizes, statistics, and source-data filenames must agree exactly.
- Default to Supplementary Figures rather than “Extended Data” unless the current target-journal instructions explicitly permit Extended Data.
- Do not hide negative or non-confirmatory results; position them appropriately and describe them accurately.

## Completion criterion

A task is not complete because prose reads well. It is complete only when:

- claims trace to evidence;
- statistical directions are verified;
- all manuscript sections agree;
- figures and source data are complete;
- limitations match the actual study status;
- data/code availability is actionable;
- voice-protected sections are author-approved;
- independent red-team reviewers find no unresolved critical issue;
- remaining administrative blockers are explicitly listed.
