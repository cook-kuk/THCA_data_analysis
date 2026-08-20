# Task 16 — Reconciliation: `20260820_data/W_LEDGER_COMPLETENESS_SCREEN.md` vs the 2026-08-20 completion-session SOURCE closures

**Date** 2026-08-20. Requested by the author after two independent sessions, run the same day,
produced conflicting closure decisions for the same five sources named in
`V3_FINAL_DECISION_MEMO.md` §G item 1 / §L (Mu 2024/HRA004166, GSE299988, GSE112202, GSE184362,
GSE221088). This task reads both in full and reports where they agree, where they conflict, and
which side's claim is corroborated by an artifact actually on disk.

**Inputs read in full**: `20260820_data/W_LEDGER_COMPLETENESS_SCREEN.md` and
`20260820_data/CONSOLIDATED_SCREENING_LEDGER_V4_2026_08_21.tsv` (the untracked drop under review);
this directory's `SOURCE_mu2024_hra004166.md`, `SOURCE_gse299988.md`, `SOURCE_gse112202.md`,
`SOURCE_gse184362.md`, `SOURCE_gse221088.md`, `TASK13_SOURCE_LEVEL_DECISIONS_rows145.md`; file
listings and timestamps for both trees.

**File-safety statement**: no existing file was edited or deleted. `20260820_data/` is untracked
and unmodified. Nothing under `audit/rai_integration_20260820/` was touched. This is a new file,
continuing this directory's own TASK-numbering (TASK13–15 already exist). **Provenance note on
this task itself**: this is a desk comparison of two prior sessions' written outputs — no new
primary-source fetch, no new molecular test, no new file download was performed here. Where a
claim below rests on an artifact "on disk," that means a file already present in the repository
before this task started (verified by `ls`/timestamp only, not re-fetched).

**Repository-layout note**: this file was authored in a git worktree (bg-session isolation), while
`audit/rai_integration_20260820_completion/` in the author's working copy is untracked. It is not
yet copied into the working copy — see the session's closing report for the exact path to copy it
from.

---

## 0. Headline

**The two sessions do not merely disagree on classification — the completion session holds
archived primary-source files and completed test results that the W memo says are impossible to
obtain.** Where they conflict, the completion session's claims are the ones corroborated by an
artifact on disk (a `.soft.gz` file, a byte-identical reproduction of a 2026-05-21 result, a raw
10x matrix already local since 2026-08-19). The W memo's `[UNVERIFIED-TRANSCRIPTION]` tier should
be treated as superseded, not as the more cautious of two equally-uncertain estimates, for four of
the five sources.

**The most consequential single finding**: W's GO/NO-GO section frames "download the GSE299988
xlsx and run the differentiation-axis test" as a bonus next step, not yet done, blocked by a 403.
The completion session shows this test was **already run in the original 2026-05-21 pipeline**,
using the series matrix (not the xlsx W was chasing), and reproduces it byte-for-byte. A sixth
test W never attempted — GSE184362 — was **also already run**, on 2026-08-19. **Both already-run
tests are directionally null or reversed, not supportive.** R4 does not carry "four completed
tests and one pending" — depending on how the two sessions' Layer/status calls are reconciled, it
may already carry five or six, and the newest ones cut against the paper's central hypothesis, not
for it.

---

## 1. The archival/access contradiction — check this first, it governs everything else

W's provenance section states: fetched pages "return a model-summarised rendering rather than the
raw record, and no fetched page was archived; file download is separately blocked by the sandbox
proxy (`curl` to `ftp.ncbi.nlm.nih.gov` returns HTTP 403 at the tunnel), so no SOFT, MINiML or
supplementary file could be stored either."

This directory contains, timestamped 2026-08-20 09:40–09:58 (i.e. the same calendar day, and
before the 15:03 timestamps on the `20260820_data/` files):

| Artifact | Source claim | On disk? |
|---|---|---|
| `_raw_GSE221088_series.txt` | Series-level SOFT record, fetched `form=text&view=full` | **Yes**, 2026-08-20 09:47 |
| `_raw_GSE221088_family.soft` + `.soft.gz` | Full Series+Sample family SOFT file, 160 GSM records, "fetched directly from the NCBI GEO FTP mirror `ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221088/soft/`" | **Yes**, 2026-08-20 09:47 |
| PMC11031230 full-text extraction (Mu 2024 Methods/Results/DAS, quoted) | Two live WebFetch passes | Not saved as a raw file, but quoted with quotation marks in `SOURCE_mu2024_hra004166.md` |
| NGDC HRA004166 accession page (title, access-status, DAC contact) | Two live WebFetch passes + 1 web search | Same — quoted, not archived as a separate raw file |
| GSE184362, GSE299988, GSE112202 raw data (10x matrices, series matrices, FPKM tables) | Already local from **earlier** pipeline runs (2026-05-21 / 2026-08-19), reused here, not re-downloaded this session | **Yes** — pre-existing at `/data/rai_atlas/raw/...` and `rai-response-genomics-atlas/data/raw/...` |

The completion session's own GSE221088 file states the fetch path explicitly: the same
`ftp.ncbi.nlm.nih.gov/geo/series/.../soft/` host and path W names as returning 403. One of these
two same-day sessions is describing an environment that does not match the other's. Two
non-exclusive explanations, neither of which resolves in W's favor as written: (a) the completion
session may not have used literal `curl` — a different fetch mechanism could reach the same URL
where a sandboxed `curl` cannot; (b) network egress may differ between sessions/sandboxes. Either
way, **the blanket claim "no SOFT/MINiML/supplementary file could be stored" is false for at least
one of the five sources, with a file on disk to prove it**, and the claim "no fetched page was
archived" undersells what the completion session actually saved (a full raw SOFT family file, not
merely a summarized quote).

**Consequence for the evidence-lock tiering**: W's `[UNVERIFIED-TRANSCRIPTION]` tag on GSE221088's
sample counts, patient counts, and field inventory should be **upgraded** — the completion
session's numbers for this source come from a full programmatic parse of an archived raw file, the
strongest tier this project's own convention recognizes, not a tool-summarized paraphrase. The same
upgrade applies, with the caveats noted per-source below, to GSE299988, GSE184362, and GSE112202,
whose core facts (sample counts, condition labels, panel-gene coverage) all trace to raw files
already on disk before either session started — not to anything fetched fresh this session by
either side.

---

## 2. Mu 2024 / NGDC HRA004166

| | W memo | Completion (`SOURCE_mu2024_hra004166.md`) |
|---|---|---|
| n=220 / n=214 / n=281 | Resolves identically: 220 enrolled, 214 with determinable RAI uptake (6 missing scan), 281 specimens | **Same resolution**, independently reached |
| Layer | **1**, INCLUDED (row 32) — "an extension of the written rule, not an application of it," explicitly flagged as the author's call to ratify or reverse | **2** — literature-only, no data product obtained or re-derived; explicitly "not promoted to Layer 1 merely because the accession exists" |
| Counts toward N/M | Yes — this is one of the two sources driving W's M: 10→12 | **No** — Layer 2 by this session's mechanical application of the same rule text W quotes |
| Assay | Not separately assessed in W beyond "no expression measurement… fails Block I-1" (§1c) | Same conclusion, plus an explicit 8-gene coverage table: **0/8 panel genes** have any expression readout on ThyroLead; 2/8 (TSHR, PAX8) exist only as DNA mutation-call targets. Classified **ASSAY-LIMITED**, independent of DAC status |
| Driver-class contingency (BRAF 56/87, RAS 1/22, Fusion 12/58, Others 9/43; OR 8.29 BRAF vs non-BRAF, P=7.7e-12) | **Central finding of W's §1d** — described as re-derived by script from prose percentages, with an explicit caveat that collapsing to "BRAF/RAS-negative vs rest" is internally incoherent | **Not present at all.** The completion session's two PMC fetches were prompted for enrollment/analysis denominators, class definitions, the DAS, and assay platform — not for the driver-mutation-by-I-RAIR contingency table. This number is **neither corroborated nor contradicted** by the completion session; it remains sourced only to W's own fetch. |
| Class breakdown 80/48/19/10 (existing manuscript drafts) | Not addressed | **Flags as UNVERIFIED** — the primary-text fetch could only locate I-RAIR (n=80) vs the I-RAIA umbrella (n=134) in aggregate; the 48/19/10 split was only locatable as *within-BRAF-mutant-subgroup percentages* (n=31), not whole-cohort counts. 48.4% is noted as suspiciously close in surface form to "n=48." |
| DAC / access | Not the focus of W's Mu 2024 section | Controlled access, DAC `HDAC002371`, contact confirmed to match this project's existing DAC-request tracking file |

**Reconciliation.** The n=214/220/281 resolution is corroborated twice, independently — treat that
part as settled. Everything past it is where the two sessions disagree, and the disagreement is
about the Layer-1 entry rule itself (which is exactly §3 of W's own memo — "the rule does not
describe the ledger it governs"). Applying W's own §3 diagnosis to this row: under the rule *as
written*, Mu 2024 fails it (no data product obtained or re-derived) — which is the completion
session's position. W's Layer-1 placement for Mu 2024 is the "screened-to-a-decision" reading W
itself says is an extension, not an application, and W explicitly asks the author to ratify or
reverse it. **This reconciliation recommends treating Mu 2024 as Layer 2**, consistent with the
completion session and with a literal reading of the rule text both sessions quote — which W
itself computes drops N to 35 and M to 11 (§9 of W's own sensitivity table).

The driver-class contingency table (OR 8.29) is the single highest-stakes number in either
document (it gates R3, pillar 3, and the abstract per W §8 item 2) and is currently **sourced to
one session's unarchived fetch only**. Neither session's artifact confirms or refutes it against a
stored primary source. This is the one open item from this reconciliation that most needs a fresh,
archived re-fetch of Mu et al. 2024 (PMC11031230) — specifically Table 1/2 or Fig. 3c — before OR
8.29 or the abstract sentence built on it are used again.

---

## 3. GSE299988

| | W memo | Completion (`SOURCE_gse299988.md`) |
|---|---|---|
| Layer / status | Layer 1, INCLUDED, "candidate fifth Category-A direct test" | Layer 1, INCLUDED, **completed** Category A test |
| Test run? | **"The test was not run"** — blocked on `GSE299988_Processed_data_files.xlsx`, HTTP 403 | **Yes — twice.** Re-ran the parent project's own 2026-05-21 pipeline (`04_score_gene_panel.py` / `05_validate_rai_labels.py`) against the **already-locally-cached raw series matrix** (`/data/rai_atlas/raw/GSE299988/...`, pulled 2026-05-21, not re-downloaded this session), byte-for-byte identical output to what was already on disk. Independently re-derived a second time with a from-scratch script, same result to 3–4 s.f. |
| Result | Not available (test not run) | Cohen's d = **+1.128** (refractory − avid), reversed from predicted polarity; MWU P = **0.1508**; exact permutation P = **0.1270**; analytic 95% CI [−0.207, 2.463] spans zero |
| Confound | RAI status = LN status by design; n=5 vs 5 | **Same confound**, independently derived from the raw `characteristics_ch1` text: N-status, M-status, stage, and ATA risk each separate the two arms perfectly (5/5 vs 5/5) |
| Data actually used | W was chasing the processed **xlsx** supplementary file | Completion used the **series matrix** (`GSE299988_series_matrix.txt.gz`), which is a different, already-obtained file — the processed xlsx was never needed |

**Reconciliation.** This is not a disagreement about the science — the confound characterization
is identical on both sides. It is a **process gap**: W did not check whether this project already
held and had already used GSE299988's raw data before concluding the test was blocked. It was not
blocked; it was already done, three months earlier, and re-verified twice on 2026-08-20 itself.
**W's GO/NO-GO §10 "download and run" action item is stale and should be replaced** with: GSE299988
is a completed, directionally-null/reversed, confounded fifth test, already in hand.

---

## 4. GSE184362 — the sharpest classification conflict

| | W memo (row 35) | Completion (`SOURCE_gse184362.md`) |
|---|---|---|
| Status | **EXCLUDED**, R1b mode 2 ("molecular data exists, RAI phenotype absent") | **INCLUDED**, R4 Category A, completed direct test |
| Basis for W's exclusion | "The words RAI / radioiodine / refractory / avid appear in **zero** sample titles and zero characteristics fields... independently event-limited at n = 2" | A four-category per-sample "tissue/treatment condition" field — `primary_untreated` / `RAI-treated_distant_met` / `RAI-refractory_LN` / `RAI-refractory_distant_met` — which the completion session states is "recovered from the GEO record's sample titles and characteristics" |
| Test status | Not run (excluded pre-test) | **Already run on 2026-08-19**, independently re-verified 2026-08-20 with a from-scratch script against the raw 10x barcode/feature/matrix triplets (already local, not re-downloaded) |
| Result | N/A | 6 samples / **5 patients** (not "n=2" as W's framing implies): direction opposite the bulk-cohort hypothesis (RAI-refractory shows *higher*, not lower, panel expression); no formal test reportable at this size; one reproducibility slip found and corrected (a pooled mean, +0.00 reported vs +0.15 re-derived — does not change the qualitative pattern) |

**This is a direct, substantive disagreement, not a documentation gap, and it needs the author's
own read of the raw GEO record to settle — not a repeat of either bot's characterization.**

What can be said now: the completion session's "condition" labels for GSE184362 are **not** all
literal GEO field values in the same sense as, say, GSE299988's `rai response: Avid`. Reading
`SOURCE_gse184362.md` closely: sample identity (patient 11, right-LN vs subcutaneous, "three
completed courses of iodine ablation") is described there as recovered from "GEO sample titles and
characteristics," but the underlying per-patient treatment-course history (1× vs 3× RAI ablation)
most plausibly comes from the **published paper's** clinical description of these specific patients
(Pu et al. 2021), matched onto GEO sample IDs by the project's own prior 2026-08-19 pipeline — not
from a `characteristics_ch1: RAI-refractory` field GEO itself carries. W's claim ("zero hits for
RAI/refractory/avid in sample titles and characteristics fields") is plausibly **true of the raw
GEO deposit taken alone** — GSE184362's actual sample titles are anatomic (`PTC11_RightLN`,
`PTC11_SC`), not RAI-status labels. If so, both sessions are right about different things: W is
right that GEO's own fields don't carry the word "RAI"; completion is right that a genuine
treatment-history phenotype for these specific patients exists and is documented **in the paper**,
already matched to sample IDs by this project's own earlier work — which is a Block II-7-adjacent
question (is the phenotype recoverable from the GEO record alone, or does it require an external,
already-performed match to the paper) that neither file states in those terms.

**This reconciliation does not adjudicate INCLUDED vs EXCLUDED for GSE184362.** It flags that: (1)
a completed test with a specific, reversed-direction result already exists in this repository,
dated 2026-08-19, independent of W's session; (2) W's "n=2, event-limited" framing undercounts the
actual analysis unit used (5 patients / 6 samples, per the existing pipeline) if the phenotype is
accepted at all; (3) the question that actually decides EXCLUDED-vs-INCLUDED is whether a
phenotype matched from the *paper* onto GEO sample IDs by this project's own earlier work satisfies
the Layer-1/R1b entry rule the same way a GEO-native field does — which is a specific instance of
W's own §3 finding (the entry rule doesn't cleanly describe the ledger) and should be decided by
the same author ruling that resolves §3, not separately per row.

---

## 5. GSE112202

| | W memo (row 34) | Completion (`SOURCE_gse112202.md`) |
|---|---|---|
| Status | EXCLUDED, R1b mode 2 — "RAI wording is in the title, not the data... a drug-exposure proxy" | EXCLUDED, for **two independent, compounding** reasons |
| Sample-count discrepancy | Flags 22 (file O) vs 25 (live record) as unresolved | **Same discrepancy**, same two numbers, also left open (UNVERIFIED) — independently converged on |
| Reason 1 | (implicit) phenotype absent from deposit | **PHENOTYPE-BLOCKED, explicitly**: the underlying paper (Tesselaar 2021) *does* report a genuine per-patient post-RAI outcome (Table 2: remission/persistent/recurrent, 3–15y follow-up) — but Table 2 uses independent sequential patient numbering that cannot be cross-walked to the GEO tissue IDs (`T08-33030`-style). No linkage exists, confirmed by a second targeted live fetch. |
| Reason 2 | Not raised | **ASSAY-LIMITED, independently of phenotype**: the series matrix's own expression block is empty, and all seven deposited FPKM files are 2-column arm-level aggregates (untreated vs digoxin-treated), never per-sample — so the pre-specified per-patient z-score method cannot be computed from this deposit under *any* phenotype choice |
| Data already held? | Not stated | **Yes** — this project has held and used these seven FPKM files since 2026-05 (`gse112202_redifferentiation_quick.py`, two existing figures); re-run today, reproduces the prior direction-of-effect (6/8 genes up, median log2FC +0.302) exactly |

**Reconciliation.** No conflict on the bottom line (both EXCLUDE it), but the completion session's
account is substantially more complete and should supersede W's shorter entry: it identifies a
real clinical RAI-outcome phenotype in the source literature (something W's "title, not data"
framing undersells), then shows independently why it still cannot be used (unlinkable, and the
molecular deposit is structurally arm-level-only regardless). It also surfaces a wording-register
caution on two of this project's own existing figure captions ("digoxin restores 8-gene panel
expression") that overstate what a two-value-per-arm log2FC can support — worth the author's
attention regardless of this ledger question.

---

## 6. GSE221088

| | W memo (row 36) | Completion (`SOURCE_gse221088.md`) |
|---|---|---|
| Status | EXCLUDED, construct mismatch, R1b mode 2 | EXCLUDED, PHENOTYPE-BLOCKED |
| RAI variable | None | None — confirmed by an **exhaustive, archived** search: 18 keywords × all fields × all 160 GSM records in the raw family SOFT file, zero hits anywhere |
| Sample arithmetic (41×4=164 vs 160 deposited) | Left **explicitly open**, flagged `[UNVERIFIED]` | **Resolved**: FTA-14 and FTA-19 are missing **both** cell-free library types (long + small RNA), i.e. 2 patients × 2 missing library types = 4; 164 − 4 = 160. Independently re-derived from per-sample library-type tallies (41 EV-long + 41 EV-small + 39 cf-long + 39 cf-small = 160), not merely asserted. |
| Why structurally excluded | Diagnostic endpoint, no RAI variable | **Same**, plus a structural argument W does not make: blood was drawn **preoperatively**, before the surgery that establishes the FTC/FTA diagnosis at all — so no patient could yet have been RAI-eligible at collection time. This is a design-level reason no RAI phenotype could exist here, independent of what fields happen to be populated. |

**Reconciliation.** No conflict on the exclusion decision. The completion session closes the one
item W left open (the 164-vs-160 arithmetic) and adds a structural (not just empirical) reason for
the exclusion. Treat the arithmetic as resolved; adopt the "preoperative, therefore structurally
RAI-ineligible" framing as the stronger version of the exclusion rationale.

---

## 7. Net effect on N/M if this reconciliation's classifications are adopted

Substituting the completion session's Layer/status calls for W's, row by row (Mu 2024 → Layer 2,
not counted; GSE299988 → INCLUDED, completed; GSE112202 → EXCLUDED, counted in N only; GSE184362 →
provisionally INCLUDED pending the author's §4 ruling; GSE221088 → EXCLUDED, counted in N only):

| | V3 (08-20) | W's V4 | This reconciliation, if adopted |
|---|---|---|---|
| Layer-1 N | 31 | 36 | **35** (Mu 2024 moves to Layer 2) |
| M (usable RAI variable) | 10 | 12 | **12** — same count, different composition: **Mu 2024 out, GSE184362 in** |
| R4 Category-A completed tests | 4 | 4 (+1 pending: GSE299988) | **6** — GSE299988 and GSE184362 are both already-completed, not pending, and both are directionally null/reversed |

The number 12 surviving under both accountings could read as reassuring, but the composition
change is the material fact: **W's two additions to M were one driver-class-heterogeneity exhibit
(Mu 2024, provisionally supportive-looking, OR 8.29) and one untested candidate (GSE299988). This
reconciliation's two additions are two *completed* tests (GSE299988, GSE184362), both directionally
null-or-reversed relative to the paper's central hypothesis.** That is a different, and less
favorable, story for R3/R4 than either session states on its own, and it should be decided before
R3, R4, or the abstract are drafted — not after, per the same reasoning W's own §8 item 2 already
gives for Mu 2024 alone.

---

## 8. What the author needs to decide, in order

1. **Mu 2024's Layer**: this reconciliation recommends Layer 2 (literature-only, DAC-gated,
   assay-limited), following a literal reading of the rule both sessions quote. If the author
   instead ratifies W's Layer-1 extension, that widening should be applied consistently to
   GSE184362 too (§4), since both rest on the same underlying question — does a phenotype
   documented in the source paper and matched to GEO sample IDs by this project's own work count
   toward Layer 1's "re-derived… data product" test.
2. **The Mu 2024 driver-class contingency table (OR 8.29)**: sourced to one unarchived fetch only,
   corroborated by neither session's saved artifact. Needs a fresh, archived re-fetch of
   PMC11031230 (ideally the actual PDF tables/figures, not a WebFetch summarization pass) before
   any manuscript sentence uses it.
3. **GSE184362 INCLUDED vs EXCLUDED**: the sharpest live conflict (§4). Whichever way it resolves,
   R4/pillar-3 language must account for the fact that a completed, reversed-direction result for
   this source already exists in the repository, dated before either of today's two sessions.
4. **W's own §3 finding** (the Layer-1 entry rule doesn't describe the ledger it governs) is
   unaffected by anything here and remains, in this task's view, the correct diagnosis of the root
   problem — items 1 and 3 above are two more instances of exactly the gap W's §3 already names,
   not new problems.
