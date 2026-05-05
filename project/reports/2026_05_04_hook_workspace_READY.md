# Hook ¶1 — Fact-Only Workspace (NOT a draft)
**Date:** 2026-05-04
**Owner:** Seungho Cook
**Status:** READY for user-keyboard Hook ¶1.
**Voice protection:** Hook is voice-protected. **Claude must not draft Hook prose.** This file holds *only* facts, structural labels, and forbidden lists — not sentences.

---

## 0. Hard rule

> Hook ¶1 is **user keyboard only**.
> Claude's role here is *zero*: this file gives the user the materials. After the user writes Hook ¶1, Claude may audit it (factual accuracy, overclaim, scope contamination, citation risk). Claude must not propose, suggest, paraphrase, or "lightly polish" Hook sentences.

---

## 1. Recommended 4-sentence skeleton (labels only — no example sentences)

| # | Label | Function |
|---|---|---|
| S1 | Clinical risk gap | Establish that PTC outcomes are heterogeneous and current risk stratification leaves a residual gap. |
| S2 | Driver paradigm utility | Acknowledge what BRAF/RAS/fusion-driver paradigm did achieve. |
| S3 | Driver-only limitation | State what the driver-only paradigm does *not* explain. |
| S4 | Transcriptional differentiation axis gap | Frame the unmet need that motivates the paper without naming the answer. |

Optional: S1+S2 may be combined if it tightens the paragraph; S3 and S4 should remain separate to keep the gap statement crisp.

---

## 2. Safe numbers (use at most TWO; do not invent precision)

- **DTC recurrence:** 15–35% over long-term follow-up. (Range; do not state a single number.)
- **ATA risk-of-recurrence range:** ~3–5% (low) to ~50–75% (high).
- **Optional: one driver-paradigm limitation number** — use only if you have a specific, citable figure ready *at the time of writing*. Do not insert a placeholder.

> Hard cap: **maximum 2 numbers** in Hook ¶1. Citation discipline applies — every number must be tied to a real reference at the moment of writing.

---

## 3. Do NOT include in Hook ¶1

- DM1 / DM2 / molecular dark matter terminology
- 8-gene composition (no gene names, no "RAI lineage", no "THYROID_NONOVERLAP")
- Hashimoto / HT / autoimmune / HLA
- TROP2 / TACSTD2
- GSE76039 / Landa 2016 / ATC / PDTC
- H&E / WSI / pathology-AI / morphology
- GigaTIME / spatial / CosMx / Visium
- Paper 2, Paper 3, Paper 4 terms (Hashimoto-overlap PTC, ICI vulnerability, Korean GD HLA)
- "Korean cohort", "K2", "Bundang", "AFND"
- Korean-specific framing (Pan-Asian, etc.)
- Survival HRs from this study
- Any reference to *our* findings — Hook ¶1 is **the gap**, not the answer

---

## 4. Tone & overclaim guard (for the user to self-check)

Avoid in Hook ¶1:
- "We show / we demonstrate / our work" — that belongs in Aim.
- "Proves", "establishes", "all risks resolved", "Cancer Cell-ready", "Nat Cancer reach".
- Monotonic-progression language (no "PTC → PDTC → ATC").
- Mechanistic causal claims (no "fusion-driven, epigenetically silenced", no "promoter hypermethylation drives X").
- Pathology-AI claims (no "H&E-inferable", "WSI-validated", "morphology-derived DM1").
- Spatial co-localisation claims (no "DM1-high spots are TROP2-high", no "spatially colocalizes").

Tense: present tense for the field state, simple past for what prior work did.
Voice: third-person field framing in S1–S3; S4 may shift to "remains underexplored / unresolved" — still no first-person.

---

## 5. After-draft Claude audit checklist (Claude runs this AFTER the user has written Hook ¶1)

Claude will check, **without rewriting**:

1. **Factual accuracy** — every number, range, and claim is verifiable from a citation already in the bibliography or specifiable on the spot.
2. **Overclaim** — none of the forbidden phrases in §4 appears.
3. **Scope contamination** — none of the §3 forbidden topics appears.
4. **Citation risk** — any non-trivial claim has a reference attached or flagged for the user to attach.
5. **Internal consistency** — Hook ¶1 does not contradict the Aim or the manuscript's framing of Paper 1 as a driver-orthogonal transcriptional differentiation axis.

Output of audit: a bulleted list of *findings only*. **No suggested rewrites.** If a fix is needed, Claude flags it; the user edits.

---

## 6. Files to keep open while writing

- `project/manuscript_v8/04_intro_1_1_hook.md` — write here.
- `project/manuscript_v8/03_introduction.md` — for paragraph-flow context only.
- `project/reports/2026_05_04_paper1_strategy_plus_hook_checklist.md` — the longer Hook checklist (if it has one).
- This file (`2026_05_04_hook_workspace_READY.md`) — the fact-only workspace.

---

**End of workspace. Claude did not draft a Hook sentence. Ready for user-written Hook.**
