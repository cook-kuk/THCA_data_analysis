# Paper 1 — LOW batch L1–L6 closure (2026-05-04)

**Author:** Seungho Cook
**Date:** 2026-05-04 (marathon scaffolding/infra; L1–L9 cosmetic batch — claude-keyboard portion only)
**Mode:** Single-token / single-line edits across 3 active manuscript files. Diff-friendly.
**Voice-risk:** ZERO. L7/L8/L9 are voice-protected — deferred to user keyboard.

---

## 1. Closures applied (L1–L6)

### L1 — TIERA67 / pan-genome ARI rounding ✓

Unify all instances to 3-decimal precise form (already used in Results prose and QA):

| File | Before | After |
|---|---|---|
| `manuscript_v8_OUTLINE.md:24` (Abstract) | `TIERA67 ARI=0.90 ≈ ... top-5000 MAD ARI=0.92` | `TIERA67 ARI=0.903 ≈ ... top-5000 MAD ARI=0.918` |
| `manuscript_v8_OUTLINE.md:146` | `TIERA67 ARI=0.90, pan-genome top-5000 MAD ARI=0.92` | `TIERA67 ARI=0.903, pan-genome top-5000 MAD ARI=0.918` |

R-PRO (`results_R1_R5_prose.md` line 41) and QA (`reviewer_QA_consolidated.md` line 47, 50) already used precise form — no change needed.

### L2 — Fisher p (Hashimoto-DM2) precision ✓

Unify shorthand `6e-10` → precise `6.4×10⁻¹⁰` (consistent with R-PRO line 61 source-of-truth):

| File | Before | After |
|---|---|---|
| `manuscript_v8_OUTLINE.md:24` (Abstract) | `OR up to 5×, p=6e-10` | `OR up to 5×, p=6.4×10⁻¹⁰` |
| `manuscript_v8_OUTLINE.md:160` (R5b) | `OR=0.20, p=6e-10` | `OR=0.20, p=6.4×10⁻¹⁰` |
| `figure_captions_all.md:139` (F5B) | `Fisher p=6e-10 to 8e-9` | `Fisher p=6.4×10⁻¹⁰ to 8e-9` |

`figure_captions_all.md:367` already used `6×10⁻¹⁰` (caption summary, narrative shorthand) — left as-is. Outreach email + internal scaffolding files (D8D / GRAND_CONSOLIDATION / code_repo_README) NOT touched (out of v8 manuscript scope).

### L3 — BRAF transcript Cohen d precision ✓

Unify rounded `−0.04` → precise `−0.044` (R-PRO line 31 + GRAND_CONSOLIDATION precise form):

| File | Before | After |
|---|---|---|
| `manuscript_v8_OUTLINE.md:24` (Abstract) | `Cohen d=−0.04 vs WT (n=273 vs 182, p=0.57)` | `Cohen d=−0.044 vs WT (n=273 vs 182, p=0.567)` |
| `manuscript_v8_OUTLINE.md:136` (R3) | `Cohen d=−0.04 (p=0.57, n=273 vs 182)` | `Cohen d=−0.044 (p=0.567, n=273 vs 182)` |

p-value also unified to 3-decimal (`0.57` → `0.567`) for consistency within the same parenthetical.

### L4 — references.bib citation count header ✓

| File | Before | After |
|---|---|---|
| `references.bib:380` | `% Citation count: ~22 papers + tools.` | `% Citation count: 29 entries (papers + tools).` |

Verified count: `grep -c '^@article' references.bib` = 29. ✓

### L5 — `Chu2018JMG` note: ★ Unicode → ASCII ✓

| File | Before | After |
|---|---|---|
| `references.bib:73` | `note = {... ★ Citation correction: ...}` | `note = {... (NOTE) Citation correction: ...}` |

Replaces the ★ Unicode with ASCII `(NOTE)` to avoid BibTeX-pipeline encoding risk.

### L6 — `Pan2025NatComm` page format ✓

Verified via PubMed PMID **40234451** (WebFetch): Nature Communications 16, **Article number 3601**, DOI 10.1038/s41467-025-58910-3.

| File | Before | After |
|---|---|---|
| `references.bib:181` | `pages = {s41467-025-58910-3}` (DOI suffix, malformed) | `number = {1}, pages = {3601}` |

The `s41467-025-58910-3` was the DOI suffix mistakenly placed in the `pages` field. Corrected to the proper Nature Communications article number (3601) and added `number = {1}` per Cell Press / BibTeX convention.

### Bonus — Outline ★ decorative marker

While editing the Abstract for L1+L2+L3, I noticed and dropped a decorative `★` prefix on result-bullet (5):

| File | Before | After |
|---|---|---|
| `manuscript_v8_OUTLINE.md:24` | `(5) **★ Autoimmune-PTC mechanism**:` | `(5) **Autoimmune-PTC mechanism**:` |

The ★ was a typographic highlight on a bullet number, NOT semantic content. Removal is consistent with the L5 ASCII rule. Other bullets (1)–(4) had no ★, so consistency is improved. **NOT a voice-protected change.**

---

## 2. Deferred (L7, L8, L9 — voice-protected)

| ID | Issue | Reason for deferral |
|---|---|---|
| L7 | 4 bib entries not yet cited (Tuttle2019 / Yi2016 / Chen2024 / Kim2014) | Per-entry decision belongs to user voice draft of Discussion §3.1 / §3.2 / §3.3 |
| L8 | "★ exceptional" verb in F2D caption — borderline subjective | User voice polish; flagged in caption text, not here |
| L9 | Possible Krishnamoorthy 2025 bib gap (already aliased to `Landa2016JCI`) | Decide at §3.1 voice draft; memory `v17_landa2016_cite_save` carries the 3-layer reverse-causality argument |

These remain on the user-keyboard track. No claude action.

---

## 3. Files modified (this commit)

```
project/reports/2026_05_03_manuscript_v8_OUTLINE.md    | 8 ++--  (L1+L2+L3 + ★ removal)
project/reports/2026_05_03_figure_captions_all.md      | 2 +-   (L2)
project/reports/2026_05_03_references.bib              | 7 ++--  (L4+L5+L6)
project/reports/2026_05_04_low_batch_L1_L6_closure.md  | NEW    (this report)
```

**Voice-protected sections untouched:** Hook ¶1, Aim ¶4, Discussion §3.1, Limitations §3.4, Cover letter ¶1, Reviewer Q9 — all confirmed by edit-target inspection.

---

## 4. Marathon safety attestation

| Constraint | Status |
|---|---|
| Voice-protected sections | ✓ untouched |
| Paper 3 design bundle | ✓ untouched |
| Paper 3 Track B | ✓ NOT started |
| Paper 4 backlog | ✓ untouched |
| New analysis | ✓ none |
| New data download | ✓ none (one read-only WebFetch to PubMed for L6) |
| H&E-DM1 retry | ✓ NOT attempted |
| TCGA WSI download | ✓ NOT attempted |
| RunPod API | ✓ NOT used |
| Background pod backflow files (`PAPER1_DM1_FULL_2026_05_04.md` modified +49 lines, etc.) | ✓ NOT included in this commit (left for user review) |

---

## 5. Remaining open issues post-L1–L6

- **L7 / L8 / L9** — voice-protected, deferred (above)
- **Background pod backflow** — `PAPER1_DM1_FULL_2026_05_04.md` modified, plus untracked Pod C/D output files. **Pod C/D need user-side stop on RunPod** before clean working tree is achievable.
- **`project_external_st/src/12_gpu/` + `post_pod_*.py`** — defer-to-user decision pending Pod outcome.

---

## 6. Recommended next action

- **`voice-hook`** — User opens Hook ¶1 in `manuscript_v8_OUTLINE.md` and writes voice-protected content directly. First voice-protected slot of marathon. Unblocks L7/L8/L9 + Discussion §3.1 (Landa 2016 cite save).
- **`pod-cleanup`** — Stop Pod C/D on RunPod console, then triage backflow files (`PAPER1_DM1_FULL_2026_05_04.md` review + 12_gpu/ commit-or-delete decision + post_pod scripts decision).

Either fits marathon. `pod-cleanup` first probably saves voice-hook from being interrupted by further backflow.

---

L1–L6 closed. Voice-protected work remains user-driven. Marathon discipline preserved.
