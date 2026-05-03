# Paper 1 — M3/M4 bib resolution + straggler cleanup (2026-05-04)

**Author:** Seungho Cook
**Date:** 2026-05-04 (marathon scaffolding/infra)
**Scope:** Close M3 (Lee2024 author placeholder) + M4 (Lim2025 title verification) + handle 2 stragglers from prior commit cycle.
**Voice-risk:** ZERO (single-token initial corrections in methods/STAR/data/ack only; no Hook/Aim/Discussion/Limitations/Cover/Q9 touched).

---

## 1. M3 outcome — Lee2024GSE213647

### 1.1 Verification source

PubMed PMID **38331894** (verified via WebFetch from `pubmed.ncbi.nlm.nih.gov/38331894/`) + GEO Series page for GSE213647 (`ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213647`).

| Field | Value (verbatim from source) |
|---|---|
| First author | **S.E. Lee** (NOT "Y") |
| Full author count | 36 authors |
| Title | "Unraveling the role of the mitochondrial one-carbon pathway in undifferentiated thyroid cancer by multi-omics analyses" |
| Journal | Nature Communications |
| Year | 2024 |
| Volume | 15 |
| Issue | 1 |
| Article number | 1163 |
| DOI | 10.1038/s41467-024-45366-0 |
| GEO accession | GSE213647 (linked dataset, Macrogen Seoul, n=632) |

### 1.2 Bib edit applied

Replaced placeholder entry with full author list (36 names) + actual published title + complete Nature Communications citation metadata (journal/year/volume/number/pages/doi/pmid). Note clarifies dataset usage scope.

### 1.3 Manuscript ripple — 5 single-token "Lee Y" → "Lee SE" corrections

The placeholder "Y" leaked into 4 prose files as the cited-author initial. The correct initial is "SE" per PubMed. All 5 occurrences are in non-voice-protected sections (data/methods/STAR/ack):

| File | Line | Change |
|---|---|---|
| `2026_05_03_STAR_Methods_KeyResources.md` | 18 | "Lee Y et al. 2024" → "Lee SE et al. 2024" |
| `2026_05_03_STAR_Methods_KeyResources.md` | 24 | "Lee Y et al. 2024" → "Lee SE et al. 2024" |
| `2026_05_03_paper_supp_tables_draft.md` | 14 | "Lee Y et al. 2024 Macrogen" → "Lee SE et al. 2024 Macrogen" |
| `2026_05_03_manuscript_v8_OUTLINE.md` | 234 | "Lee Y / Macrogen" → "Lee SE / Macrogen" |
| `2026_05_03_acknowledgments_highlights_etoc.md` | 16 | "Lee Y / Macrogen Seoul" → "Lee SE / Macrogen Seoul" |

Voice-protected sections (Hook/Aim/Discussion §3.1/Limitations/Cover/Q9) NOT touched.

**M3 status: ✓ CLOSED.**

---

## 2. M4 outcome — Lim2025GSE286332

### 2.1 Verification source

PubMed PMID **41113708** (verified via WebFetch from `pubmed.ncbi.nlm.nih.gov/41113708/`).

| Field | Bib (before) | PubMed (verified) |
|---|---|---|
| Author count | 2 (Lim, Kim) | 8 (Lim, Jeong, Lee, Choi, Fang, Wang, Kim H, Kim S-M) |
| Title | "Dysregulation of Vitamin D and Its Signaling in Hashimoto's Thyroiditis in Korean Population" | "**Dysregulated vitamin D signaling in Hashimoto's thyroiditis: an integrated transcriptomic study in a Korean cohort**" ⚠ different |
| Journal | "GEO Database" | Frontiers in Endocrinology (Lausanne) |
| Volume / pages | — | 16 / 1666115 |
| DOI | — | 10.3389/fendo.2025.1666115 |
| PMID | 41113708 | 41113708 ✓ matches |

### 2.2 Bib edit applied

Replaced 2-author shorthand with full 8-author list. Replaced paraphrased title with verbatim PubMed title. Added journal / volume / pages / DOI. Note expanded to clarify (a) the published vitamin-D / Hashimoto framing, and (b) our use of the same cohort as a PTC-vs-PTC+HT sub-analysis.

### 2.3 Vitamin-D framing concern (audit T2 §6) — resolved

The audit flagged that the Vitamin-D framing might not match our PTC-vs-PTC+HT usage. After verification: the published paper IS about vitamin-D / Hashimoto integrated transcriptomics, but the cohort design (9 PTC + 9 PTC+HT, matched) is exactly what we use. Our framing is a legitimate sub-analysis of the same dataset. The bib `note` field now explicitly states this distinction.

No manuscript prose change required for M4 — manuscript already uses "Lim DW et al. 2025" form which is correct (et al. covers all 8 authors).

**M4 status: ✓ CLOSED.**

---

## 3. references.bib diff summary

```
project/reports/2026_05_03_references.bib | 24 ++++++++++++++++-------
1 file changed, 17 insertions(+), 7 deletions(-)
```

Two entries replaced (M3 + M4), no other entries touched.

Other entries left at audit-deferred state:
- `Tuttle2019JCEM`, `Yi2016KTA`, `Chen2024EndocrConnect`, `Kim2014KoreanGraves` — 4 entries in bib but not yet cited; per audit recommendation (c), defer per-entry decision to voice-protected Discussion drafting.
- `Pan2025NatComm` — pages format unusual (L6); cosmetic batch W6.
- `Chu2018JMG` — ★ Unicode in note (L5); cosmetic batch W6.
- `Cancer2014TCGA` double-braces institutional author — already correct.

---

## 4. Straggler handling

### 4.1 `project/.runpod_pods.json` → gitignored

Content: 77-byte JSON dict mapping `{D, C, B}` to RunPod pod ID strings (e.g., `"D": "fbg10t9c5ftlcc"`). Runtime/ephemeral state — pod IDs change per session. Not sensitive (pod IDs alone don't grant access; auth is via separate API key kept out-of-tree). Not committed.

`.gitignore` entry added:
```
# RunPod live pod ID dictionary (runtime state; ephemeral pod IDs)
project/.runpod_pods.json
```

### 4.2 `project/notebooks_or_scripts/runpod_pod_B_WSI_pathology.sh` → committed

Sibling of `runpod_pod_C_AlphaFold.sh` and `runpod_pod_D_K2_STAR.sh` (both committed in `25d693c`). 7,616-byte bash setup script for a TCGA-THCA WSI pathology RunPod pod (UNI/CTransPath embedding + Hashimoto-overlap classifier). No secrets / tokens / API keys (verified by `grep -E 'token|api_key|secret|password|bearer|RUNPOD_API|HF_TOKEN|sk-|ghp_|aws_access|AKIA'` returning empty).

This script documents an experimental direction superseded by the Phase A NO-GO verdict (`PHASE_A_NOGO_REPORT_2026_05_04.md`) — committed as historical/sibling artifact, not as an active pipeline. Marathon discipline preserved: committing the script does NOT trigger any download or analysis. Per "TCGA WSI 다운로드 금지" rule, no execution of this script during marathon.

### 4.3 .gitignore diff

```
.gitignore | 3 +++
1 file changed, 3 insertions(+)
```

Single addition (3 lines, comment + entry + blank line). No other ignore rules modified.

---

## 5. Other open items found this session (NOT in scope of this commit)

While checking working tree, additional untracked items appeared post-`25d693c`. **These are NOT addressed in this commit** — flagged for the next decision turn:

| Path | Size | Likely origin | Suggested handling |
|---|---|---|---|
| `project_external_st/src/12_gpu/` | 52 KB (7 .py + 1 .sh) | New external-ST GPU pipeline scripts (g1_embed, g1_ridge_loso, g2_slide_regress, g2_tcga_he_pipeline.sh, g2_tile_embed, g3_celldart, extract_external_tiles) | Inventory only. Per `paper2-pause` decision, do not actively run. Future commit if user wants to archive code only. |
| `project_external_st/results/extra/g1_external_tile_metadata.tsv.gz` | 116 KB | Output from `g1_*` scripts above | Inventory only. Can commit as small artifact or leave untracked pending pause decision. |
| `project/results/03_pathology_poc/tile_metadata.tsv.gz` | 180 KB → 94 KB (modified) | File previously committed in `25d693c` was overwritten by a still-active local pathology process | ⚠ Investigate: a background process may be running on this VM despite paper2-pause. User should check `ps aux` / RunPod pod status. Do NOT commit the new content during marathon. |
| `project/results/03_pathology_poc/tile_metadata_resid.tsv.gz` | 578 KB → 802 KB (modified) | Same as above | ⚠ Same investigation flag. |

These are flagged in this report only — Claude has not committed, deleted, or modified them.

---

## 6. Marathon safety attestation

| Constraint | Status |
|---|---|
| Voice-protected sections (Hook / Aim / Disc §3.1 / Limitations / Cover ¶1 / Q9) | ✓ untouched |
| Paper 3 design bundle (chmod 444) | ✓ untouched |
| Paper 3 Track B start | ✓ NOT started |
| Paper 4 (Korean GD HLA backlog) | ✓ untouched |
| New analysis | ✓ none performed |
| New data download | ✓ none performed |
| H&E-DM1 closure retry | ✓ NOT attempted |
| TCGA WSI download | ✓ NOT attempted |
| WebFetch usage | ✓ minimal — 3 read-only PubMed/GEO landing-page pulls (PMID 38331894, PMID 41113708, GEO GSE213647) for citation metadata only |

---

## 7. Remaining open citation-related issues

| ID | Issue | Status | Action |
|---|---|---|---|
| L4 | references.bib header claims "~22 papers" — actual 29 entries | deferred | W6 cosmetic batch |
| L5 | ★ Unicode in `Chu2018JMG` note may break BibTeX | deferred | W6 cosmetic batch (replace with ASCII "(NOTE)") |
| L6 | `Pan2025NatComm` pages format unusual ("s41467-025-58910-3") | deferred | W6 cosmetic batch (verify Nat Comm convention) |
| L7 | 4 bib entries not yet cited (Tuttle2019, Yi2016, Chen2024, Kim2014) | deferred | Decide per-entry during voice-protected Discussion drafting |
| L9 | Possible Krishnamoorthy 2025 bib gap (already aliased to `Landa2016JCI` per memory `v17_landa2016_cite_save`) | deferred | User voice draft of Discussion §3.1 |

No remaining MEDIUM citation issues. All HIGH issues remain CLOSED.

---

## 8. Recommended next action

Two safe options:

- **`voice-hook`** — User opens Hook ¶1 in `manuscript_v8_OUTLINE.md` and writes voice-protected content directly. First voice-protected slot of marathon. Earliest unlock for Discussion §3.1 (Landa 2016 cite save) work later.
- **`low-batch`** — Pull L1–L9 cosmetic pass forward from W6 to now (~1–2 hours, single-file edits, all marathon-safe scaffolding). Includes L4/L5/L6 bib polish + L1/L2/L3 numerical rounding unification + L7 bib reserve cleanup decision + L8 ★-removal in F2D caption.

Either is valid. `voice-hook` unblocks user-keyboard work; `low-batch` continues claude-keyboard scaffolding.

---

Bib M3/M4 closed. Stragglers handled. Marathon discipline preserved.
