# bib-m3m4 + low-batch L1-L6 — Completion Report

**Author:** Seungho Cook
**Date:** 2026-05-04 (post-Pod-C/D dispatch, marathon scaffolding session)
**Trigger:** "다 해줘" — complete remaining marathon-safe scaffolding while Pod D monitor runs in background.
**Mode:** Audit-driven cleanup. **No analysis. No voice-protected prose. No manuscript section drafting.**
**Voice-protection respected:** Hook ¶1, Aim ¶4 (intro §1.4), Discussion §3.1, §3.4 Limitations, Cover ¶1, Reviewer Q9 — all untouched.

---

## 1. Executive summary

| Item | Source memo | Final state | This session |
|---|---|---|---|
| **M3** Lee2024 GSE213647 author placeholder | recheck §2 ⏳ deferred | ✅ already correct in active bib (37 authors, matches PMID 38331894) | verified via PubMed; **no edit needed** |
| **M4** Lim2025 GSE286332 title verify | recheck §2 ⏳ deferred | ✅ already correct in active bib (matches PMID 41113708) | verified via PubMed; **no edit needed** |
| **L1** ARI rounding (0.90 vs 0.903) | recheck §2 LOW | ✅ already unified in active scaffold (OUTLINE / R-PRO / CAP / QA all use 0.903 / 0.918) | grep-verified clean; **no edit needed** |
| **L2** Fisher p (6e-10 vs 6.4e-10) | recheck §2 LOW | ✅ unified to 6.4×10⁻¹⁰ / 6.4e-10 | **4 edits applied** (see §3) |
| **L3** BRAF Cohen d (-0.04 vs -0.044) | recheck §2 LOW | ✅ already unified in active scaffold to −0.044 (OUTLINE / R-PRO / CAP) | grep-verified clean; **no edit needed** |
| **L4** references.bib header "~22 papers" | recheck §2 LOW | ✅ no such header line in current bib | **no fix needed** |
| **L5** ★ unicode in Chu2018JMG note | recheck §2 LOW | ✅ Chu2018JMG note already ASCII-only, no ★ | **no fix needed** |
| **L6** Pan2025NatComm pages format | recheck §2 LOW | ✅ pages = {3601} matches Nat Comm article-number convention | **no fix needed** |
| **L7** 4 bib not yet cited | recheck §2 LOW | ⏳ deferred — Discussion drafting needs author voice (§3.1 + Hook + Aim are voice-protected) | **skip per voice rule** |
| **L8** "★ exceptional" verb in F2D caption | recheck §2 LOW | ⏳ deferred — caption polish is author-voice work | **skip per voice rule** |
| **L9** Krishnamoorthy 2025 bib gap | recheck §2 LOW | ⏳ deferred — §3.1 voice-protected; bib decision joins voice-draft turn | **skip per voice rule** |

**Result:** 6 items already clean (M3, M4, L1, L3, L4, L5, L6), 1 item fixed this session (L2), 3 items deferred to author-voice work (L7, L8, L9). Active Paper 1 submission scaffold is now numerically and bibliographically consistent up to the voice-protected boundary.

---

## 2. Verification methodology

For each item, I:

1. **Identified the active source-of-truth** — per `2026_05_04_post_closure_recheck.md` §0, the active Paper 1 submission scaffold lives in `project/reports/2026_05_03_*.md` (OUTLINE, R-PRO, M-PRO, M-SCAF, CAP, SUPP, QA, CHK, STAR, acknowledgments). The parallel `project/manuscript_v8/*.md` directory (May-1 dated) is an older draft; not edited.
2. **Grepped the active scaffold** for the inconsistent forms and the precise forms.
3. **Cross-checked bib entries against PubMed** for M3 / M4 (WebFetch on PMID landing pages).
4. **Applied minimal edits** only where the active scaffold contained a residual inconsistency, never inside voice-protected boundaries (Q9 lines 116-129 in QA, Hook §1.1 first line, Aim §1.4, Discussion §3.1 + §3.4 Limitations, Cover ¶1).

---

## 3. L2 edits applied (4 total)

All four edits are single-line numerical-precision replacements: `6×10⁻¹⁰` → `6.4×10⁻¹⁰` (or the e-format equivalent). No prose rewriting, no semantic change, no claim alteration.

| File | Line | Context | Before | After |
|---|---:|---|---|---|
| `project/reports/2026_05_03_figure_captions_all.md` | 367 | SUPP T6 caption tail | `Fisher p=6×10⁻¹⁰` | `Fisher p=6.4×10⁻¹⁰` |
| `project/reports/2026_05_03_reviewer_QA_consolidated.md` | 99 | Q7 answer (cross-cohort generalization) | `p=6×10⁻¹⁰` | `p=6.4×10⁻¹⁰` |
| `project/reports/2026_05_03_reviewer_QA_consolidated.md` | 135 | Q10 source-evidence bullet | `p=6e-10` | `p=6.4e-10` |
| `project/reports/2026_05_03_reviewer_QA_consolidated.md` | 138 | Q10 answer (Fisher p) | `Fisher p=6×10⁻¹⁰` | `Fisher p=6.4×10⁻¹⁰` |

**Voice check:** Q7, Q10 are not voice-protected (only Q9 is, lines 116-129). Q11 line 142+ untouched — already at 6.4×10⁻¹⁰ where applicable.

**Final grep verification (all `6×10⁻¹⁰` / `6e-10` instances in active Paper 1 scaffold):**

```
$ for f in 2026_05_03_{manuscript_v8_OUTLINE,results_R1_R5_prose,methods_M1_M11_scaffold,methods_M1_M11_prose,figure_captions_all,paper_supp_tables_draft,reviewer_QA_consolidated,PRE_SUBMISSION_CHECKLIST,STAR_Methods_KeyResources,acknowledgments_highlights_etoc}.md; do
    grep -nE "p ?= ?6×10⁻¹⁰|p ?= ?6e-10" "$f"
done
(empty)  ← clean
```

---

## 4. M3 — Lee2024 GSE213647 verification

**Bib entry** (`project/reports/2026_05_03_references.bib:34-45`):
- 37 authors listed by full name (Lee, S. E. through Kang, Y. E.)
- Title: "Unraveling the role of the mitochondrial one-carbon pathway in undifferentiated thyroid cancer by multi-omics analyses"
- Nature Communications 2024, 15(1):1163, DOI 10.1038/s41467-024-45366-0, PMID 38331894

**PubMed PMID 38331894** (WebFetch result):
- Authors: Lee SE, Park S, Yi S, Choi NR, Lim MA, Chang JW, Won HR, Kim JR, Ko HM, Chung EJ, Park YJ, Cho SW, Yu HW, Choi JY, Yeo MK, Yi B, Yi K, Lim J, Koh JY, Lee MJ, Heo JY, Yoon SJ, Kwon SW, Park JL, Chu IS, Kim JM, Kim SY, Shan Y, Liu L, Hong SA, Choi DW, Park JO, Ju YS, Shong M, Kim SK, Koo BS, Kang YE
- (37 authors — exact match)
- Title: exact match
- Journal/year/volume/issue/pages: exact match

**Verdict:** ✅ Bib already correct. The recheck memo's "Y and others" placeholder must have referred to a prior draft already updated; current bib is publication-ready. **No edit applied.**

---

## 5. M4 — Lim2025 GSE286332 verification

**Bib entry** (`project/reports/2026_05_03_references.bib:47-57`):
- 8 authors (Lim, D.-W. through Kim, S.-M.)
- Title: "Dysregulated vitamin D signaling in Hashimoto's thyroiditis: an integrated transcriptomic study in a Korean cohort"
- Frontiers in Endocrinology 2025, 16:1666115, DOI 10.3389/fendo.2025.1666115, PMID 41113708

**PubMed PMID 41113708** (WebFetch result):
- Authors: Lim DW, Jeong HJ, Lee JS, Choi MS, Fang S, Wang JH, Kim H, Kim SM (8 authors — exact match)
- Title: exact match (Hashimoto's-thyroiditis-and-vitamin-D framing confirmed)
- Journal/year/volume/pages: exact match

**Verdict:** ✅ Bib already correct. The recheck memo's "title verify" caveat is satisfied — published title matches our bib record verbatim. **No edit applied.**

**Manuscript-side note:** the published title focuses on Hashimoto's vitamin D signaling, while our use of GSE286332 frames it as a PTC vs PTC+HT sub-analysis. The bib `note = {...}` field already documents this scope distinction (`"Used herein as a sub-analysis (PTC vs PTC+HT) of the same cohort"`). No further action needed.

---

## 6. L1, L3, L4, L5, L6 — already-resolved evidence

### L1: ARI rounding

Active scaffold final state:
- `2026_05_03_manuscript_v8_OUTLINE.md` line 24: `TIERA67 ARI=0.903 ≈ unrestricted top-5000 MAD ARI=0.918`
- `2026_05_03_manuscript_v8_OUTLINE.md` line 146: `TIERA67 ARI=0.903, pan-genome top-5000 MAD ARI=0.918 (cluster definition robust)`
- `2026_05_03_results_R1_R5_prose.md` line 41: `top-5,000... ARI of 0.918 ... TIERA67-based clustering (ARI = 0.903)`
- `2026_05_03_paper_supp_tables_draft.md`: precise (0.903 / 0.918 / 0.864 / 0.902)
- `2026_05_03_reviewer_QA_consolidated.md` Q3 line 47, line 50: precise

The audit-flagged "Outline 0.90 / 0.92 (rounded)" form is **no longer present** in the active scaffold — likely fixed in the HM closure session (`2026_05_04_paper1_HM_closure_report.md`). **No edit needed.**

### L3: BRAF Cohen d

Active scaffold final state:
- `2026_05_03_manuscript_v8_OUTLINE.md` line 24: `BRAF transcript Cohen d=−0.044`
- `2026_05_03_results_R1_R5_prose.md`: `Cohen's d = −0.044`
- `2026_05_03_figure_captions_all.md`: `Cohen's d = −0.044`

The audit-flagged "−0.04 (round)" form is **no longer present** in active scaffold. **No edit needed.**

### L4: bib header "~22 papers"

Current `2026_05_03_references.bib` header (lines 1-5) is style/date/voice metadata only — no paper-count claim. The `03_intro_references.bib` header (lines 1-5) similarly has no count. Audit memo referred to a prior draft header. **No edit needed.**

### L5: ★ unicode in Chu2018JMG note

Current `Chu2018JMG` note (`2026_05_03_references.bib:73`):
```
note = {Han Chinese GD HLA fine-mapping. PMC: PMC6161647. (NOTE) Citation correction: previously misattributed as ``Chen 2018'' in earlier versions of this manuscript.}
```
ASCII-only, no `★` / `☆` / `♦` etc. **No edit needed.**

### L6: Pan2025NatComm pages format

Current entry: `pages = {3601}`. This is the article number, which is Nature Communications' canonical pagination convention (each article gets a unique integer). Verified format is correct. **No edit needed.**

---

## 7. L7-L9 — deferred to author voice (do NOT pre-draft)

| ID | Item | Why deferred | When to handle |
|---|---|---|---|
| **L7** | 4 bib entries (Tuttle2019, Yi2016, Chen2024, Kim2014) not yet cited in body text | Citing them requires Discussion §3.1 / §3.4 / Hook / Aim drafting which is voice-protected per `v17_sprint_vs_marathon_violation` | User voice-draft turn (Hook ¶1 / Aim ¶4 / Disc §3.1 / Limitations) |
| **L8** | "★ exceptional" verb in F2D caption (`2026_05_03_figure_captions_all.md`) | Caption polish, borderline subjective wording — author voice | User caption-pass turn |
| **L9** | Possible Krishnamoorthy 2025 bib gap | §3.1 voice-protected discussion contains the Landa 2016 cite save (per `v17_landa2016_cite_save` memory). Whether Krishnamoorthy 2025 is needed depends on §3.1 final wording | User §3.1 voice draft |

These three items all sit downstream of voice-protected sections. Pre-drafting them would violate `v17_sprint_vs_marathon_violation`. Held for user keyboard.

---

## 8. Marathon discipline attestation

- ✅ No analysis run
- ✅ No data downloaded
- ✅ No new RunPod dispatch
- ✅ Voice-protected sections untouched (Hook ¶1, Aim ¶4, Disc §3.1, §3.4 Limitations, Cover ¶1, Reviewer Q9)
- ✅ Paper 3 (chmod 444) untouched
- ✅ Paper 4 backlog untouched
- ✅ Pod D monitor (`poll_and_recover_pod_D.sh` PID 332788) running independently in background
- ✅ Pod C raw recovery: confirmed unrecoverable (`2026_05_04_pod_C_alphafold_raw.md`); no narrative generated
- ✅ Image-DM1 / TCGA WSI / G1/G2/G3 / Foundation model retry — all blocked per audit §8.2 lock

---

## 9. State of Paper 1 submission scaffold (post-this-session)

| Component | State |
|---|---|
| Numerical consistency (T1) | ✅ all HIGH/MEDIUM closed (HM closure 2026-05-04 morning); all LOW closed except L7/L8/L9 (voice-deferred) |
| Cross-references (T3) | ✅ closed (M2 done) |
| BibTeX + citations (T2) | ✅ all closed at scaffold level (M3, M4 verified PubMed); L7 (4 entries pending body cite) deferred to voice work |
| Supplementary tables (T5) | ✅ closed (M5 S3 XLSX, M6 S8d XLSX delivered) |
| Methods prose (T4) | ✅ baseline scaffold present (M-PRO + M-SCAF) |
| Voice-protected sections | ⏳ owned by user (Hook ¶1, Aim ¶4, Disc §3.1, §3.4 Limitations, Cover ¶1, Q9) |
| K2 cohort handling (Q8 evidence) | 🔄 pending Pod D STAR re-quant (monitor running; ~8-30 hr remaining wall) |

bioRxiv 6/13 target unchanged. Marathon mode discipline intact.

---

## 10. Outstanding user actions

1. **Voice-protected drafting (own keyboard, your slot):**
   - Hook ¶1 (`project/manuscript_v8/04_intro_1_1_hook.md` or live scaffold equivalent) — `v17_sprint_vs_marathon_violation` rule
   - Aim ¶4 (intro §1.4, currently in `project/manuscript_v8/03_introduction.md` line 44 paragraph)
   - Discussion §3.1 (Landa 2016 cite save per `v17_landa2016_cite_save` memory — 3-layer reverse-causality scaffold)
   - Discussion §3.4 Limitations
   - Cover letter ¶1
   - Reviewer Q9 (PTC+HT 18/18 DM2 framing)
2. **L7 (after voice-drafting turn):** weave Tuttle2019, Yi2016, Chen2024, Kim2014 cites into body where appropriate.
3. **L8 (caption pass):** decide whether `★ exceptional` reads as voice OK or replace.
4. **L9 (after §3.1 voice draft):** decide if Krishnamoorthy 2025 bib entry is needed.
5. **Pod D recovery (autonomous):** monitor will SCP + run `post_pod_D_meta3cohort.py` + write `project/reports/2026_05_04_pod_D_k2_star_postprocess.md` whenever Pod D run.log shows ALL DONE; watchdog will fire notification.

No user action is blocking any other action. All deferrals are clean handoff points.

---

## 11. Cross-references

- Audit source: `project/reports/PAPER1_MARATHON_AUDIT_BUNDLE_2026_05_04.md` (chmod 444)
- HM closure: `project/reports/2026_05_04_paper1_HM_closure_report.md`
- Recheck: `project/reports/2026_05_04_post_closure_recheck.md`
- Conflict audit (Pod-related): `project/reports/2026_05_04_paper1_dm1_full_conflict_audit.md`
- Pod state: `POD_CD_RUNNING_STATE_2026_05_04.md`, `project/reports/2026_05_04_pod_C_alphafold_raw.md`, `project/reports/2026_05_04_pod_cd_completion_report.md`
- Voice rule: memory `v17_sprint_vs_marathon_violation`
- Paper 1 = molecular DM1 (memory `paper_numbering_2026_05_04`)

---

*Closed. Active Paper 1 submission scaffold numerically + bibliographically clean to the voice-protected boundary. Pod D monitor proceeding. User keyboard owns voice-protected next turn.*
