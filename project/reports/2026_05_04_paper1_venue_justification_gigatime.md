# Paper 1 — Venue justification memo (GigaTIME publish context + reach venue argument)

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Status:** Internal venue analysis memo for Yu advisor input. NOT manuscript-facing. NOT a recommendation to over-claim. Provides factual context for the venue-ladder decision (Sci Rep base → Cell Rep Med / JCI Insight reach → Cell / Nature Cancer aspiration).

**Authority:** `STATUS_PAPER1_2026_05_02.md` venue ladder (Sci Rep base; Cell Rep Med / JCI Insight reach; Cell / Nature Cancer reach pending Yu split re-evaluation). This memo updates the reach-venue argument with the GigaTIME (Valanarasu et al., 2026 Cell) publish context and the closure battery negative-feasibility evidence.

**Voice rule:** Internal memo. Not voice-protected. Author / Yu advisor read and decide.

---

## 1. The GigaTIME publish context (objective)

| fact | detail | source |
|---|---|---|
| Journal | Cell | OpenAlex |
| Volume / pages | 189(2):386–400.e19 | Cell metadata |
| Date | January 2026 | OpenAlex |
| Citation impact (FWCI) | 10.47 (top 10%) within 4 months | OpenAlex |
| Topic | Multimodal AI for tumor microenvironment, H&E ↔ mIF translation | abstract |
| Cohort scale | 14,256 patients × 51 hospitals × 24 cancers | abstract |
| Microsoft Research × Providence Health | yes | author affiliations |

**Implication for the THCA molecular-axis space**: Cell editor accepted a high-resource, multimodal-AI, population-scale paper in tumor pathology in early 2026. The space is **actively reviewed** at Cell-tier venues. Our Paper 1 — molecular-axis-only at smaller scale but mechanistically rigorous — is not in the same submission lane (we are RNA-axis-on-molecular-evidence, they are image-translation-at-population), but the editor's interest in multimodal tumor pathology suggests adjacent angles are receivable.

## 2. Where Paper 1 stands in the reach-venue argument

### Strengths (in the manuscript-safe envelope per molecular-only lock §2)

- **16 source-verified claims** spanning multi-cohort RNA, ST, clinical outcome, mechanism, drug target, tumor-vs-normal window, dark-matter discrimination, robustness
- **Independent prognostic effect** after BRAF / RAS / age / stage adjustment (DFI multivariate HR = 1.41, p = 0.025)
- **Coherent 4-step mechanism** (TF collapse → DNMT → STAT3/AP-1 → TROP2)
- **Druggable target** (TROP2 / TACSTD2 ADC at tumor-population level)
- **Honest subtractive negatives** (N1 dark-matter NS, Q3 spot-level NEG, H&E NO-GO) — pre-empts referee Q
- **Pre-registered negative feasibility** (closure battery) — bounds what H&E alone can recover, strengthens the "RNA molecular axis is the right tool" argument
- **Venue-adjacent recent work**: Valanarasu 2026 Cell (this memo §1) demonstrates the multimodal-AI lane is publishing at Cell

### Weaknesses (honest)

- **No wet-lab functional validation** (TROP2 IHC, sacituzumab IC50, organoid, xenograft = 0 done)
- **Single primary cohort + bootstrap** for spatial bootstrap CI (n = 12 slides for a key metric)
- **Pan-cancer scope is THCA-specific outlier** — not a pan-cancer story; framing must be thyroid-specialized
- **Image-DM1 angle dropped** — we are not the multimodal-AI lane

## 3. Venue ladder re-assessment

| venue | reasoning | probability (qualitative) | Yu input needed? |
|---|---|---|---|
| **Sci Rep / Endocrine-Related Cancer / Thyroid** | base; molecular axis + clinical association + ST validation packaged cleanly | high (manuscript-safe envelope is sufficient) | confirm |
| **JCI Insight** | reach; clinical association + mechanism + drug target + honest negative + pre-registered closure all align with JCI Insight scope (mechanistic + translational, no Cell-tier scale required) | medium-high (best calibrated reach venue) | confirm |
| **Cell Rep Med** | reach; comparable scope to JCI Insight + somewhat more emphasis on translational impact (sacituzumab framing) | medium-high (alternative reach) | confirm |
| **Cancer Discov** | aspiration; would require either more cohorts OR functional validation; current package is borderline | low-medium (without wet-lab gain, more likely R&R or rejection) | yes |
| **Nat Cancer** | aspiration; would require pan-cancer relevance OR a multimodal AI bridge OR significant wet-lab; current scope is too THCA-specialized | low (without major scope expansion, unlikely) | yes |
| **Cell** | far aspiration; would require all of the above + multimodal AI integration + scale; not currently feasible without years of additional work | very low | yes |

**The honest reading**: our manuscript at submission time targets **JCI Insight or Cell Rep Med as primary reach**, with **Sci Rep / Endocrine-Related Cancer / Thyroid** as base, and **Cancer Discov** only if Yu advises and we can add 1–2 cohorts or functional layer in revision.

GigaTIME's Cell publication does NOT directly elevate our submission lane (we are not GigaTIME-lane), but it does establish that:
(a) the broader thyroid molecular axis space is editorially current
(b) reviewers will likely ask "why didn't you do GigaTIME-style multimodal AI" — for which we have the closure battery + reviewer Q scaffold + 5 re-entry conditions

## 4. Cover-letter-facing framing (already drafted, see S5 scaffold)

`2026_05_04_paper1_cover_letter_para23_scaffold.md` provides Para 2 (molecular-axis contribution) and Para 3 (orthogonal complementarity to GigaTIME). Editor reads our cover letter understanding:
- We are the **molecular axis** lane (RNA-on-mechanism)
- They are the **multimodal AI** lane (H&E-on-image-translation)
- The two are **complementary**, not competitive
- Future work bridges them; current work establishes the molecular ground-truth

## 5. Concrete venue-decision recommendation

**Default submit ladder**:

```
1. JCI Insight (or Cell Rep Med) — primary reach submission target
   ↓ if rejected without R&R
2. Sci Rep / Endocrine-Related Cancer / Thyroid — base submission
   ↓ if accepted
   publish + plan v2 with wet-lab + multimodal AI bridge
   ↓ v2 timeline: 12–18 months post-publication
   target: Cancer Discov OR Cell Rep Med
```

**Yu advisor input requested on**:
1. Should we attempt Cancer Discov first as an even higher reach (low probability but high reward)?
2. Is JCI Insight or Cell Rep Med the better primary reach venue (which has stronger mechanistic / translational fit for our package)?
3. Should Para 3 of cover letter cite GigaTIME explicitly (Option 3A in S5 scaffold), or save the reference for Discussion only?
4. Is the closure battery best published as Paper 1 Methods supplement OR as a standalone short-form companion paper (Paper 2A) at a methods-oriented venue?

## 6. What this memo is NOT

- ❌ A claim that we are at Cell/Nature/Cancer Discov tier (we are not without wet-lab + multimodal AI)
- ❌ A recommendation to write venue-probability percentages anywhere in the manuscript or cover letter
- ❌ A reason to delay submission (current envelope is JCI Insight / Cell Rep Med ready)
- ❌ A justification to revisit image-DM1 work (closure battery NO-GO + 5 re-entry conditions remain canonical)

## 7. What this memo IS for

- ✅ Yu advisor briefing input for venue decision
- ✅ Background context for cover letter Para 2/3 framing (S5 scaffold)
- ✅ Internal record of the multimodal-AI publication landscape as of 2026-05-04
- ✅ Reviewer Q anticipation reinforcement (S2 reviewer Q scaffold Q2 specifically)

---

## Cross-references

- `STATUS_PAPER1_2026_05_02.md` (existing venue ladder)
- `2026_05_04_paper1_dm1_full_molecular_only_lock.md` (manuscript-safe envelope)
- `2026_05_04_paper1_cover_letter_para23_scaffold.md` (S5 cover letter Para 2/3)
- `2026_05_04_paper2a_reviewer_q_gigatime_scaffold.md` (S2 reviewer Q)
- `2026_05_04_paper2a_gigatime_methods_comparison_table.md` (S1 comparison table)
- `2026_05_04_gigatime_thca_relevance_scan.md` (S3 scan)
- `2026_05_04_paper2a_reentry_protocol_v0.md` (S4 future-work design)
- Memory: `v17_npj_ship_status` (Sci Rep base context), `v17_dark_matter_pivot_2026_04_29` (8-gene framing)
