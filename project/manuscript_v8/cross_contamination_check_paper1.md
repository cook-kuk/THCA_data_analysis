---
title: "Paper 1 cross-contamination check — forbidden term grep + boundary audit"
date: 2026-05-04
author: Seungho Cook
scope: Cancer-only Paper 1 (BRAF/RAS-negative PTC dark matter — fusion + epigenetic)
files_audited:
  - project/manuscript_v8/00_title_candidates.md (75 lines)
  - project/manuscript_v8/01_abstract.md (73 lines)
  - project/manuscript_v8/02_outline.md (328 lines)
forbidden_terms: [Hashimoto, HT, thyroiditis, lymphocytic, Graves, GD, TSAb, hyperthyroid, autoimmune-overlap, autoimmune-PTC, Hashimoto-overlap, "DPB1*05:01", PTC+HT, "PTC ≠ GD"]
allowed: general immune microenvironment language (e.g. "BRAF V600E HIGHER HLA-I" counter-intuitive finding); "immune-overlap subtype" phrase as Paper 2 hook
status: draft for advisor review
---

# Cross-contamination check — Paper 1

## 0. Executive summary

| File | Hits | HARD violations | SOFT violations (meta only) | Pass |
|---|---|---|---|---|
| 00_title_candidates.md | 2 | 0 (titles themselves clean) | 2 (Cand 3 commentary cells) | Title Cand 1 ★ — clean |
| 01_abstract.md | 1 (+ 1 boundary risk) | 0 (abstract body clean) | 1 (review-focus note item 3) | Abstract body PASS |
| 02_outline.md | 4 + 6 boundary | **3 HARD** (§2.4b stat + §2.5 cohort + §3.3 Pan-Asian) | 1 (Risk-flag table) | §2.4b/2.5/3.3 require fix |

**Bottom line** — abstract body and Title Candidate 1 are scope-clean. **Outline §2.4b, §2.5, §3.3, §3.4(iii)** carry HARD violations that will propagate into manuscript Results / Discussion if not pruned now. Pillar 5 boundary statement at line 104 ("Paper 1 stops at 'older + immune-hot DM1 sub-B is mechanistically distinct'") is correct in *intent*, but the same row still inlines Paper 2 statistics (Hashimoto-like %, HLA n=874, BCR clonal) — the boundary needs to be enforced by *removing* the Paper 2 facts from the Paper 1 row, not just labelling them.

---

## 1. Title check

### Candidate 1 (★ recommended) — line 24
> "Fusion-driven, epigenetically silenced dark matter in BRAF/RAS-negative papillary thyroid carcinoma"

- forbidden-term count: **0**
- "dark matter" = Xing 2014 NEJM term — cancer-context, scope-clean
- ✅ **PASS** — proceed with Cand 1 as the v3 title.

### Candidate 2 — line 40
> "An 8-gene RAI panel enables reflex fusion testing and re-induction therapy in BRAF/RAS-negative thyroid cancer"

- forbidden-term count: **0**
- ✅ PASS

### Candidate 3 — line 56
> "An 8-gene panel sub-stratifies BRAF/RAS-negative thyroid cancer into fusion-driven and immune-overlap subtypes"

- title body: 0 forbidden terms ("immune-overlap" allowed per scope as Paper 2 hook)
- ⚠️ commentary lines 62 + 65 contain forbidden terms (`autoimmune-PTC`, `Hashimoto-DM2 axis`, `GSE286332`)
  - line 62 (Foregrounds cell): `Pillar 5 autoimmune-PTC layer`
  - line 65 (Pros cell): `"Immune-overlap" foreshadows Pillar 5 (Hashimoto-DM2 axis + GSE286332)`
- These are working-doc analysis cells, not destined for the rendered paper. **Status: SOFT.** Recommend sanitizing for clean working-doc archive.

**Recommended sanitization (lines 62, 65)** — replace meta-commentary text:
- L62: `Pillar 5 autoimmune-PTC layer` → `R4-2 immune-overlap teaser layer (Paper 2 backbone reserved)`
- L65: `"Immune-overlap" foreshadows Pillar 5 (Hashimoto-DM2 axis + GSE286332)` → `"Immune-overlap" phrase reserves Pillar 5 main story for Paper 2`

---

## 2. Abstract check

### Abstract body (lines 11–19)
- Background / Methods / Results / Conclusions paragraphs: forbidden-term count **0**
- Conclusions (b) verified — line 19:
  > "DM1 is a fusion-driven, epigenetically silenced subtype, providing a clinical sub-stratification algorithm and motivating evaluation of fusion-targeted therapy and **epigenetic-targeted RAI re-induction**."
  - mentions HMA + RAI re-induction *rationale* via "epigenetic-targeted RAI re-induction" phrase — autoimmune mechanism not invoked
  - ✅ **PASS** — Conclusions (b) is scope-clean and reads as cancer-mechanism rationale, not autoimmune-overlap rationale.

### Boundary risk — line 15 (Methods statement) + line 49 (numerical claims table)
- L15: `Korean cohorts (n=874)` — number itself is fine
- L49: `Korean cohorts n=874 | D4-P1: K2(235) + Lee(630) + GSE286332-PTC(9)` — **the n=874 includes 9 patients from GSE286332**, which is per memory `v17_gse286332_strong_go.md` a Korean PTC vs PTC+HT cohort = Paper 2 territory.
- **Recommendation** — confirm with advisor whether to (a) keep GSE286332-PTC(9) in n=874 with STAR Methods caveat that "only the PTC arm of GSE286332 is included; the PTC+HT arm is reserved for Paper 2", or (b) recompute n=865 (K2+Lee only) and move GSE286332 entirely to Paper 2.
  - Option (b) is cleaner for cancer-only scope but loses 9 patients of statistical power.
  - Option (a) keeps the n but adds 1 line of STAR Methods text. Lower risk if advisor accepts.

### Annotation note — line 66 (review-focus item 3)
- L66: `**Pillar 5 (autoimmune-PTC) 는 abstract 에서 의도적 보류 — Paper 2 의 backbone 으로 reserve.**`
- ⚠️ forbidden term `autoimmune-PTC` in working note. Not in abstract proper.
- **Status: SOFT.** Not rendered in final paper if the abstract block (lines 11-19) is the only excerpted content.
- **Recommended sanitization** — replace `Pillar 5 (autoimmune-PTC) 는 abstract 에서 의도적 보류 — Paper 2 의 backbone 으로 reserve` with `Paper 2 sub-stratification axis (immune-overlap) is intentionally deferred — Paper 2 backbone reserve`.

---

## 3. Outline §2.4 split check

### §2.4a — line 103
> "DM1 epigenetic silencing of differentiation machinery"
- claims: TPO d=2.30, DIO1 d=1.24, TSHR d=1.20, PAX8 d=0.97, TG d=0.86, FOXE1 d=0.84, NKX2-1 d=0.63, SLC5A5 NS; mean β 0.385 vs 0.253; fusion-independent
- forbidden-term count: **0**
- ✅ **PASS** — pure epigenetic mechanism, fusion-independence noted, scope-clean.

### §2.4b — line 104 ⚠️ HARD VIOLATION
> "Fusion-negative DM1 = immune-overlap subtype (teaser) | P5 partial | Sub-A vs sub-B: age (37 vs 51), fusion% (84.7 vs 57.9), stage III/IV (15.3% vs 44.4%), CD8/IFN-γ low vs high (d=-0.5 to -0.6), **Hashimoto-like 40% vs 67% (OR 0.34, p=0.064 trend)**; **Pillar 5 mediation/HLA n=874/BCR clonal/sub-B NBNR is Paper 2 backbone — Paper 1 stops at 'older + immune-hot DM1 sub-B is mechanistically distinct'**"

**Violations (HARD — these statistics will propagate into Results §2.4b prose):**
1. `Hashimoto-like 40% vs 67% (OR 0.34, p=0.064 trend)` — reports a Paper 2-territory stratification statistic in a Paper 1 results row
2. `Pillar 5 mediation` — references Paper 2's mediation analysis as if it might be brought into Paper 1
3. `HLA n=874` — Paper 2's HLA cohort, mentioned as an inline statistic
4. `BCR clonal` — Paper 2's BCR analysis
5. `sub-B NBNR` — NBNR is K2-cohort label per memory, Paper 2 territory

**Boundary statement check** — the row's tail does correctly say `"Paper 1 stops at 'older + immune-hot DM1 sub-B is mechanistically distinct'"`. Intent is right, **but listing the Paper 2 statistics in the same cell will encourage them to bleed into Results §2.4b prose during drafting.**

**Recommended §2.4b row content (advisor approval before edit)** — replace the "Key claims" cell with cancer-only sub-A vs sub-B contrasts:
- Drop: `Hashimoto-like 40% vs 67% (OR 0.34, p=0.064 trend)` (Paper 2)
- Drop: `Pillar 5 mediation/HLA n=874/BCR clonal/sub-B NBNR` mention list
- Keep: age (37 vs 51), fusion% (84.7 vs 57.9), stage III/IV (15.3% vs 44.4%), CD8/IFN-γ low vs high d=−0.5 to −0.6
- Keep boundary tail: `"Paper 1 stops at 'older + immune-hot DM1 sub-B is mechanistically distinct'; the autoimmune-overlap mechanism is reserved for Paper 2."` ← this last clause should be the only place "Paper 2" appears in this row, and replaces all the Paper 2 statistic listings.

**§2.4b take-home (line 112) — clean**
> "...with fusion-negative DM1 representing a mechanistically distinct immune-overlap subtype warranting separate inquiry." → Paper 2 hook
- "immune-overlap subtype" + "warranting separate inquiry" = scope-allowed framing
- "→ Paper 2 hook" annotation is meta-commentary, will not enter prose
- ✅ **PASS** — §2.4b take-home is the boundary as advisor specified.

### §2.5 — line 106 ⚠️ HARD VIOLATION
> "Cross-cohort + Reflex testing algorithm | P1 + (translational) | sc external (Pu 2021 r=0.798-0.886; Lu 2023 thyrocyte-intrinsic); FFPE robust (KS p=0.44); **Korean K2/Lee/GSE286332 n=874 + GSE213647 n=632 (28% Hashimoto-like)**; **DM1 captures 81.8% TCGA RET+ → reflex algorithm**; population estimate 48 selpercatinib-eligible per 1000 PTC | Figs 3, 5"

**Violations (HARD):**
1. `(28% Hashimoto-like)` — Hashimoto-like proportion of GSE213647 = Paper 2/3 territory
2. `GSE286332` — by memory `v17_gse286332_strong_go.md` this is the PTC vs PTC+HT cohort = Paper 2 backbone

**Recommended §2.5 row content:**
- Drop: `(28% Hashimoto-like)` parenthetical
- For `GSE286332`: align with the abstract resolution (option (a) or (b) above). If option (a), keep `GSE286332` in cohort list with a STAR Methods note "PTC arm only". If option (b), drop entirely → `Korean K2/Lee n=865 + GSE213647 n=632`.

### §3.3 Discussion — line 122 ⚠️ HARD VIOLATION
> "East-Asian generalizability + Pu 2021 cross-cohort | 200-300 w | 4,300+ East Asian PTC tumors covered; Korean DM 37.8% vs TCGA 28.4%; K2 NBNR mixed phenotype (vascular invasion); **Pan-Asian Hashimoto-like replication (Korean GSE213647 22-28% ≈ TCGA 18-20%)**; ..."

**Violation:**
- `Pan-Asian Hashimoto-like replication` is not generalizability evidence for Paper 1 (cancer-only). It is generalizability evidence for Paper 2's autoimmune-overlap thesis.

**Recommended §3.3 replacement** — drop the Hashimoto-like sentence entirely. The DM1 prevalence consistency claim (Korean DM 37.8% vs TCGA 28.4%) is sufficient by itself for East-Asian generalizability, and the cancer-only scope is preserved.

### §3.4 Limitations (iii) — line 123 ⚠️ MIXED
> "(iii) Methodological: HLA d 1.5+ partial autocorr (Δd=0.30); K2 mini-index calibration mismatch (R4-4) propagates GSE286332 reference."

- `HLA d 1.5+` — if Paper 1 does not include HLA results in §2 Results, this limitation has no §2 anchor and reads as dangling. Either remove (cleanest), or rephrase to address only the panel autocorrelation analysis (no HLA term).
- `propagates GSE286332 reference` — Paper 2 territory.

**Recommended §3.4(iii) replacement:**
- Drop the HLA clause entirely (no §2 anchor in cancer-only Paper 1).
- Replace `propagates GSE286332 reference` with a generic phrasing such as: `K2 mini-index calibration mismatch (R4-4) — same calibration is propagated through any cohort built on the K2 reference; cross-Korean cohort findings are interpreted with this caveat.` (Avoids naming GSE286332.)

### §2.4b take-home + risk-flag table (lines 112, 250)
- Line 112 (take-home) — ✅ scope-clean per §2.4b check above.
- Line 250 (risk-flag table) is internal documentation; the row "Paper 2 backbone reserve: K2 + Pan-Asian HLA n=874 + GSE286332 mediation + BCR/TLS + sub-B NBNR" is **planning-doc level**, not destined for paper. SOFT only.

---

## 4. Discussion 3.1 (Landa 2016 paired-cancer continuum) — line 120

> "L1 genetic (76.8% fusion+) · L2 mechanism heterogeneity (sub-A vs sub-B) · L3 epigenetic (TPO d=2.30 fusion-independent); ★ **Landa et al 2016 JCI 126(3):1052-1066** characterized advanced thyroid cancer (84 PDTC + 33 ATC) genomic landscape ... profound suppression of differentiation transcripts (TG, TSHR, TPO, PAX8, SLC26A4, DIO1, DUOX2); DM1 framework "**identifies an upstream signature consistent with this dedifferentiation trajectory**" — same machinery silenced at ATC end already epigenetically attenuated in DM1 PTCs. ... **Bradley 2010 contradict**: BRAF V600E HIGHER HLA-I (counter-intuitive, R3 finding 14) — alternate immune evasion mechanism in BRAF-driven PTC"

- forbidden-term count: **0** (HLA-I mention is the explicitly-allowed counter-intuitive cancer-context observation per scope)
- Landa 2016 PDTC→ATC continuum framing = cancer paired-trajectory, scope-clean
- ✅ **PASS** — Discussion 3.1 framing is scope-clean. Convergence framing (8-gene panel via Yoo 2016 + Landa 2016 ATC silenced gene list) addresses reverse-causality risk without invoking autoimmune mechanism. Bradley 2010 HLA-I claim is counter-intuitive cancer immune-microenvironment finding (BRAF V600E HIGHER HLA-I), explicitly allowed.

---

## 5. Discussion 3.2 (ATA 2015/2025 + LIBRETTO-001) — line 121

> "**ATA 2015 (Haugen 2016) — recently updated as ATA 2025 (Ringel 2025) — incorporates BRAF V600E as sole molecular risk modifier** (fusion drivers and epigenetic silencing not reflected in current risk tiers); BRAF/RAS-negative dark matter (~23%) currently defaults to clinico-pathological-only stratification; **DM1 reflex algorithm = orthogonal molecular axis** providing ATA intermediate-risk tier additional input. 57% DM1 actionable (RET 40%, NTRK 12%, ALK 5%, BRAF 6%). DM1 captures 81.8% TCGA RET+; 48 selpercatinib-eligible/1000 PTC. **Selpercatinib FDA accelerated approval May 2020 based on LIBRETTO-001 (Wirth et al 2020 NEJM, ORR 79% in RET-fusion+)** — DM1 RNA-score positivity provides upstream candidate pool for prospective evaluation. HMA + RAI re-induction rationale: TPO/DIO1/TSHR promoter hypermethylation → epigenetic-targeted RAI re-induction (decitabine retrospective ATC trials NCT00085293, NCT01065090); SLC5A5/NIS exception → combination strategies (HMA + lithium for NIS membrane trafficking). **Evaluation, not trial.**"

- forbidden-term count: **0**
- ATA 2015 / ATA 2025 / LIBRETTO-001 / Wirth 2020 NEJM — all cancer-context citations
- HMA + RAI re-induction rationale derives from epigenetic mechanism (TPO/DIO1/TSHR promoter hypermethylation) — autoimmune mechanism not invoked
- ✅ **PASS** — Discussion 3.2 framing is scope-clean. "Evaluation, not trial." dial-back phrase is preserved per advisor honest-tone policy.

---

## 6. Pillar 5 boundary clarity — Paper 1 vs Paper 2

### Boundary phrase locations (Paper 1 hook → Paper 2 reserve)
| File:Line | Phrase | Direction | Status |
|---|---|---|---|
| 02_outline.md:31 | "Pillar 5 main story (HLA n=874 mediation, BCR/TLS, sub-B NBNR) 는 Paper 2 reserve. Paper 1 에서는 fusion- DM1 = older + immune-hot 까지만 teaser. Paper 2 differentiation 보호." | meta-commentary | ✅ explicit, planning-doc only |
| 02_outline.md:56 | "2.4b Fusion- DM1 = immune-overlap teaser ... Pillar 5 hint · Paper 2 reserve" | flowchart annotation | ✅ planning-doc |
| 02_outline.md:104 | "Paper 1 stops at 'older + immune-hot DM1 sub-B is mechanistically distinct'" | inline boundary statement | ⚠️ correct wording but adjacent statistics violate (see §3 above) |
| 02_outline.md:112 | "...with fusion-negative DM1 representing a mechanistically distinct immune-overlap subtype warranting separate inquiry." → Paper 2 hook | take-home | ✅ scope-clean |
| 02_outline.md:250 | "Pillar 5 Paper 1 vs Paper 2 boundary 모호 ... 2.4b 마지막 sentence ... → Paper 2 hook 명시. Paper 2 backbone reserve: K2 + Pan-Asian HLA n=874 + GSE286332 mediation + BCR/TLS + sub-B NBNR" | risk-flag table | ✅ planning-doc only |

### Verdict
- The boundary intention (Paper 1 stops at "older + immune-hot DM1 sub-B") is documented in 4 places (line 31, 56, 104, 112).
- The boundary is **enforceable in the take-home at line 112** ("warranting separate inquiry") which reads cleanly.
- The boundary is **leaky in the §2.4b key-claims cell at line 104** because Paper 2 statistics (Hashimoto-like %, HLA n=874, BCR clonal, NBNR) are listed adjacent to the Paper 1 sub-A vs sub-B cancer contrast.
- If the leaky cell at line 104 is sanitized (as recommended in §3 above), the boundary becomes clean and §2.4b prose drafting will not pull Paper 2 facts.

---

## 7. Action items (advisor approval gate)

Before any §2 prose drafting on §2.4b / §2.5 / §3.3 / §3.4(iii), confirm with advisor:

| # | File:Line | Action | Severity |
|---|---|---|---|
| A1 | 02_outline.md:104 | Drop 4 Paper 2 inline statistics from §2.4b cell (Hashimoto-like %, Pillar 5 mediation, HLA n=874, BCR clonal, NBNR mention list); keep age / fusion% / stage / CD8 contrasts + boundary tail | HARD |
| A2 | 02_outline.md:106 | Drop `(28% Hashimoto-like)`; resolve GSE286332 inclusion question (option a STAR Methods caveat vs option b drop n=865 only) | HARD |
| A3 | 02_outline.md:122 | Drop "Pan-Asian Hashimoto-like replication" sentence from §3.3 | HARD |
| A4 | 02_outline.md:123 | Drop HLA d 1.5+ partial autocorr clause (no §2 anchor); rephrase GSE286332 reference to generic K2-calibration caveat | HARD |
| A5 | 01_abstract.md:15, 49 | Confirm GSE286332-PTC(9) inclusion in n=874 — cohere with A2 decision | MED |
| A6 | 01_abstract.md:66 | Sanitize working note item 3 — replace `Pillar 5 (autoimmune-PTC)` with `Paper 2 sub-stratification axis (immune-overlap)` | SOFT |
| A7 | 00_title_candidates.md:62, 65 | Sanitize Cand 3 commentary cells (replace `autoimmune-PTC` and `Hashimoto-DM2 axis + GSE286332`) | SOFT |

After advisor confirms A1–A4 (HARD), §2.4b / §2.5 / §3.3 / §3.4(iii) prose drafting may proceed within Paper 1 scope.

---

## 8. Files NOT touched in this check
Per session scope, the following files were NOT audited and are out-of-scope for this report:
- `04_results.md`, `06_discussion.md`, `07_star_methods.md`, `08_cover_letter.md`, `09_reviewer_qa.md` — Paper 1 prose files (advisor-protected voice)
- `_audit_HT_vs_GD_2026_05_02.md` — Paper 2/3 boundary doc, separate session
- All `_for_web_claude_*.md` and `_prep_*.md` — review prep, not source-of-truth

If A1–A7 are accepted, a follow-up grep on the prose files (04_results.md, 06_discussion.md, 07_star_methods.md) under the same forbidden-term list is recommended before the next manuscript pass.
