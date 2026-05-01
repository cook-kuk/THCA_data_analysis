# Paper 1 manuscript v8 — directory README

**Project.** "An 8-gene panel reveals fusion-driven, epigenetically silenced dark matter in BRAF/RAS-negative thyroid cancer"

**Author.** Seungho Cook (1st) · Yu Hyeong-won (corresponding)

**Target venue.** Cell Reports Medicine (1순위) → JCI Insight + Nat Commun dual stretch → npj Precision Oncology (fallback)

**Status (2026-05-01).** v3 outline + first-pass drafts of all 7 prompts (Title/Abstract/Outline/Intro/Results/Figures/Discussion/Methods/Cover Letter/Reviewer Q&A) complete. 본인 voice 적용 + cite verify (5/4 read 후) + 공저자 review (W6) → bioRxiv preprint.

---

## Files

### Manuscript draft files (00-09)

| # | File | Content | Words |
|---|---|---|---|
| 00 | `00_title_candidates.md` | 3 후보 + Hybrid Cand 1 ★ + verb options | — |
| 01 | `01_abstract.md` | Cell Press structured Abstract (Conclusions b 채택) | 153 |
| 02 | `02_outline.md` | Single-page outline (Mermaid + bullets, 5/4-6/13 schedule) | 4,805 |
| 03 | `03_introduction.md` | Section 1 Intro (1.1-1.4) | 715 |
| 03b | `03_intro_references.bib` | BibTeX references (17 cites) | — |
| 04 | `04_results.md` | Section 2 Results (5 sub-results, 대안 A ordering) | 3,250 |
| 05 | `05_figure_captions.md` | Figure 1-8 captions + Suppl S1-S9 | — |
| 06 | `06_discussion.md` | Section 3 Discussion (3.1-3.4) + Limitations | 1,380 |
| 07 | `07_star_methods.md` | STAR Methods 5 sub-sections | — |
| 08 | `08_cover_letter.md` | Cell Rep Med cover letter + 5 suggested reviewers | — |
| 09 | `09_reviewer_qa.md` | 12 reviewer Q&A pre-empt | — |
| 10 | `10_full_manuscript_v1.md` | **Compiled full manuscript** (Title + Abstract + Intro + Results + Discussion + Methods + Refs) | ~5,514 |
| 11 | `11_fig7_dm1_mechanism.py` | Figure 7 code (matplotlib, 5 panels) | — |
| 12 | `12_fig8_epigenetic.py` | Figure 8 code (matplotlib, 3 panels) | — |
| 13 | `13_supplementary_tables.md` | Suppl Tables S1-S10 index + source TSVs | — |

### Prep files (5/4 본인 read 후 verify)

| File | Content |
|---|---|
| `_prep_01_hook_sources.md` | Hook stat 검증 + 3 alt phrasings + Alt B refined |
| `_prep_02_ata2015_cheatsheet.md` | ATA 2015 risk strat + ATA 2025 (시나리오 1/2/3) + Discussion 3.2 phrasings |
| `_prep_03_krishnamoorthy_summary.md` | ★ Landa 2016 JCI cite + dial-back framing + reverse-causality 3-layer 차단 + 5 추가 cites |
| `_prep_04_reading_urls.md` | 본인 5/4 read URL 우선순위 (10 papers) |

### Web Claude review briefs

| File | Status |
|---|---|
| `_for_web_claude_prompt1_review.md` | v1 outline review (1차) |
| `_for_web_claude_prep_review.md` | prep work review (2차) |
| `_for_web_claude_full_v1_review.md` | full v1 manuscript review (3차, 다음 step) |

---

## Schedule (5/4-6/13, 6 weeks marathon)

| Week | Focus | Deliverable |
|---|---|---|
| 5/3 EOD | Citation correction + memory entry | ✅ 완료 |
| 5/4 (Mon) | 본인 read (Landa, Yoo, ATA, Pu, Bradley) + Yu professor 미팅 약속 | ATA cheatsheet 검증, cite verify |
| 5/5-5/10 W1 | v2.5 미세 조정 + Prompt 1 confirm + Prompt 2 시작 | 03_introduction v2 |
| 5/11-5/17 W2 | Intro 마무리 + Results 2.1-2.2 | 04_results v2 partial |
| 5/18-5/24 W3 | Results 2.3-2.5 + Figure plan | 04_results v2 final + 05_figures |
| 5/25-5/31 W4 | Figure code + Discussion + Limitations | 11_fig7 + 12_fig8 + 06_discussion v2 |
| 6/1-6/7 W5 | STAR Methods + Cover letter | 07_star_methods v2 + 08_cover_letter v2 |
| 6/8-6/13 W6 | 공저자 review + revision + bioRxiv | bioRxiv preprint |

---

## Marathon mode rule (5/4 부터)

- Pillar 1 STRONG = 분석 끝 (Memory: `v17_marathon_mode_post_pillar1`)
- D9, D10, D11 분석 sprint 없음
- 새 분석 제안 시 paper-blocking 만 진행, 그 외 Phase 1 (Bundang Graves' paper) 또는 Phase 2 (DIAL audit method paper) reserve
- 매일 1-2 단락 본인 voice — Claude Code = draft tool

---

## Critical findings (2026-04-30 → 05-01)

1. **Krishnamoorthy 2025 Nat Comm misattribution catch** — 진짜 cite 는 Landa 2016 JCI 126(3):1052-1066 (Krishnamoorthy GP 9번째 공저자). 8-gene 5/8 ATC silenced gene list overlap → reverse-causality 3-layer 차단 strategy. (Memory: `v17_landa2016_cite_save`)

2. **Chen 2018 → Chu X et al. 2018 정정** — Paper 2 Pillar 1 Han Chinese GD forest meta cite. (Memory: `v17_paper2_pillar1_forest_strong`)

3. **ATA 2025 published** — PMID 40844370, 본 paper Discussion 3.2 dual cite (ATA 2015 + ATA 2025).

4. **Hook stat 검증** — "5-15%" → "intermediate-risk ~20%" specific, ATA framing thread.

5. **5 추가 cite identified** — Yoo 2016, TCGA 2014, Pu 2021, Bradley 2010, Wirth 2020.

---

## 본인 voice 보호 영역 (W1-W6 본인 키보드)

- 1.1 Hook 첫 줄 (Alt B refined vs hybrid "compass" 본인 결정)
- 1.4 Aim paragraph (paper identity)
- 3.1 첫 paragraph (mechanism story)
- 3.4 Limitations (정직 disclosure 정신)
- Cover letter Para 1 (paper 진짜 의미)
- Reviewer Q9 (mechanism story)
