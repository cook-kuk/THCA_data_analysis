---
title: "Paper 1 manuscript v8 — self-verification report"
date: 2026-05-04
author: Seungho Cook
status: post-cleanup verification for Paper 1 portfolio prompt compliance
word_count: 560
---

# Self-verification report

## Files updated in this pass

- [00_title_candidates.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/00_title_candidates.md)
- [03_introduction.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/03_introduction.md)
- [04_intro_1_1_hook.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/04_intro_1_1_hook.md)
- [01_abstract.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/01_abstract.md)
- [04_results.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/04_results.md)
- [05_figure_captions.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/05_figure_captions.md)
- [06_discussion.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/06_discussion.md)
- [07_star_methods.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/07_star_methods.md)
- [08_cover_letter.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/08_cover_letter.md)
- [09_reviewer_qa.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/09_reviewer_qa.md)
- [10_full_manuscript_compiled.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/10_full_manuscript_compiled.md)
- [13_supplementary_tables.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/13_supplementary_tables.md)

## Cross-paper boundary check

Command run:

```bash
rg -n "autoimmune-PTC|autoimmune PTC|Hashimoto|Hashimoto-like|PTC\+HT|HLA-II|BCR|TLS|AICDA|Chu|Graves|GD|TSAb|TSI|hyperthy|thyrotox|exophthal|thyroid eye disease|Chen 2018" \
  project/manuscript_v8/00_title_candidates.md \
  project/manuscript_v8/01_abstract.md \
  project/manuscript_v8/03_introduction.md \
  project/manuscript_v8/04_intro_1_1_hook.md \
  project/manuscript_v8/04_results.md \
  project/manuscript_v8/05_figure_captions.md \
  project/manuscript_v8/06_discussion.md \
  project/manuscript_v8/07_star_methods.md \
  project/manuscript_v8/08_cover_letter.md \
  project/manuscript_v8/09_reviewer_qa.md \
  project/manuscript_v8/13_supplementary_tables.md
```

Result: `0 hits` in the audited Paper 1 draft set. `rg` exited with code `1`, which in this case means no matches were found.

## Citation check

- `Chen 2018`: `0 hits`
- `Chu 2018` / `Chu et al. 2018`: `0 hits` in Paper 1 primary draft files after cleanup
- No new external literature was introduced beyond citations already anticipated in the draft set
- Citations retained in edited passages remain existing manuscript citations: `Yoo et al. 2016`, `Landa et al. 2016`, `Haugen et al. 2016`, `Ringel et al. 2025`, `Wirth et al. 2020`, `Pu et al. 2021`
- Reference source file remains [03_intro_references.bib](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/03_intro_references.bib)

## Voice protection check

Protected sections were left as explicit author placeholders rather than Codex prose:
- [03_introduction.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/03_introduction.md): `1.1` hook first line
- [04_intro_1_1_hook.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/04_intro_1_1_hook.md): author hook anchor
- [06_discussion.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/06_discussion.md): `3.1` opening paragraph, `3.4` limitations
- [08_cover_letter.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/08_cover_letter.md): paragraph 1
- [09_reviewer_qa.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/09_reviewer_qa.md): `Q9`

## Open items for author review

- Citation placeholders in [03_introduction.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/03_introduction.md) still need literature-level verification for the SEER-incidence and Bethesda statements
- ~~Decide whether `GSE286332` should remain as a small Korean reference arm inside the `n=874` aggregate or be fully removed from Paper 1 cohort summaries~~ — **Resolved 2026-05-07 (audit P1-2):** GSE286332-PTC(9) dropped from Paper 1 aggregate; new Korean cohort summary statistic is `n=865 (K2 235 + Lee 630)`. GSE286332 reference arm description retained in STAR Methods with explicit boundary note. Rationale: preserve scope separation from Paper 2 (GSE286332 = Paper 2 PTC vs PTC+HT main cohort).
- Confirm the final public repository URL and archival DOI language in [07_star_methods.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/07_star_methods.md)
- Fill affiliations and corresponding-author email in [08_cover_letter.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/08_cover_letter.md)
- Final supplement numbering can be re-tuned once figure build is complete

---

## Re-verification 2026-05-13 (pre-Yu-meeting audit)

Re-ran the boundary + cohort + citation + figure-ref + key-number audit ahead of the 2026-05-14 Yu meeting. All checks GREEN.

- **Paper 2 boundary**: 0 leak. All `Hashimoto-like` / `HLA-II` / `BCR` hits in main draft files are explicitly boundary-marked as "reserved for Paper 2" (05_figure_captions.md:96/98/142, 06_discussion.md:49, 09_reviewer_qa.md:80/84).
- **Cohort number**: `n=865 (K2 235 + Lee 630)` consistent in 01_abstract.md, 03_introduction.md, 04_results.md, 07_star_methods.md. No lingering `n=874`. GSE286332-PTC(9) confirmed dropped per 2026-05-07 P1-2 audit.
- **Placeholders**: 0 inline `TBD/TODO/XXX/FIXME` in main prose. Two structured TODOs intentionally left in 07_star_methods.md:39 and :49 awaiting Zenodo DOI + public repo URL at submission. One `VERIFY` line in 01_abstract.md:18 awaits Yu-Hyeongwon institutional email.
- **Figure numbers**: Main Fig 1-8, Supp S1-S6, S5b, SX, SX_v13, SX_v14, SB2 all defined with consistent cross-refs.
- **Key effect sizes confirmed across sections**: Fusion 76.8% / OR 7.41 (5 files); pooled HR 2.53 [1.31, 4.89] (5 files); TPO d=2.30 (7 files); mean β 0.385 vs 0.253 (5 files); ARI 0.49/0.90/0.92 (4 files).

### Edits applied 2026-05-13 (post-audit)

- `04_results.md`: section order reordered 2.1 → 2.3 → 2.4a → 2.4b → 2.2 → 2.5 ⇒ **2.1 → 2.2 → 2.3 → 2.4 → 2.5a → 2.5b** per 17_results_reorder_plan.md (selection-layer paper framing, discovery-first). +14 words = header comment only; no prose added/removed/reworded. 1 internal cross-ref `(2.4a)` → `(2.5a)`.
- `07_star_methods.md`: cohort drift fixed (K2 260 → 235, Lee 632 → 630, lines 59-61 and 141); 2 placeholder TODOs added for Zenodo DOI + repo URL at submission.
- `09_reviewer_qa.md`: Q13/Q14 Paper 2 territory drift removed ("HT/B-cell route" → "MAPK-low route" / "non-MAPK route"; GSE286332 dual-role parenthetical added in Q14(ii)). Q9 untouched (voice-protected). Author confirmation requested: is GSE286332 OK to keep in v14 forest pool (n=1,287) while dropped from primary Pillar I cohort (n=865)?

### Open items still requiring author keyboard (voice-protected)

- 03_introduction.md `1.1` hook + 04_intro_1_1_hook.md anchor (Alt B refined / hybrid "compass" decision)
- 06_discussion.md `3.1` opening (Landa 2016 JCI cite, Krishnamoorthy 2025 misattribution 정정) + `3.4` limitations
- 08_cover_letter.md paragraph 1
- 09_reviewer_qa.md Q9
