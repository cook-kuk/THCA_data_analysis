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
- Decide whether `GSE286332` should remain as a small Korean reference arm inside the `n=874` aggregate or be fully removed from Paper 1 cohort summaries
- Confirm the final public repository URL and archival DOI language in [07_star_methods.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/07_star_methods.md)
- Fill affiliations and corresponding-author email in [08_cover_letter.md](/home/seungho/personal/THCA_data_analysis/project/manuscript_v8/08_cover_letter.md)
- Final supplement numbering can be re-tuned once figure build is complete
