# Terminology correction — 2026-05-04

**Trigger:** Yu advisor scope (2026-05-04) — Paper 2 isolated to **Hashimoto-overlap PTC** only.
GD phenotype/mechanism is Paper 3 territory. The cross-disease "autoimmune-thyroid continuum"
framing conflates HT (Paper 2) with GD (Paper 3) and must be removed from Paper 2 context.

## Substitution rules

| Forbidden phrase | Replacement |
|---|---|
| `autoimmune-PTC` / `autoimmune PTC` | `Hashimoto-overlap PTC` |
| `autoimmune-thyroid axis` | `Hashimoto-thyroid overlap axis` |
| `Pan-Asian autoimmune-thyroid susceptibility ...` | `Pan-Asian thyroid HLA susceptibility (HT context)` |
| `autoimmune-thyroid 분야` | `thyroid HLA susceptibility 분야` |
| Cookhla validation list `RA / T1D / Crohn's / GD` | `RA / T1D / Crohn's` (drop GD in Paper 2 context) |

## Brand expression PRESERVED (per directive)

| Phrase | Status | Reason |
|---|---|---|
| `autoimmune × thyroid 학자` (본인 자기소개) | preserved | Brand identity, not Paper 2 mechanism claim |
| `RA / T1D / Crohn's 같은 autoimmune disease 검증 history` | preserved (GD removed) | cookHLA tool's prior-validation domain |

## Files modified (this commit)

| File | Type | What changed |
|---|---|---|
| `project/manuscript_p2_brief/p2_advisor_discussion.html` | edit | 3 occurrences fixed (line 591, 658-paragraph, 1160); long discussion paragraph wrapped with deprecation banner |
| `project/manuscript_p2_brief/paper2_brief.html` | edit | 1 occurrence fixed (line 697 가설 framing) + audit notice inline |
| `project/results/p2_pillar1_forest/discussion_paragraph.md` | rewrite | DEPRECATED header + word substitutions; preserved for record only |
| `project/results/p2_pillar1_forest/PILLAR1_FOREST_SUMMARY.md` | edit | DEPRECATED header + one-line conclusion reframed; record only |
| `project/notebooks_or_scripts/v17_paper2_pillar1_forest.py` | edit | DEPRECATED docstring banner; figure title + 3 string-literal occurrences substituted (re-runs of v1 will not re-introduce conflated terms) |

## Files NOT modified (in scope but out of edit scope)

These contain conflated terminology but were left alone with rationale:

| File | Why not edited |
|---|---|
| `project/manuscript_v8/*` (Paper 1 main paper bib/sections) | Paper 1 voice-protected; user keyboard. Paper 1 is THCA driver-excluded paper, not Paper 2; "autoimmune-PTC" mentions there are intra-paper context, not Paper 2 conflation. |
| `project/reports/2026_*` (decision/audit history) | Historical record; rewriting would falsify the timeline. New reports going forward use the corrected terminology. |
| `project/notebooks_or_scripts/v17_*.py` (other analysis scripts) | Non-Paper-2 analyses (D5P6 BCR, P2 power planB, hla_analysis, etc.) — string literals there describe their own analytical context, not Paper 2 manuscript framing. |
| `project/scripts/v17_hla_autoimmune/run_autoimmune_hla.py` | Pre-Paper-2 era code; outputs are not in Paper 2 manuscript pipeline. |

## Deferred to Task A (Pillar I forest v2)

The deprecation banners on v1 outputs note that the **structural** problem (Korean
PTC vs Han Chinese GD as a `continuum`) is not just terminology — the analytical
framing must change:

- **v1**: Korean PTC pool (n=874) vs **Chu 2018 Han Chinese GD (n=1,468)** + ctrl (n=1,490)
- **v2 (Task A)**: Korean PTC pool (n=874) vs **Korean baseline allele frequency**
  - Source: Lee 2014 Tissue Antigens Korean reference, KOTRY donor cohort, or cookHLA Korean ref
  - Chu 2018 GD comparison: ONE Discussion line ("shared HLA background with GD,
    but disease distinct — see Paper 3")
  - Chu 2018 vs Chinese ctrl forest: Paper 3 reserve

Output dir for Task A: `project/results/p2_pillar1_forest_v2/`

## Cross-contamination footprint check (Paper 2 context only)

After this commit:

```bash
$ grep -rln "autoimmune-PTC\|autoimmune PTC" project/manuscript_p2_brief/ project/results/p2_pillar1_forest/
# (none — clean)

$ grep -rln "Pan-Asian autoimmune-thyroid" project/manuscript_p2_brief/ project/results/p2_pillar1_forest/
# (none — clean)

$ grep -rln "Graves\|GD\b" project/manuscript_p2_brief/p2_advisor_discussion.html
# Remaining matches: deprecation banner mentions ("see Paper 3"), prior-version
# preserved paragraph (color-greyed + banner-marked), Q12 future-work paragraph
# referencing Asian HLA literature where GD is a brand-context allele frequency
# data point. All match the directive's "shared genetic background — see Paper 3"
# allowance.
```

The remaining `Graves`/`GD` mentions in p2_advisor_discussion.html are inside the
explicit deprecation banner ("GD comparison reduced to one Discussion line") and
the v1-preserved paragraph (color-greyed, banner-marked). Net Paper 2 narrative
is HT-only; v1 paragraph is record-only.

---

*Generated 2026-05-04 by Task B of Paper 2 isolated session.*
