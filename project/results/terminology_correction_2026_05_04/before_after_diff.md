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

---

## Second sweep (2026-05-04, cross-portfolio review)

추가 발견된 occurrences — 이전 sweep 누락된 portfolio-level 자료:

| File | Line | Before | After |
|---|---|---|---|
| `project/submission/papers_overview.html` | 264 (HTML comment) | `<!-- Paper 2 — Autoimmune-overlap PTC -->` | `<!-- Paper 2 — Hashimoto-overlap PTC -->` |
| `project/submission/papers_overview.html` | 370 (Table 1 row) | `Autoimmune-overlap PTC ★` | `Hashimoto-overlap PTC ★` |
| `project/manuscript_p2_brief/p2_advisor_discussion.html` | 6 (`<title>` tag) | `Paper 2 — Autoimmune-overlap PTC \| Advisor Discussion Brief` | `Paper 2 — Hashimoto-overlap PTC \| Advisor Discussion Brief` |

## Files NOT modified — second sweep additions

| File | Why not edited (this sweep) |
|---|---|
| `project/manuscript_p3_brief/*.html` | **Paper 3 territory** (GD/autoimmune-thyroid). SCOPE 명시: "Paper 3 file touch 금지 (source file read 도 금지)". Paper 3 brief 의 occurrences 는 Paper 3 own scope. |
| `project/three_papers_index.html` | Cross-paper meta-index (3 papers all-context). "autoimmune-overlap PTC" 표현이 portfolio-level branding (advisor 지시 미포함 영역). Future revision 권장이나 본 sweep scope 제외. |
| `project/STATUS_PAPER2_2026_05_02.md`, `STATUS_PAPER3_2026_05_02.md` | Status meta-doc; "autoimmune-PTC" mentions = substitution-rule documentation context (recording the change). Active framing 아님. |
| `project/notebooks_or_scripts/v17_F5_autoimmune_PTC_mechanism.py` | Pre-isolation script (figure title "Autoimmune-PTC mechanism"). 본 figure 가 Paper 1 의 Pillar 5 figure 라서 Paper 1 voice-protected 영역. Out of Task B scope. |
| `project/scripts/v17_hla_autoimmune/run_autoimmune_hla.py:665` | Pre-isolation era; bib reference quote inside script. Historical record. |
| `project/results/audit_2026_04_30/FINAL_COMPREHENSIVE_SUMMARY_v{2,3}.md` | Historical audit summaries — frozen artifact. |

## Updated cross-contamination footprint check (after second sweep)

```bash
# Paper 2 active narrative
$ grep -rln "Autoimmune-overlap PTC\|autoimmune-PTC\|autoimmune PTC" project/manuscript_p2_brief/
# Only in deprecation banner (p2_advisor_discussion.html line 671) — recording
# the substitution context itself; not an active framing claim.

$ grep -rln "Pan-Asian autoimmune-thyroid" project/manuscript_p2_brief/
# (none — clean ✓)

$ grep -rln "Autoimmune-overlap PTC\|autoimmune-overlap PTC" project/submission/papers_overview.html
# (none — clean ✓ after second sweep)
```

Net Paper 2 active narrative (deliverables that go to advisor):
- `paper2_brief.html` (story version, 38 pages PDF) — clean ✓
- `p2_advisor_discussion.html` (extended, 50 pages PDF) — clean except 1 deprecation banner (intentional record) ✓
- `p2_pillar1_forest/` v1 outputs — DEPRECATED-banner wrapped (record only) ✓
- `p2_pillar1_forest_v2/` (Task A output) — born clean, Korean-baseline framed ✓
- `papers_overview.html` (master index) — clean after second sweep ✓

---

*Updated 2026-05-04 by Task B second sweep of Paper 2 isolated session.*
