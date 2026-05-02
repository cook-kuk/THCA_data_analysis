# Paper 2 Pillar I Forest v2 — Korean PTC vs Korean baseline (HT context)

**Date:** 2026-05-02
**Scope:** Paper 2 (Hashimoto-overlap PTC) — GD comparison reserved for Paper 3.

## ★ One-line conclusion (HT-only)

Korean PTC pool (n=874) carries DPB1*05:01 at
**53.2%** — significantly elevated above Korean baseline
**36.7%** (OR 1.96, p=1.1e-10);
3/4 focal alleles direction-consistent with
HT-context expectation.

## Forest table

| Allele | Direction | PTC freq | Baseline freq | OR (95% CI) | p |
|---|---|---|---|---|---|
| A*02:07 | risk | 8.1% (71/874) | 3.4% (n=4,613, 2 pops) | 2.54 (1.90–3.40) | 3.1e-10 |
| B*46:01 | risk | 10.3% (90/874) | 4.7% (n=4,613, 2 pops) | 2.33 (1.80–3.01) | 1.2e-10 |
| DPB1*05:01 | risk | 53.2% (465/874) | 36.7% (n=680, 3 pops) | 1.96 (1.60–2.41) | 1.1e-10 |
| DRB1*07:01 | protective | 11.4% (100/874) | 5.0% (n=201, 1 pops) | 2.46 (1.26–4.79) | 0.0084 |

## Methods

- Korean PTC pool: K2 (n=235) + Lee 2024 (n=630) + GSE286332-PTC (n=9) =
  n=874, from per-subcohort `Combined` row in v1
  `korean_PTC_pool_per_subcohort.tsv`.
- Korean baseline: sample-weighted pooled carrier frequency across
  AFND-indexed South-Korea-resident populations. For alleles with no
  AFND South Korea entry, fall back to China Harbin Korean (Korean ethnic
  minority resident in China; per-row `source` column flags this proxy).
  For alleles AFND has neither: skipped from forest, flagged in
  `korean_baseline_allele_freq.tsv` as `unavailable` — Lee 2014 Tissue
  Antigens Korean reference paper recommended for full-panel coverage.
- Effect size: 2×2 OR with continuity correction for zero cells; 95% CI
  via log-OR ± 1.96 SE.
- Direction prior: from Chu 2018 Han Chinese GD vs ctrl labels — used only
  for plot color coding, NOT as the comparator.

## Caveats

- **DRB1*07:01 direction discord**: prior label "protective" (from Chu 2018
  GD context, where ctrl 15.3% → GD 7.1%) does not generalize cleanly to
  Korean PTC vs Korean baseline (PTC 11.4% > Harbin Korean 5.0%, OR 2.46).
  Two alternative explanations: (i) the small Harbin Korean baseline (n=201,
  Korean diaspora) under-represents the true Korean DRB1*07:01 frequency;
  (ii) the GD-context "protective" effect does not transfer to PTC. Both
  hypotheses point to needing a published Korean reference (Lee 2014) for
  high-confidence direction calls; this v2 forest reports the AFND-driven
  result honestly with the proxy caveat.
- **Coverage**: 4 of 6 focal alleles have an AFND Korean baseline. C*01:02
  and DQB1*02:01 are skipped — Lee 2014 reference recommended.
- **Population stratification within Korea**: the AFND South Korea pops
  are heterogeneous studies (different ascertainment cohorts, donor
  registries vs disease cohorts). The pooled baseline is a reasonable
  population-average estimate but does not control for sub-population
  structure within Korea.

## Why v2 supersedes v1

| | v1 (deprecated) | v2 (HT-only) |
|---|---|---|
| Comparator | Chu 2018 Han Chinese ctrl + GD | Korean baseline (AFND South Korea pool) |
| Disease scope | Conflated HT/GD continuum | Hashimoto-overlap PTC only |
| GD comparison | Pillar I main claim | One Discussion line (Paper 3 cross-ref) |
| Population strat | Korean vs Han Chinese (different) | Korean vs Korean (same population) |
| Cell Press positioning | Pan-Asian autoimmune-thyroid | Pan-Asian thyroid HLA susceptibility (HT context) |

## Files

```
project/results/p2_pillar1_forest_v2/
├── korean_baseline_allele_freq.tsv       (AFND South Korea pool, sample-weighted)
├── korean_PTC_vs_korean_baseline_forest.json
├── forest_paper2_HT_only.{pdf,png}
├── discussion_paragraph_v2.md            (Cell Press paste-ready)
└── PILLAR1_FOREST_V2_SUMMARY.md (this file)
```

## Paper 3 reserve (NOT this session)

Chu 2018 Han Chinese GD (n=1,468) vs ctrl (n=1,490) forest material —
already in `project/results/p2_pillar1_forest/` v1 — should move to
`project/results/p3_graves_pillar1_forest/` for the Paper 3 brief.
That migration is Paper 3 territory, NOT this session.

---

*Generated 2026-05-02 15:43 by v17_paper2_pillar1_forest_v2.py.*
