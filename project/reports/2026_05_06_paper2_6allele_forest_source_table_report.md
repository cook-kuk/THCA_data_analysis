---
title: "Paper 2 6-allele forest source table report"
date: 2026-05-06
status: "SOURCE TABLE ONLY - no OR, no Fisher p, no forest, no meta-analysis, no manuscript prose"
output:
  - "project/results/p2_pillar1_forest_v2/paper2_6allele_forest_source_table.tsv"
inputs:
  - "project/results/p2_pillar1_forest_v2/paper2_korean_baseline_source_table.tsv"
  - "project/reports/2026_05_06_paper2_korean_baseline_extraction_report.md"
  - "project/results/p2_pillar1_forest/korean_PTC_pool_per_subcohort.tsv"
  - "project/results/p2_pillar1_forest_v2/korean_baseline_allele_freq.tsv"
---

# 1. Executive verdict

The 6-allele source/provenance table is prepared, but it is not ready for forest execution without Yu approval because the Korean PTC values are carrier frequencies while the selected Korean baseline values are published allele frequencies.

All six PTC values are available from the existing Korean PTC pool file:
`A*02:07`, `B*46:01`, `C*01:02`, `DPB1*05:01`, `DQB1*02:01`, and `DRB1*07:01`.

Primary baseline choices are source-resolved:
`In 2015` for `A*02:07`, `B*46:01`, `C*01:02`, `DQB1*02:01`, and `DRB1*07:01`; current v2 `AFND South Korea` for `DPB1*05:01`, with `Jung 2023` as a direct Korean 11-locus cross-check.

# 2. Six-allele readiness table

| allele | PTC value available? | baseline primary | baseline value | metric compatibility | ready_for_forest | recommended_use |
|---|---|---|---:|---|---|---|
| A*02:07 | yes, 0.0812 carrier frequency | In 2015 | 3.51% AF | PTC carrier vs baseline allele frequency | needs_yu_approval | primary_6allele |
| B*46:01 | yes, 0.103 carrier frequency | In 2015 | 5.06% AF | PTC carrier vs baseline allele frequency | needs_yu_approval | primary_6allele |
| C*01:02 | yes, 0.2414 carrier frequency | In 2015 | 17.81% AF | PTC carrier vs baseline allele frequency | needs_yu_approval | primary_6allele |
| DPB1*05:01 | yes, 0.532 carrier frequency | current v2 AFND South Korea | 0.3667 AF | PTC carrier vs baseline allele frequency | needs_yu_approval | primary_6allele |
| DQB1*02:01 | yes, 0.0 carrier frequency | In 2015 | 2.12% AF | PTC carrier vs baseline allele frequency | needs_yu_approval | primary_6allele |
| DRB1*07:01 | yes, 0.1144 carrier frequency | In 2015 | 6.93% AF | PTC carrier vs baseline allele frequency | needs_yu_approval | primary_6allele |

# 3. Metric compatibility check

| side | metric | source |
|---|---|---|
| Korean PTC | carrier frequency | existing Korean PTC pool Combined rows |
| Korean baseline, In 2015 | allele frequency | published Table 1 over `2n_alleles` |
| Korean baseline, AFND South Korea DPB1 | allele frequency | existing v2 baseline table |
| Korean baseline, Jung 2023 cross-check | allele frequency | abstract-level value over `2n_alleles` |

No conversion was performed. Because metric types differ, every allele is marked `needs_yu_approval`.

# 4. Recommended primary baseline per allele

| allele | recommended primary baseline | reason |
|---|---|---|
| A*02:07 | In 2015 | direct Korean healthy SBT Table 1; Choe/Chung agree closely |
| B*46:01 | In 2015 | direct Korean healthy SBT Table 1; Choe/Chung agree closely |
| C*01:02 | In 2015 | direct Korean healthy SBT Table 1; Jung/Choe/Chung support similar range |
| DPB1*05:01 | current v2 AFND South Korea | already in current v2 baseline table; Jung 2023 is a close direct Korean cross-check but full table is not archived |
| DQB1*02:01 | In 2015 | direct Korean healthy SBT Table 1; current AFND v2 has no Korean entry |
| DRB1*07:01 | In 2015 | direct Korean healthy SBT Table 1; replaces Harbin proxy as primary candidate |

# 5. Sensitivity source per allele

| allele | sensitivity / cross-check sources |
|---|---|
| A*02:07 | Choe 2021, Chung 2010, AFND South Korea current v2 |
| B*46:01 | Choe 2021, Chung 2010, AFND South Korea current v2 |
| C*01:02 | Jung 2023, Choe 2021, Chung 2010 |
| DPB1*05:01 | Jung 2023; Baek 2023/Song 2002 only after full table access |
| DQB1*02:01 | Lee 2005/Jung 2023/Baek 2023/Song 2002 only after full table access |
| DRB1*07:01 | Choe 2021, Chung 2010, Harbin Korean proxy as sensitivity only |

# 6. Items requiring Yu approval

1. Whether carrier-frequency PTC values can be compared to published allele-frequency Korean baselines in the planned forest, or whether a conversion/spec change is required.

2. Whether `In 2015` should be the primary source for `C*01:02`, `DQB1*02:01`, and `DRB1*07:01`.

3. Whether current v2 `AFND South Korea` remains primary for `DPB1*05:01`, with `Jung 2023` as cross-check, or whether `Jung 2023` should replace it after full table archiving.

4. Whether the Harbin Korean `DRB1*07:01` fallback is formally demoted to sensitivity only.

5. Whether `DQB1*02:01` with Korean PTC carrier frequency `0.0` can enter a later forest, and under what pre-specified continuity rule. No such rule was applied here.

# 7. Whether full forest can be run next

Not yet. The source/provenance table is complete enough for advisor review, but the full forest should not run until Yu approves metric compatibility and source hierarchy.

Current status: `ask Yu approval` is the clean next step. `prepare A1 forest script only` is possible after that if the script only consumes the source table and does not execute.

# 8. Next command options

```text
ask Yu approval
```

Recommended next step.

```text
prepare A1 forest script only
```

Allowed only as script preparation. Do not execute.

```text
run full A1 forest after approval
```

Blocked until Yu/source/metric approval.

```text
freeze and return to Paper 1 Hook
```

Stop HLA work here.

# Source anchors

- Korean PTC Combined values: `project/results/p2_pillar1_forest/korean_PTC_pool_per_subcohort.tsv`
- Baseline extraction: `project/results/p2_pillar1_forest_v2/paper2_korean_baseline_source_table.tsv`
- Current v2 baseline: `project/results/p2_pillar1_forest_v2/korean_baseline_allele_freq.tsv`
