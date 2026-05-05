---
title: "Paper 2 Korean baseline HLA table extraction report"
date: 2026-05-06
status: "SOURCE/PROVENANCE EXTRACTION ONLY - no OR, no Fisher p, no forest, no meta-analysis, no manuscript prose/edit"
outputs:
  - "project/results/p2_pillar1_forest_v2/paper2_korean_baseline_source_table.tsv"
---

# 1. Executive verdict

`C*01:02` and `DQB1*02:01` can be filled from `In 2015` now. `In 2015` reports both as allele frequencies in healthy Korean controls: `C*01:02 = 17.81%` and `DQB1*02:01 = 2.12%`.

`DPB1*05:01` can be cross-checked from `Jung 2023` now at the abstract level: `DPB1*05:01:01 = 35.1%`, collapsed to `DPB1*05:01`. `Song 2002` and `Baek 2023` may provide additional DPB1/DQB1 support, but the target table values were not directly accessible in this pass.

The Harbin Korean fallback can be demoted to sensitivity for `DRB1*07:01`: Korean sources now provide direct Korean baseline values for DRB1*07:01 (`In 2015`, `Choe 2021`, `Chung 2010`). No forest or OR was computed.

# 2. Source-by-source extraction summary

| source | status | extracted target values | unusable / inaccessible items |
|---|---|---|---|
| `In 2015` | strongest immediate primary candidate | A*02:07, B*46:01, C*01:02, DRB1*07:01, DQB1*02:01 | DPB1 not covered |
| `Lee 2005` | citation verified; table values not directly visible | none in this pass | A/B/C/DRB1/DQB1 target values marked `table_inaccessible`; DPB1 not covered |
| `Jung 2023` | modern 11-locus candidate; abstract values extracted | C*01:02, DPB1*05:01 | A*02:07, B*46:01, DRB1*07:01, DQB1*02:01 require full table |
| `Choe 2021` | exact PDF Table 2 extraction | A*02:07, B*46:01, C*01:02, DRB1*07:01 subtypes | DQB1 and DPB1 not covered |
| `Chung 2010` | indexed PDF snippet extraction | A*02:07, B*46:01, C*01:02, DRB1*07:01 | DQB1 and DPB1 not covered; archive full table before numeric reuse |
| `Baek 2023` | citation/landing metadata only | none in this pass | all target values marked `table_inaccessible` |
| `Song 2002` | class II family source | none in this pass | DPB1/DQB1/DRB1 table values marked `table_inaccessible`; A/B/C not covered |

# 3. Target-allele coverage matrix

Legend: `value` = extracted allele frequency; `TI` = table inaccessible; `NC` = locus not covered.

| allele | In 2015 | Lee 2005 | Jung 2023 | Choe 2021 | Chung 2010 | Baek 2023 | Song 2002 |
|---|---:|---:|---:|---:|---:|---:|---:|
| DPB1*05:01 | NC | NC | 35.1% | NC | NC | TI | TI |
| A*02:07 | 3.51% | TI | TI | 3.5% | 3.5% | TI | NC |
| B*46:01 | 5.06% | TI | TI | 5.5% | 5.1% | TI | NC |
| DRB1*07:01 | 6.93% | TI | TI | 5.1% + 0.8% as separate 8-digit subtype rows | 7.3% | TI | TI |
| C*01:02 | 17.81% | TI | 18.44% | 19.9% | 17.4% | TI | NC |
| DQB1*02:01 | 2.12% | TI | TI | NC | NC | TI | TI |

# 4. Recommended primary baseline source per allele

| allele | recommended primary candidate | rationale |
|---|---|---|
| DPB1*05:01 | Jung 2023, pending full-table archiving | only directly visible Korean healthy DPB1 value in this pass; 11-locus NGS, n=339 |
| A*02:07 | In 2015 | direct healthy Korean SBT Table 1, n=613; Choe/Chung agree closely |
| B*46:01 | In 2015 | direct healthy Korean SBT Table 1, n=613; Choe/Chung agree closely |
| DRB1*07:01 | In 2015 | direct healthy Korean SBT Table 1, n=613; enough to move Harbin proxy out of primary |
| C*01:02 | In 2015 | direct healthy Korean SBT Table 1, n=613; Jung/Choe/Chung all support similar range |
| DQB1*02:01 | In 2015 | only directly extracted Korean healthy DQB1*02:01 value in this pass |

# 5. Recommended sensitivity sources

| source | use |
|---|---|
| Choe 2021 | A/B/C/DRB1 high-resolution NGS cross-check; keep DRB1*07:01 subtypes separate unless a collapse rule is explicitly approved |
| Chung 2010 | A/B/C/DRB1 cross-check after full table is archived |
| Jung 2023 | all-locus modern NGS cross-check once full table is accessible |
| Lee 2005 | legacy 5-locus cross-check if full table is retrieved |
| Song 2002 | class II family/haplotype sensitivity only, not primary allele-frequency baseline |
| Baek 2023 | defer until full table/sample size are accessible |

# 6. Data caveats

Allele frequency vs carrier frequency: all extracted numeric values in the TSV are recorded as published allele frequencies. They were not converted to carrier frequencies.

Population/control comparability: sources are Korean healthy/control references, but donor/cord blood/family-based sampling differs by source.

Typing resolution: Choe 2021 reports 6/8-digit values and Jung 2023 abstract reports 3-field values. These are collapsed to 4-digit in `allele_4digit`, with the original allele stored in `collapsed_from`. Choe DRB1*07:01 appears as two subtype rows to avoid hidden summation.

Missing DPB1: older Korean references often omit DPB1. Lee 2005, In 2015, Choe 2021, and Chung 2010 cannot fill DPB1*05:01.

Table access: Lee 2005, Baek 2023, and Song 2002 need full-table access before values can be used. No values were inferred.

# 7. Next command options

```text
prepare 6-allele forest source table
```

Use the TSV as provenance input schema only. Do not compute a forest.

```text
ask Yu approval
```

Confirm whether `In 2015` can serve as the primary fill for `C*01:02`, `DQB1*02:01`, and `DRB1*07:01`, and whether `Jung 2023` is acceptable for DPB1 cross-check after full table archiving.

```text
freeze and return to Paper 1 Hook
```

Stop HLA extraction work here.

# Source anchors

- In 2015 PMC/PDF: `https://pmc.ncbi.nlm.nih.gov/articles/PMC4446582/`
- Choe 2021 PDF: `https://synapse.koreamed.org/upload/synapsexml/3039alm/pdf/alm-41-310.pdf`
- Jung 2023 PubMed/search metadata: `https://pubmed.ncbi.nlm.nih.gov/36719349/`
- Lee 2005 PubMed: `https://pubmed.ncbi.nlm.nih.gov/15853898/`
- Chung 2010 KCI/indexed PDF metadata: `https://www.kci.go.kr/kciportal/landing/article.kci?arti_id=ART001497939`
- Baek 2023 Ovid metadata: `https://www.ovid.com/journals/hlan/fulltext/10.1111/tan.14981~distributions-of-11loci-hla-alleles-typed-by-ampliconbased`
- Song 2002 metadata: `https://snu.elsevierpure.com/en/publications/hla-class-ii-allele-and-haplotype-frequencies-in-koreans-based-on`
