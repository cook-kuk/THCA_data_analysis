# HLA more-data registry — Paper 2 / Paper 4 only

**Date:** 2026-05-06  
**Status:** REGISTRY ONLY. No analysis execution. No forest computation. No manuscript prose.  
**Scope:** Paper 2 Korean baseline / AITD evidence pack and Paper 4 Pan-Asian Graves disease HLA backlog registry.  
**Frozen boundary:** Paper 1 is fixed as “A thyroid-lineage differentiation axis stratifies thyroid cancer beyond canonical driver mutations.” No Paper 1 edits or data expansion.

---

## 1. Executive verdict

**Paper 2:** Korean baseline frequency pack plus Korean AITD / Hashimoto-disease evidence can be strengthened now. The most useful baseline sources are `Lee 2005`, `In 2015`, `Chung 2010`, `Choe 2021`, `Jung 2023`, and `Baek 2023`. `Shin 2019` and `Cho 2011` are useful for Korean pediatric AITD / HD context, but their GD-specific rows must remain separated from Paper 2.

**Paper 4:** Pan-Asian GD HLA registry can be expanded now, but only as a backlog registry. `Chu 2018`, `Liao 2022`, `Shin 2019`, `Park 2005`, `Cho 1987`, Taiwanese / Japanese / Hong Kong studies, and Asian GD HLA reviews can be indexed for later extraction. No meta-analysis and no forest computation are allowed in this task.

**Immediate boundary:** This registry should support later source extraction decisions only. It does not compute ORs, forest plots, pooled effects, or new HLA calls.

---

## 2. Paper 2 Korean baseline registry

Target alleles:

```text
DPB1*05:01
A*02:07
B*46:01
DRB1*07:01
C*01:02
DQB1*02:01
```

| source_id | citation | year | population | sample_size | healthy/control status | loci covered | typing method | exact 4-digit data available? | DPB1*05:01 | A*02:07 | B*46:01 | DRB1*07:01 | C*01:02 | DQB1*02:01 | can fill Pillar I v2? | extraction priority | action |
|---|---|---:|---|---:|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P2_BASE_LEE2005 | Lee KW, Oh DH, Lee C, Yang SY. *Tissue Antigens* 65:437-447. doi:10.1111/j.1399-0039.2005.00386.x | 2005 | Korean | 485 | apparently unrelated healthy Koreans | A, B, C, DRB1, DQB1 | high-resolution DNA typing | yes, but full table needed | not covered | verify | verify | verify | verify | verify | yes for A/B/C/DRB1/DQB1, not DPB1 | high | use / verify table |
| P2_BASE_IN2015 | In JW et al. *Ann Lab Med* 35:429-435. doi:10.3343/alm.2015.35.4.429 | 2015 | Korean | 613 | healthy unrelated donors / cord-blood-unit reference | A, B, C, DRB1, DQB1 | sequence-based typing | yes, article table available | not covered | yes, extract | yes, extract | yes, extract | yes, extract | yes, extract | yes for 5/6, not DPB1 | top | use |
| P2_BASE_CHUNG2010 | Chung HY, Yoon JA, Han BY, Song EY, Park MH. *Korean J Lab Med* 30:685-696. doi:10.3343/kjlm.2010.30.6.685 | 2010 | Korean | 474 | healthy Koreans | A, B, C, DRB1 | high-resolution DNA typing | yes from abstract/metadata; full table needed | not covered | verify | verify | verify | verify | not covered | yes for A/B/C/DRB1 support | medium-high | use / verify |
| P2_BASE_CHOE2021 | Choe W et al. *Ann Lab Med* 41:310-317. doi:10.3343/alm.2021.41.3.310 | 2021 | Korean | 128 | healthy unrelated Korean adults | A, B, C, DRB1 | One Lambda AllType NGS, 8-digit resolution | yes, collapse 8-digit to 4-digit with note | not covered | verify / collapse | verify / collapse | verify / collapse | yes, `C*01:02:01` reported frequent | not covered | support for A/B/C/DRB1 only | medium | use as cross-check |
| P2_BASE_JUNG2023 | Jung K et al. *HLA* 101:602-612. doi:10.1111/tan.14980 | 2023 | Korean | 339 | unrelated healthy Korean subjects | A, B, C, DRB1, DRB3/4/5, DQB1, DQA1, DPB1, DPA1 | NGS, two 11-locus kits | yes, full article table needed | yes, `DPB1*05:01:01` 35.1% in abstract | verify | verify | verify | yes, `C*01:02:01` 18.44% in abstract | verify | yes for all six if full table accessed | top | use / verify table |
| P2_BASE_BAEK2023 | Baek IC et al. *HLA* 101:613-622. doi:10.1111/tan.14981 | 2023 | South Korean | not confirmed in accessible abstract | healthy South Korean donor/reference cohort; full text needed | A, B, C, DRB1/3/4/5, DQA1, DQB1, DPA1, DPB1 | amplicon-based NGS | likely yes, full table required | likely | likely | likely | likely | likely | likely | possible all-locus support | medium | verify before use |
| P2_BASE_JEKARL2021 | Jekarl DW et al. *HLA* 97:188-197. doi:10.1111/tan.14167 | 2021 | Korean | 1,293 selected from registry; 132-sample NGS performance panel | unrelated healthy donors / registry reference | A, B, C, DRB1 | NGS performance plus frequency/haplotype reference | yes if full table obtained | not covered | verify | verify | verify | verify | not covered | support for A/B/C/DRB1 only | medium | verify / optional |
| P2_BASE_SONG2002 | Song EY et al. *Tissue Antigens* 59:475-486. doi:10.1034/j.1399-0039.2002.590604.x | 2002 | Korean families | 107 families | family-based Korean reference | DRB1, DRB3/4/5, DQA1, DQB1, DPB1 | class II DNA typing / family haplotypes | yes if full table obtained | verify | not covered | not covered | verify | not covered | verify | class II background only | medium | verify / defer |

Notes:

- `In 2015` is the most immediate fill for A/B/C/DRB1/DQB1.
- `Jung 2023` is the strongest modern all-locus Korean candidate because it covers DPB1/DPA1/DQA1/DQB1 plus class I.
- `Baek 2023` is promising but should remain `verify` until sample size and table values are inspected.
- `Choe 2021` is high-resolution but lacks DQB1/DPB1, so it is a cross-check, not a full Pillar I source.

---

## 3. Paper 2 AITD / Hashimoto evidence registry

| source_id | citation | disease: HT / HD / GD / AITD / control | country / ancestry | n_case | n_control | loci | alleles relevant to 6-allele forest | effect direction | OR / p / frequency if reported | Paper 2 use |
|---|---|---|---|---:|---:|---|---|---|---|---|
| P2_AITD_SHIN2019 | Shin DH, Baek IC, Kim HJ, Choi EJ, Ahn M, Jung MH, et al. *PLoS ONE* 14:e0216941. doi:10.1371/journal.pone.0216941 | AITD split into GD and HD; controls | Korean children | 116 AITD total: 71 GD, 45 HD | 142 | A, B, C, DRB1, DQB1, DPB1 | HD: A*02:07 and DPB1*02:02; GD: B*46:01, C*01:02, DPB1*05:01; DPB1 amino-acid signatures | HD/AITD supports Korean thyroid-autoimmunity substrate; GD rows are Paper 4 only | GD: B*46:01 OR 3.96 Pc=0.008; C*01:02 OR 2.51 Pc=0.04; DPB1*05:01 OR 4.6 Pc=0.003. HD: A*02:07 OR 4.68 Pc=0.045; DPB1*02:02 OR 6.57 Pc=0.0001 | direct HD evidence + AITD background; GD extraction Paper 4 only |
| P2_AITD_CHO2011 | Cho WK, Jung MH, Choi EJ, Choi HB, Kim TG, Suh BK. *Horm Res Paediatr* 76:328-334. doi:10.1159/000331134 | AITD split into HD and GD; controls | Korean children | 73 AITD total: 32 HD, 41 GD | 159 | HLA class I and DRB1 | A*02, B*46, Cw*01, DRB1*08; DRB1*07/Cw*07 lower in GD; lower resolution than 4-digit | AITD / pediatric thyroid-autoimmunity background only | Accessible abstract reports AITD: higher A*02/B*46/Cw*01/DRB1*08 and lower A*30/B*07/Cw*07/DRB1*01 vs controls | AITD background; not quantitative Pillar I fill |
| P2_AITD_UEDA2014 | Ueda S et al. *J Clin Endocrinol Metab* 99:E379-E383. doi:10.1210/jc.2013-2841 | GD and HT; controls | Japanese | 991 AITD total: 547 GD, 444 HT | 481 | A, C, B, DRB1, DQB1, DPB1 | GD: B*46:01 and DPB1*05:01 among reported susceptible alleles; HT alleles separate | Non-Korean disease-separation support | Exact OR/p requires table extraction | AITD background only / defer |
| P2_AITD_HUH1986 | Huh KB et al. *Korean J Intern Med* 1:243-248. doi:10.3904/kjim.1986.1.2.243 | AITD / Graves-heavy historical cohort | Korean | table extraction needed | table extraction needed | serologic A, B, C, DR | no direct 4-digit mapping | historical only | old serologic associations; not compatible with 4-digit source table | drop from quantitative Paper 2; optional historical note |

Paper 2 rule:

- Use `Shin 2019` and `Cho 2011` only for **Korean AITD / HD background**.
- Do not move GD-specific rows into Paper 2 forest or Paper 2 main claim.
- Do not reactivate Chu 2018 GD inside Paper 2.

---

## 4. Paper 4 Pan-Asian GD registry

This table is a backlog-only registry. It is **not** a Paper 2 input table.

| source_id | citation | country | ancestry | n_GD | n_control | HLA loci | key alleles | amino-acid positions if reported | OR / p | summary-stat extractable? | backlog priority |
|---|---|---|---|---:|---:|---|---|---|---|---|---|
| P4_GD_CHU2018 | Chu X et al. *J Med Genet* 55:685-692. doi:10.1136/jmedgenet-2017-105146 | China | Han Chinese | 1,468 | 1,490 | HLA imputation + MHC amino-acid fine mapping | DPA1*02:02, DPB1*05:01, B*46:01, C*01:02, A*02:07, DQB1*02:01, DRB1*07:01 | HLA-DPβ1 position 205; HLA-DPα1 Met11; HLA-B positions 66/99; HLA-DRβ1 position 28 | Table 2: DPB1*05:01 OR 1.90, p=1.73e-26; B*46:01 OR 2.38, p=8.78e-21; DQB1*02:01 and DRB1*07:01 protective | yes | top |
| P4_GD_LIAO2022 | Liao WL et al. *Front Endocrinol* 13:842673. doi:10.3389/fendo.2022.842673 | Taiwan | Taiwanese | 2,998 | 29,083 EMR controls without thyroid disorder / thyroid medication / abnormal thyroid labs | HLA imputation, class I and II | A*02:07, B*46:01, C*01:02 and multi-locus genotype associations | not a primary amino-acid fine-mapping paper | genotype associations and comorbidity associations extractable; allele-frequency vs genotype distinction required | yes, with caution | high |
| P4_GD_SHIN2019 | Shin DH et al. *PLoS ONE* 14:e0216941. doi:10.1371/journal.pone.0216941 | Korea | Korean children | 71 | 142 | A, B, C, DRB1, DQB1, DPB1 | B*46:01, C*01:02, DPB1*05:01; HD comparator A*02:07 | DPB1 Leu35 and Glu55 strongly associated with GD | B*46:01 OR 3.96 Pc=0.008; C*01:02 OR 2.51 Pc=0.04; DPB1*05:01 OR 4.6 Pc=0.003; Leu35/Glu55 OR 23.38 P=0.0002 | yes | top Korean bridge |
| P4_GD_PARK2005 | Park MH et al. *Hum Immunol* 66:741-747. doi:10.1016/j.humimm.2005.03.001 | Korea | Korean | 198 | 200 | DRB1, DQB1 | DRB1*08:03, DRB1*16:02, DQB1 alleles; protective DRB1*07:01 noted in reviews | none noted | PubMed abstract reports DR/DQ associations; table extraction needed for full OR/p | yes | high |
| P4_GD_CHO1987 | Cho BY et al. *Tissue Antigens* 30:119-121. doi:10.1111/j.1399-0039.1987.tb01607.x | Korea | Korean | 128 | 220 | serologic A, B, C, DR | B13, DR5, DRw8 historical signals | none | relative risks reported; low-resolution only | partial | medium historical |
| P4_GD_JANG2011 | Jang HW et al. *Immunol Invest* 40:172-182. doi:10.3109/08820139.2010.525571 | Korea | Korean | 133 in review table | 200 in review table | DRB1 by PCR-SBT | DRB1*03:01, DRB1*08:02, DRB1*14:03 risk; DRB1*07:01, DRB1*13:02 protective per review | none expected | table extraction needed | likely | medium |
| P4_GD_CHEN2011 | Chen PL et al. *PLoS ONE* 6:e16635. doi:10.1371/journal.pone.0016635 | Taiwan | ethnic Chinese Han | 499 in review table | 504 in review table | A, B, C, DRB1, DQB1, DPB1 | B*46:01, DPB1*05:01, DQB1*03:02, DRB1*15:01, DRB1*16:02; DRB1*12:02 protective | not primary amino-acid paper | PubMed abstract reports B*46:01 OR 1.33, DPB1*05:01 OR 2.34, DRB1*16:02 OR 2.63 and corrected p-values | yes | high |
| P4_GD_UEDA2014 | Ueda S et al. *J Clin Endocrinol Metab* 99:E379-E383. doi:10.1210/jc.2013-2841 | Japan | Japanese | 547 | 481 | A, C, B, DRB1, DQB1, DPB1 | B*35:01, B*46:01, DRB1*14:03, DQB1*06:04, DPB1*05:01; protective A*24:02/A*33:03 and others in review | epistasis; exact residue positions not primary | OR/p table extractable | yes | high |
| P4_GD_ONUMA1994 | Onuma H et al. *Hum Immunol* 39:195-201. doi:10.1016/0198-8859(94)90260-7 | Japan | Japanese | needs extraction | needs extraction | DPB1 and B | DPB1*05:01, B46 | none | association extractable from paper | likely | medium |
| P4_GD_NAITO1987 | Naito S, Sasaki H, Arakawa K. *Endocrinol Jpn* 34:685-688. doi:10.1507/endocrj1954.34.5_685 | Japan | Japanese | 61 for A/B/C; 53 for DR/DQ | 1,998 | serologic A, B, C, DR, DQ | Bw46 | none | Bw46 23.0% vs 8.4%, Pc<0.003 | partial; old serology | medium historical |
| P4_GD_TSAI1989 | Tsai KS et al. *Taiwan Yi Xue Hui Za Zhi* 88:336-341. PMID:2794934 | Taiwan | Chinese in Taiwan | 93 | 106 | serologic DR/DQ | DR2, DR9, DQw1, DRw53 | none | DR2 and DR9 signals in abstract; old serology | partial | low-medium historical |
| P4_GD_CAVAN1994 | Cavan DA et al. *Clin Endocrinol* 40:63-66. doi:10.1111/j.1365-2265.1994.tb02444.x | Hong Kong | Hong Kong Chinese | needs extraction | needs extraction | A, B, DR, DQ | sex-specific HLA associations; exact allele list needs full text | none | extractable only after full text | likely partial | defer |
| P4_GD_HAWKINS1985 | Hawkins BR et al. *Clin Endocrinol* 23:245-252. doi:10.1111/j.1365-2265.1985.tb03267.x | Hong Kong | Hong Kong Chinese | needs extraction | needs extraction | serologic HLA antigens | historical Hong Kong GD HLA signals | none | table extraction needed | partial | defer |
| P4_GD_REVIEW2023 | Stasiak M et al. *Front Immunol* 14:1256922. doi:10.3389/fimmu.2023.1256922 | mixed | Asian and Caucasian review | review | review | broad HLA | Asian GD: B*46:01, DPB1*05:01, DRB1*08:02/03, DRB1*16:02, DRB1*14:03, DRB1*04:05, DQB1*05:02, DQB1*03:03; protective DRB1*07:01 etc. | review summarizes amino-acid evidence including Shin 2019 | source-discovery table only | yes | high registry map |
| P4_GD_B46_META2013 | Li Y et al. *Int J Med Sci* 10:164-170. doi:10.7150/ijms.5158 | mixed Asian | Asian meta-analysis | 1,743 | 5,689 | HLA-B | B*46 | none | pooled B*46 effect and heterogeneity | yes, review/meta context only | medium |
| P4_GD_DRB1_META2024 | Li W et al. *Horm Metab Res* 56:859-868. doi:10.1055/a-2298-4366 | mixed Asian | Asian meta-analysis | meta-analysis | meta-analysis | DRB1 | DRB1 alleles | none | pooled DRB1 summary | yes, review/meta context only | medium |

---

## 5. Allele harmonization rules

1. **Use 4-digit format.** Normalize all target alleles to two-field names: `DPB1*05:01`, `A*02:07`, `B*46:01`, `DRB1*07:01`, `C*01:02`, `DQB1*02:01`.

2. **Collapse 6/8-digit calls only with a note.** Example: `DPB1*05:01:01` -> `DPB1*05:01`; `C*01:02:01` -> `C*01:02`; `DRB1*07:01:01:01` -> `DRB1*07:01`.

3. **Carrier frequency is not allele frequency.** Carrier/sample frequency uses people as denominator; allele frequency usually uses `2n` chromosomes. Do not mix them without an explicit conversion column.

4. **Genotype frequency is separate.** Liao 2022 reports many genotype-level combinations. These cannot be treated as allele frequencies.

5. **Control ancestry label is mandatory.** Keep `Korean`, `South Korean`, `Han Chinese`, `Taiwanese`, `Japanese`, `Hong Kong Chinese`, `Southern Chinese`, and any proxy controls separate.

6. **Disease labels are strict.** Use only: `HT`, `HD`, `GD`, `AITD`, `control`, `PTC`. Mixed AITD papers must be split by disease before any future extraction.

7. **Old serology is historical only.** Terms like `Bw46`, `B13`, `DR5`, `DRw8`, `DQw1`, and `Cw*01` are not direct 4-digit alleles.

---

## 6. Immediate Paper 2 action candidates

| candidate | verdict | rationale |
|---|---|---|
| Use Lee/In/Chung/Choe/Jung/Baek as Korean baseline pack | go for extraction planning only | gives legacy high-resolution + modern NGS coverage; no OR or forest here |
| Start with In 2015 for A/B/C/DRB1/DQB1 | go | n=613, healthy unrelated Korean reference, exact 4-digit table expected |
| Use Jung 2023 for DPB1*05:01 and all-locus cross-check | go after table access | n=339, 11-locus NGS, abstract confirms DPB1*05:01:01 and C*01:02:01 frequencies |
| Use Shin 2019 and Cho 2011 as AITD/HD background | go with strict separation | HD/AITD context strengthens Paper 2; GD rows remain Paper 4 |
| Decide whether C*01:02 / DQB1*02:01 can be filled | go for source decision | In 2015 likely fills both; Jung/Baek can cross-check |
| Move DRB1*07:01 Harbin fallback to sensitivity if Korean sources fill it | go | Korean baseline references should replace proxy ancestry |

---

## 7. Paper 4 future action candidates

| candidate | verdict |
|---|---|
| Chu 2018 + Liao 2022 + Shin 2019 + Park 2005 registry | keep as Paper 4 backlog |
| Add Chen 2011 Taiwan and Ueda 2014 Japan | keep as high-priority backlog |
| Add Cho 1987 / Naito 1987 / Tsai 1989 / Cavan 1994 / Hawkins 1985 | keep as historical source registry |
| Use Stasiak 2023 systematic review as source-discovery map | useful, no computation |
| Run Paper 4 meta-analysis or forest | no-go until Paper 4 gates clear |

---

## 8. Risk notes

| risk | note |
|---|---|
| allele frequency vs carrier frequency | Do not mix `2n` allele frequency with carrier/sample frequency. |
| disease mismatch HT vs GD | Paper 2 can use HT/HD/AITD context; GD evidence is Paper 4 unless explicitly labeled as background only. |
| ancestry mismatch | Korean, Han Chinese, Taiwanese, Japanese, Hong Kong Chinese, and Southern Chinese references are not interchangeable. |
| old low-resolution HLA typing | Serology and antigen-level studies are historical context, not 4-digit source tables. |
| sample-size imbalance | Large registry references and small NGS cohorts should not be pooled without a later pre-specified method. |
| table access | Several Wiley/Karger sources need full table verification before extraction. |
| imputation vs direct typing | Chu/Liao use HLA imputation; Korean baseline sources use direct typing/NGS/SBT. Keep method labels explicit. |

---

## 9. Next command options

Allowed next commands:

```text
extract Paper2 Korean baseline tables only
```

Extract published Korean baseline table rows only. Do not compute OR. Do not run forest.

```text
prepare 6-allele forest source table
```

Create source/provenance input table only. Do not compute forest statistics.

```text
freeze Paper4 GD registry
```

Keep all GD rows as Paper 4 backlog only.

```text
stop and return to Paper1 Hook
```

Leave HLA frozen and return to Paper 1 Hook work.

Blocked commands:

```text
run meta-analysis
run forest
reactivate Chu 2018 inside Paper 2
run Bundang Graves
run cookHLA SNP integration
touch Paper 3
edit manuscript
```

---

## 10. Source anchors checked

- Lee 2005: PubMed PMID 15853898 / DOI 10.1111/j.1399-0039.2005.00386.x
- In 2015: PMC / Ann Lab Med, 613 healthy unrelated donors, SBT, A/B/C/DRB1/DQB1
- Chung 2010: SNU Pure / DOI 10.3343/kjlm.2010.30.6.685, 474 healthy Koreans
- Choe 2021: Ann Lab Med 41:310-317, PMID 33303716
- Jung 2023: PubMed PMID 36719349 / DOI 10.1111/tan.14980, 339 unrelated healthy subjects, 11 loci
- Baek 2023: PubMed PMID 36720674 / DOI 10.1111/tan.14981, 11 loci amplicon-based NGS
- Shin 2019: PLoS ONE / DOI 10.1371/journal.pone.0216941, Korean children AITD: 71 GD, 45 HD, 142 controls
- Cho 2011: Karger / DOI 10.1159/000331134, 73 AITD children and 159 controls
- Chu 2018: J Med Genet / DOI 10.1136/jmedgenet-2017-105146, Han Chinese GD n=1,468 and controls n=1,490
- Liao 2022: Front Endocrinol / DOI 10.3389/fendo.2022.842673, Taiwanese GD EMR/imputation
- Park 2005: PubMed PMID 15993720 / DOI 10.1016/j.humimm.2005.03.001
- Cho 1987: SNU Pure / DOI 10.1111/j.1399-0039.1987.tb01607.x
- Chen 2011: PubMed PMID 21307958 / DOI 10.1371/journal.pone.0016635
- Stasiak 2023 systematic review: Front Immunol / DOI 10.3389/fimmu.2023.1256922

---

Registry complete. No analysis execution. No forest computation. No manuscript edit.
