---
title: "HLA more-data evidence registry"
date: 2026-05-04
status: "REGISTRY ONLY - no analysis execution"
scope:
  paper2: "Korean baseline + Korean AITD/Hashimoto evidence pack"
  paper4: "Pan-Asian Graves disease HLA backlog registry"
hard_stops:
  - "No cookHLA SNP integration"
  - "No new raw genotype data download"
  - "No Bundang Graves work"
  - "No Chu 2018 GD forest reactivation inside Paper 2"
  - "No Paper 4 meta-analysis execution"
  - "No Paper 3 touch"
  - "No HLA LOH / neoantigen"
  - "No voice-protected prose"
  - "No manuscript edit"
  - "No large analysis execution"
allowed:
  - "Read-only PubMed / PMC / journal landing page lookup"
  - "Published table / allele frequency / OR / p-value extraction plan"
  - "Citation metadata confirmation"
  - "Registry TSV/MD writing"
---

# 1. Executive verdict

Paper 2 can be strengthened now with a Korean baseline frequency pack plus Korean AITD/HD evidence. The cleanest baseline sequence is `In 2015 + Lee 2005 + Jung 2023 + Jekarl 2021 + Chung 2010`, with `Shin 2019` and `Cho 2011` reserved for Korean AITD/HD disease context.

Paper 4 should remain a Pan-Asian GD registry only. Chu 2018, Liao 2022, Shin 2019, Park 2005, Cho 1987, and accessible Japanese/Hong Kong/Taiwanese GD sources can be indexed for future extraction, but no meta-analysis, forest, Bundang integration, or Paper 2 reactivation is allowed here.

# 2. Paper 2 Korean baseline registry

Target alleles: `DPB1*05:01`, `A*02:07`, `B*46:01`, `DRB1*07:01`, `C*01:02`, `DQB1*02:01`.

| source_id | citation | year | population | sample_size | healthy/control status | loci covered | typing method | exact 4-digit data available? | relevant alleles | can fill Pillar I v2? | extraction priority | action |
|---|---|---:|---|---:|---|---|---|---|---|---|---|---|
| P2_BASE_LEE2005 | Lee KW, Oh DH, Lee C, Yang SY. Allelic and haplotypic diversity of HLA-A, -B, -C, -DRB1, and -DQB1 genes in the Korean population. Tissue Antigens 65:437-447. doi:10.1111/j.1399-0039.2005.00386.x | 2005 | Korean | 485 | apparently unrelated healthy individuals | A, B, C, DRB1, DQB1 | high-resolution DNA typing | yes, but older notation must be normalized | A*02:07 likely; B*46:01 yes; DRB1*07:01 yes; C*01:02 yes; DQB1*02:01/02:02 may need ambiguity check; DPB1 absent | yes for 5/6; not DPB1 | high | use / verify table |
| P2_BASE_IN2015 | In JW et al. Allele and haplotype frequencies of HLA-A, -B, -C, -DRB1, and -DQB1 from sequence-based DNA typing data in Koreans. Ann Lab Med 35:429-435. doi:10.3343/alm.2015.35.4.429 | 2015 | Korean | 613 | healthy unrelated donors and cord blood units | A, B, C, DRB1, DQB1 | SBT | yes, Table 1 | A*02:07 yes; B*46:01 yes; DRB1*07:01 yes; C*01:02 yes; DQB1*02:01 yes; DPB1 absent | yes for 5/6; strongest immediate C/DQB1 fill | top | use |
| P2_BASE_JUNG2023 | Jung K, Kim JG, Shin S, Roh EY, Hong YJ, Song EY. Allele and haplotype frequencies of 11 HLA loci in Koreans by next-generation sequencing. HLA 101:602-612. doi:10.1111/tan.14980 | 2023 | Korean | 339 | unrelated healthy subjects | A, B, C, DRB1, DRB3/4/5, DQB1, DQA1, DPB1, DPA1 | NGS, two kits: NGSgo-MX11-3 and AllType NGS 11 loci | yes if full tables obtained | all six expected, including DPB1*05:01 and DQB1*02:01 | yes, best 11-locus modern baseline candidate | top | use / verify full table |
| P2_BASE_JEKARL2021 | Jekarl DW et al. HLA-A, -B, -C, -DRB1 allele and haplotype frequencies of the Korean population and performance characteristics of HLA typing by next-generation sequencing. HLA 97:188-197. doi:10.1111/tan.14167 | 2021 | Korean | 1,293 | unrelated healthy donors | A, B, C, DRB1 | NGS performance / allele frequency from registry data | yes if full table obtained | A*02:07, B*46:01, C*01:02, DRB1*07:01; DPB1/DQB1 absent | yes for class I + DRB1 cross-check; cannot fill DQB1/DPB1 | high for A/B/C/DRB1 | use / verify |
| P2_BASE_CHUNG2010 | Chung HY, Yoon JA, Han BY, Song EY, Park MH. Allelic and haplotypic diversity of HLA-A, -B, -C, and -DRB1 genes in Koreans defined by high-resolution DNA typing. Korean J Lab Med 30:685-696. doi:10.3343/kjlm.2010.30.6.685 | 2010 | Korean | 474 | healthy Koreans | A, B, C, DRB1 | high-resolution DNA typing | yes if KCI/full table obtained | A*02:07, B*46:01, C*01:02, DRB1*07:01; DQB1/DPB1 absent | yes for A/B/C/DRB1 support only | medium-high | use / verify |
| P2_BASE_BAEK2021 | Baek IC et al. Allele and haplotype frequencies of HLA-A, -B, -C, -DRB1, -DRB3/4/5, -DQA1, -DQB1, -DPA1, and -DPB1 by NGS-based typing in Koreans in South Korea. PLoS ONE 16:e0253619. doi:10.1371/journal.pone.0253619 | 2021 | South Korean | 173 | healthy Koreans | 11 loci | TruSight HLA v2 NGS | yes, Table 1 plus supplements | all six expected | yes as 11-locus support, but smaller n than Jung 2023 | medium | use as support |
| P2_BASE_BAEK2023 | Baek IC et al. Distributions of 11-loci HLA alleles typed by amplicon-based NGS in South Koreans. HLA 101:613-622. doi:10.1111/tan.14981 | 2023 | South Korean | not confirmed in abstract | healthy donors | 11 loci | amplicon-based NGS | likely yes | all six likely | possible but redundant if Jung 2023 + Baek 2021 suffice | low-medium | defer |
| P2_BASE_SONG2002 | Song EY et al. HLA class II allele and haplotype frequencies in Koreans based on 107 families. Tissue Antigens 59:475-486. doi:10.1034/j.1399-0039.2002.590604.x | 2002 | Korean families | 107 families; 207 parents + 291 children | family-based population reference | DRB1, DRB3/4/5, DQA1, DQB1, DPB1 | family SBT / class II haplotyping | yes if full table obtained | DPB1*05:01, DRB1*07:01, DQB1*02:01 possible; no A/B/C | class II background only | medium | use as background |

# 3. Paper 2 AITD / Hashimoto evidence registry

| source_id | citation | disease: HT / HD / GD / AITD / control | country / ancestry | n_case | n_control | loci | alleles relevant to 6-allele forest | effect direction | OR / p / frequency if reported | Paper 2 use |
|---|---|---|---|---:|---:|---|---|---|---|---|
| P2_AITD_SHIN2019 | Shin DH et al. HLA alleles, especially amino-acid signatures of HLA-DPB1, might contribute to the molecular pathogenesis of early-onset autoimmune thyroid disease. PLoS ONE 14:e0216941. doi:10.1371/journal.pone.0216941 | AITD split into GD and HD; controls | Korean children | 116 AITD total: 71 GD, 45 HD | 142 | A, B, C, DRB1, DQB1, DPB1 | GD: B*46:01, C*01:02, DPB1*05:01; HD: A*02:07 and DPB1*02:02; shared haplotype context includes A*02:07-B*46:01-C*01:02-DPB1*05:01 | HD/AITD background supports Korean thyroid-autoimmunity HLA substrate; GD rows must stay Paper 4/backlog | GD associations include B*46:01 OR 3.96 Pc=0.008, C*01:02 OR 2.51 Pc=0.04, DPB1*05:01 OR 4.6 Pc=0.003; HD includes A*02:07 OR 4.68 Pc=0.045; DPB1*02:02 OR 6.57 Pc=0.0001 | direct HD/AITD background; GD-specific rows Paper 4 only |
| P2_AITD_CHO2011 | Cho WK et al. Association of HLA alleles with autoimmune thyroid disease in Korean children. Horm Res Paediatr 76:328-334. doi:10.1159/000331134 | AITD split into HD and GD; controls | Korean children | 73 AITD total: 32 HD, 41 GD | 159 | class I and DRB1 | B*46, Cw*01, A*02, DRB1*07/08; lower resolution than 4-digit | HD and GD susceptibility overlap noted; use only as AITD/HD background in Paper 2 | reported higher B*46 and Cw*01 in HD; A*02/B*46/Cw*01/DRB1*08 higher in GD; DRB1*07 lower in GD | AITD background; not Pillar I quantitative fill |
| P2_AITD_UEDA2014_BACKGROUND | Ueda S et al. Identification of independent susceptible and protective HLA alleles in Japanese autoimmune thyroid disease and their epistasis. J Clin Endocrinol Metab 99:E379-E383. doi:10.1210/jc.2013-2841 | GD and HT; controls | Japanese | 991 AITD total: 547 GD, 444 HT | 481 | A, C, B, DRB1, DQB1, DPB1 | B*46:01 and DPB1*05:01 for GD; A*02:07 and DRB4 for HT per secondary citation | non-Korean support for disease split; not primary Paper 2 Korean pack | exact OR/p needs table extraction; PubMed verifies n and loci | citation background only; not direct Korean evidence |
| P2_AITD_ITO1989_BACKGROUND | Ito M et al. Association of HLA antigen and TCR beta-chain RFLP with Graves' disease and Hashimoto's thyroiditis. J Clin Endocrinol Metab 69:100-104. doi:10.1210/jcem-69-1-100 | GD and HT | Japanese | 61 GD, 50 HT | normal population comparator | serologic HLA antigen + TCR beta RFLP | Bw46 maps only broadly to B*46; no 4-digit resolution | supports old Japanese GD/HT shared HLA-Bw46 signal | Bw46 frequency 23.0% in GD, 24.0% in HT vs 8.0% normal population; corrected P reported | background only; drop from quantitative Paper 2 |
| P2_AITD_HUH1986 | Huh KB et al. Human leukocyte antigen in Korean patients with autoimmune thyroid diseases. Korean J Intern Med 1:243-248. doi:10.3904/kjim.1986.1.2.243 | AITD / Graves-heavy | Korean | not fully parsed | not fully parsed | serologic HLA | A11, DRw8, A10/B12; not 4-digit six-allele set | old Korean AITD/GD context only | no direct 4-digit extraction | drop from Paper 2 quantitative; optional historical background |

# 4. Paper 4 Pan-Asian GD registry

This table is backlog only. It is not a Paper 2 input table.

| source_id | citation | country | ancestry | n_GD | n_control | HLA loci | key alleles | amino-acid positions if reported | OR / p | summary-stat extractable? | backlog priority |
|---|---|---|---|---:|---:|---|---|---|---|---|---|
| P4_GD_CHU2018 | Chu X et al. Fine mapping MHC associations in Graves' disease and its clinical subtypes in Han Chinese. J Med Genet 55:685-692. doi:10.1136/jmedgenet-2017-105146 | China | Han Chinese | 1,468 | 1,490 | classical HLA imputation + amino-acid fine mapping | DPB1*05:01, B*46:01, DPA1*02:02, HLA-B and DR/DQ alleles | HLA-DP alpha Met11; HLA-B Lys66-Arg69-Val76 in pTRAb-negative GD; other residue signals in supplementary tables | published OR/p in tables; DPB1*05:01 and amino acid results extractable | yes | top |
| P4_GD_LIAO2022 | Liao WL et al. Analysis of HLA variants and Graves' disease and its comorbidities using a high resolution imputation system to examine electronic medical health records. Front Endocrinol 13:842673. doi:10.3389/fendo.2022.842673 | Taiwan | Taiwanese | 2,998 | control design needs extraction | HLA imputation, class I and II | A*02:07, B*46:01, C*01:02, DPA1/DPB1/DQA1/DQB1/DRB1 signals | not primary amino-acid paper | tables report genotype associations and comorbidity associations | yes, but phenotype/control definition needs careful extraction | high |
| P4_GD_SHIN2019 | Shin DH et al. PLoS ONE 14:e0216941. doi:10.1371/journal.pone.0216941 | Korea | Korean children | 71 | 142 | A, B, C, DRB1, DQB1, DPB1 | B*46:01, C*01:02, DPB1*05:01; DPB1*02:02; A*02:07 in HD | DPB1 Leu35 and Glu55 strong GD amino-acid signatures | B*46:01 OR 3.96 Pc=0.008; C*01:02 OR 2.51 Pc=0.04; DPB1*05:01 OR 4.6 Pc=0.003; Leu35/Glu55 OR 23.38 P=0.0002 | yes | top Korean bridge |
| P4_GD_PARK2005 | Park MH et al. Association of HLA-DR and -DQ genes with Graves disease in Koreans. Hum Immunol 66:741-747. doi:10.1016/j.humimm.2005.03.001 | Korea | Korean | 198 | 200 | DRB1, DQB1 | DRB1*08:03, DRB1*16:02, DQB1 alleles | none noted | abstract reports DRB1*08:03 OR 2.27 Pc=0.03 and other DR/DQ results in tables | yes | high |
| P4_GD_CHO1987 | Cho BY et al. HLA and Graves' disease in Koreans. Tissue Antigens 30:119-121. doi:10.1111/j.1399-0039.1987.tb01607.x | Korea | Korean | 128 | 220 | serologic A, B, C, DR | B13, DR5, DRw8 | none | relative risks reported: B13, DR5, DRw8 increased | partial; old serology | medium historical |
| P4_GD_JANG2011 | Jang HW et al. Identification of HLA-DRB1 alleles associated with Graves' disease in Koreans by sequence-based typing. Immunol Invest 40:172-182. doi:10.3109/08820139.2010.525571 | Korea | Korean | needs full-text extraction | controls included | DRB1 | DRB1 alleles | none expected | extract from table only | likely yes | medium |
| P4_GD_CHEN2011 | Chen PL et al. Comprehensive genotyping in two homogeneous Graves' disease samples reveals major and novel HLA association alleles. PLoS ONE 6:e16635. doi:10.1371/journal.pone.0016635 | Taiwan | ethnic Chinese Han | two samples; exact n by stage needs table extraction | population controls | A, B, C, DQB1, DRB1, DPB1 | DPB1*05:01, B*46:01, DQB1*03:02, DRB1*15:01, DRB1*16:02, DRB1*12:02 | no amino-acid fine mapping focus | published combined model: DPB1*05:01 OR 2.34 P(Bc)=2.58e-10; B*46:01 OR 1.33 P(Bc)=1.17e-2; other alleles reported | yes | high |
| P4_GD_UEDA2014 | Ueda S et al. J Clin Endocrinol Metab 99:E379-E383. doi:10.1210/jc.2013-2841 | Japan | Japanese | 547 | 481 | A, C, B, DRB1, DQB1, DPB1 | GD: B*35:01, B*46:01, DRB1*14:03, DQB1*06:04, DPB1*05:01; HT: A*02:07, DRB4 | epistasis; exact residue positions not primary | OR/p table extractable from paper | yes | high |
| P4_GD_ONUMA1994 | Onuma H et al. Association of HLA-DPB1*0501 with early-onset Graves' disease in Japanese. Hum Immunol 39:195-201. doi:10.1016/0198-8859(94)90260-7 | Japan | Japanese | needs table extraction | needs table extraction | DPB1 and B | DPB1*05:01, B46 | none | published association | likely yes | medium |
| P4_GD_NAITO1987 | Naito S, Sasaki H, Arakawa K. Japanese Graves' disease: association with HLA-Bw46. Endocrinol Jpn 34:685-688. doi:10.1507/endocrj1954.34.685 | Japan | Japanese | 61 for A/B/C; 53 for DR/DQ | 1,998 | serologic A, B, C, DR, DQ | Bw46; CX46 | none | Bw46 23.0% vs 8.4%, Pc<0.003 | partial; old serology | medium historical |
| P4_GD_TSAI1989 | Tsai KS et al. Association of HLA-DR tissue types with Graves' disease in Taiwan. Taiwan Yi Xue Hui Za Zhi 88:336-341. PMID:2794934 | Taiwan | Chinese residing in Taiwan | 93 | 106 | serologic DR/DQ | DR2, DR9, DQw1, DRw53; DR3/DRw52 lower | none | DR2/GD 40.9% vs 21.7% Pc=0.028; DR9 29.0% vs 13.2% Pc=0.048; DQw1/DRw53 increased | partial; old serology | low-medium historical |
| P4_GD_CAVAN1994 | Cavan DA et al. The HLA association with Graves' disease is sex-specific in Hong Kong Chinese subjects. Clin Endocrinol 40:63-66. doi:10.1111/j.1365-2265.1994.tb02444.x | Hong Kong | Chinese | needs abstract/table extraction | needs extraction | A, B, DR, DQ | likely Bw46/DR/DQ signals | none | sex-specific results in full text | likely partial | defer |
| P4_GD_HUANG2021_GO | Huang X et al. Human leucocyte antigen alleles confer susceptibility and progression to Graves' ophthalmopathy in a Southern Chinese population. Br J Ophthalmol 105:1462-1468. doi:10.1136/bjophthalmol-2020-317091 | China | Southern Han Chinese | GD/GO split needs extraction | controls within n=683 total | A, B, C, DRB1, DQB1, DQA1, DPA1, DPB1 | GO/GD alleles; exact six-allele overlap needs table extraction | not primary Chu-style fine mapping | four-digit SBT association tables | yes | medium for GO backlog |
| P4_GD_B46_META2013 | Li Y et al. Association between HLA-B*46 allele and Graves disease in Asian populations: a meta-analysis. Int J Med Sci 10:164-170. doi:10.7150/ijms.5158 | Asian populations | mixed Asian | 1,743 | 5,689 | HLA-B | B*46 | none | pooled OR reported; heterogeneity reported | yes, review/meta context only | medium |
| P4_GD_DRB1_META2024 | Li W et al. Association between HLA-DRB1 alleles and Graves' disease in Asian populations: a meta-analysis. Horm Metab Res 56:859-868. doi:10.1055/a-2298-4366 | Asian populations | mixed Asian | meta-analysis | meta-analysis | DRB1 | DRB1 alleles | none | pooled DRB1 association summary | yes, review/meta context only | medium |
| P4_GD_REVIEW2023 | Stasiak M et al. Significance of HLA in Graves' disease and Graves' orbitopathy in Asian and Caucasian populations: a systematic review. Front Immunol 14:1256922. doi:10.3389/fimmu.2023.1256922 | Asian/Caucasian | mixed | review | review | broad HLA | country-specific GD and GO HLA alleles | review table only | Table 1 indexes Asian GD studies and sample sizes | yes for source discovery | high registry map |

# 5. Allele harmonization rules

1. Normalize every extracted allele to 4-digit two-field form: `A*02:07`, `B*46:01`, `C*01:02`, `DRB1*07:01`, `DQB1*02:01`, `DPB1*05:01`.

2. Collapse 6-digit or 8-digit calls to 4-digit only with a `collapsed_from` note. Example: `C*01:02:01` -> `C*01:02`.

3. Do not merge carrier frequency and allele frequency. Published population references usually report allele frequency across `2n`; arcasHLA cohort summaries may be carrier/sample frequency.

4. Keep genotype frequency separate from allele frequency. Liao 2022 and some older GD papers report genotype-level combinations; these must not be treated as per-allele frequencies.

5. Preserve control ancestry labels exactly: Korean, South Korean, Han Chinese, Taiwanese, Japanese, Hong Kong Chinese, Southern Han Chinese, Singaporean Chinese, Harbin Korean proxy.

6. Disease labels are strict: `HT`, `HD`, `GD`, `AITD`, `control`, `PTC`. A mixed AITD paper must be split into disease-specific rows before any later extraction.

7. Serologic results such as `Bw46`, `B13`, `DRw8`, `DQw1` are historical context only unless a paper gives allele-level conversion. Do not merge serologic data into a 4-digit forest source table.

# 6. Immediate Paper 2 action candidates

| candidate | verdict | reason |
|---|---|---|
| Use Lee/In/Jung/Jekarl as Korean baseline pack | go for extraction planning | collectively covers old and modern Korean baseline, including large A/B/C/DRB1 and 11-locus DP/DQ support |
| Use In 2015 to fill C*01:02 / DQB1*02:01 | go for table extraction planning | exact Table 1 includes both alleles and n=613; no computation yet |
| Use Jung 2023 to fill DPB1*05:01 and cross-check all six | go after full-table access | n=339 11-locus NGS is the best modern Korean all-locus baseline candidate |
| Use Jekarl 2021 for A/B/C/DRB1 cross-check | go after table access | n=1,293 healthy donors; no DQB1/DPB1 |
| Use Shin 2019 as AITD/HT/GD background with strict separation | go with restriction | HD/AITD rows can support Paper 2 context; GD rows must remain Paper 4/backlog |
| Decide whether DRB1*07:01 Harbin fallback moves to sensitivity | go for decision | Korean baseline sources now exist for DRB1*07:01, so Harbin proxy should not drive primary Paper 2 evidence |

# 7. Paper 4 future action candidates

| candidate | verdict |
|---|---|
| Chu 2018 + Liao 2022 + Shin 2019 + Park 2005 registry | go as backlog registry only |
| Add Japanese Ueda 2014 / Onuma 1994 / Naito 1987 to backlog | go as source discovery |
| Add Taiwanese Chen 2011 / Tsai 1989 and Hong Kong Cavan 1994 to backlog | go as source discovery |
| Use Stasiak 2023 systematic review as source map | go as registry index |
| Run meta-analysis or forest | no-go until Paper 4 gates clear |

# 8. Risk notes

| risk | note |
|---|---|
| allele frequency vs carrier frequency | never mix AF across 2n chromosomes with carrier/sample frequency without explicit conversion plan |
| disease mismatch HT vs GD | Paper 2 is Hashimoto/HD/AITD context; GD evidence belongs to Paper 4 unless used only as disease-separation background |
| ancestry mismatch | Korean, Han Chinese, Taiwanese, Japanese, Hong Kong Chinese, and Southern Han Chinese are separate labels |
| old low-resolution typing | serologic studies are useful for history but not for 4-digit allele-level extraction |
| sample-size imbalance | Jekarl 2021 n=1,293, In 2015 n=613, Lee 2005 n=485, Jung 2023 n=339, Baek 2021 n=173; do not over-weight small modern NGS tables without disclosure |
| table access | some Wiley/HLA articles require full-text access for exact table values; registry status should remain `verify` until table-level values are inspected |

# 9. Next command options

Allowed next commands:

```text
extract Paper2 Korean baseline tables only
```

Published Korean baseline frequency extraction only. No OR computation. No forest.

```text
prepare 6-allele forest source table
```

Create source/provenance table schema only. Do not compute forest statistics.

```text
freeze Paper4 GD registry
```

Keep all GD sources as backlog only.

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

# Source anchors checked

- Lee 2005: `https://pubmed.ncbi.nlm.nih.gov/15853898/`
- In 2015: `https://synapse.koreamed.org/upload/synapsexml/3039alm/pdf/alm-35-429.pdf`
- Jung 2023: `https://pubmed.ncbi.nlm.nih.gov/36719349/`
- Jekarl 2021: `https://pubmed.ncbi.nlm.nih.gov/33314756/`
- Chung 2010: `https://www.kci.go.kr/kciportal/ci/sereArticleSearch/ciSereArtiView.kci?sereArticleSearchBean.artiId=ART001497939`
- Shin 2019: `https://journals.plos.org/plosone/article/file?id=10.1371%2Fjournal.pone.0216941&type=printable`
- Cho 2011: `https://cuk.elsevierpure.com/en/publications/association-of-hla-alleles-with-autoimmune-thyroid-disease-in-kor`
- Chu 2018: `https://pubmed.ncbi.nlm.nih.gov/29987165/`
- Liao 2022: `https://www.frontiersin.org/articles/10.3389/fendo.2022.842673/full`
- Park 2005: `https://pubmed.ncbi.nlm.nih.gov/15993720/`
- Cho 1987: `https://snu.elsevierpure.com/en/publications/hla-and-graves-disease-in-koreans`
- Ueda 2014: `https://pubmed.ncbi.nlm.nih.gov/24285682/`
- Chen 2011: `https://pubmed.ncbi.nlm.nih.gov/21307958/`
- Naito 1987: `https://www.jstage.jst.go.jp/article/endocrj1954/34/5/34_5_685/_article`
- Tsai 1989: `https://pubmed.ncbi.nlm.nih.gov/2794934/`
- Cavan 1994: `https://pubmed.ncbi.nlm.nih.gov/8306482/`
- Stasiak 2023 systematic review: `https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2023.1256922/pdf`
