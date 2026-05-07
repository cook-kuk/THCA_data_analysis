---
title: "HLA external evidence acquisition registry"
date: 2026-05-04
status: "REGISTRY / PLAN ONLY - no analysis execution"
scope:
  paper2: "Hashimoto-overlap PTC HLA susceptibility + HLA-II-mediated dedifferentiation"
  paper4: "Korean GD HLA / Pan-Asian Graves backlog"
hard_boundaries:
  - "GD / Graves / Bundang are forbidden in Paper 2 execution."
  - "Paper 4 sources are registry backlog only."
  - "Paper 3 ICI is frozen."
forbidden:
  - "new large data download"
  - "cookHLA SNP integration"
  - "Bundang Graves work"
  - "Chu 2018 GD meta-analysis execution"
  - "Paper 4 main analysis execution"
  - "Paper 3 touch"
  - "manuscript prose"
  - "voice-protected sections"
  - "H&E / WSI / pathology image analysis"
  - "TCGA WSI"
  - "RunPod / GPU"
allowed:
  - "PubMed / PMC / journal landing page read-only lookup"
  - "citation metadata verification"
  - "HLA frequency / published OR / p-value / sample-size extraction plan"
  - "registry table"
  - "execution-ready plan without execution"
---

# 1. Executive verdict

Paper 2 HLA can be strengthened now by a Korean baseline plus HT/AITD reference pack. The immediate evidence pack should use Korean general-population HLA references for baseline frequency support and Korean early-onset AITD/HD references for disease-background support.

Paper 4 GD HLA can be strengthened only as a backlog registry. Chu 2018, Taiwan GD, Korean Graves historical studies, Japanese GD, and Asian GD meta-analyses should be indexed for future extraction, but not executed or merged into Paper 2.

Critical boundary: Paper 2 may cite HT/HD/AITD substrate evidence and Korean baseline frequency evidence. It must not run GD/Graves/Bundang work, must not re-activate Chu 2018 as a Paper 2 comparator, and must not touch Paper 3.

# 2. Paper 2 evidence registry

Target allele set: DPB1*05:01, A*02:07, B*46:01, DRB1*07:01, C*01:02, DQB1*02:01.

| source_id | citation | population | sample_size | disease/control status | HLA loci | alleles relevant to our 6-allele forest | exact frequency table available | can fill Pillar I v2? | priority | action |
|---|---|---:|---:|---|---|---|---|---|---|---|
| P2-KR-BASE-IN2015 | In JW et al. Ann Lab Med 2015;35:429-435. doi:10.3343/alm.2015.35.4.429 | Korean | 613 | healthy unrelated donors / cord blood units | A, B, C, DRB1, DQB1 | A*02:07 yes; B*46:01 yes; DRB1*07:01 yes; C*01:02 yes; DQB1*02:01 yes; DPB1 absent | yes, Table 1 | yes for 5/6; not DPB1 | high | use |
| P2-KR-BASE-LEE2005 | Lee KW et al. Tissue Antigens 2005;65:437-447. doi:10.1111/j.1399-0039.2005.00386.x | Korean | 485 | apparently unrelated healthy individuals | A, B, C, DRB1, DQB1 | covers A/B/C/DRB1/DQB1 set; DPB1 absent; some legacy notation/ambiguity possible | likely yes, published allele/haplotype frequencies | yes as cross-check for 5/6; weaker than In 2015 for clean fill | medium | verify |
| P2-KR-BASE-BAEK2021 | Baek IC et al. PLoS ONE 2021;16:e0253619. doi:10.1371/journal.pone.0253619 | South Korean | 173 | healthy Koreans | A, B, C, DRB1/3/4/5, DQA1, DQB1, DPA1, DPB1 | all six covered, including DPB1 and DQB1; Table 1 lists >1%, S9/S10 full tables | yes, Table 1 plus supplements | yes, but small n; best as 11-locus NGS cross-check | high for DPB1/DQ support; medium as primary | use as support, not sole primary |
| P2-KR-BASE-BAEK2023 | Baek IC et al. HLA 2023;101:613-622. doi:10.1111/tan.14981 | South Korean | not confirmed from abstract | healthy donors | 11 loci by amplicon NGS | likely all six covered; exact sample size/table access needs full-text check | likely yes | possible, but not needed before Baek 2021/In 2015 | medium | defer |
| P2-KR-BASE-JEKARL2021 | Jekarl DW et al. HLA 2021;97:??. doi needs full verification; PMID 33314756 | Korean | 1,293 | unrelated healthy donors | A, B, C, DRB1 | A*02:07, B*46:01, C*01:02, DRB1*07:01; no DQB1/DPB1 | yes if full table obtained | yes for class I + DRB1 only; cannot fill DQB1/DPB1 | medium | verify |
| P2-KR-BASE-CHOE2021 | Choe W et al. Ann Lab Med 2021;41:310-317. doi:10.3343/alm.2021.41.3.310 | Korean | 128 | healthy unrelated adults | A, B, C, DRB1 | A*02:07 possible; B*46:01 possible; C*01:02 yes; DRB1*07:01 possible; no DQB1/DPB1 | yes | no for full Pillar I; useful for high-resolution C/DR cross-check | low-medium | citation background / verify |
| P2-KR-BASE-SONG2002 | Song EY et al. Tissue Antigens 2002;59:475-486. doi:10.1034/j.1399-0039.2002.590604.x | Korean families | 107 families; 207 parents + 291 children | family-based population reference | DRB1, DRB3/4/5, DQA1, DQB1, DPB1 | DPB1*05:01, DRB1*07:01, DQB1*02:01 possible; no A/B/C | yes, class II allele/haplotype tables | yes for class II context; not full 6-allele | medium | use as class II haplotype background |
| P2-AITD-SHIN2019 | Shin DH et al. PLoS ONE 2019;14:e0216941. doi:10.1371/journal.pone.0216941 | Korean children | 116 AITD: 71 GD, 45 HD; 142 controls | disease + control; includes HD | A, B, C, DRB1, DQB1, DPB1 | GD: B*46:01, C*01:02, DPB1*05:01; HD: A*02:07, DPB1*02:02; DPB1*05:01 differs GD vs HD | yes, manuscript + supplements | yes for HT/HD/AITD disease substrate, not baseline primary | high | use for Paper 2 HT/HD context only; do not use GD arm for Paper 2 analysis |
| P2-AITD-CHO2011 | Cho WK et al. Horm Res Paediatr 2011;76:328-334. doi:10.1159/000331134 | Korean children | 73 AITD: 32 HD, 41 GD; 159 controls | disease + control; includes HD | HLA class I and DRB1 | B*46 and Cw*01 reported in HD/AITD; A*02/B*46/Cw*01/DRB1*08 in GD; DRB1*07 lower in GD | yes, likely allele frequencies | supportive for early-onset AITD overlap; lower resolution | medium | use as citation background / verify |
| P2-LEE2014-KOTRY | "Lee 2014 Tissue Antigens / KOTRY donor reference" per local prior notes | Korean | not verified in this pass | presumed donor/reference | presumed multi-locus HLA | reported locally as needed for C*01:02 / DQB1*02:01 fill and DPB1*05:01 baseline cross-check | not verified in this pass | do not use until citation identity and table are verified | high but blocked | verify before any extraction |

# 3. Paper 2 specific outputs

## Korean baseline panel candidate sources

Primary candidate pack:

| role | source_id | rationale | caveat |
|---|---|---|---|
| clean 5-locus baseline | P2-KR-BASE-IN2015 | high-resolution SBT, n=613, exact table includes A*02:07, B*46:01, C*01:02, DRB1*07:01, DQB1*02:01 | no DPB1 |
| 11-locus NGS baseline | P2-KR-BASE-BAEK2021 | includes DPB1 and DQB1 with supplements; clean nomenclature | n=173; use as support/cross-check |
| legacy 5-locus baseline | P2-KR-BASE-LEE2005 | n=485, historically cited Korean reference | older typing/notation; verify exact allele resolution |
| class II family/haplotype context | P2-KR-BASE-SONG2002 | DPB1/DQB1/DRB1 haplotype substrate | family-based; not direct unrelated-control frequency panel |
| large class I/DRB1 donor panel | P2-KR-BASE-JEKARL2021 | n=1,293 and NGS, useful for A/B/C/DRB1 | no DQB1/DPB1 |

The "Baek 2021 n=339" candidate did not verify as written. The verified Baek 2021 PLOS ONE 11-locus NGS paper is n=173. A separate 2023 Baek HLA paper exists, but sample size was not confirmed from accessible metadata in this pass.

## HT/HD/AITD evidence sources

| role | source_id | Paper 2 use |
|---|---|---|
| Korean early-onset HD/AITD reference | P2-AITD-SHIN2019 | supports Korean AITD HLA background and HD-specific allele context; use HD/AITD parts only for Paper 2 |
| Korean pediatric AITD background | P2-AITD-CHO2011 | supports overlap between HD and GD susceptibility alleles in Korean children; use as background only |
| Korean class II haplotype substrate | P2-KR-BASE-SONG2002 | supports DPB1/DQB1/DRB1 Korean haplotype interpretation |

## Which sources support Lee 2014 fill

The current local plan names "Lee 2014 Tissue Antigens / KOTRY" as the fill source, but this registry pass did not verify the exact bibliographic identity or table. Do not extract from "Lee 2014" until citation, DOI/PMID, sample size, locus coverage, and frequency-table location are verified.

If Lee 2014 verification stalls, the nearest clean fill alternatives are:

| fill target | immediate candidate | reason |
|---|---|---|
| C*01:02 | In 2015 | exact high-resolution Korean AF table; n=613 |
| DQB1*02:01 | In 2015 | exact high-resolution Korean AF table; n=613 |
| DPB1*05:01 cross-check | Baek 2021 + Song 2002 | 11-locus NGS and class II haplotype context |
| A*02:07 / B*46:01 / DRB1*07:01 cross-check | In 2015 + Lee 2005 + Jekarl 2021 | independent Korean baseline references |

## Citation background only

Choe 2021 is useful for 8-digit A/B/C/DRB1 background, but it cannot fill the full 6-allele forest because DQB1 and DPB1 are absent. Baek 2023 should stay "defer" unless the full text/table is needed after the In 2015/Baek 2021 pack.

# 4. Paper 4 GD evidence registry

This section is backlog only. No Paper 4 main analysis, no Chu 2018 execution, no Bundang work, and no GD/Graves integration into Paper 2.

| source_id | citation | country / ancestry | sample_size_case | sample_size_control | disease subtype | HLA loci | relevant alleles | OR / p / amino acid positions if reported | summary-stat extractable | paper4 priority | action |
|---|---|---|---:|---:|---|---|---|---|---|---|---|
| P4-GD-CHU2018 | Chu X et al. J Med Genet 2018;55:685-692. doi:10.1136/jmedgenet-2017-105146; PMID 29987165 | Han Chinese | 1,468 | 1,490 | GD, clinical subtypes | classical HLA imputation and amino acid fine mapping | DPB1*05:01, B*46:01, A*02:07, C*01:02, DRB1/DQB1 alleles | published allele OR/p and amino acid residue associations; B residues Lys66-Arg69-Val76 noted for pTRAb-negative GD | yes | top | backlog use |
| P4-GD-SHIN2019 | Shin DH et al. PLoS ONE 2019;14:e0216941 | Korea | 71 GD; 45 HD also present | 142 | early-onset AITD; GD/HD split | A, B, C, DRB1, DQB1, DPB1 | GD: B*46:01, C*01:02, DPB1*05:01; HD: A*02:07, DPB1*02:02; DPB1*05:01 GD vs HD difference | published OR/Pc for selected alleles and DPB1 Leu35/Glu55 amino acid signal | yes | top for Korean GD/HD bridge | backlog use; Paper 2 may use HD context only |
| P4-GD-LIAO2022 | Liao WL et al. Front Endocrinol 2022;13:842673. doi:10.3389/fendo.2022.842673 | Taiwanese | 2,998 | EMR/control design needs verification | GD + comorbidities | HLA imputation, class I and II | A*02:07, B*46:01, C*01:02, DPB1*05:01, DQB1/DRB1/DPA1/DQA1 variants | published genotype associations; tables include GD and comorbidity associations | yes, but phenotype/control design needs careful extraction | high | backlog use |
| P4-GD-PARK2005 | Park MH et al. Hum Immunol 2005;66:741-747. doi:10.1016/j.humimm.2005.03.001 | Korea | 198 | 200 | GD | DRB1, DQB1 | DRB1*08:03, DRB1*16:02; DQB1 associations | published OR/corrected p in abstract/table | yes | high for Korean historical GD | backlog use |
| P4-GD-CHO1987 | Cho BY et al. Tissue Antigens 1987;30:119-121. doi:10.1111/j.1399-0039.1987.tb01607.x | Korea | 128 | 220 | GD | serologic HLA-A, B, C, DR | B13, DR5, DRw8; not 4-digit | relative risks reported for serologic antigens | partial; not 4-digit harmonized | medium | backlog historical background |
| P4-GD-JANG2011 | Jang HW et al. Immunol Invest 2011;40:172-182. doi:10.3109/08820139.2010.525571 | Korea | full case/control split needs table verification | likely controls included | GD | DRB1 by SBT | DRB1 alleles | published DRB1 association table | likely yes | medium | defer/verify |
| P4-GD-CHEN2011 | Chen PL et al. PLoS ONE 2011;6:e16635. doi:10.1371/journal.pone.0016635 | Taiwan | two homogeneous GD samples; exact per-stage n needs extraction | controls in study | GD | comprehensive HLA genotyping | DPB1*05:01, B*46:01, DQB1*03:02, DRB1*15:01, DRB1*16:02, DRB1*12:02 | published allele model and population-attributable risk; DPB1*05:01 noted as major | yes | high for Taiwan historical GD | backlog use |
| P4-GD-UEDA2014 | Ueda S et al. J Clin Endocrinol Metab 2014;99:E379-E383. doi:10.1210/jc.2013-2841 | Japan | exact n needs extraction | exact n needs extraction | Japanese AITD, GD/HD | HLA alleles | B*35:01, B*46:01, DRB1*14:03, DQB1*06:04, DPB1*05:01 for GD; A*02:07 and DRB4 for HD per Shin 2019 citation | independent susceptible/protective alleles and epistasis | likely yes | high for Japan backlog | defer/verify |
| P4-GD-ONUMA1994 | Onuma H et al. Hum Immunol 1994;39:195-201. PMID 8026987 | Japan | exact n needs extraction | exact n needs extraction | early-onset GD | DPB1 and B | DPB1*05:01, B46 | published association | likely yes | medium | defer/verify |
| P4-GD-LI2013 | Li Y et al. Int J Med Sci 2013;10:164-170. doi:10.7150/ijms.5158 | Asian populations | 1,743 | 5,689 | GD | HLA-B*46 | B*46 | published pooled OR and heterogeneity | yes | medium | backlog use for context only |
| P4-GD-LI2024 | Li W et al. Horm Metab Res 2024;56:859-868. doi:10.1055/a-2298-4366 | Asian populations | meta-analysis | meta-analysis | GD | DRB1 | DRB1 alleles | published meta-analysis summary | yes | medium | backlog use |
| P4-GD-SYSREV2023 | Significance of HLA in Graves' disease and Graves' orbitopathy in Asian and Caucasian populations. PMC10568027 | Asian and Caucasian | systematic review | systematic review | GD/GO | broad HLA | country-specific GD/GO alleles | table of Asian studies including Cho 1987, Park 2005, Jang 2011 | yes for registry triage | medium | backlog citation map |

# 5. Allele harmonization plan

1. Standard allele naming: normalize to 4-digit two-field form for the registry, e.g. `DPB1*05:01`, `A*02:07`, `B*46:01`, `DRB1*07:01`, `C*01:02`, `DQB1*02:01`.

2. Preserve higher-resolution originals: keep original 3-field/4-field values in a separate `reported_allele` column during extraction, then map to `allele_4digit`.

3. Carrier frequency vs allele frequency: do not mix without an explicit flag. PTC arcasHLA outputs often behave as carrier/genotype-derived sample frequency; population references often report allele frequency across 2n chromosomes. Extraction tables must include `frequency_type = carrier | allele | genotype | phenotype | unclear`.

4. Phenotype frequency warning: older serologic studies report antigen/phenotype frequencies such as B13, DR5, or DRw8. These are not equivalent to 4-digit allele frequencies and should not be merged with 4-digit allele-level forests.

5. Ancestry label: use explicit labels, e.g. `Korean`, `South Korean`, `Han Chinese`, `Taiwanese`, `Japanese`, `Korean children`, `Korean families`, `Harbin Korean proxy`. Do not collapse to "Asian" in extraction tables.

6. Disease label: use controlled labels only: `HT`, `HD`, `GD`, `AITD`, `control`, `PTC`. For Paper 2 extraction, `GD` rows may be present in registry metadata but must be excluded from Paper 2 analysis execution.

7. Table provenance: every extracted row should carry `source_id`, `table_or_supplement`, `page_or_table_label`, `n_case`, `n_control`, `frequency_type`, `typing_method`, and `extraction_status`.

# 6. Paper 2 go/no-go decisions

| decision | verdict | reason |
|---|---|---|
| Proceed with Lee/In/Baek/Jekarl reference pack? | go for registry and extraction prep only | In 2015 and Baek 2021 are verified enough for planning; Jekarl and Lee 2014 need citation/table verification before extraction |
| Keep AFND + KOTRY/Lee as primary? | conditional go | AFND remains primary existing baseline; KOTRY/Lee can join only after Lee 2014 identity/table verification |
| Use In 2015 as immediate backup fill? | go for plan | In 2015 has exact C*01:02 and DQB1*02:01 frequency table and is cleaner than an unverified Lee 2014 citation |
| Move DRB1*07:01 to sensitivity? | go | DRB1*07:01 is ancestry/haplotype sensitive and current baseline has Harbin Korean proxy in existing v2; primary claim should not depend on it |
| Expand 4-allele to 6-allele? | conditional go | expand only after baseline source decision and exact frequency extraction; do not run forest in this registry step |
| Use Shin 2019 in Paper 2? | go with restriction | use HD/AITD context; do not use GD arm for Paper 2 analysis |
| Use Chu 2018 in Paper 2? | no-go | Paper 4 backlog only |

# 7. Paper 4 future plan

Pan-Asian GD HLA meta-analysis design, future only:

1. Build a source registry from Chu 2018, Park 2005, Cho 1987, Jang 2011, Shin 2019, Chen 2011, Liao 2022, Ueda 2014, Onuma 1994, Li 2013, Li 2024, and the 2023 systematic review.

2. Extract only published summary-stat evidence: allele/genotype counts, reported OR, p/Pc, confidence interval, typing method, and ancestry.

3. Separate evidence strata before any future analysis: Korean GD, Han Chinese GD, Taiwanese GD, Japanese GD, pediatric GD, adult GD, Graves orbitopathy, GD comorbidities.

4. Harmonize allele resolution and frequency type before cross-study comparison.

5. Do not execute until Paper 4 gates clear.

Bundang cohort integration plan, future only:

| step | status |
|---|---|
| Bundang sample arrival and n check | blocked |
| Korean GD phenotype harmonization | blocked |
| HLA typing/imputation plan | blocked |
| Korean GD vs Korean control comparison | blocked |
| Pan-Asian GD meta-analysis | blocked |

# 8. Immediate next commands

Allowed next command set:

```text
run Paper2 baseline extraction only
```

Meaning: extract published Korean baseline frequencies only from verified baseline references. Do not compute OR. Do not run a forest.

```text
prepare 6-allele forest script only
```

Meaning: prepare script skeleton and input schema only. Do not execute the script.

```text
freeze Paper4 registry
```

Meaning: keep Paper 4 sources in backlog and do not touch analysis.

```text
stop and return to Paper1 Hook
```

Meaning: leave HLA work frozen and return to Paper 1 Hook tasks.

Blocked commands:

```text
run Chu meta-analysis
run Paper4 analysis
run Bundang Graves
run cookHLA SNP integration
run forest
edit manuscript
touch Paper3
```

# Verification notes and source anchors

Read-only sources checked in this pass:

- In 2015, Annals of Laboratory Medicine PDF: `https://synapse.koreamed.org/upload/synapsexml/3039alm/pdf/alm-35-429.pdf`
- Lee 2005 PubMed: `https://pubmed.ncbi.nlm.nih.gov/15853898/`
- Baek 2021 PLOS ONE PDF: `https://journals.plos.org/plosone/article/file?id=10.1371%2Fjournal.pone.0253619&type=printable`
- Baek 2023 PubMed/Ovid metadata: `https://pubmed.ncbi.nlm.nih.gov/36720674/`
- Jekarl 2021 PubMed metadata: `https://pubmed.ncbi.nlm.nih.gov/33314756/`
- Choe 2021 KCI / Ann Lab Med metadata: `https://www.kci.go.kr/kciportal/ci/sereArticleSearch/ciSereArtiView.kci?sereArticleSearchBean.artiId=ART002718219`
- Song 2002 SNU metadata: `https://snu.elsevierpure.com/en/publications/hla-class-ii-allele-and-haplotype-frequencies-in-koreans-based-on`
- Shin 2019 PLOS ONE PDF: `https://journals.plos.org/plosone/article/file?id=10.1371%2Fjournal.pone.0216941&type=printable`
- Cho 2011 institutional metadata: `https://cuk.elsevierpure.com/en/publications/association-of-hla-alleles-with-autoimmune-thyroid-disease-in-kor`
- Chu 2018 PubMed/PMC metadata: `https://pubmed.ncbi.nlm.nih.gov/29987165/`
- Liao 2022 Frontiers: `https://www.frontiersin.org/articles/10.3389/fendo.2022.842673/full`
- Park 2005 PubMed: `https://pubmed.ncbi.nlm.nih.gov/15993720/`
- Cho 1987 institutional metadata: `https://snu.elsevierpure.com/en/publications/hla-and-graves-disease-in-koreans`
- Chen 2011 PubMed/PMC metadata: `https://pubmed.ncbi.nlm.nih.gov/21307958/`
- Asian HLA-B*46 GD meta-analysis: `https://pmc.ncbi.nlm.nih.gov/articles/PMC3547214/`
- Asian DRB1 GD meta-analysis PubMed: `https://pubmed.ncbi.nlm.nih.gov/38698581/`
- GD/GO systematic review: `https://pmc.ncbi.nlm.nih.gov/articles/PMC10568027/`
