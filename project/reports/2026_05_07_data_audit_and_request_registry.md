# Data audit + request registry (2026-05-07)

**Purpose:** (1) 현재 hub data dictionary v1/v2 가 빠뜨린 cohort 점검, (2) 더 가져올 public 데이터 후보, (3) **공개됐지만 다운로드 불가 / application/PI request 필요한 데이터** 전체 정리. Paper 1–4 + 5–12 portfolio 통합 view.

**원칙:** 분석 명령 아니라 **registry 만**. 마라톤 모드 무엇도 위반 안 함. 신청 결정은 사용자 권한.

---

## 0. TL;DR — 빠진 것 + 신청 후보 한눈에

| 카테고리 | 누락 / 차단 cohort |
|---|---|
| **사전 누락 (이미 manuscript 사용 중)** | GSE213647 (Lee 2024 Korean n=632), MSK-IMPACT thyroid (Landa 2016 n=117), GSE184362 (Pu 2021 sc n=6), GSE193581 (Lu 2023 n=23), GSE241184 (Phase 1 sc n=1), cBioPortal `thca_tcga_pub`, HM450 TCGA-THCA (n=503) |
| **추가 가져올 가치 있는 public** | GSE78220 (Hugo), GSE91061 (Riaz), GSE60542 (PTC nodal n=92), GSE151179 (post-RAI n=52), GSE27155, GSE97466, Pozdeyev 2018 CCR (n=779 advanced) supp table |
| **dbGaP controlled — application 필요** | TCGA-THCA paired tumor-normal BAM, Liu phs000452 (멜라노마 ICI WES+RNA), TCGA-CDR phs000178 |
| **EGA — application 필요** | EGAS00001003540, EGAS00001002556, EGAD00001004845, Yoo SK 2019 Korean ATC (EGA likely) |
| **PI / pharma request — 영원히 비공개 가능** | KEYNOTE-158 thyroid raw, KEYNOTE-028 thyroid raw, Spartalizumab ATC (NCT02404441 Novartis), Sehgal lenva+pembro, Dierks lenva+pembro, Cabanillas thyroid ICI, NCT03246958 nivo+ipi, NCT03181100 atezo |
| **Korean institutional outreach** | Bundang SNUH Graves' BTC, 연세 / 서울대 / 삼성 thyroid, KARE / KoGES SNP cohort |
| **Population HLA — open but Korean-specific aggregation needed** | AFND South Korea pool (이미 사용), Chen 2018 Han Chinese GD HLA fine-mapping (PMC supp 또는 dbGaP), Liao 2022 / Shin 2019 / Park 2005 / Cho 1987 supplementary tables |

---

## 1. 현재 v1/v2 dictionary 누락 cohort (이미 manuscript / reports에서 사용 중인데 표에 안 들어감)

| Accession / Name | Source | Modality | n | 사용처 | 추가해야 함 |
|---|---|---|---|---|---|
| **GSE213647** | Lee 2024 Korean PTC | RNA-seq | **632** | Paper 1 Korean external validation main; 22-28% Hashimoto-like ≈ TCGA 18-20%; Suppl S6/S7 | ★★★ critical — manuscript 04_results.md / 06_discussion.md / 05_supp_figure_captions에서 자주 등장 |
| **MSK-IMPACT thyroid (Landa 2016 Cell)** | MSK-IMPACT | targeted DNA + clinical | **117** (84 PDTC + 33 ATC) | Paper 1 advanced-disease anchor; Cox HR 2.67; cBioPortal `thca_mskcc_2016` | ★★★ — GSE76039 만 사전에 있고 정작 본 cohort은 빠짐 |
| **GSE184362** | Pu 2021 Nat Commun, Fudan | scRNA | **6 PTC patients, 158K cells** | Paper 1 single-cell external; per-patient r=0.798–0.886 thyrocyte-intrinsic; Fig SB1/SC1 | ★★★ — series_matrix 다운로드 됨 (`paper3_ici_public_data/GSE184362/`), RAW.tar (970 MB) skipped |
| **GSE193581** | Lu 2023 | scRNA + bulk cell line | **23 samples / 32 GSM** | Paper 1 Lu 2023 thyrocyte; SA4/SC2; cell line bulk count 다운로드됨 | ★★★ |
| **GSE241184** | Phase 1 single-cell | scRNA | **1 patient** | Paper 1 author-independence check; SA4 | ★★ |
| **cBioPortal `thca_tcga_pub`** | cBioPortal (TCGA THCA) | mutation + clinical | TERT promoter n=36 recovered | Paper 1 TERT⁺ Cox HR=4.33 p=4.9e-6 (memory `v17_tert_recovery_v2`) | ★★★ |
| **HM450 methylation TCGA** | TCGA Illumina HM450 | methylation array | **503** | Paper 1 promoter hypermethylation; manuscript abstract | ★★ |
| **GSE60542** | GEO PTC nodal mets | GPL570 microarray | **92** (17 PTC N+ / 11 PTC N0 / 17 LN-met / 24 paired-N / 4 normal LN) | Paper 3 ICI public download (already pulled); held from external_expression sweep | ★★ |
| **GSE78220 (Hugo)** | Hugo 2016 melanoma ICI | bulk RNA-seq + clinical | **28** | Paper 3 Track B-lite anchor (Hugo+Riaz pooled n=415) | ★★ — 이미 pooled 으로 사전에 있지만 개별 분리 명시 권장 |
| **GSE91061 (Riaz)** | Riaz 2017 nivolumab | bulk RNA-seq + WES + clinical | **109** | Paper 3 Track B-lite anchor | ★★ |
| **GSE151179** | post-RAI thyroid | gene expression | **52** | Paper 3 RAI-dediff axis; thyroid_diff d=−1.01 p=1e-4 | ★★ — 사전 v2에 ID 있지만 details 부족 |
| **GSE27155** | thyroid benchmark series | array | (verify) | manuscript / reports에 2,348회 등장 — 용도 verify 필요 | ★ verify needed |
| **GSE97466** | (verify) | — | — | 518회 등장 | ★ verify needed |
| **GSE39582 / GSE31210** | colon (39582), lung (31210) | array | — | pan-cancer comparison? Paper 11/12 검토 필요 | ★ verify |

**Action:** 위 1행~7행은 즉시 사전에 추가. 8~14행은 사용처 verify 후 추가.

---

## 2. 추가로 가져올 가치 있는 public 데이터 (분석 명령 아님, registry 만)

| Source | 가져오기 방법 | 예상 사이즈 | 가치 / 용도 |
|---|---|---|---|
| **Pozdeyev 2018 CCR (n=779 advanced DTC+ATC, 4-cluster ATC)** | publication supplementary table; raw 일부 dbGaP | supp <1 MB | Paper 1 Discussion 3.1 Landa 보강; advanced cohort 광범위 |
| **GSE60542 (PTC nodal mets)** | GEO Series Matrix | ~27 MB | external sweep에서 hold 였지만 lymph-node confounding 분석시 가치 — Paper 1 supp |
| **GSE151179 (post-RAI)** | GEO | ~7.5 MB (이미 다운로드됨) | Paper 3 RAI-dediff axis 추가 분석 |
| **PRJEB23709 (Gide)** | ENA | medium | Paper 3 Track B melanoma ICI 추가 |
| **PRJEB25780 (Kim GC)** | ENA | medium | Paper 3 pan-cancer Tier 3 |
| **GSE184362 RAW.tar (Pu 2021 scRNA full)** | GEO Supp | 970 MB | Paper 1 single-cell deep dive; Track B Wk5 |
| **GSE193581 RAW.tar (Lu 2023 scRNA full)** | GEO Supp | medium | Paper 1 single-cell deep dive |
| **GSE232237 / GSE191288 / GSE148673 RAW** | GEO Supp | medium | scRNA atlas integration |
| **TCGA-THCA processed matrices** | TCGAbiolinks / recount3 R 패키지 | medium | 재처리 검증; 현재는 Paper 1 ETL output 재사용 |
| **AFND South Korea full pool** | http://www.allelefrequencies.net | <10 MB | Paper 2/4 Korean baseline 정밀화 |

---

## 3. Controlled access / application 필요 데이터 (★★★ 핵심 차단 항목)

### 3.1 dbGaP

| Accession | PI / cohort | Modality | n | 사용 목적 | Application 시간 | 결정 |
|---|---|---|---|---|---|---|
| **TCGA-THCA paired tumor-normal BAM** | NCI TCGA | WES BAM | ~496 | LOHHLA / HLA LOH / neoantigen presentability (Paper 3 Module C) | **4–8 weeks** | 사용자 결정 대기 |
| **phs000452** | Liu et al. melanoma ICI | WES + RNA-seq | — | DIAL audit Tier 2 anchor (Paper 3) | 4–8 weeks | 사용자 결정 대기 |
| **phs000178** (TCGA-CDR) | TCGA Clinical Data Resource | clinical | — | survival downstream | 표준 application | low priority |

**현재 상태:** 사용자 명시 결정 전까지 신청 안 함 (`paper3_ici_data_access_blockers.md` §6).

### 3.2 EGA (European Genome-phenome Archive)

| Accession | Source | 추정 사용 | Application 시간 | 비고 |
|---|---|---|---|---|
| **EGAS00001003540** | (manuscript에 등장) | TBD verify | EGA standard | DAC (Data Access Committee) approval 필요 |
| **EGAS00001002556** | (manuscript에 등장) | TBD verify | EGA standard | DAC approval |
| **EGAD00001004845** | (manuscript에 등장) | TBD verify | EGA standard | dataset-level access |
| **Yoo SK 2019 Korean ATC** | Korean ATC RNA-seq | Asian dedifferentiation generalization | likely 4–8 weeks if controlled | EGA 또는 Korean repository — accession verify 필요 |

**Action:** EGAS / EGAD accession 3건 사용 목적 verify (manuscript에서 어디 cite 되는지) → 사용자 application 결정.

### 3.3 KoGES / KARE — Korean public-but-controlled

| Cohort | Authority | Modality | 사용 | Application |
|---|---|---|---|---|
| **KARE/KoGES** | KCDC Korean Genome Project | ~10K SNP genotype, ImmunoChip-like | Korean GD HLA imputation reference (Paper 4) | 한국 Korean Genome 신청 절차 — 영문 폼 + IRB |
| **Chen 2018 Han Chinese GD HLA fine-mapping** | PMC6161647 supplementary 또는 dbGaP | summary stats | Paper 4 Pan-Asian | supp PDF 추출 또는 dbGaP application |

---

## 4. PI / pharma request 필요 — **영원히 비공개일 수 있음**

이 데이터는 publication 만 있고 raw 가 publicly 비공개. **K1 (response predictor 클레임 금지)** 의 단일 근본 원인. 비공개 기관 코호트와 협력만이 K1 해제 경로.

| Trial / Paper | NCT | Sponsor / PI | 데이터 종류 | 현재 사용 | Request 경로 |
|---|---|---|---|---|---|
| **KEYNOTE-158 thyroid** | NCT02628067 | Merck | clinical + IHC; raw RNA-seq 비공개 | citation only | Merck investigator-initiated submission |
| **KEYNOTE-028 thyroid** | NCT02054806 | Merck | 동일 | citation only | Merck IIS |
| **Spartalizumab ATC** | NCT02404441 | Novartis | 동일 | citation only | Novartis ETOP |
| **NCT03246958 nivo+ipi aggressive thyroid** | NCT03246958 | — | clinical | citation only | trial PI |
| **NCT03181100 atezolizumab matched TT ATC** | NCT03181100 | Roche/Genentech | 동일 | citation only | Genentech IIS |
| **Dierks 2021 lenva+pembro** | — | publication only | 임상 데이터 | citation only | corresponding author |
| **Sehgal lenva+pembro** | various | publication only | clinical | citation only | corresponding author |
| **Cabanillas thyroid ICI** | various | publication only | clinical | citation only | corresponding author |
| **Pozdeyev 2018 CCR (n=779)** | — | Pozdeyev / Schweppe | mutation + clinical | supp table only | corresponding author for raw |

**현재 사용:** 모두 *clinical evidence cite* 만 가능. response prediction claim 차단됨 (K1).

---

## 5. Korean institutional outreach (Paper 4 Pan-Asian + Paper 2 Korean retrospective FFPE)

| 기관 / 코호트 | 데이터 종류 | 현재 상태 | 차단 해제 경로 |
|---|---|---|---|
| **Bundang SNUH Graves' BTC** | Korean GD + thyroid Ca | outreach 단계, 0% | Yu 교수 → Bundang outreach 응답 (6+ wk wait) |
| **연세 의료원 thyroid** | retrospective FFPE PTC ± Hashimoto | 미접촉 | Yu 교수 협력 |
| **서울대 본원 thyroid** | 동일 | 미접촉 | Yu 교수 협력 |
| **삼성 의료원 thyroid** | 동일 | 미접촉 | Yu 교수 협력 |
| **Bundang prospective Korean cohort** | 새 prospective FFPE | outreach 단계 0% | Bundang IRB |

**현재 결정 (memory v17_2026_04_30_pivot):** Bundang outreach **non-responsive 6+ weeks 또는 모든 access 거부** 시 Plan B trigger.

---

## 6. Population HLA — open이지만 country-specific aggregation 필요

| Source | 형태 | 사용 |
|---|---|---|
| **AFND South Korea pool** | Allele Freq Net DB online | Paper 2 Korean baseline (Chu 2018 GD comparator 대체) — 이미 사용 중 |
| **In 2015** | publication PDF | Korean baseline anchor |
| **Lee 2005 / Jung 2023 / Jekarl 2021 / Chung 2010** | publications | Korean baseline pack |
| **Shin 2019 / Cho 2011** | publications | Korean AITD/HD context |
| **Liao 2022** | publication | Pan-Asian GD registry (Paper 4 backlog) |
| **Park 2005 / Cho 1987** | publications | Pan-Asian GD registry |
| **Ueda 2014 / Chen 2011 / Naito 1987 / Tsai 1989 / Cavan 1994** | publications | Pan-Asian extended |
| **Stasiak 2023 systematic review** | publication | meta-context |

**현재 사용:** AFND South Korea pool 만 직접 통합. 나머지는 registry 단계 (`2026_05_04_hla_more_data_registry.md`).

---

## 7. 결정 행렬 (사용자 결정 대기 항목)

| 옵션 | 데이터 | 영향 | 비용 (시간) | 비용 (금전) |
|---|---|---|---|---|
| ☐ TCGA-THCA paired BAM dbGaP 신청 | LOHHLA Module C | Paper 3 Track B 진입 가능 | 4–8 weeks lag | 0 |
| ☐ Liu phs000452 dbGaP 신청 | DIAL Tier 2 melanoma | Paper 3 Track B 보강 | 4–8 weeks lag | 0 |
| ☐ Yoo SK 2019 EGA 신청 | Asian ATC anchor | Paper 1 generalization 강화 | 4–8 weeks lag | 0 |
| ☐ EGAS00001003540 / 002556 / EGAD00001004845 verify + 신청 | TBD | 사용처 verify 후 결정 | 4–8 weeks lag | 0 |
| ☐ Yu 교수 → 비공개 thyroid ICI 코호트 협상 | KEYNOTE-158/028 / Spartalizumab raw | **K1 해제** (response prediction claim 가능) | months | possibly MTA fee |
| ☐ Bundang SNUH outreach follow-up | Korean GD BTC | Paper 4 backlog 진입 | weeks | 0 |
| ☐ KoGES / KARE Korean Genome 신청 | Korean SNP, HLA imputation reference | Paper 4 Korean GD HLA strength 강화 | months | KoGES fee 가능 |
| ☐ GSE184362 RAW.tar pull (970 MB) | Pu 2021 full scRNA | Paper 1 SB/SC figures 본격화 | hours | 0 |
| ☐ GSE193581 RAW.tar pull | Lu 2023 full scRNA | 동일 | hours | 0 |
| ☐ Pozdeyev 2018 CCR supp 추출 | n=779 advanced DTC+ATC | Discussion 3.1 보강 | hours | 0 |
| ☐ Chen 2018 Han Chinese GD HLA supp 추출 | Pan-Asian fine-mapping | Paper 4 backlog 강화 | hours | 0 |

---

## 8. Sprint 권한 한계 (이번 audit 에서 *하지 않은* 것)

- ✅ 새로운 GEO 검색 안 함 (이미 manuscript / reports에서 등장하는 accession 만 inventory).
- ✅ 새로운 다운로드 안 함.
- ✅ dbGaP / EGA / KoGES application 안 함.
- ✅ Bundang outreach 메일 안 보냄.
- ✅ Pharma / PI request 안 보냄.
- ✅ Manuscript prose 안 만짐.
- ✅ Voice-protected section 안 적음.
- ✅ Paper 3 Track B 시작 안 함.
- ✅ Commit 안 함.

**registry only.** 의사결정은 사용자 권한.

---

## 9. 다음 권장 행보 (우선순위 순)

1. **사전 v3 갱신** — §1 의 ★★★ 7개 cohort 즉시 추가 (`_inject_dict.py` v3 + 47 page 재주입).
2. **Pozdeyev 2018 supp 추출** — public, hours, Paper 1 Discussion 3.1 강화. 사용자 GO 시 시작.
3. **EGAS / EGAD accession 사용처 verify** — manuscript grep 으로 확인 후 application 결정.
4. **Yu 교수 미팅에서 결정 받기:**
   - dbGaP TCGA paired-BAM 신청 여부
   - Yoo SK 2019 EGA 신청 여부  
   - Bundang outreach status
5. **Plan B trigger 조건 재확인** — Bundang non-responsive 6+ weeks 도달 시 Paper 4 backlog 정상화 / Plan B 진입.
