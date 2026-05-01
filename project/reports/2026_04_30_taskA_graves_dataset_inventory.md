# Task A — Graves' / Autoimmune Thyroid Open Dataset Inventory

**Date:** 2026-04-30 (Day 1 PM)
**Status:** Initial sweep complete (4 hr)

---

## ★ TL;DR (1 줄)

**한국 PTC+Hashimoto overlap RNA-seq (GSE286332, n=18) + Han Chinese GD (n=2,958 + 1,490 HLA imputation panel) 가 Top 2.** 분당 cohort 에 Graves' sample 있으면 **STRONG GO 시나리오 2 (multi-omic Graves' subtype paper)** — Cell Rep Med / JCI Insight reach.

---

## 1. Top dataset candidates (15 sweep, top 5 deep dive)

| Source | Accession | n (case/ctrl) | Modality | Population | Accessibility | HLA-typeable | ★ Tier |
|---|---|---|---|---|---|---|---|
| **GEO** | **GSE286332** ★★ | **9 PTC + 9 PTC+HT** | RNA-seq (Illumina NovaSeq X) | **Korean (Dongguk Univ)** | Open + SRA PRJNA1208932 | ✅ arcasHLA | ★★★ |
| Nat Genet GWAS | Cooper 2012 | ~3,000 / 7,500 | ImmunoChip | EUR | Controlled | ✅ direct | ★★★ |
| **Han Chinese fine-map** | Chen 2018 | **1,468 GD / 1,490 ctrl** | SNP array + HLA imputation | Han Chinese | Open summary | ✅ pan-Asian panel | ★★★ |
| Taiwan CMUH | 2024 (PMC8936090) | **2,998 GD / 29,083 ctrl** | EMR + HIBAG SNP imputation | Taiwanese | Controlled access | ✅ done | ★★ |
| GEO | GSE71956 | ~10 / ~10 | Microarray | EUR | Open | ⚠ small | ★ |
| GEO | GSE29315 | 6 HT / 8 TPH | Microarray | EUR | Open | ⚠ small | ○ |
| GEO | GSE138198 | 13 HT / 3 ctrl | Microarray | EUR | Open | ⚠ small | ○ |
| GEO | GSE308553 | TED orbital adipose | RNA-seq | EUR | Open | × tissue-specific | ○ |
| GEO | GSE105149 / GSE58331 | TED lacrimal | Microarray | EUR | Open | × | ○ |
| Japanese | Okada 2015 | ~2,000 GD | HLA imputation | Japanese | Controlled | ✅ Pan-Asian ref | ★★ |
| Korean lit | Park 2005 | ~200 | HLA-DR/DQ typing | Korean | Closed | ✅ done | ★ |
| Han 2010 | KJIM | ~5 / ~5 | Microarray | Korean | Closed | ✅ small | ○ |
| KARE/KoGES | TBD | ~10K SNP genotype | ImmunoChip-like | Korean | dbGaP / KoGES | ✅ if access | ★★ |
| 분당 cohort | TBD outreach | TBD | TBD | Korean | Pending | ✅ if SNP/RNAseq | ★★★ if Graves' |

---

## 2. Top 3 deep dive

### ★★★ #1 — GSE286332 (Korean PTC+Hashimoto RNA-seq, Macrogen Seoul)
- **Why critical:** Direct molecular phenotype of Korean PTC/Hashimoto overlap — exactly the autoimmune-PTC sub-axis hypothesis we've been chasing
- **Method fit:** RNA-seq → arcasHLA imputation possible (already validated pipeline)
- **Limitation:** Only n=18 (9+9) — too small for HLA association alone
- **Use case:** **Validation cohort** for Hashimoto-like signature (Q12 GSE213647 17% finding) + arcasHLA HLA imputation cross-cohort meta
- **Action:** 즉시 download (PRJNA1208932 SRA), 5M-read partial → arcasHLA → integrate with 우리 K2+Lee meta (n=865+18=883)

### ★★★ #2 — Han Chinese GD HLA fine-mapping (Chen 2018, n=1,468 + 1,490)
- **Why critical:** Largest published Asian Graves' HLA association study with public summary statistics
- **Findings overlap with우리:** HLA-DPB1*05:01 + B*46:01 둘 다 confirmed Asian Graves' risk
- **Method fit:** Pan-Asian HLA reference panel applicable to Korean
- **Limitation:** Han Chinese ≠ Korean (subtle allele frequency 차이)
- **Use case:** **External validation reference** — paper Discussion 의 "Pan-Asian replication" 직접 비교
- **Action:** Summary statistics download → 우리 K2+Lee 결과 와 forest plot 메타분석

### ★★ #3 — Taiwan CMUH (n=2,998 + 29,083)
- **Why critical:** Massive sample size, EMR-linked, HIBAG imputation already done
- **Method fit:** SNP-based HLA imputation comparable to arcasHLA at population level
- **Limitation:** Controlled access (need formal request), Taiwanese ≠ Korean
- **Use case:** **Largest reference for HLA-disease association meta** if access granted
- **Action:** Defer to Bundang prospective cohort approval (parallel)

---

## 3. 분당 cohort Graves' query (Day 1 PM action)

분당 outreach email 에 추가:
> "분당서울대 cohort 에 Graves' / Hashimoto / autoimmune thyroid sample 이 있는지 — 있으면 (1) modality (microarray / RNA-seq / WES), (2) sample n, (3) HLA typing/imputation 가능 수준, (4) 임상 metadata (TRAb, anti-TPO, anti-Tg titer, 치료력)"

---

## 4. Paper feasibility — 3 시나리오

| 시나리오 | 조건 | Timeline | Target Venue |
|---|---|---|---|
| **시나리오 1 — Pure HLA association** | GSE286332 + Han Chinese summary stats | 4-6 mo | Front Immunol / Clin Exp Immunol (IF 5-7); reach J Autoimmun (IF 12-14) |
| **★ 시나리오 2 — Multi-omic Graves' subtype** | + 분당 Graves' cohort 받음 | 6-9 mo | J Clin Endocrinol Metab / Thyroid (IF 5-7); reach **Nat Commun** (본인 메모리상 비슷한 angle 게재 경험) |
| 시나리오 3 — Method paper (cookHLA application) | n < 50 만 가능 시 | 3-4 mo | Brief Bioinform / Bioinformatics (IF 4-7) |

---

## 5. Decision tree (EOD Day 2)

```
분당 cohort 에 Graves' sample 있나?
├── YES (n > 50): 시나리오 2 STRONG GO → Western discovery (GSE286332 + Han) + Korean validation (분당)
├── YES (n < 50): 시나리오 1 + 분당 부분 cite → 6 mo
├── NO (Graves' 없음): 시나리오 1 (Western only) → 4-6 mo, IF 5-7 venue
└── 분당 outreach 6+ wk 응답 없음: 시나리오 3 (method paper) backup
```

---

## 6. 본인 unique angle (5-7 년 plan 정렬)

- 본인 = autoimmune (cookHLA Nat Comm 류마티스+T1D+CD) + thyroid (현 PTC paper 진행 중) **두 분야 교차**
- 두 분 (주영석/최정균) 못 가는 영역
- 한국 medical AI 첫 mover identity 와 정렬 — autoimmune × thyroid 횡단

---

## 7. Sources

Sources:
- [Analysis of HLA Variants and Graves' Disease (Frontiers)](https://www.frontiersin.org/journals/endocrinology/articles/10.3389/fendo.2022.842673/full)
- [Construction of HLA imputation reference (Nat Genet)](https://www.nature.com/articles/ng.3310)
- [Han Chinese GD HLA fine-mapping (PMC6161647)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6161647/)
- [Korean GD HLA-DR/-DQ (PMID 15993720)](https://pubmed.ncbi.nlm.nih.gov/15993720/)
- [GSE286332 Korean PTC+HT RNA-seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286332)
- [Genome-Wide Association Studies of AITD (E-ENM)](https://e-enm.org/upload/pdf/enm-33-175.pdf)
- [Korean cohort Hashimoto vitamin D (Front Endocrinol 2025)](https://www.frontiersin.org/journals/endocrinology/articles/10.3389/fendo.2025.1666115/full)
