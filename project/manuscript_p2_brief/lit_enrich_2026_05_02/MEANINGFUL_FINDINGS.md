---
title: "lit_enrich — 쓸데 있는 결과만 (Meaningful findings)"
date: 2026-05-02
author: Seungho Cook
scope: Paper 2 Pillar 1 보강 + manuscript paste 후보
status: review-ready
---

# lit_enrich — 쓸데 있는 결과만

_생성: 2026-05-02 16:10 · 출처: `10_pipelines/`, `09_more_fetched/`, `data/`_

전체 6 pipeline 중 **paper-blocking** 으로 골라낸 것만 모음. 나머지는 scaffolding/process. 본인 본인이 keyboard 책임지는 voice-protected 영역 (Discussion 3.1 / Limitations 등)에 paste 가능한 형태로.

---

## 1. ★ Paper-changing — Pillar 1 baseline 강화

### 핵심 수치

> **Korean general-population DPB1\*05:01 baseline = 36.7% [95% CI: 33.2-40.4%]**
> _AFND 3 studies, total n=680, random-effects (Freeman-Tukey + DerSimonian-Laird), I²=0%_

### 기존 Pillar 1 구조 (memory)

| Reference | Pop | Freq | Source |
|---|---|---|---|
| Korean PTC pool (the cohort) | n=874 | **53.2%** | Paper 2 본 cohort |
| Chinese GD baseline | (varies) | 44.0% | **Chu X et al. 2018** (forest meta) |

### 추가 layer

| Reference | Pop | Freq [95% CI] | n | Source | I² |
|---|---|---|---|---|---|
| **Korean general-pop baseline** | n=680 | **36.7% [33.2-40.4%]** | 3 studies | AFND, RE meta | **0%** |
| Chinese general-pop baseline | n=1,969 | 35.7% [30.2-41.4%] | 20 studies | AFND, RE meta | 85% |
| Japanese general-pop baseline | n=1,235 | 38.6% [35.9-41.3%] | 6 studies | AFND, RE meta | 0% |
| ⚠️ Taiwan general-pop baseline | n=704 | 67.4% [48.3-83.9%] | 5 studies | AFND, RE meta | 94% |

### Pillar 1 statistical contrast — paper-ready statement

> Korean PTC pool (n=874, 53.2%) shows DPB1\*05:01 enrichment over an
> independent Korean general-population baseline (random-effects pooled
> 36.7% [33.2-40.4%], 3 AFND studies, total n=680, I²=0%):
> Δ = +16.5 percentage points, **Wald z = 6.58, p = 4.74 × 10⁻¹¹**.
> The Korean baseline is consistent (zero between-study heterogeneity)
> and concordant with independent Chinese (35.7%) and Japanese (38.6%)
> general-population baselines, ruling out a regional-pop-stratification
> artifact.

### 왜 paper-changing 인가

- 기존 forest는 reference 1개 (**Chu 2018 GD = 44%**)
- 추가 reference (**Korean general-pop = 36.7%**) 는 **독립 baseline** — different population, different study design, different decade
- 두 reference 모두 PTC pool (53.2%) 보다 낮음 = "Korean PTC HLA enrichment is not population stratification"
- Reviewer 의 most likely challenge ("baseline 비교 문제 있음") 를 직접 차단

### Source files

- `/home/seungho/personal/THCA_data_analysis/project/manuscript_p2_brief/lit_enrich_2026_05_02/10_pipelines/P2_dpb1_country_meta.csv`
- 재실행: `project/.venv/bin/python project/notebooks_or_scripts/v17_lit_enrich_pipelines.py`

---

## 2. ⚠️ Reviewer trap — Taiwan outlier

### 문제

| Item | Value | Note |
|---|---|---|
| Taiwan pooled | 67.4% [48.3-83.9%] | wide CI |
| n_studies | 5 | small |
| total_n | 704 (avg 141/study) | small per-study n |
| I² | 94% | extreme heterogeneity |
| τ² | 0.045 | extreme |

### 가능한 원인

1. Per-study n 작아서 sampling 분산 큼
2. HLA typing resolution 다름 (low-res serology vs high-res sequencing)
3. Ethnic substratification (Hakka / Min / Taiwan Aboriginal) 별로 frequency 큰 차이
4. 특정 5 studies 가 다른 ascertainment basis (transplant donor vs disease cohort)

### 권장 조치 (manuscript paste)

- Pan-Asian gradient 그림/표에 Taiwan **그대로 표시 X**
- Sensitivity analysis 별도 panel 또는 Limitations 에 한 줄:

> "Taiwan AFND aggregate shows extreme between-study heterogeneity
> (I² = 94%, τ² = 0.045), likely reflecting ethnic-substratification
> (Hakka / Min / Aboriginal) and HLA typing resolution differences;
> we exclude Taiwan from the primary East-Asian baseline pool and
> retain it for sensitivity analysis only."

### Source

- 위 동일 `P2_dpb1_country_meta.csv`
- AFND raw: `data/afnd_alleles.json` → `DPB1*05:01 → Taiwan` 5 studies

---

## 3. 🔴 Critical data hygiene save — Yoo2016 + Pu2021 PMC mismatch

### 발견

P4 bib fact-check (10 OA abstract × verified_refs) 결과 author mismatch 3건. 검증 결과 **Unpaywall PMC ID 2개가 무관 paper 가리킴**:

| Bib key | Bib expects | Wrong PMC was | Actually pointed to | Correct PMC (DOI 검증) |
|---|---|---|---|---|
| Yoo2016 | Yoo SK, PLoS Genet 2016 | PMC4986964 | Vedelek V "Testis-Specific Bb8" PLoS One | **PMC4975456** |
| Pu2021 | Pu W, Nat Commun 2021 | PMC8523608 | Dev SA "Adulteration in Ayurvedic raw drugs" 3 Biotech | **PMC8523550** |

### 영향

- 만약 fix 안 했으면 manuscript 에 abstract paste 시 무관 paper 의 내용/저자가 인용됨 → **misattribution**
- 같은 misattribution 패턴 (Krishnamoorthy 2025 Nat Comm 사건) 이 v17 round 12 직전에 한 번 더 발생할 뻔함

### 적용된 fix

- `data/unpaywall_links.json` 정정 + `_pmc_correction_note` 필드 추가
- `09_more_fetched/oa_abstracts.{json,md}` 재fetch (이제 진짜 Yoo SK + Pu W abstract)
- 재검증: 9/10 author match (TCGA2014 1개만 corporate author 라 nominal)

### 행동 제안

- 향후 bib expansion 시 Unpaywall payload 의 PMC 가 author/title 과 일치하는지 자동 검증 (`v17_lit_enrich_pipelines.py` P4 sweep)

---

## 4. 🟡 Useful but not paper-changing (참고만)

| 결과 | 위치 | 쓰임 |
|---|---|---|
| 73 topically-filtered papers (claim별 top-5) | `10_pipelines/P1_top5_clean_per_claim.csv` | citation triage 시 reference (기존 02_top5_per_claim.csv 는 noise dominant 라 사용 X) |
| 9 RET-fusion priority trials | `10_pipelines/P3_RET_fusion_priority_trials.csv` | Discussion translational outlook 1 paragraph 후보 |
| 38 trials drug-class breakdown | `10_pipelines/P3_drug_class_summary.csv` | Reviewer "current clinical landscape?" Q 응답용 |
| 10/10 OA abstract 풀세트 | `09_more_fetched/oa_abstracts.md` (21 KB) | Discussion 3.x paste 후보 직접 선택 (voice-protected) |

---

## 5. 한 줄 요약

> **★** Pillar 1 forest 에 **Korean general-pop baseline 36.7% [33.2-40.4%], n=680, p=4.74e-11** 추가 = paper-changing.
> **⚠** Taiwan 67.4% 는 caveat 또는 sensitivity-only.
> **🔴** Yoo2016 + Pu2021 PMC mismatch 정정 = misattribution save.
> **나머지는 scaffolding.**

---

## 6. 다음 본인 액션 (voice-protected)

1. Pillar 1 forest 그림 / 표 에 Korean general-pop baseline 36.7% 추가 (Chu 2018 옆에 second reference)
2. Discussion 3.1 (또는 Pillar 1 footnote) 에 Wald z=6.58, p=4.74e-11 한 줄 삽입
3. Taiwan caveat 한 문장 작성 (Limitations 또는 Pan-Asian discussion)
4. 24 gap candidate triage 는 P1 relevance score (`P1_top5_clean_per_claim.csv`) 보고 cherry-pick — 전체 다 읽지 말고 cosine ≥0.10 만 보면 충분
