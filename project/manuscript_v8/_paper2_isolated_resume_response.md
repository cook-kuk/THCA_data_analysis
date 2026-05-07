---
title: "Paper 2 isolated session — resume response (deliverable options A-F)"
date: 2026-05-04
purpose: Truncation 정정 후 Paper 2 isolated scope 재확정 + deliverable 후보 6개 + marathon safety score + 본인 선택 요청
status: Spec/plan only — edit/scan 실행 0. 본인 선택 받기 전 까지 대기.
related_files:
  - _session_origin_paper2_isolated.md (origin context)
  - _audit_HT_vs_GD_2026_05_02.md (HT vs GD conflation audit)
  - memory/v18_paper2_HT_isolated.md (Yu professor 5/4 결정)
  - memory/paper_numbering_2026_05_04.md (canonical numbering)
---

# Paper 2 isolated session — resume after truncation

## [1] Truncation acknowledged

이전 5/4 session 활성 메시지가 "(2) Cancer-only terms — DM" 에서 cut off. 본 resume message 로 scope + forbidden list + deliverable 후보 모두 명시 받음. 추측 없이 confirm 된 instruction 만 따름.

---

## [2] Paper 2 isolated scope

| 항목 | 내용 |
|---|---|
| **Topic** | Hashimoto-overlap PTC 의 분자 축 — antigen-driven HLA-II–mediated dedifferentiation |
| **Substrate** | Hashimoto-overlap PTC ONLY (GSE286332 + Korean PTC HLA pool n=874 + GSE213647 Hashimoto-like 22-28%) |
| **Mechanism focus** | HT 의 destructive lymphocytic infiltration → thyroid-lineage dedifferentiation + RAI loss; HLA-II antigen presentation; mediation model |
| **GD reference 허용 한도** | "shared genetic background" 정도 1-line; disease equivalence 주장 X; mechanism 직접 사용 X |
| **Paper 1/3/4 boundary** | DM1/DM2 cluster mechanism deep dive X · ICI/HLA LOH/neoantigen X · GD/TSAb/hyperthyroid X |

---

## [3] Forbidden vs allowed language

### Forbidden (3 categories, ~30 terms)

**Paper 4 territory (GD)**:
- Graves, GD, TSAb, TSI, thyroid-stimulating antibody
- hyperthyroidism, thyrotoxicosis, exophthalmos
- Graves' ophthalmopathy, TED, thyroid eye disease
- Bundang Graves cohort

**Paper 1 territory**:
- DM1/DM2 cluster mechanism deep dive
- BRAF/RAS driver mechanism deep dive
- TERT promoter kinetics
- kinase fusion mechanism
- dark matter subtype mechanism
- TROP2 / sacituzumab
- TCGA WSI / H&E-DM1
- image-DM1

**Paper 3 territory**:
- ICI response prediction
- HLA LOH
- neoantigen
- DIAL audit
- pan-cancer ICI
- C5AR1 anti-PD-1 synergy

### Allowed Paper 2 vocabulary

- Hashimoto-overlap PTC
- destructive lymphocytic thyroiditis
- antigen presentation
- HLA-II
- thyroid-lineage dedifferentiation
- RAI/thyroid differentiation loss
- BCR/TLS (overexpansion 금지 — 이미 Paper 2 scope 에 포함된 부분만)
- mediation model
- Korean PTC vs Korean baseline
- AFND South Korea baseline
- HT-specific tissue context

---

## [4] Deliverable options A–F

| ID | Deliverable | Output type |
|---|---|---|
| A | Paper 2 HT-isolated outline skeleton | heading + bullet skeleton, voice-protected 표시, prose X |
| B | Pillar I v2 forest meta spec | input/output/schema/statistical plan only, 분석 X |
| C | Terminology correction plan | sed-level edit plan, 실행 X |
| D | Paper 2 5-pillar 재정의 | HT-only substrate 기준 5-pillar 구조 |
| E | Cross-contamination audit plan | Paper 1/2/3/4 forbidden terms scan plan, 실행 전 승인 |
| F | Advisor-facing one-page brief | structured bullets only, prose X |

---

## [5] Marathon safety score

| ID | Score | Reasoning |
|---|---|---|
| A | ★★ | Heading+bullet 만이면 safe, 그러나 skeleton 의 bullet 이 prose 로 expand 할 risk (특히 voice-protected sections — Hook, Aim, mechanism story). 본인 키보드 영역 명시 + bullet 길이 cap 필요 |
| B | ★★★ | Spec only (input/output/schema). 분석 실행 없음. 본인 5/2 audit 의 Conflation 1 정정 직결 (Korean PTC vs AFND South Korea baseline). Yu professor 미팅 confirmation 후 실행 가능 |
| C | ★★★ | sed-level plan, dry-run output 만. 실행 전 본인 승인 절차 포함. manuscript_v8 + memory 양쪽 file list + 정정 전후 diff preview 가능 |
| D | ★★ | 5-pillar 구조 변경은 paper 정체성 결정 — 본인 voice 영역. 5-pillar 의 categorical structure 만 plan 으로 가능, 각 pillar 의 narrative 는 본인 |
| E | ★★★ | Audit plan 작성 → 본인 승인 후 scan 실행. forbidden terms grep + 위치 list 만, 자동 정정 X |
| F | ★★ | Yu professor 외부 share 용 — 본인 voice 가 advisor 에게 전달되는 것 → 본인이 직접 본인 voice 로 작성하는 게 진짜 권장. Bullets 만이면 OK 하지만 advisor 가 "draft 톤" 으로 인식할 risk |

**범례**:
- ★★★ = full marathon-safe (analysis-free, plan-only, no voice creep)
- ★★ = marathon-safe with discipline (skeleton/structure 만 유지 시 OK, prose creep risk)

---

## [6] Recommended next choice

본인 default 추천 (**B + C**) 에 동의. 근거 정렬:

- **B 먼저** — Conflation 1 (Korean PTC vs Chu 2018 GD direct comparison) 의 정정 spec. AFND South Korea baseline 결정이 Paper 2 의 Pillar I framing 전체 결정. Spec lock 후 다른 deliverable 가 그 위에 build 가능.
- **C 동시** — Terminology contamination 차단이 outline/prose 작성 전에 base layer. "autoimmune-PTC" + "GD" mentions 가 manuscript_v8 + memory 에 잔존 (이전 sprint 잔여물). Scan plan + diff preview → 본인 승인 → 실행.

추가 제안 (선택):
- **E 보조 흡수** — C 의 scan plan 의 일부로 E 통합 가능. Forbidden words grep 이 C 의 "정정 대상 list" 와 동일 작업. C 안에 E 의 audit plan 흡수 가능.

순서: **B (spec) → C (scan plan + 본인 승인 후 실행) → 그 다음 다른 deliverable**.

---

## [7] 본인 선택 요청

다음 중 1개 (또는 조합) 선택하시면 그것만 작성:

- **(i)** B + C 동시 진행 (default 권장)
- **(ii)** B 먼저, C 는 다음 round
- **(iii)** C 먼저, B 는 다음 round (terminology 깨끗한 base 위에 spec 작성하고 싶으시면)
- **(iv)** B + C + E 통합 (E 가 C 안에 자연스럽게 흡수)
- **(v)** 다른 조합 (D / F / A 추가)

**본인 선택 받기 전 까지 — edit 또는 scan 실행 0. Spec/plan 작성도 본인 confirm 후 시작.**

---

## Marathon mode reaffirmation

- 5/4–6/13 active
- 새 분석 실행 금지
- 새 데이터 다운로드 금지
- 5,000w sprint generation 금지
- voice-protected prose generation 금지
- Claude = scaffolding/audit/spec only
