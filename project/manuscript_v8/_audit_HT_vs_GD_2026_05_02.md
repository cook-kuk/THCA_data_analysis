---
title: "★ Critical audit — HT (Hashimoto) vs GD (Graves) conflation + 3-topic isolation"
date: 2026-05-02
trigger: User catch — "그레이브스병과 하시모토갑상선염을 혼용해서 사용한다는 것"
scope: 거기까지만 — audit + 3-topic spec + remaining analysis 명시. 분석 실행 X.
status: ★★★ Paper 2 framing 의 fundamental error 발견. Paper 1 은 영향 없음.
---

# 1. 의학적 정확 정의 (사용자 제공, 본인 검증)

| 항목 | Hashimoto's Thyroiditis (HT) | Graves' Disease (GD) |
|---|---|---|
| **기전** | 만성 항원 노출 + 림프구 침윤 → 갑상선 조직 파괴 | TSH receptor 자극 자가항체 (TSI) → 갑상선 호르몬 과다 분비 |
| **결과** | Hypothyroidism (기능 저하) | Hyperthyroidism (기능 항진) |
| **조직학** | 림프구 침윤 + TLS (3차 림프 구조) 형성 + thyrocyte 파괴 | 미만성 갑상선 비대 + follicle hyperplasia |
| **PTC 와의 동반율 (Korean)** | 20-30% (Korean PTC + HT) | 매우 드묾 (PTC 와 GD 동반은 임상적으로 흔치 않음) |
| **HLA 위험 allele (East Asian)** | DPB1\*05:01 등 (HT-specific 연구 빈약) | DPB1\*05:01 OR=1.90 (Chu 2018, well-published) |
| **본 project 데이터 substrate** | GSE286332 (n=18, PTC vs PTC+HT, Korean) | 별도 cohort 없음 — Chu 2018 published summary stats only |

**핵심 차이**: 두 질환은 같은 thyroid autoimmunity 우산 아래 있지만 **서로 다른 유전 병태생리**. HT-specific HLA 연구는 GD 보다 적음.

---

# 2. ★ Conflation 발견 — 3 핵심 location

## Conflation 1 (★★★ critical — Paper 2 framing fundamental error)

**위치**: `v17_paper2_pillar1_forest_strong.md` line 26
> "Korean PTC pool DPB1\*05:01 (53.2%) **HIGHER than Chu Han Chinese GD itself (44.0%)** — Korean PTC has more pronounced autoimmune-thyroid risk allele profile than Chinese GD."

**문제점**:
- Korean PTC pool (n=874) 의 autoimmune comorbidity 는 **HT (20-30%)**, GD 가 아님
- Chu 2018 cohort 는 **GD (n=1,468)**, HT 가 아님
- 직접 비교 (Korean PTC HLA vs Chinese GD HLA) 는 **HT cohort vs GD cohort** 의 비교 → 구조적 mismatch
- 현재 framing: "Korean PTC HLA pattern shares GD susceptibility" → reviewer 가 "왜 HT cohort 를 GD 와 비교?" 잡힘

**올바른 framing**:
- DPB1\*05:01 은 **East Asian thyroid autoimmunity 일반 risk allele** (HT + GD 양쪽 모두 share)
- Korean PTC 53.2% 는 일반 Korean baseline (~38-42%, Lee 2014) 위 — Korean PTC 의 HT comorbidity 가 이 elevation 의 source 일 가능성 높음
- Chu 2018 GD 44% 는 Chinese GD 의 elevation, Chinese baseline 38% 위
- 직접 "Korean PTC > Chinese GD" 비교 ❌ — 두 cohort 가 다른 disease entity
- 올바른 비교: Korean PTC (HT-enriched) vs **Korean general population** (Lee 2014) + Chinese GD (Chu 2018) vs **Chinese general population** (Chu controls) — 각 cohort 의 within-population elevation 비교

## Conflation 2 (★★ medium — Paper trajectory contradiction)

**위치 1**: `v17_2026_04_30_pivot.md` line 9
> "Task A — Graves'/autoimmune thyroid open dataset 탐색"
> "Graves' paper trajectory 는 cancer paper 와 별개 별도 paper"

**위치 2**: `v17_graves_pivot.md`
> "Graves' disease autoimmune molecular phenotype paper" (별도 paper trajectory)

**위치 3**: `v17_paper2_pillar1_forest_strong.md` (위 Conflation 1) — Paper 2 에서 GD cohort 사용

**문제점**:
- 4/29 + 4/30 결정: **GD = 별도 paper** (J Autoimmun / Front Immunol target)
- 5/3 결정: **Paper 2 (PTC + HT)** 의 Pillar 1 forest 에 GD cohort 사용 → 4/29-4/30 결정 위반
- Paper 2 의 substrate 는 HT, GD 가 아님 → GD comparison 은 Paper 3 (별도) 으로 이동해야

## Conflation 3 (★ minor — terminology slippage)

**위치**: 다수 memory 파일 + manuscript_v8 파일 — "autoimmune-PTC" 라는 generic term 사용

**문제점**:
- "autoimmune-PTC" 가 HT-PTC 와 GD-PTC 양쪽 다 포괄하는 듯 모호
- 본 project 의 실제 데이터: GSE286332 = PTC+HT, 즉 "Hashimoto-PTC" 또는 "PTC with concurrent HT"
- "autoimmune-PTC" 사용 시 reviewer 가 GD-PTC 까지 포함하는지 의심

**올바른 표현**:
- ✅ "PTC with concurrent Hashimoto's thyroiditis (PTC+HT)"
- ✅ "Hashimoto-overlap PTC"
- ❌ "autoimmune-PTC" (모호)
- ❌ "autoimmune-thyroid-PTC" (HT 와 GD 짬뽕)

---

# 3. ★ 3-Topic 격리 spec (strict)

## Topic 1 — Paper 1 (Cell Reports Medicine, 2026 6월 bioRxiv)

| 항목 | 내용 |
|---|---|
| **Title** | "An 8-gene panel reveals fusion-driven, epigenetically silenced dark matter in BRAF/RAS-negative thyroid cancer" |
| **Substrate** | TCGA-THCA n=504 + MSK n=117 + Korean PTC n=874 (with HT comorbidity 20-30% baseline) + sc external |
| **Findings** | 8-gene DM1/DM2 + DM1 76.8% fusion+ + DM1 epigenetic hypermethylation + meta HR 2.53 |
| **HT 언급 정도** | minimal — Section 2.4b sub-B 의 "older + immune-hot, possibly Hashimoto-overlap" teaser 만. **상세 mechanism 다루지 않음** |
| **GD 언급** | ❌ **완전 제외** |
| **Paper 1 영향** | ✅ **Paper 1 은 conflation 영향 없음** — 이미 sub-B = "Paper 2 reserve" 로 분리해 둠 |

## Topic 2 — Paper 2 (Korean Hashimoto-PTC paper, J Autoimmun / JCI Insight target)

| 항목 | 내용 |
|---|---|
| **Working Title** | "HLA susceptibility and B cell mechanism in Korean papillary thyroid carcinoma with concurrent Hashimoto's thyroiditis" (★ HT 명시) |
| **Substrate** | GSE286332 (n=9 PTC + n=9 PTC+HT, Korean) + Korean PTC HLA pool n=874 + GSE213647 Hashimoto-like 22-28% |
| **Findings** | HT-PTC 의 HLA-II d=+3.65, IGHV clonal d>0.5, TLS d=+1.96, AICDA up, mediation 140% (HLA-II), DM1 sub-B = NBNR cluster |
| **HLA framing (정정)** | "DPB1\*05:01 53.2% in Korean PTC pool exceeds Korean general population baseline (~38-42%, Lee 2014), consistent with thyroid autoimmunity comorbidity. **The allele is shared with Graves' disease (Chu 2018, Chinese GD 44% with OR=1.90 vs Chinese controls 38%) as a general East-Asian thyroid autoimmunity susceptibility marker, but our cohort substrate is HT not GD.**" |
| **GD 언급** | ⚠️ **Background reference only** (DPB1\*05:01 의 East-Asian thyroid autoimmunity 일반 위험 맥락 1줄). Forest meta 직접 비교 ❌ → Korean PTC vs Korean baseline + Chinese GD vs Chinese baseline (separate panels) |
| **Forest meta 정정** | 현재 "Korean PTC vs Chu GD direct" → "Within-population elevation: Korean PTC vs Korean baseline (Lee 2014); Chinese GD vs Chinese baseline (Chu 2018 controls)" 두 separate panel |
| **Paper 1 timeline 영향** | ❌ 없음 — Paper 2 별도 trajectory (Phase 1, 6-18 개월) |

## Topic 3 — Paper 3 (Graves' disease paper, J Autoimmun target, BACKBURNER)

| 항목 | 내용 |
|---|---|
| **Working Title** | "Graves' disease autoimmune molecular phenotype" (별도 paper trajectory, `v17_graves_pivot.md` 참조) |
| **Substrate** | TBD — Bundang Graves' BTC cohort (없으면 Chu 2018 published summary stats only) |
| **Status** | ★ **Backburner — Bundang outreach 결과 wait** |
| **HLA framing** | DPB1\*05:01 OR=1.90 in Chinese GD (Chu 2018) — 본 paper 의 main pillar |
| **HT 언급** | ❌ **완전 제외** (paper substrate 가 GD only) |

**3-topic 격리 strict rule**: Topic 1/2/3 substrate (cohort + disease) 정확히 다름. Cohort comparison 시 same-disease 만 직접 비교, cross-disease (HT vs GD) 는 background reference 만.

---

# 4. ★ 더 해야 할 분석 (limited, 거기까지만)

본인 catch 후 진짜 paper-blocking 분석:

## 4.1 Paper 2 forest meta 재구성 (★★★ critical)

**목표**: Conflation 1 정정.

**현재 (잘못된)**: Korean PTC vs Chu GD 직접 forest meta
**정정 필요**:
- (a) Korean PTC pool n=874 vs **Korean general population baseline** (Lee et al. Tissue Antigens 2014, KOTRY donor) — **same population, within-elevation**
- (b) Chu 2018 GD vs Chu 2018 Chinese controls (already published) — separate panel, background context
- (c) DPB1\*05:01 의 "East Asian thyroid autoimmunity allele" 일반 framing — HT + GD 양쪽 risk

**소요 시간**: 1-2 hr (Lee 2014 Korean baseline lookup + 새 forest plot 2-panel)

**Paper-blocking?**: Paper 2 는 yes (현재 framing 잘못). Paper 1 은 no.

## 4.2 GSE286332 의 HT 명시화 (★★ medium)

**목표**: Conflation 3 정정 — "autoimmune-PTC" → "PTC with concurrent HT" 일관 표현.

**작업**: manuscript_v8 + memory 파일들에서 "autoimmune-PTC" → "Hashimoto-overlap PTC (PTC+HT)" 일괄 정정. Paper 1 의 Section 2.4b + Discussion 3.3 mention.

**소요 시간**: 30 min (sed 일괄 + 본인 verify)

## 4.3 Memory 정리 (★ minor)

**목표**: Conflation 2 정정 — 4/29-4/30 결정 (GD = 별도 paper) 와 5/3 결정 (Paper 2 에 GD 사용) 모순 해결.

**작업**: 
- `v17_paper2_pillar1_forest_strong.md` 의 framing 정정 (HT-substrate 명시 + Chu GD 는 background reference only)
- `v17_graves_pivot.md` cross-reference 추가 (Paper 3 = backburner)

**소요 시간**: 15 min

## 4.4 Paper 2 outline 시작 (★ optional, marathon-compliant)

Paper 2 본격 outline 은 Paper 1 6월 bioRxiv 후 (Phase 1, 6-18 개월). 지금은 conflation 정정만 + Paper 2 의 substrate (HT-PTC) explicit 명시까지만.

---

# 5. Paper 1 영향 평가 (★ 좋은 소식)

| Paper 1 component | Conflation 영향 | 정정 필요 |
|---|---|---|
| Title | ❌ 영향 없음 (HT/GD 단어 없음) | 없음 |
| Abstract | ❌ 영향 없음 | 없음 |
| Section 2.1-2.3 (panel + fusion + epigenetic) | ❌ 영향 없음 | 없음 |
| Section 2.4b (immune-overlap teaser) | ⚠️ minor — "autoimmune-PTC" → "Hashimoto-overlap PTC" 표현 정정 | 1-line edit |
| Section 2.5 (cross-cohort) | ⚠️ minor — Hashimoto-like 22-28% 표현은 OK (이미 HT-specific) | 없음 |
| Discussion 3.3 (East-Asian) | ⚠️ minor — "autoimmune-PTC" 표현 정정 | 1-line edit |
| Methods | ❌ 영향 없음 | 없음 |
| Cover letter | ❌ 영향 없음 | 없음 |

**결론**: Paper 1 은 **conflation 영향 거의 없음** (sub-B teaser + 1-2 line 표현 정정만). 본인이 4/29 부터 "Pillar 5 main story = Paper 2 reserve" 분리를 정확히 해 둔 덕분.

**Paper 2 는 framing 재구성 필요** (Pillar 1 forest meta 가 잘못된 cross-disease comparison).

---

# 6. Paperclip (https://gxl.ai/blog/adding-arxiv-and-abstracts) 적용성

**Tool**: agent-native 색인 시스템 — 3M arXiv + 150M+ OpenAlex abstracts. `paperclip cat / search / grep / sql` 명령어로 full-text section structure 보존하며 검색.

**핵심 기능**:
- `paperclip search "..." -s abstracts` — 150M abstract hybrid 검색 (BM25 + vector)
- `paperclip cat arxiv_2501.12948 --section Methods` — section-aware full-text fetch
- `paperclip grep` — 3M papers full-text grep (80초 demo)

**본 project 적용성**:
- ★★★ **Paper 2 forest meta 재구성**: Korean general population DPB1\*05:01 baseline (Lee 2014 Tissue Antigens) + 다른 East Asian HT-specific HLA studies — `paperclip search "DPB1*05:01 Korean Hashimoto" -s abstracts` 로 빠르게 reference 확보
- ★★ **Paper 1 Bradley 2010 cite verify**: 본인 5/4 prep item — `paperclip search "BRAF V600E HLA-I downregulation thyroid 2010"` 로 정확 citation 확인
- ★ **Paper 2/3 cite expansion**: HT-specific vs GD-specific HLA literature 분리 검색
- ★ **specific paragraph 막힐 때**: 본인 voice 작성 중 cite 한 줄 verify 시 1-shot

**비고**: Paperclip 은 marathon-compliant tool — sprint generate 가 아니라 fact-check + cite verification. 본인 W1 voice-first 작업에 보조 도구로 OK.

**가입 필요?**: gxl.ai 계정 + API 또는 CLI install 필요. 본인이 cost-benefit 판단 (paid service 가능성).

---

# 7. 거기까지만 — 다음 step 명시

**오늘 (5/2 EOD 까지)**: 본 audit 까지 완료. 분석 실행 X.

**5/4 또는 5/5 (본인 prep day)**:
- 4.1 Paper 2 forest meta 재구성 (1-2 hr) — 본인이 Lee 2014 Korean baseline 확보 후 새 forest plot 2-panel
- 4.2 "autoimmune-PTC" → "Hashimoto-overlap PTC" 일괄 정정 (30 min)
- 4.3 memory 정리 (15 min)
- (선택) Paperclip 가입 + Bradley 2010 cite verify 사용

**5/5-5/10 W1 (Paper 1 voice-first)**:
- Paper 1 Section 2.4b + Discussion 3.3 의 1-2 line 표현 정정 (HT 명시)
- Hook 첫 줄 본인 voice
- Section 1.1 + 1.2 본인 voice

**Paper 2 본격 outline**: Paper 1 bioRxiv (6월) 후 (Phase 1, Q3-Q4 2026)

**Paper 3 (Graves)**: Bundang outreach 응답 wait → backburner

---

# 8. 한 줄 정리

**본인 catch 정확함. Paper 1 은 거의 영향 없음 (sub-B teaser + 1-line 표현 정정), Paper 2 forest meta framing 은 fundamental error (HT cohort vs GD cohort 직접 비교 → within-population elevation 두 separate panel 로 정정 필요), Paper 3 (Graves) 는 별도 trajectory backburner. 거기까지만 — 분석 실행은 5/4 본인 prep day 에. Paperclip 은 cite verify 보조 도구로 marathon-compliant.**
