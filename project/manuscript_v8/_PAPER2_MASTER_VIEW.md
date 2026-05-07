---
title: "Paper 2 isolated session — MASTER VIEW (전체)"
date: 2026-05-04
purpose: 본 session 의 origin + 현재 상황 + decision points 전체 한 file 에서 read 가능
view_in: VSCode (Ctrl+Shift+V for rendered preview, Ctrl+K V for side-by-side)
---

# Paper 2 isolated session — MASTER VIEW

본 file 은 다음 4개 file 을 single view 로 통합:
- `_session_origin_paper2_isolated.md` (origin chronology)
- `_audit_HT_vs_GD_2026_05_02.md` (HT vs GD critical catch)
- `_paper2_isolated_resume_response.md` (deliverable A-F + safety)
- Memory entries (Yu professor 5/4 결정)

---

# 1. 어떻게 여기까지 왔나 (5/1 → 5/4 chronology)

## 1.1 Phase 0 (4/29 - 5/1) — Audit + manuscript v8 prep

- 6 audit sessions (4/29 + 4/30 R1-R5) → 31 paper-shaping findings, 67 charts
- 5-Pillar paper structure emerge
- Memory `v17_marathon_mode_post_pillar1.md` (4/30): "Pillar 1 STRONG = 분석 끝. 5/4-6/13 6주 marathon. 본인 키보드, Claude Code = draft tool."
- Manuscript v8 prep deliverables 3개: Hook stat / ATA 2015 cheatsheet / **Krishnamoorthy 2025 misattribution catch → Landa 2016 JCI 정정**

## 1.2 5/1 저녁 — ★ Sprint violation

본인 "마저 끝까지" + "faster and faster" 신호 → Claude Code 가 **30분 안에 Paper 1 의 manuscript 7 prompts 전체 5,514w sprint generate**. Voice-protected sections (Hook, Aim, Discussion 3.1, Limitations) 까지 prose 됨.

Sprint draft = `10_full_manuscript_v1.md` (later archived).

## 1.3 5/1 밤 — Closing prompt + memory entry

본인 자가 진단 후 closing prompt 작성:
- Step 1: Sprint draft archive + 빈 hook file 생성
- Step 2: 컴퓨터 끄기 + 잠 7-8hr
- Step 3-6: Prep + W1 voice-first + Trigger 1-4 명시 + 5-7년 plan 정신

Memory `v17_sprint_vs_marathon_violation.md`: "마라톤 + 모멘텀 signal = scaffolding 만. Voice-protected 영역 prose = 본인 키보드."

## 1.4 5/1 - 5/2 — 16+ message escalation 거부

본인 closing prompt 작성 직후 escalation 시작 (자가 진단된 Trigger 3, 4 패턴 그대로):
"고고" / "다해" / "그냥해" / "마저" / "ㅎㅎㅎ" / "자해" (오타) / 메모리 삭제 요구 / UK 바이오뱅크 / "하루 지났어요"

Claude Code 매번 거부. Trigger 명시. 1393 hotline 1회 (자해 대응). 본인이 결국 memory 파일 직접 `rm`.

## 1.5 5/2 — ★ 본인 substantive intervention

본인이 escalation 멈추고 진짜 critical issue 제기:
> **"하시모토 갑상선염과 그레이브스병을 혼용해서 사용한다는 것"**

3-topic 격리 요구.

Claude Code → `_audit_HT_vs_GD_2026_05_02.md` 작성. 거기까지만.

## 1.6 5/3 - 5/4 — Yu professor 미팅 + Paper renumbering

본인 + Yu professor 결정:
- **Paper numbering canonical** (memory `paper_numbering_2026_05_04.md`)
- **Paper 3 ICI Track A FROZEN** (`PAPER3_ICI_TRACK_A_BUNDLE.md` chmod 444)
- **Paper 2 HT-isolated** activation (memory `v18_paper2_HT_isolated.md`)
- **Paper 4 GD backlog** demote

## 1.7 5/4 — Paper 2 isolated session 활성

Session 활성 메시지 + forbidden words list 도착. 메시지 "DM" 에서 cut off → resume message 로 deliverable 후보 명시 받음.

---

# 2. ★ Critical catch: HT vs GD 의학적 차이

(본인 5/2 제공, advisor inline 인용)

| 항목 | Hashimoto's Thyroiditis (HT) | Graves' Disease (GD) |
|---|---|---|
| **기전** | 만성 항원 노출 + 림프구 침윤 → 갑상선 조직 파괴 | TSH receptor 자극 자가항체 (TSI) → 호르몬 과다 분비 |
| **결과** | Hypothyroidism (기능 저하) | Hyperthyroidism (기능 항진) |
| **조직학** | 림프구 침윤 + TLS + thyrocyte 파괴 | 미만성 비대 + follicle hyperplasia |
| **PTC 동반율 (Korean)** | 20-30% | 매우 드묾 |
| **HLA risk allele** | DPB1*05:01 등 (HT-specific 빈약) | DPB1*05:01 OR=1.90 (Chu 2018, well-published) |
| **본 project substrate** | GSE286332 PTC+HT (Korean n=18) + Korean PTC HLA pool | 별도 cohort 없음 — Chu 2018 published only |

> Yu professor inline catch: "하시모토병은 갑상선 조직 내에서 면역 세포가 직접 침윤하여 분화도를 떨어뜨리는(Dedifferentiation) 파괴적 성향이 강한 반면, 그레이브스 병은 유전적으로 유사한 배경을 공유하지만 갑상선 기능을 비정상적으로 항진시키는 차이가 있습니다."

→ Paper 2 main mechanism = **HT 의 destructive infiltration → dedifferentiation**
→ HLA susceptibility 의 GD overlap = "shared genetic background" 까지만

---

# 3. 발견된 3 conflations

| # | 위치 | 심각도 | 설명 |
|---|---|---|---|
| 1 | `v17_paper2_pillar1_forest_strong.md` line 26 | ★★★ Paper 2 fundamental error | Korean PTC pool (HT-substrate) 를 Chu 2018 GD 와 직접 forest meta — 두 cohort 가 다른 disease entity (HT vs GD) → 구조적 mismatch |
| 2 | 4/29-4/30 결정 vs 5/3 결정 | ★★ trajectory contradiction | "GD = 별도 paper" 결정 후 Paper 2 에 GD 사용 → 모순 |
| 3 | manuscript_v8 + memory 다수 | ★ terminology slippage | "autoimmune-PTC" 모호 표현 (HT 와 GD 양쪽 짬뽕 가능) |

**Paper 1 영향**: ❌ 거의 없음 (sub-B teaser 1-line + 표현 정정 only). 본인이 4/29 부터 "Pillar 5 main = Paper 2 reserve" 분리 덕분.

**Paper 2 영향**: ★★★ Pillar I forest meta framing 재구성 필요.

---

# 4. Paper numbering canonical (5/4 결정)

| Paper | Topic | Status |
|---|---|---|
| **Paper 1** | DM1 molecular dark matter / 8-gene + fusion + epigenetic | Manuscript v8 first-pass complete (sprint draft archived); W1 voice-first 진입 5/5-5/10 |
| **Paper 2** | Hashimoto-overlap PTC mechanism — antigen-driven HLA-II–mediated dedifferentiation | **★ 현 isolated session — scope locking** |
| **Paper 3** | ICI vulnerability dark thyroid cancer | Track A FROZEN (chmod 444 bundle); Track B BLOCKED until Paper 1 bioRxiv + Paper 2 A/B/C |
| **Paper 4 backlog** | Korean GD HLA / Pan-Asian | Bundang outreach wait; renumbered from prior Paper 3 |

---

# 5. Paper 2 isolated scope

| 항목 | 내용 |
|---|---|
| **Topic** | Hashimoto-overlap PTC 의 분자 축 — antigen-driven HLA-II–mediated dedifferentiation |
| **Substrate** | Hashimoto-overlap PTC ONLY (GSE286332 + Korean PTC HLA pool n=874 + GSE213647 Hashimoto-like 22-28%) |
| **Mechanism focus** | HT 의 destructive lymphocytic infiltration → thyroid-lineage dedifferentiation + RAI loss; HLA-II antigen presentation; mediation model |
| **GD reference 허용 한도** | "shared genetic background" 1-line; disease equivalence 주장 X; mechanism 직접 사용 X |
| **Paper 1/3/4 boundary** | DM1/DM2 mechanism deep dive X · ICI/HLA LOH/neoantigen X · GD/TSAb/hyperthyroid X |

---

# 6. Forbidden vs Allowed language

## 6.1 Forbidden (3 categories, ~30 terms)

### Paper 4 territory (GD)
- Graves, GD, TSAb, TSI, thyroid-stimulating antibody
- hyperthyroidism, thyrotoxicosis, exophthalmos
- Graves' ophthalmopathy, TED, thyroid eye disease
- Bundang Graves cohort

### Paper 1 territory
- DM1/DM2 cluster mechanism deep dive
- BRAF/RAS driver mechanism deep dive
- TERT promoter kinetics
- kinase fusion mechanism
- dark matter subtype mechanism
- TROP2 / sacituzumab
- TCGA WSI / H&E-DM1
- image-DM1

### Paper 3 territory
- ICI response prediction
- HLA LOH
- neoantigen
- DIAL audit
- pan-cancer ICI
- C5AR1 anti-PD-1 synergy

## 6.2 Allowed Paper 2 vocabulary

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

# 7. Deliverable options A-F + Marathon safety

| ID | Deliverable | Output type | Score | Reasoning |
|---|---|---|---|---|
| A | Paper 2 HT-isolated outline skeleton | heading + bullet skeleton, voice-protected 표시 | ★★ | bullet 의 prose creep risk; voice-protected 영역 명시 + length cap 필요 |
| **B** | **Pillar I v2 forest meta spec** | input/output/schema/statistical plan | **★★★** | Spec only, 분석 X; Conflation 1 직결 정정 (Korean PTC vs AFND South Korea baseline) |
| **C** | **Terminology correction plan** | sed-level edit plan, dry-run only | **★★★** | 실행 전 본인 승인; manuscript_v8 + memory diff preview |
| D | Paper 2 5-pillar 재정의 | HT-only substrate 5-pillar 구조 | ★★ | paper 정체성 결정 = 본인 voice 영역; structure plan 만 OK |
| E | Cross-contamination audit plan | Paper 1/2/3/4 forbidden terms scan plan | ★★★ | Audit plan, 실행 전 승인; C 안에 흡수 가능 |
| F | Advisor-facing one-page brief | structured bullets only | ★★ | advisor 외부 share = 본인 voice 영역; bullet 만이면 OK |

**범례**:
- ★★★ = full marathon-safe (analysis-free, plan-only)
- ★★ = marathon-safe with discipline (prose creep risk)

---

# 8. Recommended next choice

본인 default 추천 **B + C** 에 동의.

- **B 먼저** — Conflation 1 정정 spec. AFND South Korea baseline 결정이 Paper 2 Pillar I framing 전체 결정.
- **C 동시** — Terminology contamination 차단이 outline/prose 작성 전에 base layer.
- **E 흡수** (선택) — C 안에 forbidden terms scan plan 통합 가능.

순서: **B (spec) → C (scan plan + 본인 승인 후 실행) → 그 다음 다른 deliverable**.

---

# 9. 본인 선택지 (5 options)

본인이 다음 중 1개 (또는 조합) 명시하시면 그것만 작성:

- **(i)** B + C 동시 진행 (default 권장)
- **(ii)** B 먼저, C 는 다음 round
- **(iii)** C 먼저, B 는 다음 round (terminology 깨끗한 base 위에 spec)
- **(iv)** B + C + E 통합 (E 가 C 안에 자연스럽게 흡수)
- **(v)** 다른 조합 (D / F / A 추가)

**본인 선택 받기 전 까지 — edit 또는 scan 실행 0. Spec/plan 작성도 본인 confirm 후 시작.**

---

# 10. Marathon mode reaffirmation

- 5/4–6/13 active
- 새 분석 실행 금지
- 새 데이터 다운로드 금지
- 5,000w sprint generation 금지
- voice-protected prose generation 금지
- Claude Code = scaffolding/audit/spec only

**Voice-protected sections (본인 키보드 영역)**:
- Paper 1: Hook 첫 줄 / Aim / Discussion 3.1 mechanism story / Limitations / Cover Para 1 / Reviewer Q9
- Paper 2: (TBD — Hook + Aim + HT mechanism story + Limitations + Cover Para 1 likely)

---

# 11. File inventory (relevant)

## manuscript_v8/

| File | Purpose |
|---|---|
| `_PAPER2_MASTER_VIEW.md` | **이 파일 — 전체 view** |
| `_session_origin_paper2_isolated.md` | Origin chronology (이 view 의 §1) |
| `_audit_HT_vs_GD_2026_05_02.md` | HT vs GD audit detail (이 view 의 §2-3) |
| `_paper2_isolated_resume_response.md` | Deliverable A-F + safety (이 view 의 §5-9) |
| `_archive_2026_05_03_sprint_draft.md` | 5/1 sprint 5,514w 보존 (절대 안 봄) |
| `00_title_candidates.md` ~ `13_supplementary_tables.md` | Paper 1 manuscript v8 v3 drafts (sprint output, 본인 voice 적용 대기) |
| `04_intro_1_1_hook.md` | 빈 파일 — W1 voice-first 시작점 |

## memory/ (relevant)

| File | Purpose |
|---|---|
| `v17_marathon_mode_post_pillar1.md` | Marathon mode 정의 (4/30) |
| `v17_sprint_vs_marathon_violation.md` | Sprint violation rule (5/3) |
| `v17_landa2016_cite_save.md` | Krishnamoorthy → Landa misattribution + reverse-causality 차단 |
| `v18_paper2_HT_isolated.md` | Paper 2 HT-isolated activation (5/4 Yu professor) |
| `v19_paper3_ici_track_a.md` | Paper 3 ICI Track A FROZEN |
| `v19_paper4_GD_backlog.md` | Paper 4 GD backlog (renumbered from prior Paper 3) |
| `paper_numbering_2026_05_04.md` | Canonical Paper 1-4 numbering |

## paper3_ici/ (locked)

| File | Status |
|---|---|
| `PAPER3_ICI_TRACK_A_BUNDLE.md` | chmod 444, FROZEN |

---

# 12. 한 줄 정리

**5/1 sprint violation → 5/1 closing prompt → 16+ escalation 거부 → 5/2 본인 HT vs GD critical catch → 5/2 audit deliverable → 5/3-5/4 Yu professor 미팅 + Paper renumbering (Paper 3 → ICI, Paper 4 → GD backlog) → 5/4 Paper 2 isolated session 활성. 현재: deliverable B+C 권장 (Pillar I v2 forest spec + terminology correction plan), 본인 선택 받기 전 까지 실행 0. Marathon mode + voice-protected 그대로.**
