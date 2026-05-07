---
title: "Paper 2 HLA β — Yu professor send/meeting package"
date: 2026-05-04
audience: 본인 → Yu Hyeong-won 전달용
parent_packet: project/reports/2026_05_04_hla_beta_advisor_packet.md
parent_plan: project/reports/2026_05_04_hla_strengthening_beta_plan.md
status: COMMUNICATION DRAFT ONLY — no execution, no analysis, no lookup, no forest run
forbidden: Lee 2014 lookup execution · 6-allele forest execution · B5/B3/B4/C4 execution · Paper 3/4 touch · GD/Graves analysis · manuscript prose · voice-protected sections · new data download
---

# Paper 2 HLA β — Yu professor send/meeting package

---

## [1] KakaoTalk message (5 lines, Korean, polite)

```
교수님 안녕하세요.
Paper 2 HLA 작업 정리해서 advisor packet 만들어 두었습니다.
v1 (Korean PTC vs Chu 2018 GD 직접 비교) 은 HT/GD 혼용 문제로 폐기하고,
v2 = Korean PTC pool n=874 vs Korean baseline (HT-only scope) 으로 재구성했습니다.
Q1–Q7 결정만 미팅에서 함께 보면 좋을 것 같아 편하신 시간 알려주시면 감사하겠습니다.
```

---

## [2] Email version (Korean)

### Subject

`[Paper 2 HLA] β plan + advisor packet 검토 요청 — Q1–Q7 decisions`

### Body

```
교수님,

Paper 2 HLA 작업 정리를 마무리했습니다. 짧게 요약드립니다:

1. Scope: Hashimoto-overlap PTC HLA susceptibility + HLA-II–mediated dedifferentiation mechanism (HT-only).
2. v1 (Korean PTC pool n=874 vs Chu 2018 Han Chinese GD 직접 forest meta) 은
   HT-substrate 와 GD-substrate 의 cross-disease mismatch 로 폐기했습니다.
3. v2 = Korean PTC pool n=874 vs Korean general population baseline
   (AFND South Korea pool 또는 + Lee 2014 KOTRY donor reference) 으로 재구성.
   v2 4-allele forest 는 이미 실행되어 있습니다 (project/results/p2_pillar1_forest_v2/).
   Lee 2014 fill 시 6-allele (C*01:02, DQB1*02:01 추가) 으로 확장 가능합니다.
4. GD / Graves / Bundang 관련 작업은 Paper 4 backlog 로 분리했습니다
   (4/4 gating: Paper 1 bioRxiv + Paper 2 A/B/C + Bundang n>50 + 교수님 합의).

미팅에서 결정 부탁드릴 사항 (Q1–Q7):
  Q1. Pillar I v2 (Korean PTC vs Korean baseline) framing 승인 여부
  Q2. Lee 2014 fill (C*01:02 / DQB1*02:01) 사용 여부 → 6-allele 확장
  Q3. 6-allele forest 를 primary 로 할지 supplementary 로 할지
  Q4. Baseline pool 구성 — AFND only / AFND + KOTRY(Lee 2014) / 다른 옵션
  Q5. Harbin Korean fallback (있을 경우) 처리 — 제외 / sensitivity / primary 포함
  Q6. Tier 2 spec 우선순위 — B5(antigen presentation) / B3(sc HLA-II) / B4(BCR/TLS×HLA-II) / C4(TCGA HLA QC)
  Q7. 실행 시점 — 지금 / Paper 1 bioRxiv 이후

각 Q 에 대한 추천 default 는 advisor packet §5 에 정리되어 있습니다.

첨부 / 링크:
  - Advisor packet: project/reports/2026_05_04_hla_beta_advisor_packet.md
  - β plan (full spec): project/reports/2026_05_04_hla_strengthening_beta_plan.md
  - v2 4-allele forest 기존 결과: project/results/p2_pillar1_forest_v2/PILLAR1_FOREST_V2_SUMMARY.md
  - HLA strengthening situation overview: project/reports/2026_05_04_HLA_strengthening_situation.md

새 분석 / lookup / forest 실행은 교수님 결정 후 시작하겠습니다.

미팅 가능하신 시간 알려주시면 감사하겠습니다.

감사합니다,
Seungho Cook
```

### Decision requested

Yu professor 의 Q1–Q7 답변. 미팅 시 §[3] decision sheet 직접 fill-in 또는 email reply 모두 가능.

---

## [3] Q1–Q7 decision sheet

| Q | Decision needed | Recommended default | Yu decision | Note |
|---|---|---|---|---|
| **Q1** | Pillar I v2 framing (Korean PTC vs Korean baseline) approve? | ✅ **YES** approve v2 | ⏸ | v1 의 HT/GD cross-disease conflation 해소; reviewer-defensible |
| **Q2** | Lee 2014 fill for C\*01:02 + DQB1\*02:01 approve? | ✅ **YES** Lee 2014 fill | ⏸ | 4-allele → 6-allele 확장. G3 = YES (user 5/4 lock); advisor confirm |
| **Q3** | 6-allele forest primary OR supplementary? | **6-allele primary IF Lee 2014 추출 깨끗; otherwise 4-allele primary + 6-allele supp** | ⏸ | Quality-conditional default |
| **Q4** | AFND South Korea baseline pool 구성? | **AFND + KOTRY (Lee 2014)** — G2-2 | ⏸ | 두-source 풀, 적당한 n, KOTRY = 인정된 한국 reference; heterogeneity 공개 manageable |
| **Q5** | Harbin Korean fallback (AFND 에 있을 경우) 처리? | **Sensitivity panel only, not primary** | ⏸ | Lee 2014 가 primary 알릴 cover 시 cleaner Korean baseline |
| **Q6** | Tier 2 priority — B5 / B3 / B4 / C4? | **B5 first** (antigen presentation gene score) | ⏸ | 기존 v17_hla data 활용; lowest new-execution cost; Paper 2 mechanism core |
| **Q7** | Execute now OR after Paper 1 bioRxiv? | **Paper 1 user-voice cleanup 후 ~ Paper 1 bioRxiv 제출 사이** | ⏸ | Paper 2 Pillar I lock 을 Paper 1 ship 전에 확보 |

---

## [4] If-then execution mapping

| Condition | Next step | Command (per β plan §7) |
|---|---|---|
| **Q1 YES + Q2 YES + Q4 chosen + Q5 chosen** | Tier 1 unblocked → A2 lookup 시작 | `run A2 lookup only` |
| **A2 Lee 2014 fill output 깨끗 (post-A2 verify)** | A1 6-allele forest script 작성 (실행 X) | `prepare A1 forest script only` |
| **A1 script 검토 OK + 다시 사용자 confirm** | A1 6-allele forest pipeline 실행 | `run full A1 forest` (사용자 explicit) |
| **Q6 = B5 first + Tier 1 완료** | Tier 2 B5 spec → execution | `execute B5` |
| **Advisor says "defer Pillar I"** | Pillar I freeze, Paper 1 W1-W6 marathon 집중 | `stop and bring to Yu meeting` (이미 한 상태) → freeze |
| **Advisor says "move HLA to Paper 4"** | Paper 2 scope memo 업데이트만, HLA 분석 X | `v18_paper2_HT_isolated.md` 메모리 업데이트, 분석 0 |
| **Q7 = after Paper 1 bioRxiv** | 모든 HLA execution 을 Phase 1 (Q3-Q4 2026) 으로 deferred | Freeze, W1-W6 Paper 1 집중 |
| **Q3 = 4-allele primary** | 기존 v2 4-allele 유지, 6-allele = supplementary 만 | A2 + A1 still useful for supp; revised script |
| **Q1 NO (v2 reject)** | β plan 전체 freeze, Paper 2 framing 재논의 | 모든 execution 0; Paper 2 outline 재작성 결정 필요 |

---

## [5] 60-second meeting script (bullets only)

**Slide-less spoken version, ~60초:**

- **Paper 2 HLA scope (5초)**: Hashimoto-overlap PTC HLA susceptibility + HLA-II–mediated dedifferentiation. **HT-only**. GD 는 Paper 4 backlog 분리.
- **v1 problem (10초)**: Korean PTC pool (HT-substrate) 을 Chu 2018 Chinese GD cohort 와 직접 forest meta. **Cross-disease mismatch** — reviewer 가 1차 fact-check 에서 잡힘.
- **v2 fix (10초)**: Korean PTC pool n=874 (K2 235 + Lee 630 + GSE286332-PTC 9) **vs Korean general population baseline** (AFND + Lee 2014 KOTRY). **Same population, within-elevation**. 4-allele forest 이미 실행 완료.
- **3 key decisions (15초)**:
  1. **Q1** — v2 framing 승인 (recommended YES)
  2. **Q2** — Lee 2014 fill 시 6-allele 확장 (recommended YES, G3 user-lock 됨)
  3. **Q4** — baseline pool 구성 (recommended AFND + KOTRY/Lee 2014)
- **추가 (10초)**: Q3 primary/supp · Q5 Harbin handling · Q6 Tier 2 priority (B5 추천) · Q7 timing (Paper 1 bioRxiv 전)
- **Marathon mode (10초)**: 새 분석 / Lee 2014 lookup / forest 실행 — **교수님 결정 후 시작**. 지금까지 완료된 건 spec/plan + 기존 4-allele forest + terminology cleanup.

---

## [6] Files to attach / open

### 핵심 첨부 (advisor packet 우선)

| Priority | File | Purpose |
|---|---|---|
| ★★★ | `project/reports/2026_05_04_hla_beta_advisor_packet.md` | **1-page packet, 미팅 메인 문서** |
| ★★★ | `project/reports/2026_05_04_hla_strengthening_beta_plan.md` | β plan full spec (Tier 1 + Tier 2) |
| ★★ | `project/results/p2_pillar1_forest_v2/PILLAR1_FOREST_V2_SUMMARY.md` | 기존 v2 4-allele forest 결과 |
| ★★ | `project/reports/2026_05_04_HLA_strengthening_situation.md` | Cross-paper situational overview (24 options + scoring) |

### 보조 cross-reference (필요 시)

| File | Purpose |
|---|---|
| `project/manuscript_v8/_audit_HT_vs_GD_2026_05_02.md` | HT vs GD 의학적 차이 + 3 conflations 발견 (v1 deprecation 근거) |
| `project/manuscript_v8/PAPER2_PILLAR1_V2_FOREST_META_SPEC_2026_05_04.md` | Pillar I v2 8 decision gates G1-G6 spec |
| `project/manuscript_v8/_LEE_2014_FILL_DECISION_LOCK_2026_05_04.md` | G3 = YES user lock (5/4) |
| `project/manuscript_v8/_PAPER2_MASTER_VIEW.md` | Paper 2 isolated session origin + scope full chronology |

### Memory cross-reference (text-based context)

- `v18_paper2_HT_isolated.md` — Paper 2 HT-isolation 5/4 decision (canonical)
- `v17_paper2_pillar1_forest_strong.md` — v1 history (SUPERSEDED note)
- `v19_paper4_GD_backlog.md` — Paper 4 GD gating
- `paper_numbering_2026_05_04.md` — Paper 1-4 canonical numbering

### 미팅 자리에서 펴 둘 것

1. Advisor packet (§[3] decision sheet 인쇄해서 직접 fill-in 권장)
2. β plan (Q 답변에 따라 §6 execution gates 참조)
3. v2 forest existing result (Q1 framing 검토 시 reviewer 시각 sanity check)

---

**Yu send package ready. No execution performed.**
