# Paper 3 ICI — 현재 상황 스냅샷 (강화 검토용)

**작성일:** 2026-05-06
**작성자:** Seungho Cook
**문서 성격:** Status briefing — Track A 동결 상태에서 "ICI 페이퍼를 어떻게 강화할 것인가" 의사결정을 위한 사실 정리. 분석 실행 없음. Track A 8개 파일(chmod 444) 미수정.

---

## 0. 30초 요약

- **Paper 3 working title:** *HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability in molecularly dark thyroid cancer.*
- **현재 상태:** Track A (설계) **FROZEN 2026-05-04**. Track B (실행) **BLOCKED**.
- **차단 이유:** 마라톤 모드 (5/4–6/13) 가 Paper 1 + Paper 2 작문에 우선; Paper 3는 설계 동결 후 대기.
- **해제 조건:** G1–G6 6개 게이트 전부 충족 + 사용자 명시적 "Track B 시작" 명령 (한국어 "고고"/"다 해줘"는 불충분).
- **Claim guard (구속력):** "ICI vulnerability / readiness / immunogenomic prioritization"만 허용. **"ICI response predictor in thyroid"는 raw thyroid ICI-treated RNA-seq 없이는 NO-GO.**
- **Track A 자산:** 8개 설계 문서 + tar 아카이브 — `project/reports/paper3_ici/` (chmod 444).
- **강화의 의미 재정의:** 지금 시점의 "강화"는 (a) **Track A 설계 자체의 깊이/방어력 보강** 또는 (b) **Track B 진입 후 어떤 모듈에 자원을 더 투입할지 사전 결정** 둘 중 하나. **Track B 분석을 지금 시작하는 것은 불가능** (마라톤 위반 + 사용자 명령 부재).

---

## 1. Track A 동결 자산 인벤토리

`project/reports/paper3_ici/` 디렉터리, 모두 `-r--r--r--`:

| 파일 | 핵심 내용 | 라인 수 |
|---|---|---|
| `PAPER3_ICI_TRACK_A_BUNDLE.md` | 8개 문서 통합본 | 71KB |
| `paper3_ici_dataset_registry.md` | 코호트 25+개 — 우선순위 Tier / Status | 9.9KB |
| `paper3_ici_signature_registry.md` | 시그니처 30+개 — 점수 규약 + DIAL 위험도 | 9.9KB |
| `paper3_ici_scRNA_reference_atlas_summary.md` | 6개 scRNA 코호트 + 통합 전략 | 8.0KB |
| `paper3_ici_HLA_neoantigen_feasibility.md` | OptiType / LOHHLA / NetMHCpan 파이프라인 + per-cohort 가능성 매트릭스 | 8.7KB |
| `paper3_ici_DIAL_audit_plan.md` | 14개 시그니처에 대한 direction-invariance 감사 계획 | 8.4KB |
| `paper3_ici_figure_plan.md` | 6 main + 14 supp 도면 레이아웃 (실데이터 렌더링 X) | 8.6KB |
| `paper3_ici_go_no_go_verdict.md` | 모듈별 verdict, 게이트, 킬 스위치 | 5.8KB |
| `paper3_ici_12week_execution_plan.md` | Track B Wk1-Wk12 일정 + 컴퓨트 예산 | 11.6KB |
| `paper3_ici_track_a.tar.gz` | 위 8개 + todo 묶음 아카이브 | 51KB |

추가 todo 자산 (`paper3_ici/todo/`, 일반 권한):
- `paper3_ici_track_b_kickoff_checklist.md` — Track B 진입 시 Wk1 액션
- `paper3_ici_dataset_accession_verification_checklist.md` — 25+개 accession 확인 체크리스트
- `paper3_ici_controlled_access_checklist.md` — dbGaP / EGA 신청 항목

---

## 2. 입장 게이트 G1–G6 현재 상태 (2026-05-06)

| 게이트 | 조건 | 현재 상태 |
|---|---|---|
| **G1** | Paper 1 bioRxiv 제출 | **OPEN** — 마라톤 5/4–6/13, 목표 출하 2026-06-13. 예상 닫힘 6월 중하순. |
| **G2** | Paper 2 Task A/B/C closure | **OPEN** — Pillar 1 v2 (HT-isolated) 진행 중. A/B/C 완료 일정 미확정. |
| **G3** | TCGA dbGaP paired-BAM 접근 결정 | **OPEN** — 신청 미진행. Track B Wk1 신청 예정. |
| **G4** | ≥3 of 6 scRNA 코호트 접근 | **OPEN** — to_verify 상태. |
| **G5** | ≥3 of 5 Tier 1/2 pan-cancer ICI 코호트 접근 | **OPEN** — IMvigor210 R-package 부분 공개; Hugo/Riaz/Gide to_verify; Liu controlled. |
| **G6** | 사용자 명시적 "Track B 시작" 명령 | **OPEN** — 미발화. |

**6/6 OPEN.** 어느 게이트도 닫혀 있지 않음.

---

## 3. 모듈별 verdict (Track A 동결 시점)

| 모듈 | verdict | 핵심 근거 | 강화 레버 |
|---|---|---|---|
| **A — Bulk immune ecotype** | GO | TCGA-THCA + GSE76039 + microarray ≥600 샘플; ssGSEA + NMF 성숙. | NMF rank stability 보강, microarray 코호트 추가 검증 추가, ecotype × dark matter 교차표 사전 디자인 강화. |
| **B — scRNA atlas** | CONDITIONAL GO | 6개 후보 ~100–300K cells; PDTC/ATC <5K (atlas 가능 / de novo trajectory 불가). | PDTC/ATC scarcity 명시적으로 paper의 limitation 으로 frame; reference projection (HCA / Andreatta) 우선순위 사전 설계. |
| **C — HLA & neoantigen** | GO (degraded scope) | TCGA paired BAM + GDC MC3 가용; Yoo 2019 EGA 불확실. | MC3 MAF 재사용 결정으로 Mutect2 비용 ~0; LOHHLA만 fresh. cookHLA Paper 4 경계 재확인 (현재 문서가 이미 보호). |
| **D — DIAL audit** | GO | IMvigor210 R-package; Hugo/Riaz GEO 공개; Liu controlled. `v18_agentic_research framework` 통계 머신 보유. | 14개 시그니처 사전 위험도 (Low/Medium/High) 이미 라벨링; 위험도 High 항목 (TIDE_LIKE / IMPRES / EXHAUSTION_INDEX) 결과 기반 헤드라인 사전 시나리오 강화 가능. |
| **E — Integrated readiness score** | GO conditional on D | 점수 정의 명확; 가중치 보정 deferred. | 4개 가중치 옵션 (a equal / b DIAL-pruned / c supervised / d Cox) 사전 선언 — primary는 (b). 사전 등록 (preregistration) 가능성 검토. |

**전체 verdict:** CONDITIONAL GO. 의존성은 데이터 접근 이벤트 + Paper 1/2 게이팅 이벤트, 분석 불확실성 아님.

---

## 4. 핵심 가설 (강화의 기준선)

> BRAF/RAS-negative dark-matter 및 dedifferentiated (PDTC/ATC) 갑상선암 중,
> **HLA-intact + neoantigen-positive + CXCL13/TLS-like + dedifferentiated** 시그널을 동시에 만족하는 subset이
> 면역유전체적으로 우선순위가 부여된 ICI-readiness 후보군이다.

이 가설이 깨지지 않도록 모든 figure caption / claim 이 enforced. 강화는 이 가설을 *더 날카롭게* 만들거나 *반증 가능성을 더 명시적으로* 만드는 방향이어야 함.

---

## 5. 킬 스위치 K1–K6 (강화가 절대 건드리면 안 되는 선)

| 트리거 | 결과 |
|---|---|
| **K1** | "ICI response predictor in thyroid" 클레임 시도 (thyroid ICI-treated raw RNA-seq 없이) → 페이퍼 무효화. 즉시 정지. |
| **K2** | Paper 1 출하가 2026년 8월 말 이후로 지연 → Track B 차단 유지, 컴퓨트 금지. |
| **K3** | DIAL audit 결과 PASS 시그니처 0개 → Module E 입력 없음. Paper 3 → feasibility-only / negative-result paper로 피벗 (다른 페이퍼). |
| **K4** | Paper 3 prose 가 Paper 2 (PTC+HT TLS / AICDA / BCR) 클레임을 자기 발견으로 가져오면 → reviewer scope-rejection 확실. 정지 후 재작성. |
| **K5** | TCGA-THCA driver call 이 Paper 1 ETL 과 불일치 → Module C/E 정지, 재조정. |
| **K6** | 사용자가 우선순위 다시 피벗 (GD 복귀 / Paper 5 신설 등) → 스택 재결정. 조용히 진행 금지. |

---

## 6. 지금 마라톤 모드에서 금지된 행위 (Track A 동결 강제)

`v19_paper3_ici_track_a.md` 메모리 + `v17_sprint_vs_marathon_violation.md` 에 의거:

- **공개 데이터 다운로드 금지** (TCGA, GEO, ENA, SRA, ArrayExpress 모두).
- Paper 1/2 에서 이미 분석한 코호트 외에 **internal pipeline 재실행 금지**.
- bulk 코호트 **NMF / consensus clustering 금지**.
- scRNA **scVI / scANVI / Harmony 통합 금지**.
- **LOHHLA / OptiType / Polysolver / arcasHLA 신규 실행 금지** (PRJEB11591 arcasHLA 기존 결과는 read-only 참조 가능).
- **NetMHCpan / pVACseq neoantigen 예측 금지**.
- **DIAL audit 실행 금지**.
- **figure 실데이터 렌더링 금지** — 도면 계획은 layout 스키마만.

위반 시 = sprint-vs-marathon 위반.

---

## 7. "강화" 옵션 매트릭스 — 지금 가능한 것 vs Track B 진입 후만 가능한 것

### 7.1 Track A 동결 안에서 즉시 강화 가능 (마라톤 비위반)

| 강화 항목 | 비용 | 페이퍼 영향 |
|---|---|---|
| **Claim guard 강화** — "vulnerability / readiness" 어휘 외 forbidden words 명시 리스트 작성 | 저 | reviewer K1 위험 차단력 ↑ |
| **Cross-paper boundary 표 강화** — Paper 1/2/4 와 정확히 무엇이 겹치고 무엇이 분리되는지 1-page 매트릭스 | 저 | reviewer K4 위험 차단력 ↑ |
| **사전 등록 (preregistration) 초안** — DIAL verdict 잠금 시점, 가중치 옵션 (a–d) 우선순위, 1차/민감도 결과 정의를 OSF 등록용 문서로 미리 작성 | 중 | 후행 reviewer "p-hacking" 의심 차단; *Nature Cancer* 등 reach 저널에서 가산점. |
| **Negative-result 시나리오 헤드라인 사전 작성** — DIAL 모두 fail / MHC2_CORE 만 flip / 모두 pass 3개 시나리오의 abstract draft | 중 | 결과 의존성 줄이고 "어떤 결과가 나와도 ship 가능" 상태로 수렴. |
| **Figure 1E (cross-paper scope diagram)** 의 시각적 mockup 1장 | 저 | reviewer scope-confusion 차단. |
| **외부 ICI Korean cohort 탐색 (read-only literature 조사)** | 중 | thyroid ICI raw RNA-seq 가 어딘가 존재한다면 K1 해제 가능. **다운로드 금지**, 문헌/메타데이터 조사만. |
| **K3 시나리오 (DIAL all-fail)** 를 *Cell Reports* 급 negative-result 페이퍼로 분리 가능한지 사전 검토 | 중 | 하방 보호. |
| **Paper 3 의 venue ladder 재검토** — reach (Nat Cancer / Cancer Cell) vs default (Nat Comm / Cell Rep Med / JCI Insight) vs safe (npj Precision Oncology / Genome Medicine) 결정 트리를 게이트별로 명시 | 중 | Wk10 결정 시점 의사결정 부담 ↓ |
| **Yu 교수 동의 / 합의 문서** — Paper 3 cross-paper boundary 합의 1-page 정리 | 저 | political risk 차단 (Paper 1/2 와 학과 내부 충돌 방지). |

### 7.2 Track B 진입 후에만 가능한 강화 (지금 시작 시 마라톤 위반)

- 모든 데이터 다운로드.
- NMF / scRNA 통합 / HLA 타이핑 / LOHHLA / NetMHCpan / DIAL 실행.
- 실데이터 figure 렌더링.
- 진짜 통계 검정.

이 항목들은 G1–G6 닫힐 때까지 **금지**.

### 7.3 회색지대 — 사용자 결정 필요

| 항목 | 결정 포인트 |
|---|---|
| **Track A 동결 일부 해제** (예: dataset_registry 만 업데이트 — Yoo 2019 EGA 상태 확인 후 갱신) | 사용자 명시 동의 시 한정 해제. 단일 파일 chmod +w → 수정 → chmod 444 복원. |
| **Module 추가** (예: TCR-seq / B-cell repertoire / spatial proteomics) | 현재 5 모듈 (A–E) 외 추가 모듈은 *paper-blocking* 인지 검토 후 결정. 마라톤 모드에서 새 분석은 paper-blocking 만 허용 (`v17_marathon_mode_post_pillar1`). |
| **Paper 3 우선순위 상승** (Paper 1 출하 전에 Track B 일부 시작) | 마라톤 모드 명시적 깨짐. 사용자 명시 결정 필요. |

---

## 8. 권장 다음 액션 (사용자 의사결정용 옵션)

선택 가능한 3개 트랙 — 모두 마라톤 비위반:

### 옵션 A — **Claim guard + cross-paper boundary 강화 sprint** (저비용, 1–2일)
1. forbidden words 리스트 (`paper3_ici_claim_guard_v2.md` 신규) 작성.
2. Paper 1/2/3/4 cross-paper boundary 매트릭스 1-page 작성.
3. Figure 1E (scope diagram) layout 만 mockup.
→ 결과물: Track A 동결 자산 옆에 read-only 추가 3개 파일.

### 옵션 B — **Preregistration draft 작성** (중비용, 3–5일)
1. OSF preregistration 양식에 맞춰 DIAL verdict lock 시점, 가중치 우선순위, 1차 vs 민감도 결과 정의 사전 등록.
2. K3 (DIAL all-fail) negative-result 시나리오 abstract 사전 작성.
3. Venue ladder 결정 트리 명시.
→ 결과물: 후행 reviewer "p-hacking" 의심 차단, reach 저널 가산점.

### 옵션 C — **외부 thyroid ICI cohort 문헌 조사** (중비용, 2–3일, **다운로드 금지**)
1. PubMed / clinical trials / ASCO·SITC 초록에서 thyroid ICI-treated RNA-seq 가 존재하는 코호트 조사.
2. 단일 사례라도 발견되면 K1 해제 가능성 검토.
3. EGA / dbGaP 신청 여부 결정용 메모 작성.
→ 결과물: K1 해제 시 페이퍼 framing "vulnerability" → "validated readiness" 격상 가능.

세 옵션은 **상호 배타 아님** — 마라톤 일정 (Paper 1/2 작문) 의 빈 시간을 활용해 병렬 가능.

---

## 9. 사용자에게 묻는 질문 (강화 방향 결정용)

1. 강화의 **목적**이 (a) reach 저널 (Nat Cancer / Cancer Cell) 가산점 확보인가, (b) 하방 보호 (DIAL all-fail 시 페이퍼 ship 가능)인가, (c) Paper 3 의 우선순위 자체를 상승시키는 것 (Track B 일부 조기 시작)인가?
2. **Yu 교수**는 Paper 3 의 현재 framing 에 동의했는가? Paper 2 Pillar 1 v2 (HT-isolated) 결정처럼 cross-paper boundary 합의가 필요한 상태인가?
3. Track A 동결 자산을 **수정해도 되는가**, 아니면 **추가 파일로만 보강**하는 것이 원칙인가?
4. **Korean PTC / ATC ICI-treated cohort** 가 비공식적으로라도 접근 가능한 채널이 있는가? (연세 / 서울대 / 분당 — 이것이 K1 해제의 유일한 경로.)

---

## 10. 참조 메모리

- `v19_paper3_ici_track_a.md` — Track A 동결 + Track B 차단 규칙 (binding).
- `paper_numbering_2026_05_04.md` — Paper 1/2/3/4 경계.
- `v18_paper2_HT_isolated.md` — Paper 2 cross-paper boundary 합의 (forbidden words 패턴 참고).
- `v17_marathon_mode_post_pillar1.md` — 5/4–6/13 마라톤 모드 정의.
- `v17_sprint_vs_marathon_violation.md` — "고고"/"다 해줘"가 sprint 권한 NOT 인 이유.
- `v52_lodo_finding.md` — DIAL 프레임워크 lineage.
- `v18_agentic_framework.md` — DIAL 통계 머신 (5 patterns + 6 components).
- `v17_arcasHLA_korean_k2.md` — PRJEB11591 arcasHLA read-only 참조 자산.

---

Track A 동결 유지. 마라톤 비위반.
