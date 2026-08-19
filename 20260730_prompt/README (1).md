# DM1 Nature Communications — Claude Code Submission Kit

이 폴더는 Paper 1 (DM1) 저장소를 **투고 직전 수준으로 감사·수정·검증·패키징**하기 위한 Claude Code 프로젝트 설정입니다.

단순한 영문 교정 프롬프트가 아닙니다. 다음을 분리된 reviewer agent로 수행합니다.

- 주장–근거–수치 추적
- 통계 모델과 reference direction 감사
- 임상적 과장 및 RAI 관련 인과 비약 차단
- Figure/legend/source-data 완결성 검사
- Methods 재현성 및 데이터·코드 가용성 검사
- 인용 정확성 및 허위 citation 차단
- 적대적 peer review와 편집자 desk-reject simulation
- Nature Communications 제출 폴더 빌드
- 저자 고유 voice가 필요한 6개 구역 보호

## 설치

이 키트의 **내용물**을 DM1 repository root에 복사합니다. 기존 `CLAUDE.md`가 있으면 먼저 백업하고 두 파일의 규칙을 수동 병합하십시오.

예시 구조:

```text
<DM1_REPO>/
├── CLAUDE.md
├── .claude/
│   ├── agents/
│   ├── rules/
│   └── skills/
├── manuscript_v8/
│   ├── NATURE_COMMUNICATIONS_FULL_DRAFT_v2_2026_07_08.md
│   └── SUBMISSION_CONTROL/
└── ...
```

Claude Code를 repository root에서 시작합니다.

```bash
cd <DM1_REPO>
claude
```

첫 세션에서 아래 순서로 실행합니다.

```text
/context
/manuscript-audit
/figure-lock
/voice-interview
/nc-submission repair
/submission-build
```

`/context`에서 root `CLAUDE.md`와 `.claude/rules/*.md`가 로드됐는지 확인하십시오.

## 권장 운영 순서

### 1. Audit only

```text
/manuscript-audit
```

파일을 고치기 전에 현재 원고·수치·그림·코드·인용의 불일치를 보고서로 냅니다. Audit agent는 원칙적으로 원고를 수정하지 않습니다.

### 2. Figure and number lock

```text
/figure-lock
```

Figure 번호, 패널, legend, source data, manuscript citation을 하나의 manifest로 묶습니다. 숫자나 reference group이 불명확하면 수정하지 않고 blocker로 남깁니다.

### 3. Author voice capture

```text
/voice-interview
```

Claude가 저자 대신 임상 동기나 해석을 발명하지 않도록, 6개 voice-protected section에 대해 국승호 저자의 원문 답변을 먼저 받습니다.

### 4. Repair

```text
/nc-submission repair
```

독립 reviewer 결과를 합친 뒤 수정합니다. 수정 전후 claim ledger와 git diff를 남깁니다.

### 5. Final build

```text
/submission-build
```

업로드 가능한 파일 세트와 누락 행정 목록을 만듭니다. **submission portal 업로드나 외부 공개는 자동 수행하지 않습니다.**

## 중요 원칙

1. Audit-locked 숫자는 원본 output을 확인하지 않고 변경하지 않습니다.
2. 숫자는 맞아도 reference group 또는 endpoint 방향이 불명확하면 문장화하지 않습니다.
3. RAI 흡수 또는 치료 반응을 직접 측정하지 않은 분석을 “predicts RAI response”라고 쓰지 않습니다.
4. `DM1 × BRAF` PFI interaction은 치료효과 interaction이 아니므로 기본적으로 “predictive biomarker”라고 부르지 않습니다.
5. methylation은 인과가 아니라 epigenetic correlate입니다.
6. IHC 3-plex는 in-silico proxy이며 임상 검증 assay가 아닙니다.
7. 내부 실험은 2026-07-30 중단 확정입니다. 원고에서 “validation pending”처럼 곧 완료될 것처럼 쓰지 않습니다.
8. 확인되지 않은 DOI, PMID, cohort size, p-value, software version을 생성하지 않습니다.

## 생성되는 작업물

Claude는 실행 시 `manuscript_v8/review_runs/<timestamp>/` 아래에 다음을 남기도록 설계되어 있습니다.

```text
00_inventory.md
01_claim_evidence_ledger.tsv
02_statistics_audit.md
03_clinical_claims_audit.md
04_figure_audit.md
05_citation_audit.md
06_reproducibility_audit.md
07_editorial_red_team.md
08_change_log.md
09_submission_readiness.md
BLOCKERS.md
```

## 사용 모델

정교한 repair와 red-team은 가능한 가장 강한 reasoning model로 실행하십시오. 단순 inventory나 filename check는 저비용 model로 분리해도 됩니다. 모델 성능보다 중요한 것은 **서로 독립된 context에서 감사한 뒤, 마지막에 합치는 구조**입니다.
