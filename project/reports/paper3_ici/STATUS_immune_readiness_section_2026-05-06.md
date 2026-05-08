# Immune Readiness ≠ Response Prediction — section status

**작성일:** 2026-05-06 KST
**대상 페이지:** THCA Hub `papers_hub_2026_05_04/index.html`
**Evidence layer:** Paper 3 ICI — Track B-lite Analysis Report

## URL

라이브 (포트 8012 · `serve_secure.py`):

```
http://40.82.129.113:8012/papers_hub_2026_05_04/index.html#immune-readiness
```

내부 LAN:

```
http://127.0.0.1:8012/papers_hub_2026_05_04/index.html#immune-readiness
```

응답: `HTTP/1.0 200 OK · 84,024 bytes` (편집 전 63,564 → 약 +20KB 섹션 추가).

> 참고: `http://40.82.129.113:80/papers_hub_2026_05_04/index.html` (papers-web.service 의 별도 복사본, `/var/www/papers/...`) 은 의도적으로 손대지 않았습니다. 사용자 정책상 hands-off. 동기화가 필요하면 별도 deploy step에서 처리.

## 편집한 파일 (단 1개)

```
/home/seungho/personal/THCA_data_analysis/project/papers_hub_2026_05_04/index.html
```

- 백업: `index.html.bak-immune-readiness-2026-05-06` (같은 디렉토리, 63,564 bytes)
- 삽입 위치: "Today · 2026-05-06 · professor-facing strengthening" 카드 직후 + IDEKER-style subplot 직전
- 앵커: `#immune-readiness`

## 검증

| 체크 | 결과 |
|---|---|
| `curl -I` 8012 | `HTTP 200 · 84024B` |
| HTML tag balance (section / div / style / script / head / body) | 모두 delta=0 |
| 새 섹션 핵심 string 포함 (`Cohen's d  +2.53`, `69.4%`, `3 / 7 modules`, `ATC Immune Vulnerability Map`, `DIAL — Direction-Invariance`) | 10/10 hit |
| `/var/www/papers/...` 복사본 md5 | 변경 없음 (b0beed5c…) |
| 기존 페이지 sections | 보존 (hero, dashboard, today's updates, IDEKER, portfolio cards, drug platform, process docs) |

## 새 섹션 구조 (id `immune-readiness`)

```
[eyebrow]   Paper 3 ICI · Track B-lite · 2026-05-06 · ThyroidVax / Lumenix evidence layer
[title]     Immune Readiness ≠ Response Prediction
[subtitle]  4 thyroid cohorts · n=415 · DIAL 3/7 sign · readiness ≠ prediction
[meta-grid] evidence · comparison · audit · stack

▣ Visual module 1 — ATC Immune Vulnerability Map  (3-col grid, 7 cards)
  ├─ Myeloid suppression ↑   d = +2.53  (4/4 sign · large)
  ├─ Checkpoint exhaustion ↑ d = +1.01  (4/4)
  ├─ IFNG / T-cell inflamed ↑ d = +0.55  (4/4)
  ├─ HLA-II ↑                 d = +0.68  (3/4)
  ├─ HLA-I ↑                  d = +0.58  (4/4)
  ├─ Thyroid differentiation ↓ d = −2.51 (4/4 · large)
  └─ TLS / CXCL13 ~ not uniformly ↑  d = −0.15  (3/4)  · Paper 2 cross-boundary guard

▣ Two-axis immune map (2-col)
  ├─ PC1  69.4%  pan-immune activation (purple)
  └─ PC2  14.2%  dedifferentiation (blue) ·  83.6% 분산 설명

▣ Visual module 2 — DIAL · Direction-Invariance Audit Layer  (5-step flow)
  Pan-cancer ICI cohorts
    → Module direction test (7 modules × responder vs non-responder)
    → Bootstrap (1000) + permutation (1000) sign agreement
    → Thyroid transfer eligibility (≥5/7 · LODO)
    → Hypothesis only / blocked claim  ← orange terminal step
  현재 audit: 3/7 sign (HLA-I / HLA-II / TLS only) · 4 모듈 sign-flip · binomial p=0.77

▣ Claim boundary  (orange "Do NOT claim" | teal "DO claim · research-grade")
  NOT: response prediction / treatment recommendation / 개별화 결정 / K1 kill
  DO : ICI-readiness · immunogenomic vulnerability prioritization
       · vaccine-combination hypothesis · DIAL methodological necessity

▣ Pull quote (teal/orange gradient · Cormorant Garamond italic)
  "Every immune score is treated as a hypothesis until it passes
   cross-cohort direction-invariance, HLA/neoantigen architecture,
   and thyroid-specific validation."

[foot] source path · Track A frozen · Track B unlock 조건 · paper3 + cancer_vaccine_agent 링크
```

## 디자인 토큰

자체적으로 격리된 `.iv-` 프리픽스 CSS (기존 페이지 스타일과 충돌 없음):

| 토큰 | 값 | 용도 |
|---|---|---|
| `--iv-bg` / `--iv-bg-2` / `--iv-bg-3` | `#07101a` → `#0d1422` → `#122036` | 다크 biotech 대시보드 (페이지 hero `Cancer vaccine agent` 버튼과 동일 톤) |
| `--iv-teal` | `#35d39d` | validated analysis signal · DO claim · up arrows |
| `--iv-orange` | `#ff8a3d` | risk / guardrail · DO NOT claim · DIAL block · down arrows |
| `--iv-purple` | `#8c7bff` | immune modules (myeloid / checkpoint / IFNG) |
| `--iv-blue` | `#5fb1ff` | HLA modules · PC2 axis |
| `--iv-ink` / `--iv-ink-soft` / `--iv-ink-dim` | `#e6f1ee` / `#93a8b8` / `#6c7e91` | 본문 / sub / metadata |

## 반응형

- ≤ 1100px → vmap 3→2열, DIAL flow 5→2열, claim 2→1열
- ≤ 720px → 모두 1열, 패딩 축소, 타이틀 26px
- 인쇄 (Ctrl+P) 모드: 기존 페이지의 print-css 가 그대로 적용됨 (다크 배경 → 흰 배경 자동 전환은 별도 처리 안 함; 다크 톤 유지가 의도)

## Claim guard (binding · 페이지에 명시됨)

**Do NOT claim** — thyroid 환자에서 ICI response prediction · 임상 treatment recommendation · 개별화 면역치료 결정 · K1 kill 영역 ("ICI response predictor" 류).

**DO claim · research-grade** — ICI-readiness framing · immunogenomic vulnerability prioritization · vaccine-combination hypothesis generation · DIAL audit methodological necessity.

## Source 정렬

- 본 분석 보고서: `reports/paper3_ici/paper3_ici_track_B_lite_analysis_report.md`
- Track A 디자인 번들: `reports/paper3_ici/PAPER3_ICI_TRACK_A_BUNDLE.md` (chmod 444 frozen + sha256)
- Figure plan: `reports/paper3_ici/paper3_ici_figure_plan.md`
- Dataset registry: `reports/paper3_ici/paper3_ici_dataset_registry.md`

## Track B 게이트 상태 (페이지에 명시됨)

| 게이트 | 상태 |
|---|---|
| G1 — Paper 1 bioRxiv 제출 | OPEN (target 2026-06-13) |
| G2 — Paper 2B Pillar II/III closure | OPEN |
| G3 — literal "Paper 3 Track B 시작" 명령 | ❌ 미발화 |

본 섹션은 Track A frozen 의 **evidence layer 가시화**일 뿐이며 Track B 의 본 분석 (NMF / LOHHLA / NetMHCpan / 5+ cohort 완전 DIAL audit / scRNA atlas) 을 대체하지 않음.

## 다음에 할 만한 것 (선택)

- **`/var/www/papers/...` 미러 동기화** — 현재 papers-web.service 가 서빙하는 80번 포트 페이지에는 이 섹션이 없습니다. 의도적이지만, 동기화하려면 수동 cp 또는 별도 deploy step 필요. 정책상 systemd 유닛은 손대지 않음.
- **Plotly로 PCA biplot 추가** — PC1 vs PC2 위에 8 cohort × n=415 를 그리고 ATC 만 강조. 페이지에 이미 `https://cdn.plot.ly/plotly-2.35.2.min.js` 가 로드돼 있으니 의존성 추가 없음.
- **Paper 3 hub 페이지 (paper3.html) 에서도 동일 섹션 link** — 현재 `Track A frozen` 이라는 status badge 와 자연스러운 짝.
- **인쇄 시 라이트 톤 자동 전환** — `@media print` 안에 `.iv-section{background:#fff;color:#000}` 정도 한 줄.
- **PNG figure embed** — 보고서 옆 `figures_png/` 의 PCA / hclust / DIAL forest 도 섹션 안 figgrid 로 직접 임베드 가능.
