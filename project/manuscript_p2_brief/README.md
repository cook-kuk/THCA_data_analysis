---
title: "Paper 2 advisor brief — directory README"
date: 2026-05-01
purpose: "advisor 미팅 single-page brief 의 모든 deliverable 와 source 검증 인덱스."
---

# Paper 2 advisor brief — directory README

`Paper 2 (Autoimmune-overlap PTC)` 의 advisor 미팅용 single-page discussion brief + 보조 자료. **HLA / autoimmune-overlap (DM2) 축만** 다룸. Paper 1 (Dark Matter, DM1) 은 별도 — 본 디렉토리 scope 밖.

## 1. Deliverables (10 files)

### Primary (advisor 미팅 메인 자료)

| File | Size | Purpose |
|---|---|---|
| **`paper2_brief.html`** ★ NEW | 102 kB | **Story-driven version** — 5 chapter (Introduction / Evidence / Synthesis / Translation / Decision), 14 sections, 11 figures (6 paper-defining), Mermaid mediation chain, 균일한 Pillar 구조 (Q→방법→결과→해석→한계), 약어 사전 § 0 + inline tooltips. Advisor 미팅 권장. |
| **`paper2_brief.pdf`** ★ NEW | 2.8 MB | 38-page printable / email-able PDF. Chapter 별 divider page 로 시각적 분리. |

### Extended reference (kitchen-sink, deep-dive)

| File | Size | Purpose |
|---|---|---|
| `p2_advisor_discussion.html` | 130 kB | **Extended kitchen-sink version** — 21 figures, 5 tables, WHY 박스, 4-scenario sensitivity 등 모든 detail. Reviewer 답변 / supplementary 작성 시 reference. |
| `p2_advisor_discussion.pdf` | 2.6 MB | 50-page extended PDF. 8 sub-section anchors (s311–s318). |
| `DATA_SOURCES_INDEX.md` | 11 kB | 21 figures + 5 tables → underlying TSV/JSON path mapping. **26/26 source paths verified ✅**. § C.2 verbatim ratio (Methods 99.78%, Discussion 97.81%). § E.0 Figure# ↔ Plotly id mapping table. |
| `COHORT_ACCESS_GUIDE.md` | 7 kB | 5 cohorts (TCGA/GSE286332/GSE213647/K2/Chu 2018) raw data accession + download method + processing script + IRB. ~600 MB 처리된 자료로 brief 전체 reproduce 가능. |
| `PROMPT_DECISION_LOG.md` | 12 kB | Audit log — v4 audit 통합 작업 중 scope 오판단 → HLA-only 정정 trace + 16-row propagation checklist |
| `CHANGELOG.md` | 5 kB | Version history — pre-2026-05-01 baseline + 2026-05-01 자율 infra batch + companion infra changes |
| `POST_COMMIT_STATUS.md` | 5 kB | Post-commit infra audit — 4 commits 자율 작업 통합 + pillar 폴더 reorganization 발견 + working tree clean verification |
| `REFERENCES_BIB_AUDIT.md` | 3 kB | manuscript_v8/03_intro_references.bib 의 6/17 incomplete entries 보고. Chu 2018 정확 entry 제안 (voice-protected 영역, content 변경 X) |
| `README.md` | this file | Directory navigation |

## 2. Quick start

### Browser (interactive Plotly)
```bash
open http://40.82.129.113:8012/manuscript_p2_brief/p2_advisor_discussion.html
```

### PDF (offline / email)
```bash
xdg-open project/manuscript_p2_brief/p2_advisor_discussion.pdf
```

### PDF 재생성 (브리프 갱신 후)
```bash
CHROME=/home/seungho/.cache/ms-playwright/chromium-1208/chrome-linux64/chrome
$CHROME --headless=new --disable-gpu --no-sandbox --hide-scrollbars \
  --virtual-time-budget=15000 --run-all-compositor-stages-before-draw \
  --print-to-pdf="project/manuscript_p2_brief/p2_advisor_discussion.pdf" \
  --print-to-pdf-no-header --no-pdf-header-footer \
  http://40.82.129.113:8012/manuscript_p2_brief/p2_advisor_discussion.html
```

## 3. Brief 구성 — 5 Pillars + Mediation

| § | Pillar | Status | 핵심 figure |
|---|---|---|---|
| 3.1 | I — HLA enrichment hypothesis | ★ STRONG (격상 2026-05-03) | Fig 2 ★★★ 6-allele forest (Korean PTC vs Chu ctrl/GD), Fig 2b/c/d 신규 |
| 3.2 | II — In-cohort mechanism (GSE286332) | STRONG | Fig 6 ★★★ HLA-II Cohen d=+3.65 |
| 3.3 | III — Cross-cohort generalization | STRONG | Fig 9 ★★★ TCGA OR=0.20, p=6.4e−10 |
| 3.4 | IV — Antigen-driven specificity | STRONG | Fig 12 ★★★ TLS d=+1.96, IGHV clonality d=+2.09 |
| 3.5 | V — DM1 sub-B = NBNR cluster | PARTIAL | Fig 14 ★★★ 96% mutation-negative |
| 3.6 | Mediation (causal arrow) | STRONG | Fig 16 ★★★ HLA-II 140% mediation |

**현재 status**: 5 pillars 중 4 STRONG + 1 PARTIAL. Cell Rep Med (1순위) → JCI Insight (2순위) → Nat Commun (reach) realistic.

## 4. 2026-05-03 핵심 신규 사항 (Pillar I 격상)

- **Citation correction**: "Chen 2018" → **Chu X et al. 2018** *J Med Genet* 55(10):685–692 (PMC 6161647, doi:10.1136/jmedgenet-2017-105146)
- **Pillar I 정량화**: Korean PTC pool n=874 vs Chu 2018 Han Chinese GD n=1,468 / ctrl n=1,490 random-effects DerSimonian-Laird forest meta
- **6 alleles direction-consistent**: 4 risk (A*02:07, B*46:01, C*01:02, **DPB1*05:01**) + 2 protective (DQB1*02:01, DRB1*07:01)
- **Top signal**: DPB1*05:01 53.2% Korean vs 31.3% Chu ctrl, OR 2.50 [2.10, 2.97], pooled OR 2.16 [1.65, 2.83]
- **Robustness**: 4-scenario sensitivity (전체 / GSE286332 제외 / Lee only / K2 only) 모두 52–56% 안정
- **Sub-cohort heterogeneity**: DPB1*05:01 Cochran I²=0% across K2 / Lee / GSE286332-PTC (가장 robust)

→ Cell Rep Med reach 의 reviewer 저항 가장 큰 위험 (Pillar I 정량 부재) 가 제거됨.

## 5. 다른 source 와 boundary

| Boundary | Paper 1 (Dark Matter) | Paper 2 (이 brief) |
|---|---|---|
| Cluster | DM1 (BRAF/RAS-negative) | DM2 (autoimmune-overlap) |
| Mechanism | Fusion (RET/NTRK/ALK) + epigenetic silencing | HLA-II / IFN-γ / antigen-driven B cell |
| Key finding | DM1 76.8% fusion+, TPO methylation d=2.30 | HLA-II d=+3.65, mediation 140%, TCGA OR=0.20 |
| Manuscript dir | `project/manuscript_v8/` | (TBD — Paper 2 writing 6/11 시작) |
| Brief dir | (없음) | **`project/manuscript_p2_brief/`** ← here |

→ v4 audit 의 R5 fusion / methylation / 3-layer DM1 model 은 **Paper 1 영역, 본 brief 에 들어오지 않음**. 자세한 분류는 `PROMPT_DECISION_LOG.md` § 3 참조.

## 6. Companion 자료

- **Pillar I forest analysis**: `../results/p2_pillar1_forest/` — `PILLAR1_FOREST_SUMMARY.md`, `methods_paragraph.md`, `discussion_paragraph.md` (Cell Press paste-ready)
- **Deprecated 자료**: `../results/d4p1_panasian_meta/DEPRECATED.md` — Chen 2018 wrong citation 사용 금지 표
- **Memory**: `~/.claude/.../memory/v17_paper2_pillar1_forest_strong.md` — 2026-05-03 STRONG 격상 fact

## 7. Next steps (advisor 미팅 후)

1. § 7 Discussion Q01–Q12 advisor 답변 받음
2. Q3 (mini-index calibration), Q9 (citation propagation) 의 advisor 결정 → action items
3. 분당 outreach 응답 대기 (Scenario A vs B 분기)
4. 6/11 W6 end → Paper 2 본격 writing 시작 (`project/manuscript_p2/`)

---

*Generated 2026-05-01. 마지막 brief 갱신: Pillar I STRONG (2026-05-03 자료 통합). Author: Seungho Cook (kukshomr@gmail.com).*
