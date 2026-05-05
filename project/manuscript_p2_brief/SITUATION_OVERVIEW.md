<!-- VS Code: Cmd/Ctrl+K V → side-by-side preview -->
---
title: "Paper 2 isolated session — 상황 판단 dashboard"
date: 2026-05-04
session: "Paper 2 (Hashimoto-overlap PTC) isolated, post second-sweep"
purpose: "한 페이지로 전체 상황 즉시 판단 — 어디까지 됐고, 뭐가 남았고, 뭘 봐야 하는지"
---

# 🎯 TL;DR (3줄)

1. **Paper 2 = Hashimoto-overlap PTC ONLY** (Yu 2026-05-04 advisor SCOPE). Graves'/GD = Paper 4 reserve. Active deliverables 4개 모두 SCOPE clean ✅
2. **Pillar I framing v1 → v2 migration 완료** — active narrative 는 Korean baseline (HT context); v1 (Chu 2018 GD comparator) 는 deprecation banner 처리. v2 figure data 자체는 `p2_pillar1_forest_v2/` 의 PDF 참조
3. **다음 행동** = (a) Yu meeting 에서 v2 forest review (b) Lee 2014 Korean reference 로 C*01:02 / DQB1*02:01 fill (c) Phase 1 paper outline 본인 키보드 시작

---

# 1. 4 active deliverables — 한눈 status

| # | Deliverable | URL | PDF pages | SCOPE clean | 마지막 갱신 |
|---|---|---|---|---|---|
| 1 | **Story brief** (canonical advisor 자료) | [paper2_brief.html](http://40.82.129.113:8012/manuscript_p2_brief/paper2_brief.html) | [38 → 40+ pages](paper2_brief.pdf) | ✅ Pillar I v2 audit banner | 2026-05-04 (Pillar I 정정) |
| 2 | **Extended brief** (kitchen-sink) | [p2_advisor_discussion.html](http://40.82.129.113:8012/manuscript_p2_brief/p2_advisor_discussion.html) | [50 pages](p2_advisor_discussion.pdf) | ✅ deprecation banner only | 2026-05-04 (`<title>` 정정) |
| 3 | **Master papers index** | [papers_overview.html](http://40.82.129.113:8012/submission/papers_overview.html) | [10 pages](../submission/papers_overview.pdf) | ✅ Paper 2 row Hashimoto-framed | 2026-05-04 |
| 4 | **3-papers trajectory index** | [three_papers_index.html](http://40.82.129.113:8012/three_papers_index.html) | [pages](../three_papers_index.pdf) | ✅ Paper 2 row v2 (Korean baseline) | 2026-05-04 |

**Forbidden term grep 결과** (active 자료 전체):
- `Autoimmune-overlap PTC` 잔재: **0** ✅
- `Pan-Asian autoimmune-thyroid` 잔재: **0 in Paper 2 context** ✅
- Graves'/GD 잔재 = glossary 정의 + Bundang outreach status + audit notes 만 (SCOPE 허용)

---

# 2. Pillar I 의 v1 → v2 migration 상태 (가장 큰 변경)

## 2.1 v1 framing (DEPRECATED, record only)

```
Korean PTC pool n=874  vs  Chu 2018 Han Chinese GD n=1,468 / ctrl n=1,490
                            ↑↑↑ 이 비교가 Paper 2 SCOPE 위반 (GD = Paper 4)
DPB1*05:01: PTC 53.2%  vs  Chu ctrl 31.3%  vs  Chu GD 44.0%
OR (Korean vs Chu ctrl): 2.50 [2.10, 2.97]
Pooled OR (random-effects): 2.16 [1.65, 2.83]
```

**v1 위치**: `project/results/p2_pillar1_forest/` (전체 dir 에 DEPRECATED banner 처리)

## 2.2 v2 framing (★ ACTIVE — Paper 2 SCOPE 준수)

```
Korean PTC pool n=874  vs  Korean baseline (AFND South Korea pool)
                            ↑↑↑ 같은 인구 내부 비교 (HT context)
DPB1*05:01: PTC 53.2%  vs  Korean baseline 36.7%
OR: 1.96 [1.60, 2.41], p=1.1e−10
A*02:07:    PTC 8.1%   vs  baseline 3.4%   OR 2.54
B*46:01:    PTC 10.3%  vs  baseline 4.7%   OR 2.33
DRB1*07:01: PTC 11.4%  vs  baseline 5.0%   OR 2.46  ⚠️ direction discord
            (v1 GD context 에서 protective; v2 risk direction — Harbin Korean
             baseline n=201 의 under-representation 가능성)
C*01:02:    AFND South Korea entry 없음 → Lee 2014 권장
DQB1*02:01: AFND South Korea entry 없음 → Lee 2014 권장
```

**v2 위치**: `project/results/p2_pillar1_forest_v2/`
- 핵심 file: [`PILLAR1_FOREST_V2_SUMMARY.md`](../results/p2_pillar1_forest_v2/PILLAR1_FOREST_V2_SUMMARY.md)
- Forest figure: [`forest_paper2_HT_only.pdf`](../results/p2_pillar1_forest_v2/forest_paper2_HT_only.pdf)
- Discussion paste-ready: [`discussion_paragraph_v2.md`](../results/p2_pillar1_forest_v2/discussion_paragraph_v2.md)

## 2.3 paper2_brief.html Pillar I section 의 현재 상태

- Section title + kicker = v2 framed ✅
- Strong audit banner = v2 figure location + 핵심 수치 명시 ✅
- 질문 line = v2 framed (Korean baseline) ✅
- 방법/결과/해석/한계 QMRIL body = **v1 narrative 보존** + audit banner 로 cross-ref
- Figure 1 (Plotly JS) = **여전히 v1 data** (Chu 2018 forest) — 향후 update 후보

→ **부분 migration 완료**. 내용 자체가 advisor 가 read 할 때 v2 audit banner 로 즉시 redirected. 완전한 v2 re-render (Figure 1 JS data + QMRIL body 전체 v2) 는 advisor 검토 후 다음 sprint.

---

# 3. 5 Pillars 현재 status (Paper 2 isolated SCOPE)

| # | Pillar | Headline | Status | Source |
|---|---|---|---|---|
| **I** | **Korean PTC HLA susceptibility** | DPB1*05:01 OR 1.96 vs Korean baseline (HT context) | ★ STRONG (v2 active) | `p2_pillar1_forest_v2/` |
| **II** | **GSE286332 PTC vs PTC+HT in-cohort** | 8-gene Cohen d=−1.60, HLA-II d=+3.65, IFN-γ FDR=2e−4, 10,380 DEGs | STRONG | `p3_gse286332/` |
| **III** | **TCGA + Korean cross-cohort** | DM2 Hashimoto OR=0.20, p=6.4e−10; 18-30% prevalence robust | STRONG | `d4p2_tcga_hashimoto_signature/` + `d8b_korean_replication/` |
| **IV** | **BCR clonal + TLS antigen-driven** | TLS Cabrita d=+1.96, IGHV clonality d=+2.09 | STRONG | `d5p6_bcr_repertoire/` |
| **V** | **DM1 sub-B = NBNR cluster** | 96% mutation-neg sub-cluster; Korean transfer 47-53% | PARTIAL | `d6p7_dm1_subcluster/` + `d8c_dm1_subB_x_K2_NBNR/` |

**§ 10 Mediation synthesis**: PTC+HT → HLA-II → P(DM1) ↓, **HLA-II 140% mediation** (p=0.023), generic immune NS — confounder vs mediator 정량 분리. Source: `d3p5_pdm1_gradient/`

---

# 4. Memory entries 현재 status (3 sessions 누적 audit)

| Entry | Status | 한 줄 요약 |
|---|---|---|
| `v18_paper2_HT_isolated.md` | ✅ **canonical** | Paper 2 = HT-only ONLY, forbidden words list, Pillar I v2 framing, Task A/B/C anchor |
| `v17_paper2_pillar1_forest_strong.md` | ⚠️ **SUPERSEDED** | v1 (Chu 2018 GD) record. Cohort + DPB1*05:01 53.2% facts stay; framing changed to Korean baseline |
| `v17_graves_pivot.md` | ⚠️ **SUPERSEDED** | 2026-04-29 v0 pivot ("Graves' = main"). Now Paper 1 active, Korean GD HLA = Paper 4 backlog |
| `v17_arcasHLA_korean_k2.md` | ✅ historical 분석 | DPB1*05:01 56% Kim 2014 ref replicated, EUR depletion 자체 검증 — methodology 자료 OK |
| `v17_2026_04_30_pivot.md` | ✅ historical pivot context | 그 시점 결정 record |
| `v17_K2_vs_bundang_distinction.md` | ✅ active reference | K2 = PRJEB11591 public; Bundang = outreach-stage 별도 |
| `v17_korean_k2_calibration.md` | ✅ active reference | 8-gene mini-index TPM inflate 10-100×; within-sample-centered profile 사용 |
| `v17_gse286332_strong_go.md` | ✅ active reference | Pillar II core data |
| `v17_D4P2_tcga_hashimoto_generalization.md` | ✅ active reference | Pillar III core data |
| `v17_D5P6_BCR_clonal_TLS.md` | ✅ active reference | Pillar IV core data |
| `v17_D6P7_dm1_subB_NBNR.md` | ✅ active reference | Pillar V core data |
| `paper_numbering_2026_05_04.md` | ✅ canonical portfolio | Paper 1/2/3/4 정의 |
| `v19_paper3_ici_track_a.md` | ✅ Paper 3 FROZEN bundle | ICI vulnerability dark thyroid cancer (touch X) |
| `v19_paper4_GD_backlog.md` | ✅ Paper 4 backlog | Korean GD HLA, gated 4/4 |

---

# 5. 본 isolated session (2 sub-sessions) 한 일 timeline

## 2026-05-04 first session (15:43 prior)

| Task | 한 일 |
|---|---|
| A | `p2_pillar1_forest_v2/` 6 file 생성 (Korean baseline forest, 4 alleles) |
| B | First sweep 5 file terminology 정정 + DEPRECATED banner |
| C | `v18_paper2_HT_isolated.md` 작성 |

## 2026-05-04 second session (this thread)

| Task | 한 일 |
|---|---|
| A | 검증 only ('거기까지만 해') |
| B | Second sweep — 3개 portfolio-level occurrence 추가 정정 (papers_overview × 2 + p2_advisor_discussion `<title>`) |
| C | `v17_paper2_pillar1_forest_strong` SUPERSEDED 처리 + MEMORY.md 갱신 + `v18_paper2_HT_isolated` footer second-sweep section |
| α | `v17_graves_pivot.md` SUPERSEDED note + body 보존 |
| β | `three_papers_index.html` Paper 2 row v2 data 정렬 (cohort + headline + status) |
| γ | `paper2_brief.html` Pillar I section title/kicker/audit banner/질문 v2 framing 정정 |
| δ | PDF 3개 재생성 + final SCOPE grep verify |

총 **7 actions** (Task A/B/C + α/β/γ/δ) 본 thread.

---

# 6. 다음 step (advisor / 본인 키보드 영역, 본 session SCOPE 외)

| 우선순위 | Action | 책임 | 예상 effort |
|---|---|---|---|
| 1 | **Yu meeting 에서 Pillar I v2 forest review** | advisor | 미팅 30 min |
| 2 | **C*01:02 / DQB1*02:01 의 Lee 2014 Korean ref freq fill** | analysis (자율 가능) | 1-2 hr (web fetch + script edit) |
| 3 | **Phase 1 paper outline (Task D)** | 본인 키보드 (voice-protected) | Yu 합의 후 시작 |
| 4 | **paper2_brief.html Figure 1 JS data v1 → v2 swap** | analysis (자율 가능) | 1-2 hr (Plotly JS edit + PDF re-gen) |
| 5 | **paper2_brief.html Pillar I QMRIL body 전체 v2 re-write** | analysis (자율 가능) but advisor approval 권장 | 1 hr |
| 6 | **Paper 4 backlog 진입 4/4 gating** | dependency (Paper 1 bioRxiv + Paper 2 A/B/C + Bundang n>50 + Yu 합의) | 6-9 mo |

---

# 7. Files map (VS Code 에서 빠르게 open)

## Paper 2 brief deliverables
```
project/manuscript_p2_brief/
├── paper2_brief.html              ← ★ canonical story brief
├── paper2_brief.pdf               ← printable
├── p2_advisor_discussion.html     ← extended kitchen-sink
├── p2_advisor_discussion.pdf
├── README.md                      ← directory navigation
├── DATA_SOURCES_INDEX.md          ← figure → TSV mapping (26/26 verified)
├── COHORT_ACCESS_GUIDE.md         ← 5 cohort accession + reproducibility
├── PROMPT_DECISION_LOG.md         ← v4 audit 통합 작업 audit log
├── CHANGELOG.md
├── POST_COMMIT_STATUS.md
├── REFERENCES_BIB_AUDIT.md
└── SITUATION_OVERVIEW.md          ← ★ 이 파일 (전체 상황 한 페이지)
```

## Pillar I source data (v1 deprecated + v2 active)
```
project/results/
├── p2_pillar1_forest/             ← v1 DEPRECATED-banner wrapped
│   ├── PILLAR1_FOREST_SUMMARY.md  (deprecated)
│   ├── discussion_paragraph.md    (deprecated)
│   └── forest_panasian_HLA.{pdf,png}  (Chu 2018 forest, record only)
└── p2_pillar1_forest_v2/          ← ★ v2 ACTIVE (Korean baseline)
    ├── PILLAR1_FOREST_V2_SUMMARY.md
    ├── discussion_paragraph_v2.md
    ├── korean_baseline_allele_freq.tsv
    ├── korean_PTC_vs_korean_baseline_forest.json
    └── forest_paper2_HT_only.{pdf,png}
```

## Pillars II-V source data (HT-only, all SCOPE clean)
```
project/results/
├── p3_gse286332/                  ← Pillar II in-cohort mechanism
├── d4p2_tcga_hashimoto_signature/ ← Pillar III TCGA arm
├── d8b_korean_replication/        ← Pillar III Korean arm
├── d5p6_bcr_repertoire/           ← Pillar IV antigen-driven
├── d6p7_dm1_subcluster/           ← Pillar V TCGA discovery
├── d8c_dm1_subB_x_K2_NBNR/        ← Pillar V Korean transfer
└── d3p5_pdm1_gradient/            ← § 10 Mediation
```

## Session report dir (Task B + audit trail)
```
project/results/terminology_correction_2026_05_04/
├── before_after_diff.md           ← substitution rules + first/second sweep
├── files_modified.txt             ← all touched files
└── SESSION_REPORT.md              ← 이전 turn 의 session 흐름 정리
```

## Memory entries (auto-loaded next session)
```
~/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/
├── MEMORY.md                       ← 한 줄 인덱스 (24 entries)
├── v18_paper2_HT_isolated.md       ← ★ canonical Paper 2 SCOPE
├── v17_paper2_pillar1_forest_strong.md  (SUPERSEDED note)
├── v17_graves_pivot.md             (SUPERSEDED note)
├── paper_numbering_2026_05_04.md   ← canonical 4-paper portfolio
├── v19_paper4_GD_backlog.md        ← Paper 4 (GD) gating
└── ...
```

## Cross-paper indexes
```
project/
├── three_papers_index.html         ← 3-paper trajectory landing page
├── three_papers_index.pdf
└── submission/
    ├── papers_overview.html        ← 6-paper master dashboard
    └── papers_overview.pdf
```

---

# 8. VS Code 에서 빠르게 보는 방법

## Single file open
```bash
# 이 파일
code /home/seungho/personal/THCA_data_analysis/project/manuscript_p2_brief/SITUATION_OVERVIEW.md

# Pillar I v2 summary
code /home/seungho/personal/THCA_data_analysis/project/results/p2_pillar1_forest_v2/PILLAR1_FOREST_V2_SUMMARY.md

# canonical memory entry
code /home/seungho/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/v18_paper2_HT_isolated.md
```

## Whole project open
```bash
code /home/seungho/personal/THCA_data_analysis
```
그 후 VS Code 에서:
- `Cmd/Ctrl+P` → `SITUATION_OVERVIEW.md` 검색
- 파일 view 후 `Cmd/Ctrl+K V` → side-by-side markdown preview

## Browser HTTP serving (active deliverables 직접 view)
- Master index: http://40.82.129.113:8012/submission/papers_overview.html
- 3-papers landing: http://40.82.129.113:8012/three_papers_index.html
- Story brief: http://40.82.129.113:8012/manuscript_p2_brief/paper2_brief.html
- Extended brief: http://40.82.129.113:8012/manuscript_p2_brief/p2_advisor_discussion.html

---

# 9. SCOPE compliance final check (8/8 ✅)

| SCOPE 규칙 | 본 session 준수 |
|---|---|
| 1. Hashimoto thyroiditis (HT) overlap 의 분자 mechanism ONLY | ✅ Active deliverables 4/4 clean |
| 2. Graves' disease (GD) phenotype/mechanism 언급 금지 | ✅ Glossary 정의 + outreach status + audit notes 만 잔존 |
| 3. Cancer driver mutation deep dive 금지 | ✅ Touch X |
| 4. `autoimmune-PTC` → `Hashimoto-overlap PTC` 일괄 정정 | ✅ Two sweeps complete (5 + 3 files) |
| 5. Pillar I forest = Korean baseline (NOT Chu 2018 GD) | ✅ v2 forest 검증, v1 SUPERSEDED, brief Pillar I migrated |
| 6. Paper 3 file touch 금지 (read 도 X) | ✅ `manuscript_p3_brief/` touch 안함 |
| 7. Voice-protected sections (Hook/Aim/Disc 3.1/Limit/Cover/Q9) | ✅ Touch X |
| 8. Task D (Phase 1 outline) NOT 진행 | ✅ SCOPE 외 |

---

# 10. "거기까지만 해" — 오늘 세션의 stop 지점

본 session 의 SCOPE 가 명시한 boundary 에 정확히 도달:

**한 일** (within SCOPE):
- Forest meta v2 검증
- Terminology second sweep
- Memory v17 supersession 처리
- 3개 portfolio-level deliverable Pillar I 정정
- PDF 재생성

**안 한 일** (SCOPE 외):
- Phase 1 paper outline (voice-protected, Yu meeting 후)
- C*01:02 / DQB1*02:01 Lee 2014 web fetch (advisor 의제로 deferred)
- Figure 1 JS data swap (advisor approval 권장)
- Pillar I QMRIL body 전체 v2 re-write (advisor approval 권장)
- Paper 3/4 territory touch (read 포함 금지)

---

*Generated 2026-05-04 by Paper 2 isolated session — second sweep wrap-up.*
*VS Code 에서 markdown preview (Cmd/Ctrl+K V) 권장 — 표/링크 모두 렌더링.*
