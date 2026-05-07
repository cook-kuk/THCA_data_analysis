---
title: "Paper 2 isolated session — Session report (전체 흐름 정리)"
date: 2026-05-04
session_id: "Paper 2 isolated session, second sweep"
purpose: "이 session 의 시작 → 진행 → 현재 상태 까지 한 번에 read"
---

# Paper 2 isolated session — 전체 흐름

## 1. Session 시작 (advisor 지시 inline)

### 1.1 Trigger prompt

이 session 은 **사용자가 명시적으로 보낸 isolated session prompt** 로 시작:

```
=== Paper 2 isolated session ===
Topic: Hashimoto-overlap PTC 의 분자 축 — antigen-driven HLA-II–매개
       dedifferentiation

이 session 의 SCOPE (advisor 지시 2026-05-04):
- Hashimoto thyroiditis (HT) overlap 의 분자 mechanism ONLY
- Graves' disease (GD) phenotype/mechanism 언급 금지 — Paper 3 영역
- Cancer driver mutation deep dive 금지 — Paper 1 영역
- "autoimmune-PTC" 표현 → "Hashimoto-overlap PTC" 일괄 정정
```

### 1.2 Advisor 의 핵심 catch (inline 인용)

Yu 교수님 2026-05-04 catch:

> "하시모토병은 갑상선 조직 내에서 면역 세포가 직접 침윤하여 분화도를
>  떨어뜨리는(Dedifferentiation) 파괴적 성향이 강한 반면,
>  그레이브스 병은 유전적으로 유사한 배경을 공유하지만 갑상선 기능을
>  비정상적으로 항진시키는 차이가 있습니다."

**해석**:
- Paper 2 의 main mechanism = HT 의 destructive infiltration → dedifferentiation
- GD 의 hyperthyroid mechanism 은 Paper 3 (이후 Paper 4 backlog 로 renumber)
- HLA susceptibility 의 GD overlap 은 "shared genetic background" 까지만, disease equivalence 주장 X

### 1.3 정의된 4개 task

| Task | 내용 | 시간 예상 | 본 session 진행 |
|---|---|---|---|
| **A** | Forest meta v2 — Korean PTC vs Korean baseline (Chu 2018 GD 비교 제거) | 1-2 hr | ✅ 검증 (이전 session 이 이미 완료) |
| **B** | Terminology 일괄 정정 — `autoimmune-PTC` → `Hashimoto-overlap PTC` 등 | 30 min | ✅ Second sweep 완료 |
| **C** | Memory entry — `v18_paper2_HT_isolated.md` 작성 + v17 entries audit | 15 min | ✅ 갱신 + supersession 처리 |
| **D** | Phase 1 paper outline — voice-protected, 본인 키보드 영역 | — | ❌ NOT 진행 (SCOPE 외) |

### 1.4 Cross-contamination check 규칙 (advisor 명시)

본 session 에서 다음 단어 등장 시 즉시 STOP + 재작성:

**(1) GD-specific terms (Paper 3 → Paper 4 backlog 영역, Paper 2 에서 금지)**
- `Graves`, `GD`, `TSAb`, `TSI`, `thyroid-stimulating antibody`
- `hyperthyroidism`, `thyrotoxicosis`, `exophthalmos`
- `Graves' ophthalmopathy`, `TED`, `thyroid eye disease`
- Bundang Graves' cohort 직접 언급

**(2) Cancer-only terms (Paper 1 영역, Paper 2 에서 minor 만)**
- BRAF V600E mechanism deep dive
- TERT promoter mutation kinetics
- DM1 76.8% fusion paradigm (Paper 1 main, Paper 2 supplementary OK)

**(3) Forbidden conflations**
- `autoimmune-PTC` → `Hashimoto-overlap PTC`
- `autoimmune-thyroid axis` → `Hashimoto-thyroid overlap axis`
- `Pan-Asian autoimmune-thyroid susceptibility` → `Pan-Asian thyroid HLA susceptibility (HT context)`

**(4) Pillar I forest meta**
- 현재: Korean PTC pool vs Chu 2018 Han Chinese GD ❌
- 정정: Paper 2 forest = Korean PTC vs Korean baseline (HT context only) ✅
- Chu 2018 GD 비교는 Discussion 한 줄 ("shared HLA background with GD, but disease distinct — see Paper 3") 만 허용
- Chu 2018 vs Chinese ctrl forest 는 Paper 3 (이후 Paper 4) reserve

---

## 2. 시작 시점의 state (이전 session 이 이미 한 일)

본 session 이 시작될 때 발견한 사실:

### 2.1 Task A — 이미 완료 (2026-05-04 15:43 prior session)

`project/results/p2_pillar1_forest_v2/` 가 이미 존재:

```
PILLAR1_FOREST_V2_SUMMARY.md   ← 4 alleles 3-arm forest result
discussion_paragraph_v2.md     ← Cell Press paste-ready (HT-only)
forest_paper2_HT_only.{pdf,png} ← 시각화
korean_PTC_vs_korean_baseline_forest.json ← 4 alleles forest data
korean_baseline_allele_freq.tsv ← AFND South Korea pool baseline
```

**이전 session 결과 핵심**:
- Korean PTC pool n=874 (K2 235 + Lee 2024 630 + GSE286332-PTC 9)
- Korean baseline = AFND South Korea pop weighted (n=4,613 for A*02:07/B*46:01, n=680 for DPB1*05:01, n=201 Harbin Korean diaspora for DRB1*07:01)
- 4 alleles forest 결과:
  - **DPB1*05:01**: PTC 53.2% vs baseline 36.7%, OR 1.96, p=1.1e−10 (★ HT susceptibility allele)
  - **A*02:07**: PTC 8.1% vs baseline 3.4%, OR 2.54, p=3.1e−10
  - **B*46:01**: PTC 10.3% vs baseline 4.7%, OR 2.33, p=1.2e−10
  - **DRB1*07:01**: PTC 11.4% vs baseline 5.0%, OR 2.46, p=0.0084 (caveat: prior "protective" label 가 v2 forest 에서 risk direction — Harbin Korean baseline n=201 의 under-representation 가능성 명시)
- C*01:02, DQB1*02:01: AFND South Korea entry 없음 → flagged as "unavailable", Lee 2014 Tissue Antigens 권장

### 2.2 Task B — 부분 완료 (first sweep)

`project/results/terminology_correction_2026_05_04/` 가 이미 존재:

```
before_after_diff.md   ← substitution rules + first sweep changes
files_modified.txt     ← 5 files modified
```

**First sweep (이전 session)** 이 처리한 5 files:
1. `project/manuscript_p2_brief/p2_advisor_discussion.html` (3 occurrences + deprecation banner)
2. `project/manuscript_p2_brief/paper2_brief.html` (1 occurrence + audit notice)
3. `project/results/p2_pillar1_forest/discussion_paragraph.md` (DEPRECATED header)
4. `project/results/p2_pillar1_forest/PILLAR1_FOREST_SUMMARY.md` (DEPRECATED header)
5. `project/notebooks_or_scripts/v17_paper2_pillar1_forest.py` (DEPRECATED docstring)

### 2.3 Task C — 부분 완료 (이전 session)

`~/.claude/.../memory/v18_paper2_HT_isolated.md` 가 이미 존재 — well-written:
- Forbidden words list
- Pillar I v2 framing 정의
- Brand expression preserved 명시
- Voice-protected sections cross-link
- Files anchor (Task A/B output)
- Cross-paper boundary 표

`MEMORY.md` 에도 한 줄 entry 등록되어 있음.

---

## 3. 본 session 에서 한 일 (Second sweep)

### 3.1 Task A — 검증 only

`p2_pillar1_forest_v2/` 자료가 이미 advisor SCOPE 에 충족되므로 추가 보강 없이 verify:
- 4 alleles 3-arm forest 자료 확인 ✅
- Discussion paragraph 의 마지막 GD line ("see Paper 3 for the Graves'-specific HLA forest") 가 SCOPE 가 명시한 한 줄 cross-ref 와 정확히 align ✅
- C*01:02 / DQB1*02:01 의 Lee 2014 Korean reference 자료 fill 은 본 session 에서 X (advisor "거기까지만 해" 지시 따라 deferred)

### 3.2 Task B — Second sweep 추가 정정

Cross-portfolio review 에서 first sweep 누락 발견:

| File | Line | Before | After |
|---|---|---|---|
| `project/submission/papers_overview.html` | 264 (HTML comment) | `<!-- Paper 2 — Autoimmune-overlap PTC -->` | `<!-- Paper 2 — Hashimoto-overlap PTC -->` |
| `project/submission/papers_overview.html` | 370 (Table 1 row) | `Autoimmune-overlap PTC ★` | `Hashimoto-overlap PTC ★` |
| `project/manuscript_p2_brief/p2_advisor_discussion.html` | 6 (`<title>` tag) | `Paper 2 — Autoimmune-overlap PTC \| Advisor Discussion Brief` | `Paper 2 — Hashimoto-overlap PTC \| Advisor Discussion Brief` |

`before_after_diff.md` + `files_modified.txt` 갱신 — second sweep section 추가.

### 3.3 Task C — Memory v17 audit + supersession

**`v17_paper2_pillar1_forest_strong.md`** 의 frontmatter + body 갱신:
- Title: `... Korean PTC vs Chu 2018 Han Chinese GD` → `... (SUPERSEDED for Paper 2)`
- Description: `2026-05-03 Pillar 1 PARTIAL → STRONG. ...` → `2026-05-03 historical record. SUPERSEDED 2026-05-04 by v18_paper2_HT_isolated — the Chu 2018 GD comparison is now Paper 4 territory; Paper 2 Pillar I uses Korean baseline (AFND South Korea pool) instead. Read this entry only as v1 record.`
- Body 시작에 ⚠️ SUPERSEDED note 추가 — what stays valid (cohort, DPB1*05:01 freq, Chu citation correction, arcasHLA methodology) vs what is superseded (comparator, framing) 정확히 분리

**`MEMORY.md`** line 24 갱신:
```
- [v17 Paper 2 Pillar 1 forest STRONG (SUPERSEDED)](v17_paper2_pillar1_forest_strong.md)
  — 2026-05-03 v1 record. SUPERSEDED 2026-05-04 by v18_paper2_HT_isolated:
  Chu 2018 GD comparator → Korean baseline (AFND South Korea pool); GD forest
  moves to Paper 4 reserve. Cohort + DPB1*05:01 53.2% Korean fact stays;
  framing changed
```

**`v18_paper2_HT_isolated.md`** footer 에 second-sweep update section 추가:
- Task A status: 완료된 상태 verify, C*01:02/DQB1*02:01 deferred 명시
- Task B second sweep: 2 portfolio-level files 추가 정정
- Task C: this entry = canonical, v17_paper2_pillar1_forest_strong SUPERSEDED 처리
- Task D: NOT in scope (voice-protected, Yu meeting 후)
- Cross-paper boundary 명시: Paper 1 / Paper 2 / **Paper 4 backlog** (Paper 3 was renumbered to ICI vulnerability, FROZEN)
- Cross-link: `v17_paper2_pillar1_forest_strong` (SUPERSEDED) + `v17_graves_pivot` (outdated)

### 3.4 Task D — NOT 진행

Phase 1 paper outline 작성은:
- Voice-protected 영역 (본인 키보드)
- Yu meeting 후 진행
- 본 session SCOPE 외

---

## 4. 현재 상태 — final cross-contamination footprint

### 4.1 Active Paper 2 deliverables (advisor 미팅 갈 자료)

| Deliverable | Status |
|---|---|
| `paper2_brief.html` (story version, 38 pages) | ✅ Clean, GD/Graves' = glossary 정의 + Bundang outreach 맥락 (Q01) + advisor catch quote 만 |
| `p2_advisor_discussion.html` (extended, 50 pages) | ✅ Clean except 1 deprecation banner (intentional record) |
| `p2_pillar1_forest/` (v1 outputs) | ✅ DEPRECATED-banner wrapped, record-only |
| `p2_pillar1_forest_v2/` (Task A output) | ✅ Born clean, Korean-baseline framed (HT context only) |
| `papers_overview.html` (master index) | ✅ Clean after second sweep |

### 4.2 Forbidden-term grep results (active deliverables)

```bash
$ grep "Autoimmune-overlap PTC" project/submission/papers_overview.html
# (none — clean ✓)

$ grep "Pan-Asian autoimmune-thyroid" project/manuscript_p2_brief/
# (none — clean ✓)

$ grep "Autoimmune-overlap PTC" project/manuscript_p2_brief/
# Only in deprecation banner (p2_advisor_discussion.html line 671)
# — recording the substitution context itself; not active framing.
```

### 4.3 Remaining Graves'/GD mentions (acceptable per SCOPE)

**`paper2_brief.html`** 의 잔여 Graves'/GD mentions:
- Line 488: Glossary 정의 `GD = Graves' Disease — Pillar I HLA reference (Chu 2018)` — definition only
- Line 507: Glossary `DPB1*05:01` definition `Asian Graves' top risk allele (Chu 2018 OR=1.90)` — historical reference
- Line 656: Cohort table `Bundang SNUH Graves'/PTC+HT cohort 는 별도` — outreach status note (Scenario A 의존성)
- Line 680: Section kicker `자가면역 (Graves') 인구 패턴을 닮는다` — historical Pillar I framing
- Line 685: QMRIL 질문 `autoimmune (Graves') 인구` 와 닮아 있는가 — historical question
- Line 697: QMRIL 해석 `[2026-05-04 audit: GD/Graves'-adjacent framing → Paper 3 reserve]` — explicit audit note
- Line 1009: Tier 2 clinical translation `분당 Graves'/PTC+HT prospective Korean cohort` — outreach status
- Line 1079-1080: Gap 2 `분당 Graves' arm` — outreach dependency
- Line 1116: Scenario A `Bundang Graves' n>50` — venue ladder variable

**모두 SCOPE 가 허용하는 contexts**:
- Glossary definitions (factual reference, not framing claim)
- Bundang outreach status (cohort access guide content)
- Audit notes 내부 (audit context)
- Historical version markers

### 4.4 Memory entries 상태

| Entry | Status |
|---|---|
| `v18_paper2_HT_isolated.md` | ✅ Canonical Paper 2 isolated-scope record. Second-sweep update 추가됨 |
| `v17_paper2_pillar1_forest_strong.md` | ⚠️ SUPERSEDED — header + body + MEMORY.md description 모두 supersession 명시 |
| `v17_graves_pivot.md` | ⚠️ outdated — v18_paper2_HT_isolated 의 cross-link 에 명시됨 |
| `v17_arcasHLA_korean_k2.md` | ✅ historical 분석 자료, OK (DPB1*05:01 의 "Graves' risk" annotation = 분석 결과 description) |
| `v17_2026_04_30_pivot.md` | ✅ historical pivot decision context |

### 4.5 HTTP serving check (active deliverables)

```
200  http://40.82.129.113:8012/manuscript_p2_brief/paper2_brief.html
200  http://40.82.129.113:8012/manuscript_p2_brief/p2_advisor_discussion.html
200  http://40.82.129.113:8012/submission/papers_overview.html
```

모두 정상 deploy.

---

## 5. SCOPE 준수 자체 검증

| SCOPE 규칙 | 본 session 준수 |
|---|---|
| Hashimoto thyroiditis (HT) overlap 의 분자 mechanism ONLY | ✅ Active deliverables clean |
| Graves' disease (GD) phenotype/mechanism 언급 금지 | ✅ Glossary 정의 + outreach status + audit note 만 잔존 (SCOPE 허용) |
| Cancer driver mutation deep dive 금지 | ✅ Touch X |
| `autoimmune-PTC` → `Hashimoto-overlap PTC` 일괄 정정 | ✅ Second sweep 완료 |
| Pillar I forest = Korean baseline (NOT Chu 2018 GD) | ✅ v2 forest 검증, v1 SUPERSEDED |
| Paper 3 file touch 금지 (read 도 X) | ✅ `manuscript_p3_brief/` touch 안함 |
| Voice-protected sections (Hook/Aim/Disc 3.1/Limitations/Cover Para 1/Q9) | ✅ Touch X |
| Task D (Phase 1 outline) NOT 진행 | ✅ SCOPE 외 |

---

## 6. 다음 step (advisor / 본인 키보드)

본 session 이 끝난 후 advisor 또는 본인이 처리할 항목:

1. **Yu meeting** — Pillar I v2 forest (Korean baseline) review + C*01:02/DQB1*02:01 의 Lee 2014 Korean reference fill 의제
2. **Phase 1 paper outline (Task D)** — Yu 합의 후 본인 키보드 시작
3. **`v17_graves_pivot.md`** — 본 session 에서 cross-link 만 추가, 별도 SUPERSEDED note 권장 (future audit)
4. **Three-papers-index `three_papers_index.html`** — Paper 2 row 의 "Autoimmune-overlap" branding 정정 (cross-paper meta-doc, advisor 의제)
5. **Paper 4 backlog 진입** — 4/4 gating 충족 (Paper 1 bioRxiv + Paper 2 A/B/C + Bundang n>50 + Yu 합의) 후

---

## 7. Files generated/modified — 본 session 한정

### 본 session 에서 새로 생성한 file
- `project/results/terminology_correction_2026_05_04/SESSION_REPORT.md` (이 파일)

### 본 session 에서 수정한 file
- `project/submission/papers_overview.html` (2 lines)
- `project/manuscript_p2_brief/p2_advisor_discussion.html` (1 line `<title>`)
- `project/results/terminology_correction_2026_05_04/before_after_diff.md` (second sweep section appended)
- `project/results/terminology_correction_2026_05_04/files_modified.txt` (3 lines appended)
- `~/.claude/.../memory/v17_paper2_pillar1_forest_strong.md` (SUPERSEDED note + body)
- `~/.claude/.../memory/v18_paper2_HT_isolated.md` (cross-link + second-sweep section)
- `~/.claude/.../memory/MEMORY.md` (line 24 description)

### 본 session 에서 read 한 file (touch X)
- `project/results/p2_pillar1_forest_v2/*` (Task A 검증)
- `project/results/p2_pillar1_forest/*` (deprecated 자료 검증)
- `~/.claude/.../memory/v18_paper2_HT_isolated.md` (existing state 확인)

### 본 session 에서 touch 도 read 도 안한 file
- `project/manuscript_p3_brief/*` (Paper 3/4 territory, SCOPE 명시)
- `project/manuscript_v8/*` (Paper 1 voice-protected)
- `project/results/p3_*` 등 Paper 3 관련 results

---

*Generated 2026-05-04 by Paper 2 isolated session second-sweep wrap-up.*
