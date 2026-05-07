<!-- VS Code: Cmd/Ctrl+K V → side-by-side preview -->
---
title: "Paper 2 — Yu professor meeting review packet"
date: 2026-05-04
session: "Pre-meeting packet, no new analysis executed"
purpose: "한 장으로 advisor 가 결정해야 할 것만 정리. Bullet/table only."
---

# Paper 2 — Yu Meeting Review Packet (2026-05-04)

## 1. Paper 2 one-line scope

> **Paper 2 = Hashimoto-overlap PTC 의 분자 mechanism — antigen-driven HLA-II–매개 dedifferentiation. ONLY HT, not GD.**

- HT = destructive immune infiltration → dedifferentiation
- GD (TSAb-stimulated hyperthyroidism) = Paper 4 backlog
- HLA shared genetic background acknowledged Discussion 한 줄, disease equivalence 주장 X

---

## 2. What changed from v1 to v2

| Layer | v1 (DEPRECATED) | v2 (★ ACTIVE 2026-05-04) |
|---|---|---|
| Comparator | Chu 2018 Han Chinese GD (n=1,468) + ctrl (n=1,490) | AFND South Korea pool (n=4,613 / 680 / 201 per allele) |
| Population frame | Korean PTC vs Han Chinese (cross-population) | Korean PTC vs Korean baseline (same population) |
| Disease scope | "Pan-Asian autoimmune-thyroid continuum" (HT+GD conflated) | "Hashimoto-thyroid overlap (HT context)" — GD = Paper 4 |
| Forest claim | DPB1*05:01 OR 2.50 vs Chu ctrl, pooled OR 2.16 [1.65, 2.83] | DPB1*05:01 OR 1.96 [1.60, 2.41] vs Korean baseline |
| Direction count | 6/6 alleles direction-consistent (Chu 2018 GD priors) | 3/4 risk direction-consistent (DRB1*07:01 caveat) |
| 2 alleles status | C*01:02, DQB1*02:01 included via Chu 2018 | C*01:02, DQB1*02:01 unavailable in AFND South Korea — flagged |
| Discussion GD mention | Pan-Asian framing as main claim | One Discussion line cross-ref to Paper 4 |
| Cell Press positioning | Pan-Asian autoimmune-thyroid susceptibility | Pan-Asian thyroid HLA susceptibility (HT context) |

---

## 3. Pillar I v2 key numbers (Korean PTC vs Korean baseline)

| Allele | Direction | PTC freq (n=874) | Baseline freq | Baseline n / pops | OR (95% CI) | p |
|---|---|---|---|---|---|---|
| **DPB1*05:01** ★ | risk | 53.2% (465/874) | 36.7% | 680 / 3 AFND-SK pops | **1.96 [1.60, 2.41]** | 1.1e−10 |
| A*02:07 | risk | 8.1% (71/874) | 3.4% | 4,613 / 2 AFND-SK pops | 2.54 [1.90, 3.40] | 3.1e−10 |
| B*46:01 | risk | 10.3% (90/874) | 4.7% | 4,613 / 2 AFND-SK pops | 2.33 [1.80, 3.01] | 1.2e−10 |
| DRB1*07:01 | ⚠️ discord | 11.4% (100/874) | 5.0% | 201 / 1 China-Harbin Korean (proxy) | 2.46 [1.26, 4.79] | 0.0084 |
| C*01:02 | n/a | 24.1% (211/874) | — | AFND-SK unavailable | — | — |
| DQB1*02:01 | n/a | 0% (0/874) | — | AFND-SK unavailable | — | — |

**Korean PTC pool composition** (n=874 unchanged from v1):
- K2 / PRJEB11591 (Yoo 2016 SNU-GMI) — n=235 typeable
- Lee 2024 / GSE213647 — n=630
- GSE286332-PTC arm (Lim 2025 Dongguk) — n=9

---

## 4. Why Chu 2018 GD comparator was deprecated

- **Yu advisor catch (2026-05-04)**:
  > "하시모토병은 ... Dedifferentiation 의 파괴적 성향이 강한 반면,
  >  그레이브스 병은 ... 갑상선 기능을 비정상적으로 항진시키는 차이가 있습니다."
- HT and GD are **mechanistically distinct** despite shared HLA background
- v1 "Pan-Asian autoimmune-thyroid susceptibility allele continuum" = HT + GD conflated
- Korean PTC vs Han Chinese GD = **cross-population stratification confound** (different populations)
- v2 Korean baseline = same population, cleaner comparator
- GD-related forest = Paper 4 backlog (gated 4/4 — Paper 1 bioRxiv + Paper 2 A/B/C + Bundang n>50 + Yu approval)

---

## 5. What remains uncertain

### 5.1 C*01:02 / DQB1*02:01 missing AFND entries
- AFND has no South Korea population entry for these 2 alleles
- v1 Chu 2018 reference covered them; v2 currently leaves them blank
- Lee 2014 *Tissue Antigens* Korean reference 가 표준 후보 (web fetch 필요, 본 session 금지)
- 2 alleles 부재 → Pillar I 4 alleles 만 main forest, 6 alleles 전체 framing 약화

### 5.2 DRB1*07:01 direction discord
- v1 (Chu 2018 GD context): protective (ctrl 15.3% → GD 7.1%)
- v2 (Korean baseline): **risk direction** — PTC 11.4% > Harbin Korean baseline 5.0%, OR 2.46
- Hypothesis 1: Harbin Korean baseline n=201 (China-Korean diaspora) 가 true Korean DRB1*07:01 freq under-representation
- Hypothesis 2: Chu 2018 GD-context "protective" 효과가 PTC 로 transfer 안 됨
- Lee 2014 Korean reference 로 hypothesis 1 직접 검정 가능

### 5.3 paper2_brief.html Figure 1 still v1 data
- Plotly JS data = Chu 2018 forest (Korean PTC vs Chu ctrl vs Chu GD)
- v2 audit banner 가 redirect 명시했지만 실제 figure 는 v1 visualization
- v2 figure 는 `p2_pillar1_forest_v2/forest_paper2_HT_only.{pdf,png}` 별도 파일

### 5.4 paper2_brief.html Pillar I QMRIL body partial v2
- Section title + kicker + audit banner + 질문 line = v2 framed ✅
- 방법/결과/해석/한계 (4 QMRIL rows) = v1 narrative 보존
- Brief 가 "audit banner + body 보존" 부분-migration 상태

---

## 6. Three advisor questions

| Q | Question | 본인 default position |
|---|---|---|
| **Q1** | Korean baseline AFND pool (sample-weighted; Harbin Korean fallback for some alleles) 이 Pillar I v2 reference 로 acceptable 한가? | Acceptable, single-population stratification advantage. Harbin Korean fallback 명시 + Lee 2014 future-fill 권장 |
| **Q2** | C*01:02 / DQB1*02:01 의 Lee 2014 Korean reference fill 진행할 것인가? | Yes 권장. 2 alleles 추가 시 Pillar I 6-allele 전체 forest 가능, framing 강화 |
| **Q3** | paper2_brief.html Figure 1 (JS data) 와 Pillar I QMRIL body 를 advisor approval 후 fully v2 re-render 할 것인가? | Yes 권장. 현재 audit banner redirect 만으로는 reviewer 혼동 위험 |

---

## 7. Decision table (advisor 선택지)

| Decision | 의미 | 결과 / next action |
|---|---|---|
| **A. Approve v2 as active** | Korean baseline framing 그대로 ship-ready 처리 | paper2_brief.html Figure 1 + QMRIL fully v2 re-render → PDF re-gen |
| **B. Request Lee 2014 fill (C*01:02 + DQB1*02:01)** | 2 alleles 추가 web fetch + analysis | 1-2 hr 작업; v2 forest 4 → 6 alleles 확장 + DRB1*07:01 caveat 검정 |
| **C. Defer Pillar I to next sprint** | 현 audit banner 상태 유지 | Pillar II–V 와 Mediation 만으로 advisor 미팅 진행, Pillar I 후속 |
| **D. Move more HLA material to Paper 4** | DPB1*05:01 / B*46:01 같은 strong Asian risk allele 도 Paper 4 reserve | Paper 2 Pillar I 더 축소; HT-specific allele 만 남김 |
| **E. A + B 결합** | Approve v2 + request Lee 2014 fill | Default 권장 path |

**본인 default 추천 = E** (approve v2 + Lee 2014 fill).

---

## 8. After-meeting action plan (조건부)

### If Decision A or E approved
1. paper2_brief.html Figure 1 JS data v1 → v2 swap (Plotly bar 4 alleles)
2. paper2_brief.html Pillar I QMRIL body 전체 v2 re-write (4 rows)
3. paper2_brief.pdf 재생성
4. p2_advisor_discussion.html 도 동일 migration (extended brief)
5. CHANGELOG + SITUATION_OVERVIEW 갱신
6. ETA: 2-3 hr

### If Decision B (Lee 2014 fill) added
1. Web fetch: Lee 2014 *Tissue Antigens* Korean HLA reference paper, table 추출
2. C*01:02 + DQB1*02:01 Korean baseline freq 채움
3. DRB1*07:01 hypothesis 1 검정 (Lee 2014 freq 가 Harbin Korean 5.0% 와 일치하는가)
4. v2 forest 6-allele 확장 + figure regen
5. PILLAR1_FOREST_V2_SUMMARY.md update
6. ETA: 2-4 hr

### If Decision C (defer)
1. paper2_brief.html 현 audit banner 유지
2. advisor 미팅에서는 Pillar II–V + Mediation 위주
3. Pillar I 후속 sprint 시점 결정

### If Decision D (move to Paper 4)
1. Pillar I scope 재정의 — HT-specific allele 만 남김
2. SCOPE 문서 (`v18_paper2_HT_isolated.md`) 갱신
3. Paper 4 backlog entry 에 추가 HLA material 재분류

---

## 9. Files to bring to meeting

- **★ 본 packet**: `project/manuscript_p2_brief/YU_MEETING_REVIEW_PACKET_2026_05_04.md`
- **Situation overview**: `project/manuscript_p2_brief/SITUATION_OVERVIEW.md`
- **Story brief PDF** (38 pages, Pillar I audit banner 포함): `project/manuscript_p2_brief/paper2_brief.pdf`
- **v2 forest PDF**: `project/results/p2_pillar1_forest_v2/forest_paper2_HT_only.pdf`
- **v2 forest summary**: `project/results/p2_pillar1_forest_v2/PILLAR1_FOREST_V2_SUMMARY.md`
- **v2 Discussion paragraph**: `project/results/p2_pillar1_forest_v2/discussion_paragraph_v2.md`

---

## 10. Cross-paper boundaries (advisor 미팅 재확인)

- **Paper 1** (Cancer / Dark Matter v8) — manuscript_v8/, 8-gene + DM1 fusion + epigenetic, Cell Rep Med target
- **Paper 2** (★ this packet) — Hashimoto-overlap PTC, Korean PTC HLA + GSE286332 mechanism + TCGA/Korean generalization + BCR + DM1 sub-B
- **Paper 3** (FROZEN) — ICI vulnerability dark thyroid cancer (`v19_paper3_ici_track_a` bundle, chmod 444)
- **Paper 4 backlog** (gated 4/4) — Korean GD HLA / Pan-Asian autoimmune-thyroid (`v19_paper4_GD_backlog`), Chu 2018 forest 자료 reserve

---

Yu meeting packet ready. No analysis executed.
