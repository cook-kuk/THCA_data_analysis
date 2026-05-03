# THCA 3-Paper 마라톤 — 전체 상황 스냅샷 (2026-05-04)

**Author:** Seungho Cook
**작성 시각:** 2026-05-04 (직전 Claude Code 세션 크래시 복구용)
**목적:** 오프라인 상태에서 다음 행동 결정용 단일 문서. 별도 자료 없이 이 파일만으로 판단 가능하도록 구성.
**모드:** 마라톤 (5/4 ~ 6/13). Paper 1 bioRxiv 6/13 ship target.

---

## 0. 한 줄 요약 + 결정해야 하는 질문

**상태:** 직전 세션은 새벽 00:04~00:57 사이 (1) Paper 3 ICI Track A 동결 → (2) Paper 1 마라톤 audit 5종 + 통합 bundle → (3) HIGH/MEDIUM 6개 닫는 Closure 작업까지 완료한 뒤 크래시. **데이터 손실 없음. 모든 산출물 디스크 보존.**

**지금 결정해야 할 질문 3개 (병렬):**

| Q | 질문 | 옵션 | 추천 |
|---|---|---|---|
| **Q1** | Paper 3 ICI 다음 단계 | (A) Track B 시작 / (B) Track-A scaffolding 추가 / (C) 마라톤 복귀 | **(C)** |
| **Q2** | 직전 세션 미커밋 변경 처리 | 커밋 / 폐기 / 보류 | **커밋** (위험 없음) |
| **Q3** | Paper 1 다음 행동 | LOW 9건 cosmetic / Discussion §3.1 본인 키보드 / Methods M3-M9 inline cite / Voice-protected draft | 본인 판단 |

답은 끝의 §10 결정 메뉴에서 한 줄씩 골라서 회신하시면 됩니다.

---

## 1. 직전 세션 timeline (2026-05-04 새벽)

| 시각 | 산출물 | 상태 |
|---|---|---|
| 00:04 | `paper3_ici_dataset_registry.md` | Paper 3 Deliverable 1 |
| 00:05 | `paper3_ici_signature_registry.md` | D2 |
| 00:06 | `paper3_ici_scRNA_reference_atlas_summary.md` | D3 |
| 00:07 | `paper3_ici_HLA_neoantigen_feasibility.md` | D4 |
| 00:08 | `paper3_ici_DIAL_audit_plan.md` | D5 |
| 00:09 | `paper3_ici_figure_plan.md` | D6 |
| 00:10 | `paper3_ici_go_no_go_verdict.md` | D7 (Conditional GO) |
| 00:11 | `paper3_ici_12week_execution_plan.md` | D8 |
| 00:16 | `PAPER3_ICI_TRACK_A_BUNDLE.md` (72 KB) + `paper3_ici_track_a.tar.gz` (51 KB) | Track A 통합 + chmod 444 동결 |
| 00:17–00:18 | todo/ 3종 checklist | editable (644) |
| 00:22 | `PAPER3_ICI_TRACK_A_FREEZE_2026_05_04.md` (root) | 동결 통지서 |
| 00:24 | `PAPER3_ICI_DECISION_BRIEF_2026_05_04.md` (root) | A/B/C 결정 요청 brief |
| 00:38 | `pathology_dm1_phaseA_cpu_verdict_2026_05_04.md` | Paper 2 GPU/CPU 평가 |
| 00:39 | `2026_05_04_paper1_numerical_consistency_audit.md` | Paper 1 T1 |
| 00:41 | `2026_05_04_paper1_bib_citation_audit.md` | T2 |
| 00:42 | `2026_05_04_paper1_xref_audit.md` | T3 |
| 00:43 | `2026_05_04_paper1_methods_prose_gap_report.md` | T4 |
| 00:43 | `2026_05_04_GPU_needs_runpod_assessment.md` | infra |
| 00:44 | `2026_05_04_paper1_supp_tables_completeness.md` | T5 |
| 00:46 | `PAPER1_MARATHON_AUDIT_BUNDLE_2026_05_04.md` (51 KB, 955 lines, **chmod 444**) | T1–T5 통합 + 의사결정 1-pager |
| 00:52 | `2026_05_03_results_R1_R5_prose.md` 수정 | **HM closure 적용** |
| 00:52 | `2026_05_03_manuscript_v8_OUTLINE.md` 수정 | HM closure |
| 00:53 | `2026_05_03_figure_captions_all.md` 수정 | HM closure |
| 00:55 | `2026_05_03_paper_supp_tables_draft.md` 수정 | HM closure |
| 00:56 | `2026_05_04_paper1_HM_closure_report.md` | Closure 리포트 (HIGH 2/2 + MEDIUM 4/6) |
| 00:57 | `2026_05_03_reviewer_QA_consolidated.md` 수정 | HM closure (마지막 변경) |

**해석:** "ici 하던게 꺼짐" 기억은 00:24 시점 ICI Decision Brief 작성 직후를 가리키는 것으로 보이지만, 실제 세션은 그 뒤로도 ~33분 더 진행됨 → Paper 1 마라톤 audit + HM closure까지 완료된 뒤 크래시. **잃어버린 작업 없음.**

---

## 2. 미커밋 인벤토리 (git status)

### 2.1 수정된 파일 (M) — 5개

모두 Paper 1 HM closure가 적용된 manuscript 파일. **HM closure report와 정확히 일치**. 안전하게 커밋 가능.

```
M  project/reports/2026_05_03_figure_captions_all.md       (F5D + F7A 패널)
M  project/reports/2026_05_03_manuscript_v8_OUTLINE.md     (Abstract + R5d 2 lines)
M  project/reports/2026_05_03_paper_supp_tables_draft.md   (S1/S3/S8 source paths)
M  project/reports/2026_05_03_results_R1_R5_prose.md       (R5b/R5d 2 paragraphs)
M  project/reports/2026_05_03_reviewer_QA_consolidated.md  (Q12 2 places)
M  .gitignore                                              (작은 변경, 확인 필요)
```

`git diff --stat HEAD` (closure report 기록): **5 files, +25/-17, net +8 lines**.

### 2.2 신규 파일 (??) — 주요 카테고리별

**Paper 3 ICI (디자인 동결물, root + paper3_ici/):**
- `PAPER3_ICI_TRACK_A_FREEZE_2026_05_04.md`
- `PAPER3_ICI_DECISION_BRIEF_2026_05_04.md`
- `project/reports/paper3_ici/` (전체 디렉터리, chmod 444 디자인 + tar.gz + todo/)

**Paper 1 audit + closure (5/3~5/4 신규):**
- `project/reports/PAPER1_MARATHON_AUDIT_BUNDLE_2026_05_04.md` (chmod 444)
- `project/reports/2026_05_04_paper1_*.md` × 6 (audit T1–T5 + GPU + closure)
- `project/reports/2026_05_04_paper1_HM_closure_report.md`

**Paper 2 / Pathology DM1 (이전 세션 산출물):**
- `PATHOLOGY_DM1_FEASIBILITY_2026_05_03.md`, `PHASE_A_NOGO_REPORT_2026_05_04.md`
- `azure_gpu_phaseA/`, `project/results/03_pathology_poc/`
- `project/reports/pathology_dm1_phaseA_cpu_verdict_2026_05_04.md`

**Paper 1 spatial / external ST 검증 (이전 세션):**
- `EXTERNAL_ST_TRIAGE_ALL_2026_05_03.md`, `GSE250521_*.md`
- `PAPER1_SPATIAL_SUPPLEMENT_FREEZE_2026_05_03.md`
- `project_external_st/`, `project/results/01_spatial_score/`, `02_stage_trend/`, `figures_for_advisor/`

**Supplementary 데이터 (HM closure로 새로 생성된 표):**
- `project/results/p2_pillar1_forest/cohort_assembly.tsv` (H2 build, 1.8 KB)
- `project/results/p3_gse286332/SuppTable_S3_GSE286332_DEGs.xlsx` (M5, 2.2 MB)
- `project/results/d6p7_dm1_subcluster/SuppTable_S8_DM1_subcluster.xlsx` (M6, 3.1 MB)

**기타:**
- `project/src/`, `project/supplementary/`, `project/results/00_qc/`

### 2.3 권장 커밋 그룹화 (3 commit)

1. **Paper 3 ICI Track A FREEZE bundle** — Paper 3 디자인 동결 (모두 chmod 444 안전)
2. **Paper 1 marathon audit + HM closure** — bundle + audit T1–T5 + HM closure + 5 manuscript file edits + 3 supp tables
3. **이전 세션 잔여물 (spatial, pathology, external_st)** — 별도 또는 분할

위험: `.gitignore` 변경 1줄은 diff 확인 후 결정.

---

## 3. 3-Paper 1-screen 현황표

| Paper | 제목 | 상태 | 다음 마일스톤 | 차단? |
|---|---|---|---|---|
| **Paper 1** | DM1 molecular dark matter (8-gene + Pillar 1 forest + GSE286332 + DM1 sub-A/B + TCGA Hashi-like) | manuscript v8 scaffold + audit done + HM closure done | bioRxiv submit **6/13** | LOW 9건 cosmetic + voice-protected sections (본인 키보드) |
| **Paper 2** | H&E → DM1 / pathology projection / TCGA validation (HT-only scope) | Pillar I v2 forest done; Pathology Phase A NO-GO 결정 | Task A/B/C 마무리 | Phase A NO-GO 후 재설계 필요 여부 |
| **Paper 3** | HLA loss · neoantigen · ICI vulnerability in dark thyroid cancer | **Track A FROZEN** (8 deliverables + bundle, chmod 444) | Track B는 (1) Paper 1 bioRxiv + (2) Paper 2 A/B/C 종료 + (3) "Paper 3 Track B 시작" 명시 명령 | **3 조건 모두 미충족 → BLOCKED** |
| **Paper 4** | Korean GD HLA / Pan-Asian (구 Paper 3) | backlog, 4/4 entry condition gated | Paper 1·2·3 정리 후 | Paper 3 freeze 동안 작업 금지 |

---

## 4. Paper 1 상세

### 4.1 디자인 (5/3 base 완료)

manuscript v8 scaffold = OUTLINE / Methods M1–M11 prose / Results R1–R5 prose / Figure & Supp captions / Supp Tables S1–S8 spec / Reviewer Q&A consolidated / Pre-submission checklist / STAR Methods Key Resources / cover letter / data availability / code repo README / GitHub Actions CI / acknowledgments / submission workflow / venue comparison (Cell Rep Med vs JCI Insight). 9 manuscript 핵심 파일.

### 4.2 Cohort sizes (audit 통과 — 모두 ✓ CONSISTENT)

TCGA-THCA 500 (DM1 140 / DM2 360); BRAF-V600E 273 vs WT 182; K2 PRJEB11591 raw 260 / valid HLA 235; Lee2024 GSE213647 raw 632 / valid HLA 630; GSE286332 18 (PTC 9 + PTC+HT 9); **Korean PTC pool 874** (235+630+9 산술 ✓); Chu2018 GD 1,468 / ctrl 1,490 / total 2,958; DM1 sub-A 84 / sub-B 56.

### 4.3 직전 세션 audit 결과 (PAPER1_MARATHON_AUDIT_BUNDLE 0:46 freeze)

5종 audit 통합 → 의사결정 1-pager (issue inventory by severity). HIGH 2건 + MEDIUM 6건 + LOW 9건.

### 4.4 직전 세션 HM closure (00:52~00:57, 5 file edit + 3 supp build)

**닫음:**
- **H1** DM1 sub-B "96%" mut-neg 산술 inconsistency → 후보 (a) 채택, **53/56 = 94.6%** 통일 (5 파일 7 위치)
- **H2** S1 cohort_assembly.tsv "build TBD" → 5×10 컬럼 TSV 빌드
- **M1** sub-A 69% vs 61% denominator → "51/74 mut-tested (69%); 61% of n=84 sub-A" 분리 표기
- **M2** F6/F7 R-PRO 미참조 → "(Figure 6)" / "(Figure 7)" 명시
- **M5** S3 DEG XLSX 미패키지 → 29,672행 XLSX 빌드 (padj 정렬)
- **M6** S8d sub-A/B DEG → S8a/b/c/d 4 시트 XLSX 빌드

**보류 (M3, M4):** Lee2024 GSE213647 bib 저자 placeholder + Lim2025 GSE286332 bib title 검증 → GEO/PubMed lookup 필요, 마라톤 displacement 안전 (이번 주~다음 주 수동).

**LOW 9건 (cosmetic, W6 6/8–6/13 한 번에 처리):**
- L1 ARI 반올림 (0.90 vs 0.903) / L2 Fisher p (6e-10 vs 6.4e-10) / L3 BRAF Cohen d (-0.04 vs -0.044) / L4 ref bib header "~22 papers" → 실제 29 / L5 ★ Unicode in Chu2018JMG / L6 Pan2025NatComm pages 형식 / L7 4 bib (Tuttle2019, Yi2016, Chen2024, Kim2014) 미인용 / L8 "★ exceptional" F2D caption 정서적 / L9 Krishnamoorthy 2025 bib 결손 가능성

### 4.5 Voice-protected sections (본인 키보드 strict)

`v17_sprint_vs_marathon_violation` 메모리 규칙: 마라톤 모드에서 "고고"/"faster"/"다 해줘" 가 voice-protected sprint generate 권한이 NOT 됨. Default = scaffolding/infra. **본인 키보드 의무:**
- Hook ¶1 (introduction 첫 문단)
- Aim ¶4
- Discussion §3.1 (Krishnamoorthy 2025 → Landa 2016 cite save 포함; 메모리 `v17_landa2016_cite_save`에 3-layer 차단 논리 보존)
- Discussion §3.4 Limitations
- Cover letter ¶1
- Reviewer Q9

이 섹션들은 manuscript v8에 아직 비어있거나 scaffold 상태로만 있음 → **본인 작업 슬롯 6/13까지 확보 필요**.

### 4.6 6/13 ship 까지 남은 일 (체크리스트)

- [ ] HM closure 결과 git commit (위험 없음, 5 file edit + 3 신규 supp)
- [ ] M3 + M4 bib 검증 (GEO/PubMed lookup, ~20 분)
- [ ] LOW L1–L9 cosmetic 일괄 (W6, 1~2 시간)
- [ ] Voice-protected 6 sections 본인 키보드 작성 (Hook / Aim / §3.1 / Limitations / Cover ¶1 / Q9)
- [ ] Discussion §3.1 Landa 2016 cite save 적용
- [ ] Pre-submission checklist 최종 sweep
- [ ] bioRxiv 메타데이터 + 첨부 + cover letter 제출 (6/13)

---

## 5. Paper 2 상세

### 5.1 Scope (`v18_paper2_HT_isolated` 메모리, 2026-05-04 Yu professor 결정)

- **Paper 2 = Hashimoto-overlap PTC ONLY** (HT-isolated)
- Pillar I v2 = **Korean PTC vs Korean baseline** (AFND South Korea pool)
- GD / Chu 2018 forest는 **Paper 4로 이동**
- 금지어 리스트 + Task A/B/C anchor 파일 존재

### 5.2 Pillar 1 forest STRONG (`v17_paper2_pillar1_forest_strong` 메모리, 2026-05-03)

Korean PTC pool n=874 vs **Chu X et al. 2018** (NOT Chen) random-effects DerSimonian-Laird forest meta. DPB1*05:01 53.2% Korean > 44.0% Chu GD; 6 alleles direction-consistent. **단, Yu 결정 후 Pillar I v2 = Korean baseline 비교로 변경, Chu 2018 자체는 Paper 4로.**

### 5.3 Pathology DM1 Phase A (00:38 verdict)

`pathology_dm1_phaseA_cpu_verdict_2026_05_04.md` + `PHASE_A_NOGO_REPORT_2026_05_04.md` (root). Phase A는 CPU/GPU 둘 다 평가 → **NO-GO 결정**. 즉 H&E → DM1 projection 1차 디자인은 일단 보류. `azure_gpu_phaseA/` 디렉터리에 작업물 보존됨. `2026_05_04_GPU_needs_runpod_assessment.md`로 GPU 재시작 옵션 평가 가능.

### 5.4 Task A/B/C 상태 (커밋 메시지 기준)

- **Task A:** Pillar I forest v2 (HT-only) — 8e2c827 커밋에서 done
- **Task B:** terminology correction (HT-only scope) — c060d08 커밋
- **Task C:** memory entry — 8e2c827 커밋

→ Task A/B/C 형식상 closed. **다만 Phase A NO-GO 결정 이후 Paper 2 다음 단계(다른 pathology 데이터셋? TCGA validation?)가 미정.**

### 5.5 외부 ST / GSE250521 (이전 세션 산출물, 미커밋)

- `EXTERNAL_ST_TRIAGE_ALL_2026_05_03.md` — 외부 spatial transcriptomics 데이터셋 triage
- `GSE250521_ALL_2026_05_03.md` + 4종 보고서 (validation / inclusion verdict / meaningful or not / FOR_WEB_CLAUDE)
- `paper1_spatial_supplement_freeze_2026_05_03.md`
- `PAPER1_SPATIAL_SUPPLEMENT_FREEZE_2026_05_03.md` (root)
- `project/results/01_spatial_score/`, `02_stage_trend/`, `figures_for_advisor/`

→ **Paper 1에 GSE250521이 들어가는지 결정** (verdict 파일 존재) + **외부 ST 결과를 Paper 1 supplement로 freeze** 했는지 재확인 필요.

---

## 6. Paper 3 ICI 상세 (FROZEN)

### 6.1 정체성

- **Working title:** "HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability in molecularly dark thyroid cancer"
- **Allowed claim:** "ICI vulnerability" / "ICI-readiness" / "immunogenomic prioritization"
- **FORBIDDEN claim:** "ICI response predictor" (thyroid ICI raw RNA-seq 부재 → kill switch K1)

### 6.2 Track A 8 deliverables (00:04~00:11, 모두 chmod 444)

1. **Dataset Registry** — bulk ≤9 / scRNA ≤6 / spatial ≤4 / pan-cancer ICI ≤7. 모두 `to_verify` 상태. cross-paper boundary 주석.
2. **Signature Registry** — IFNγ / TIS / Cytolytic / MHC-I / MHC-II / TLS / CXCL13 / Myeloid / Treg / TIDE-like / IMPRES + tumor-intrinsic (BRS/ERK/TDS/lineage) + sc-derived placeholders. 통합 score 가중치는 Module E로 deferred.
3. **scRNA Atlas Plan** — scVI/scANVI/Harmony, Level-1/2 taxonomy, sc-derived signature lock for bulk projection. PDTC/ATC <5K caveat.
4. **HLA & Neoantigen Feasibility** — OptiType / Polysolver / arcasHLA / xHLA + LOHHLA + NetMHCpan / NetMHCIIpan. **cookHLA은 Paper 4 reserved.** dbGaP unauthorized 시 degraded scope.
5. **DIAL Audit Plan** — sign-consistency × effect-size × tissue-transfer → PASS/FLIP/COLLAPSE/AMBIGUOUS verdict per signature. 통합 readiness score는 PASS 시그니처만.
6. **Figure Plan** — F1 cohort/dark matter / F2 ecotype NMF / F3 sc atlas / F4 HLA+neoantigen / F5 DIAL / F6 integrated score; 14 supp.
7. **Go/No-Go Verdict** — Conditional GO. G1–G6 hard gates + K1–K6 kill switches.
8. **12-Week Execution Plan** — Wk 1 verify+ETL → Wk 10 integrated score → Wk 11–12 manuscript draft.

### 6.3 동결 패키지

- `project/reports/paper3_ici/PAPER3_ICI_TRACK_A_BUNDLE.md` — 71.7 KB, 1236 lines, 8 deliverables 통합 (chmod 444)
- `project/reports/paper3_ici/paper3_ici_track_a.tar.gz` — 51 KB, sha256 `290349f6ebed7de64c448afbebf633baa05f1845e0685f7e7a2df9fb3d8936aa`
- `project/reports/PAPER3_ICI_TRACK_A_FREEZE_2026_05_04.md` — root freeze notice
- `project/reports/PAPER3_ICI_DECISION_BRIEF_2026_05_04.md` — A/B/C 결정 brief

### 6.4 Editable todo 3종 (chmod 644, freeze 외부)

- `paper3_ici/todo/paper3_ici_track_b_kickoff_checklist.md` — Track B kickoff verification (G1–G6 gate)
- `paper3_ici/todo/paper3_ici_controlled_access_checklist.md` — dbGaP/EGA 신청 (≤30분/세션, 마라톤-displacement-safe)
- `paper3_ici/todo/paper3_ici_dataset_accession_verification_checklist.md` — Wk 1 default deferred

### 6.5 Track B unlock 조건 (3개 모두 필요)

1. Paper 1 bioRxiv 제출 완료
2. Paper 2 Task A/B/C 종료 (현재 형식상 closed지만 Phase A NO-GO 후 재정의 필요할 수 있음)
3. **명시적 명령:** "Paper 3 Track B 시작"
   - "고고" / "다 해줘" / "faster" / "추가 분석 더 해줘" → **NOT count**

### 6.6 Decision Brief A/B/C (재게시)

| 옵션 | 액션 | 비용 | gate | verdict |
|---|---|---|---|---|
| **(A)** | Track B 즉시 시작 (12주 plan Wk 1) | Paper 1 bioRxiv 6/13 슬립 + Paper 2 정지 + 컴퓨트 1500–2500 CPU-hr + 2 GPU-day | **FAILS** | NO. Yu 명령 같은 override 사유 필요 |
| **(B)** | Track-A 호환 scaffolding 1개 추가 (① reviewer Q anticipation / ② methods M1–Mxx scaffold / ③ SAP / ④ pre-reg template / ⑤ cover letter / ⑥ code-repo README / ⑦ ethics 보일러) | Paper 1/2 시간 displacement 2–6 hr/항목, 컴퓨트 0 | Track B 여전히 BLOCKED, 마라톤 displaced 하지만 strict 위반 아님 | 항목 명시 + ≤4시간 timebox 시 acceptable |
| **(C)** | 마라톤 복귀 (default) | 0 | full discipline | **RECOMMENDED** |

---

## 7. Paper 4 backlog

`v19_paper4_GD_backlog` 메모리 (2026-05-04 renumber). Korean GD HLA / Pan-Asian. 4/4 entry condition gated, scope rules는 prior `v19_paper3_GD_gating` 동일. Paper 3 freeze 동안 작업 금지.

Pillar 1 forest의 Chu 2018 비교 분석은 **Paper 4로 옮겨감** (Paper 2가 HT-only로 좁아진 후속).

---

## 8. 마라톤 모드 규칙 (위반 회피)

`v17_marathon_mode_post_pillar1` + `v17_sprint_vs_marathon_violation` 메모리:

- **Window:** 5/4 ~ 6/13 (6주). 새 분석은 **paper-blocking only**.
- **Default 권한:** scaffolding / infra만.
- **Voice-protected 6 sections** (위 §4.5)는 **본인 키보드 strict**. "고고"/"faster"/"다 해줘"로 Claude가 generate 권한 NOT.
- **Track B unlock 명령 ≠ "고고"** — Paper 3에 한해 **literal "Paper 3 Track B 시작"** 만 인정.
- **컴퓨트 spend** = paper-blocking 만. 일반 audit/scaffolding 은 read-only OK.

---

## 9. 외부 의존성 / 미해결 결정

| 항목 | 상태 | 결정 필요? |
|---|---|---|
| Bundang SNUH outreach (이메일 5/3 final draft) | 발송 완료 여부 미확인. `v17_K2_vs_bundang_distinction` 메모리: Bundang은 outreach-stage, 데이터 없음 | 발송 상태 확인 |
| dbGaP / EGA 신청 (Paper 3 controlled access) | 미신청. todo/ checklist 생성됨 | Track B 시작 결정 후 |
| AFND South Korea pool (Paper 2 Pillar I v2 baseline) | 데이터 풀 정의 상태 미확인 | Task A 재확인 |
| GSE250521 Paper 1 inclusion | verdict 파일 있음 (`gse250521_paper1_inclusion_verdict_2026_05_03.md`), 결과 미확인 | 본문 read 필요 |
| Lee2024 / Lim2025 bib 저자·제목 (M3, M4) | GEO/PubMed lookup 미실행 | 20분 작업 |
| Voice-protected 6 sections | scaffold만, 본인 작성 미시작 | 본인 일정 |

---

## 10. 결정 메뉴 (이걸 답변으로 회신)

### Q1 — Paper 3 ICI

- `(A)` Track B 시작 (override 사유 명시)
- `(B) ① / ② / ③ / ④ / ⑤ / ⑥ / ⑦` 중 하나, ≤4시간 timebox
- `(C)` 마라톤 복귀 ✅ default

### Q2 — 미커밋 변경

- `commit-3groups` (P3 freeze / P1 audit+closure / 잔여물 3개로 분할 커밋)
- `commit-2groups` (P3 freeze 1개 + 나머지 1개)
- `commit-all-one` (전부 1커밋)
- `commit-skip` (보류 — 추천 안 함, 작업 손실 위험은 없지만 inventory 흐림)

### Q3 — Paper 1 다음 행동 (마라톤 복귀 시)

- `voice-hook` Hook ¶1 본인 키보드 시작
- `voice-aim` Aim ¶4
- `voice-disc31` Discussion §3.1 (Landa 2016 cite save 적용)
- `voice-lim` §3.4 Limitations
- `voice-cover` Cover letter ¶1
- `voice-q9` Reviewer Q9
- `bib-m3m4` Lee2024 / Lim2025 bib lookup (~20분, 마라톤-safe scaffolding)
- `low-batch` LOW L1–L9 cosmetic 일괄 (1~2시간; 정상은 W6에 두지만 당겨 가능)
- `audit-recheck` HM closure 산출물 재검증 (XLSX 시트 / TSV 행수 확인)
- `gse250521-decide` Paper 1 GSE250521 inclusion 결정 read-through
- `paper2-redesign` Phase A NO-GO 후 Paper 2 다음 디자인 (마라톤-blocking 가능성 있음 → 신중)

### Q4 (선택) — Paper 2

- `paper2-task-status-recheck` Task A/B/C 결과물 보존 상태 점검
- `paper2-phaseB-design` Phase B 데이터셋 후보 brainstorm (scaffolding only)
- `paper2-pause` Paper 1 끝낼 때까지 보류

---

## 11. 핵심 파일 인덱스

### 마라톤 audit
- `project/reports/PAPER1_MARATHON_AUDIT_BUNDLE_2026_05_04.md` (chmod 444)
- `project/reports/2026_05_04_paper1_HM_closure_report.md`

### Paper 1 manuscript (HM closure 적용된 5 file)
- `project/reports/2026_05_03_manuscript_v8_OUTLINE.md`
- `project/reports/2026_05_03_results_R1_R5_prose.md`
- `project/reports/2026_05_03_methods_M1_M11_prose.md` (수정 없음)
- `project/reports/2026_05_03_figure_captions_all.md`
- `project/reports/2026_05_03_paper_supp_tables_draft.md`
- `project/reports/2026_05_03_reviewer_QA_consolidated.md`
- `project/reports/2026_05_03_PRE_SUBMISSION_CHECKLIST.md`
- `project/reports/2026_05_03_STAR_Methods_KeyResources.md`
- `project/reports/2026_05_03_cover_letter_assembly.md`

### Paper 1 supp data (HM closure 신규)
- `project/results/p2_pillar1_forest/cohort_assembly.tsv`
- `project/results/p3_gse286332/SuppTable_S3_GSE286332_DEGs.xlsx`
- `project/results/d6p7_dm1_subcluster/SuppTable_S8_DM1_subcluster.xlsx`

### Paper 3 ICI freeze
- `project/reports/PAPER3_ICI_TRACK_A_FREEZE_2026_05_04.md`
- `project/reports/PAPER3_ICI_DECISION_BRIEF_2026_05_04.md`
- `project/reports/paper3_ici/PAPER3_ICI_TRACK_A_BUNDLE.md` (chmod 444, 1236 lines)
- `project/reports/paper3_ici/paper3_ici_track_a.tar.gz` (51 KB)
- `project/reports/paper3_ici/paper3_ici_*.md` × 8 (chmod 444)
- `project/reports/paper3_ici/todo/*.md` × 3 (chmod 644)

### Paper 2 / pathology
- `PHASE_A_NOGO_REPORT_2026_05_04.md` (root)
- `PATHOLOGY_DM1_FEASIBILITY_2026_05_03.md` (root)
- `project/reports/pathology_dm1_phaseA_cpu_verdict_2026_05_04.md`
- `project/reports/2026_05_04_GPU_needs_runpod_assessment.md`
- `azure_gpu_phaseA/`

### Paper 1 spatial / external ST
- `EXTERNAL_ST_TRIAGE_ALL_2026_05_03.md`, `GSE250521_ALL_2026_05_03.md`
- `PAPER1_SPATIAL_SUPPLEMENT_FREEZE_2026_05_03.md`
- `project/reports/gse250521_*.md` × 4
- `project/reports/paper1_spatial_supplement_freeze_2026_05_03.md`
- `project_external_st/`, `project/results/01_spatial_score/`, `02_stage_trend/`, `figures_for_advisor/`

### 메모리 인덱스 (`~/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/`)
- `paper_numbering_2026_05_04.md` — canonical 4-paper numbering
- `v19_paper3_ici_track_a.md` — Paper 3 Track A FROZEN
- `v18_paper2_HT_isolated.md` — Paper 2 = HT-only scope
- `v19_paper4_GD_backlog.md` — Paper 4 backlog 게이팅
- `v17_marathon_mode_post_pillar1.md` — 마라톤 정의
- `v17_sprint_vs_marathon_violation.md` — 본인 키보드 strict
- `v17_landa2016_cite_save.md` — Discussion §3.1 cite 정정
- `v17_paper2_pillar1_forest_strong.md` — Korean PTC pool 874 vs Chu 2018
- `v17_K2_vs_bundang_distinction.md` — K2 vs Bundang
- `v17_D6P7_dm1_subB_NBNR.md` — DM1 sub-B 정체성

---

**다음 행동:** 위 §10에서 Q1 + Q2 + Q3 (그리고 Q4 선택) 한 줄씩 회신하시면 즉시 진행. 회신 없이 **`continue`** 한 단어만 보내면 default 조합 = `(C)` + `commit-3groups` + `voice-hook` 으로 시작.
