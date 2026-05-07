# THCA 5-Paper 상황 스냅샷 (공유용)

**Author:** Seungho Cook
**작성:** 2026-05-04
**모드:** 마라톤 (5/4 → 6/13). Paper 1 bioRxiv 6/13 ship target.
**목적:** 5개 paper trajectory를 한 문서로. 한 개 실패(image-DM1) 포함 전체 상황 공유.

---

## 1. 한 줄 요약 표

| # | Paper | 한 줄 정체 | 상태 | 다음 마일스톤 |
|---|---|---|---|---|
| **1** | **DM1 molecular dark matter** | 8-gene driver-excluded sub-stratifier (RAI/DM1) | 🟢 **Active marathon** — manuscript v8 + audit + HM closure done | **bioRxiv 6/13** |
| **2A** | **H&E → DM1 image projection** *(원래 Paper 2 디자인)* | spot-aligned H&E tile로 depth-residualized DM1 score 예측 | 🔴 **실패 / 폐기 (NO-GO 2026-05-04)** | **재시도 금지** (5 entry conditions 충족 시까지) |
| **2B** | **Hashimoto-overlap PTC mechanism** *(Paper 2 reframe, 2026-05-04 Yu)* | HT-overlap PTC molecular-only (Pillar I forest + GSE286332 + TCGA HM-like + BCR/TLS + DM1 sub-B = NBNR) | 🟢 **Active marathon** (Task A/B/C 형식 closed; Pillar II/III 정의 TBD) | Pillar II/III 정의 → manuscript |
| **3** | **ICI vulnerability dark thyroid cancer** | HLA loss · neoantigen · ICI-readiness in dark TC | 🟡 **Track A FROZEN (chmod 444)** | Track B = (Paper 1 bioRxiv + Paper 2 A/B/C + 명시 명령) 모두 충족 시 |
| **4** | **Korean GD HLA / Pan-Asian** *(backlog)* | Korean Graves' DPB1*05:01 + pan-Asian forest | ⚪ **Backlog (4/4 entry gated)** | Paper 1·2·3 정리 후 |

> 5개 = "active 3개 (1, 2B, 3 frozen) + backlog 1개 (4) + 실패 1개 (2A)". 한 개 실패 = **2A image-DM1**.

---

## 2. Paper 별 상세

### Paper 1 — DM1 molecular dark matter 🟢

- **Status:** manuscript v8 scaffold + 5종 audit + HIGH/MEDIUM closure 완료. LOW 9건 cosmetic + voice-protected 6 sections 남음.
- **Cohort:** TCGA-THCA 500 (DM1 140 / DM2 360); BRAF-V600E 273 vs WT 182; Korean PTC pool 874; Chu 2018 GD 1,468 / ctrl 1,490; DM1 sub-A 84 / sub-B 56.
- **Key result:** 8-gene driver-excluded RAI score → DM1 sub-B (n=56, 96% mut-neg) = TCGA NBNR equivalent + 4× HM-like rate. TERT⁺ Cox HR=4.33 (p=4.9e-6).
- **Voice-protected (본인 키보드 strict):** Hook ¶1 / Aim ¶4 / Discussion §3.1 (Landa 2016 cite save) / §3.4 Limitations / Cover ¶1 / Reviewer Q9.
- **남은 일:** bib M3/M4 GEO/PubMed lookup (~20분) → LOW L1–L9 일괄 (1–2시간) → voice-protected 6 sections 본인 작성 → bioRxiv 제출 (6/13).
- **Venue:** Sci Rep base / Cell Rep Med / JCI Insight reach.
- **앵커 파일:** `project/reports/PAPER1_MARATHON_AUDIT_BUNDLE_2026_05_04.md` (chmod 444), `project/reports/2026_05_04_paper1_HM_closure_report.md`, `project/reports/2026_05_03_manuscript_v8_OUTLINE.md`.

---

### Paper 2A — H&E → DM1 image projection 🔴 **실패**

**원래 Paper 2의 1차 디자인이었음. 2026-05-04 closure battery로 NO-GO 확정.**

- **Hypothesis:** Paper 1의 depth-residualized 8-gene `DM1_like_score_resid`가 spot-aligned H&E tile로부터 예측 가능.
- **Result:** **모든 model class / target preprocessing / tile size / aggregation에서 reject.**
  - DM1_resid LOSO 224 px ResNet50 ImageNet → Spearman r = **0.022**, AUROC = **0.511** (chance level).
  - 30 random 8-gene panels: real 0.022는 random max 0.031보다 낮음.
  - Raw 신호 (RGB 0.224)는 ResNet50 (0.061)을 이김 → **stain ↔ sequencing-depth artifact**, 생물학적 신호 아님.
  - Tile size 448 px → 0.017 (224보다 worse), 672 cancelled.
  - Foundation model (UNI/CONCH/Virchow2) 최대 +0.15 gain 가정해도 BORDERLINE 0.20 floor 미달.
- **Dominant failure cause:** spot-level ST molecular label과 tile-level H&E의 **fundamental mismatch** (stage ordinal r = −0.05 포함, 모든 continuous molecular target r ≤ 0.08).
- **Dropped:** image-DM1 Paper 2 active pillar / TCGA WSI Phase C ~500 GB / multi-cohort LODO Phase B / RunPod·Azure GPU spend / HF UNI/CONCH/Virchow2 access / Paper 2 H&E IP claim.
- **Re-entry conditions (5/5 모두 필요):** (1) full-res scanner WSI + paired bulk RNA-seq, (2) external molecular label, (3) post-marathon (≥ Paper 1 bioRxiv 후), (4) foundation model access + non-thyroid 벤치 검증, (5) pre-registered 가설.
- **앵커 파일:** `project/reports/2026_05_04_image_dm1_final_nogo_decision.md`, `project/reports/pathology_dm1_closure_battery_2026_05_04.md`, `CLOSURE_BATTERY_2026_05_04.md`, `PHASE_A_NOGO_REPORT_2026_05_04.md`.
- **Archived data (보존, 미삭제):** `project/results/03_pathology_poc/closure_battery_metrics.tsv` 외 8 파일.

---

### Paper 2B — Hashimoto-overlap PTC mechanism 🟢

**2A 실패 후 Yu professor 결정 (2026-05-04)으로 Paper 2 = HT-overlap PTC molecular-only로 좁힘.**

- **Pillar I v2:** Korean PTC pool n=874 vs **AFND South Korea pool** baseline. (이전 Chu 2018 GD 비교는 → **Paper 4로 이동**.)
- **GSE286332 STRONG GO:** Korean PTC vs PTC+HT n=18, 10,380 DEGs, 8-gene d=−1.6, HLA-II d=+3.65, IFN-γ FDR=2e-4.
- **TCGA Hashimoto-like generalization (D4-P2):** GSE286332 PTC+HT signature → TCGA에서 DM2 enrichment OR 최대 5×, p=6e-10. **18/18 DM2 paradox 해소.**
- **BCR clonal + TLS (D5-P6):** PTC+HT TLS d=+1.96, IGHV clonality d>0.5, AICDA up. Antigen-driven B cell response 확정.
- **DM1 sub-B = NBNR (D6-P7):** DM1 sub-B (n=56, 96% mut-neg) = TCGA equivalent of K2 NBNR + 4× HM-like.
- **Task 상태:** Task A (Pillar I v2 forest) + Task B (terminology) + Task C (memory) 형식상 closed. **Pillar II/III 정의 미정** (Yu 다음 미팅 input 대기).
- **금지어:** GD / TSAb / exophthalmos / image-DM1 / pathology projection 표현은 manuscript에서 제거 또는 reframe (audit `2026_05_04_paper2_post_image_dm1_nogo_status.md` §4).
- **Voice-protected:** Hook / Aim / Discussion 3.1 / Limitations / Cover ¶1 / Q9 (본인 키보드 strict).
- **Venue 후보:** JCI Insight / Cell Rep Med (HT-overlap mechanism + TLS + IGHV + TCGA generalization 패키지).
- **앵커 파일:** `project/reports/2026_05_03_paper2_pillar1_for_web_claude.md`, GSE286332/TCGA-Hashimoto/BCR 결과는 `v17_GSE286332_strong_go`, `v17_D4P2_tcga_hashimoto_generalization`, `v17_D5P6_BCR_clonal_TLS`, `v17_D6P7_dm1_subB_NBNR` 메모리 참조.

---

### Paper 3 — ICI vulnerability dark thyroid cancer 🟡 **FROZEN**

- **Working title:** "HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability in molecularly dark thyroid cancer".
- **Allowed claim:** "ICI vulnerability" / "ICI-readiness" / "immunogenomic prioritization".
- **FORBIDDEN claim:** "ICI response predictor" — thyroid ICI raw RNA-seq 부재 (kill switch K1).
- **Track A 8 deliverables (chmod 444 freeze):**
  1. Dataset Registry (bulk ≤9 / scRNA ≤6 / spatial ≤4 / pan-cancer ICI ≤7)
  2. Signature Registry (IFNγ / TIS / Cytolytic / MHC-I/II / TLS / CXCL13 / Myeloid / Treg / TIDE-like / IMPRES + tumor-intrinsic + sc-derived)
  3. scRNA Atlas Plan (scVI/scANVI/Harmony, Level-1/2 taxonomy)
  4. HLA & Neoantigen Feasibility (OptiType / Polysolver / arcasHLA / xHLA + LOHHLA + NetMHCpan/NetMHCIIpan; **cookHLA은 Paper 4 reserved**)
  5. DIAL Audit Plan (PASS/FLIP/COLLAPSE/AMBIGUOUS verdict)
  6. Figure Plan (F1–F6 + 14 supp)
  7. Go/No-Go Verdict (Conditional GO; G1–G6 gates + K1–K6 kill switches)
  8. 12-Week Execution Plan
- **Track B unlock 조건 (3개 모두 필요):**
  1. Paper 1 bioRxiv 제출 완료
  2. Paper 2 Task A/B/C 종료 (Phase A NO-GO 후 재정의 가능성 있음)
  3. **명시적 명령:** "Paper 3 Track B 시작" — "고고"/"다 해줘"/"faster" → **NOT count**
- **앵커 파일:** `project/reports/paper3_ici/PAPER3_ICI_TRACK_A_BUNDLE.md` (chmod 444, 1236 lines), `project/reports/paper3_ici/paper3_ici_track_a.tar.gz` (51 KB, sha256 `290349f6…d8936aa`), `project/reports/PAPER3_ICI_TRACK_A_FREEZE_2026_05_04.md`, `project/reports/PAPER3_ICI_DECISION_BRIEF_2026_05_04.md`.
- **Editable todo (chmod 644):** Track B kickoff checklist / Controlled access (dbGaP/EGA) / Dataset accession verification.

---

### Paper 4 — Korean GD HLA / Pan-Asian ⚪ **Backlog**

- **이전 Paper 3였음 → 2026-05-04 ICI에 Paper 3 슬롯 양보, Paper 4로 renumber.**
- **Scope:** Korean Graves' disease HLA + pan-Asian forest. Chu 2018 GD 비교는 Paper 2가 HT-only로 좁아진 후 **Paper 4로 이동**.
- **Entry condition (4/4 필요):** Paper 1 bioRxiv submitted + Paper 2 Pillar I+II+III committed + Bundang FFPE n>50 + Yu green light.
- **Status:** **No touch.** Forbidden words active. Paper 3 freeze 동안 작업 금지.
- **앵커 메모리:** `v19_paper4_GD_backlog`, `v19_paper3_GD_gating` (legacy).

---

## 3. 마라톤 모드 규칙 (위반 회피)

- **Window:** 2026-05-04 → 2026-06-13 (6주). 새 분석은 **paper-blocking only**.
- **Default 권한:** scaffolding / infra / audit. **Voice-protected 6 sections는 본인 키보드 strict.**
- **금지 키워드:** "고고" / "faster" / "다 해줘" → Claude의 voice-protected sprint generate 권한 **NOT**.
- **Track B unlock:** Paper 3에 한해 literal **"Paper 3 Track B 시작"** 만 인정.
- **컴퓨트 spend:** paper-blocking 만. 일반 audit/scaffolding은 read-only OK.
- **마라톤 displaced 시간:** image-DM1 NO-GO로 ~4–6주 freed → Paper 2B Pillar II/III 정의에 사용 가능.

---

## 4. 외부 의존성 / 미해결 결정

| 항목 | 상태 | 결정 필요? |
|---|---|---|
| Bundang SNUH outreach (이메일 5/3 final draft) | 발송 완료 여부 미확인 | 발송 상태 확인 |
| dbGaP / EGA 신청 (Paper 3 controlled access) | 미신청. todo/ checklist 생성됨 | Track B 시작 결정 후 |
| AFND South Korea pool (Paper 2B Pillar I v2 baseline) | 데이터 풀 정의 상태 미확인 | Task A 재확인 |
| GSE250521 Paper 1 inclusion | verdict 파일 존재, 결과 미확인 | 본문 read 필요 |
| Lee2024 / Lim2025 bib 저자·제목 (Paper 1 M3, M4) | GEO/PubMed lookup 부분 진행 | ~20분 잔여 |
| Voice-protected 6 sections (Paper 1 + Paper 2B) | scaffold만, 본인 작성 미시작 | 본인 일정 |
| RunPod Pod C (AlphaFold) + Pod D (K2 STAR) | 살아있음, dispatch 대기 | paper-blocking 여부 확인 후 stop/let |

---

## 5. 다음 행동 추천 우선순위

1. **Paper 1 bib M3/M4 마무리** (~20분, 마라톤-safe scaffolding).
2. **Paper 1 voice-protected sections** 본인 키보드 (Hook → Aim → Discussion §3.1 [Landa 2016 cite save] → Limitations → Cover ¶1 → Q9).
3. **Paper 2B brief audit** (`2026_05_04_paper2_post_image_dm1_nogo_status.md` §4 grep + image-DM1 wording 정리).
4. **Paper 2B Pillar II/III 정의** (Yu 다음 미팅 input 후).
5. **Paper 1 supplement methods scaffolding** for closure battery negative-feasibility entry (Claude OK, 마라톤-safe).

**Do NOT during marathon:** image-DM1 재분석 / foundation model retry / TCGA WSI 다운로드 / Paper 3 Track B 시작 / Paper 4 touch / voice-protected prose Claude generate.

---

## 6. 핵심 파일 인덱스

### Paper 1 manuscript
- `project/reports/PAPER1_MARATHON_AUDIT_BUNDLE_2026_05_04.md` (chmod 444)
- `project/reports/2026_05_04_paper1_HM_closure_report.md`
- `project/reports/2026_05_03_manuscript_v8_OUTLINE.md`
- `project/reports/2026_05_03_results_R1_R5_prose.md`
- `project/reports/2026_05_03_methods_M1_M11_prose.md`
- `project/reports/2026_05_03_figure_captions_all.md`
- `project/reports/2026_05_03_paper_supp_tables_draft.md`
- `project/reports/2026_05_03_reviewer_QA_consolidated.md`
- `project/reports/2026_05_03_PRE_SUBMISSION_CHECKLIST.md`
- `project/reports/2026_05_03_cover_letter_assembly.md`

### Paper 2A 실패 기록 (preserved)
- `project/reports/2026_05_04_image_dm1_final_nogo_decision.md` ← **decision authority**
- `project/reports/pathology_dm1_closure_battery_2026_05_04.md`
- `project/reports/pathology_dm1_phaseA_cpu_verdict_2026_05_04.md`
- `CLOSURE_BATTERY_2026_05_04.md` (root, Korean summary)
- `PHASE_A_NOGO_REPORT_2026_05_04.md` (root)
- `project/results/03_pathology_poc/` (8 결과 파일 아카이브)

### Paper 2B
- `project/reports/2026_05_03_paper2_pillar1_for_web_claude.md`
- `project/reports/2026_05_04_paper2_post_image_dm1_nogo_status.md`
- 메모리: `v18_paper2_HT_isolated`, `v17_paper2_pillar1_forest_strong`, `v17_GSE286332_strong_go`, `v17_D4P2_tcga_hashimoto_generalization`, `v17_D5P6_BCR_clonal_TLS`, `v17_D6P7_dm1_subB_NBNR`

### Paper 3 ICI freeze
- `project/reports/PAPER3_ICI_TRACK_A_FREEZE_2026_05_04.md`
- `project/reports/PAPER3_ICI_DECISION_BRIEF_2026_05_04.md`
- `project/reports/paper3_ici/PAPER3_ICI_TRACK_A_BUNDLE.md` (chmod 444, 1236 lines)
- `project/reports/paper3_ici/paper3_ici_track_a.tar.gz` (51 KB)
- `project/reports/paper3_ici/paper3_ici_*.md` × 8 (chmod 444)
- `project/reports/paper3_ici/todo/*.md` × 3 (chmod 644)

### Paper 4 backlog
- 메모리: `v19_paper4_GD_backlog`, `v19_paper3_GD_gating` (legacy reference)

### 메모리 인덱스
`~/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/MEMORY.md`

---

**End of overview.** 5개 카운트 = Paper 1 / Paper 2A (실패) / Paper 2B / Paper 3 (frozen) / Paper 4 (backlog). 한 개 실패 = **Paper 2A image-DM1**.
