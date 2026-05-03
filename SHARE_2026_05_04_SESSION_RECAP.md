# THCA 3-Paper Marathon — 2026-05-04 Session Recap

> Single-file shareable summary. VS Code: open in preview (Ctrl+Shift+V).
> **Period:** 2026-05-04 dawn → afternoon. **Mode:** Marathon 5/4–6/13 (Paper 1 bioRxiv 6/13).
> **Author:** Seungho Cook.

---

## TL;DR

오늘 세션은 (1) 새벽 크래시 복구 → (2) 4단계 audit·closure → (3) image-DM1 pivot 결정 → (4) Paper 1 manuscript-safe scope lock 으로 끝남. **Paper 1 molecular axis 는 ship-ready, image axis 는 NO-GO 폐기.** 8커밋 추가, marathon 규칙 무위반.

---

## 1. 오늘 commits (시간순)

| # | Hash | 메시지 | 핵심 |
|---|---|---|---|
| 1 | `c55d582` | freeze: paper3 ici track a design bundle | Paper 3 ICI Track A 디자인 8 deliverables 동결 (chmod 444) |
| 2 | `ef023dd` | fix: close paper1 marathon audit high medium issues | HIGH 2/2 + MEDIUM 4/6 closed (sub-B 96%→94.6%, S1/S3/S8 빌드) |
| 3 | `25d693c` | chore: archive spatial and pathology feasibility artifacts | 205 파일 archive (reports/code/small results, 대용량 .gitignore) |
| 4 | `e4bb222` | chore: resolve bib medium issues and runpod stragglers | M3 Lee2024 (Lee SE 36 author + Nat Comm 2024) + M4 Lim2025 (full title + 8 author) + .runpod_pods.json gitignore + pod_B sh |
| 5 | `778b255` | decision: archive image-dm1 closure battery no-go | H&E→DM1 closure battery (D 폐기) + dispatch shutdown 리포트 |
| 6 | `1d16a4e` | decision: drop image-dm1 pathology angle after closure battery | (사용자 또는 다른 세션에서 실행) image-DM1 angle 폐기 |
| 7 | `c42bd6a` | fix: close paper1 marathon audit low batch L1-L6 | L1 ARI / L2 Fisher p / L3 BRAF Cohen d / L4 bib header / L5 ★→ASCII / L6 Pan2025 article# 3601 |
| 8 | `153084e` | chore: tag paper1 dm1 bundle for marathon scope (image-axis dropped) | PAPER1_DM1_FULL 재태깅 (manuscript-safe / internal-only / DEPRECATED) + §12 candidate inventory |

---

## 2. 4-Paper 현황 (2026-05-04 EOD)

| Paper | 제목 | 상태 | 다음 |
|---|---|---|---|
| **Paper 1** | DM1 molecular dark matter | manuscript v8 scaffold + 모든 audit + HIGH/MEDIUM/LOW 마무리. **image-axis dropped.** Voice-protected 6 sections 본인 키보드 대기 | Hook ¶1 → Aim ¶4 → §3.1 (Landa 2016 cite save) → §3.4 Limitations → Cover ¶1 → Q9 → 6/13 bioRxiv |
| **Paper 2** | H&E → DM1 / pathology projection / TCGA validation (HT-only) | **paused.** Phase A NO-GO + closure battery NO-GO 둘 다 확정. 새 디자인 보류 | Paper 1 ship 후 재정의 |
| **Paper 3** | ICI vulnerability dark thyroid cancer | Track A FROZEN (chmod 444 bundle + tar.gz). Track B BLOCKED | Paper 1 bioRxiv + Paper 2 close + 명시적 "Paper 3 Track B 시작" 명령 |
| **Paper 4** | Korean GD HLA / Pan-Asian | backlog | 4/4 entry condition gated |

---

## 3. Paper 1 manuscript-safe scope (post-image-axis-drop)

### Main figure 후보 3개 (molecular axis only)

1. **DM1/RAI cross-validation** — TCGA bulk r=−0.885 + spatial replicate 2-platform consilience
2. **Clinical outcome (Cox)** — dichotomized PFI HR=2.04 (p=0.015) + multivariate DFI HR=1.41 (p=0.025) [+ N1 dark matter NS limitation 명시]
3. **TF-driven 4-step mechanism** — NKX2-1/FOXE1 collapse → DNMT activation → STAT3/AP-1 activation → (TROP2 elevation, tumor-level only); GSEA Hallmark + TF coordination 으로 cross-checked

### Supplementary 후보 8개

drop-one-out 8 RAI / Moran's I + bivariate / pan-cancer THCA outlier / tumor vs normal (**tumor-population vulnerability** — DM1-subset-specific 아님) / DM1-DM2 cluster bridge / drug-target volcano (tumor-level만) / Q3 TROP2 spot-level NEGATIVE (정직 보고) / GSEA Hallmark IL6_JAK_STAT3+EMT+IFN-γ↑/OXPHOS↓

### Already-frozen spatial supplement

`project/supplementary/spatial_freeze_2026_05_03/` SuppFig X1·X2·X3 + SuppTable SX·SX1 (모두 molecular score, image-axis 아님 → 그대로 유지)

### Manuscript에서 금지

- "DM1-high spot = TROP2-high spot" 류 spot-level subset claim (Q3 ρ=−0.016 negative)
- "DM1-axis defines sacituzumab candidate population" 류 (wet-lab 검증 전)
- Image-axis 어떤 figure / supp / Methods 도 포함 금지 (closure battery NO-GO)
- venue probability 문장 (cover letter / prose / submission package 어디에도)

---

## 4. Voice-protected 6 sections (본인 키보드 슬롯 대기)

| Section | 메모리 anchor | 비고 |
|---|---|---|
| **Hook ¶1** | — | manuscript_v8_OUTLINE 첫 단락 |
| **Aim ¶4** | — | OUTLINE Aim 단락 |
| **Discussion §3.1** | `v17_landa2016_cite_save` | Krishnamoorthy 2025 → Landa 2016 cite 정정 + 3-layer reverse-causality 차단 |
| **Discussion §3.4 Limitations** | `PAPER1_DM1_FULL §0` 6항목 + §2 N1 + §7-Q3 negative | seed 풍부 |
| **Cover letter ¶1** | venue 결정 시점 | venue probability 표는 §3 internal-only |
| **Reviewer Q9** | — | manuscript_v8 Q&A consolidated |

규칙: "고고" / "faster" / "다 해줘" → voice-protected sprint 권한 NOT (per `v17_sprint_vs_marathon_violation`).

---

## 5. Marathon 위반 0건 attestation

- Voice-protected 섹션 미수정 ✓
- Paper 3 design bundle (chmod 444) 미수정 ✓
- Paper 3 Track B 미시작 ✓
- Paper 4 미터치 ✓
- 새 분석 0건 ✓
- 새 데이터 download 0건 (WebFetch 4회 read-only: PMID 38331894 / 41113708 / 40234451 + GEO GSE213647 — 모두 citation 검증용)
- H&E-DM1 retry 미시도 ✓ (closure battery NO-GO 인용만)
- TCGA WSI download 미시도 ✓
- RunPod API Claude 사용 0건 ✓ (token 미보유)
- `.runpod_pods.json` 미수정 ✓ (gitignored)

---

## 6. ⚠ 사용자 action 대기

### RunPod Pod C / Pod D stop

| Pod | ID | Role | SSH |
|---|---|---|---|
| C | `158ic3wrtf2j8l` | AlphaFold | 216.81.151.3:10960 |
| D | `8kdsltl8s2dqbo` | K2 STAR | 157.157.221.29:21123 |

**24시간 auto-shutdown forced** 됐지만 즉시 stop 권장 (marathon mode 정신):

```bash
RUNPOD_KEY=$(grep RUNPOD_API_KEY ~/.runpod/config | cut -d= -f2)
for ID in 158ic3wrtf2j8l 8kdsltl8s2dqbo; do
  curl -s -X POST https://api.runpod.io/graphql \
    -H "Authorization: Bearer $RUNPOD_KEY" \
    -H "Content-Type: application/json" \
    -d "{\"query\":\"mutation { podStop(input: {podId: \\\"$ID\\\"}) { id desiredStatus } }\"}"
done
```

### Working tree 잔여물 (decision pending)

```
?? SITUATION_BACKGROUND_ACTIVITY_2026_05_04.md            (이 세션 정황 분석 MD)
?? SHARE_2026_05_04_SESSION_RECAP.md                       (이 파일)
?? project/notebooks_or_scripts/post_pod_B/C/D_*.py        (3 post-processing scripts)
?? project/reports/2026_05_04_paper1_dm1_full_conflict_audit.md  (사용자 작성)
?? project_external_st/results/extra/g1_external_tile_metadata.tsv.gz
?? project_external_st/src/12_gpu/                         (8 GPU pipeline files, image-axis NO-GO 후 deferred)
```

→ `12_gpu/` + `post_pod_*` 은 image-axis NO-GO 와 paper2-pause 둘 다 적용 → commit 보류 또는 archive 결정 필요.

---

## 7. 핵심 reference 파일

### 마라톤 audit 산출물 (오늘)

- `project/reports/PAPER1_MARATHON_AUDIT_BUNDLE_2026_05_04.md` (chmod 444, 51KB)
- `project/reports/2026_05_04_paper1_HM_closure_report.md`
- `project/reports/2026_05_04_post_closure_recheck.md`
- `project/reports/2026_05_04_bib_m3m4_straggler_cleanup.md`
- `project/reports/2026_05_04_low_batch_L1_L6_closure.md`
- `project/reports/2026_05_04_background_dispatch_shutdown.md`
- `project/reports/pathology_dm1_closure_battery_2026_05_04.md`
- `CLOSURE_BATTERY_2026_05_04.md` (root, final 1-page decision)
- `PAPER1_DM1_FULL_2026_05_04.md` (root, audited + tagged)

### Paper 3 ICI freeze

- `project/reports/PAPER3_ICI_TRACK_A_FREEZE_2026_05_04.md`
- `project/reports/PAPER3_ICI_DECISION_BRIEF_2026_05_04.md`
- `project/reports/paper3_ici/PAPER3_ICI_TRACK_A_BUNDLE.md` (chmod 444, 1236 lines)
- `project/reports/paper3_ici/paper3_ici_track_a.tar.gz` (51 KB)

### Paper 1 manuscript v8 scaffold

- `project/reports/2026_05_03_manuscript_v8_OUTLINE.md`
- `project/reports/2026_05_03_results_R1_R5_prose.md`
- `project/reports/2026_05_03_methods_M1_M11_prose.md`
- `project/reports/2026_05_03_figure_captions_all.md`
- `project/reports/2026_05_03_paper_supp_tables_draft.md`
- `project/reports/2026_05_03_reviewer_QA_consolidated.md`
- `project/reports/2026_05_03_PRE_SUBMISSION_CHECKLIST.md`
- `project/reports/2026_05_03_STAR_Methods_KeyResources.md`
- `project/reports/2026_05_03_cover_letter_assembly.md`
- `project/reports/2026_05_03_references.bib` (29 entries, M3+M4 fixed, L4/L5/L6 polished)

### Paper 1 supp data

- `project/results/p2_pillar1_forest/cohort_assembly.tsv` (5 cohorts × 10 cols)
- `project/results/p3_gse286332/SuppTable_S3_GSE286332_DEGs.xlsx` (29,672 DEGs)
- `project/results/d6p7_dm1_subcluster/SuppTable_S8_DM1_subcluster.xlsx` (4 sheets)

---

## 8. 다음 행동 (recommended)

| 우선순위 | 작업 | 비고 |
|---|---|---|
| 🔴 **즉시** | RunPod Pod C/D stop | 위 §6 명령 |
| 🟡 next | `voice-hook` Hook ¶1 작성 | 본인 키보드, 마라톤 첫 voice slot |
| 🟢 backlog | conflict_audit.md 검토 + 12_gpu/ 처분 결정 | 사용자 작성 audit 본 후 |
| 🟢 backlog | Voice §3.1 (Landa 2016 cite save) | Hook 후 |
| 🟢 backlog | Voice §3.4 Limitations | seed = §0 6 risks + §2 + Q3 |

---

## 9. 메모리 anchor (참고)

```
v17_marathon_mode_post_pillar1   — 마라톤 정의 (5/4–6/13)
v17_sprint_vs_marathon_violation — 고고/faster/다 해줘 ≠ voice-protected sprint
v17_landa2016_cite_save          — Discussion §3.1 cite 정정 3-layer
v18_paper2_HT_isolated           — Paper 2 = HT-only scope
v19_paper3_ici_track_a           — Paper 3 Track A FROZEN
v19_paper4_GD_backlog            — Paper 4 gating
paper_numbering_2026_05_04       — canonical 4-paper numbering
v17_D6P7_dm1_subB_NBNR           — DM1 sub-B = NBNR cluster
v17_K2_vs_bundang_distinction    — K2 / Bundang 분리
```

---

**End of recap.** 8 commits, marathon discipline preserved, Paper 1 ship-ready (molecular axis), image-DM1 dropped, Paper 3 frozen, Pod C/D awaits user stop.
