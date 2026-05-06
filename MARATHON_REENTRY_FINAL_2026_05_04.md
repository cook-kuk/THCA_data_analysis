# Marathon re-entry FINAL — 2026-05-04 session-end consolidated view

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Latest commit:** `742d644` · **Session phase:** Marathon re-entry final cleanup → ready for voice-hook
**Single-file shareable summary** for VS Code (Ctrl+Shift+V preview).

---

## 0. 한 줄 요약

> Paper 1 molecular axis는 ship-ready. Image-DM1 angle은 closure battery (D 폐기) + decision memo + lock memo로 완전히 차단. 이번 session에서 RunPod stop reminder + final classification memo commit 완료. **다음 step = 본인 키보드로 voice-hook ¶1.** Pod C/D는 user manual stop 필요. Claude는 RunPod API 안 쓴다.

---

## 1. 오늘 commit chain (최신 → 과거)

| # | hash | 메시지 | 핵심 |
|---|---|---|---|
| 9 | **`742d644`** | docs: finalize marathon reentry status and runpod stop reminder | **이 session** — 2 memos (runpod stop reminder + final classification) |
| 8 | `95676c2` | docs: record marathon reentry cleanup after image-axis drop | parallel — SHARE_2026_05_04_SESSION_RECAP commit |
| 7 | `f66a6bb` | docs: lock paper1 molecular-only interpretation after dm1 full audit | 직전 session — molecular_only_lock + leftover_inventory + conflict_audit + SITUATION_BACKGROUND_ACTIVITY |
| 6 | `153084e` | chore: tag paper1 dm1 bundle for marathon scope (image-axis dropped) | parallel — PAPER1_DM1_FULL banner+tags |
| 5 | `c42bd6a` | fix: close paper1 marathon audit low batch L1-L6 | parallel — L1 ARI / L2 Fisher p / L3-L6 |
| 4 | `1d16a4e` | decision: drop image-dm1 pathology angle after closure battery | image-DM1 angle 정식 폐기 + Paper 2 scope cleanup + marathon snapshot |
| 3 | `778b255` | decision: archive image-dm1 closure battery no-go | parallel — closure battery 결과 보존 |
| 2 | `e4bb222` | chore: resolve bib medium issues and runpod stragglers | parallel — bib M3/M4 |
| 1 | `25d693c` | chore: archive spatial and pathology feasibility artifacts | parallel — 205 파일 archive |

오늘 추가된 commit 9개. Marathon 규칙 무위반.

---

## 2. 4-Paper 현재 상태 (2026-05-04 EOD)

| paper | 상태 | image-DM1 영향 | 다음 |
|---|---|---|---|
| **Paper 1** — DM1 molecular dark matter | **active / ship-ready.** manuscript v8 scaffold + audit + HM closure + bib M3/M4/L1-L6 모두 done. Voice-protected 6 sections 본인 키보드 대기 | image-axis dropped — Paper 1 contribution 아니었음 | Hook ¶1 → Aim ¶4 → §3.1 (Landa 2016 cite) → §3.4 Limitations → Cover ¶1 → Q9 → bioRxiv 6/13 |
| **Paper 2** — H&E → DM1 / pathology projection | **paused.** Phase A NO-GO + closure battery NO-GO 두 layer 폐기 확정 | dropped (이게 그 paper) | Paper 1 ship 후 재정의 |
| **Paper 3** — ICI vulnerability dark thyroid cancer | **Track A frozen** (chmod 444 bundle + tar.gz). Track B blocked | none | Paper 1 bioRxiv + Paper 2 close + 명시적 "Paper 3 Track B 시작" 명령 |
| **Paper 4** — Korean GD HLA / Pan-Asian | **backlog** | none | 4/4 entry condition gated |

(5-paper schema는 `PAPERS_5_OVERVIEW_2026_05_04.md` 참조 — parallel agent file, root에 untracked)

---

## 3. 이 session 산출물 (commit `742d644` 안)

```
project/reports/2026_05_04_runpod_user_stop_reminder.md           (93 lines)
project/reports/2026_05_04_final_untracked_classification.md      (99 lines)
```

### `2026_05_04_runpod_user_stop_reminder.md` — 핵심
- Claude는 RunPod API 사용 금지 (이 session + 이후 marathon)
- Pod C / Pod D 둘 다 user 직접 stop 필요
- 24h auto-shutdown은 safety net, 즉시 stop 권장
- Marathon 동안 추가 compute 조건 4개 (모두 충족 시만): paper-blocking + cached 안 됨 + user 명시 승인 + closed angle 아님

### `2026_05_04_final_untracked_classification.md` — 핵심
- A 3건 commit (이미 처리됨)
- B 13건 leave-untracked deprecated (12_gpu/, post_pod_*, g1 metadata)
- C 14+건 user decision (P2 brief, papers_overview, parallel agent files)
- D never commit (data/embeddings/runtime/raw)

---

## 4. RunPod Pod C / Pod D — USER ACTION 필요

| pod | id | role | endpoint | shutdown |
|---|---|---|---|---|
| **Pod C** | `158ic3wrtf2j8l` | AlphaFold | 216.81.151.3:10960 | 24h auto-shutdown 설정됨 |
| **Pod D** | `8kdsltl8s2dqbo` | K2 STAR meta-3-cohort | 157.157.221.29:21123 | 24h auto-shutdown 설정됨 |

### Stop 명령 (user 직접 실행)

**Option 1 — RunPod web console (권장)**
- https://www.runpod.io/console/pods 접속
- Pod C, Pod D 각각 선택 → Stop

**Option 2 — runpod CLI**
```bash
runpod stop pod 158ic3wrtf2j8l
runpod stop pod 8kdsltl8s2dqbo
```

**Option 3 — GraphQL via curl**
```bash
curl -X POST "https://api.runpod.io/graphql?api_key=$RUNPOD_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"query":"mutation { podStop(input: {podId: \"158ic3wrtf2j8l\"}) { id desiredStatus } }"}'
curl -X POST "https://api.runpod.io/graphql?api_key=$RUNPOD_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"query":"mutation { podStop(input: {podId: \"8kdsltl8s2dqbo\"}) { id desiredStatus } }"}'
```

**Pod B**는 dispatch log에서 `null` (한 번도 시작 안 됨 — image-DM1이라 closure battery로 차단됨). cleanup 불필요.

---

## 5. Voice-hook 진입 준비 — `voice_hook_fact_brief.md`

**위치:** `project/reports/2026_05_04_voice_hook_fact_brief.md` (uncommitted)
**작성자:** parallel agent (Claude 작성한 버전보다 더 풍부)
**용도:** 본인 키보드로 Hook ¶1 작성할 때 참조용 fact-only briefing
**구조 (8 sections, bullets only)**:
1. Clinical problem (3 bullets) — DTC 15-35% recurrence, ATA 2015 anatomic-only stratification
2. Why BRAF/RAS exhausted (3) — Cohen d / single-feature AUC / driver_anchor ARI=−0.007
3. Dark-matter gap (3) — N=156 BRAF-neg ∩ RAS-neg, DM2 unexplained, prior frameworks
4. **What this paper resolves (5 bullets)** — 8-gene panel + Pan-Asian HLA + GSE286332 + TCGA Hashimoto-like + DM1 sub-B + Survival HR
5. What NOT to overclaim (5 bullets) — verbs / "all risks resolved" / image / TROP2 / venue %
6. **Candidate opening angles A/B/C/D** (label-only, no prose)
7. **Exact numbers SAFE in Hook** (table, ~22 rows with sources)
8. **Exact numbers NOT safe in Hook** (table, 11 rows with reason + proper home)

VSCode 명령:
```bash
code project/reports/2026_05_04_voice_hook_fact_brief.md
# Ctrl+Shift+V → markdown preview
```

---

## 6. 다음 사용자 액션 우선순위

### A. 즉시 (5분)
- **Pod C, Pod D stop** (위 §4 명령)

### B. Voice-hook 작성 진입 (본인 키보드, Claude 도움 X)
- `2026_05_04_voice_hook_fact_brief.md` open + preview
- §1-§5 bullets + §7 safe number table 보면서 Hook ¶1 prose 작성
- §6 4개 angle 중 하나 선택 (혹은 hybrid)
- §8 unsafe number table에 있는 표현은 **STOP**

### C. Voice-hook 끝나면 (다음 marathon block)
- Aim ¶4 → Discussion §3.1 (Landa 2016) → §3.4 Limitations → Cover Para 1 → Q9 (전부 본인 키보드)

### D. Cosmetic cleanup (언제든 OK, 마라톤 우선순위 낮음)
- C class user-decision items review:
  - Paper 2 brief HTML/PDF (3 files) — P2 paused이라 commit 보류 권장
  - papers_overview HTML/PDF — parallel agent의 update
  - three_papers_index HTML/PDF
  - manuscript_v8/_session_origin_paper2_isolated.md
  - root memos: POD_CD_RUNNING_STATE / PAPERS_5_OVERVIEW / POD_CLEANUP_AND_PAPER1_LOCK / 이 파일 (MARATHON_REENTRY_FINAL)

### E. Claude 추가 가능 (voice-protected 외, 마라톤 호환)
- Paper 1 supplementary Methods scaffolding (closure battery negative-feasibility entry)
- Figure caption scaffolding (lock §6 envelope에서)
- Reviewer Q1-Q8 scaffolding (Q9 voice-protected 제외)

### F. Claude 절대 안 함
- Hook / Aim / Discussion 3.1 / Limitations / Cover Para 1 / Q9 prose
- RunPod API · GPU spend · 새 분석 · 새 데이터 다운로드 · image-DM1 retry · foundation model · TCGA WSI · Paper 3/4 touch

---

## 7. 모든 관련 파일 한눈에

### 이 session commit `742d644` (commit 안에 있음)
- `project/reports/2026_05_04_runpod_user_stop_reminder.md`
- `project/reports/2026_05_04_final_untracked_classification.md`

### 직전 session commit `f66a6bb` (이미 commit 됨)
- `project/reports/2026_05_04_paper1_dm1_full_molecular_only_lock.md` ⭐ Paper 1 manuscript-safe envelope authority
- `project/reports/2026_05_04_paper1_dm1_full_conflict_audit.md` ⭐ claim-by-claim source verification (16/16)
- `project/reports/2026_05_04_leftover_untracked_inventory.md`
- `SITUATION_BACKGROUND_ACTIVITY_2026_05_04.md`

### 이전 session commit `1d16a4e` (이미 commit 됨)
- `project/reports/2026_05_04_image_dm1_final_nogo_decision.md` ⭐ closure decision authority (5 re-entry conditions)
- `project/reports/2026_05_04_paper2_post_image_dm1_nogo_status.md`
- `project/reports/2026_05_04_marathon_state_after_image_dm1_nogo.md`
- closure battery outputs (24 MB total, all in `project/results/03_pathology_poc/`)

### Parallel agent (이미 commit 됨)
- commit `95676c2`: `SHARE_2026_05_04_SESSION_RECAP.md` (root)
- commit `153084e`: `PAPER1_DM1_FULL_2026_05_04.md` (root) with banner + per-section [MANUSCRIPT-SAFE] / [INTERNAL-ONLY] / [DEPRECATED] tags

### Uncommitted, ready (user 결정 대기)
- `project/reports/2026_05_04_voice_hook_fact_brief.md` ⭐ **voice-hook fact briefing — 본인 키보드용**
- `MARATHON_REENTRY_FINAL_2026_05_04.md` (이 파일)
- `POD_CLEANUP_AND_PAPER1_LOCK_2026_05_04.md` (root, 이전 session)
- `POD_CD_RUNNING_STATE_2026_05_04.md` (root, parallel agent)
- `PAPERS_5_OVERVIEW_2026_05_04.md` (root, parallel agent)
- Paper 2 brief HTML/PDF + papers_overview + three_papers_index + manuscript_v8 origin

### Closure battery 핵심 결과 (`project/results/03_pathology_poc/`, 이미 commit)
- `closure_battery_metrics.tsv` (26 rows A/B/D/E/F)
- `closure_battery_predictions.tsv.gz` (per-tile preds)
- `closure_battery_summary.png` (bar plot)
- `loso_metrics_resnet50.tsv` (224 + 448)
- `loso_predictions_resnet50.tsv.gz`
- `negative_controls_summary.tsv` (resid mode)
- `negative_controls_raw_summary.tsv` (raw mode)
- `pred_vs_obs_resnet50.png`, `pred_vs_obs_resnet50_per_slide.png`
- `tile_metadata*.tsv.gz` (8,521 tiles 224/448/672)

### Top-level summary memos (root, scp-ready)
- `MARATHON_REENTRY_FINAL_2026_05_04.md` ⭐ **이 파일**
- `POD_CLEANUP_AND_PAPER1_LOCK_2026_05_04.md`
- `SHARE_2026_05_04_SESSION_RECAP.md` (commit 95676c2)
- `SITUATION_BACKGROUND_ACTIVITY_2026_05_04.md` (commit f66a6bb)
- `CLOSURE_BATTERY_2026_05_04.md` (commit 778b255)
- `PHASE_A_NOGO_REPORT_2026_05_04.md` (commit 25d693c)
- `PATHOLOGY_DM1_FEASIBILITY_2026_05_03.md` (commit 25d693c)
- `PAPER1_DM1_FULL_2026_05_04.md` (commit 153084e, banner+tags)

---

## 8. scp 명령 (Windows PowerShell)

```powershell
# 이 단일 요약 파일 (모든 것을 보여주는 한 파일):
scp seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/MARATHON_REENTRY_FINAL_2026_05_04.md .

# Voice-hook fact briefing (본인 키보드 작성용):
scp seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/project/reports/2026_05_04_voice_hook_fact_brief.md .

# 이번 session commit의 2 memos:
scp seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/project/reports/2026_05_04_runpod_user_stop_reminder.md .
scp seungho@40.82.129.913:/home/seungho/personal/THCA_data_analysis/project/reports/2026_05_04_final_untracked_classification.md .

# Paper 1 lock + audit (직전 session):
scp seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/project/reports/2026_05_04_paper1_dm1_full_molecular_only_lock.md .
scp seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/project/reports/2026_05_04_paper1_dm1_full_conflict_audit.md .

# Image-DM1 closure decision authority:
scp seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/project/reports/2026_05_04_image_dm1_final_nogo_decision.md .

# 모든 reports 한 번에:
scp -r seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/project/reports .
```

VSCode 열기:
```bash
code MARATHON_REENTRY_FINAL_2026_05_04.md
# Ctrl+Shift+V → markdown preview
```

---

## 9. 한 줄 next action

**Pod C/D 즉시 stop (위 §4) → `2026_05_04_voice_hook_fact_brief.md` 열기 → 본인 키보드로 Hook ¶1 작성.** Claude는 prose 한 줄도 안 씀.

---

*Generated 2026-05-04 by Claude (Opus 4.7) under marathon-mode discipline. No new analysis. No RunPod API. No Pod start/stop. No image-DM1 retry. No new data download. No voice-protected prose. Pod C/D observed read-only and left to user manual stop.*
