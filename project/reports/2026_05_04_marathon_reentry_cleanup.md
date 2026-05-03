# Marathon re-entry cleanup (2026-05-04)

**Author:** Seungho Cook
**Date:** 2026-05-04 (post-image-axis drop, marathon mode 5/4–6/13)
**Mode:** read-only triage + small docs commit. Diff-friendly. No analysis.
**Voice-risk:** ZERO.

---

## 1. Paper status (active scope only)

| Paper | Status | Active scope |
|---|---|---|
| **Paper 1** — DM1 molecular dark matter | 🟢 ACTIVE | Molecular axis (TCGA bulk r=−0.885 + spatial replicate, Cox PFI HR=2.04, 4-step TF mechanism, GSEA Hallmark, TF coordination, drop-one-out, Moran's I, pan-cancer specificity, tumor-population TROP2, DM1/DM2 cluster bridge) |
| **Paper 2** — H&E → DM1 | 🔴 NO-GO 폐기 | Phase A NO-GO + closure battery NO-GO. Manuscript에서 제외. RunPod G1/G2/G3 폐기. |
| **Paper 3** — ICI vulnerability | 🔒 FROZEN | Track A bundle chmod 444. Track B BLOCKED (3 unlock 조건). |
| **Paper 4** — Korean GD HLA | 🔒 BACKLOG | 4/4 entry condition gated. Touch 금지. |

**Image-axis 어떤 결과물도 Paper 1 manuscript / supplement / Methods 에 포함하지 않음.** TROP2 framing은 tumor-population vulnerability only (Q3 spot-level ρ=−0.016 negative).

---

## 2. Working tree 잔여물 분류

### A. Commit as report (이 commit 포함)

| Path | Size | Why |
|---|---|---|
| `SHARE_2026_05_04_SESSION_RECAP.md` | 9.9 KB | 오늘 세션 8 commit + 4-paper 현황 1-page 공유본 |
| `SITUATION_BACKGROUND_ACTIVITY_2026_05_04.md` | 11 KB | 03:00–01:25 KST 백그라운드 dispatch 상황 분석 (kill-dispatch 결정 근거) |
| `project/reports/2026_05_04_marathon_reentry_cleanup.md` | (this file) | 본 cleanup report |

### B. Already committed elsewhere

| Path | Commit | Note |
|---|---|---|
| `project/reports/2026_05_04_paper1_dm1_full_conflict_audit.md` | `f66a6bb` | 평행 세션이 이미 commit. 추가 action 불필요. |

### C. Leave untracked (image-axis / RunPod 영역 — user decision pending)

| Path | Reason |
|---|---|
| `project/notebooks_or_scripts/post_pod_B_wsi_analysis.py` | image-axis WSI post-processing (Pod B null) |
| `project/notebooks_or_scripts/post_pod_C_alphafold_analysis.py` | RunPod C 후처리 (사용자 stop 권한 영역) |
| `project/notebooks_or_scripts/post_pod_D_meta3cohort.py` | RunPod D 후처리 (사용자 stop 권한 영역) |
| `project_external_st/src/12_gpu/` (8 파일, 56 KB) | image-axis GPU pipeline (NO-GO 후 deferred) |
| `project_external_st/results/extra/g1_external_tile_metadata.tsv.gz` (116 KB) | image-axis result artifact (NO-GO) |

### D. Leave untracked (user-side parallel session 산출물)

| Path | Reason |
|---|---|
| `PAPERS_5_OVERVIEW_2026_05_04.md` (root) | 사용자/평행 세션 작성 5-paper status snapshot — user commit list 외 |
| `POD_CD_RUNNING_STATE_2026_05_04.md` (root) | Pod C/D 상태 + execution plan — user commit list 외 |
| `POD_CLEANUP_AND_PAPER1_LOCK_2026_05_04.md` (root) | Pod cleanup + Paper 1 lock session report — user commit list 외 |
| `project/manuscript_v8/_session_origin_paper2_isolated.md` | 평행 세션 origin record |
| `project/results/terminology_correction_2026_05_04/SESSION_REPORT.md` | 평행 세션 산출물 |

이 D 그룹은 사용자가 직접 commit / 정리 권한. 본 cleanup commit 에 포함하지 않음.

### E. Modified files (M) — 평행 세션 변경 (이 commit 미포함)

| Path | Change | Reason for skip |
|---|---|---|
| `project/manuscript_p2_brief/p2_advisor_discussion.html` | 2 lines | Paper 2 brief HTML rebuild — paper2-pause 상태에서 평행 세션이 실행. 사용자 검토 후 |
| `project/results/terminology_correction_2026_05_04/before_after_diff.md` | +49 lines | 평행 세션 산출물 |
| `project/results/terminology_correction_2026_05_04/files_modified.txt` | +5 lines | 평행 세션 산출물 |
| `project/submission/papers_overview.html` | 4 lines | overview rebuild |
| `project/three_papers_index.html` | 14 lines | 3-paper index rebuild |

이 5개 M 파일은 본 cleanup 의 scope 외 (image-axis drop / Paper 1 lock 와 직접 연관 없음). 사용자가 검토 후 별도 commit.

### F. Permanent gitignore (이미 적용)

| Path | Where ignored |
|---|---|
| `project/.runpod_pods.json` | commit `e4bb222` |
| `project/.runpod_ssh.json` | commit `1d16a4e` (사용자) |
| `azure_gpu_phaseA/phaseA_gpu_pkg.tar.gz` | commit `25d693c` |
| `embeddings_*.npz`, raw data, `*.h5ad`, `*.fastq.gz`, `*.bam` | repo .gitignore 기존 rule |

---

## 3. NOT touched this cleanup

- ❌ RunPod API (Bearer token 미보유; Pod C/D stop 사용자 직접)
- ❌ H&E-DM1 retry (closure battery NO-GO 인용만)
- ❌ G1/G2/G3 RunPod sprint 실행 (모두 폐기)
- ❌ TCGA WSI download
- ❌ Paper 3 design bundle (chmod 444)
- ❌ Paper 4 backlog
- ❌ 새 분석 / 새 데이터 download / 새 sprint
- ❌ Voice-protected 6 sections (Hook ¶1, Aim ¶4, Discussion §3.1, §3.4 Limitations, Cover ¶1, Reviewer Q9)
- ❌ Modified HTML files / terminology correction (평행 세션 영역)
- ❌ User-written D-group MDs (PAPERS_5_OVERVIEW, POD_CD_RUNNING_STATE, POD_CLEANUP_AND_PAPER1_LOCK)

WebFetch 0회. 모든 데이터는 disk 기존 산출물 / 메모리 reference.

---

## 4. Marathon discipline summary (post-cleanup)

| Constraint | Status |
|---|---|
| Voice-protected sections | ✓ untouched this session |
| Paper 3 design bundle (chmod 444) | ✓ untouched |
| Paper 3 Track B | ✓ NOT started |
| Paper 4 (Korean GD HLA backlog) | ✓ untouched |
| New analysis launched | ✓ none |
| New data download | ✓ none |
| H&E-DM1 retry | ✓ NOT attempted (closure battery NO-GO 인용만) |
| TCGA WSI download | ✓ NOT attempted |
| RunPod API used by Claude | ✓ NOT used |
| `.runpod_pods.json` modified | ✓ NOT modified |
| `12_gpu/` pipeline committed | ✓ NOT committed (deferred) |
| post_pod_*.py committed | ✓ NOT committed (deferred) |

---

## 5. Paper 1 voice-hook readiness check

전제 조건:
- ✅ Manuscript v8 scaffold (9 files) 모두 disk 에 존재 + audit 통과
- ✅ HIGH 2/2 + MEDIUM 4/6 + LOW 6/9 closed (M3/M4/L1–L6)
- ✅ M3 (Lee SE 36 author Nat Comm 2024) + M4 (Lim 8 author Front Endocrinol 2025) bib 검증 완료
- ✅ L1–L6 cosmetic 완료 (ARI 3-decimal / Fisher p / BRAF Cohen d / bib header / ★→ASCII / Pan2025 article# 3601)
- ✅ Image-axis dropped (PAPER1_DM1_FULL §9–11 DEPRECATED 태그 + §12 manuscript-safe candidate inventory)
- ✅ TROP2 reframed (tumor-population vulnerability binding rule)
- ✅ Limitations seed 충분 (§0 6 risks + §2 N1 + §7-Q3 negative)
- ✅ Discussion §3.1 cite save anchor (memory `v17_landa2016_cite_save`)
- ⏳ Voice-protected 6 sections 본인 키보드 슬롯 대기

**voice-hook 진입 가능.** Hook ¶1 부터 본인 키보드.

---

## 6. Next recommended action

**`voice-hook`** ⭐

User opens `project/reports/2026_05_03_manuscript_v8_OUTLINE.md` and writes Hook ¶1 directly. Marathon 첫 voice-protected slot. Claude 는 이 시점 default = scaffolding/infra만; "고고"/"faster"/"다 해줘" 가 voice 권한 아님 (per `v17_sprint_vs_marathon_violation`).

Claude 가 도울 수 있는 것:
- voice draft 후 typo / fact-check / xref 정합성 audit
- supp table / figure caption 의 voice-after edits
- Discussion §3.1 cite save 적용 시 bib + xref 검증

Claude 가 못 하는 것:
- Hook / Aim / Discussion / Limitations / Cover / Q9 의 sprint generate
- Pod C/D stop (user RunPod console)
- `12_gpu/`, `post_pod_*` 처분 결정 (user)
- 평행 세션 modified HTML 검토 (user)

---

Cleanup closed. Marathon re-entry 가능.
