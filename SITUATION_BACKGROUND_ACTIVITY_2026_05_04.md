# 백그라운드 활동 상황 판단용 MD (2026-05-04 01:25)

**Author:** Seungho Cook
**작성:** 2026-05-04 01:25 (4번째 commit `e4bb222` 직후)
**목적:** 마라톤 모드 / paper2-pause 와 충돌 가능한 백그라운드 활동을 사용자 판단할 수 있도록 한 파일에 정리.
**Claude 처리:** **아래 모든 항목 read-only 관찰. 어떤 파일도 수정·삭제·kill 안 함.**

---

## 0. 한 줄 요약 + 결정 질문

**상황:** 이 세션 시작 직전 (01:14)에 사용자가 직접 시작한 RunPod 디스패치 루프가 실행 중이며, **마라톤 모드 안에서 컴퓨트 spend가 일어나고 있다.** 일부는 사용자 본인의 closure-of-closure 작업 (✅ 정당), 일부는 paper2-pause / H&E-DM1 금지 규칙과 모서리에 닿음 (⚠️ 판단 필요).

**결정 질문 3개:**

| Q | 질문 | 옵션 |
|---|---|---|
| **Q1** | `/tmp/runpod_dispatch_v2.sh` 디스패치 루프 (240 iter × 10s = 40분 max, etime 03:11) | **kill / let finish / selective** |
| **Q2** | RunPod Pod C (AlphaFold) + Pod D (K2 STAR) 실행 권한 | **계속 / 중단** |
| **Q3** | 새로 생긴 working tree 파일 9개 처리 | **commit / gitignore / leave-untracked / delete** |

답: §6 결정 메뉴에서 한 줄씩.

---

## 1. 핵심 발견 — `/opt/thyroid-dash` = 본 repo symlink

```
$ readlink /opt/thyroid-dash/project
/home/seungho/personal/THCA_data_analysis/project
```

→ 디스패치 스크립트가 `/opt/thyroid-dash/...` 경로로 쓰는 모든 결과가 **본 repo 작업 트리에 직접 흘러들어옴**. 이게 working tree에 untracked 파일이 계속 추가되는 메커니즘.

---

## 2. 활성 프로세스 (현재 시각 기준)

```
PID 185803  etime 03:11  /bin/bash /tmp/runpod_dispatch_v2.sh
PID 186051  etime 02:56  tail -F /tmp/dispatch.log
```

### 2.1 디스패치 스크립트 동작 (`/tmp/runpod_dispatch_v2.sh`)

```bash
# 240 × 10s = 40 분 max. 10초마다 RunPod GraphQL API 폴링.
# Pod D / C / B 가 SSH-ready 되는 즉시 dispatch.
# Bearer 토큰: ~/.runpod/config 에서 읽음.
```

### 2.2 디스패치 로그 (가장 최근)

```
[01:21:37] Starting dispatch loop
[01:22:38] iter=6,  dispatched=0/3 (D=0 C=0 B=0)
[01:23:46] iter=12, dispatched=0/3 (D=0 C=0 B=0)
```

**현재 상태:** 폴링 중, 아무 pod도 dispatch 안 됨. Pod 들이 SSH-ready 상태가 아님 (cold/booting/없음).

### 2.3 pod ID 변화 (.runpod_pods.json)

| Pod | 이전 (이 세션 초) | 현재 (01:22) | 의미 |
|---|---|---|---|
| D | `fbg10t9c5ftlcc` | `8kdsltl8s2dqbo` | **새 pod로 교체됨** (K2 STAR) |
| C | `eq1hcajjwa41wb` | `158ic3wrtf2j8l` | **새 pod로 교체됨** (AlphaFold) |
| B | `lrhd42ciql44wp` | `null` | **kill됨** (WSI — 다행) |

→ 사용자가 이 세션 동안 .runpod_pods.json을 직접 업데이트했음. **C, D 는 의도적으로 살아있음.** B는 의도적으로 죽어있음.

### 2.4 Cron / systemd 없음

```
$ crontab -l → no crontab
```

→ 스케줄링 없음. `runpod_dispatch_v2.sh`는 **사용자가 직접 foreground로 시작**한 것. 종료 = 자연스러운 240 iter 끝나거나 pod 모두 dispatch되거나 사용자가 kill.

---

## 3. 새로 생긴 9개 파일 — 시간순

| 시각 | 경로 | 출처 / 추정 의도 |
|---|---|---|
| 01:12–01:16 | `project_external_st/src/12_gpu/` (8 파일) | external_st GPU 파이프라인 새 디자인 (extract_external_tiles, g1_embed/g1_ridge_loso, g2_slide_regress/g2_tile_embed, g2_tcga_he_pipeline.sh, g3_celldart, setup.sh, README.md). **사용자 직접 작성 추정** (Pod C/D dispatch 대상 코드). |
| 01:13 | `project_external_st/results/extra/g1_external_tile_metadata.tsv.gz` | g1_embed 출력. 116 KB. |
| 01:16 | `project/results/03_pathology_poc/negative_controls_raw_summary.tsv` | closure battery output |
| 01:16 | `project/results/03_pathology_poc/tile_metadata_size672.tsv.gz` | closure battery output (tile size 변경) |
| 01:16 | `project/results/03_pathology_poc/tile_metadata_resid.tsv.gz` (M) | closure battery 재계산 (이전 버전 group 3에 commit됨) |
| 01:17 | `project/results/03_pathology_poc/tile_metadata.tsv.gz` (M) | 동일 |
| 01:16:53 | `project/notebooks_or_scripts/post_pod_D_meta3cohort.py` | Pod D (K2 STAR) 후처리 코드 |
| 01:17:25 | `project/notebooks_or_scripts/post_pod_C_alphafold_analysis.py` | Pod C (AlphaFold) 후처리 코드 |
| 01:17:57 | `project/notebooks_or_scripts/post_pod_B_wsi_analysis.py` | Pod B 후처리 코드 (Pod B는 null인데 코드만 작성) |
| **01:21:05** | `project/reports/pathology_dm1_closure_battery_2026_05_04.md` | **closure battery 결과 보고서** |

---

## 4. `pathology_dm1_closure_battery_2026_05_04.md` 핵심 (01:21 작성)

```markdown
# Closure battery — H&E → DM1_like_score_resid (failure-mode dissection)
**Run:** 2026-05-04 · **Owner:** Seungho Cook · **Wall time:** ~80 min CPU
**Goal:** Phase A NO-GO를 인정하고 폐기 전, 5가지 가능 원인을 구분.

> TL;DR — **D (full 폐기) 권장.** ResNet50 ImageNet 임베딩은 raw 신호조차 random panel을 못 이김.
> ... (3) depth-residualization이 morphology-visible signal을 제거하고 (4) spot-level ST label이 H&E와
> fundamentally mismatched가 dominant cause. Foundation model 업그레이드도 (3)+(4) 못 풂.
> **원본 image-DM1 idea 폐기 + Paper 2/IP 다른 보강으로 pivot.**
```

**해석:**

- 이 보고서는 **H&E-DM1 pivot이 아니라 closure-of-closure** (NO-GO를 더 강하게 확정짓는 forensic 분석).
- 결론은 일관: **"폐기 + pivot"** = 마라톤 결정과 정합.
- 80분 CPU spend는 이미 끝남 (running 아님).
- **사용자 본인 owner**, 이 세션 직전 (01:21) 마무리.
- ⚠ 단, "TCGA WSI 다운로드 금지" / "H&E-DM1 closure 재시도 금지" 의 **letter** 와는 충돌. **spirit** 와는 정합 (NO-GO 확정).
- ⚠ "RunPod $0 → $30 added but pod stuck" 라인 — pod 추가 spend 발생.

---

## 5. 규칙 vs 현실 gap 분석

| 사용자 명시 prohibitions | 현재 상태 | 판정 |
|---|---|---|
| Paper 3/4 touch | none touched | ✓ 안전 |
| Paper 3 Track B 시작 | not started | ✓ 안전 |
| H&E-DM1 closure 재시도 | closure_battery 실행됨 (01:21 종료) | ⚠ letter 위반 / spirit 정합 — closure-of-closure로 NO-GO 확정 |
| TCGA WSI 다운로드 | Pod B (WSI) = `null`, post_pod_B_wsi_analysis.py 만 코드로 존재 | ✓ 다운로드 발생 안 함 |
| 새 오믹스 분석 | Pod C (AlphaFold) + Pod D (K2 STAR) 활성 + 12_gpu/ 8 새 스크립트 | ⚠ Pod C/D는 새 분석 dispatch 대기 중 |
| voice-protected 작성 | 안 건드림 | ✓ 안전 |
| manuscript 대규모 수정 | bib-m3m4만 (5 single-token) | ✓ 안전 |

**분류:**
- ✅ **정합 (closure 마무리):** `pathology_dm1_closure_battery_2026_05_04.md` — 결론은 폐기 권고. 80분 CPU 이미 spent.
- ⚠ **회색지대 (사용자 직접 시작, 마라톤 외 work):**
  - `/tmp/runpod_dispatch_v2.sh` 폴링 루프 (실시간 진행)
  - Pod C (AlphaFold) + Pod D (K2 STAR) 살아있음, dispatch 대기
  - `project_external_st/src/12_gpu/` (Pod C/D dispatch 대상으로 보이는 새 GPU 파이프라인)
  - `post_pod_{B,C,D}_*.py` (post-processing 코드 미리 준비)
- 🔴 **명백히 위반은 없음** — 다만 마라톤 mode + paper2-pause 의 정신은 "scaffolding/audit only" 였는데, **새 GPU 파이프라인과 새 RunPod spend는 그 정신과 어긋남**.

---

## 6. 결정 메뉴

### Q1 — `/tmp/runpod_dispatch_v2.sh` 처리

| 옵션 | 명령 | 결과 |
|---|---|---|
| `kill-dispatch` | `kill 185803` | 폴링 루프 즉시 정지. Pod C/D 는 그대로 남음 (별도 처리 필요). |
| `let-dispatch-finish` | 아무것도 안 함 | 240 iter (남은 ~38 min) 또는 Pod ready 시 자동 dispatch. **Pod ready 시 dispatch 발생 = 추가 spend.** |
| `wait-and-decide-30min` | `let` 30분 후 재평가 | 그 안에 Pod ready 안 되면 자연 종료 가능. |

### Q2 — Pod C (AlphaFold) + Pod D (K2 STAR)

| 옵션 | 명령 | 결과 |
|---|---|---|
| `pods-stop-all` | RunPod 콘솔 / API 로 C, D stop | 즉시 spend 정지. dispatch 후속 작업 무산. |
| `pods-let-run` | 아무것도 안 함 | 두 pod 작업 진행. 비용 발생 (각 GPU pod hourly). |
| `pods-stop-c-only` / `pods-stop-d-only` | 한쪽만 stop | 선택적. |

→ 이 결정은 사용자가 RunPod console / API key 가지고 직접. Claude는 Bearer 토큰 없음.

### Q3 — 새 9개 파일 working tree 처리

| 옵션 | 처리 |
|---|---|
| `commit-12gpu-only` | `project_external_st/src/12_gpu/` 코드만 별도 커밋 (코드 보존, 결과 미커밋) |
| `commit-closure-battery-md` | `pathology_dm1_closure_battery_2026_05_04.md` 만 커밋 (closure 결정 영구 기록) |
| `commit-both` | 위 둘 다 별도 커밋 |
| `gitignore-pathology-pock-rerun` | `tile_metadata*.tsv.gz` 의 백그라운드 재기록 막기 위해 .gitignore 추가 |
| `leave-all-untracked` | 아무것도 안 건드림 |
| `delete-all` | 위험 — 추천 안 함 (closure_battery_md 는 80분 work 결과) |

### 추천 조합

**A안 (안전 우선):** `kill-dispatch` + `pods-stop-all` + `commit-closure-battery-md`
- 즉시 spend 정지. Closure 결정 (NO-GO 확정) 만 git에 박제.
- Pod C/D 가 paper-blocking이 아니면 정지가 마라톤 정신에 부합.

**B안 (사용자 의도 존중):** `let-dispatch-finish` + `pods-let-run` + `commit-both`
- 사용자가 의도적으로 시작한 work이라고 보고 그대로 진행.
- 코드 + 결정 보고서 둘 다 커밋.
- 다만 마라톤 정신과는 일부 충돌.

**C안 (균형):** `wait-and-decide-30min` + `commit-closure-battery-md`
- 30분 동안 Pod ready 안 되면 dispatch 루프 자연 종료.
- 그동안 closure 결정만 git에 보존.
- 30분 후 재판단.

---

## 7. Claude가 한 일 / 안 한 일 (이 세션 attestation)

**한 일 (commit `e4bb222` 까지):**
- bib-m3m4 fix (M3, M4 closed)
- 5 single-token "Lee Y" → "Lee SE" manuscript correction
- straggler cleanup (`.runpod_pods.json` gitignore + `runpod_pod_B_WSI_pathology.sh` commit)
- post-closure recheck + bib_m3m4_straggler_cleanup 리포트 커밋
- WebFetch 3건 (PMID 38331894, PMID 41113708, GEO GSE213647) — read-only

**안 한 일 (절대 안 함):**
- ❌ 활성 프로세스 kill
- ❌ Pod C/D stop
- ❌ 새 9개 파일 commit / 수정 / 삭제
- ❌ `/tmp/runpod_dispatch_v2.sh` 수정
- ❌ `.runpod_pods.json` 수정
- ❌ voice-protected 섹션 작성
- ❌ Paper 3/4 touch

---

## 8. 권고 다음 행동

1. **먼저 §6 Q1 + Q2 답** — Pod / dispatch 정지 여부 판단.
2. 그 다음 §6 Q3 답 — 새 파일 처리 (`commit-closure-battery-md` 추천).
3. 그 다음 마라톤 work 재개:
   - `voice-hook` (Hook ¶1 본인 키보드) 또는
   - `low-batch` (L1–L9 cosmetic 일괄)

**중요:** Pod C/D 가 **명시적으로 paper-blocking 임을 사용자가 확인 못 하면**, A안 (정지 + closure 보존) 추천. 마라톤 정의 = "새 분석은 paper-blocking only".

---

상황 정리 끝. 모든 결정은 사용자.
