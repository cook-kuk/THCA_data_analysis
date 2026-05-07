# RunPod Pod C/D 현재 상태 + "오래 걸리더라도 철저하게 고고" 실행계획

**Author:** Seungho Cook · **Date:** 2026-05-04 (post-closure marathon mode)
**작성 시점 컨텍스트:** PAPER1_DM1_FULL × Closure Battery conflict audit 완료 직후, 사용자 "오래 걸리더라도 철저하게 고고" 신호.
**모드:** 백그라운드 GPU 작업 + 결과 회수 + 후처리. 새 dispatch 없음.

---

## 0. 한 줄 요약

**Pod C (AlphaFold) + Pod D (K2 STAR) 두 GPU pod 현재 RunPod 상에서 running 중.** Pod B (WSI) 는 closure battery NO-GO 후 의도적 kill 됨. "분석 안겹치게 고고" 의 가장 자연스러운 해석 = **새 dispatch 없이, 두 pod 결과 기다려서 SCP + 후처리 fully execute**. 단 Pod C 는 Paper 4 (Korean GD HLA backlog) territory 에 닿는 결과물을 만들어내므로 사용자 confirmation 필요 (§3.3).

---

## 1. Live state (확인 시각 ~ 02:?? KST 기준)

### 1.1 Pod 상태

| Pod | RunPod ID | SSH | Dispatched | 작업 | 예상 runtime | Auto-shutdown | 상태 |
|---|---|---|---|---|---|---|---|
| **D** (K2 STAR) | `xgk5y4oqwvo7te` | `157.157.221.29:21123` (dispatch.log 기준; `.runpod_ssh.json` 의 `209.170.80.132:12043` 은 stale) | ✅ 01:27:27 | K2 PRJEB11591 30-sample PoC STAR re-quant + featureCounts 8-gene | **~30 hr on 16-core / 8-12 hr on 32-core** (Pod 사양에 따라) | 24 hr forced | **running** |
| **C** (AlphaFold) | (pods.json 에 `null` 로 기록됐으나 dispatch log 는 dispatched 표시) | `216.81.151.3:10960` | ✅ 01:26:53 | LocalColabFold 4 structures: HLA-DPB1\*05:01 × {TSHR p1, Tg p1, TPO p1} + HLA-B\*46:01 × TSHR p1 | **8-12 hr** | 24 hr forced | **running** |
| **B** (WSI) | `null` | — | — | — | — | — | **intentionally killed** (closure NO-GO 후) |

### 1.2 Dispatch 인프라

- `/tmp/runpod_dispatch_v2.sh` 폴링 루프 — 이미 종료 (`ALL_3_DISPATCHED` 도 아님; D + C 두 pod 만 dispatched, B null 이라 skip)
- 새 polling/dispatch 프로세스 없음 (`ps aux | grep runpod` 비어있음)
- `/tmp/dispatch.log` 마지막 줄: `[01:27:38] iter=6, dispatched=2/3 (D=1 C=1 B=0)`
- RunPod API key: `~/.runpod/config` (Bearer 사용 가능)
- SSH key: `~/.ssh/id_ed25519` (dispatch 가 `accept-new` 로 known_hosts 등록 완료)

### 1.3 본 VM 의 GPU 상태

- `nvidia-smi: command not found` — main VM (Standard_D8as_v5) 에 GPU **없음**
- 모든 GPU 작업은 RunPod 측에서만 진행됨 (= 사용자 의도와 일치 = "GPU 있어 분석 안겹치게")

---

## 2. 사용자 신호 해석

직전 turn 메시지 흐름:

1. **First turn (audit 요청):** "지금 바로 RunPod G1/G2/G3 실행하지 말고, 먼저 conflict audit 을 수행한다" + 절대 금지 8개 (RunPod / G1/G2/G3 / TCGA WSI / H&E-DM1 / Foundation model / Paper 3/4 / 새 분석 / manuscript prose 대규모 수정)
2. **Audit 완료** → image-DM1 dropped, RunPod G1/G2/G3 blocked, 5-condition re-entry gate 0/5
3. **Second turn:** "고고 gpu 있어 분석 안겹치게"
4. **Third turn:** "오래 걸리더라도 철저하게 고고"

**해석 (best read):**

- "GPU 있어" = Pod C/D 가 이미 RunPod 상에서 LIVE — 새 GPU 확보가 아니라 **이미 running 중인 pod 자원을 가리킴**
- "분석 안겹치게" = (a) 새 pod dispatch 안 하고 (b) Pod C/D 와 동일 작업 중복 안 하고 (c) image-DM1 (closure NO-GO 된 line) 재시도 안 함
- "오래 걸리더라도" = Pod C 8-12hr + Pod D 8-30hr 가 끝날 때까지 대기 OK
- "철저하게" = SCP + 후처리 + 결과 보고서 까지 fully execute, 중간에 짜르지 말고

**해석에 따른 "고고"의 scope:**

- ✅ Pod C/D 결과 기다림 (새 dispatch 0)
- ✅ SCP back to main VM
- ✅ `post_pod_D_meta3cohort.py` 실행 → Paper 1/2 K2 calibration
- ⚠ `post_pod_C_alphafold_analysis.py` 실행 → **Paper 4 (Korean GD HLA backlog) territory**, §3.3 별도 confirmation 필요
- ❌ Pod B 재시작 (closure NO-GO)
- ❌ G1/G2/G3 image-DM1 retry (closure NO-GO + audit lock)
- ❌ Foundation model UNI/CONCH (5-condition gate 0/5)

---

## 3. 작업별 scope 판정

### 3.1 Pod D (K2 STAR re-quantification) — ✅ Paper 1/2 직접 관련, GO

**스크립트:** `project/notebooks_or_scripts/runpod_pod_D_K2_STAR.sh`

**작업:**
- K2 PRJEB11591 (Korean PTC, n=260) 중 **30 sample PoC** ENA 5M-read partial download
- GENCODE v44 primary assembly + basic annotation
- STAR sparse index → STAR alignment → featureCounts on **8-gene panel** (`SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1`)

**Scope:** Paper 1 (8-gene DM1 dark-matter) 및 Paper 2 (Korean PTC vs HT-overlap) 의 K2 cohort 처리 quality 직접 영향. Memory `v17_korean_K2_calibration` 의 "K2 TPM 10-100× inflated, TCGA-trained classifier wrong-direction" 문제가 STAR-기반 raw count → CPM normalize 로 해결되는지 검증.

**Output (post_pod_D 가 만들 파일):**
```
project/results/d7p3_k2_starred/
  k2_counts.tsv                       (raw counts)
  k2_starred_8gene_logCPM.tsv         (log2 CPM)
  k2_starred_8gene_centered.tsv       (within-sample-centered profile)
  3cohort_meta_updated.json
```

**판정:** **GO**. Paper 1/2 directly relevant. 사용자 prohibition 리스트의 "Paper 3/4 touch" 룰 영향 없음.

### 3.2 Pod B (WSI pathology) — ❌ NO-GO (이미 kill 됨)

**상태:** `.runpod_pods.json` 에서 `null`. Pod B 는 Phase A NO-GO + closure battery 결정 후 의도적으로 kill 됨.

`post_pod_B_wsi_analysis.py` 코드는 디스크에 존재하지만 **실행 트리거 없음** (입력 파일이 영원히 안 옴 → 스크립트 line 21-24 의 fallback 으로 종료).

**판정:** **NO-GO 유지**. closure battery 결정 + audit §8.2 lock 그대로.

### 3.3 Pod C (AlphaFold HLA-peptide) — ⚠ Paper 4 territory, **사용자 confirmation 필요**

**스크립트:** `project/notebooks_or_scripts/runpod_pod_C_AlphaFold.sh`

**작업:**
- LocalColabFold 4 structure 예측:
  1. HLA-DPB1\*05:01 + TSHR p1 (residues 252-261, autoreactive 9-mer)
  2. HLA-DPB1\*05:01 + Thyroglobulin p1 (residues 2549-2570)
  3. HLA-DPB1\*05:01 + TPO p1 (residues 535-552)
  4. HLA-B\*46:01 + TSHR p1 (Asian-specific HLA-I)

**post_pod_C 가 만드는 narrative (스크립트 line 64-90 발췌):**
- *"asian_pacific_islander HLA reference"*
- *"Korean Asian-specific Graves' risk allele"* — Graves' disease (GD) 는 **Paper 4** scope
- *"Pillar 5 Supplement"* — Pillar 5 는 어느 Paper 의 Pillar 인지 코드만으론 모호. Paper 2 가 Pillar I/II/III (HT-overlap), Paper 4 도 자체 pillar 구조 가능. 코드 내 "Pillar 5" 는 cross-paper extension 표시 가능.
- TSHR / Tg / TPO autoantigen + DPB1\*05:01 = Korean autoimmune thyroid panel — **Paper 4 (GD HLA Pan-Asian) 의 핵심 mechanism leg**

**Scope 충돌 분석:**

| 룰 / lock | Pod C 후처리 영향 |
|---|---|
| 사용자 first-turn 절대 금지 "Paper 3/4 touch" | ⚠ post_pod_C 결과물 (`pillar5_supplement.md`) = Paper 4 territory output |
| Memory `v19_paper4_GD_backlog` "Track B 시작 전까지 Paper 4 frozen" | ⚠ 직접 위반 |
| Memory `v18_paper2_HT_isolated` "GD/Chu 2018 forest → Paper 4" | TSHR/Tg/TPO autoantigen + Korean HLA = GD core, Paper 4 영역 |
| Audit §8.2 "Paper 3/4 not touched" | ⚠ post_pod_C 가 결과 파일 만들면 위반 |
| 단, Pod C 자체는 RunPod 상에서 이미 사용자가 직접 dispatch 한 작업 | "dispatch" ≠ "touch local files"; pod 가 알아서 돌고 있는 건 ok |

**3가지 옵션 (사용자 confirmation 항목):**

| 옵션 | 의미 | 결과 파일 만듦 | Paper 4 룰 영향 |
|---|---|---|---|
| **C1: Pod C 결과 회수만, 후처리 보류** | SCP `SUMMARY.json` + per-structure PDB 만 회수, `pillar5_supplement.md` narrative 생성 안 함, raw 파일만 보존 | raw 파일만 (`structure_summary.tsv` 가능, `pillar5_supplement.md` 보류) | 회색지대 — pod 가 만든 raw output 보존은 "touch" 보다 약함 |
| **C2: Pod C 후처리 fully execute** | post_pod_C 그대로 실행, Pillar 5 supplement 생성 | full | ⚠ Paper 4 룰 위반 — 명시적 lift 필요 |
| **C3: Pod C 결과 폐기** | RunPod 자동 24hr shutdown 까지 두고 회수 안 함 | 0 | clean — 룰 그대로 유지 |

**기본값 (사용자 명령 없으면):** **C1** — raw 결과만 보존 (8-12hr 컴퓨트 + RunPod 비용 이미 지출되었으므로 폐기 비합리), 후처리/narrative 는 Paper 4 unfreeze 시점까지 보류.

**사용자 명시 확인 필요:** "Pod C 결과는 C1 / C2 / C3 중 어느 것으로?" 답변 안 오면 C1 진행.

---

## 4. 실행 계획 (오래 걸리더라도 철저하게)

### Step 1 — Pod D 모니터링 (지금 시작)

- 폴링 인터벌: 30분 (Pod D 8-30hr runtime, 너무 자주 안 봐도 됨)
- 완료 마커: `/workspace/k2_star/results/k2_counts.tsv` 존재 + `/workspace/run.log` 끝에 `=== ALL DONE ===` (D 스크립트는 `=== Phase 7 complete ===` 등 stage marker 가능; tail 로 확인)
- SSH/SCP target: `root@157.157.221.29:21123` (dispatch log 기준)

### Step 2 — Pod C 모니터링 (병렬)

- 폴링 인터벌: 30분
- 완료 마커: `/workspace/alphafold/outputs/SUMMARY.json` 존재 + `=== ALL DONE ===` (run.sh tail)
- SSH/SCP target: `root@216.81.151.3:10960`
- 회수 후 처리는 §3.3 옵션 결정에 따름

### Step 3 — Pod D 완료 시점

1. SCP `k2_counts.tsv` → `project/results/d7p3_k2_starred/`
2. `python3 post_pod_D_meta3cohort.py` 실행
3. Output 검증:
   - 8-gene log CPM range 합리적인가
   - Within-sample-centered profile 분포
   - TCGA-trained classifier 적용 시 K2 sample distribution 이 memory `v17_korean_K2_calibration` 의 "TPM 10-100× inflated, classifier wrong-direction" 문제 해결되는가
4. 결과 보고서 `project/reports/2026_05_04_pod_D_k2_star_postprocess.md` 작성
5. Paper 1/2 manuscript scaffold (`project/results/p2_pillar1_forest/cohort_assembly.tsv`) 에 K2 cohort 처리 update 권고 (실제 manuscript prose 수정 안 함; 권고만)

### Step 4 — Pod C 완료 시점 (옵션 C1 default)

1. SCP `SUMMARY.json` + `*_scores_rank_001*.json` 4 structures × 3 models = 12 file → `project/results/c_alphafold/`
2. `structure_summary.tsv` 생성 (post_pod_C line 22-26 부분만; narrative 부분은 skip)
3. raw plDDT / iPTM 통계만 단순 dump
4. `pillar5_supplement.md` narrative 생성 보류 (사용자 명시적 GO 시 진행)
5. 단순 raw report `project/reports/2026_05_04_pod_C_alphafold_raw.md` 작성: "raw scores landed; Paper 4 narrative deferred per backlog freeze"

### Step 5 — 모든 pod 완료 후 (수 시간 ~ 1일)

1. Cumulative report `project/reports/2026_05_04_pod_cd_completion_report.md`
2. Marathon 복귀 — manuscript writing (voice-protected 사용자 키보드)

---

## 5. Hard locks (이번 세션 안 건드림)

audit `2026_05_04_paper1_dm1_full_conflict_audit.md` §8.2 그대로:

- ❌ H&E → DM1 / image-DM1 retry (any tile size, any embedder)
- ❌ TCGA WSI Phase C download
- ❌ Image-DM1 IP claim
- ❌ RunPod G1/G2/G3 image stack
- ❌ Foundation model UNI/CONCH/Virchow2 retry
- ❌ Wet-lab roadmap
- ❌ Paper 3 Track A/B (chmod 444 freeze 그대로)
- ❌ 새 RunPod pod dispatch (B 부활 포함)
- ⚠ Paper 4 (Pod C 후처리만 §3.3 옵션 따라 — narrative 생성은 보류 default)

---

## 6. 사용자 confirmation 요청 1개

**Q1 — Pod C 후처리 어떤 옵션으로?**

- **C1 (default, 권장):** raw `SUMMARY.json` + plDDT/iPTM 통계만 회수, `pillar5_supplement.md` narrative 생성 보류 → Paper 4 freeze 유지
- **C2:** Pod C 후처리 fully execute (Pillar 5 narrative 생성 포함) → Paper 4 룰 명시적 lift 필요
- **C3:** Pod C 결과 폐기, 24hr auto-shutdown 까지 두고 회수 안 함 (8-12hr 컴퓨트 매몰비)

답변 안 오면 default **C1** 진행.

기타 모든 액션은 §4 계획대로 자동 진행. Pod 완료까지 8-30시간 대기 동안 main VM 에서는 manuscript bib-m3m4 / scaffold-only 작업 OR 사용자 voice-protected work 가능.

---

## 7. 첨부 파일 + 경로

**Live state inputs:**
```
/tmp/dispatch.log                                    dispatch trace
/tmp/runpod_dispatch_v2.sh                           dispatch script (이미 종료)
project/.runpod_pods.json                            {"D": "xgk5y4oqwvo7te", "C": null, "B": null}
project/.runpod_ssh.json                             stale C/D ports — 사용 시 검증 필요
~/.runpod/config                                     Bearer key
~/.ssh/id_ed25519                                    SSH key
~/.ssh/known_hosts                                   pod hosts already registered (line 01:26-01:27)
```

**Pod scripts (이미 dispatched, pod 측 /workspace/run.sh 로 업로드됨):**
```
project/notebooks_or_scripts/runpod_pod_C_AlphaFold.sh
project/notebooks_or_scripts/runpod_pod_D_K2_STAR.sh
project/notebooks_or_scripts/runpod_pod_B_WSI_pathology.sh   (not dispatched)
```

**Post-pod scripts (main VM 에서 실행):**
```
project/notebooks_or_scripts/post_pod_C_alphafold_analysis.py
project/notebooks_or_scripts/post_pod_D_meta3cohort.py
project/notebooks_or_scripts/post_pod_B_wsi_analysis.py       (won't run — pod B 죽음)
```

**Reference reports:**
```
CLOSURE_BATTERY_2026_05_04.md
PHASE_A_NOGO_REPORT_2026_05_04.md
SITUATION_BACKGROUND_ACTIVITY_2026_05_04.md
project/reports/2026_05_04_image_dm1_final_nogo_decision.md
project/reports/2026_05_04_paper1_dm1_full_conflict_audit.md     ← 직전 audit
project/reports/2026_05_04_post_closure_recheck.md
```

---

*VS Code 에서 보고 결정. C1 / C2 / C3 답주시면 그대로 §4 계획 실행. 답 없으면 C1 default. Pod 완료 모니터 중 main VM 작업 병행 가능 (bib-m3m4, low-batch, voice-protected user work).*
