# Session Full Report — 2026-05-07 (Disk · Pod · Hub Section 0)

**Branch:** `paper9-perturbation-extension-20260506`
**Date:** 2026-05-07 (Korea time)
**Mode:** Marathon (5/4-6/13), default scaffolding/infra
**Model/Author:** Claude Opus 4.7 (1M) · Seungho Cook

---

## Executive summary

이번 세션은 처음에 paper 1 사고 흐름 추가 요청에서 출발해, 디스크 ENOSPC, RunPod SSH wedge, /data 경로 정리, 12개 paper 사고 흐름 풀세트, hub badge까지 자율 사이클로 확장됐다. 끝나고 보면 4가지 큰 변화가 한 자리에 떨어졌다:

1. **VM 디스크 4.7G → 52G free** (97% → 58%) — caches + arcasHLA 32G 이전 + project/results·data bind mount
2. **RunPod L40S → A6000 교체** ($0.79/hr → $0.33/hr, **−58%**) — 기존 pod SSH wedged, L40S host 자리 없음 → 새 pod
3. **`CLAUDE.md` 신규 작성** + 메모리 `v19_runpod_a6000_2026_05_07` — 다음 세션 첫 페이지에서 disk layout + pod state + 5 rules 자동 인지
4. **12 paper 사고 흐름 (#s0)** + 9 hub badges — 모든 paper hub 페이지 첫 화면에서 "왜 이 결정?" 4가지 질문이 닫힘

---

## 타임라인 (phase by phase)

### Phase 1 — Paper 1 section 0 (시작)

> 사용자: "어떻게 8개 유전자가 되었는지... BRAF/RAS/TERT 왜 제외... TROP2 왜 갑자기... 표 그림 엄청난 설명으로 최초에 떄려야해"

- `rerun_v2.py:121-200` 직접 확인 — TIERA67 7 categories (67 entries) + `Driver_anchor` 12개 의도적 제외 → 55-gene clean pool → RandomForest top-8
- 살아남은 8개 모두 `TDS_core` 출신 (SLC5A5/TPO/TG/TSHR/PAX8/NKX2-1/FOXE1/DIO1)
- TROP2/TACSTD2는 별도 surface-marker exploratory track, GPL570 + spatial 외부 검증 실패로 demoted
- `paper1.html`에 `<section id="s0">` 추가 (158 lines, 5 tables + 2 figure refs + Q1-4 story-map)

### Phase 2 — Paper 2B / 3 / 4 (다른 패턴 적용)

> 사용자: "다른 패턴도 적용해줘"

- **Paper 2B (HT-overlap PTC mechanism)** — cohort 선택 (GSE286332 n=18), HT-isolated framing (Yu 2026-05-04 lock), 5 modules (HLA-II / IFN-γ / TLS / BCR / 8-gene), mediation analysis (Baron-Kenny + 5,000-iter bootstrap), TCGA n=500 generalization
- **Paper 3 (ICI vulnerability)** — vulnerability vs predictor (K1 kill), chmod 444 freeze, 4-tool HLA (cookHLA → Paper 4 reserve), DIAL audit 4-verdict, Track B-lite findings
- **Paper 4 GD HLA backlog** — renumber Paper 3→4 (Yu 2026-05-04), 4-gate gating, registry-only, disease firewall (shared allele ≠ shared disease)

### Phase 3 — 디스크 ENOSPC

> 사용자: "ENOSPC: no space left on device, write 디바이스 더 사야하나?"

- 현황: VM 123G 중 119G 차서 4.7G만 남음 (97%)
- 큰 범인 발견: `arcasHLA_GSE213647/*.alignment.p` 24G + `arcasHLA/*.alignment.p` 8.4G = **32G** (intermediate pickle, 최종 .genes.json은 git-tracked)

> 사용자: "일단 정리하고 조금 더 늘리긴 하자"

- 캐시류 5G 회수 (apt clean + journalctl --vacuum-time=7d + npm cache + ms-playwright + pip + torch)
- 결과: 4.7G → 9.8G free

### Phase 4 — 512G data disk 발견

> 사용자: "디스크 늘려"

- `lsblk` 확인 시 **`/dev/sdb` 512G xfs `/data`에 이미 mount됨** (FSAVAIL 438.2G, 14% 사용중)
- `/data/thca/` 이미 73G 사용중 (raw FASTQ, kallisto refs, etc.)
- **디스크 살 필요 전혀 없었음 — 이미 있던 게 그냥 놀고 있던 것**
- arcasHLA 32G를 `/data/thca/_repo_offload/`로 rsync 후 source 삭제 + symlink
- 결과: 4.7G → 42G free

### Phase 5 — RunPod SSH wedge

> 사용자: "럿팟도 부족해" / "넵 고고고 럿팟도 해결해줘"

- API/SSH 정보 발견: `~/.runpod/config.toml` + `~/.runpod/ssh/RunPod-Key-Go`
- Pod `thca-spark-dm-l40s-v3` (id `50uem6t82i0t9g`) SSH wedged (disk-full로 추정)
- `runpodctl exec`도 SSH 통해 작동 → 같이 hang
- Pod 상태 stop (EXITED, $0/hr)

> 사용자: "중요 데이터 다 /runpod-volume에 있음 — restart 고고"

- `podResume` 시도 → "There are not enough free GPUs on the host machine" (L40S 같은 host에 자리 없음)
- 60초 간격 retry loop 백그라운드 가동 (PID 2465461)
- 2시간/126번 retry 후에도 같은 메시지

> 사용자: "없어 딴거 가성비 좋은거로 할까 ㅠㅠ"

- L40S `$0.79` vs RTX A6000 `$0.33` (−58%) vs A40 `$0.35`
- 50GB volume = **pod-specific persistent** (NOT networkVolume) → terminate 시 손실
- 사용자 confirm: volume 손실 감수, RTX A6000 새 pod

> 사용자: "RTX A6000 새 pod — 경제적 충돌 회피 (volume 손실)"

- Loop kill, old pod terminate
- GraphQL `podFindAndDeployOnDemand` → 403 Forbidden
- `runpodctl create pod` 사용 (deprecated alias이지만 작동) → 새 pod `thca-spark-dm-a6000-v4` (id `uvp9i2r9s6l85y`)
- SSH `135.84.176.142:20788`, RTX A6000 검증, /runpod-volume 50GB 빈 상태
- 새 pod scaffolding: `/workspace/_tmp` + `/runpod-volume/{results,cache,refs,logs,scratch}` + persistent `TMPDIR=/workspace/_tmp` in `~/.bashrc`

### Phase 6 — `/data` 경로 정리 (symlink → bind mount)

> 사용자: "/data로 결과물 떨어뜨리는 경로 정리해줘" / "응 만들어줘"

- `project/results` (3.8G) + `project/data` (6.4G) rsync to `/data/thca/repo_results` + `/data/thca/repo_data`
- 처음 symlink로 했더니 **git이 3595개 tracked file을 "deleted"로 인식** — symlink는 git 호환 안 됨
- 즉시 **bind mount로 전환** (`mount --bind /data/thca/repo_* /home/seungho/.../project/*`)
- `/etc/fstab` 영구화 (`bind,nofail` 옵션)
- `findmnt` 검증: `/dev/sdb[/thca/repo_results]` xfs
- **부수효과 발견:** rsync source/destination 둘 다에 1782개 tracked file이 missing — **pre-existing** stale tracked 상태 (이번 마이그레이션 무관)

### Phase 7 — Documentation

- **`/home/seungho/personal/THCA_data_analysis/CLAUDE.md`** 신규 작성:
  - Disk layout (root vs /data) + bind mount 설명
  - 5 rules (>100MB write, new dir 만들 때, /tmp, df 75%, /data mount 검증)
  - RunPod conventions (current pod info, GraphQL 403 → runpodctl, volume types)
  - Marathon mode notes (voice-protected sections)
- **메모리 `v19_runpod_a6000_2026_05_07.md`** 신규 — pod 교체 사실 + GraphQL 403 → runpodctl 사용 패턴
- **`MEMORY.md`** index 업데이트 (33줄 → 34줄)

### Phase 8 — Paper 5 / 6 / 7 / 9 / 10 / 11 / 12 / 2A section 0

> 사용자: "더더더더" / "6 7 9 10 12 다 박아줘" / "다음 해요"

순차 추가 (각각 100-160 lines, 4-5 tables + Q1-4 + closing callout):

| Paper | Hero theme | Section 0 angle |
|---|---|---|
| **Paper 5** (npj 8-gene RAI biomarker) | crimson #7B1F2A | venue (npj Prec Onc) + framing (RAI biomarker vs DM1 dark matter 시계열 관계) + TERT recovery (cBioPortal `thca_tcga_pub` n=36 → Cox HR=4.33) + paused 이유 (옵션 보존) |
| **Paper 6** (v15 NeurIPS DIAL · v52 LODO) | green #5C7B6E | NeurIPS submission 정체 + v52 self-audit (DIAL 0.494 → 0.000 collapse) + algorithm vs application 분리 + downstream halt rule (4 paper에서 cite 차단) |
| **Paper 7** (v18 agentic_research framework) | teal #3F7A8A | framework + Python pkg 결정 + 5 patterns 추출 (case study에서) + 6 components codify (test pass) + 22-slide Korean talk portable |
| **Paper 9** (Metabolic Vulnerability Map) | teal | "map" framing (SL paper와 분리) + Paper 1 anchor + 9/10/12 trajectory 분리 + research-use boundary (clinical actionability 금지) |
| **Paper 10** (DM1 SL Atlas) | teal | Paper 9의 visualization-first companion + 17 figs + 9 sortable tables + DepMap deferred 4 reason + live-inspection UX |
| **Paper 11** (Pan-Cancer DM1 Transfer Atlas) | teal | pancan transfer motivation + substrate (TCGA pancan 11,069 × 33 lineage) + per-lineage anchor activity caveat + Paper 9/10/12 boundary 공유 |
| **Paper 12** (Co-expression Network) | teal | 36 candidate 입력 lineage (Paper 9/10) + network 추출 절차 (n=572 → 137 edges → 7 modules → 5 hubs) + state-conditional Δr 의미 + co-expression vs perturbation boundary |
| **Paper 2A** (H&E NO-GO) | green #3C6B4F | hypothesis motivation + NO-GO 판정 강도 (random보다 낮음) + publish 가치 (actionable upper bound) + 5/5 re-entry binding |

**Index hub badges 9개** (paper2a + paper6/7는 hub 카드 없음):
- paper1#s0 / paper2b#s0 / paper3#s0 / paper4_gd_hla#s0 / paper5#s0 / paper9#s0 / paper10_atlas#s0 / paper11_pancancer#s0 / paper12_network#s0

**HTML balance 검증**: 12/12 모두 unclosed=0, extra-close=0

---

## 확보한 artifacts (전체 list)

### 새로 생성된 파일
- `CLAUDE.md` (project root) — disk layout + RunPod conventions + 5 rules + marathon notes
- `~/.claude/projects/.../memory/v19_runpod_a6000_2026_05_07.md` — pod 교체 메모리
- `SESSION_FULL_2026_05_07_DISK_POD_HUB.md` (this file)

### 수정된 paper hub 페이지 (12개)
- `paper1.html` (section 0 — 유전자 사고 흐름)
- `paper2a.html` (section 0 — NO-GO 사고 흐름)
- `paper2b.html` (section 0 — 방법론 사고 흐름)
- `paper3.html` (section 0 — 전략 사고 흐름)
- `paper4_gd_hla.html` (section 0 — backlog 사고 흐름)
- `paper5.html` (section 0 — npj 사고 흐름)
- `paper6.html` (section 0 — forensic 사고 흐름)
- `paper7.html` (section 0 — framework 사고 흐름)
- `paper9.html` (section 0 — map 사고 흐름)
- `paper10_atlas.html` (section 0 — atlas 사고 흐름)
- `paper11_pancancer.html` (section 0 — transfer 사고 흐름)
- `paper12_network.html` (section 0 — network 사고 흐름)
- `index.html` (9 hub badges 추가)

### 시스템 변경
- `/etc/fstab` (sudo) — bind mount 영구화 2줄 추가
- `~/.bashrc` (pod) — TMPDIR=/workspace/_tmp persistent

### 디스크 이전 (~ 42G total to /data)
- `arcasHLA/` 8.4G → `/data/thca/_repo_offload/arcasHLA` + symlink
- `arcasHLA_GSE213647/` 24G → `/data/thca/_repo_offload/arcasHLA_GSE213647` + symlink
- `project/results/` 3.8G → `/data/thca/repo_results/` + bind mount
- `project/data/` 6.4G → `/data/thca/repo_data/` + bind mount

### RunPod 변경
- Pod `thca-spark-dm-l40s-v3` (id `50uem6t82i0t9g`) **terminated** — 50GB volume 손실 감수
- Pod `thca-spark-dm-a6000-v4` (id `uvp9i2r9s6l85y`) 신규 — RTX A6000 $0.33/hr, SSH `135.84.176.142:20788`
- Container disk 100GB + /runpod-volume 50GB (빈 상태)

---

## 현 상태 (post-session snapshot)

### Azure VM
| metric | before | after |
|---|---|---|
| Root disk used | 119G (97%) | 72G (58%) |
| Root disk free | 4.7G | 52G |
| `/data` (512G) used | 73G (15%) | 113G (22%) |
| `/data` free | 438G | 399G |
| Bind mounts | none | results + data |

### RunPod
| metric | value |
|---|---|
| Active pod | `thca-spark-dm-a6000-v4` (uvp9i2r9s6l85y) |
| GPU | RTX A6000 48GB |
| Cost | $0.33/hr (community cloud) |
| SSH | `135.84.176.142:20788` |
| /runpod-volume | 50GB **빈 상태** (이전 50GB 손실) |
| /workspace/_tmp | 생성됨 (TMPDIR 영구화) |
| /runpod-volume layout | results/ cache/ refs/ logs/ scratch/ |
| Bioinfo tools | 미설치 (필요 시 apt-get) |

### Git
| state | count |
|---|---|
| Modified (이번 세션 변경) | 14 (12 paper + index.html + CLAUDE.md) |
| Modified (auto-injected data dictionary) | ~50 (paper hub HTML 자동 주입) |
| Deleted (pre-existing stale) | 1782 |
| Untracked | 145 (paper3 ICI 결과 + bak + figures + TSV + CLAUDE.md + this MD) |

---

## 마지막 — 판단 필요한 것 (User decisions)

이 세션이 자율로 처리한 것은 모두 reversible / non-destructive scaffolding이지만, 아래 항목들은 user 결정이 필요하다.

### D1. Git 1782 deleted entries 처리 방향
- 상황: rsync 마이그레이션 무관하게 **session 시작 전부터** 1782개 tracked file이 filesystem에 없는 stale 상태. `/data/thca/repo_results`에도 없음.
- 옵션 A: `git rm --cached <files>` → 정식 deletion으로 commit (가장 빠름, 정직하게 정리)
- 옵션 B: 일부는 git history에서 복원, 일부는 cached 제거 (선별 검토 필요, 시간 소요)
- 옵션 C: 그대로 두고 다음 commit에서 자연스럽게 정리 (덜 깨끗하지만 disruption 적음)
- **추천: A** — 어차피 disk에 없는 파일들이고, 필요하면 git history에서 retrieve 가능. 단 commit message에 "session 시작 전 stale tracked deletion 일괄 정리" 명시.

### D2. RunPod /runpod-volume 재구축
- 상황: 이전 pod의 50GB volume 데이터 손실. 새 pod의 /runpod-volume은 빈 상태.
- 옵션 A: 필요한 ref/cache 데이터를 local /data/thca에서 새 pod로 scp upload (재현 가능한 정확한 list 필요)
- 옵션 B: 당장 안 쓸 거면 그대로 두고 필요 시 on-demand
- **추천: B** (marathon mode에서 manuscript writing 우선이라 pod GPU 당장 안 씀). 6/13 Paper 1 ship 후 Track B 시작 시 다시 셋업.

### D3. arcasHLA `.alignment.p` 32G 운명
- 상황: intermediate pickle, regeneratable. 현재 `/data/thca/_repo_offload/`에 있음 (root disk 영향 0). 최종 outputs (.genes.json/.tsv)은 별개 파일로 살아있음.
- 옵션 A: 그대로 보존 — re-analysis 시 재사용 (디스크 32G 사용 유지)
- 옵션 B: 즉시 삭제 — `/data/thca/repo_results/v17_korean/` 안의 .genes.json만 남기고 32G 회수 → 주1회 재분석 비용 vs 32G 절약
- **추천: A** (현재 /data 358G 여유, 압박 없음). 6/13 ship 후 정리 검토.

### D4. CLAUDE.md 분량/scope
- 상황: 신규 작성 (5 sections — disk layout / rules / RunPod / marathon). 처음이므로 reference 분량이 어느 정도가 적절한지 user 판단.
- 옵션 A: 현재 분량 유지 — 다음 Claude 세션이 disk + pod context까지 한 번에 인지
- 옵션 B: Disk layout만 남기고 RunPod/Marathon은 메모리로 분리 — CLAUDE.md를 더 lean하게
- **추천: A** (CLAUDE.md는 자동 로드되니 가장 critical한 운영 지식을 모아두는 게 효율)

### D5. Section 0 12-paper 일관성 유지
- 상황: 12개 paper 모두 #s0 박혀있음. 향후 paper 변경 시 (e.g., Paper 1 manuscript revision) 사고 흐름도 같이 update해야 일관성 유지.
- 옵션 A: section 0를 paper hub의 영구 convention으로 박기 (template 포함)
- 옵션 B: 일회성으로 두고 자연스럽게 stale되도록 허용 (낮은 priority)
- 옵션 C: section 0 자동생성 script 작성 (data dict v2처럼 auto-injection)
- **추천: B** (marathon 끝까지 manuscript v8 안정화에 집중, section 0은 frozen). Marathon 후 C 검토.

### D6. Hub badges 추가 여부
- 상황: 12개 paper 중 9개에 hub 카드 + badge. paper2a (NO-GO) / paper6 / paper7는 hub 카드 자체가 없어 footer text link로만 접근.
- 옵션 A: paper2a/6/7도 hub 카드 추가 — 가시성 ↑ but 카드 수 ↑ (busy)
- 옵션 B: 그대로 — sub-link / footer text link 유지
- 옵션 C: paper2a만 추가 (NO-GO finding은 reviewer 자주 묻는 영역)
- **추천: B** (현재 hub가 9 카드면 시각적으로 적절. NO-GO/forensic/framework는 cross-reference로 도달).

### D7. Untracked 145개 파일 처리
- 큰 카테고리:
  - `index.html.bak-immune-readiness-2026-05-06` — 백업 (delete 가능)
  - `paper3_ici_track_b_lite/figures_png/fig{7,8,9}_*.png + fig_montage_8panel.png` — Track B-lite 결과 figures
  - `paper3_ici_track_b_lite/*.tsv + .json` — Track B-lite 분석 산물
  - `scores_per_cohort/*.tsv` — module score outputs
  - `STATUS_immune_readiness_section_2026-05-06.md` — status note
  - `fix_rai_split.py + track_b_lite_phase4_continue.py` — Track B-lite scripts
  - `CLAUDE.md` (이번 세션 신규) + `SESSION_FULL_2026_05_07_DISK_POD_HUB.md` (이 파일)
- 옵션 A: Track B-lite results은 staging area에 commit, bak/.md는 .gitignore
- 옵션 B: Marathon 안에서는 안 commit, 일괄 처리는 6/13 ship 후
- **추천: A 부분 적용** — bak 파일 삭제 + CLAUDE.md/이 MD는 commit 후보 + Track B-lite는 user 검토.

### D8. 다음 자연스러운 작업
세션 마지막 시점에서 가능한 진행:
1. **Paper 1 manuscript v8 voice-protected 작성** (Hook/Aim/Disc 3.1/Limitations/Cover Para 1/Q9) — **user 키보드만**, Claude 진입 금지 (per `v17_sprint_vs_marathon_violation`)
2. **Paper 2B Pillar II/III closure** (Yu input pending)
3. **Paper 1 bioRxiv 제출 준비** (target 6/13)
4. **D1-D7 결정 후 git commit** — 단일 commit "session: 2026-05-07 disk migration + 12-paper section 0 + CLAUDE.md"

---

## Reference — 변경된 파일 전체 목록

### 신규
```
/home/seungho/personal/THCA_data_analysis/CLAUDE.md
/home/seungho/personal/THCA_data_analysis/SESSION_FULL_2026_05_07_DISK_POD_HUB.md
~/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/v19_runpod_a6000_2026_05_07.md
```

### 직접 수정 (이번 세션)
```
project/papers_hub_2026_05_04/paper1.html         (+158 lines section 0)
project/papers_hub_2026_05_04/paper2a.html        (+section 0)
project/papers_hub_2026_05_04/paper2b.html        (+section 0)
project/papers_hub_2026_05_04/paper3.html         (+section 0)
project/papers_hub_2026_05_04/paper4_gd_hla.html  (+section 0)
project/papers_hub_2026_05_04/paper5.html         (+section 0)
project/papers_hub_2026_05_04/paper6.html         (+section 0)
project/papers_hub_2026_05_04/paper7.html         (+section 0)
project/papers_hub_2026_05_04/paper9.html         (+section 0)
project/papers_hub_2026_05_04/paper10_atlas.html  (+section 0)
project/papers_hub_2026_05_04/paper11_pancancer.html (+section 0)
project/papers_hub_2026_05_04/paper12_network.html (+section 0)
project/papers_hub_2026_05_04/index.html          (+9 NEW badges)
~/.claude/projects/-home-seungho-personal-THCA-data-analysis/memory/MEMORY.md (+1 entry)
```

### 시스템 (sudo)
```
/etc/fstab — 2 bind mount lines added
~/.bashrc (pod) — TMPDIR exports
```

### 이전 (rsync + symlink/bind mount)
```
project/results        → /data/thca/repo_results       (bind mount)
project/data           → /data/thca/repo_data          (bind mount)
project/results/v17_korean/arcasHLA           → /data/thca/_repo_offload/arcasHLA           (symlink)
project/results/v17_korean/arcasHLA_GSE213647 → /data/thca/_repo_offload/arcasHLA_GSE213647 (symlink)
```

### Pod
```
TERMINATED: thca-spark-dm-l40s-v3 (50uem6t82i0t9g) — 50GB volume lost
NEW:        thca-spark-dm-a6000-v4 (uvp9i2r9s6l85y), RTX A6000, $0.33/hr
SSH:        135.84.176.142:20788
```

---

## Sanity check commands (다음 세션이 검증용)

```bash
# 디스크 상태
df -h / /data

# Bind mounts 살아있나
findmnt /home/seungho/personal/THCA_data_analysis/project/results
findmnt /home/seungho/personal/THCA_data_analysis/project/data

# Section 0 모두 박혔나
cd /home/seungho/personal/THCA_data_analysis/project/papers_hub_2026_05_04
for f in paper*.html; do grep -q 'id="s0"' "$f" 2>/dev/null && echo "✓ $f"; done | wc -l   # → 12

# Pod 상태
curl -s -X POST https://api.runpod.io/graphql \
  -H "Authorization: Bearer $(grep apikey ~/.runpod/config.toml | cut -d'"' -f2)" \
  -H "Content-Type: application/json" \
  -d '{"query":"query { myself { pods { id name desiredStatus runtime { ports { ip publicPort privatePort type } } } } }"}'

# Pod SSH 검증
ssh -i ~/.runpod/ssh/RunPod-Key-Go -p 20788 root@135.84.176.142 "df -h && nvidia-smi -L"
```

---

*Generated 2026-05-07 by Claude Opus 4.7 (1M) under marathon mode.*
*Original session: `de1b36b4-2b56-46bf-8948-06af3822451d`.*
