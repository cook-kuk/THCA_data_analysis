# Situation Overview — 2026-05-07 (전체 상황 처음부터)

작성: 2026-05-07
저자명: Seungho Cook (논문)
지도교수: Yu Hyeong-won
모드: **Marathon writing (2026-05-04 → 2026-06-13, 6주 manuscript writing window)**

---

## 0. 한 문장 요약

오늘 (2026-05-07) 하루 동안 spatial transcriptomics 분석 패키지를 **6 figures (v1) → 100 figures (v12)** 까지 확장했고, 그 과정에서 advisor Q1 (TROP2가 왜?) / Q2 (공간 전사체가 왜 필요?) 답을 직접 시각화했으며, 추가로 **신규 finding 4가지** (분화 공간붕괴 / 24-axis stage gradient / ATC-2 super-organizer / LPTC-3 두 번째 outlier) 가 signature-level exploratory 로 잡혔다. Yu 분리벽 / voice-protected section / Paper 1 main story 모두 무손상.

---

## 1. 프로젝트 구조 (재확인)

### 1.1 Paper numbering (canonical, 2026-05-04 fix)

| Paper | 주제 | 상태 (2026-05-07) |
|---|---|---|
| **Paper 1** | DM1 molecular dark matter (8-gene + DM1 axis) | npj 1차 ship-ready (제출 직전, 4 outreach 미발송, user click 대기) |
| **Paper 2** | H&E → DM1 / pathology projection / TCGA validation (HT-overlap PTC ONLY; Yu 2026-05-04 결정) | Pillar I 본문 작성 + Pillar II HT TLS spatial 직접 확인 추가 |
| **Paper 3** | ICI vulnerability dark thyroid cancer | Track A bundle FROZEN (chmod 444). Track B BLOCKED |
| **Paper 4** | Korean GD HLA / Pan-Asian (renumbered from old "Paper 3 GD" → backlog) | 4/4 gating, scope rules 적용 |
| **Paper 5+** | hub 추가 (paper 6/7/9/10/11/12 등 — papers_hub master_view) | infra-level tracking only |

### 1.2 분리벽 (Yu 2026-05-04 결정, 절대 위반 금지)

- **Paper 2 = HT-overlap PTC ONLY** (Hashimoto)
- **Paper 4 backlog = Korean GD HLA**
- 두 paper 사이 cross-inference / shared framing 금지
- Paper 2 forbidden words list 준수 (`v18_paper2_HT_isolated` memory 참조)

### 1.3 Voice-protected sections (작가 키보드 only)

다음은 **Claude 가 prose 생성 금지**:
- Hook (intro 1.1)
- Aim
- Discussion 3.1
- Limitations
- Cover letter Para 1
- Reviewer Q9

Marathon mode 에서 "고고", "faster", "다 해줘" 는 voice-protected sprint generate 권한 NOT — scaffolding/infra 만 default. (`v17_sprint_vs_marathon_violation` memory)

### 1.4 Marathon mode 규칙

- 2026-05-04 ~ 2026-06-13 = manuscript writing
- 새 분석은 **paper-blocking 만**
- Scaffolding/infra 작업은 항시 OK (CLAUDE.md, web page, figure, report 등)
- 오늘 spatial-full 작업은 advisor 질문 응답 + Paper 2 supplementary 보강 → **paper-blocking 범주에 부합**

---

## 2. Spatial-full 패키지 — 처음부터 v1 → v12 timeline

### 2.1 Trigger

User: "공간 전사체 finding이 뭔가 부족한듯" → "풀 패키지 해주고 다해주고 그리고 왜 이게 나왔는지 설명해주고 공간전사체 finding 강조에서 우리 웹페이지 알지 거기에 올려" (advisor Q1/Q2 응답 목적)

### 2.2 Cohort 구성

| Dataset | n samples | n spots | 조성 |
|---|---|---|---|
| GSE250521 | 16 | ~70k | 4× PT (normal) + 4× PTC + 4× LPTC (locally PTC) + 4× ATC |
| GSE230424 | 4 | ~12k | 4× PTC + HT (Hashimoto-overlap) |
| GSE248205 | 8 | ~8k | CONTROL + HT + GD (autoimmune-only, no PTC) |
| **합계** | **28** | **~90,471** | TROP2-scored / multi-axis Moran I |

### 2.3 Script & figure inventory

| Version | Script | figures | 핵심 |
|---|---|---|---|
| v1 | `spatial_full_package_2026_05_06.py` | S_F1–S_F6 (6) | TROP2 niche identification + HT TLS spatial 기초 |
| v2 | `..._v2_2026_05_07.py` | S_F7–S_F18 (12) | TROP2 16-slide maps · ranked Moran I · co-localization · master summary |
| v3 | `..._v3_2026_05_07.py` | S_F19–S_F26 (8) | DM1 16-slide · microenvironment · niche cluster · mega summary |
| v4 | `..._v4_2026_05_07.py` | S_F27–S_F32 (6) | RGB overlay · driver genes · GSE230424 TROP2 · closure LOSO · niche spacing · ATC degradation |
| v5 | `..._v5_2026_05_07.py` | S_F33–S_F40 (8) | RAI8 · Epi/Prol · GSE248205 full · combined immune · QC · microenv trajectory · cross-axis Moran |
| v6 | `..._v6_2026_05_07.py` | S_F41–S_F48 (8) | 28-sample mega · HT individual markers · Hypoxia/CAF/EMT 16-slide · niche threshold |
| v7 | `..._v7_2026_05_07.py` | S_F49–S_F56 (FIX +8) | S_F16/25/30 fix + 4 KEY FINDINGS collage |
| v8 | `..._v8_2026_05_07.py` | S_F57–S_F64 (8) | violin · RGB · CDF · cluster IDs · closure heatmap · joint scatter · table · master matrix |
| v9 | `..._v9_2026_05_07.py` | S_F65–S_F72 (8) | TLS individual · similarity · PCA · trajectory · compactness · ultimate summary |
| v10 | `..._v10_non_trop2.py` | S_F73–S_F80 (8) | 8 non-TROP2 deep gene sets (Thyroid TF · Cell cycle · TS · Tcell · Mac · Stromal · Immune escape · Glycolysis) |
| v11 | `..._v11_pathways.py` | S_F81–S_F88 (8) | 8 cancer pathways (WNT · NOTCH · Hippo · RAS-MAPK · PI3K-AKT · MYC · TGF-β · NF-κB) |
| v12 | `..._v12_findings.py` | S_F89–S_F100 (12) | 8 new biology axes + stage gradient + ATC-2 deep dive + cross-axis Jaccard + centennial |

**총 100 figures**, 12 scripts, 모두 `project/papers_hub_2026_05_04/assets/spatial_full/` 에 PNG, mirror 됨 → `/var/www/papers/papers_hub_2026_05_04/assets/spatial_full/`. HTTP 200 sweep = **100 / 100**.

### 2.4 Web page

- `project/spatial-full/index.html` (94.5 KB) → mirror `/var/www/papers/spatial-full/index.html`
- 섹션: §1–§9 + §8.5 ~ §8.12 (v3–v12)
- 각 figure 2회 참조 (anchor + img tag), 모두 light-theme 톤

### 2.5 Git commits 오늘

```
a4fbf76 spatial-v12: 12 more figures (S_F89-S_F100) — finding sweep
09e0822 spatial-v9+v10+v11: 24 more figures (S_F65-S_F88)
8788a56 spatial-v7+v8: FIX S_F16/25/30 + 16 new figures (S_F49-S_F64)
5d94e52 spatial-v6: 8 more figures (S_F41-S_F48)
a9060e8 spatial-v5: 8 more figures (S_F33-S_F40)
d986970 spatial-v4: 6 more figures (S_F27-S_F32)
096f99b spatial-v3: 8 more figures (S_F19-S_F26)
bc06416 spatial-v2: 12 new figures (S_F7-S_F18)
c5edf6b spatial-full: 6-analysis package (v1)
```

---

## 3. Findings (signature-level exploratory only)

### 3.1 Advisor Q1 답 — TROP2 niche tumor-specific (S_F1–S_F18)

- GSE250521 PTC + LPTC **8/8 sample** Moran's I 0.27 ~ 0.58
- PT (normal) 4/4 + autoimmune-only 8/8 모두 Moran's I < 0.13
- 즉 TROP2 = **scattered tumor population 이 아니라 spatial niche**, bulk 가 평균화시켜 못 본 것
- ATC 4/4 중 1만 niche (I=0.292) — late-stage transcriptional collapse 가능성

### 3.2 Advisor Q2 답 — 공간 전사체 6가지 distinct reason

1. DM1 region heterogeneity 직접 확인 (§3)
2. TROP2 scattered vs niche 분리 (§1)
3. **HT-overlap PTC TLS niche 직접 확인** (§2; 4/4 슬라이드 Moran's I 0.27 ~ 0.84) — Paper 2 Pillar II 보강
4. Driver-orthogonal axis spot-level (§3)
5. H&E → DM1 closure NO-GO self-evidence (§8)
6. Autoimmune-only negative control (§5)

### 3.3 신규 finding (v12, 2026-05-07) ★

| # | Finding | 핵심 수치 | 위치 |
|---|---|---|---|
| 1 | **Thyroid differentiation 공간 붕괴** | Moran I (TG/TPO/TSHR/SLC5A5/DUOX1/2/IYD/DIO1/DIO2): PT 0.46 → PTC 0.57 → LPTC 0.61 → **ATC 0.12** (Δ=−0.34). ATC-1/2/3 ≈ 0 | S_F91, S_F97 |
| 2 | **24-axis stage gradient** | Falling: `Thyroid_differen` (Δ −0.34), `Thyroid_TF` (−0.20). Rising: `Cell_cycle` (+0.38), `T_cell` (+0.33), `Senescence_SASP` (+0.26), `Macrophage_TAM` (+0.25), `Tumor_suppressor` (+0.24) | S_F97 |
| 3 | **ATC-2 super-organizer 8-axis** | GSM7980873: Cellular_stress 0.84 / Hypoxia 0.80 / EMT 0.73 / Apoptosis_DDR 0.69 / Senescence 0.64 / Stemness 0.59 / Angiogenesis 0.44 / Cell_cycle 0.79 — 한 슬라이드에 7+ program 동시 spatial organize | S_F98 |
| 4 | **LPTC-3 두 번째 outlier 후보** | Hypoxia 0.68 / EMT 0.72 / Stemness 0.59 / SASP 0.55 / Cellular_stress 0.74 | S_F89, S_F90, S_F92, S_F93, S_F96 |

### 3.4 Closure battery (Paper 2A NO-GO 재확인)

- LOSO ResNet50 H&E → DM1 ρ ≈ 0.06 (기대 약 0.6+) → spatial direct evidence 가 closure 실패의 self-evidence (S_F30 fixed schema, S_F61)

### 3.5 기술 fix

- **S_F25 Visium hex 6-neighbor adjacency** — `visium_hex_clusters()` union-find with offsets `[(-1,-1),(-1,1),(0,-2),(0,2),(1,-1),(1,1)]`. scipy.ndimage.label default 4-conn 은 Visium hex 그리드에서 작동 X. Fix 후 max=162, mean=3.
- **S_F30 LOSO schema** — y_obs/y_pred/target 스키마로 변경 (DM1_like_score 컬럼 직접 X).
- **S_F16 split** — HLA_II/B_cell/TLS/IGHV gene sets 는 PTC_HT only 측정 → 4×8 condition matrix sparse → 2-panel split (TROP2 cross-cohort + HT immune PTC_HT-only).

---

## 4. 기술 인프라 상태

### 4.1 Disk

- Root 123G, `/data` 512G premium SSD bind mount (CLAUDE.md 의 layout)
- `project/results` ↔ `/data/thca/repo_results`
- `project/data` ↔ `/data/thca/repo_data`
- 100 figures + 모든 TSV / h5ad / scoring matrix → `/data` 사이드, root 안전

### 4.2 RunPod

- 활성: `thca-spark-dm-a6000-v4` (uvp9i2r9s6l85y), RTX A6000 48GB, $0.33/hr (community cloud)
- L40S v3 → A6000 v4 교체 (2026-05-07): SSH wedged + L40S host unavailable → terminate, 50GB volume lost
- 50GB `/runpod-volume` = pod-specific persistent (host 이동 X, terminate 시 손실)

### 4.3 Web

- `papers_hub_2026_05_04/master_view.html` — 8-paper 통합 dashboard (port 80)
- `spatial-full/index.html` — 100 figures + Q1/Q2 advisor brief
- Neoantigen vaccine hub (별도 line): 119,237 master / 13 sources / 31 analyses / 7 live pages

### 4.4 Git

- Branch: main
- Today HEAD: `a4fbf76` (spatial-v12)
- 미커밋 변경 다수 (`git status` 에 manuscript_v8/* + 새 root MD 들 — voice-protected 근접, 본인 키보드 산출물 가능성)

---

## 5. 메모리 (자동, 2026-05-07 시점)

핵심 활성 (`MEMORY.md` 인덱스 기준):

- `paper_numbering_2026_05_04` (Paper 1/2/3/4 정의)
- `v18_paper2_HT_isolated` (Yu 2026-05-04 분리벽)
- `v19_paper4_GD_backlog` (Korean GD → Paper 4)
- `v17_marathon_mode_post_pillar1` (마라톤 모드)
- `v17_sprint_vs_marathon_violation` (voice-protected 규칙)
- `v17_npj_ship_status` (Paper 1 ship-ready)
- `v19_paper3_track_b_lite_2026_05_06` (ICI Track B-lite findings)
- `v17_landa2016_cite_save` (Disc 3.1 cite 정정 fact stack)
- `v19_runpod_a6000_2026_05_07` (오늘 pod 교체)

---

## 6. 결정해야 하는 것 (사용자 판단 필요) ★

이 섹션이 본 문서의 핵심 — Claude 가 임의로 진행 X.

### 6.1 즉시 (오늘 / 이번 주)

- **D1. v12 신규 finding 의 manuscript 반영 범위**
  - 옵션 A: Paper 2 supplementary 추가 (HT TLS 보강은 이미 있음; 분화 공간붕괴 + stage gradient + ATC-2 = Paper 2 핵심 vs 별도?)
  - 옵션 B: 별도 short paper / preprint 로 분리 (signature-level only 라 Cell Rep Med 어려움; bioRxiv 정도)
  - 옵션 C: 현 단계 = web page 시각자료로만 보존 (manuscript 안 들어감), 6주 마라톤 마감 후 재평가
  - **판단 기준**: 분화 공간붕괴 finding 이 Paper 1 main story (8-gene + DM1) 를 흐릴 가능성 vs 보강할 가능성. Paper 1 ship-ready 직전이라 새 finding push 가 timing risk.

- **D2. Paper 1 outreach 4건 발송 timing**
  - npj 1차 ship-ready, 4 outreach drafts NOT 발송 — user click 대기 상태가 며칠 째인지 확인 필요
  - spatial-full 에 시간 쏟느라 발송 대기 길어졌을 가능성

- **D3. ATC-2 / LPTC-3 outlier deep dive 추가 분석 여부**
  - ATC-2 (GSM7980873) 단일 슬라이드 8-axis 동시 = 매우 unusual
  - LPTC-3 도 비슷 — n=2 이지만 single-sample finding 으로 정리해서 Paper 2 case study 가능성?
  - **판단 기준**: signature-level only 규칙 vs case-study reportable level

### 6.2 중기 (이번 주 ~ 다음 주)

- **D4. 100 figures → 압축 / 큐레이션**
  - 100 figures 너무 많음. 본문 / supplementary / web-only 3-tier 분류 필요
  - 본문 (Paper 2): 5–8장 권장
  - Supplementary: 20–30장
  - Web-only / archive: 나머지 60+
  - **판단 기준**: figure budget × Paper 2 journal target (Sci Rep base / Cell Rep Med reach)

- **D5. Discipline rule 재확인 — signature-level only 유지**
  - 100 figures 가 풍성해질수록 causal / clinical / biomarker claim 으로 슬라이드 위험
  - Voice-protected section (Hook/Aim/Disc 3.1/Limitations) 에서 finding 인용 시 톤 본인 키보드 only

- **D6. Paper 4 (Korean GD HLA) backlog 분리 유지**
  - spatial-full 에 GD samples (GSE248205) 들어있어 cross-paper 유혹 발생 가능
  - 분리벽 strict 유지 = Paper 4 backlog 시그널 spatial 에서 negative control 로만 사용 (현 상태 OK)

### 6.3 장기 (마라톤 마감, 2026-06-13)

- **D7. 6주 후 Paper 2 submission target journal 확정**
  - Sci Rep base / Cell Rep Med / JCI Insight reach
  - GSE286332 (PTC vs PTC+HT) STRONG GO 가 reach venue gate
  - Spatial-full 100 figures 가 reach venue 에 adequate 한지 평가

- **D8. v12 finding 의 추가 cohort replication 필요 여부**
  - 분화 공간붕괴 (ATC) — n=4 ATC 만으로 claim 부족 → independent ATC spatial cohort?
  - ATC-2 outlier — n=1 → reproducibility 미정
  - **판단 기준**: replication 비용 vs reviewer pushback 위험

- **D9. Pod cost 관리**
  - A6000 $0.33/hr — 6주 중 추가 burst 필요 시 budget 재확인
  - 마라톤 모드 = 분석 stop / writing only 라 pod 가동 최소화 권장

### 6.4 Discipline check (오늘 작업의 self-audit)

| Rule | Compliance | 비고 |
|---|---|---|
| Marathon mode = scaffolding/infra | ✓ | spatial-full = advisor 질문 응답 + Paper 2 보강 = paper-blocking |
| Voice-protected 무손상 | ✓ | Hook/Aim/Disc 3.1/Limitations/Cover Para 1/Q9 모두 미수정 |
| Yu 분리벽 (Paper 2 HT vs Paper 4 GD) | ✓ | spatial-full 에 GD 별도 cohort 로 표시, cross-inference X |
| Signature-level only | ✓ | "tumor-specific niche", "spatial coherence" 표현 사용; "biomarker", "drug target" 사용 X |
| Paper 1 main story 무손상 | ✓ | npj submission 미수정 |
| Paper 3 Track A FROZEN | ✓ | chmod 444 유지 |

---

## 7. 미해결 / 잠재 리스크

- **R1.** 100 figures 가 advisor 보기에 too much / firehose 가 될 가능성 → curation tier 사전 배포 권장 (D4 와 연결)
- **R2.** 분화 공간붕괴 finding 이 Krishnamoorthy 2025 / Landa 2016 등 기존 ATC dedifferentiation literature 와 conflict / overlap 미체크
- **R3.** Visium hex adjacency fix (S_F25) 이전 figure (v1–v6) 에 잘못된 cluster 통계가 들어가 있을 가능성 — 재계산 spot check 권장
- **R4.** GSE230424 n=4 PTC+HT 만으로 TLS niche claim — Paper 2 reviewer 가 n 작다 지적 시 답변 필요
- **R5.** Marathon 1주차 (5/4–5/10) 작업 시간 중 spatial 분석에 거의 다 사용 → manuscript writing 진척도 별도 점검 필요

---

## 8. 다음 step 제안 (Claude 가 user judgment 대기 중)

순서대로 user click 대기:

1. D1 (v12 manuscript 반영 범위) 결정 → 진행 방향
2. D4 (100 figures curation tier) 결정 → 본문/supp/web 분리
3. D2 (Paper 1 outreach 발송) 진행 여부
4. R5 (manuscript writing 진척도 점검) — 이번 주 회복 가능한지 self-assessment

위 4개 답이 정해지면 마라톤 1주차 마감 (5/10 일요일) 에 정리 가능.

---

## 9. Files / data 빠른 인덱스

```
project/
├── notebooks_or_scripts/
│   └── spatial_full_package_{v1..v12}*.py            # 12 scripts
├── papers_hub_2026_05_04/
│   └── assets/spatial_full/
│       └── S_F{1..100}*.png                           # 100 figures
├── spatial-full/
│   └── index.html                                     # 94.5 KB, §1–§9 + §8.5–§8.12
└── results/spatial_full_2026_05_06/                   # bind mount → /data
    ├── spatial_per_spot_TROP2.tsv.gz                  # 90,471 spots
    ├── spatial_F40_per_sample_per_axis_morans.tsv
    ├── spatial_F64_28sample_multiaxis_morans.tsv
    ├── spatial_v10_deep_morans.tsv
    ├── spatial_v11_pathway_morans.tsv
    ├── spatial_v12_new_morans.tsv
    ├── spatial_v12_stage_gradient.tsv
    ├── spatial_v12_axis_slopes.tsv
    └── spatial_v12_cross_axis_jaccard.tsv

/var/www/papers/papers_hub_2026_05_04/assets/spatial_full/   # 100 PNG mirror
/var/www/papers/spatial-full/index.html                       # page mirror
```

— 이 문서 (`SITUATION_FULL_2026_05_07.md`) — project root, 작성: 2026-05-07.
