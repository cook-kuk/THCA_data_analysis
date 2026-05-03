# H&E → 8-gene DM1 axis 예측 — feasibility audit & staged plan

**Run:** 2026-05-03 · **Owner:** Seungho Cook · **Path:** `/home/seungho/personal/THCA_data_analysis/PATHOLOGY_DM1_FEASIBILITY_2026_05_03.md`
**For:** scp 로컬 검토 · senior scientist (compu-path + ST + cancer genomics) framing

---

## 0. 한 줄 결론

> 풀 파이프라인을 한 번에 돌리는 건 **GPU 부재 + TCGA WSI 500 GB + HF gated 모델 access** 때문에 비현실적. **Phase A (GSE250521 단독 ResNet50 LOSO)** 가 첫 gate. 거기서 LOSO Spearman ≥ 0.30 / DM1_high AUROC ≥ 0.70이 안 나오면 Paper 2 / IP 양쪽 다 즉시 no-go. 통과 시에만 B (멀티 코호트) → C (TCGA WSI subsample) → D (go/no-go report) 순차 진행.

---

## 1. 가설 & 목표 (사용자 정의)

**Hypothesis.** ST-aligned H&E morphology can predict a depth-residualized DM1_like score, and the resulting image-DM1 score recovers the RNA-defined DM1 subtype in TCGA-THCA.

**Targets**
- Primary: `DM1_like_score_resid`
- Secondary: `RAI_8_score_resid`, `THYROID_NONOVERLAP_score_resid`
- Binary: `DM1_high` = top quartile (within sample / cohort)

**Strict rules** (사용자 spec)
1. depth-residualized label만 사용
2. random tile split 금지 (LOSO / LODO 만)
3. Required validation: leave-one-slide-out, leave-one-dataset-out, PTC-only sensitivity (가능 시)
4. RET fusion 진단 주장 X — molecular triage / DM1-like program prediction만
5. random 8-gene panel + housekeeping panel으로 negative control

---

## 2. 현재 상태 inventory (정확히 뭐가 있고 뭐가 없는지)

### 2.1 이미 완료된 것 ✅

| 자산 | 경로 | 비고 |
|---|---|---|
| GSE250521 16 slide parse | `project/data/processed/GSE250521/*.raw.h5ad` | 55,873 spots post-QC |
| 8-gene + 보조 score 계산 | `project/results/01_spatial_score/all_spots_scored.tsv.gz` | RAI_8 / DM1_like / TDS_like / CAF_ECM / EMT / Hypoxia / Proliferation / Epithelial |
| Per-sample spot scores | `project/results/01_spatial_score/per_sample/GSM*.tsv.gz` (16 files) | |
| Stage trend + depth-corrected | `project/results/02_stage_trend/depth_corrected_trend.csv` | epi top-50% RAI ρ=−0.35, p=0.18 (n=16 underpowered) |
| **Tile metadata (ready for embedding)** | `project/results/03_pathology_poc/tile_metadata.tsv.gz` | **3,200 tiles** balanced 800/stage; columns: spot_id, sample_id, stage, tile_path, x_hires, y_hires, RAI_8_score, DM1_like_score, TDS_like_score |
| Tile extraction script | `project/src/05_pathology_poc/extract_tiles.py` | hires PNG에서 spot-centered crop |
| TCGA-THCA bulk RNA-seq (이미 있음) | `/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv` | external molecular validation에 즉시 사용 가능 |
| TCGA-THCA bulk z-score | `/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv` | |
| Pathology POC scaffold doc | `project/src/05_pathology_poc/README.md` | embedding/training **skeleton (미구현)** |
| Advisor brief | `project/reports/gse250521_dm1_spatial_validation_report.md` | "Pathology POC is GO… do not run on this single-dataset POC alone — wait for HRA003537/GSE230424 ingestion" 게이트 명시 |

### 2.2 빠진 것 ❌

| 자산 | 비고 |
|---|---|
| GSE230424 ST raw (4 PTC+HT) | 미다운로드 |
| GSE248205 ST raw (8 autoimmune control) | 미다운로드 |
| HRA003537 (~30 Chinese ST, advisor brief에서 언급) | 미다운로드 — 권한 issue 가능 (CNGB/HRA China access) |
| TCGA-THCA WSI (~500 slides × ~1 GB ≈ **500 GB**) | `/data/thca/data_raw/gdc/TCGA-THCA/` 디렉터리만 존재, WSI 자체 없음 |
| Pathology foundation model 가중치 | UNI / CONCH / Virchow2 / Prov-GigaPath 모두 HF gated, **token + acceptance 필요** |
| GPU | `nvidia-smi` 미설치 — **현재 머신 GPU 없음** |
| Embedding code (실구현) | scaffold만, 학습 루프 미작성 |
| TCGA WSI tile inference 파이프라인 | 미구현 |
| Slide-level aggregation (mean / top-k mean / high-risk fraction) | 미구현 |

### 2.3 자원 제약 (가장 중요)

| 자원 | 상태 | 영향 |
|---|---|---|
| `/data` 가용 디스크 | **454 GB** | TCGA WSI 500 GB 다 받으면 거의 풀; subsample 필요 |
| `/` (root) | 9.5 GB | 모델 weight (~3-5 GB)는 `/data`에 둬야 함 |
| GPU | **없음** | UNI/CONCH (ViT-L) 추론: GSE250521 3,200 tiles → CPU 30-60분 OK; TCGA WSI 25M tiles → CPU 사실상 불가 (수 주) |
| HF token | 미확인 | UNI 등 access 신청 후 수 시간~며칠 |
| 외부 GPU burst (Azure $4.80 패턴, memory에 success) | 옵션 가능 | UNI TCGA-WSI 추론에는 사실상 필수 |

---

## 3. Marathon mode 제약 (memory 기반)

- 2026-04-30 결정: Pillar 1 STRONG = **분석 끝**. 5/4-6/13 6주 manuscript writing marathon. **새 분석은 paper-blocking만**.
- 본 프로젝트는 **신규 paper feasibility audit** 성격 → 풀 파이프라인 한 번에 진행은 마라톤 위반.
- Voice-protected section sprint는 별개 — 코드/임베딩/그림 생성은 마라톤 모드에서도 OK, 단 **gate를 통과한 후의 단계만**.

**해석**: Phase A는 "Paper 2 / IP feasibility를 결정짓는 1회 분석"으로 정당화 가능. Phase B/C는 A의 결과로 ROI가 정당화될 때만 진입.

---

## 4. Paper 경계 (memory v18 Paper 2 isolation)

- Paper 2 = **Hashimoto-overlap PTC ONLY** (Yu 2026-05-04 결정). Pillar I v2 = Korean PTC vs Korean baseline (AFND).
- 본 H&E → DM1 axis는 **Paper 1의 8-gene driver-excluded sub-stratifier** territory.
- **Decision**: 본 분석이 strong이면 → Paper 1 reach venue (Cell Rep Med / JCI Insight) 보강 figure가 1순위. "Paper 2 viable"는 **현재 정의된 Paper 2 (HT-overlap)**가 아니라 **신규 별도 paper (image-DM1 triage)** 의미.
- **금지어 가드 (Paper 2)**: HT/PTC overlap / TLS / BCR / AICDA — 본 분석에 포함 X
- **금지어 가드 (Paper 3)**: GD / HLA / Graves' — 포함 X

---

## 5. Staged 실행 plan

### Phase A — GSE250521 단독, ResNet50 baseline (gate 1)

**시간** ~30-60분 CPU · **디스크** 무시할 수준 · **외부 access** 없음

**할 것**
1. ResNet50 ImageNet (torchvision) frozen → 3,200 tile (224 px) 임베딩 (2048-dim)
2. Ridge / ElasticNet → `DM1_like_score_resid` 회귀, **leave-one-slide-out** (16 fold)
3. Metrics: per-fold Spearman ρ, Pearson r, R², `DM1_high` (top quartile) AUROC
4. **Negative controls** (사용자 spec):
   - Random 8-gene panel (genome-wide, balanced detection rate matched) → 같은 파이프라인
   - Housekeeping panel (ACTB, GAPDH, B2M, HPRT1, PPIA, RPL13A, RPLP0, TBP) → 같은 파이프라인
5. Multi-resolution: 224 / 448 / 672 px tile 비교 (사용자 spec)
6. PTC-only sensitivity (4 slides) — 데이터 크기 작아 LOSO ≈ 4-fold만 가능, 보고로

**Gate (Phase B 진입 조건)**
- LOSO mean Spearman ρ(`DM1_like_score_resid`) ≥ **0.30**
- LOSO `DM1_high` AUROC ≥ **0.70**
- Random panel ρ ≤ 0.10 (signal specificity)
- Housekeeping panel ρ ≤ 0.10

**산출물**
- `project/results/03_pathology_poc/embeddings_resnet50_224.npz`
- `project/results/03_pathology_poc/loso_metrics_resnet50.tsv`
- `project/results/03_pathology_poc/loso_pred_vs_obs.png`
- `project/results/03_pathology_poc/negative_controls_summary.tsv`

**Risk**
- ResNet50 ImageNet은 weak baseline (UNI/CONCH 대비 보통 2-3× 약함). 0.30 임계가 너무 빡셀 수 있음 → **조기 false-no-go 위험**
- Mitigation: ρ 0.20-0.30 borderline이면 UNI access 신청 + 재시도 1회 허용

### Phase B — 멀티 코호트 LODO (gate 2)

**시간** ~2-4시간 (다운로드 + parse + score + retrain) · **디스크** ~10-20 GB · **외부** GEO 다운로드만

**할 것**
1. GSE230424 (4 PTC+HT) + GSE248205 (8 autoimmune control) 다운로드 + parse (이미 있는 `01_parse_visium` 패턴 재사용)
2. 같은 8-gene + housekeeping + random scoring + depth residualization
3. Combined dataset = GSE250521 (16) + GSE230424 (4) + GSE248205 (8) = **28 slides, ~3 datasets**
4. **leave-one-dataset-out**: train 2 datasets / test 1 dataset, 3 fold
5. PTC-only sensitivity: GSE250521 PTC+LPTC (8) + GSE230424 PTC+HT (4) = 12 slides LOSO

**Gate (Phase C 진입 조건)**
- LODO Spearman ρ ≥ **0.25** (cross-dataset 약간 완화)
- 적어도 2/3 dataset에서 양의 신호

**산출물**
- `project/results/04_pathology_lodo/lodo_metrics.tsv`
- `project/results/04_pathology_lodo/per_dataset_pred_vs_obs.png`

**Paper 2 isolation 주의**
- GSE230424 (PTC+HT) parse하면 자연스럽게 HT 신호 있음. 본 분석에서는 **DM1_like / RAI_8 score만 회귀 target**. HT signature / TLS / BCR 점수 계산 X. 이는 Paper 2 영역.

### Phase C — TCGA-THCA WSI external validation (gate 3)

**시간** ~4-8시간 (subsample) ~ 수일 (전체) · **디스크** subsample 50 GB / 전체 500 GB · **외부** GDC manifest

**할 것 (subsample 우선)**
1. TCGA-THCA WSI manifest에서 **DM1 high vs low로 stratified subsample 50 slide** (bulk 8-gene RNA score 기준)
2. GDC client로 subsample 다운로드 (~50 GB)
3. WSI → tile (5x or 10x at 224 px non-overlap, tissue mask 적용) — `openslide` + 기본 thresholding
4. Phase A/B에서 학습된 모델로 tile-level DM1_like 예측
5. Slide-level aggregation 3가지 (사용자 spec):
   - mean
   - top-k mean (k = top 10% tile)
   - high-risk fraction (tile score > slide-internal threshold)
6. **External validation**: image-DM1 (slide level) ↔ bulk RNA 8-gene score (이미 `/data/thca/data_processed/bulk_rnaseq/`에 있음) Spearman + AUROC for DM1/DM2 subtype
7. **Enrichment test**: BRAF/RAS-negative dark matter (memory v17 dark_matter_phase2 활용 가능) + fusion-positive (`/data/thca/data_raw/v3_ext/tcga_fusions`) 에서 image-DM1 score 분포 비교

**Gate (Phase D 진입 조건)**
- Image-DM1 vs bulk 8-gene score Spearman ρ ≥ **0.30**
- DM1 vs DM2 subtype AUROC ≥ **0.65**

**산출물**
- `project/results/05_tcga_wsi_external/slide_level_image_dm1.tsv`
- `project/results/05_tcga_wsi_external/image_vs_rna_correlation.png`
- `project/results/05_tcga_wsi_external/dark_matter_enrichment.tsv`

**Risk**
- TCGA-THCA WSI는 frozen section + FFPE 혼재. ST training (GSE250521도 frozen vs FFPE 확인 필요)과 도메인 갭 큼. UNI/CONCH가 더 robust하지만 ResNet50은 여기서 무너질 가능성 큼 → Phase C는 **UNI/CONCH로 다시 학습한 후 진입 권장**

### Phase D — 통합 go/no-go report

**할 것**
1. A/B/C 결과를 single report에 통합
2. **Paper 2 viability**:
   - 신규 image-DM1 paper로 viable한가? (Phase C가 strong이어야)
   - 현재 정의된 Paper 2 (HT-overlap)와의 boundary 명확화
3. **Business / IP viability**:
   - prior art 검색 (UNI/CONCH 학습기반 thyroid 예측 연구 — `v17_lit_enrich` 패턴 재사용)
   - claim defensibility: tile-level DM1 prediction + multi-cohort LODO + external bulk RNA validation = patentable?
4. **External validation gap**:
   - HRA003537 (Chinese ST) access 시도 가치
   - Bundang FFPE cohort (memory: outreach-stage, no data yet) → 본 모델의 prospective FFPE pilot으로 매력도 ↑

**산출물**
- `project/reports/PATHOLOGY_DM1_GO_NOGO_2026_05_XX.md`
- 결정: GO (어느 venue로) / PIVOT / NO-GO

---

## 6. 4가지 진행 옵션 (user 선택)

| Option | 설명 | 예상 wall time | 디스크 | 외부 access |
|---|---|---|---|---|
| **1. A만 지금** | ResNet50 baseline → gate 결과 보고 | 30-60분 | ~0 | 없음 |
| **2. A→B 자동** | A 통과 시 GSE230424/248205 다운로드 → LODO까지 | 3-5시간 | +20 GB | GEO만 |
| **3. 풀 강행** | A→B→C(TCGA WSI subsample 50 slide) | 8-12시간 | +70 GB | GDC token |
| **4. UNI 우선** | HF token + UNI access 신청 → A부터 UNI로 시작 | 신청 대기 + A 60분 | ~0 | HF gated |

**senior scientist 추천**: **Option 1 → 결과 보고 후 2 또는 4 결정**. 이유:
- ResNet50 baseline이 Phase A에서 ρ ≥ 0.30 못 만들면 UNI도 Paper 2/IP 관통할 강도 안 나올 가능성 높음 (UNI gain ≈ 2× 수준)
- 만약 0.20-0.30 borderline이면 UNI access 신청 + 재시도 (Option 4 진입)
- 만약 ≥ 0.40 strong이면 바로 Option 2 (멀티 코호트)

---

## 7. User 확인 필요 항목

1. **HF token / UNI·CONCH·Virchow2 access** 신청 상태 (없으면 ResNet50로 Phase A 강행)
2. **Azure GPU burst** ($4.80 success 패턴, memory에 있음) 사용 의향 — TCGA WSI Phase C에 사실상 필수
3. **TCGA WSI subsample 50 slides vs 전체 500 GB** — subsample 권장
4. **HRA003537 (Chinese ST)** access 시도 의향 — 시간/외교 비용 큼
5. **Marathon mode 재확인** — 본 분석이 "paper-blocking 신규"로 분류되는지, 사용자 OK 사인 필요

---

## 8. Paper isolation 가드 (본 plan 전체)

- ✅ Paper 1 territory (8-gene driver-excluded sub-stratifier; image triage)
- 🚫 Paper 2 (HT-overlap) 침범 X — TLS / BCR / AICDA / IGHV / HT signature 회귀 target X
- 🚫 Paper 3 (GD / HLA) 침범 X — HLA / Graves / 면역 점수 X
- ⚠️ GSE230424 다운로드는 Paper 2 dataset이지만, 본 plan에서는 **DM1_like / RAI_8 score regression target만 사용**

---

## 9. Reproducibility (Phase A 실행 명령 — gate 통과 후 추가 단계 문서화)

```bash
cd /home/seungho/personal/THCA_data_analysis

# Phase A (Option 1)
python3 project/src/05_pathology_poc/extract_tiles.py --tile-sizes 224,448,672  # multi-res
python3 project/src/05_pathology_poc/embed_resnet50.py                          # NEW
python3 project/src/05_pathology_poc/train_loso_ridge.py                        # NEW
python3 project/src/05_pathology_poc/negative_controls.py                       # NEW
python3 project/src/05_pathology_poc/plot_pred_vs_obs.py                        # NEW
```

---

## 10. 한 줄 next action

**user**: Option 1/2/3/4 중 선택 + HF token / Azure GPU 의향만 답변 → 즉시 실행. 본 문서가 그 선택의 근거.

---

*generated 2026-05-03 by Claude (Opus 4.7) under marathon-mode discipline*
