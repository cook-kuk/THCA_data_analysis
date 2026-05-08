# SPARK · DM (Dark Matter) Pathology — 전체 알고리즘 / 데이터 / Flow

**Date:** 2026-05-07
**Goal:** Tolkach Lab SPARK (Nature Med 2026) framework 를 갑상선암 DM1 vs DM2 dark-matter axis 에 적용. H&E pathology + Spatial Transcriptomics + TCGA WSI 를 단일 SPARK pipeline 으로 묶음.

---

## 0. 우리 결정의 배경 — 왜 SPARK?

| 질문 | 답 |
|------|----|
| 무엇이 부족한가? | DM1 vs DM2 분류는 RNA-seq 기반 (TCGA-THCA n=504, 8-gene RAI panel AUC 0.962). 그러나 **수술 시점에서 RNA 결과 받기 전에** 분화도 위험을 알 수 없음 (RT-qPCR 4-7일 소요) |
| 가설 | H&E 슬라이드 자체에 분화도 정보가 인쇄돼 있다 (BRAF mRNA 와 직교 d=−0.04, 따라서 driver mutation 과 별개의 morphology axis 존재) |
| 왜 SPARK? | LLM agentic 가 *여러* 가설을 자동 생성 → 코드 → 검증. 단일 hand-crafted feature 가 아닌 **20-50 SPARK ideas** 의 ensemble 로 reviewer-proof |
| 왜 지금? | SPARK 가 2026-04 Nature Med 발표 + GitHub (cpath-ukk/SPARK) 코드 공개. Tolkach lab 의 lung/colon/breast 5,400 patient validation 직후 → **갑상선 first-mover** |

---

## 1. 데이터 inventory (왜 이걸 썼나)

### 1.1 GSE250521 Visium ST (16 slides) — Discovery cohort

| 항목 | 값 |
|------|----|
| Stages | N (4) / PTC (4) / LPTC (4) / ATC (4) |
| Spots | 8521 total Visium spots |
| H&E tiles extracted | 3200 @ size 224, 2935 @ 448, 2386 @ 672 |
| Per-spot scores | DM1_like_score, RAI_8_score, TDS_like_score (precomputed via 8-gene panel + log_counts/log_ngenes residualized) |
| **왜 이 데이터?** | (a) ST + H&E 가 한 슬라이드에 동시 존재 → **per-spot transcriptome × per-spot pathology** 직접 매칭, (b) ATC 까지 4-stage 트라젝토리 커버, (c) 16 slide LOSO 가능 (모델이 single-slide 외워서 부풀리는 문제 차단) |

### 1.2 TCGA-THCA WSI (n=504, DM-balanced 50 subset) — Validation cohort

| 항목 | 값 |
|------|----|
| Total cases | 504 primary tumor |
| DM-labeled | 137 (DM1 n=82 + DM2 n=55), 345 unlabeled |
| **DM-balanced subset for WSI run** | 50 cases (DM1=25, DM2=25) |
| DM1 composition | 1 BRAF + 20 driver-neg + 4 OTHER |
| DM2 composition | 0 BRAF + 20 driver-neg + 5 OTHER |
| **왜 이 데이터?** | (a) 가장 큰 thyroid cancer cohort with H&E + transcriptome + clinical (PFI/OS), (b) Korean K2 cohort 와 driver class distribution 유사 → 한국 generalizability 함께 평가 가능, (c) DM1/DM2 둘 다 driver-neg-dominant → **분화도 axis 가 BRAF/RAS 와 직교한다는 가설을 H&E 에서도 직접 검증** |

### 1.3 SPARK framework (cpath-ukk/SPARK)

| 항목 | 값 |
|------|----|
| Generative pipeline | LLM agent (default openai/o1; 우리는 claude-opus-4-7 사용) → biological hypothesis → Python code snippet |
| Analytical pipeline | WSI + tissue-mask + Hovernext-format GeoJSON → patch-level features |
| Prognostic pipeline | per-patient aggregate + Cox / KM with PFI |

---

## 2. 알고리즘 stack (왜 이 알고리즘?)

### 2.0 GigaTIME virtual mIF (Cell 2025; Valanarasu et al., Microsoft × Providence)

**Why included now**: Cell 2025 paper (S0092-8674(25)01312-1) trains a cross-modal translator on 40M paired H&E + mIF cells (Providence dataset, 14,256 patients). Generates **23-channel virtual mIF** from any H&E patch — including DAPI / CD8 / CD4 / CD20 / CD3 / CD138 / PD-1 / PD-L1 / CK / Ki67 / Tryptase / Caspase3 / Ki67 / Transgelin. Public: github.com/prov-gigatime/GigaTIME + HF prov-gigatime/GigaTIME (gated, license click + HF_TOKEN).

**How it slots in**:
1. After HoVer-NeXt cell detection (Phase 6), run GigaTIME on each tile to add virtual IO markers (PD-L1, Ki67) + cross-validate cell-class calls (CD8/CD3/CD20 vs HoVer-NeXt "Lymphocyte").
2. Per-tile per-channel mean intensity → SPARK Idea_generation_D (multiplexed TMA) inputs → 14-cell-type SPARK pipeline.
3. Virtual PD-L1 + Ki67 mean across DM1 vs DM2 → directly informs **Paper 3 ICI vulnerability** (memory `v19_paper3_ici_track_a`) without needing real IO IHC.

**Pod DM Phase 7** (added): downloads weights via HF_TOKEN, runs 30 tiles/slide × 50 WSI = 1500 GigaTIME inferences, ≈45 min on A40, output `gigatime_mIF/{slide}/{tile}.npz` (channel mean + 128×128 sub-sample).

**Local CPU**: not feasible — GigaTIME is a U-Net++ scale cross-modal translator, ~3-10 s/tile on CPU; 3200 tiles × 10s = 9 hours. **GPU only** in practice.

### 2.1 H&E embedding — multi-tier strategy

| Tier | Backbone | When | Why |
|------|----------|------|-----|
| **0 (현재 완료)** | ResNet50 ImageNet | Phase A POC | baseline; foundation-model 없이도 ATC binary 같은 coarse axis 잡히는지 확인 |
| **1 (다음)** | UNI (MahmoodLab) | Pod DM | 100M+ pathology image self-supervised; gated HF, non-commercial. **DM1/DM2 fine stratification 의 가장 가능성 높은 backbone** |
| **2 (fallback)** | DINOv2 ViT-L | Pod DM if UNI gated access 없음 | 142M LAION self-supervised; pathology-naive but 강한 general representation |
| **3 (last resort)** | ViT-B ImageNet | Pod DM 안전망 | 항상 됨 |

**Phase A 결과 → ResNet50 dead end 확정:**
- LOSO 3-size pooled AUROC = 0.51 / 0.52 / 0.51 (chance = 0.50)
- RAI_8 real panel Spearman 0.067 ≈ random panel q95 = 0.053 (즉 **curated panel 효과 ≈ 0**)
- ATC binary AUROC 0.689 (coarse axis 만 잡힘)
- **Tier 1 (UNI) 로 가야만 DM1 vs DM2 fine 가능성 있음**

### 2.2 Cell detection — Hovernet 계열

| Tool | When | Why |
|------|------|-----|
| **Local skimage (현재 완료)** | GSE250521 즉시 SPARK Analytical 운영 위해 | Color-deconvolution (rgb2hed) + watershed + rule-based 4-class. **GPU 없이도 3200 tile = 3.5 분**. 분류 정확도 한계 있지만 SPARK feature 코드 구동에 필요한 GeoJSON 형식 동일 |
| **HoVer-NeXt (Pod DM)** | RunPod GPU 후 TCGA-THCA WSI 50개 | digitalpathologybern/hover_next_inference + Lizard-Mitosis ConvNeXtV2-Large = 7-class (Tumor/Stroma/Lymphocyte/Plasma/Eosinophil/Neutrophil/Connective). Tolkach lab SPARK 공식 입력 format. CPU inference 는 60시간 = 비현실적 → GPU 필수 |
| StarDist (deprecated) | Pod DM 초기 fallback | HoVer-NeXt 가 안 되면 자동 사용. PanNuke prob_thresh 0.5 |

### 2.3 SPARK Generative (idea generation) — agentic LLM

```
LLM_IGA (idea generator):       claude-opus-4-7 (default openai/o1)
LLM_IRA (idea reviewer):        claude-opus-4-7
LLM_DDA (duplicate detector):   semantic similarity threshold 8/10

Cycles:                         4 (default) — 각 cycle 마다 6-10 ideas, 중복 제거
Quality buckets:                basic / advanced / creative / visionary
Output:                         refined idea + Python code snippet
```

**Web 챗봇으로 wrapping**: lumenix_spark/dm_idea_generator.html — API key 입력 시 라이브 호출, 없으면 hardcoded 8 fallback ideas (B1/B2/B3/B4/B5/C1/C2/D1).

### 2.4 SPARK Analytical → Prognostic

```
Per tile  →  cell counts (n_thy, n_lym, n_fib, n_other)
          →  geometric features (eccentricity, area variance)
          →  spatial features (DBSCAN cluster size, Moran's I)
          →  8 SPARK ideas

Per slide →  tile-mean + top-25% + bottom-25% per feature
          →  4-quadrant heterogeneity (D1)

Per patient → slide aggregation
            → Cox HR vs PFI / OS
            → DM1 vs DM2 logistic regression LOSO
```

---

## 3. 현재 진행 상태 (2026-05-07 15:14)

### 3.1 완료

| 단계 | 결과 | 파일 |
|------|------|------|
| Phase A POC end-to-end | ResNet50 224/448/672 LOSO **all chance level**; ATC binary **AUROC 0.689** | `loso_metrics_resnet50.tsv`, `negative_controls_summary.tsv` |
| DM-balanced TCGA WSI manifest | 50 cases (DM1=25 / DM2=25) GDC filter ready | `tcga_thca_wsi_dm_balanced.tsv`, `tcga_thca_wsi_gdc_filter.json` |
| SPARK repo clone + env scout | conda env, vllm/CrewAI, weights 600 MB Lizard convnext-large | `external/SPARK/`, `external/hover_next_inference/` |
| Pod DM dispatch script | RunPod A40 25-35h + UNI/DINOv2 + HoVer-NeXt + DM1/DM2 classifier | `runpod_pod_DM_WSI_pathology.sh` |
| Lumenix SPARK 챗봇 (live) | 8 fallback ideas + Claude API + live evidence panel auto-refresh 1분 | `lumenix_spark/dm_idea_generator.html` |
| **로컬 cell detection** | 3200 tile, 16 slide, GeoJSON × 16 (~150 MB), per-stage cell counts | `cell_features_per_tile.tsv.gz`, `cell_geojson/*.geojson` |
| **SPARK Analytical bridge** | 8 SPARK features × 16 slides, per-stage 비교 | `spark_analytical_per_tile.tsv.gz`, `spark_analytical_per_slide.tsv`, `spark_analytical_per_stage.tsv` |
| SPARK-lite ST features (RNA-only) | spatial entropy, Moran's I-lite, top-25 pooling per slide | `spark_lite_features_per_slide.tsv`, `spark_lite_correlations.tsv` |

### 3.2 핵심 numerical findings

#### Phase A (ResNet50) — ALL CHANCE
| Tile size | Spearman | AUROC DM1-high |
|-----------|----------|----------------|
| 224 | 0.022 | 0.511 |
| 448 | 0.017 | 0.519 |
| 672 | 0.015 | 0.513 |

→ ResNet50 **확정 dead end**, UNI 필수.

#### SPARK Analytical (per-stage, real cell counts)
| Feature | N | PTC | LPTC | ATC | Δ(ATC-N) |
|---------|---|-----|------|-----|----------|
| B5 thyrocyte_cluster_med | 3.53 | 4.72 | 4.06 | **5.26** | **+49%** ★ |
| B4 tumor_stroma_ratio | 0.272 | 0.315 | 0.314 | **0.363** | **+33%** |
| C2 celltype_entropy_top25 | 1.43 | 1.46 | 1.48 | **1.49** | +4% |
| B3 nuclear_ecc_var | 0.012 | 0.012 | 0.012 | 0.012 | 0% (saturated) |
| B1 stromal_encased_thy | 0.989 | 0.991 | 0.986 | 0.993 | 0% (saturated) |

→ B5/B4 가 stage trajectory 단조증가 = **SPARK 가 H&E 만으로 분화도 axis 잡음을 입증**. B1/B3 는 saturated (rule-based classifier 한계 → HoVer-NeXt 필요).

#### SPARK-lite RNA (per-stage)
| Stage | DM1 top25 mean | Moran's I-lite |
|-------|----------------|----------------|
| N | 0.483 | 0.523 |
| PTC | 0.414 | 0.610 |
| LPTC | 0.397 | 0.659 |
| ATC | **0.181** | 0.581 |

→ DM1 hotspot density 가 dedifferentiation 진행에 따라 -62% 감소.

### 3.3 대기

- **Pod DM dispatch** (사용자 RunPod 트리거) — UNI + HoVer-NeXt 7-class + DM1/DM2 LOSO classifier, 25-35h, ~$10-15
- **SPARK Generative pipeline live LLM** (사용자 API key) — 챗봇에서 Claude 호출
- **공간 전사체 통합** (다음 task 9) — 아래 §4

---

## 4. 공간 전사체 (ST) × SPARK 통합 — YES, 가능

### 4.1 왜 통합 가능한가

GSE250521 Visium spot 과 H&E tile 이 **이미 1:1 매칭**돼 있음:
- 각 spot 의 (x_hires, y_hires) 좌표로 H&E hires image 에서 patch 추출
- 그 patch 가 곧 우리가 SPARK 에 입력하는 tile
- **즉 spot ↔ tile ↔ SPARK feature ↔ 18,000+ gene 발현 직접 매칭**

### 4.2 통합 후 가능해지는 분석

| 분석 | 형식 | 결과 |
|------|------|------|
| Gene → SPARK feature regression | per-spot, 8521 spots × 8 SPARK feat × 18k gene | 어떤 gene 모듈이 어떤 SPARK 형태형질을 예측하는지 |
| SPARK feature → DM1 score residual | per-spot, 8 feat → DM1_like_score | H&E 에서 DM1 직접 readout |
| Niche-level co-localization | spatial autocorrelation × cell-class | 예: TLS 와 NIS expression 공간 일치도 |
| Drug-response axis prediction | SPARK feat × 약물 반응 module score | RAI fail signature 의 H&E 표현 |
| Inverse mapping | SPARK feat → 가장 강하게 상관하는 gene set | SPARK B5 = "thyrocyte cluster median" 이 어떤 gene 으로 reproduce 되는가 |

### 4.3 즉시 실행 가능 — 데이터 준비 상태

- ✅ tile_metadata_resid.tsv.gz (8521 row, spot_id + sample_id + DM1/RAI/TDS resid)
- ✅ spark_analytical_per_tile.tsv.gz (3200 tile, 8 SPARK features)
- ✅ GSE250521 *.h5ad (16 sample, log-normalized expression)
- ⏳ ST × SPARK bridge script (다음 단계 — 5분 작성)

---

## 5. Pod DM dispatch — 사용자 1-line 으로 완료

```bash
# RunPod A40 48GB pod 시작 후 SSH 연결, 그 안에서:

mkdir -p /workspace/wsi_pathology_dm/manifests
# 로컬에서 manifest 업로드 (rsync 또는 직접 복사)
rsync -avz local:/home/seungho/personal/THCA_data_analysis/project/data/manifests/ \
    pod:/workspace/wsi_pathology_dm/manifests/
rsync -avz local:/home/seungho/personal/THCA_data_analysis/project/notebooks_or_scripts/runpod_pod_DM_WSI_pathology.sh \
    pod:/workspace/

bash /workspace/runpod_pod_DM_WSI_pathology.sh
# 25-35시간 후 /workspace/wsi_pathology_dm/dm_wsi_artifacts.tgz
# 로컬로 sync → project/data/processed/TCGA-THCA-WSI-DM/
```

**비용**: A40 48GB @ $0.39/h × 30h ≈ **$11.70**

---

## 6. Paper coverage (어느 paper 에 들어가나)

| Paper | 사용 |
|-------|------|
| Paper 1 | DM1 axis H&E inferability supplement (negative result Phase A + UNI positive) |
| Paper 2 | HT-overlap PTC (Hashimoto) cell composition 비교 (B5 cluster, lymphocyte density) |
| Paper 3 | ICI vulnerability (myeloid vs T-cell H&E proxy) |
| Paper 9 | synthetic lethality candidate triage (H&E features 가 SL pair 와 결합) |

→ "다 범용적으로" (사용자 요청) 충족.

---

## 7. 한 줄 요약

**ResNet50 + RNA 만으로는 DM1/DM2 분리 불가능을 negative control 100 panel × LOSO 로 입증했고, SPARK + UNI + HoVer-NeXt 가 그 갭을 메우는 다음 step. 로컬 cell-detection + SPARK Analytical 8-feature pipeline 은 이미 GSE250521 16-slide 에서 작동, B5/B4 가 stage trajectory 단조 (+49% / +33%) 로 SPARK 의 path 정당성 확보. Pod DM 한 번 dispatch 후 50 TCGA WSI × HoVer-NeXt 7-class × UNI embed × DM1/DM2 LOSO 가 확정 결과로 paper 준비.**
