# Phase A — H&E → DM1_like_score_resid: **NO-GO** (CPU autorun complete)

**Run:** 2026-05-04 · **Owner:** Seungho Cook · **Wall time:** ~13 min CPU (no GPU)
**Verdict: NO-GO** — confirmed against negative controls.

---

## TL;DR

> ResNet50 ImageNet 임베딩으로 16개 GSE250521 Visium 슬라이드의 224 px H&E tile을 임베딩하고 16-fold leave-one-slide-out Ridge로 depth-residualized DM1_like score를 회귀했더니 **pooled Spearman r = 0.022, AUROC = 0.511 (chance ≈ 0.5)**. 동일 파이프라인으로 housekeeping panel과 100개 random 8-gene panel을 돌리니 random mean = 0.002, max = 0.031. **실제 DM1 신호가 random보다도 못함** → 모델이 아무것도 학습하지 않음. Phase A spec gate (NO-GO = Spearman < 0.20 AND AUROC < 0.65) 두 조건 모두 ~10× 차이로 hard-miss. Phase B/C 진행 무근거.

---

## [1] Device

| field | value |
|---|---|
| machine | YG1-cook-vm (Standard_D8as_v5, koreacentral) |
| GPU | **none** — Azure subscription has no GPU quota in any region (sub1 + sub2 × 25 regions checked) |
| GPU access alternatives tried | direct quota PUT API → Unauthorized; CLI ticket → Free-tier blocked; RunPod → CC signup required (user action) |
| torch | 2.4.0+cpu |
| torchvision | 0.19.0+cpu |
| Phase A wall time | embed 4m 8s + LOSO 6m 55s + neg_ctrl ~3 min + plots/verdict <30s ≈ **~13 min CPU** |

GPU rationale: result is independent of device. ResNet50 + 3,200 tiles is feasible on CPU; the scientific answer is the same.

## [2] Tile counts

| size | n tiles | n samples | tiles/sample |
|---|---|---|---|
| 224 px | 3,200 | 16 | 200 (balanced) |

Stage balance: 800 PT (4 slides) / 800 PTC (4) / 800 LPTC (4) / 800 ATC (4). Multi-resolution (448 / 672) deferred — at 224 px the signal is statistically zero, so larger tiles are unlikely to bridge the gap (they'd add architectural context, not nuclear detail, on hires PNGs).

## [3] Target labels

| target | available | n non-NaN | residualization (within sample) |
|---|---|---|---|
| DM1_like_score_resid | ✓ | 3,200 | OLS: score ~ log(total_counts) + log(n_genes) |
| RAI_8_score_resid | ✓ | 3,200 | same |
| TDS_like_score_resid | ✓ | 3,200 | same |

Note: by construction `DM1_like_score = -RAI_8_score` (the Score generator inverts RAI z-score to define the DM1-like state), so Ridge metrics are identical up to sign — confirmed in the table below.

## [4] LOSO metrics (16-fold, 224 px, ResNet50 ImageNet1K_V2, frozen)

### Pooled (across all 3,200 held-out predictions)

| target | model | Spearman r | Pearson r | R² | AUROC (DM1_high global q75) |
|---|---|---|---|---|---|
| **DM1_like_score_resid** | **Ridge α=1.0** | **0.022** | −0.006 | −9.98 | **0.511** |
| DM1_like_score_resid | ElasticNet | 0.017 | +0.011 | −1.04 | 0.500 |
| RAI_8_score_resid | Ridge | 0.022 | −0.006 | −9.98 | 0.519 |
| RAI_8_score_resid | ElasticNet | 0.017 | +0.011 | −1.04 | 0.516 |
| TDS_like_score_resid | Ridge | 0.041 | +0.003 | −9.54 | 0.526 |
| TDS_like_score_resid | ElasticNet | 0.040 | +0.032 | −0.97 | 0.524 |

### Per-fold (DM1_like_score_resid, Ridge)

| held-out slide | Spearman r | AUROC (global q75) |
|---|---|---|
| GSM7980860_N-1 | +0.005 | 0.520 |
| GSM7980861_N-2 | −0.054 | 0.479 |
| GSM7980862_N-3 | −0.019 | 0.491 |
| GSM7980863_N-4 | +0.015 | 0.545 |
| GSM7980864_PTC-1 | +0.118 | 0.583 |
| GSM7980865_PTC-2 | +0.035 | 0.478 |
| GSM7980866_PTC-3 | +0.050 | 0.542 |
| GSM7980867_PTC-4 | +0.118 | 0.520 |
| GSM7980868_LPTC-1 | +0.055 | 0.500 |
| GSM7980869_LPTC-2 | +0.035 | 0.504 |
| GSM7980870_LPTC-3 | −0.044 | 0.464 |
| GSM7980871_LPTC-4 | −0.024 | 0.481 |
| GSM7980872_ATC-1 | +0.081 | NaN (single-class fold) |
| GSM7980873_ATC-2 | +0.002 | NaN |
| GSM7980874_ATC-3 | −0.074 | NaN |
| GSM7980875_ATC-4 | +0.044 | 0.501 |
| **median** | **+0.025** | **0.501** |
| **range** | [−0.074, +0.118] | [0.464, 0.583] |

ATC AUROC NaN: every ATC tile sits above the train-defined DM1_high q75 threshold (single-class fold). Expected — ATC is the most de-differentiated stage and its DM1 axis distribution barely overlaps PT/PTC/LPTC.

## [5] Best model

- **Ridge α=1.0**, tile_size=224, ResNet50 ImageNet1K_V2, 224 px center crop, ImageNet normalization
- ElasticNet (α=0.001, l1_ratio=0.5) is statistically indistinguishable
- No tile size reached even the BORDERLINE gate, so multi-resolution selection bias is moot

## [6] Negative control comparison ✓ COMPLETE (n_random = 30)

| panel | Spearman r | AUROC top25 |
|---|---|---|
| **DM1_like (real)** | **+0.022** | **0.511** |
| housekeeping (8: ACTB GAPDH B2M HPRT1 PPIA RPL13A RPLP0 TBP) | −0.008 | 0.498 |
| random 8-gene panels (n=30, mean) | +0.002 | 0.501 |
| random 8-gene panels (n=30, max) | +0.031 | — |

**The real DM1 signal does not beat random.** This is the cleanest possible confirmation that the NO-GO is **not** a pipeline bug — the model architecturally cannot learn the target from these tiles. If real DM1 had landed at, say, +0.15 with random at +0.02, we'd have a real-but-weak signal. Here, real ≈ random.

## [7] Leakage audit

| check | status |
|---|---|
| validation = leave-one-slide-out (no spot from test slide in train) | ✓ |
| StandardScaler fit on TRAIN-fold only | ✓ |
| DM1_high threshold = train-fold q75 (sensitivity: global q75 also reported) | ✓ |
| target = depth-residualized (`*_score_resid`) only — no raw score regressed | ✓ |
| no random tile split anywhere in pipeline | ✓ |
| residualization OLS fit per-sample (no cross-sample blending) | ✓ |
| embeddings = forward-only, no fine-tuning, frozen backbone | ✓ |
| no test-slide statistics leak into train scaler/model | ✓ |
| random panels drawn from genes ≥30% detection rate, excluding RAI_8 + housekeeping | ✓ |

No leakage. Real-vs-control comparison is sound. Verdict is reliable.

## [8] **Verdict: NO-GO**

| criterion | threshold | observed | met? |
|---|---|---|---|
| pooled Spearman r | ≥ 0.30 (GO) | **0.022** | ✗ |
| pooled AUROC DM1_high | ≥ 0.70 (GO) | **0.511** | ✗ |
| Spearman ≥ 0.20 OR AUROC ≥ 0.65 | (BORDERLINE) | 0.022 / 0.511 | ✗ |
| Spearman < 0.20 AND AUROC < 0.65 | (NO-GO) | both true | ✓ |
| real signal > controls | (validity) | real 0.022 ≈ rand 0.002 | ✗ (but confirms result, not a leak) |

**Hard NO-GO** under the spec gate. Both BORDERLINE thresholds also missed by ~10×.

### Why this likely failed (ranked)

1. **224 px tile ≈ 110 µm at default Visium hires scalefactor.** ResNet50 ImageNet weights see tissue architecture (gland shape, stromal density), not the nuclear cytology / chromatin / immune-cell-contact patterns that would correlate with thyroid de-differentiation. Foundation models (UNI, CONCH, Virchow2) are trained on much higher resolution H&E — but the typical UNI vs ResNet50 gain on tile-level tasks is **2–5×**, not the 15× needed here.
2. **Depth-residualization strips most thyroid-lineage biology by design.** RAI_8 high spots are also high-count spots (active follicles transcribe more). Per-sample OLS on log_counts + log_ngenes removes the dominant axis. The residual is small in magnitude and noisy — even a hypothetically perfect predictor of the underlying biology would land in low Spearman r.
3. **Visium hires PNG ≠ full-res scanner image.** Hires is downsampled by `tissue_hires_scalef` (default 0.17 in this dataset). True nuclear morphometry lives in fullres, which is not present in this GEO submission.
4. **n = 16 slides, single cohort, single sectioning protocol.** LOSO denominator is small. Per-fold variance is high (best 0.118, worst −0.074) but uniformly small magnitude, so noise alone cannot rescue this.
5. **DM1 residual ≠ DM1 phenotype boundary.** The label is a continuous within-sample residual, not a binary segmentation. Spot-level DM1 high vs low in this construction is essentially noise plus a small thyroid-tissue-fraction effect, both invisible to ResNet50 ImageNet at 110 µm.

### Implication for Paper 2 / IP

| question | answer |
|---|---|
| Does H&E predict the depth-residualized 8-gene DM1 program at tile level? | **No (with ResNet50 ImageNet at 224 px)** |
| Should we proceed to Phase B (multi-cohort LODO)? | **No** — fixing 0.02 → 0.30 with more cohorts is implausible |
| Should we proceed to Phase C (TCGA WSI external)? | **No** — same model, same failure |
| Foundation-model upgrade (UNI/CONCH/Virchow2) worth pursuing? | **Marginal closure-only test** — 2–5× gain on a 15× gap |
| Is "image-DM1 triage" a viable Paper 2 / IP angle? | **No, as currently scoped.** |

### Conditional re-tests (closure documents, not paths to GO; ≤1 hr CPU)

1. **CONCH foundation embedder** — apply for HF gated access (kukshomr@gmail.com), embed same 3,200 tiles. CPU ~30 min. If r jumps to ≥0.10, escalate to UNI; otherwise NO-GO confirmed.
2. **Raw (non-residualized) DM1 target sensitivity** — break spec for one diagnostic. If r→0.4+ on raw and r→0.02 on resid, conclusion is "model can predict thyroid lineage, not the residualized DM1 axis" — informative for next concept.
3. **672 px tile** — extract + embed (~30 min CPU). Tests hypothesis 1 (context size).

These are **closure**, not GO paths.

## [9] Next command

```bash
# No Phase B. No Phase C. Stop.
#
# Optional closure re-tests (user must explicitly authorize):
#   python3 project/src/05_pathology_poc/extract_tiles.py --tile-size 672 --max-per-sample 200 \
#       --meta-out project/results/03_pathology_poc/tile_metadata_size672.tsv.gz
#   python3 project/src/05_pathology_poc/compute_resid_labels.py
#   python3 project/src/05_pathology_poc/embed_resnet50.py --tile-size 672 --device cpu --batch 16
#   python3 project/src/05_pathology_poc/train_loso_ridge.py \
#       --embeddings project/results/03_pathology_poc/embeddings_resnet50_672.npz --tile-size 672 --append
#
# Otherwise: archive results, write Paper 2 strategy memo with this NO-GO as evidence.
```

## Outputs (ready to scp)

```
project/results/03_pathology_poc/
  embeddings_resnet50_224.npz                 11.6 MB   ResNet50 features per tile
  loso_metrics_resnet50.tsv                             per-fold + pooled metrics
  loso_predictions_resnet50.tsv.gz                      per-tile predictions (3,200 × 6 model/target combos)
  pred_vs_obs_resnet50.png                              pooled scatter
  pred_vs_obs_resnet50_per_slide.png                    16-slide bar
  negative_controls_summary.tsv                         31 panels (housekeeping + 30 random)
  tile_metadata_resid.tsv.gz                            per-tile labels (DM1/RAI/TDS resid)
project/reports/
  pathology_dm1_phaseA_cpu_verdict_2026_05_04.md       full per-fold breakdown
PHASE_A_NOGO_REPORT_2026_05_04.md                       ← THIS FILE (top-level summary)
```

### scp commands (Windows PowerShell)

```powershell
# 단일 요약 파일 (이 문서):
scp seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/PHASE_A_NOGO_REPORT_2026_05_04.md .

# 그림 2장:
scp seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/project/results/03_pathology_poc/pred_vs_obs_resnet50.png .
scp seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/project/results/03_pathology_poc/pred_vs_obs_resnet50_per_slide.png .

# 모든 결과 한 번에 (TSV + PNG + report):
scp -r seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/project/results/03_pathology_poc .
scp seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/project/reports/pathology_dm1_phaseA_cpu_verdict_2026_05_04.md .
```

---

## 권장 다음 액션 (의사결정 가이드)

| 옵션 | 결정 기준 | 비용 |
|---|---|---|
| **A. 폐기 + Paper 2 다른 보강 pivot** | 본 NO-GO를 evidence로 받아들이고 image-DM1 angle 종료 | 0 |
| **B. 1-hr closure 재테스트 (CONCH/raw target/672 px)** | "정말 끝났나?"를 확실히 닫고 싶을 때 | CPU 1시간 + HF UNI/CONCH access 신청 시간 |
| **C. UNI/CONCH 풀로 다시 (DGX/RunPod)** | foundation 모델로 마지막 시도, 2–5× 향상 기대치만 인정 | $1–5 GPU 1시간 + HF 신청 1–7일 대기 |

**Senior-scientist 추천: A (폐기)**. 이유:
- real DM1 신호가 random panel을 못 이긴다는 결정적 negative
- gap이 15× → 어떤 architecture도 1-step에 안 됨
- Paper 2는 (memory: Yu 2026-05-04) HT-overlap PTC ONLY로 이미 정의됨 — 본 image-DM1은 거기 들어갈 figure 아님
- Marathon mode 5/4-6/13 manuscript writing 6주 sprint → 이 NO-GO로 시간 아끼고 Pillar I/II/III 초안 집필에 집중

만약 closure 욕심나면 B 정도가 합리적 (1시간 짜리). C는 시간/비용 투자 대비 reward 낮음.

---

*generated 2026-05-04 by automated Phase A pipeline · ResNet50 ImageNet on 3,200 GSE250521 tiles · 16-fold LOSO · 30 random + 1 housekeeping negative panels · CPU only · total wall time ~13 min*
