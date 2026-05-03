# Phase A verdict — H&E → DM1_like_score_resid (ResNet50 LOSO)

**Run:** 2026-05-04 · **Owner:** Seungho Cook · **Path:** `project/reports/pathology_dm1_phaseA_cpu_verdict_2026_05_04.md`
**Verdict: NO-GO**

> ResNet50 ImageNet embedding cannot recover the depth-residualized DM1_like signal from 224 px H&E tiles in GSE250521. Pooled LOSO Spearman = 0.022, AUROC = 0.511 (chance ≈ 0.5). All 16 per-fold Spearman r ∈ [−0.07, +0.12]; no single slide reaches the GO threshold. R² strongly negative confirms predictions are worse than the per-fold mean.

---

## [1] Device

| field | value |
|---|---|
| machine | YG1-cook-vm (Standard_D8as_v5, koreacentral) |
| GPU | none — Azure subscription has no GPU quota in any region |
| torch | 2.4.0+cpu |
| torchvision | 0.19.0+cpu |
| Phase A wall time so far | embed 4m 8s + LOSO 6m 55s + plots <2s ≈ **11 min CPU** |

GPU rationale: both 전산신약개발연구실 subscriptions checked across 25+ regions show GPU family quota = 0/0. Free Support plan blocks programmatic ticket creation. Decision: run on CPU since ResNet50 + 3,200 tiles is feasible (~10 min embed) and the scientific answer is independent of device.

## [2] Tile counts

| size | n tiles | n samples | tiles/sample |
|---|---|---|---|
| 224 px | 3,200 | 16 | 200 (balanced) |

Stage balance: 800 PT (4 slides) / 800 PTC (4) / 800 LPTC (4) / 800 ATC (4). Multi-resolution (448 / 672) deferred — Phase A NO-GO at 224 px makes upscaling unlikely to bridge the gap.

## [3] Target labels

| target | available | n non-NaN | within-sample residualization |
|---|---|---|---|
| DM1_like_score_resid | ✓ | 3,200 | OLS on log_counts + log_ngenes |
| RAI_8_score_resid | ✓ | 3,200 | same |
| TDS_like_score_resid | ✓ | 3,200 | same |

Note: by construction `DM1_like_score = -RAI_8_score` (the score generator inverts RAI z-score to define DM1-like state). Hence the two targets give identical Ridge regression coefficients up to sign — observed in metrics below.

## [4] LOSO metrics (16-fold, 224 px, ResNet50 ImageNet)

| target | model | pooled Spearman r | pooled Pearson r | pooled R² | pooled AUROC (DM1_high global q75) | n |
|---|---|---|---|---|---|---|
| **DM1_like_score_resid** | **Ridge** | **0.022** | −0.006 | −9.98 | **0.511** | 3,200 |
| DM1_like_score_resid | ElasticNet | 0.017 | +0.011 | −1.04 | 0.500 | 3,200 |
| RAI_8_score_resid | Ridge | 0.022 | −0.006 | −9.98 | 0.519 | 3,200 |
| RAI_8_score_resid | ElasticNet | 0.017 | +0.011 | −1.04 | 0.516 | 3,200 |
| TDS_like_score_resid | Ridge | 0.041 | +0.003 | −9.54 | 0.526 | 3,200 |
| TDS_like_score_resid | ElasticNet | 0.040 | +0.032 | −0.97 | 0.524 | 3,200 |

### Per-fold detail (DM1_like_score_resid, Ridge)

| held-out slide | Spearman r | AUROC (global q75) |
|---|---|---|
| GSM7980860_N-1 | +0.005 | 0.520 |
| GSM7980861_N-2 | −0.054 | 0.479 |
| GSM7980862_N-3 | −0.019 | 0.491 |
| GSM7980863_N-4 | +0.015 | 0.545 |
| GSM7980864_PTC-1 | **+0.118** | **0.583** |
| GSM7980865_PTC-2 | +0.035 | 0.478 |
| GSM7980866_PTC-3 | +0.050 | 0.542 |
| GSM7980867_PTC-4 | +0.118 | 0.520 |
| GSM7980868_LPTC-1 | +0.055 | 0.500 |
| GSM7980869_LPTC-2 | +0.035 | 0.504 |
| GSM7980870_LPTC-3 | −0.044 | 0.464 |
| GSM7980871_LPTC-4 | −0.024 | 0.481 |
| GSM7980872_ATC-1 | +0.081 | NaN (single class) |
| GSM7980873_ATC-2 | +0.002 | NaN (single class) |
| GSM7980874_ATC-3 | −0.074 | NaN (single class) |
| GSM7980875_ATC-4 | +0.044 | 0.501 |
| **median** | **+0.025** | **+0.501** |
| **range** | [−0.074, +0.118] | [0.464, 0.583] |

ATC AUROC = NaN: every ATC tile sits above the train-defined DM1_high threshold (single-class held-out fold). Expected — ATC is the most de-differentiated stage and its DM1 axis distribution barely overlaps PT/PTC/LPTC.

## [5] Best model

- **Ridge alpha=1.0**, tile_size=224, ResNet50 ImageNet1K_V2, 224 px center crop, ImageNet normalization
- ElasticNet (alpha=0.001, l1_ratio=0.5) is statistically indistinguishable
- No tile size reached the GO gate, so multi-resolution selection bias is moot

## [6] Negative control comparison

> n_random=30 panels in progress (background, ETA ~15 min on CPU). Will be appended when complete.

Preliminary expectation: random 8-gene panels and housekeeping panel should be near-zero, **confirming** that the DM1_like signal absence is not a pipeline bug but a real failure of ResNet50 ImageNet embeddings to encode the relevant histological features. If, instead, a control accidentally beats DM1_like, that flags either a label-leak or a pre-processing artifact.

| panel | Spearman r | AUROC top25 |
|---|---|---|
| **DM1_like (real)** | **0.022** | **0.511** |
| housekeeping (8 genes) | _pending_ | _pending_ |
| random 8-gene (n=30, mean) | _pending_ | _pending_ |
| random 8-gene (max) | _pending_ | — |

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

No leakage found. The near-zero signal is genuine.

## [8] **Verdict: NO-GO**

| criterion | threshold | observed | met? |
|---|---|---|---|
| pooled Spearman r | ≥ 0.30 (GO) | **0.022** | ✗ |
| pooled AUROC DM1_high | ≥ 0.70 (GO) | **0.511** | ✗ |
| Spearman ≥ 0.20 OR AUROC ≥ 0.65 | (BORDERLINE) | 0.022 / 0.511 | ✗ |
| Spearman < 0.20 AND AUROC < 0.65 | (NO-GO) | both true | ✓ |

**Hard NO-GO** under the spec gate. Both BORDERLINE thresholds also missed by ~10×.

### Why this likely failed (ranked)

1. **224 px tile = ~110 µm at default Visium hires scalefactor**. ResNet50 ImageNet weights see tissue architecture (gland shape, stromal density), not the nuclear cytology / chromatin / immune-cell contact patterns that would correlate with thyroid de-differentiation. Foundation models (UNI, CONCH, Virchow2) are trained on much higher resolution H&E and might recover some — but going from r = 0.02 to r = 0.30 is a **15× jump**; published UNI vs ResNet50 gains on tile-level tasks are typically **2–5×**, not 15×.
2. **Depth-residualized target by design strips most thyroid lineage biology**. RAI_8 high spots are also high-count spots (active follicles transcribe more). Per-sample OLS on log_counts + log_ngenes removes the dominant axis. What remains is a residual that is small in magnitude (and noisy) — even a perfect predictor of the underlying biology would land in low Spearman r if the residual is mostly noise.
3. **Visium hires PNG ≠ full-res scanner image**. Hires is downsampled by `tissue_hires_scalef` (default 0.17). True nuclear morphometry lives in fullres which is not present in this dataset's GEO submission. This is a hard ceiling on tile-level resolution.
4. **n = 16 slides, all from one cohort, single sectioning protocol**. Even with 200 tiles/slide, the LOSO denominator is 16. Per-fold Spearman variance is high (best 0.118, worst −0.074) — but the magnitude is uniformly small, so noise alone cannot rescue this.
5. **DM1 residual ≠ DM1 phenotype boundary**. The label is a continuous program score, not a binary segmentation of "DM1 zone" within tissue. Spot-level DM1 high vs low is a morphology-blind bulk readout — H&E may simply not contain the information.

### Implication for Paper 2 / IP

| question | answer |
|---|---|
| Does H&E predict the depth-residualized 8-gene DM1 program at tile level? | **No (with ResNet50 ImageNet at 224 px)** |
| Should we proceed to Phase B (multi-cohort LODO)? | **No** — fixing 0.02 → 0.30 with more cohorts is implausible |
| Should we proceed to Phase C (TCGA WSI external)? | **No** — same model would be applied; same failure |
| Is foundation-model upgrade (UNI/CONCH/Virchow2) worth pursuing? | **Marginal** — see "Conditional re-test" below; will not bridge a 15× gap |
| Is "image-DM1 triage" a viable Paper 2 / IP angle? | **No, as currently scoped.** |

### Conditional re-test (only if user explicitly authorizes)

If the user wants to be **certain** before closing the concept entirely, three small follow-ups can be done locally on CPU within 1 hour total:

1. **Foundation embedder (CONCH, smaller than UNI)** — apply for HF gated access (kukshomr@gmail.com), embed the same 3,200 tiles. CPU runtime ~30 min. If r jumps to ≥ 0.10, escalate to UNI; otherwise, confirmed NO-GO.
2. **Raw (non-residualized) DM1 target sensitivity** — break the spec rule for one diagnostic to see whether any signal lives in the part of biology we residualized away. If r → 0.4+ on raw target and r → 0.02 on residualized, the conclusion is "model can predict thyroid lineage, not the residualized DM1 axis" — informative but not the original concept.
3. **Larger context tile (672 px)** — extract + embed (~30 min CPU). Tests hypothesis 1 above; tile-architecture vs nuclear-cytology question.

These are *not* a path to GO. They're a closure document.

## [9] Next command

```bash
# No Phase B. No Phase C. Stop.
#
# If user authorizes the conditional re-tests above, run:
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
  embeddings_resnet50_224.npz                 11.6 MB
  loso_metrics_resnet50.tsv                   per-fold metrics
  loso_predictions_resnet50.tsv.gz            per-tile predictions (3,200 × 6 model/target combos)
  pred_vs_obs_resnet50.png                    pooled scatter
  pred_vs_obs_resnet50_per_slide.png          16-slide bar chart
  negative_controls_summary.tsv               (in progress)
project/reports/
  pathology_dm1_phaseA_cpu_verdict_2026_05_04.md   ← this file
```

---

*Generated 2026-05-04 by automated Phase A pipeline. ResNet50 ImageNet on 3,200 GSE250521 tiles, 16-fold LOSO. No Azure GPU quota; CPU run. Final wall time including this report: ~12 min.*
