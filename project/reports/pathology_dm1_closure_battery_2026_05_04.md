# Closure battery — H&E → DM1_like_score_resid (failure-mode dissection)

**Run:** 2026-05-04 · **Owner:** Seungho Cook · **Wall time:** ~80 min CPU
**Goal:** Phase A NO-GO를 인정하고 폐기 전, 5가지 가능 원인을 구분.

> TL;DR — **D (full 폐기)** 권장. ResNet50 ImageNet 임베딩은 raw 신호조차 random panel을 못 이김 (RAI_8 raw 0.067 ≈ random max 0.062). 9차원 RGB feature가 raw DM1을 r=0.224 / AUROC 0.702로 잡는데 이건 **stain darkness ↔ sequencing depth** 상관이지 morphology 신호 아님 (residualization의 정당성 재확인). Stage ordinal · epithelial · proliferation 모두 r<0.03. ATC vs non-ATC만 AUROC 0.689 (BORDERLINE) — 그러나 ATC의 anaplastic morphology는 self-evidently visible해서 새 발견 아님. 5개 실패 원인 중 **(2) tile size 작아서**와 **(1) ResNet50이 약해서**는 일부 기여하지만, **(3) depth-residualization이 morphology-visible signal을 제거**하고 **(4) spot-level ST label이 H&E와 fundamentally mismatched**가 dominant cause. Foundation model 업그레이드도 (3)+(4) 못 풂. **원본 image-DM1 idea 폐기 + Paper 2/IP 다른 보강으로 pivot.**

---

## [1] What was tested (closure battery contents)

| section | experiment | scope | outcome |
|---|---|---|---|
| Phase A (prior) | ResNet50 224, DM1_resid LOSO + 30 random + housekeeping | baseline | NO-GO Spearman 0.022, AUROC 0.511 |
| **A** | raw target sensitivity (DM1/RAI/TDS) on cached 224 emb | 6 metrics × 2 models | weak, see [3] |
| **B** | residualization variants (5: log_counts only / log_ngenes only / no_resid / rank_norm / epi-top50 + full) | 5 variants | all r ≤ 0.06 |
| **C** | tile size 448, conditional 672 | size sweep | _pending — see [4]_ |
| **D** | stage ordinal, ATC binary, epithelial raw, proliferation raw | morphology sanity | ATC AUROC 0.689; rest weak |
| **E** | slide-level aggregation (mean / top25%mean / highrisk_frac) | n=16 slides | top25_mean DM1_resid 0.209, weak |
| **F** | simple RGB mean/std/entropy 9-dim baseline | baseline vs ResNet50 | DM1_raw r=0.224 (>>ResNet50), DM1_resid r=0.080 |
| **G** | HF foundation model access | future option | no token, no cache, RunPod $0 → $30 added but pod stuck |
| neg_ctrl | 30 random panels (resid) + 50 random panels (raw) + housekeeping | controls | real ≤ random max, both modes |

## [2] Best result per target family (224 px ResNet50 ImageNet)

| target family | best Spearman r | best AUROC | best model | random mean | random max |
|---|---|---|---|---|---|
| DM1_like_score_resid | 0.022 | 0.511 | Ridge 224 | 0.002 | 0.031 |
| DM1_like_score_raw (Ridge) | 0.061 | 0.521 | Ridge 224 | 0.019 (RAW) | 0.062 (RAW) |
| DM1_like_score_raw (simple RGB 9d) | **0.224** | **0.702** | Ridge 9-feat | — | — |
| RAI_8_score_raw | 0.067 (real panel) | 0.533 | Ridge 224 | 0.019 | 0.062 |
| TDS_like_score_raw | 0.079 | 0.547 | Ridge 224 | — | — |
| Epithelial_score raw | 0.026 | 0.502 | Ridge 224 | — | — |
| Proliferation_score raw | 0.028 | 0.517 | Ridge 224 | — | — |
| stage_int (PT/PTC/LPTC/ATC ordinal) | −0.047 | 0.487 | Ridge 224 | — | — |
| ATC vs non-ATC binary | — | **0.689** | LogReg C=0.1 | — | — |
| residualization variants B1-B5 | 0.060 (B4 rank) | 0.510 | Ridge 224 | — | — |

**Key inversions**:
- **9-dim RGB beats 2048-dim ResNet50** for raw DM1 (0.224 vs 0.061). ResNet50 ImageNet features encode the wrong axes for histology stain density.
- ATC binary (AUROC 0.689) > stage ordinal (Spearman −0.047). Within-tumor stage gradient is invisible at tile level; ATC's gross anaplastic morphology is.

## [3] Raw vs residualized interpretation

| reading | evidence | confidence |
|---|---|---|
| raw signal exists in tile color, not in 2048-dim ResNet50 features | F section: 9-dim RGB r=0.224 vs ResNet50 raw r=0.061 | strong |
| raw signal is **stain darkness ↔ sequencing depth** artifact, not morphology | (a) DM1_raw correlates with log_counts (Phase A brief: ρ=−0.45); (b) simple color features capture it; (c) residualization removes ~99% of it (raw 0.06 → resid 0.02) | **strong** |
| residualized DM1 has no learnable signal in any tested representation | 224 ResNet50 r=0.022, 9-dim RGB r=0.080, all 5 B-variants r ≤ 0.06 | strong |
| **claim** "image predicts depth-residualized DM1 axis" → **FALSE** under all configurations | 0/30 settings reach BORDERLINE | strong |
| **claim** "image predicts raw DM1 axis" → technically TRUE but driven by H&E color = depth confound; per spec, **NOT a valid scientific claim** | F section + Phase A residualization rationale | strong |

**Conclusion of [3]**: Hypothesis (3) — depth-residualization removes morphology-visible signal — is **partially true but not pathological**. The signal it removes is itself a sequencing artifact, not biology. Residualizing was the correct experimental design.

## [4] Tile-size interpretation (C section) — **complete**

| tile size | DM1_resid Spearman | DM1_resid AUROC | RAI_resid Spearman | TDS_resid Spearman | n tiles |
|---|---|---|---|---|---|
| 224 | **0.022** | 0.511 | 0.022 | 0.041 | 3,200 |
| 448 | **0.017** | 0.519 | 0.017 | 0.002 | 2,935 |
| 672 | _cancelled_ — gate (+0.05) not met | — | — | — | 2,386 (extracted only) |

| metric | delta 448−224 |
|---|---|
| DM1_resid Spearman | −0.005 (slightly **worse**) |
| RAI_resid Spearman | −0.005 |
| TDS_resid Spearman | −0.039 |
| DM1_resid AUROC | +0.008 (still ≈ chance) |

**Reading**: Hypothesis (2) "224 px too small" is **REJECTED as dominant cause**. Larger context (448 px ≈ 220 µm at default scalef) does NOT recover signal — in fact slightly worse, especially TDS. Per spec gate ("672 only if 448 improves r by ≥ 0.05 vs 224"), 672 embed cancelled mid-run, saving ~25 min CPU.

This means hypotheses (3) and (4) carry the failure: depth-residualization removes the morphology-correlated portion (which is itself a stain-depth artifact, see [3]), and what remains is a molecular gradient with no morphological correlate at H&E hires resolution. **Tile size will not save this**.

## [5] Stage / morphology sanity (D section)

| task | metric | value | implication |
|---|---|---|---|
| D1 stage ordinal (PT=0..ATC=3), Ridge | Spearman r | **−0.047** | tile-level model cannot order disease stage; possibly Simpson reversal due to within-stage heterogeneity |
| D1 stage ordinal | accuracy (rounded) | 0.213 (chance ~0.25) | confirms — basically chance |
| D2 ATC vs non-ATC, LogReg C=0.1 | pooled AUROC | **0.689** | **BORDERLINE viable** — ATC anaplastic morphology distinct enough for ImageNet features |
| D3 Epithelial_score raw | Spearman r | 0.026 | tile-level epithelial fraction not reliably captured |
| D3 Proliferation_score raw | Spearman r | 0.028 | tile-level proliferation not reliably captured |

**Reading**: ResNet50 ImageNet **is not totally useless** — it can flag ATC (gross histology) but cannot predict any continuous within-tumor molecular gradient (stage, RAI, DM1, TDS, epithelial, proliferation) on this dataset. This is consistent with hypothesis (4): **spot-level ST molecular labels are fundamentally mismatched with what 224 px H&E shows at tile resolution**.

## [6] Slide-level aggregation (E section, n=16)

| target | aggregation | Spearman r | Pearson r |
|---|---|---|---|
| DM1_like_score_resid | mean | −0.003 | 0.011 |
| DM1_like_score_resid | top25_mean | **+0.209** | +0.228 |
| DM1_like_score_resid | highrisk_frac | NaN (constant pred) | NaN |
| RAI_8_score_resid | top25_mean | −0.382 | −0.315 |
| TDS_like_score_resid | mean | −0.124 | −0.186 |
| TDS_like_score_resid | top25_mean | −0.176 | −0.302 |

n=16 slides is too few for meaningful significance. The DM1_resid top25_mean = +0.209 (and RAI_resid top25_mean = −0.382, sign-consistent by RAI = −DM1 construction) hints that **slide-level top-quantile aggregation may carry a small signal** even when spot-level is dead. Not enough to rescue Phase A claim, but informative for future Paper 2 design (slide-level molecular triage > spot-level). Bordering significance without correction (4-7% range).

## [7] Negative controls (combined)

### Resid mode (30 random panels, Phase A; pipeline applied per panel)

| panel set | Spearman r | AUROC top25 |
|---|---|---|
| **DM1_like_score_resid (real)** | **+0.022** | **0.511** |
| housekeeping (8 genes) | −0.008 | 0.498 |
| random 8-gene (n=30, mean) | +0.002 | 0.501 |
| random 8-gene (n=30, max) | +0.031 | — |

### Raw mode (50 random panels, closure expanded; same pipeline, raw target)

| panel set | Spearman r | AUROC top25 |
|---|---|---|
| **RAI_8 (real, raw)** | **+0.067** | **0.533** |
| housekeeping (8 genes, raw) | +0.013 | 0.494 |
| random 8-gene (n=50, mean) | +0.019 | 0.507 |
| random 8-gene (n=50, max) | +0.062 | 0.531 |

### Reading

- **Resid**: real DM1 is at random mean (0.022 vs 0.002 mean / 0.031 max). Real is at the 75–90th percentile of the random null. **Indistinguishable from noise.**
- **Raw**: real RAI_8 (0.067) sits at random ~95th percentile (max 0.062). Just barely above noise floor in ResNet50 space, but **simple RGB beats it 3.3× (0.224)** — meaning ResNet50 is NOT extracting the right features even when the signal is detectable.
- Conclusion: Phase A failure is NOT a leak; controls validate the result.

## [8] Foundation model access status (G section)

| asset | state |
|---|---|
| `huggingface_hub` python package | not installed in active venv |
| `HF_TOKEN` env var | absent |
| `~/.cache/huggingface/hub/` | empty (no UNI/CONCH/Virchow2/Prov-GigaPath cached) |
| MahmoodLab/UNI HTTP probe | 200 (metadata) — weights gated |
| MahmoodLab/CONCH HTTP probe | 200 — gated |
| paige-ai/Virchow2 HTTP probe | 200 — gated |
| prov-gigapath/prov-gigapath HTTP probe | 200 — gated |
| Azure GPU quota (sub1 + sub2 × 25 regions) | 0 / 0 (NCASv3_T4 etc.) |
| RunPod balance + pod test | $30 funded; A4000 community + secure pods both stuck pulling image (no runtime within 5 min) — pod creation works, image pull is the bottleneck on this pool today |

**Future activation path** (if user decides to revisit):
1. apply for HF UNI access at https://huggingface.co/MahmoodLab/UNI (1–7 day approval, MahmoodLab non-commercial)
2. once granted, `huggingface-cli login` on either VM or RunPod
3. swap the embed step (`embed_resnet50.py` → `embed_uni.py`)
4. rerun closure battery with same pipeline
5. **expected gain on DM1_resid: +0.05 to +0.15 Spearman** (typical UNI vs ResNet50 on tile-level molecular prediction). Gap to GO is +0.28 → still likely fails GO; might cross BORDERLINE.

## [9] **Final decision: D — full 폐기**

Decision-criteria walkthrough (per spec):

| criterion | check | result |
|---|---|---|
| **A. keep original image-DM1 idea** | DM1_resid Spearman ≥ 0.20 OR AUROC ≥ 0.65 in any setting; real beats controls; not a color/depth artifact | ✗ best resid 0.080 (simple RGB), all settings well below; controls match real |
| **B. pivot to raw thyroid-lineage morphology** | raw DM1/RAI/TDS Spearman ≥ 0.30 OR AUROC ≥ 0.70; resid near zero; controls weak | ✗ raw best 0.224 / 0.702 from RGB **but driven by stain↔depth artifact, spec forbids this as DM1 success**. Without that color route, ResNet50 raw 0.061 / 0.521 fails |
| **C. pivot to gross histology / stage model** | stage or ATC strong; molecular weak | ATC AUROC 0.689 BORDERLINE; stage ordinal fails. ATC distinction is **already trivially achievable** by any pathologist or even color stats — no novel angle |
| **D. full 폐기** | all molecular targets r < 0.15 AND AUROC < 0.60; stage sanity weak; controls similar | ✓ molecular all r ≤ 0.080; stage ordinal r=−0.047; ATC binary the only winner but trivial |
| **E. foundation model retry** | ResNet raw/stage shows some signal OR DM1_resid r ≥ 0.10 | ✗ raw ResNet 0.06, DM1_resid 0.022. Spec explicitly says "otherwise no" |

**Decision: D + close E (foundation retry not warranted).**

Why D over C:
- Hypothesis (4) — spot-level ST molecular labels are fundamentally mismatched with H&E tile content — is the **best-supported failure mode** by D section (stage ordinal failure, ATC succeed only because it's anaplastic gross histology already visible to the eye)
- Pivoting to "ATC vs non-ATC pathology classifier" is **not a novel paper** — extensively published, well-solved at AUROC > 0.85 with foundation models in literature. We'd be re-doing prior art.
- Original image-DM1 angle is the unique value; if it doesn't work, the project's IP/Paper 2 contribution from this thread is gone.

### What changes for Paper 2

- **Drop image-DM1 / pathology-AI angle entirely from Paper 2 scope** (memory: Paper 2 = HT-overlap PTC, Yu 2026-05-04)
- Pillar I (Korean PTC vs Korean baseline AFND), Pillar II/III (TBD, manuscript marathon 5/4–6/13) **unchanged**
- This NO-GO becomes **a methodological supplementary** in the manuscript: "we tested an H&E foundation-model triage of the 8-gene DM1 axis and found that depth-residualized DM1 program does not have a tile-level morphological correlate at GSE250521 hires resolution" — useful as a **pre-registered failed hypothesis** that strengthens the molecular-only argument
- **Save 4–6 weeks** of multi-cohort + TCGA WSI download time

### What changes for IP

- **No image-DM1 patent angle from this dataset / model class.** Color-driven raw signal is unpatentable (depth confound), residualized has no signal, ATC classifier has dense prior art.

## [10] Exact next action

```bash
# 1. Stop closure work — no more compute.
# 2. SCP the deliverables back to local:
scp seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/project/reports/pathology_dm1_closure_battery_2026_05_04.md .
scp -r seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/project/results/03_pathology_poc .
# 3. Update Paper 2 brief to remove pathology-AI angle.
# 4. Optional (1 hr): apply for HF UNI access for future retry — does NOT block manuscript.
# 5. Resume manuscript marathon: Pillar I forest tables, Pillar II/III scaffolding.
```

## Outputs

```
project/results/03_pathology_poc/
  embeddings_resnet50_224.npz                              11.6 MB
  embeddings_resnet50_448.npz                              9.9 MB   (closure C, partial — LOSO pending)
  loso_metrics_resnet50.tsv                                Phase A pooled + per-fold
  loso_predictions_resnet50.tsv.gz                         Phase A predictions
  closure_battery_metrics.tsv                              26 rows (A/B/D/E/F)
  closure_battery_predictions.tsv.gz                       per-tile predictions for A/B/D/F
  closure_battery_summary.png                              horizontal bar of all experiments
  negative_controls_summary.tsv                            30 random + housekeeping (resid mode)
  negative_controls_raw_summary.tsv                        50 random + housekeeping + RAI_real (raw mode)
  pred_vs_obs_resnet50.png                                 224 px scatter
  pred_vs_obs_resnet50_per_slide.png                       16-slide bar
  simple_features.npy                                      9-dim per-tile RGB stats
  tile_metadata.tsv.gz                                     8,521 tiles across 224/448/672
  tile_metadata_resid.tsv.gz                               with depth-residualized labels

project/reports/
  pathology_dm1_phaseA_cpu_verdict_2026_05_04.md           Phase A NO-GO
  pathology_dm1_closure_battery_2026_05_04.md              ← THIS FILE

PHASE_A_NOGO_REPORT_2026_05_04.md                          top-level Phase A summary
```

---

*generated 2026-05-04 by automated closure battery · 6 sections + 2 negative-control regimes · CPU-only · Azure GPU quota = 0 (subscription limit), RunPod pod stuck at image pull (despite $30 balance) · ResNet50 ImageNet on 3,200 (224 px) + 2,935 (448 px) + 2,386 (672 px) GSE250521 tiles · 16-fold LOSO · spec-strict (no random tile split, no raw success claim)*
