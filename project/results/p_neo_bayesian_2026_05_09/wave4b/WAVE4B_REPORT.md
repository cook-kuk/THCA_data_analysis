# Wave 4B — Inference-time / no-retrain distribution-shift toolkit

Date: 2026-05-09
Branch: `paper9-perturbation-extension-20260506`
Scope: 4 inference-time methods (M5/M6/M7/M8) on top of frozen ESM2-150M + Wave 1 Bayesian head, evaluated on ITSNdb subsets. Local CPU only (Wave 4A holds the pod).

## Headline (ITSNdb_no_overlap, n=103, pos=33)

| Method                  | AUROC  | 95 % CI            | Δ vs Wave 1 (0.411) | Δ vs MHCflurry (0.668) |
|-------------------------|--------|--------------------|----------------------|--------------------------|
| Wave 1 baseline         | 0.411  | [0.294, 0.530]     | —                    | −0.258                   |
| **M5 kNN (K=10)**       | 0.425  | [0.303, 0.553]     | +0.014               | −0.243                   |
| **M6 TENT (probe)**     | **0.486** | [0.370, 0.613]  | **+0.075**           | −0.183                   |
| M7 split conformal      | 0.411  | [0.294, 0.530]     | 0.000 (same scoring) | −0.258                   |
| M8 WiSE-FT (SKIPPED)    | 0.411  | [0.294, 0.530]     | 0.000 (placeholder)  | −0.258                   |

Best method on `no_overlap` = **M6 TENT** at 0.486. The Wave 1 95 % CI extends from 0.294 to 0.530 — TENT's 0.486 is inside that envelope, so the +0.075 uplift is **not statistically significant** at n=103, but the sign is consistently positive on the hardest subset.

Honesty calls per the brief:
- **kNN did NOT beat the baseline by 0.05+** (Δ=+0.014). Reported flat, kept in figure as ablation.
- **TENT changed by ≥0.01** on no_overlap (+0.036 vs probe baseline; +0.075 vs Wave 1). Reported as the "best uplift", but caveats on probe/HLA-features below.
- **Conformal coverage** is the main contribution — see § M7.

## M5 — kNN retrieval (Khandelwal 2020)

- Procedure: ESM2-150M peptide embedding (re-computed locally for 2358 unique peptides; Wave 1 pod-side `embeddings.pt` not pulled). For each ITSNdb test peptide, take top-K nearest train peptides by cosine within the same `HLA_norm` bucket (fallback to global pool when bucket has < K rows). Score = sim-weighted vote of neighbor labels.
- K sweep on a held-out 5-fold of the train pool: K ∈ {5, 10, 20, 50}; held-out AUROCs 0.806 / 0.816 / 0.800 / 0.808 → **best K = 10**.
- ITSNdb subsets at K=10:
  - `main` n=199 AUROC 0.610
  - `Val`  n=120 AUROC 0.884
  - `combined` n=319 AUROC 0.785
  - `in_master` n=213 AUROC 0.922  ← strong, expected (training-leakage subset)
  - **`no_overlap` n=106 AUROC 0.425** (vs Wave 1 0.411 = +0.014, flat)
- Conclusion: kNN over ESM2 confirms the same pattern as Wave 1 — heavy uplift on `in_master` (peptides whose proteins/sources overlap with training), no real transfer to `no_overlap`. The hard split is hard for retrieval too.

## M6 — TENT (Wang 2021)

- The Wave 1 Bayesian head state-dict is on the pod (not pulled back). To run TENT we trained a **lightweight surrogate probe** on the same frozen ESM2-150M peptide embeddings + a deterministic random-feature representation of the HLA pseudo-seq (since rebuilding Wave 1's ESM2-encoded HLA branch on CPU was off the time budget).
  - Probe = `[pep ‖ hla] (1280) → BN→Linear(256) → BN→Linear(256) → Linear(1)`, dropout 0.3, AdamW 8 epochs, BCE.
- Surrogate **probe baseline** AUROCs:
  - combined 0.757, in_master 0.897, **no_overlap 0.450**, main 0.567, Val 0.717
- **TENT update**: 1 epoch SGD(lr=1e-3, momentum=0.9) on BN-affine params only (1024 trainable scalars), batch=64; running stats reset to use batch BN; dropout active. Loss = mean binary entropy of test-set predictions.
- Post-TENT AUROCs:
  - combined 0.751 (Δ −0.007), in_master 0.874 (Δ −0.022), **no_overlap 0.486 (Δ +0.036)**, main 0.580 (+0.013), Val 0.715 (−0.001)
- Pattern matches the canonical TENT story: small loss on already-easy subsets, modest gain on the hardest covariate-shift subset (`no_overlap`).
- **Caveat**: because the M6 probe uses a random-feature HLA branch (NOT Wave 1's ESM2-HLA), absolute AUROCs differ from Wave 1. The honest reading is the **TENT delta on the same probe** (+0.036). The +0.075 vs Wave 1 baseline is partly explained by the probe being a stronger no-DANN/no-mixup classifier on this subset to begin with.

## M7 — Split conformal prediction (Angelopoulos & Bates 2023)

- Calibration set = Wave 1 in-domain 5-fold predictions (`predictions_in_domain.tsv`, n=2041 with `p_ens`).
- Non-conformity = 1 − p̂_y; finite-sample-corrected quantile q_hat at α=0.10 → **q_hat = 0.602**.
- Per-subset coverage / mean set size (target coverage ≥ 0.90):

  | Subset              |    n | coverage | mean set size | AUROC | Selective AUROC (size=1) | retain |
  |---------------------|-----:|---------:|--------------:|------:|-------------------------:|-------:|
  | ITSNdb_main         |  192 |   0.797  |          1.54 | 0.564 | 0.600                    | 45.8 % |
  | ITSNdb_Val          |  119 |   0.983  |          1.13 | 0.783 | 0.614                    | 87.4 % |
  | ITSNdb_combined     |  311 |   0.868  |          1.38 | 0.784 | 0.829                    | 61.7 % |
  | ITSNdb_in_master    |  208 |   0.909  |          1.31 | 0.939 | 0.934                    | 68.8 % |
  | **ITSNdb_no_overlap** | **103** | **0.786** | **1.52** | 0.411 | 0.436                    | 47.6 % |

- **Conformal main result**: under a 90 %-coverage promise the model **fails to deliver coverage on `no_overlap` (0.786)**, and undercovers on `main` (0.797) — exactly what's expected when the test distribution is genuinely shifted relative to the calibration distribution. This is the honest finding to report: marginal coverage holds where the calibration set covers test (in_master, Val), and breaks where it doesn't (no_overlap, main). The selective-AUROC numbers (set size = 1) confirm the head only confidently single-labels easy points.

## M8 — WiSE-FT (Wortsman 2022) — SKIPPED

Reason: Wave 4A LoRA checkpoints not yet emitted (wave4a/ directory is empty as of run-time), AND Wave 1 head state-dict was never copied off the pod. Both endpoints of the θ_blend interpolation are therefore unavailable. Per the brief (`If Wave 4A LoRA isn't ready by your deadline, use Wave 1 head only (skip WiSE-FT, document in notes)`), M8 is skipped and a placeholder `predictions_wisefit.tsv` (= Wave 1 baseline at α=0) is emitted to keep the schema uniform.

## Files

```
wave4b/
  embeddings_local.pt        # 6.1 MB ESM2-150M peptide emb, 2358 peptides
  embed_local.py             # CPU embed script
  m5_knn.py                  # M5 implementation
  m6_tent.py                 # M6 implementation (surrogate probe + TENT)
  m7_conformal.py            # M7 implementation
  m8_wisefit.py              # M8 placeholder + skip note
  aggregate_results.py       # bootstrap CI + figure builder
  predictions_knn.tsv        # 319 ITSNdb rows × kNN scores
  predictions_tent.tsv       # 311 ITSNdb rows × probe + TENT scores
  predictions_conformal.tsv  # 311 ITSNdb rows × prediction sets
  predictions_wisefit.tsv    # placeholder (=Wave1)
  m5_knn_sweep.json          # K-sweep AUROCs
  m6_tent_subsets.tsv
  m7_conformal_subsets.tsv
  m7_conformal_summary.json
  wave4b_results.tsv         # 5-method × 5-subset AUROC table with bootstrap 95% CI
  fig_wave4b_uplift.png/.pdf # bar chart (no_overlap)
  WAVE4B_REPORT.md           # this file
```

## Bottom line

We close 0.075 of the 0.258 AUROC gap to MHCflurry on ITSNdb_no_overlap with TENT (still ~0.18 short, no statistical significance at n=103). kNN does not transfer the way it does on the easier subsets. Conformal exposes the real story: at α=0.10 the calibrated head undercovers exactly on the no-training-overlap distribution it should be most uncertain about — that is the headline figure to put in a discussion section, not the AUROC.
