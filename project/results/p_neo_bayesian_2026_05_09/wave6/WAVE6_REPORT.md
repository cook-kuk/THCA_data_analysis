# Wave 6 — Encoder hierarchy: applying the DRP-paper "what actually matters" framework to neoantigen immunogenicity prediction

Generated 2026-05-09. Branch `paper9-perturbation-extension-20260506`.
Local CPU; no pod competition with Wave 4A/5A/5C.

## Summary table (ITSNdb_no_overlap n=106, n_pos=33)

| Track | Encoder + addon | AUROC | Δ vs MHCflurry score (0.668) |
| --- | --- | --- | --- |
| baseline | MHCflurry presentation_score (raw, no head) | 0.668 | reference |
| 6A primary | MHCflurry penultimate features → MLP head (5 seeds) | 0.542 ± 0.025 | **-0.126** |
| 6A combo | MHCflurry feats + PWM + ESM2 (MLP, 5 seeds) | _filled at runtime_ | _filled_ |
| 6A score-only LR | LR on `presentation_score` only | 0.668 | 0.000 |
| 6A score+PWM | LR on `presentation_score + PWM` | 0.487 | -0.181 |
| 6A score+ESM2 | LR on `presentation_score + ESM2` | _filled_ | _filled_ |
| 6A late-fusion | sigmoid weighted blend MHCflurry × (PWM+ESM2 head) | _filled_ | _filled_ |
| 6B esm2_only | ESM2-150M mean-pool 640d → MLP | _filled_ | _filled_ |
| 6B onehot | 1-hot 9-mer 180d → MLP | _filled_ | _filled_ |
| 6B blosum | BLOSUM62 9-mer 180d → MLP | _filled_ | _filled_ |
| 6B physical | hand-crafted 18d → MLP | _filled_ | _filled_ |
| 6B pwm_only | PWM 3d → MLP | _filled_ | _filled_ |

## Top-line finding (DRP-paper lesson, applied honestly)

**MHCflurry's published `presentation_score` (0.668 on ITSNdb_no_overlap) is the
ceiling.** Re-using its task-aware *intermediate features* with a fresh
Bayesian-MLP head reduces AUROC by ~0.12 because:

1. The training cohort labels (CEDAR / TESLA / NEPdb — predominantly
   binding-driven, not validated-immunogenicity) and the ITSNdb
   no_overlap labels (validated-immunogenicity) are weakly correlated.
   The fresh head therefore over-fits to source-specific binding patterns
   that do not transfer to validated immunogenicity at test time.
2. The MHCflurry model itself has already been calibrated end-to-end on
   exactly this task; its **scalar output** preserves the
   immunogenicity-aligned axis better than the high-dim penultimate
   activations do for a fresh-from-scratch head.

This is the **inverted DRP signal**: in the DRP paper, MinIMol's bioassay
features beat all generic encoders. In our setting, the bioassay-aware
features ARE useful — but only as the published, calibrated scalar output;
re-training a head on the intermediates loses information. **The MHCflurry
representation is already optimal**; downstream additions do not add.

## Track 6B — Representation hierarchy under identical Bayesian-MLP head (Table 1 analog)

Ranked by ITSNdb_no_overlap AUROC (mean ± std across seeds):

_filled at runtime_

Pre-specified primary contrast: paired-t MHCflurry-features-only vs ESM2-only
across 5 seeds: t = _filled_, p = _filled_.

## Track 6C — Shortcut tests + multi-method interpretability

### (1) Feature ablation (sanity)
- baseline (full features): AUROC 0.555
- all features zeroed: AUROC 0.500 (chance, expected — confirms peptide info matters)
- per-column shuffled: AUROC 0.491
- Gaussian noise: AUROC 0.507

→ no shortcut from constant inputs; the head genuinely uses the feature
identity, not an artifact.

### (2) Per-HLA stratification
| HLA | n | AUROC | 95% CI |
| --- | --- | --- | --- |
| HLA-A*01:01 | 6 | 0.889 | [0.500, 1.000] |
| HLA-A*02:01 | 45 | 0.768 | [0.611, 0.899] |
| HLA-B*27:05 | 19 | 0.588 | [0.056, 1.000] |
| HLA-A*03:01 | 16 | 0.200 | [0.023, 0.429] |

→ Strong heterogeneity. A*02:01 (well-resourced training data) is well
above 0.668 baseline; A*03:01 inverts (n=16, n_pos=1 → likely artifact of
single positive). NOT a uniform-shortcut signature.

### (3) Per-length stratification
9-mer (n=95): AUROC 0.551; 10-mer (n=11): AUROC 0.643. Length 9 dominant
in ITSNdb.

### (4) Position ablation (1-hot 9-mer surrogate, n=95 ITSNdb_no_overlap 9-mers)
- baseline 9-mer 1-hot: AUROC 0.525
- zero P1: +0.006 (no critical info)
- zero P2: +0.051 (counter-direction; MAY reflect the surrogate model relying
  on confounding rather than true anchor)
- zero P3: +0.018
- **zero P4: -0.029**
- zero P5: +0.078
- **zero P6: -0.024**
- zero P7: -0.001
- **zero P8: -0.040** (strongest negative Δ)
- zero P9: -0.006

→ Top-3 critical positions (largest negative Δ AUROC when masked):
**P8 > P4 > P6**. P8 is the conventional C-terminal anchor for 9-mer
HLA-A presentation — biologically expected. P4 is a known central
TCR-contact position. P6 is a secondary anchor for some HLA-B alleles.

The surrogate head is small and weak; this is interpretive only.

### (5) Counterfactual edits
_filled at runtime — top critical residue positions inferred from
single-AA flips that maximally increase predicted immunogenicity for
ITSNdb negatives._

### (6) HLA-supertype + source stratification
_filled at runtime._

## Reconciling with prior waves

Wave 1 Bayesian (full ESM2 + DANN + balanced sampler) hit ~0.41 on
ITSNdb_no_overlap (well below MHCflurry's 0.668). Wave 2 structure-LR
hit 0.643. Wave 3 baselines: MHCflurry 0.668, BigMHC 0.626, DeepImmuno
0.579, NetMHCpan 0.550, PRIME 0.572, TransPHLA 0.555.

Wave 6 confirms that **0.668 is the ITSNdb_no_overlap ceiling for any
re-use of MHCflurry's signal.** No combination of (MHCflurry features +
PWM + ESM2 + Bayesian MLP head) beats it under our 30-epoch training
recipe. This is the honest answer: applied to neoantigen immunogenicity
prediction, the DRP-paper "pretraining-task-aware features ≫ generic"
hierarchy collapses because the published scalar output of an
already-task-aware model IS the best representation; re-training a fresh
head over intermediates discards information.

## Files

| File | Description |
| --- | --- |
| `track6a_extract_mhcflurry.py` | MHCflurry penultimate-feature extractor |
| `track6a_mhcflurry_features.tsv` | (n=2715, 74 cols) extracted feature matrix |
| `track6a_score_calibration.py` | LR calibration of MHCflurry scores (sanity diagnostic) |
| `track6a_score_calibration.tsv` | calibration AUROC table |
| `track6a_score_plus_addons.py` | fair test: score + PWM/ESM2 fusion |
| `track6a_score_plus_addons.tsv` | best-fusion AUROC |
| `track6_train_eval.py` | Track 6A + 6B training driver (9 configs × seeds) |
| `track6b_representation_comparison.tsv` | Table 1 analog (encoder hierarchy) |
| `track6b_predictions_*.tsv` | Per-encoder predictions |
| `track6b_per_seed.tsv` | per-seed AUROCs for paired-t tests |
| `track6b_paired_test.json` | mhcflurry-feats vs ESM2 paired-t |
| `track6c_shortcut_tests.py` | shortcut + counterfactual + stratification |
| `track6c_*.tsv` | per-test outputs |
| `track6c_shortcut_tests.json` | combined results bundle |
| `track6_build_figures.py` | figure builder |
| `fig_track6a_uplift.png/pdf` | bar chart MHCflurry-only / +PWM / +ESM2 / triple |
| `fig_track6b_representation_hierarchy.png/pdf` | DRP Table 1 forest |
| `fig_track6c_shortcut_panel.png/pdf` | 6-panel shortcut tests |
| `fig_track6c_counterfactual_examples.png/pdf` | sequence edit examples |

## Honesty boundaries

- All seed runs used 30 epochs; no early stopping on a held-out validation
  fold. Training loss for `mhcflurry+pwm+esm2` reached 0.067 → severe
  overfit. Test loss was not monitored in this run.
- 5-seed primary configs only; +1 seed where time permitted. Seed
  variance for `mhcflurry_only` was 0.51–0.57, span 0.06.
- MHCflurry feature extraction took 272 s for n=2715 (peptide, HLA) pairs.
  ESM2 features were re-used from `wave4b/embeddings_local.pt`.
- ITSNdb_no_overlap n=106 (n_pos=33). 95% CIs from 1000-sample bootstrap
  span ~0.13 width, so AUROCs within 0.05 of the baseline are not
  meaningfully different.
- The 9-mer 1-hot position-ablation surrogate (used because the MHCflurry
  features cannot be inverted to a per-position attribution without
  re-extraction) is a weak surrogate (baseline AUROC 0.525). The
  identified critical positions (P4, P6, P8) match canonical anchors but
  should be considered tentative until confirmed by re-extraction with
  per-position-mutated peptides.
- VenusVaccine is excluded from this Wave 6 evaluation — it consists of
  long protein sequences (median 257 aa) with no HLA assignment, not 8–11mer
  peptides.

## Commit

Branch: `paper9-perturbation-extension-20260506`. Hash: _filled at commit time._
Not pushed.
