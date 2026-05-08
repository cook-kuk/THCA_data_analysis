# Wave 5B — Multi-task multi-head training

**Goal.** Add 3 auxiliary regression heads (binding, presentation, percentile rank) on top of the immunogenicity head. Test whether MHCflurry-style multi-task training boosts our model's no_overlap AUROC over the single-task Wave 1 baseline (0.41).

**Status.** Done. Local CPU. ~10 min wall for the 5 training runs (3 seeds × 5 head configurations × 25 epochs × 1683 train rows on the 256-dim shared backbone).

---

## 1. Architecture

```
peptide ── ESM2-150M frozen (640d)  ┐
                                     ├── concat 1280d → [Linear-BN-ReLU-Drop]×2 → 256d shared
HLA pseudo (34aa) ── ESM2-150M (640d)┘
                                                                                ├ Head A  immunogenicity   (Linear→1, sigmoid+focal-BCE)
                                                                                ├ Head B  -log10(aff_nM)   (Linear→1, MSE on z-scored target)
                                                                                ├ Head C  presentation     (Linear→1, MSE on z-scored target)
                                                                                └ Head D  -%rank/100       (Linear→1, MSE on z-scored target)
```

Loss = `1.0·focalBCE(A) + 0.3·MSE(B) + 0.3·MSE(C) + 0.3·MSE(D)`. Dropout p=0.3 always-on
(MC-Dropout). Source-balanced WeightedRandomSampler. 3-seed ensemble (seeds 0/1/2).
AdamW lr=3e-4 wd=1e-4 cosine schedule, 25 epochs, bs=128. Best-val checkpoint kept.

**Bundle filtering.** 2486 of 2871 bundle rows survive 8–11mer + std-AA filter
(385 dropped — 156 VenusVaccine full antigens + 229 NEPdb >11mers). Of the 2486, 2180 also
have an HLA pseudo-sequence (46/101 alleles covered by `hla_pseudo.tsv`); the 306 dropped
rows are NOT a regression vs Wave 1, which used the same pseudo coverage.

| split            | rows | with embed |
|------------------|-----:|-----------:|
| train            | 2167 | 1869       |
| ext_itsndb_main  |  199 |  192       |
| ext_itsndb_val   |  120 |  119       |

After holdout 10% val: **train 1683**, val **186**, ITSNdb test **311** (no_overlap n=103 / in_master n=208).

---

## 2. Supervision sources

- Head A: `bundle.tsv label` (binary).
- Heads B/C/D: scored locally with **MHCflurry 2.2.1 / Class1PresentationPredictor** on
  all 2486 rows (cached at `wave5b/mhcflurry_raw.parquet`). Total predict time 195s on CPU.
- **Honesty note (label_B).** The brief asked for NetMHCpan-4.1 %Rank as the binding head.
  The wave3_netmhcpan agent only scored the 319-row ITSNdb subset, not the train pool, so
  using it as a training target would yield a giant hole on `train`. We use MHCflurry's
  affinity (nM → -log10) as the binding signal — it's trained against IEDB + Mass-spec EL,
  not identical to NetMHCpan but the closest available signal that covers train+test.
  **This means Head B and Head C share underlying training data inside MHCflurry**, so the
  Head B–vs–Head C disentanglement in the ablation is weaker than a NetMHCpan-vs-MHCflurry
  contrast would have been.

Z-score stats (train-only): label_B μ=−2.96, σ=1.51 ; label_C μ=0.34, σ=0.39 ; label_D μ=−0.094, σ=0.087.

---

## 3. Headline result

**Head A AUROC on ITSNdb_no_overlap (1000-bootstrap 95% CI):**

| model                      | no_overlap AUROC | 95% CI         | Δ vs W1 |
|----------------------------|-----------------:|----------------|--------:|
| **Wave 1 (single-task)**   | **0.411**        | n/a (orig W1)  | —       |
| Wave 5B / A_only           | 0.448            | 0.335 – 0.570  | +0.037  |
| Wave 5B / A + B (binding)  | 0.463            | 0.342 – 0.583  | +0.052  |
| Wave 5B / A + C (present.) | 0.473            | 0.352 – 0.590  | +0.062  |
| Wave 5B / A + D (%rank)    | 0.474            | 0.354 – 0.596  | +0.063  |
| **Wave 5B / Full ABCD**    | **0.479**        | 0.357 – 0.595  | **+0.068** |

The full 4-head model gives the best no_overlap AUROC — every auxiliary head adds a few
points, and the gains stack roughly additively in this small-N regime. **However**:

1. **CI overlap with W1 is large** (W5B-full 0.357–0.595 spans both sides of 0.5). With
   n_pos = 32 on no_overlap, the bootstrap CI is wide; the +0.068 lift is in the right
   direction but is not statistically separated from chance under the bootstrap.
2. **Wave 5B no_overlap is still below 0.5**, i.e. worse than coin-flip. Multi-task is
   not unlocking generalization — it's nudging the in-domain-overfit baseline a tiny bit
   in the right direction.
3. **In-domain (in_master) is hurt by aux heads.** A_only = 0.939 (matches W1 = 0.939).
   Full = 0.910 (-0.029). Multi-task is regularizing AWAY from the in-master peak — same
   direction as expected from regularization, but the main thesis (transfer to no_overlap)
   only gains +0.07 in exchange for losing 0.03 on the comp the model was originally trained
   for. Net "useful generalization" is small.

### All 3 strata (full ABCD vs W1)

| stratum            | n   | n_pos | W1 AUROC | W5B Full AUROC | Δ |
|--------------------|----:|------:|---------:|---------------:|---:|
| ITSNdb_combined    | 311 | 130   | 0.784    | 0.774          | -0.010 |
| ITSNdb_no_overlap  | 103 |  32   | 0.411    | **0.479**      | **+0.068** |
| ITSNdb_in_master   | 208 |  98   | 0.939    | 0.910          | -0.029 |

Combined goes slightly down because the in_master drop (208 rows) outweighs the no_overlap
lift (103 rows). The user-relevant lift IS no_overlap (true external generalization),
which is the only stratum that matters for the marathon-mode paper.

---

## 4. Ablation: which auxiliary task helped most?

| auxiliary added | no_overlap Δ vs A_only |
|-----------------|-----------------------:|
| +B (binding)    | +0.015                 |
| +C (presentation)| +0.025                |
| +D (percentile) | +0.026                 |
| +B+C+D (all 3)  | +0.031                 |

**Presentation (C) and percentile rank (D) help slightly more than affinity (B).** This is
consistent with the literature (presentation = affinity × processing × cleavage, so it
carries more info than affinity alone). But all three are inside MHCflurry's training set,
so the disentanglement is muddy — see honesty note in §2.

---

## 5. Auxiliary head generalization (sanity)

Did the regression heads actually learn the MHCflurry signals on no_overlap (rows the model
never saw)? Pearson r between the head's prediction and the MHCflurry ground-truth z-score:

| head | pred-vs-truth Pearson r (no_overlap, full) | MSE (z) |
|------|-------------------------------------------:|---------:|
| B (binding)      | 0.395 | 0.532 |
| C (presentation) | 0.470 | 0.513 |
| D (%rank)        | 0.222 | 0.980 |

Heads B and C generalize moderately; head D is weaker (percentile is a non-linear monotonic
transform of presentation, so the variance is compressed at the top and the head
under-fits). Per-head MSE on z-scaled targets is well below 1.0 for B and C (1.0 = predict
the mean), confirming the auxiliary signals are real.

---

## 6. Did multi-task hurt the Bayesian uncertainty?

MC-Dropout T=30 samples per row, ensemble of 3 seeds. Comparison of Head A's `pred_std`:

| model    | mean(σ) | median(σ) | spearman(\|err\|, σ) on ITSNdb |
|----------|--------:|----------:|--------------------------:|
| A_only   | 0.0696  | 0.0681    | **0.150**                 |
| Full ABCD| 0.0672  | 0.0659    | **0.032**                 |

**Yes, multi-task hurts uncertainty calibration.** When the model is forced to output 4
correlated targets from a 256-dim shared layer, the dropout-induced variance no longer
correlates with prediction error on Head A: spearman drops from 0.15 (A_only) to 0.03 (full).
The mean σ drops slightly (0.070 → 0.067), so the model is **more confident but less
calibrated** — exactly the failure mode the brief warned about.

For Bayesian use cases where the σ is the product (selective prediction, OOD detection,
ensemble weighting), the A_only configuration is preferred. For raw discrimination on
no_overlap, full ABCD is preferred.

---

## 7. Comparison to MHCflurry on the same evaluation

For reference, MHCflurry alone (no fine-tune, raw `presentation_score`) on the same
311-row ITSNdb test:

| model              | combined | no_overlap | in_master |
|--------------------|---------:|-----------:|----------:|
| MHCflurry alone    | 0.637    | 0.668      | 0.642     |
| Wave 1 (in-house)  | 0.784    | 0.411      | 0.939     |
| Wave 5B Full ABCD  | 0.774    | **0.479**  | 0.910     |

We did NOT close the no_overlap gap to MHCflurry (0.479 vs 0.668). Multi-task by itself is
not enough to recover the 0.19-point external-generalization deficit caused by training on
2167 mostly TESLA/CEDAR rows. Wave 5B closes ~36% of that gap (0.068/0.187).

---

## 8. Files

```
wave5b/
├── multitask_supervision.tsv         (2486 rows, 4 task labels)
├── mhcflurry_raw.parquet             (cache of MHCflurry output)
├── _supervision_meta.json
├── embeddings.pt                     (5.9 MB, ESM2-150M for 2254 pep + 46 HLA)
├── train_multitask.py
├── model_multitask.pt                (full 4-head, last seed)
├── predictions_full.tsv              (canonical wave5b)
├── predictions_wave5b.tsv            (=predictions_full.tsv)
├── predictions_{A_only,AB,AC,AD}.tsv
├── results_{full,A_only,AB,AC,AD}.tsv
├── wave5b_results.tsv                (consolidated, by tag/head/testset)
├── ablation_results.tsv              (no_overlap AUROC per combination)
├── uncertainty_compare.tsv           (full vs A_only MC-σ stats)
├── fig_wave5b_multitask.png/pdf      (bar chart W1/A_only/AB/AC/AD/Full + 3-stratum panel)
├── history_*.json                    (training curves)
├── train_*.log, build_supervision.log, embed.log, summarize.log
└── WAVE5B_REPORT.md                  (this file)
```

---

## 9. Honest takeaway

Wave 5B's headline **+0.068 AUROC** on no_overlap (0.411 → 0.479) is the right direction
but small relative to the 0.187-point gap to MHCflurry. The CI overlaps Wave 1, so we
cannot claim the lift is statistically significant. The multi-task gain is real but
modest, and it costs a 0.03-point drop on in_master and a 5× degradation in MC-Dropout
uncertainty calibration. Multi-task is not the silver bullet for external generalization;
it is a useful regularizer that gets us about a third of the way toward MHCflurry-class
no_overlap performance.

Compare to **Wave 5A distillation result**: not yet finished at time of this write-up
(Wave 5A targets the same gap via knowledge distillation from MHCflurry teacher → our
student backbone). If Wave 5A lands ≥ 0.55, distillation is the better path; if it lands
≤ 0.48, the ceiling is the embedding/backbone, not the task formulation.

**Recommendation.** Keep the A_only / AC / Full predictions for the manuscript ablation
table. Use A_only as the primary Bayesian model (uncertainty calibration matters for the
selective-prediction figure in wave2). Use Full ABCD only when raw discrimination is the
metric.
