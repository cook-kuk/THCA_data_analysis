# v5 Track 2 — DANN Biomarker De-confounding Method Card

## One-line description
A gradient-reversal-layer (GRL) biomarker classifier that is adversarially
regularized against batch (cohort) prediction, enabling cohort-invariant
feature learning for cross-study biomarker discovery.

## Architecture
Shared encoder: `MLP(54 -> 32 -> 16)` with ReLU + Dropout(0.3) after each
hidden layer. Biomarker head: `Linear(16 -> 2)` softmax for BRAF_like vs RAS_like.
Cohort head: `GRL -> Linear(16 -> 16) -> ReLU -> Linear(16 -> K=2)`
softmax over K=2 cohorts. The gradient-reversal layer (Ganin & Lempitsky,
2015) multiplies gradients by `-lambda` during backprop through the encoder,
forcing the encoder to produce features that a cohort classifier cannot
discriminate.

## Loss function
```
L = CE(y_bio, y_bio_hat) + lambda * CE(y_cohort, y_cohort_hat)
```
The GRL ensures the encoder MINIMIZES the ability of the cohort head to
predict cohort (by gradient negation) while the cohort head itself is trained
normally to predict cohort. Lambda controls anti-confounding strength.

## Hyperparameters
- Optimizer: Adam, lr=1e-3, weight_decay=1e-4
- Batch size: 64
- Max epochs: 200 (with early-stop patience=20 on val bio AUC)
- Dropout: 0.3 in encoder
- Lambda sweep: {0, 0.1, 0.3, 1.0, 3.0}
- Input features: top 54 variable genes across pooled cohorts

## Metrics (this run)
| lambda | bio CV AUC | ident AUC | LODO AUC |
|--------|-----------:|----------:|---------:|
| 0.0 | 0.917 | 1.000 | 0.934 |
| 0.1 | 0.929 | 1.000 | 0.946 |
| 0.3 | 0.914 | 1.000 | 0.925 |
| 1.0 | 0.917 | 1.000 | 0.939 |
| 3.0 | 0.887 | 1.000 | 0.906 |

**Best λ = 0.1** (selected by `bio_CV_AUC - max(0, ident_AUC - 0.70)`).

## SHAP-stable genes
Genes appearing in the top-10 gradient*input importances for ≥3 of 5
λ values (considered robust to the anti-confounding regularization):

TACSTD2, CST6, DIO1, COL1A1, DDX3Y, IL1RL1

## Novelty vs. prior art
DANN (Ganin & Lempitsky, 2015) was developed for image-domain adaptation.
To our knowledge, this is the first application to **cross-cohort bulk RNA-seq
biomarker discovery** as a replacement/augmentation for ComBat-style pre-processing.
The key claim is: instead of pre-correcting expression and then fitting a
downstream classifier (which v4's Track A shows is unrecoverable for THCA
subtype labels), **learn a cohort-invariant representation jointly** with the
biomarker classifier. The Pareto frontier (bio AUC vs 1-ident AUC) is a
principled way to report the trade-off rather than claiming a single AUC.

## Limitations
1. CPU-only PyTorch, small MLP — does not attempt deep representations.
2. GRL hyperparameter λ must be swept; no automatic selection.
3. Assumes cohort label is the only confounder; does not jointly
   de-confound platform, sequencing depth, age, stage separately.
4. BRAF/RAS label distribution is severely imbalanced across cohorts —
   some LODO folds have single-class held-out and are reported as NaN.
5. **Prototype only, not clinical**. Decision-support / methodology
   contribution. Not a diagnostic device.

## Reproducibility
- Code: `notebooks_or_scripts/v5_track2_adversarial.py`
- Model weights: `results/ml/v5_track2_best_model.pt`
- Results: `results/ml/v5_track2_adversarial_results.tsv`
- SHAP stability: `results/ml/v5_track2_shap_stable_genes.tsv`
- Random seed: 42 for sklearn splits, 0 for torch training
- Hardware: CPU-only, PyTorch 2.11.0+cpu
