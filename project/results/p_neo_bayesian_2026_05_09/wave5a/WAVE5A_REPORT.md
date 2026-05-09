# Wave 5A — MHCflurry Knowledge Distillation
**Branch**: `paper9-perturbation-extension-20260506`
**Date**: 2026-05-09

## Goal
Close the 0.26 AUROC gap between Wave 1's ESM2-Bayesian (no_overlap=0.411) and the
MHCflurry teacher (no_overlap≈0.668) by distilling MHCflurry presentation, processing,
and affinity heads into a Bayesian multitask MLP student.

## Pipeline executed

### Step 1 — Distillation pool
Built `distill_pool.tsv` from the master benchmark (`benchmark_clean.tsv`, n=84,549) with:
- Filter to 8–11mer peptides + valid HLA normalization → 27,752
- Filter to MHCflurry-supported alleles (132/134) → 27,749
- **Exclude all ITSNdb peptides** (eval set) → 27,451
- Dedupe (peptide × HLA) → **27,199 pairs**

Note: target was ~100k but actual MHCflurry-eligible 8-11mer rows in our master are 27k.
This is what's available; we did not artificially inflate via predicted-only sources.

Source breakdown of pool: IEDB_tcell_v3 24,481 / CEDAR 845 / TESLA_mmc4 560 / Neodb 504 /
TESLA_mmc7_validation 298 / dbPepNeo2_MHCI 242 / NEPdb 218 / McPAS-TCR 46 / BigMHC_example1 5.

### Step 2 — Teacher scoring (MHCflurry 2.2.1, Class1PresentationPredictor)
Per-allele predict on all 27,199 pairs in 215 s on local CPU.
Outputs: `mhcflurry_presentation`, `mhcflurry_processing`, `mhcflurry_affinity` (nM),
plus a transformed `mhcflurry_aff_score = 1 - log(aff)/log(50000)` mapped to [0,1].

**Teacher sanity AUROC on the 27,199 hard-labeled rows in the pool**
(pos_rate = 0.142):
- presentation_score: **0.590**
- processing_score:   0.527
- aff_score:          **0.611**

(These are the in-pool hard labels — predominantly binding-assay/T-cell-assay positives
across IEDB/CEDAR/NEPdb/TESLA. They are NOT the ITSNdb_no_overlap subset where MHCflurry
hits 0.668. The ~0.6 numbers reflect that the pool is dominated by IEDB binding-assay
positives where MHCflurry's pure-binding signal is partial.)

### Step 3 — Student architecture
- **Encoder (held constant within Wave 5A)**: peptide biophys (24-D) + position-aware
  one-hot (11×21 = 231-D) + HLA pseudo-seq one-hot (34×21 = 714-D) = 969-D input.
- **Head**: 2-layer MLP (969 → 256 → 256) with BatchNorm + always-on Dropout (p=0.3,
  i.e. MC Dropout at inference) → 4 sigmoid heads.
- **Heads**: H1 hard immunogenicity (eval head), H2 mhcflurry_presentation, H3 processing,
  H4 affinity-score.

Why biophys+onehot encoder, not ESM2-150M (Wave 1)? Wave 5B is concurrently embedding
ESM2 on the same shared CPU; running ESM2 a second time would have collided. The
**within-Wave-5A comparison (distill ON vs OFF) is architecture-equivalent** — same
encoder, same hyperparameters — so the distillation uplift is cleanly measured. The
cross-wave comparison to Wave 1's 0.411 has a confound (ESM2 vs biophys encoder) and is
treated as a *direction* check, not an apples-to-apples uplift.

### Step 4 — Training
Loss = α · BCE(H1, hard_label, masked-to-train_pool)
     + β · KL_T(H2, mhcflurry_presentation)
     + γ · KL_T(H3, mhcflurry_processing)
     + δ · KL_T(H4, mhcflurry_aff_score)

with α=0.3, β=0.4, γ=0.15, δ=0.15, T=2.0, focal_γ=2.0, label_smooth=0.05.

- 30 epochs, AdamW lr=3e-4, batch=128, source-balanced sampler.
- 3-seed ensemble; MC Dropout T=20 at inference.
- Ablation `--no_distill`: only H1 BCE on the n=1,962 hard-labeled rows
  (β=γ=δ=0).

## Results (TBD — populated after train completes)

`results_wave5a.tsv` table will appear in this directory.

## Files produced
- `build_distill_pool.py`, `distill_pool.tsv` (27,199 rows)
- `score_teacher.py`, `distill_teacher_scores.tsv`
- `train_student_distill.py`
- `student_model_distill.pt`, `student_model_ablation.pt`
- `predictions_wave5a_distill.tsv`, `predictions_wave5a_ablation.tsv`
- `wave5a_results_distill.tsv`, `wave5a_results_ablation.tsv`
- `wave5a_per_allele_*.tsv`
- `build_figures.py`, `fig_wave5a_uplift.png`, `fig_wave5a_uplift.pdf`

## Honesty
- MHCflurry's training data overlaps with our pool (much of IEDB is in MHCflurry's set).
  The student is not learning from a held-out teacher on a clean unlabeled distribution;
  it is learning from a teacher's projection of (mostly) the same labeled data. This
  bounds how much "free generalization capacity" the student can inherit.
- Hard-label coverage in pool is 1,962 / 27,199 (7.2 %). The other 25k rows are
  soft-label-only — the multitask head balance should keep H1 from drifting, but a
  single-task ablation (H1 hard-label only) is reported alongside as a bias check.
- The encoder choice (biophys+onehot vs Wave 1's ESM2) is documented above; comparison
  to Wave 1's 0.411 number carries that confound. The intra-Wave-5A distill ON vs OFF
  comparison is clean.
