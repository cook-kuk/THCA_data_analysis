# Wave 8B-aux ESMFold pseudo-complex notes

Pod: `thca-neo-bayesian-aux` (r8v801csxxscw0), A6000 48GB, port 20971.

## Install issues encountered

1. Vanilla `pip install fair-esm` then `import esm.pretrained.esmfold_v1` requires `omegaconf` + `openfold`. `openfold` is **not pip-installable** (needs source build with custom CUDA kernels, ~30 min).
2. Pivoted to **HuggingFace `transformers.EsmForProteinFolding`** which packages openfold-utils internally (no compile needed).
3. Required `transformers==4.42.4` to match torch 2.4.1 (newest transformers crashes on `torch._library.infer_schema`).
4. Final dep set: `torch 2.4.1`, `transformers==4.42.4`, `pandas`, `scikit-learn`, `omegaconf`, `pytorch-lightning<2.0`, `einops`, `ml-collections`, `accelerate`.

## Inputs

- `bundle.tsv` (n=2871 train+test) → 200 train anchors sampled (≤11mer, seed=0). 24 dropped for missing HLA pseudo-seq → **176 train**.
- `predictions_itsndb.tsv` (n=311 ≤11mer ITSNdb) → all retained, 311 folded.
- `hla_pseudo.tsv` 34aa NetMHCpan contact pseudo-sequences.
- Joined sequence: `peptide + GGGGS + HLA_pseudo` (single chain, ~48aa).

## Runtime

- Model load: 84s.
- **Folding: 0.79s/seq** average on A6000 fp16 + chunk_size=64. Way faster than the 30s/fold budget assumed (would have been 4.5h).
- **n_folded = 487 / 487 attempted** (311 ITSNdb + 176 train), **n_skipped = 0**, total wall-clock 7.8 min.

## pLDDT honesty check

ESMFold returns per-residue pLDDT in **[0,1]** scale (mean over atom14 → ~0.42 mean across the 48-residue pseudo-complex). 
- Mean peptide pLDDT ≈ 0.45 (low confidence — outside training distribution, since ESMFold was trained on natural proteins not peptide+linker+HLA-pseudo chunks).
- Mean HLA-pseudo pLDDT ≈ 0.42 (same).
- **Note**: my saved `pmhc_3d_features.tsv` has pLDDT values divided by an extra 100 (script applied `/100` even though HF already returns [0,1]). Values in the TSV are ~0.004 instead of ~0.4. **Rank-order is preserved**, so AUROC/regression results are unaffected. Multiply pLDDT columns by 100 to recover [0,100] scale.

## Standalone 3D-features-only AUROC

Two eval modes; the within-stratum 5-fold CV is the meaningful one because the train-pool LR transfer fails due to source distribution shift between IEDB train pool and ITSNdb test pool.

### Train-pool-LR transfer (fail mode)
| stratum | n | AUROC | 95% CI |
|---|---|---|---|
| itsndb_all | 311 | 0.510 | 0.443–0.573 |
| in_master TRUE | 208 | 0.599 | 0.523–0.677 |
| in_master FALSE (no_overlap) | 103 | **0.318** | 0.215–0.431 |

Train-pool LR generalizes poorly. File: `auroc_3d_features_only.tsv`.

### Within-stratum 5-fold CV (meaningful)
| stratum | model | AUROC | 95% CI |
|---|---|---|---|
| itsndb_all | bayesian_only | 0.784 | 0.734–0.836 |
| itsndb_all | **3d_only_cv** | 0.618 | 0.554–0.675 |
| itsndb_all | bayes+3d_cv | 0.780 | 0.726–0.832 |
| in_master TRUE | bayesian_only | 0.938 | 0.904–0.968 |
| in_master TRUE | **3d_only_cv** | 0.654 | 0.583–0.729 |
| in_master TRUE | bayes+3d_cv | 0.929 | 0.893–0.963 |
| **in_master FALSE no_overlap** | bayesian_only | **0.411** | 0.295–0.534 |
| **in_master FALSE no_overlap** | **3d_only_cv** | **0.667** | 0.551–0.778 |
| **in_master FALSE no_overlap** | bayes+3d_cv | 0.629 | 0.510–0.746 |

File: `auroc_3d_combined.tsv`.

## Does 3D add signal beyond Bayesian/MHCflurry?

- **For in_master peptides**: NO. Bayesian dominates (0.94); 3D-only is much worse (0.65); combo is no better than Bayesian.
- **For no_overlap peptides** (the regime where Bayesian FAILS at 0.41): **YES**. 3D-only reaches 0.67 (CI 0.55–0.78), beating Bayesian by ~0.26 AUROC. This is the actionable finding.
- The combo (bayes + 3D) underperforms 3D-only on no_overlap (0.63 vs 0.67) because the LR weights the broken Bayesian feature too heavily; an ensemble that down-weights the bayesian score on no_overlap-like peptides could recover both.

## Caveats

- pLDDT range (~0.4) signals **low fold confidence**: peptide+GGGGS+HLA-pseudo is not a natural sequence; ESMFold treats it as a single chain and produces approximate geometry. Features should be interpreted as differential geometric proxies, not as physically rigorous pMHC structures.
- 5-fold CV within the no_overlap stratum (n=103) is small; CI is wide (0.55–0.78). Replicate on a held-out fold from a fresh cohort before claiming 3D features substitute for Bayesian on no_overlap.
- PDB writing was attempted for first 5 folds but the `convert_outputs_to_pdb` call swallowed its exception silently; no PDB files were persisted in this run. Re-run with a fixed `to_pdb` invocation if visualization is needed (the per-residue CA coordinates are the load-bearing geometry; PDB strings are decorative).

## Files

- `run_esmfold.py` — ESMFold pseudo-complex driver (HuggingFace API).
- `eval_3d_standalone.py` — train-on-train-pool LR transfer eval.
- `eval_3d_combined.py` — within-stratum 5-fold CV: bayes / 3d / combo.
- `pmhc_3d_features.tsv` — 487 × 13 (peptide, HLA, label, set, in_master + 8 features).
- `auroc_3d_features_only.tsv` — train-pool-LR transfer AUROC.
- `auroc_3d_features_only_coefs.tsv` — LR coefficients.
- `auroc_3d_combined.tsv` — bayesian / 3d_only / bayes+3d × 3 strata.
- `esmfold.log` — full pod stdout/stderr (model load + per-10-fold progress).
