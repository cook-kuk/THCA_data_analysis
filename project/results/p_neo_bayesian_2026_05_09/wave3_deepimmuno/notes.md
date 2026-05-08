# DeepImmuno-CNN — wave3 ITSNdb leakage-stratified scoring

**Algorithm:** DeepImmuno-CNN (Li et al., *Briefings in Bioinformatics* 2021,
doi:10.1093/bib/bbab160). Lightweight CNN over peptide AAindex-PCA + HLA
pseudo-sequence (paratope-residue) embeddings, sigmoid output 0-1.

**Repo:** https://github.com/frankligy/DeepImmuno (commit cloned 2026-05-09 at
`/tmp/DeepImmuno`).

## Install

`git clone https://github.com/frankligy/DeepImmuno /tmp/DeepImmuno` — cleanly
provides `deepimmuno-cnn.py` plus `data/after_pca.txt`,
`data/hla2paratopeTable_aligned.txt` (62 HLAs) and the trained TF checkpoint at
`models/cnn_model_331_3_7/`. No `requirements.txt` install needed — the
project venv already has `tensorflow==2.21.0`, `numpy`, `pandas`, sklearn.

## Issue 1 — Keras 3 cannot read the TF2 object-graph checkpoint directly

The shipped `deepimmuno-cnn.py` calls
`cnn_model.load_weights('./models/cnn_model_331_3_7/')`. Under
TF 2.21 / Keras 3 this fails:

> `ValueError: File format not supported: filepath=./models/cnn_model_331_3_7/.
> Keras 3 only supports V3 .keras files, .weights.h5 files, legacy H5 format
> files, or Orbax checkpoints.`

The shipped checkpoint is the legacy TF2 trackable-object-graph format with
hidden-file names (`.data-00000-of-00001`, `.index`) and a `checkpoint`
manifest pointing to prefix `.`.

**Workaround:** read the checkpoint with `tf.train.load_checkpoint`, enumerate
`layer_with_weights-N/{kernel,bias,gamma,beta,moving_mean,moving_variance}`
tensors, and assign them into a freshly built `seperateCNN()` Keras model by
shape-matching layer-by-layer (Conv2D / BatchNormalization / Dense). Verified
the assignment is unique by shape (the two pep-side and two HLA-side Conv
kernels have different shapes — 2x12, 2x1, 15x12, 9x1). Sigmoid output is in a
non-degenerate range [0.15, 0.99] post-load (random init would collapse near
0.5), and downstream AUROC is positive on both strata, confirming load
succeeded.

The dot-prefixed checkpoint files were copied to `ckpt.data-00000-of-00001` /
`ckpt.index` in `/tmp/DeepImmuno/models/cnn_model_331_3_7/` so
`tf.train.load_checkpoint('.../ckpt')` resolves them. Originals untouched.

## Issue 2 — peptide length restriction

DeepImmuno-CNN supports **9-mers and 10-mers only** (FAQ #1 in upstream
README). Encoding pads 9-mers to 10 by inserting `-` between positions 5 and
6. All 319 ITSNdb rows in the bundle are 9- or 10-mers, so **no peptides were
skipped for length**.

## Issue 3 — HLA format and rescue

Bundle uses `HLA-A*02:01` (with colons); DeepImmuno expects `HLA-A*0201` (no
colons). Stripped colons before lookup. DeepImmuno's pseudo-sequence table
covers only 62 alleles; for any HLA not in the table, DeepImmuno's built-in
`rescue_unknown_hla` substitutes the nearest-numerically subtype (same locus,
nearest field-1, then nearest field-2). **47 / 319 rows** were rescued; the
largest single rescue group is `HLA-A*2601` -> `HLA-A*2402` (n=30, both
A*-locus, no A*26 in DeepImmuno's table). Full rescue map in
`_run_meta.json`. **No rows were dropped** for unsupported HLA — the rescue
behaviour is part of the original tool, but it does mean the A*26:01 cohort
(n=30, half of `in_master` False's positives spread across alleles) is being
scored under an A*24:02 model.

## Test set (n=319)

| stratum | n | n_pos |
|---|---:|---:|
| ITSNdb_combined | 319 | 136 |
| ITSNdb_no_overlap (in_master False) | 106 | 33 |
| ITSNdb_in_master (in_master True) | 213 | 103 |

This matches the wave3 MHCflurry sibling run on the same bundle. (Spec said
"n=311 / 208 / 103"; live bundle file has 319 / 213 / 106 — same slight
mismatch the MHCflurry run also encountered. Used the live bundle counts.)

## Results

```
testset             in_master   n   n_pos   AUROC   95% CI
ITSNdb_combined     any         319 136     0.702   [0.643, 0.763]
ITSNdb_no_overlap   False       106  33     0.579   [0.464, 0.698]
ITSNdb_in_master    True        213 103     0.803   [0.741, 0.856]
```

**Δ (in_master - no_overlap) = +0.224 AUROC.** This is a large leakage gap:
DeepImmuno performs near-chance (0.58, CI lower bound 0.46) on the truly
external ITSNdb subset and only achieves its published-style ~0.80 AUROC on
peptides whose pHLA pairs overlap the public training pool — exactly the
pattern the leakage-aware re-benchmarking paper is designed to expose.

## Training-overlap caveat

DeepImmuno was trained on IEDB-derived immunogenicity calls (`assign_posterior`
over IEDB respond/test counts). ITSNdb-Main was specifically constructed to
test on peptides outside common training pools, but ITSNdb-Val + the Bayesian
master pool overlap heavily with IEDB. The `in_master` AUROC of 0.80 should
therefore be read as in-distribution / partially in-training, not as a fair
external generalization estimate. The `no_overlap` 0.58 (CI lower bound below
0.5) is the headline number for the paper.

## Runtime

Total 8.5 s on CPU (no CUDA available); pure CNN forward pass on 319 pairs is
0.32 s. Lightweight model (one of the smallest in the wave3 sweep).

## Files

- `run_score.py` — scoring script
- `predictions.tsv` — peptide / hla / label / in_master / score_deepimmuno / source / split (n=319)
- `auroc_summary.tsv` — three strata with bootstrap 95% CI
- `_run_meta.json` — runtime, filter counts, HLA rescue map
