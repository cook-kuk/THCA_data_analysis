# TransPHLA-AOMP scoring on ITSNdb (leakage-stratified)

## Install path
- **TransPHLA-AOMP** (Wu 2024 SJTU/A96123155 lab) cloned from
  https://github.com/a96123155/TransPHLA-AOMP into `/tmp/TransPHLA`.
- Used the released checkpoint `TransPHLA-AOMP/pHLAIformer.pkl` (n_layers=1, n_heads=9, fold=4 — matches model.py defaults).
- No fallback needed (TransPHLA-AOMP loaded and ran on first attempt; sibling
  fallbacks MixMHCpred / MHCnuggets were not used).

## Inference path
- We did **not** invoke `pHLAIformer.py` CLI — its top-level `from scipy import interp`
  is broken under modern SciPy 1.17.x (interp removed) and its CLI requires
  paired peptide/HLA fasta files. Instead, we imported `Transformer`,
  `read_predict_data` from `model.py` and ran inference on a dataframe of
  (HLA, peptide) pairs via the model's own `MyDataSet` / DataLoader path.
- We monkey-patched `scipy.interp = numpy.interp` before the import to satisfy
  any transitive imports.
- We replaced `model.eval_step` with a small `eval_step_safe` (run script) that
  skips their `transfer()` post-processing — `transfer()` has a
  `[0, 1][bool_array_element]` indexing bug under modern NumPy. Probabilities
  (`softmax(...)[:,1]`) are unchanged; we only avoided the buggy thresholding
  step, which we don't need for AUROC anyway.

## HLA format conversion
- Input bundle uses `HLA-A*02:01` style.
- TransPHLA-AOMP's `common_hla_sequence.csv` uses the same format → no conversion
  needed (`read_predict_data` merges by `HLA` column directly).

## Allele coverage
- Unique alleles in ITSNdb subset: 35
- Unsupported by TransPHLA's bundled HLA pseudo-sequence DB (112 alleles): 1
  - **HLA-B*41:02** (1 peptide: REFDKIELAY, label=1, ITSNdb_main, in_master=True)
- n_peptides skipped due to unsupported HLA: 1 / 319

## Scoring
- `score_transphla` = `softmax(model_logits, dim=1)[:, 1]` — class-1 probability
  (binder probability), higher = more likely to bind/be presented.
- Peptide lengths in subset: 9-mer (n=270), 10-mer (n=49) — both within
  TransPHLA's trained range (8-15-mer); cut_peptide=False semantically since
  our peptides are already ≤ 15.
- Batch size 1024 (single batch covered all 318 pairs).

## Runtime
- 9.3 s wall-clock (CPU, single batch) for 318 (peptide, HLA) pairs.

## Results
- **ITSNdb_combined** (n=318, n_pos=135): AUROC = 0.580 [0.515, 0.642]
- **ITSNdb_no_overlap** (n=106, n_pos=33): AUROC = 0.555 [0.440, 0.673]
- **ITSNdb_in_master** (n=212, n_pos=102): AUROC = 0.583 [0.508, 0.653]
- Δ (in_master − no_overlap) = +0.028 AUROC.

## Training-overlap caveat
- TransPHLA-AOMP was trained on a curated peptide-HLA-I binding dataset
  derived from IEDB, IEDB-Anthem, and other publicly available sets
  (Wu et al. 2024). The `in_master` flag in our bundle reflects overlap with
  the assembled neoantigen master training pool used elsewhere in this
  benchmark; it does NOT directly correspond to the TransPHLA training set,
  but a substantial fraction of `in_master` peptides (PRIME / NetMHCpan
  training data lineage) are likely also represented in the TransPHLA
  training corpus, so the small +0.028 AUROC gap on `in_master` should be
  read as an upper bound on novel-peptide generalization, not a clean
  hold-out estimate. Honest no-overlap performance (0.555) is essentially
  at chance for this immunogenicity-graded ITSNdb test, consistent with
  the broader pattern that pHLA binding predictors do not separate
  immunogenic from non-immunogenic binders.
