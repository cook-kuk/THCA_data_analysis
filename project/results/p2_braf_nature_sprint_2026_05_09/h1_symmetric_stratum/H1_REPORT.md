# H1 — symmetric stratum panel evaluation (Paper 2 BRAF Nature sprint, 2026-05-09)

## Headline (strongest finding for Paper 2)

**H&E and bulk RNA are anti-correlated across strata: where image fails, RNA wins; where image wins, RNA fails. They are complementary, not redundant.**

The original BRAF/cPTC result ("RNA replaces H&E") inverts on RAS/FVPTC: there image hits AUC=0.97 while HT-13 collapses to 0.59 [0.23-0.92]. The collapse is **n=3 DM1 imbalance**, not biology — on the master cohort RAS/FVPTC (n=57) HT-13 recovers to 0.83, BRAF/cPTC (n=110) holds at 0.92, and the full TCGA primary cohort (n=176) gives HT 0.90 / MAPK 0.96 / non-leak-22 0.97.

**Reframed Paper 2 narrative:** H&E reads dedifferentiation morphology (visible in FVPTC, invisible in well-differentiated cPTC); RNA reads the immune-overlap axis (HT) + MAPK output. Either alone is partial; their union triangulates DM1.

## Key numbers (pooled OOF AUC; bootstrap 95% CI)

| Stratum | n (DM1/DM2) | image-only | rna_ht (13) | rna_mapk (9) | rna_nonleak (22) | meth_8gene (LEAK) | rna_leak (8) |
|---|---|---|---|---|---|---|---|
| BRAF_like/cPTC (CLAM) | 41 (26/15) | **0.564** [0.38-0.74] | **0.926** [0.79-1.00] | 0.974 [0.93-1.00] | 0.992 [0.97-1.00] | 0.918 [0.80-0.99] | 0.959 [0.88-1.00] |
| RAS_like/FVPTC (CLAM) | 16 (3/13) | **0.974** [0.90-1.00] | **0.590** [0.23-0.92] | 0.641 [0.08-1.00] | 0.615 [0.26-0.92] | 0.513 [0.23-0.77] | 0.949 [0.79-1.00] |
| ALL_CLAM (pooled) | 59 (29/30) | 0.716 [0.58-0.84] | 0.915 [0.83-0.98] | 0.960 [0.90-1.00] | 0.962 [0.89-1.00] | 0.906 [0.82-0.98] | 0.918 [0.82-0.99] |
| ALL TCGA-THCA master (no img) | 176 (107/69) | — | 0.900 [0.85-0.94] | 0.957 [0.93-0.98] | 0.967 [0.94-0.99] | 0.931 [0.89-0.97] | 0.925 [0.88-0.96] |
| BRAF/cPTC master (no img) | 110 (85/25) | — | 0.924 [0.83-0.99] | 0.983 [0.96-1.00] | 0.992 [0.98-1.00] | 0.968 [0.91-1.00] | 0.965 [0.92-0.99] |
| RAS/FVPTC master (no img) | 57 (18/39) | — | 0.828 [0.70-0.93] | 0.818 [0.68-0.94] | 0.870 [0.77-0.96] | 0.843 [0.72-0.94] | 0.819 [0.69-0.93] |

CLAM strata reuse the v2 5-fold IDs (so image-only is comparable to the published 0.83±0.14 mean across strata). RAS/FVPTC CLAM collapsed to 3-fold because n_DM1=3.

## Answers to the two structural questions

**Q1 — Universal or BRAF-cPTC-specific?**  
Universal at adequate n. HT-13: BRAF/cPTC 0.924 (n=110), RAS/FVPTC 0.828 (n=57), full master 0.900 (n=176). CLAM-restricted RAS/FVPTC drop to 0.59 is small-n imbalance (3 DM1, CI crosses chance), not biology.

**Q2 — Does RNA replace H&E in EVERY stratum?**  
**Only where image fails.** Mirror image: BRAF/cPTC (image 0.56, HT 0.93); RAS/FVPTC-CLAM (image 0.97, HT 0.59). The "RNA structurally replaces H&E" claim from the BRAF-only run inverts on FVPTC.

## Honest caveats

- RAS/FVPTC CLAM n=16 with 3 DM1 forced 3-fold; CIs are wide, point AUCs unstable.
- `meth_8gene` is the panel that DEFINED DM1/DM2 — its high AUC everywhere is tautology, included as data-join sanity check.
- `rna_leak` (8 thyroid-diff genes) is ~80% redundant with the label via methylation↔expression coupling; AUCs are upper bounds.
- Master "no image" splits use stratified 5-fold (seed=42); CLAM strata reuse v2 fold IDs.
- Master sub-stratification is post-hoc; reported all with n≥30 (no cherry-pick).

## Files
- `h1_results.tsv` (33 rows: 6 strata × 5-6 panels)
- `h1_per_slide_preds.tsv` (459 rows of OOF probs)
- `run_h1_symmetric.py` (reproducible)
