# v5.2 Critical Assessment

**Date:** 2026-04-25
**Author:** Seungho Cook (논문 저자)
**Triggered by:** v5.2 self-review identified two SUBMIT-BLOCKERs in v5.1.

## A1 verdict — does the central claim survive proper LODO?

**NO. Pattern broken.**

Under v5.1 (ComBat fit on full pooled data, then LODO split):
> "THCA exhibits batch_entangled in 4/5 classifiers (DIAL > 0.3) while 4 other
> cancers show 20/20 non-flip pattern."

Under v5.2 (ComBat fit on training cohort within each LODO fold; held-out
cohort centered locally — see `notebooks_or_scripts/v5p2_combat_lodo.py`):

| Cancer | v5.1 batch_entangled | v5.2 batch_entangled | Status |
| --- | ---: | ---: | --- |
| **THCA**  | **4/5** | **0/5** | claim collapsed |
| SKCM  | 0/5 | 0/5 | unchanged |
| LGG   | 0/5 | 0/5 | unchanged |
| LUAD  | 0/5 | 0/5 | unchanged |
| COAD  | 0/5 | 0/5 | unchanged |

THCA AUC_post moves from **0.006** (v5.1, complete inversion) to **0.995**
(v5.2, near-perfect generalisation) for the LogReg_l2 line that anchored
the paper. ΔDIAL = −0.494 on the same row.

The specificity evidence (4 cancers × 5 classifiers = 20 non-flips) is
preserved. But the specificity evidence was *only meaningful as a contrast to
THCA's 4/5 flip*. With THCA also non-flipping, there is no THCA-vs-rest
asymmetry to report.

## What broke

`v5p1_common.compute_dial` calls `_combat_preserve(X_raw, Y, B)` on the full
matrix **before** entering the LODO loop (`v5p1_common.py:167`). The corrected
matrix `X_post` is then split into train/test by LODO. Even though the
classifier never sees the test rows during fit, the ComBat step that produced
the rows already saw them — its γ_{gb}, δ_{gb}, α_g, β_g were estimated using
samples that the LODO fold would later hold out. With `par_prior=True` and the
covariate matrix passed in, the empirical-Bayes priors borrow strength across
all samples, including the held-out cohort.

This is a textbook information leak. It does not invalidate ComBat itself; it
invalidates the v5.1 *evaluation protocol* of ComBat under LODO.

## A2 verdict — Table 1 schema

**Fixed.** New schema separates n_rnaseq, n_maf, n_class_A, n_class_B, n_usable
per (cancer, cohort). LGG row that previously read
`n=260, n_class_A=414, n_class_B=95` (impossible) is now:

| Cohort | RNA-seq | MAF | Class A | Class B | Usable |
| --- | ---: | ---: | ---: | ---: | ---: |
| TCGA-LGG | 260 | 509 | 414 | 95 | 257 |
| DIAL cohort (post-harmonize) | 142 | — | 116 | 26 | 142 |

`n_class_A + n_class_B = 509` corresponds to the MAF universe (509 patients
have a derivable IDH-mut/IDH-wt label), which is bigger than 260 (RNA-seq
patient cap). 257 of those 509 also have RNA-seq → that is the "usable"
intersection. Harmonization further drops 257 → 142 because we keep only the
four largest TSS sub-cohorts (TSS-HT 49, TSS-DU 46, TSS-S9 25, TSS-DB 22) for
multi-batch LODO.

Files:
- `results/v5p2_fix/cohort_availability_detailed.tsv`
- `results/v5p2_fix/table1_redesigned.tsv`
- `reports/v5p2/table1_redesigned.tex`

The same audit fix was independently applied earlier in `reports/v5/v5p1_paper.tex`
under "audit fix F2" — that table already uses the corrected schema. The v5.2
redesign extends it to TSV/TeX artefacts that downstream tooling (dashboards,
supplementary materials) consumes.

## Project status

**CRITICAL.** The central paper claim does not survive the proper LODO fix.

## Recommendation

**HALT all v15 / v10 / v13 downstream work that builds on the v5.1 THCA
batch_entangled finding** (notably anything titled `*_AAAI`, `*_NeurIPS`,
`*_drug_discovery_v13` that quotes DIAL = 0.494 or "4/5 batch_entangled").
Re-frame is not sufficient because the magnitude *and direction* of the
flagship result changed.

Two viable paths forward, ordered by integrity:

1. **Retract & rewrite v5.x as a methodological paper on the leak itself.**
   The core lesson — *"batch correction must be fit per LODO fold, not
   pooled, or the audit metric will report inverted AUC as evidence of
   biology being entangled with batch"* — is publishable. Quantify
   ΔDIAL = −0.49 on THCA, show specificity evidence (no other cancer changes),
   and present the v5.2 LinearComBat fit/transform decomposition as the fix.
   This is the honest path and it is the one I recommend. The current code
   already supports this (`v5p2_combat_lodo.py`, `v5p2_dial_proper.py`).

2. **Investigate whether the leak inflates DIAL on synthetic data with a
   *known* batch confound.** If so, v5.x can be salvaged as a calibration
   study: "DIAL with leaky ComBat reports false positives at rate X." THCA
   is then a worked example, not the headline. This requires a controlled
   simulation that I have not yet run.

Do **not** attempt to keep the THCA finding by tweaking the LinearComBat
implementation (e.g., using nearest-train-batch transform) — the auc_pre is
already 0.99 *without any correction*, so the flip in v5.1 was not measuring
"ComBat removed real biology", it was measuring "leaky ComBat hallucinated a
flip". A different transform will not change auc_pre's verdict.

## Downstream files that need updating

Citations to "DIAL = 0.494" or "batch_entangled in 4/5" or "AUC_post = 0.006"
should either be retracted or annotated as v5.1-leak-artifact. Best-effort
list (verify each before editing):

```
reports/v5/v5p1_paper.tex                         (Section 4 + Table 3 numbers)
reports/v5/v5_dial_paper.tex                      (legacy DIAL paper)
reports/v8/...supplementary*.md                   (meta-analysis citing v5.1)
reports/bioinformatics_submission/...             (Section 5 if it cites v5.1 THCA)
```

I have not auto-edited these — they are outside the v5.2 scope and the
re-framing decision needs a human author.

## Reproducibility

```bash
cd /opt/thyroid-dash/project
.venv/bin/python notebooks_or_scripts/v5p2_table1_redesign.py
.venv/bin/python notebooks_or_scripts/v5p2_dial_proper.py
```

Outputs:
- `results/v5p2_fix/v5p2_dial_proper_lodo.tsv` — 25 rows (5 cancers × 5 clf)
- `results/v5p2_fix/v5p2_lodo_comparison.tsv` — side-by-side with v5.1
- `results/v5p2_fix/v5p2_impact_analysis.md` — per-row deltas
- `results/v5p2_fix/cohort_availability_detailed.tsv`
- `results/v5p2_fix/table1_redesigned.tsv`
- `reports/v5p2/table1_redesigned.tex`
- `logs/v5p2_critical_fix.log`

Total runtime: ~16 min on the project venv (numpy/sklearn/pandas only;
no inmoose dependency in v5.2 path).
