# v5.2 — Critical Fix: ComBat-on-train-only LODO + Table 1 schema

> ⚠️ **OBSOLETE — read `reports/v5p2_morning_followup_prompt.md` first.**
> The fix this prompt was written to drive **already ran overnight**
> (results: `results/v5p2_fix/v5p2_dial_proper_lodo.tsv`,
> verdict: `reports/v5p2/v5p2_critical_assessment.md`). The THCA flip
> collapsed entirely (DIAL = 0.000 across all 5 classifiers, AUC_post
> ∈ [0.75, 0.99]). The "what's left tomorrow morning" checklist is
> in the followup prompt linked above. The body of this file is
> retained for traceability of the original Stage A1–A4 spec.

**For tomorrow morning's clean session.** Paste the entire spec below into
a fresh agent / tmux. Do not interleave with v6 / v10 / v12 / v13 / v15
work — this is the single highest-priority item before the Bioinformatics
paper can submit.

> **NOTE (2026-04-25 evening).** A `v5p2_dial_proper.py` is already in
> flight in tmux session `v5p2_fix` (8 python procs, ≈ 1 hr in). Tomorrow
> morning, **first** check whether that run finished and produced output
> at `results/v5p2_proper_lodo/` — if yes, jump straight to A2/A3/A4. If
> no, kill it, read the script, fix any bugs, and restart from A1.

---

## Why this matters

The v5.1 LODO pipeline applies `pycombat_norm` to the **entire** harmonised
matrix `X` (including the held-out cohort), then performs LODO CV on the
ComBat-corrected output. This is **leakage-equivalent** to using test-fold
information during preprocessing: ComBat's empirical-Bayes priors absorb
the held-out batch's mean and dispersion before the train fold ever sees
the data.

Proper LODO requires:

1. ComBat **fit on the training fold only** (one or more cohorts).
2. ComBat **transform applied to the held-out cohort** at inference time
   using only the training-fold parameters.
3. AUC measured on the held-out cohort exclusively.

**The risk**: under proper LODO, the THCA `batch_entangled` 4/5 flip may
shrink, disappear, or move to fewer classifiers. **Until A4 is run, do
not submit the Bioinformatics paper.**

---

## Stage A1 — ComBat fit-on-train-only

`notebooks_or_scripts/v5p2_dial_proper.py` (already exists; verify it
implements the spec below or rewrite if not):

```python
# Pseudocode for the critical inner loop.
from inmoose.pycombat import pycombat_norm

def compute_dial_proper(X, Y, B, clf_factory):
    cv = LeaveOneGroupOut().split(X, Y_bin, groups=B)
    pre, post = [], []
    for tr, te in cv:
        # Pre: classifier on raw train, eval on raw test
        c = clf_factory().fit(X[tr], Y_bin[tr])
        pre.append(roc_auc_score(Y_bin[te], c.predict_proba(X[te])[:, 1]))

        # Post: ComBat fit on TRAIN ONLY, transform applied to test
        # via the training reference distribution.
        # `pycombat_norm` returns a corrected matrix; for proper LODO
        # we need the per-batch shift parameters (gamma_star, delta_star)
        # from the training fit and apply them to the test batch.
        # If pycombat_norm exposes train/transform separately use that;
        # otherwise build the corrected train matrix, learn the residual
        # cohort shift between train and test in the *uncorrected* space,
        # and apply.
        X_tr_post = pycombat_norm(X[tr].T, batch=B[tr],
                                  covar_mod=one_hot(Y[tr])).T
        # Apply train-fit ComBat parameters to held-out cohort:
        X_te_post = apply_pretrained_combat(X[te], B[te],
                                            train_params=...)
        c = clf_factory().fit(X_tr_post, Y_bin[tr])
        post.append(roc_auc_score(Y_bin[te], c.predict_proba(X_te_post)[:, 1]))
    return ...
```

`apply_pretrained_combat` is the work item. Three options, in
preference order:

1. **Use `inmoose.pycombat`'s `pycombat_apply` if it exposes one**
   (check the v0.4+ API; `pycombat_norm` may now return parameters).
2. **Implement manually**: extract `gamma_star`, `delta_star` from a
   patched `pycombat_norm` run on train, then compute
   `X_te_corrected = (X_te - alpha - gamma_star[B_te]) / delta_star[B_te] + alpha`
   per the Johnson 2007 formulation. About 30 lines of numpy.
3. **Train-only ComBat (no covariate)**: degenerate fallback if (1) and
   (2) fail. Acceptable for a sanity check but not for the published
   number — document if used.

Run on all 5 cancers × 5 classifiers × 2-fold LODO (THCA) and the
existing fold counts for SKCM/LGG/LUAD/COAD.

**Output**:
- `results/v5p2_proper_lodo/v5p2_dial_all_cancers.tsv`
- `results/v5p2_proper_lodo/v5p2_dial_THCA.tsv` (focal)
- `logs/v5p2_dial_proper.log`

## Stage A2 — Table 1 schema fix

The current Table 1 in `reports/v5/v5p1_paper.tex` (lines 19–35) and in
the dashboard reports `n` (samples in expression matrix) but conflates
three distinct counts:

- **n_RNA**: samples in raw RNA-seq / microarray matrix
- **n_labels**: samples with BRAF/RAS or equivalent label
- **n_usable**: intersection of n_RNA × n_labels × passing v5.1 QC

For example THCA TCGA: n_RNA = 568, n_labels = 351, n_usable = 351.
For GSE27155: n_RNA = 99, n_labels = 41, n_usable = 41 (in v5.1 only
the labelled 41 enter LODO).

**Action**: split the column in Table 1 of `v5p1_paper.tex` into three
columns (n_RNA / n_labels / n_usable) and propagate the same schema to:

- `reports/html/pages/v8_statgen_supplement.html` artefact-index card
- `reports/v8/v8_supplementary.md` S3/S6 cohort tables
- Any dashboard tile that reports a per-cohort `n`.

Generate the canonical numbers from
`results/v5/v5p1_cohort_availability.tsv` (already exists per
`reports/v5/v5p1_table1_summary.tsv`); cross-check with
`results/v8p1_rigor/c_expanded_lodo/cohort_label_attempt.tsv` for the
failed-cohort rows.

## Stage A3 — downstream propagation

After A1 and A2 produce the new numbers, update three documents:

1. **`reports/v5/v5p1_paper.tex`** — replace DIAL numbers in the abstract,
   Section "DIAL results (REAL, LODO)", and Table 1. Add a one-line note
   in Limitations that v5.1's ComBat-on-full-matrix protocol is replaced
   by v5.2's train-only protocol.
2. **`reports/v8/v8_supplementary.md`** — S1 (ComBat-seq), S2 (meta-
   analysis), S3 (pathway DIAL), S6A (FastRNA centering): every one of
   these reports DIAL numbers that depend on the v5.1 protocol. Re-run
   with v5.2 protocol if affordable; if not, add a one-paragraph
   "v5.2 protocol caveat" stating that the qualitative direction is
   unchanged but absolute DIAL values may shift by up to ~Δ. The
   `v5p1_dial_all_cancers.tsv` output of Stage A1 is the canonical
   source.
3. **Bioinformatics short-paper draft** at
   `reports/v4_minipaper_skeleton.tex` (or wherever the live draft is —
   check `reports/v15_neurips/` or `reports/v10_aaai/` if v4 is stale).
   The draft must reflect v5.2 numbers before submission.

## Stage A4 — THCA 4/5 batch_entangled claim sanity check

**This is the gate.** Without this check the paper is unsafe to submit.

Compare v5.1 vs v5.2 on THCA only:

| Classifier        | v5.1 DIAL | v5.2 DIAL | Δ | Still `batch_entangled` (DIAL ≥ 0.3)? |
|-------------------|----------:|----------:|---|---------------------------------------|
| LogReg_l2         |     0.494 | ?         | ? | ?                                     |
| LogReg_elasticnet |     0.492 | ?         | ? | ?                                     |
| RandomForest      |     0.323 | ?         | ? | ?                                     |
| GradientBoosting  |     0.334 | ?         | ? | ?                                     |
| XGBoost           |     0.013 | ?         | ? | (already null)                         |

**Decision tree:**
- If 4 of 5 still `batch_entangled` under v5.2 → **paper claim survives.
  Update numbers, ship to Bioinformatics.**
- If 2–3 still `batch_entangled` → **partial survival.** Reframe the paper:
  the THCA flip is real but is partly an artifact of v5.1's
  ComBat-on-full-matrix protocol, and the proper-LODO flip is observed
  in a smaller subset. Still publishable, but the abstract changes.
- If 0–1 still `batch_entangled` → **claim collapses.** The v5.1 result
  is largely a methodological artifact. Pull the paper and rewrite as a
  methodology critique.

In all three branches, the v8/v8.1 robustness scaffolding **still
matters** — pathway-level analysis, BRS validation, druggable retention,
quantum paradigm — because those are independent of the LODO ComBat
protocol. The DESeq2 + LMM + 3-cohort BRS direction validation results
do **not** depend on v5.1 vs v5.2 ComBat; they remain ground truth for
the BRAF-vs-RAS biomarker slate.

## Wall-time estimate

- Stage A1: 10 min if `v5p2_dial_proper.py` finishes cleanly tonight,
  otherwise 30–60 min to debug + rerun.
- Stage A2: 30 min (mechanical schema fix across 3–4 documents).
- Stage A3: 90 min (re-running v8/v8.1 affected pieces with v5.2 numbers
  takes the most time; the markdown / tex edits are 15–20 min).
- Stage A4: 5 min (one comparison table).

**Total budget: 2 to 3 hours of clean morning work.**

## What NOT to do tomorrow morning

- Do **not** start v6 scGPT / VEGA GPU runs at the same time — these
  contend for the same GPU and v5.2 has no GPU need.
- Do **not** touch v15_neurips work — it's running `v15_stress_test.py`
  in a separate stream.
- Do **not** rewrite v8 / v8.1 supplementary text proactively — only
  update if v5.2 numbers force a change. Most of v8 is independent of
  the LODO ComBat protocol.
- Do **not** submit to Bioinformatics until A4 is decided.

## Files to touch (canonical list)

- `notebooks_or_scripts/v5p2_dial_proper.py` (verify or rewrite)
- `results/v5p2_proper_lodo/v5p2_dial_all_cancers.tsv` (output)
- `results/v5p2_proper_lodo/v5p2_dial_THCA.tsv` (output, focal)
- `reports/v5/v5p1_paper.tex` (abstract, Table 1, DIAL section, Limitations)
- `reports/v8/v8_supplementary.md` (S1/S2/S3/S6A — only if v5.2 changes
  numbers)
- `reports/v5/v5p1_table1_summary.tsv` (regenerate with three-column schema)
- `reports/html/pages/v8_statgen_supplement.html` (cohort tile labels if any)

## Stage 99 (end-of-morning)

```
=== v5.2 CRITICAL FIX COMPLETE ===
A1 protocol fix      : results/v5p2_proper_lodo/v5p2_dial_all_cancers.tsv
A2 schema fix        : reports/v5/v5p1_table1_summary.tsv (3-col schema)
A3 propagation       : v5p1_paper.tex / v8_supplementary.md / minipaper
A4 THCA flip survival:
  LogReg_l2          : v5.1 0.494  →  v5.2 ___  →  batch_entangled? ___
  LogReg_elasticnet  : v5.1 0.492  →  v5.2 ___  →  batch_entangled? ___
  RandomForest       : v5.1 0.323  →  v5.2 ___  →  batch_entangled? ___
  GradientBoosting   : v5.1 0.334  →  v5.2 ___  →  batch_entangled? ___
  XGBoost            : v5.1 0.013  →  v5.2 ___  →  (null control)
  PAPER VERDICT      : survives / partial / collapses
```

Once Stage 99 prints, the Bioinformatics submission decision is binary
and informed. Do not improvise the paper before A4.
