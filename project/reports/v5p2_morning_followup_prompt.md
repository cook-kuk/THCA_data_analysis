# v5.2 Morning Follow-up — What's Left After the Retraction Landed

> **Companion to** `reports/v5p2_critical_fix_prompt.md`. That earlier
> prompt was a *pre-flight* checklist for running the v5.2 ComBat-LODO
> fix. It is now mostly obsolete: the fix ran overnight, the verdict
> printed, and 16 documents were annotated with v5.2 retraction
> markers. **Read this file instead** for what's actually still open
> tomorrow morning.

## Verdict (one line)

THCA "4/5 batch_entangled" claim **collapsed** under leak-safe LODO.
All 5 THCA classifiers → `true_biology`, AUC_post ∈ [0.75, 0.99],
DIAL = 0.000. ΔDIAL on the headline LogReg_l2 row = −0.494.
Detail: `reports/v5p2/v5p2_critical_assessment.md`,
`results/v5p2_fix/v5p2_impact_analysis.md`.

## What's already done (no action needed)

- ✅ **v5.2 fix ran**: `results/v5p2_fix/v5p2_dial_proper_lodo.tsv`
  (25 rows × 5 cancers × 5 classifiers).
- ✅ **v5p1_paper.tex retraction header + abstract pivot** to
  *"covariate-aware ComBat fit on pooled train+test in LODO produces
  flip artefacts"* methodology framing.
- ✅ **16 / 16 v8 / v8.1 documents** carry SURVIVES / FALLS / SUPERSEDED
  markers (run `grep -lE 'v5\.2 (SUPERSEDED|retraction|SURVIVES|FALLS)'
  reports/v8/* results/v8p1_rigor/*/*.md results/v8_statgen/*.md` to
  verify).
- ✅ **Public dashboard v4a page** has a v5.2 retraction banner at the
  top linking to `research_reasoning.html` (the v5.2 reasoning
  timeline, 39 KB, exists).
- ✅ **Topnav** has *Reasoning* and *v12 Biology* tabs surfaced.

## What's open and gated on a single morning of clean work

### Track 1 — paper 1 (DIAL methodology, formerly Bioinformatics)

The retraction header + abstract are in `reports/v5/v5p1_paper.tex`.
**The body still has v5.1 numbers.** Three options, in preference
order:

1. **Strip and rewrite the body around the v5.2 narrative** (1 working
   day). Section "DIAL results (REAL, LODO)" rewritten to lead with
   the v5.1 → v5.2 ΔDIAL = −0.494 leak-fix table, supported by the
   leak-safe per-fold LinearComBat fit/transform decomposition (`v5p2_combat_lodo.py`).
   Section "Cancer-specificity" deleted (no asymmetry to report).
   Section "Robustness (v8/v8.1)" trimmed to the SURVIVES axes only
   (raw STAR DESeq2, full LMM, BRS 3-cohort, gold slate).
2. **Mark every v5.1 number in the body as `[v5.1 leak-artefact]`**
   inline and add a `Section{Postscript: v5.2 leak-safe re-evaluation}`
   at the end (4 hours). Keeps the original argument legible as a
   historical document; cleaner for archival than for submission.
3. **Pull the paper entirely** and start a fresh manuscript from the
   v5.2 critical assessment as the seed (1–2 working days). Cleanest
   final product but loses the v5.1 body's cohort-availability,
   harmonisation, and theory sections that are still valid.

**Recommendation: option 1.** It's the honest path that maximises
existing infrastructure. The v8.1 SURVIVES axes provide methodology
strength and the v5.2 leak-fix is a publishable methodological
contribution.

### Track 2 — paper 2 (TROP2 / biomarker / drug-repurposing standalone)

Draft already exists at `reports/v13/v13_trop2_standalone_draft.md`.

The biomarker spine that survives v5.2:
- **8/8 druggable targets** retained under both raw-count DESeq2 (v8.1
  Task A) and full-gene LMM (Task B)
- **1 786-gene gold slate** (Task D, dual-validated)
- **3-cohort cross-platform direction validation** (Task E:
  TCGA × GSE27155 × GSE33630 × GSE29265, 67.3 % full / 88.3 % ≥2/3)
- **CCLE negative-result** (Task F) honestly reported

This is publishable as a clinical / drug-discovery paper independent
of the DIAL methodology paper. It does **not** depend on the v5.1
LODO ComBat protocol; v5.2 changes nothing here.

**Recommendation**: read the existing draft, align it with the v8.1
biomarker numbers, target a clinical journal. ~1 working day.

### Track 3 — downstream halt audit (per v5p2 critical assessment §Recommendation)

The v5p2 assessment explicitly recommends halting v15 / v10 / v13
work that quotes the leak claim. Quick audit:

```bash
# Scan v15/v10/v13/v9 papers and reports for stale v5.1 quotes:
grep -rEln 'DIAL.{0,5}=.{0,5}0\.49|batch_entangled.*4/5|auc.post.{0,5}=.{0,5}0\.006' \
  reports/v15* reports/v10* reports/v13* reports/v9* 2>/dev/null
```

If the grep returns hits without v5.2 markers, those need
- a retraction header (option 2 from Track 1, half-hour each), or
- removal of the offending citation.

### Track 4 — v8 S5 quantum re-run under v5.2 LODO (optional)

The S5 quantum DIAL grid was computed using the same leaky ComBat as
v5.1 (`combat_preserve_pcs` in `v8_quantum_dial.py` calls
`_combat_preserve` on the full pooled PC matrix). Under v5.2 the
THCA/VQC DIAL = 0.386 result is probably also a leak artefact.

If you want a v5.2-clean quantum result (because reviewers will ask):
1. Patch `v8_quantum_dial.py` to use `LinearComBat.fit/.transform` per
   LODO fold instead of the pooled `_combat_preserve`.
2. Re-run on THCA only (the only cell where the v5.1 → v5.2 result
   is expected to change; non-THCA cells should be unchanged at
   DIAL = 0).
3. Update `v8_quantum_comparison.tsv` and S5 in supplementary.

Wall time: ~30 min. Skip unless paper 1 reviewers specifically
request quantum-paradigm-invariance under v5.2.

## Definitely-do-not list

- Do **not** quote `DIAL = 0.494` or `4/5 batch_entangled` in any new
  draft, slide, blog post, or grant.
- Do **not** rebuild figures showing the v5.1 flip without a
  "v5.1-leak-artefact" caption.
- Do **not** chase the v5.1 result by tweaking ComBat parameters or
  the LinearComBat transform — `auc_pre = 0.99 without correction`
  on THCA is the diagnostic: the original "flip" was leak-driven, not
  biology being entangled with batch.
- Do **not** start v6 scGPT / v15 NeurIPS / v9 Ideker work until paper
  1 body is rewritten (or at least until you've decided which option
  from Track 1).

## Decision gate for the morning

Pick ONE track to start on, in the morning, after coffee:

- [ ] **Track 1 option 1**: rewrite paper 1 body (1 day, highest
      reviewer-credibility outcome)
- [ ] **Track 1 option 2**: mark every v5.1 number inline (4 hours,
      preserves the historical paper)
- [ ] **Track 2**: align paper 2 draft with v8.1 biomarker numbers
      (1 day, parallel-publishable to clinical journal)
- [ ] **Track 3**: downstream halt audit (half day, prevents future
      mis-citation)

Track 4 is opt-in after one of 1–3 is in flight.

## File map (canonical reference)

- `results/v5p2_fix/v5p2_dial_proper_lodo.tsv` — 25-row leak-safe DIAL
- `results/v5p2_fix/v5p2_lodo_comparison.tsv` — v5.1 vs v5.2 side-by-side
- `results/v5p2_fix/v5p2_impact_analysis.md` — per-row deltas
- `reports/v5p2/v5p2_critical_assessment.md` — verdict + recommendation
- `reports/v5p2/table1_redesigned.tex` — Table 1 schema fix (audit fix F2)
- `reports/v5/v5p1_paper.tex` — retraction header + abstract pivot in;
  body still v5.1
- `reports/v8/v8p1_FINAL_STATUS.md` — survives/falls scorecard + repair plan
- `reports/v8/v8p1_rigor_summary.md` — v5.2 retraction header
- `reports/html/pages/v4a_dial_cross_cancer.html` — public-page v5.2 banner
- `reports/html/pages/research_reasoning.html` — v5.2 reasoning timeline
  (39 KB, exists)
- `notebooks_or_scripts/v5p2_combat_lodo.py` — leak-safe `LinearComBat`
  class (`fit` / `transform` semantics)
- `notebooks_or_scripts/v5p2_dial_proper.py` — leak-safe DIAL pipeline

Pick a track. Sleep well.
