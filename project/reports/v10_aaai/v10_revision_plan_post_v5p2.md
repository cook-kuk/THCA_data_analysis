# v10 AAAI paper — revision plan after v5.2 critical assessment

**Date:** 2026-04-25
**Trigger:** `reports/v5p2/v5p2_critical_assessment.md` finds the v5.1 THCA flip
(AUC_post = 0.006, DIAL = 0.494, batch_entangled in 4/5 classifiers) does
**not** survive proper per-fold LODO ComBat. Re-run numbers (THCA only):

| Classifier | AUC_post v5.1 | AUC_post v5.2 | DIAL v5.1 | DIAL v5.2 |
|---|---:|---:|---:|---:|
| LogReg_l2          | 0.006 | **0.995** | 0.494 | **0.000** |
| LogReg_elasticnet  | 0.008 | **0.995** | 0.492 | **0.000** |
| RandomForest       | 0.177 | **0.878** | 0.323 | **0.000** |

All 25 (5 cancer × 5 classifier) cells now non-flip. The "THCA-vs-rest"
asymmetry is gone. Source: `results/v5p2_fix/v5p2_lodo_comparison.tsv`.

## What this means for v10's AAAI paper

| Component | Status | Why |
|---|---|---|
| Theorem 1 (phase transition at ρ²=½) | **Survives** | Pure theory; doesn't depend on real data. |
| Lemma 1 (K=2 retention identity, K≥3 generalisation) | **Survives** | Same. |
| Proposition 2 (quantum 2^−(d−1) bound) | **Survives** | Pure theory + Qiskit numerical (no real data). |
| Synthetic benchmark (4,950 runs, ρ_YY mechanism) | **Survives & strengthens** | The flip we showed (0.764 → 0.339) was driven by *platform divergence*, not by ComBat leakage. Reproducible. |
| §5.2 "Real cancer data" — citing AUC = 0.006, DIAL = 0.494, "only THCA fires" | **BROKEN** | Falsifies under proper LODO. |
| Abstract sentence "only subtype-cohort-confounded settings trigger DIAL" | **BROKEN** | Hinges on §5.2. |
| §5.4 meta-analysis ("m_THCA ≥ 0.9999") | **BROKEN** | Cites the same v5.1 numbers. |
| §5.3 quantum empirical (VQC reproduces flip on THCA) | **BROKEN-DEPENDENT** | The "VQC also flips on THCA" finding inherits the leak. |

About 20% of the paper text needs surgery. The theory + synthetic + Qiskit
core (~75% of the contribution) is intact.

## Three viable paths

### Path A — pivot to "leakage-induced false-positive" (recommended)

**New headline:** *DIAL detects the leakage; the leakage was the bug.*

The v5.1→v5.2 contrast becomes the empirical hook. Sequence:

1. ComBat fit on full pooled data, then LODO split → AUC_post 0.006, DIAL 0.494 (THCA, LogReg-L2).
2. ComBat fit per LODO fold (proper protocol) → AUC_post 0.995, DIAL 0.000.
3. Theorem 1 + Lemma 1 explain *why* leakage produces the flip:
   the held-out cohort's distribution is implicitly preserved in the
   empirical-Bayes priors, manufacturing the platform-divergence
   condition that drives the flip.
4. The synthetic benchmark (already in §5.1) is the controlled
   demonstration. We add a "leakage knob" axis: fit ComBat on (i) full
   pooled, (ii) train-only, (iii) train-plus-leak fraction λ ∈ [0,1].
   Predict: DIAL grows monotonically in λ.
5. Proposition 2 stands as the quantum mitigation lemma.

**Pros:** stronger AAAI fit (methodology > biology); honest about the bug;
real-data anchor is concrete (the fix); we already have all the numbers.
**Cons:** requires running the leakage-knob experiment (~3–4 hr of
compute) and rewriting Abstract / §1 / §5.2 / §5.4 / §6.

### Path B — trim real-data section

Keep theorem + Lemma 1 + synthetic + Qiskit. Drop §5.2, §5.4 entirely.
Add a one-paragraph note acknowledging that the real-cancer demonstration
in the companion paper has been retracted pending v5.2 reanalysis.

**Pros:** smallest-edit option; honest. **Cons:** the empirical hook is
now synthetic-only; reviewers may ask "where's the real-world impact?";
the paper feels truncated.

### Path C — pure calibration / simulation paper

Build a controlled simulation where the per-fold-vs-pooled leak rate is
tunable, characterise DIAL's false-positive rate as a function of leak.
THCA becomes a worked example, not the headline.

**Pros:** matches the v5.2 critical-assessment author's recommendation
(option 2); cleanest scientifically. **Cons:** longest revision; loses
the cancer-genomics narrative entirely; risks making the paper feel
narrow for AAAI.

## Recommendation: **Path A**

It's the right balance for AAAI: the methodology framing is stronger than
the original biology framing, the v5.1→v5.2 contrast is a concrete and
reproducible empirical anchor, and we keep the cancer-data story for
context without depending on a falsified claim.

## Section-by-section edit map (Path A)

| Section | Edit | Source |
|---|---|---|
| Title | Replace with: *"Information leakage in batch-corrected LODO evaluation: a phase-transition characterisation and quantum mitigation."* | new |
| Abstract | Rewrite. Lead with the leak; cite the 0.006 → 0.995 contrast on THCA-LogReg-L2; keep the synthetic + Qiskit numbers. | new + v5p2_lodo_comparison.tsv |
| §1 Intro | Reframe motivation: ComBat-with-LODO is widely used; we identify a leakage failure mode that DIAL flags. | new |
| §2 Preliminaries (the flip symptom) | Cite both v5.1 (leaky) and v5.2 (proper) numbers as the contrast that motivates the analysis. | v5p2_lodo_comparison.tsv |
| §3–§4 Theory | Unchanged. | — |
| §5.1 Synthetic | Add the leakage-knob axis (λ ∈ [0,1]); show DIAL monotone in λ. ~3–4 hr of new compute. | new script v10_leakage_knob.py |
| §5.2 Real cancer | Rewrite as "v5.1 vs v5.2 contrast under LODO ComBat". Both numbers shown. The flip is a *protocol artifact*, not a biology finding. | v5p2_lodo_comparison.tsv |
| §5.3 Quantum empirical | Trim — the VQC-on-THCA-flips claim was inherited; replace with classical Qiskit synthetic. | existing |
| §5.4 Meta-analysis | Drop. The m-value claim depends on the v5.1 numbers. | — |
| §6 Discussion / §7 Conclusion | Reframe as: leakage-aware audit metric. Mention proper-LODO ComBat is the prescribed fix. | new |
| §8 Limitations | Add: "the v5.1-style protocol is not a strawman — it appears in published code (`v5p1_common.py:167`); the failure mode is real and common." | — |

Estimated effort: **8–12 hr writing + 3–4 hr compute** (mostly the
leakage-knob synthetic). The Theorem-1 / Proposition-2 / quantum-Qiskit
core needs no edits.

## New title options under Path A

| | Title | Posture |
|---|---|---|
| **A.1 (recommended)** | *Information-leakage-induced label flip in LODO-evaluated batch-corrected biomedical ML: a phase-transition characterisation and quantum mitigation bound.* | Methodology |
| **A.2** | *DIAL detects the leak: a direction-invariant audit for batch-correction protocols, with synthetic and quantum validation.* | Methodology, more accessible |
| **A.3** | *When batch correction creates the bug it claims to detect: a phase-transition theory of LODO leakage.* | Provocative |

## Companion-paper coordination

The Bioinformatics companion (`reports/v5/v5p1_paper.tex`) cites the
same v5.1 THCA numbers and is now also broken. Coordinate revisions:
either (a) retract that submission and resubmit aligned with v5.2, or
(b) split — let the AAAI paper carry the methodology story, the
Bioinformatics paper carries the v5.2 protocol-fix story.

## ★ UPDATE 2026-04-25 22:50 — leakage-knob 2-D sweep results

After running `notebooks_or_scripts/v10_leakage_knob.py` (264 LODO runs
on a 2-cohort synthetic with `inmoose.pycombat_norm` + Y covariate),
we have a much sharper picture:

**Mean DIAL by (ρ_YY, λ):**

| ρ_YY \ λ | 0.0 | 0.2 | 0.4 | 0.6 | 0.8 | 1.0 |
|---:|---:|---:|---:|---:|---:|---:|
| **−0.7** | **0.35** | 0.29 | 0.30 | **0.37** | 0.17 | 0.33 |
| **−0.3** | 0.19 | 0.20 | 0.23 | 0.12 | 0.16 | 0.14 |
|  0.0 | 0.02 | 0.14 | 0.07 | 0.07 | 0.05 | 0.02 |
|  +0.3 | 0.01 | 0.04 | 0.01 | 0.00 | 0.00 | 0.00 |
|  +0.7 | 0.00 | 0.01 | 0.00 | 0.00 | 0.00 | 0.00 |
|  +1.0 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |

**What the table says, in plain English:**

1. **DIAL is dominated by ρ_YY (platform divergence), not by λ (leak).**
   Reading down each column at fixed λ, DIAL grows monotonically as
   ρ_YY decreases below 0. Reading across each row at fixed ρ_YY,
   DIAL is roughly flat.
2. **The leak alone does not create flips.** At ρ_YY = +0.7 or +1.0,
   DIAL = 0 across all λ — even with full v5.1-style pooled fit.
3. **The leak modulates flips that already exist due to divergence.**
   At ρ_YY = 0, DIAL bumps from 0.02 at λ=0 to 0.14 at λ=0.2, then
   returns to ~0.02. Mild, non-monotonic.

### Implication for the v5.1 → v5.2 THCA contrast

The v5.2 critical assessment attributes the entire 0.006 → 0.995 swing
to the leak. The synthetic suggests this is **partially correct but
incomplete**: the THCA flip likely had two co-acting causes —

| Cause | Evidence | v5.1 | v5.2 |
|---|---|---|---|
| Real platform divergence (RNA-seq vs old microarray) | THCA: TCGA-THCA + GSE27155, very different platform vintages | active | active |
| Leak amplification | v5.1 fits ComBat on full pooled before LODO | active | removed |
| v5.2's specific protocol (no β·Y reinjection on held-out) | v5p2_combat_lodo.py:81 | n/a | shields the test cohort from the divergent β |

So the v5.2 fix is doing two things at once: (a) removing the leak,
(b) refusing to project the divergent β onto the test cohort.
**(b) is what actually rescues the AUC**, not (a).

### Revised paper story (still recommended Path A)

- **Headline:** DIAL fires on platform-divergent cohorts. Pooled-fit
  ComBat with covariate β leakage *modulates* the flip; the
  v5.2 protocol fix demonstrates the modulator's role.
- **Theorem 1 / Lemma 1** unchanged: ρ²_eff > threshold ⇒ flip.
- **§5.1 synthetic** unchanged (4 950-run ρ_YY sweep, real ComBat).
- **§5.2 leakage knob (NEW):** the 2-D ρ_YY × λ sweep above; argues
  that the leak modulates but does not cause flips.
- **§5.3 real cancer (REWRITTEN):** present BOTH v5.1 and v5.2 numbers
  on THCA. State that v5.2's protocol shields against both leakage and
  divergent β reinjection on held-out cohorts. Cite v5p2_combat_lodo.py.
- **§5.4 quantum** (rename from §5.3 in current draft) unchanged.
- **§5.5 meta-analysis (DROPPED).**

This framing is *more interesting and harder to scoop* than the
original "THCA flips, others don't" story. Net effect on the AAAI
contribution: **stronger.**

### Files

- Script: `notebooks_or_scripts/v10_leakage_knob.py`
- TSV:    `results/v10_aaai/synthetic_benchmark/leakage_knob.tsv` (264 rows)
- HTML:   `reports/html/figs_interactive/v10/leakage_knob.html` (heatmap)

## Immediate next actions (in order)

1. Add `\textbf{[v5.2 RETRACTED]}` margin notes on every line of the v10
   `aaai_paper_v0.tex` that quotes the broken numbers. (Done in the
   header notice; per-section markers still missing.)
2. Build `notebooks_or_scripts/v10_leakage_knob.py` — the
   λ-parameterised pooled-vs-per-fold ComBat sweep.
3. Rewrite Abstract + §2 (the flip symptom) + §5.2 against the v5.2
   numbers. Send a complete redline to me before global edits.
4. Drop §5.4 and tighten §5.3 to remove the inherited VQC-on-THCA claim.
5. Confirm with companion-paper coordination.

## TODO markers added by this plan

- [ ] Decide between Path A (recommended), B, or C.
- [ ] If Path A: write `v10_leakage_knob.py` and run.
- [ ] Update the dashboard supplementary
  (`reports/html/pages/v10_aaai_technical.html`) §4 cross-cancer table
  — currently shows the broken v5.1 row.
- [ ] Coordinate with `reports/v5/v5p1_paper.tex` revision.

Bottom line: the paper isn't dead; the empirical anchor moved. Path A
*improves* the AAAI fit. Recommend committing to it within the next 48 h
so the revision and the companion-paper retraction stay in sync.
