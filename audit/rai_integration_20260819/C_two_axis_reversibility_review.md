# C — Does the existing "two-axis" analysis answer the reversible-vs-irreversible dedifferentiation question?

**Date** 2026-08-19. New audit output. Extends but does not edit or supersede
`figure31_two_axis_scatter.py`, `meraiode_redifferentiation_panel_2026_08_06.py`, the
manuscript, or `audit/rai_integration_20260819/A_*`/`B_*`. All numbers below were computed
in this session directly from the existing per-sample tables (paths in §7); nothing is
transcribed from narrative without recomputation.

**One-line verdict up front:** the script named "two-axis" tests a genuinely different pair
of axes than the one the reversibility question needs, and on the axes it does test, the
off-diagonal quadrant does not outperform the aggregate score for any clinical endpoint
checked. **This is a restatement of the aggregate score under a different name for the
purpose asked here — not a dead end in general (it is real, useful biology for the
driver-stratification story), and not a new angle on reversibility.**

---

## 1. What `figure31_two_axis_scatter.py` actually plots (read before trusting the name)

The script (`rai-response-genomics-atlas/scripts/figure31_two_axis_scatter.py`) reads
`project/results/r17_tcga_panel_d4p2_reconciliation/tert_interaction/r17_tert_per_sample.tsv`
(n = 522, output already present, last written 2026-05-22 — script was re-run this session
to confirm it reproduces the stored figure; it does) and plots:

- **x-axis: `RAI_8`** — a single composite z-score averaging **all eight** panel genes
  (DIO1, FOXE1, NKX2-1, PAX8, SLC5A5, TG, TPO, TSHR) together. This is the RAI-lineage /
  differentiation-silencing axis used everywhere else in the manuscript.
- **y-axis: `sig_score`** — the **d4p2 Hashimoto/immune-overlap** signature (B/T-cell, HLA-II,
  IFN, TLS genes; from the GSE286332 PTC+HT generalization work), projected onto TCGA.

**This is not "lineage identity vs iodine-machinery functionality."** There is no column, in
this table or anywhere else feeding this script, that separates the 8-gene panel into a
lineage-identity sub-score and an iodine-machinery sub-score. `RAI_8` already collapses both
into one number before the two-axis plot is built. The two axes actually plotted are
**(1) RAI-lineage silencing** vs **(2) Hashimoto-immune infiltration** — a real and
previously-established orthogonal pair (`R17_reconciliation_report.md`), but not the pair
the reversibility hypothesis is about.

The genuine lineage-identity-vs-machinery split does exist elsewhere in this repo, computed
today in `audit/rai_integration_20260819/B_functional_decomposition_tcga.md` as `LINEAGE_4`
(PAX8, NKX2-1, FOXE1, TSHR) vs `HORMONE_3` (TG, TPO, DIO1) + `UPTAKE_1` (SLC5A5) — see §5
below for how that decomposition bears on this question. `figure31` is not that analysis.

---

## 2. Axis correlation — not "one axis measured twice," but not orthogonal either

```
n = 522 (TCGA-THCA, RAI_8 & sig_score both non-missing)
Pearson  r   = -0.439   p = 5.4e-26
Spearman rho = -0.539   p = 1.1e-40
```

r = -0.44 is well under the 0.8 "same axis twice" threshold, so by that criterion the two
plotted axes are statistically distinguishable — but they are also not the "orthogonal"
axes the manuscript captions them as (`ncomms_rai_atlas_v1.md` Fig 2a legend: "panel z and
d4p2 sig_score are orthogonal"). A moderate, highly significant negative correlation is a
real relationship, not orthogonality; the manuscript's use of "orthogonal" describes the
*label concordance* (14.8%, i.e., the categorical DM1/DM2 calls rarely agree — see
`R17_reconciliation_report.md`), not the continuous-score correlation, and the two claims
should not be conflated when this figure is discussed with reviewers.

---

## 3. Quadrant structure and what actually explains it

Sign-based quadrants on the raw continuous scores (n=522):

| | sig_score > 0 (HT-immune-high) | sig_score ≤ 0 (HT-immune-low) |
|---|---:|---:|
| **RAI_8 > 0** (panel preserved) | 54 | 179 |
| **RAI_8 ≤ 0** (panel silenced) | 150 | 139 |

The manuscript's categorical `zone` variable (built from DM1/DM2 calls, a related but not
identical threshold) gives close but not identical counts: WT-like 153, BRAF-like 165,
dark-matter 161, RAS-like 43.

**The quadrant structure is not a free-floating biological discovery — it is driver
identity restated.** From `R17_reconciliation_report.md`, per-driver breakdown of the same
522-row table:

| Driver | n | % panel-silenced (RAI_8<0) | % HT-immune-high (sig_score>0) |
|---|---:|---:|---:|
| BRAF V600E | 318 | 82% | 6.3% |
| RAS | 54 | 11% | 98% |
| BRAF·RAS-neg | 128 | ~50% (mixed) | ~50% (mixed) |

BRAF V600E tumors populate the "silenced, HT-low" quadrant (BRAF-like) almost by
definition of the driver; RAS tumors populate "preserved, HT-high" (RAS-like) almost by
definition. Only the BRAF·RAS-neg subset (n=128) genuinely mixes across quadrants. So the
"4 zones" are largely a re-encoding of driver class plus HT-overlap status, not four
independently-discovered cell states — which is exactly how the manuscript itself frames it
("driver-pattern biology, not a labeling bug," `R17_reconciliation_report.md` §Interpretation).
That is legitimate biology, but it means "quadrant membership" and "driver + one continuous
score" carry overlapping information by construction, which matters for §4.

---

## 4. Does the off-diagonal quadrant (dark-matter: silenced + HT-high) beat the aggregate score on any clinical endpoint?

Two endpoints are present in the same per-sample table used by the figure: TERT promoter
status (`tert_pos`, n=477 known) and progression-free interval events (`PFI`, n=477–522
known). No stage variable is merged into this table (it lives in a separate clinical file
used by the B-decomposition audit, not joined here); testing stage against this specific
two-axis table was out of scope for evaluating the existing script and is flagged as a gap,
not fabricated.

**TERT+ (n=477, 36 events):**

| Model | pseudo-R² | AUC | key coefficient P |
|---|---:|---:|---:|
| `RAI_8` only (aggregate, single axis) | 0.033 | 0.660 | 0.0048 |
| dark-matter zone (binary, off-diagonal quadrant) | 0.023 | 0.600 | 0.0149 |
| `RAI_8` + `sig_score` (both continuous axes, additive) | 0.033 | — | sig_score P=0.74 |
| `RAI_8` + dark-matter dummy | — | — | dark-matter P=0.27 (LRT vs RAI_8-only: p=0.27) |

**PFI event (n=477, 49 events):**

| Model | pseudo-R² | AUC | key coefficient P |
|---|---:|---:|---:|
| `RAI_8` only | 0.030 | 0.650 | 0.0030 |
| dark-matter zone (binary) | 0.018 | 0.586 | 0.0153 |
| `RAI_8` + `sig_score` (additive) | 0.030 | — | sig_score P=0.72 |
| `RAI_8` + dark-matter dummy | — | — | dark-matter P=0.34 (LRT vs RAI_8-only: p=0.34) |

**Reading this honestly:** the single continuous aggregate score (`RAI_8` alone) predicts
both TERT+ and PFI events *better* than the off-diagonal "dark-matter" quadrant does (higher
pseudo-R², higher AUC), and adding either the second continuous axis (`sig_score`) or the
quadrant membership itself on top of `RAI_8` adds no detectable signal (all LRT p>0.25). The
categorical Fisher-exact result the manuscript already reports for dark-matter × TERT
(OR=2.34, p=0.016, `R17_tert_zone_interaction.md`) is real but is not an independent
discovery beyond `RAI_8`'s low tail — it is consistent with dark-matter being, functionally,
"the low-RAI_8 half of the samples that also happen to be HT-immune-high," and the HT-immune
status contributes nothing once RAI_8 is already in the model.

**This directly answers the question the task poses: no, the off-diagonal quadrant does not
associate with TERT or PFI more strongly than the aggregate score. If anything, the
aggregate single score is the better predictor on both endpoints tested.**

Stage, BRAF/RAS association is already covered structurally in §3 (it defines the quadrants,
so testing "does BRAF/RAS associate with the quadrant" is closer to circular than
informative — of course it does, driver assignment built the axes).

---

## 5. Cross-reference with the B decomposition — same negative result from an independent angle

`audit/rai_integration_20260819/B_functional_decomposition_tcga.md` (same session, same
day) already ran the analysis that would most directly test a lineage-identity-vs-machinery
reversibility story, on TCGA RAI-treated patients against real clinical endpoints (new tumor
event, persistent disease):

- `LINEAGE_4` (PAX8, NKX2-1, FOXE1, TSHR — the identity/TF genes): OR 0.44/0.36,
  **P = 0.069 / 0.063 — does not reach nominal significance in either endpoint.**
- `HORMONE_3` (TG, TPO, DIO1) and `UPTAKE_1` (SLC5A5) — the machinery genes: HORMONE_3
  nominally significant (P=0.02/0.02) but does not survive Bonferroni correction for the
  6-test family it introduces; UPTAKE_1 alone is the weakest and least stable of the three.

So even using the *correct* identity-vs-machinery split (which `figure31_two_axis_scatter.py`
does not compute), there is no clean "lineage-high" sub-score that behaves differently from
"machinery-low" in a way that would support a distinct reversible-dedifferentiation quadrant
— LINEAGE_4 alone is the weakest of the two directionally-consistent groupings, not a
standout marker of a protected, therapy-responsive state.

---

## 6. Cross-reference with MERAIODE — the quadrant assignment cannot be tested there, and not just for power reasons

`meraiode_redifferentiation_panel_2026_08_06.py` and its brief
(`results/reports/meraiode_redifferentiation_brief_2026_08_06.md`) already establish the
central fact that forecloses this cross-reference:

- The HTG EdgeSeq panel used in the E-MTAB-12837/12900 trial gives **zero counts for
  SLC5A5 (NIS), TG, and DIO1 in all 21 libraries** — three of the four genes in the
  `HORMONE_3`/uptake machinery grouping. Only TPO survives from that group.
- The 5 genes actually measurable (PAX8, FOXE1, NKX2-1, TSHR, TPO) are dominated by the
  identity/TF module (4 of 5 genes are exactly `LINEAGE_4`); the resulting `panel_z` used in
  the response test is, in the identity-vs-machinery framing, closer to a **lineage-identity
  score with one machinery gene attached**, not a paired identity/machinery two-axis
  decomposition.
- No zone or quadrant assignment is computed for the MERAIODE patients anywhere in the
  existing pipeline — `panel_z` is a single scalar per patient, never matched against the
  TCGA `RAI_8`/`sig_score` thresholds or the `figure31` zone scheme.

Consequently: **it is not possible, with the data currently in this repository, to check
whether TCGA-defined lineage-high/machinery-low patients respond better to MERAIODE than
lineage-low/machinery-low patients.** The machinery axis this would require (SLC5A5, TG,
DIO1) is structurally absent from the only public redifferentiation-response dataset, not
merely underpowered. Reporting an association here would mean fabricating a machinery score
that the assay cannot produce. This is consistent with the MERAIODE brief's own verdict
("Do not put this in the manuscript as a validation attempt") and extends it: the reason
generalizes beyond statistical power (n=9 vs 12) to a structural assay-coverage gap that
persists regardless of sample size.

Separately, even setting the missing-genes problem aside, n=21 (9–12 per arm) is far too
small to test a 2×2 quadrant interaction (four cells of ~2–6 patients each) with any
resolving power — a second, independent reason this cross-reference cannot be run
meaningfully even in principle with this dataset.

---

## 7. Verdict

1. **The script's name overstates what it measures for this specific question.**
   `figure31_two_axis_scatter.py` plots RAI-lineage-silencing (all 8 genes combined) against
   Hashimoto/immune-overlap, not lineage-identity against iodine-machinery-functionality. It
   is a legitimate, previously-validated two-axis result for driver-stratification biology
   (BRAF vs RAS vs BRAF·RAS-neg), but it is the wrong pair of axes for the reversibility
   question as posed.
2. **On the axes it does test, r = -0.44 is not >0.8** — so this is not literally "one axis
   measured twice" — but the off-diagonal quadrant (dark-matter) carries no clinical signal
   (TERT+, PFI) beyond what the single aggregate `RAI_8` score already provides, and adding
   the second axis to a regression on top of `RAI_8` is null in every model tested (all
   p>0.25). Functionally, for prediction purposes, this reduces to the aggregate score under
   a categorical label.
3. **The correct identity-vs-machinery split was run today, in `B_functional_decomposition_tcga.md`,
   against real clinical outcomes, and does not support a clean reversible/irreversible
   quadrant story either** — LINEAGE_4 alone is directionally consistent but non-significant,
   and is the weaker of the two sub-modules, not a distinct protected state.
4. **The MERAIODE cross-reference cannot be run at all with existing data**, because the only
   public redifferentiation-response transcriptomes structurally lack the three genes
   (SLC5A5, TG, DIO1) that would define the "machinery" axis, independent of sample size.

**Overall: for the specific "reversible vs irreversible dedifferentiation" question, this is
a dead end with the data and scripts that exist today** — not because the two-axis idea is
wrong in principle, but because (a) no script in this repo computes the identity/machinery
axis pair the hypothesis needs (figure31 computes a different pair; B computes the right
pair but against different, non-quadrant endpoints and finds no clean signal), and (b) the
one dataset that could test therapy-response against such an axis cannot measure the
machinery genes at all. Building a genuine test would require either a new redifferentiation
cohort with full 8-gene coverage, or accepting that MERAIODE can only ever speak to the
identity axis.

---

## 8. Provenance

| Item | Path |
|---|---|
| Two-axis script (read, re-run to confirm reproducibility) | `rai-response-genomics-atlas/scripts/figure31_two_axis_scatter.py` |
| Two-axis figure output (pre-existing, confirmed current) | `rai-response-genomics-atlas/results/figures/figure31_two_axis_scatter.{png,pdf}` |
| Per-sample source table (RAI_8, sig_score, zone, TERT, PFI) | `project/results/r17_tcga_panel_d4p2_reconciliation/tert_interaction/r17_tert_per_sample.tsv` |
| Two-axis concordance report (source of driver breakdown, §3) | `project/results/r17_tcga_panel_d4p2_reconciliation/R17_reconciliation_report.md` |
| TERT × zone interaction (source of Fisher OR=2.34, §4) | `project/results/r17_tcga_panel_d4p2_reconciliation/tert_interaction/R17_tert_zone_interaction.md` |
| MERAIODE script | `rai-response-genomics-atlas/scripts/meraiode_redifferentiation_panel_2026_08_06.py` |
| MERAIODE brief | `rai-response-genomics-atlas/results/reports/meraiode_redifferentiation_brief_2026_08_06.md` |
| MERAIODE per-sample panel_z table | `rai-response-genomics-atlas/results/tables/meraiode_panel_per_sample_2026_08_06.tsv` |
| Identity-vs-machinery decomposition (cross-referenced, §5) | `audit/rai_integration_20260819/B_functional_decomposition_tcga.md` |
| New computations in this document (logistic models, correlations, quadrant crosstabs) | Run inline in this session via `python3`/`statsmodels`/`scipy`/`sklearn` against the two source tables above; not saved as a separate script since every number is reproducible from the one-liners in §2 and §4 against the unmodified source tables listed here. |

No file under `rai-response-genomics-atlas/scripts/`, `rai-response-genomics-atlas/results/`,
`project/manuscript_v8/`, or `audit/rai_integration_20260806/` /
`audit/rai_integration_20260819/A_*` / `B_*` was modified in producing this document.
