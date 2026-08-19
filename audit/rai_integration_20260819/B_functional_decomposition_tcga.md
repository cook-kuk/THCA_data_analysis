# B — Functional decomposition of the 8-gene panel: which sub-module drives the post-RAI structural-disease association?

Generated 2026-08-19. New analysis — extends, does not replace or restate as-is, the audit
completed 2026-08-06 in `audit/rai_integration_20260806/` (`00_EXECUTIVE_DECISION.md` §3–4,
`04_STATISTICAL_REANALYSIS.md`, `02_EVIDENCE_LEDGER.tsv`).

**Verification status of every number below: VERIFIED — reproduced by a script run in this
session, not transcribed from narrative.** Script, data, and output paths are given in
§7. No existing manuscript, audit, or decision file was edited.

---

## 0. What this does and does not change

The prior audit established, in TCGA-THCA, within radioiodine-treated patients, using
Firth-penalised logistic regression on the aggregate 8-gene panel score `RAI_8`:

| endpoint | n | events | OR (score-only) | 95% CI | P |
|---|---|---|---|---|---|
| new tumour event after initial treatment | 145 | 9 | 0.221 | 0.062–0.784 | 0.019 |
| persistent disease within 3 months of surgery | 72 | 13 | 0.224 | 0.061–0.826 | 0.025 |

This document asks a question the prior audit did not ask: **does this association come from
the panel as a whole, or is it carried by one of its three functionally distinct components?**
The panel score is decomposed into three sub-modules, defined on biological grounds set by the
task, not by data-driven gene selection:

| Sub-module | Genes | Biology |
|---|---|---|
| `LINEAGE_4` | PAX8, NKX2-1, FOXE1, TSHR | thyrocyte-identity master transcription factors + TSH receptor |
| `UPTAKE_1` | SLC5A5 | iodide trapping (NIS) — the single gene radioiodine uptake most directly depends on |
| `HORMONE_3` | TG, TPO, DIO1 | hormone synthesis / organification machinery |

Same TCGA-THCA RAI-treated cohort, same two endpoints, same covariate sources, same Firth
tiers, same leave-one-event-out and patient-bootstrap stability checks as the 2026-08-06 audit.
Nothing here re-derives or restates the aggregate result — the aggregate numbers are reproduced
only as an exactness check (§2) that the sub-module scores are built on the identical cohort-z
basis as the original `RAI_8`.

---

## 1. Method — exact reuse of the original pipeline

**Cohort construction** reproduces `tcga_rai_best_response_2026_08_06.py::build()` line for
line (radioiodine-course aggregation from the GDC BCR Biotab radiation file, purity/leukocyte
merge, clinical merge), substituting the panel-score columns for the three sub-module scores.
Source files, unchanged from the original:

- Radiation courses: `rai-response-genomics-atlas/data/raw/TCGA_THCA_biotab/nationwidechildrens.org_clinical_radiation_thca.txt`
- Purity/leukocyte fraction: `/data/thca/repo_results/p2_braf_nature_sprint_2026_05_09/r6_purity_audit/r6_purity_per_sample.tsv`
- Clinical (stage, TUMOR_STATUS, CLINICAL_STATUS_WITHIN_3_MTHS_SURGERY, NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT): `/data/thca/repo_results/p2_braf_nature_sprint_2026_05_09/h24_survival_sensitivity/cache/cbio_thca_tcga_patient_clinical.tsv`

**Sub-module scoring** reuses the exact formula in
`project/results/ncomm_push_2026_05_08/_common.py::score_panel` — per-gene **cohort z-score**
`(x − cohort_mean) / cohort_std`, computed over the same 513-sample TCGA-PTC GDC study cohort
(`panel_expression_thpa_tcga_gdc.tsv`) that `RAI_8` itself was z-scored against — then averaged
within each sub-module (not re-centered on the 145/72-patient RAI-treated subset, which would
change the scoring basis and break comparability with the published `RAI_8`).

**Statistical model** reuses the Firth IRLS implementation (`firth_logit`) and the Cohen's-d /
bootstrap-CI / Mann–Whitney AUC helpers (`cohens_d`, `boot_ci_d`, `auc_mw`) imported directly
from `rai-response-genomics-atlas/scripts/audit/03_rare_event_structural_outcomes.py` and
`rai-response-genomics-atlas/scripts/tcga_rai_best_response_2026_08_06.py` — not reimplemented,
so results are numerically identical to the original where the original applies. Same four
Firth tiers per score × endpoint (score-only; +stage; +purity; +stage+purity, exploratory),
same leave-one-event-out refits, same 2,000-draw patient bootstrap, same seed (`20260806`).

## 2. Exactness check — sub-module scores reproduce the published `RAI_8`

Before testing anything, the recomputed 8-gene aggregate was checked against the stored
`RAI_8` used throughout the manuscript:

```
max |RAI_8 − recomputed 8-gene mean-z| across 513 samples = 6.7e-16   (floating-point noise)
max |recomputed RAI_8 − (4·LINEAGE_4 + 1·UPTAKE_1 + 3·HORMONE_3)/8|   = 4.4e-16
```

Both are exact to machine precision. The sub-module scores are on the identical scale and
cohort-z basis as `RAI_8`, and `RAI_8 = (4·LINEAGE_4 + 1·UPTAKE_1 + 3·HORMONE_3)/8` exactly —
the weighted mean of the three sub-modules **is** the published score, so any discrepancy in
which module associates with outcome is a real decomposition, not a scoring artifact.

Cohort sizes after the identical merge also reproduce exactly: n=145/events=9 (new tumour
event) and n=72/events=13 (persistent disease), matching `04_STATISTICAL_REANALYSIS.md` row
for row.

## 3. Primary result — Firth, score-only, one covariate family (the sub-module) at a time

| Endpoint | Score | n | events | OR | 95% CI | P | q (BH, family of 6) | Bonferroni (α=0.0083) |
|---|---|---|---|---|---|---|---|---|
| new tumour event | **RAI_8 (reproduced)** | 145 | 9 | 0.221 | 0.062–0.784 | 0.0194 | — | — |
| new tumour event | LINEAGE_4 | 145 | 9 | 0.440 | 0.182–1.065 | 0.0688 | 0.103 | not sig. |
| new tumour event | UPTAKE_1 | 145 | 9 | 0.259 | 0.026–2.591 | 0.2502 | 0.250 | not sig. |
| new tumour event | HORMONE_3 | 145 | 9 | **0.231** | 0.066–0.801 | **0.0209** | 0.125 | not sig. |
| persistent disease | **RAI_8 (reproduced)** | 72 | 13 | 0.224 | 0.061–0.826 | 0.0245 | — | — |
| persistent disease | LINEAGE_4 | 72 | 13 | 0.360 | 0.123–1.056 | 0.0628 | 0.126 | not sig. |
| persistent disease | UPTAKE_1 | 72 | 13 | 0.309 | 0.060–1.591 | 0.1602 | 0.192 | not sig. |
| persistent disease | HORMONE_3 | 72 | 13 | **0.279** | 0.092–0.848 | **0.0244** | 0.073 | not sig. |

Source: `rai-response-genomics-atlas/results/tables/audit12_module_firth_main_2026_08_19.tsv`
(Tier 2 rows). The `RAI_8 (reproduced)` rows are the exactness check from §2, included for
comparison — they are not a new claim.

**Reading this table honestly:**

- Point estimates for **all three** sub-modules sit on the same side of 1 (odds of the adverse
  outcome lower with a higher differentiation score) in both endpoints — directionally
  consistent with the aggregate, as expected since each is a component of it.
- **HORMONE_3** (TG, TPO, DIO1) is the only sub-module whose confidence interval excludes 1 in
  both endpoints at the nominal (uncorrected) α=0.05, and its point estimate (OR 0.23 / 0.28)
  is closest to the aggregate's (OR 0.22 / 0.22).
- **LINEAGE_4** (PAX8, NKX2-1, FOXE1, TSHR) trends the same direction but does not reach
  nominal significance in either endpoint (P=0.069, P=0.063); its point estimate is weaker
  (OR 0.44 / 0.36, roughly half-way to the null compared with the aggregate).
- **UPTAKE_1** (SLC5A5 alone) is the weakest and least precise: OR 0.26–0.31, P=0.16–0.25, and
  a 95% CI spanning more than two orders of magnitude for the new-tumour-event endpoint
  (0.026–2.59) — the data cannot distinguish this single-gene score from the null.
- **None of the six sub-module tests survives Bonferroni correction for the 6-test family**
  (α = 0.05/6 = 0.0083); the smallest raw P (HORMONE_3, new tumour event, P=0.0209) is roughly
  2.5× the corrected threshold. Benjamini–Hochberg q-values are 0.07–0.25 across the family.

## 4. Sensitivity — same purity/stage adjustment as the original

The original audit's key robustness check was that the aggregate association survived
adjustment for tumour purity (leukocyte-fraction proxy). The same check was run for each
sub-module (Tier 3, `rai-response-genomics-atlas/results/tables/audit12_module_firth_main_2026_08_19.tsv`):

| Endpoint | Score | +stage OR (P) | +purity OR (P) | +stage+purity OR (P) |
|---|---|---|---|---|
| new tumour event | LINEAGE_4 | 0.449 (0.073) | 0.406 (0.059) | 0.417 (0.064) |
| new tumour event | UPTAKE_1 | 0.283 (0.259) | 0.230 (0.228) | 0.251 (0.238) |
| new tumour event | HORMONE_3 | 0.251 (0.025) | 0.232 (0.019) | 0.253 (0.024) |
| persistent disease | LINEAGE_4 | 0.376 (0.074) | 0.333 (0.059) | 0.350 (0.071) |
| persistent disease | UPTAKE_1 | 0.312 (0.162) | 0.255 (0.148) | 0.264 (0.152) |
| persistent disease | HORMONE_3 | 0.292 (0.027) | 0.289 (0.026) | 0.302 (0.029) |

Purity adjustment does not attenuate any sub-module's estimate (mirroring the aggregate
finding); HORMONE_3 remains the only sub-module nominally significant after adjustment in
both endpoints. This is a sensitivity check on an already-non-primary comparison — it is not
treated as an independent confirmatory result.

## 5. Stability — leave-one-event-out and patient bootstrap (score-only model)

`rai-response-genomics-atlas/results/tables/audit12_module_stability_2026_08_19.tsv`

| Endpoint | Score | LOO OR range | direction stable on LOO | Boot median OR | Boot 95% CI | % boot OR<1 |
|---|---|---|---|---|---|---|
| new tumour event | RAI_8 | 0.121–0.282 | yes | 0.208 | 0.033–0.774 | 98.8% |
| new tumour event | LINEAGE_4 | 0.372–0.522 | yes | 0.434 | 0.132–0.894 | 98.6% |
| new tumour event | UPTAKE_1 | **0.000**–0.359 | yes | 0.255 | **0.000**–0.954 | 98.3% |
| new tumour event | HORMONE_3 | 0.040–0.291 | yes | 0.216 | 0.005–1.140 | 96.6% |
| persistent disease | RAI_8 | 0.156–0.274 | yes | 0.214 | 0.032–0.611 | 100.0% |
| persistent disease | LINEAGE_4 | 0.276–0.440 | yes | 0.348 | 0.078–0.926 | 98.5% |
| persistent disease | UPTAKE_1 | 0.189–0.357 | yes | 0.306 | 0.027–0.759 | 99.9% |
| persistent disease | HORMONE_3 | 0.219–0.326 | yes | 0.274 | 0.060–0.602 | 100.0% |

`loo_all_below_1` holds for every sub-module and endpoint — the direction of every point
estimate is not an artifact of any single patient. But for **UPTAKE_1**, one leave-one-out
refit collapses the odds ratio to 0.000 (near-complete separation on 9 events with a single
noisy predictor) and the bootstrap 2.5th percentile also touches 0.000 — a sign of numerical
fragility in the single-gene score, not merely wide uncertainty. HORMONE_3's bootstrap upper
CI (1.140, new-tumour-event endpoint) crosses 1 even though the point estimate and Firth CI
did not, another sign that stability is driven by very few events rather than a settled
estimate.

## 6. Why HORMONE_3 tracks the aggregate more closely than gene-count alone would predict

`LINEAGE_4` averages 4 genes, `HORMONE_3` averages 3 — if sub-module strength were only a
function of how many genes are averaged (more genes → less measurement noise → tighter
estimate), LINEAGE_4 should be at least as strong as HORMONE_3. It is not: Spearman
correlation with the full 8-gene score is **0.905 for HORMONE_3 vs 0.878 for LINEAGE_4**
despite HORMONE_3 having one fewer component gene
(`rai-response-genomics-atlas/results/tables/audit12_module_scores_2026_08_19.tsv`, n=505).
`UPTAKE_1` correlates only 0.512 with the aggregate. This argues (does not prove) that the
difference in association strength between sub-modules reflects something beyond a trivial
averaging effect — but with three data points (module-vs-aggregate correlations) this is
descriptive context, not a statistical test, and is not offered as a claim on its own.

`SLC5A5` (the sole `UPTAKE_1` gene) is markedly right-skewed in this cohort (median cohort log2
expression 0.12, IQR 0.04–0.73, max 9.28 across 513 samples) — most tumours sit near a floor of
near-complete NIS silencing with a long tail of preserved expression. A single skewed gene, not
averaged with anything, is the most plausible mechanical reason `UPTAKE_1` produces the widest
intervals and the LOO/bootstrap instability in §5, independent of any biological claim about
NIS itself.

## 7. Multiplicity — explicit statement

The 2026-08-06 audit tested **one** aggregate score against these two endpoints (plus a third,
null, endpoint — `structural disease at last follow-up` — that this document does not revisit).
This document tests **three** sub-module scores against the same two endpoints: a 6-test
primary family where the original had 1. That inflation of the number of comparisons, run on
outcomes with only 9 and 13 events, is exactly the situation in which nominal P-values around
0.02–0.07 are the least trustworthy. Two multiplicity corrections were applied and are reported
alongside the raw values in §3, not substituted for them:　

- **Bonferroni** (conservative, assumes independence): α = 0.05/6 = 0.0083. **Zero of six
  sub-module tests reach this threshold** — including HORMONE_3's smallest P (0.0209), which
  is 2.5× the corrected α.
- **Benjamini–Hochberg** (less conservative, controls false-discovery rate): q-values 0.07–0.25
  across the family. None reach the conventional q<0.05 threshold either, though HORMONE_3's
  persistent-disease q (0.073) is the closest.

The sub-modules are also not independent tests in the Bonferroni sense (Spearman correlations
0.30–0.65 between module pairs, §6), so a naive Bonferroni is itself somewhat conservative here
— but going the other direction (treating 6 correlated tests as if they were 1) is not
defensible either. Reporting both, uncorrected alongside both corrections, is the honest
position given the ambiguity.

## 8. Verdict

**This decomposition is descriptive/exploratory, not a refinement that upgrades the aggregate
claim, and not evidence the aggregate association is noise.**

1. **Not a refinement in the sense of "we found the causal gene(s)."** No sub-module
   individually survives correction for the multiplicity this analysis itself introduces.
   HORMONE_3 is nominally associated with both endpoints at the uncorrected α=0.05 and has the
   point estimate closest to the aggregate's, but that is a suggestive pattern across two
   correlated endpoints on 9 and 13 events, not a confirmed independent finding.
2. **Not a restatement either** — the aggregate `RAI_8` result already existed and is not
   changed by this analysis; this document adds new information (relative behavior of the
   three sub-modules) that was not previously known and could have come out differently
   (e.g., UPTAKE_1 alone driving everything, or all three splitting the signal evenly — neither
   happened).
3. **The honest reading:** the association is not obviously concentrated in the single gene
   (`SLC5A5`) most mechanistically tied to radioiodine uptake — if anything that sub-module is
   the noisiest and least stable of the three. It is somewhat more consistent with the hormone
   synthesis/organification module (TG/TPO/DIO1) than with the lineage/master-TF module, but
   this is a gradient across three correlated, underpowered comparisons, not a clean
   dichotomy, and does not clear a defensible multiplicity-corrected significance bar.
4. **Recommended manuscript-facing language, if this is used at all:** "A post-hoc
   decomposition of the eight-gene score into lineage, uptake, and hormone-synthesis
   sub-modules found the same direction of association for all three, with the hormone-
   synthesis sub-module (TG, TPO, DIO1) nominally significant in both endpoints (P=0.02–0.02)
   and the single-gene uptake sub-module (SLC5A5) the weakest and least stable; none of the
   three sub-module tests survived Bonferroni correction for the six comparisons this analysis
   performed, so this decomposition is reported as exploratory and does not establish which
   component gene(s) are causally responsible." Anything stronger than that is not supported by
   these numbers.
5. **Do not use this analysis to name a "driver gene."** Especially not `SLC5A5` — the data
   here argue mildly against, not for, an NIS-centric mechanistic story, and even that reading
   rests on 9 and 13 events split three further ways.

## 9. Limitations specific to this decomposition (beyond those already logged in the 2026-08-06 audit)

- Splitting one 8-gene composite into three shorter composites (and one singleton) increases
  the standard error of each per-patient score relative to the aggregate; some of the
  attenuation in LINEAGE_4 and UPTAKE_1's estimates could reflect added measurement noise from
  averaging fewer genes rather than true absence of association in those genes' underlying
  biology. This analysis cannot separate "genuinely weaker association" from "noisier
  estimate of a similar true association" with 9–13 events.
- The module boundaries (lineage / uptake / hormone-synthesis) were specified by the task
  as a biological hypothesis before this script ran; they were not derived from the data
  (e.g., not by clustering genes by their individual associations with outcome), which is the
  correct order of operations for a pre-specified decomposition, but the choice of exactly
  these three groupings among several biologically plausible alternatives (e.g., TSHR could be
  argued into a "regulatory" module rather than "lineage") is itself a modeling decision worth
  disclosing.
- Endpoints, cohort, and covariate handling reuse every limitation already logged for the
  aggregate result in `audit/rai_integration_20260806/00_EXECUTIVE_DECISION.md` §3 (neither
  TCGA proxy is an adjudicated ATA response category; `persistent disease within 3 months of
  surgery` overlaps the interval to radioiodine and is not a pure post-treatment endpoint) —
  those caveats apply unchanged to every sub-module result in this document.

---

## 10. Provenance

| Item | Path |
|---|---|
| New script (this analysis) | `rai-response-genomics-atlas/scripts/audit/12_functional_decomposition_2026_08_19.py` |
| Module scores (per-sample, n=505) | `rai-response-genomics-atlas/results/tables/audit12_module_scores_2026_08_19.tsv` |
| Firth main results (all tiers, both endpoints, 4 scores) | `rai-response-genomics-atlas/results/tables/audit12_module_firth_main_2026_08_19.tsv` |
| LOO / bootstrap stability | `rai-response-genomics-atlas/results/tables/audit12_module_stability_2026_08_19.tsv` |
| Reused (unmodified) — cohort builder, Firth stats helpers | `rai-response-genomics-atlas/scripts/tcga_rai_best_response_2026_08_06.py` |
| Reused (unmodified) — Firth IRLS implementation | `rai-response-genomics-atlas/scripts/audit/03_rare_event_structural_outcomes.py` |
| Reused (unmodified) — original structural-disease pipeline (not re-run, cited for comparability only) | `rai-response-genomics-atlas/scripts/tcga_rai_structural_disease_2026_08_06.py` |
| Raw 8-gene expression + published `RAI_8`, cohort-z base (n=513) | `project/results/ncomm_push_2026_05_08/cbioportal_sweep/panel_expression_thpa_tcga_gdc.tsv` |
| Radiation courses (GDC BCR Biotab) | `rai-response-genomics-atlas/data/raw/TCGA_THCA_biotab/nationwidechildrens.org_clinical_radiation_thca.txt` |
| Purity / leukocyte fraction | `/data/thca/repo_results/p2_braf_nature_sprint_2026_05_09/r6_purity_audit/r6_purity_per_sample.tsv` |
| Clinical (stage, tumour status, new-tumour-event) | `/data/thca/repo_results/p2_braf_nature_sprint_2026_05_09/h24_survival_sensitivity/cache/cbio_thca_tcga_patient_clinical.tsv` |
| Prior audit this extends | `audit/rai_integration_20260806/00_EXECUTIVE_DECISION.md`, `04_STATISTICAL_REANALYSIS.md`, `02_EVIDENCE_LEDGER.tsv` |

No file under `audit/rai_integration_20260806/`, `project/manuscript_v8/`, or any other
existing manuscript/audit/decision file was modified in producing this document.
