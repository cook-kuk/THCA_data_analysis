# TCGA-THCA — 8-gene panel vs structural disease after radioiodine

**Date** 2026-08-06 · **Script** `scripts/tcga_rai_structural_disease_2026_08_06.py`
**Cohort** 235 RAI-treated TCGA-THCA patients (mCi-dosed courses) with an 8-gene panel score.

## Headline

Within patients who actually received radioiodine, a **lower** 8-gene differentiation score
is associated with post-treatment structural disease events, and the association survives
adjustment for the thyrocyte-content confounder that destroyed the equivalent signal in
GSE151179.

| Endpoint | n (event / no event) | Cohen's d | 95% CI | P | q (BH) |
|---|---|---|---|---|---|
| **New tumour event after initial treatment** | 9 / 136 | **−0.84** | −1.37 to −0.25 | **0.0049** | 0.015 |
| **Persistent disease within 3 mo of surgery** | 13 / 59 | **−0.78** | −1.29 to −0.31 | **0.012** | 0.018 |
| Structural disease at last follow-up | 34 / 182 | +0.06 | −0.33 to +0.48 | 0.88 | 0.88 |

Adjusted logistic regression, odds of the event per 1-unit increase in panel z:

| Endpoint | Model | OR | 95% CI | P |
|---|---|---|---|---|
| New tumour event | panel z | 0.203 | 0.054–0.761 | 0.018 |
| New tumour event | + purity | **0.177** | 0.044–0.717 | **0.015** |
| New tumour event | + purity + stage | **0.187** | 0.046–0.764 | **0.020** |
| Persistent disease | panel z | 0.196 | 0.050–0.767 | 0.019 |
| Persistent disease | + purity | **0.186** | 0.046–0.754 | **0.019** |
| Persistent disease | + purity + stage | **0.191** | 0.046–0.790 | **0.022** |

## Why the purity result matters

The eight panel genes are thyrocyte-specific, so panel score tracks tumour cellularity.
In this cohort panel z is strongly correlated with the leukocyte-fraction proxy
(Spearman ρ = −0.383, P = 1.3 × 10⁻⁹) — but leukocyte fraction is **not** associated with
either outcome (persistent disease P = 0.73; new tumour event P = 0.82). Purity therefore
predicts the exposure but not the outcome, which by definition is not a confounder here,
and the adjusted odds ratios move *away from* the null rather than toward it.

This is the mirror image of GSE151179, where the same covariate collapsed a d = +0.37
association to β = +0.035 (`gse151179_uptake_at_met_site_brief_2026_08_06.md`). Running
both analyses with the identical covariate set is what makes the contrast interpretable.

## Three things this does NOT show

1. **It is not RAI-specific — predictive versus prognostic is unresolved.** In the 265
   patients with no recorded mCi course there are only 2 new-tumour events and 1 persistent
   -disease event, so the treatment × panel interaction cannot be estimated
   (new tumour event: β = −2.23, P = 0.077; persistent disease: not estimable). The
   association is established *within* RAI-treated patients only. A prognostic signature
   would produce exactly this pattern, and the data cannot currently separate the two.
   Any claim that the panel identifies who benefits from RAI requires this interaction and
   must not be made until it is powered.
2. **The initial response to RAI is unrelated to the panel.** The RECIST-flavoured
   `treatment_best_response` endpoint on the same patients is a flat null
   (d = −0.03, P = 0.73; `tcga_rai_best_response_brief_2026_08_06.md`). The dissociation —
   null on first-treatment response, positive on subsequent structural events — is
   biologically coherent (86% of this low-risk cohort achieve complete response, leaving
   no room to discriminate) but it is a dissociation, not a single consistent effect.
3. **Event counts are small and one endpoint straddles the treatment.** Nine and thirteen
   events give wide intervals. `CLINICAL_STATUS_WITHIN_3_MTHS_SURGERY` is recorded within
   three months of surgery, a window that overlaps the usual 4–12 week interval to RAI, so
   it partly reflects post-surgical residual disease rather than response to radioiodine.
   `NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT` is the cleaner post-treatment endpoint.
   Neither is an adjudicated ATA response-to-therapy category — no public dataset carries
   those.

## Defensible wording

> Among TCGA-THCA patients who received radioiodine (n = 235, median cumulative activity
> 102.8 mCi), a lower 8-gene thyroid-differentiation score was associated with a new tumour
> event after initial treatment (OR 0.19 per unit, 95% CI 0.05–0.76, P = 0.02, adjusted for
> tumour cellularity and stage), while showing no association with the documented best
> response to the initial radioiodine course (Cohen's d = −0.03, P = 0.73). Because only
> two events occurred among patients with no recorded radioiodine course, the
> treatment-by-signature interaction could not be estimated, and the association is
> therefore reported as prognostic within a radioiodine-treated population rather than as
> evidence of differential benefit from radioiodine.

## Files

- `results/tables/tcga_rai_structural_main_2026_08_06.tsv`
- `results/tables/tcga_rai_structural_adjusted_2026_08_06.tsv`
- `results/figures/figure_tcga_rai_structural_2026_08_06.{png,pdf}`

Seed 20260806 · bootstrap 5,000 draws · statsmodels logit · covariates dropped to keep
events-per-variable ≥ 5, with the purity-inclusive models reported separately as the
decisive sensitivity analysis.
