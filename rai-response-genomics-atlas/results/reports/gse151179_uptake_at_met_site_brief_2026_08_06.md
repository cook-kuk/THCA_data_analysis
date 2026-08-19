# GSE151179 — 8-gene panel vs RAI uptake at the metastatic site

**Date** 2026-08-06 · **Script** `scripts/gse151179_uptake_at_met_site_2026_08_06.py`
**Trigger** Kang Minsu (2026-08-06): the inferential step from public survival data to
radioiodine treatment effect is the weakest link and will be attacked as a logical leap.

## What was new here

All prior atlas work tested the **patient-level** GEO field `patient rai responce`
(Refractory 46 / Avid 6) — a 7.7:1 imbalance that cannot support a test. This analysis
instead used a field that had never been tested:

> `rai uptake at the metastatic site: Yes / No` — **Yes 25 / No 27** (39 tumour
> specimens: Yes 19 / No 20, across 32 patients)

This is the closest available direct measurement of *RAI effect* — whether radioiodine
actually reached and was retained by disease — as opposed to *RAI receipt*, which every
sample in this series shares (all 52 patients were treated with RAI after thyroidectomy).

## Result: directionally correct, statistically null, and confounded

| Analysis | Unit | n (Yes / No) | Effect | 95% CI | P |
|---|---|---|---|---|---|
| **PRIMARY** patient-level, tumours | patient | 16 / 16 | d = **+0.37** | −0.34 to +0.98 | 0.53 |
| sample-level, tumours (non-independent) | sample | 19 / 20 | d = +0.27 | −0.39 to +0.84 | 0.68 |
| primary tumours only | sample | 11 / 6 | d = −0.02 | −1.51 to +1.09 | 0.52 |
| pre-RAI specimens only | sample | 14 / 8 | d = +0.34 | −0.79 to +1.26 | 0.97 |
| high-purity tumours only | sample | 12 / 10 | d = +0.03 | −1.06 to +0.74 | 0.67 |
| **NEGATIVE CONTROL** non-neoplastic thyroid | sample | 6 / 7 | d = +0.65 | −0.33 to +2.29 | 0.18 |

Direction is as hypothesised throughout (lesions that took up RAI have *higher*
differentiation-panel expression), but no comparison reaches significance, and two
findings make the nominal effect uninterpretable:

1. **The negative control is larger than the signal.** Non-neoplastic thyroid from the
   same patients shows d = +0.65 versus d = +0.37 in tumour. An association that is
   stronger in normal tissue than in tumour is a patient-level artefact, not tumour biology.

2. **Tumour purity explains the whole effect.** In the adjusted model
   `panel_z ~ uptake + purity_class + driver_class` (39 tumour specimens):

   | Term | β | 95% CI | P |
   |---|---|---|---|
   | CIBERSORT purity: low vs high | **−0.880** | −1.280 to −0.480 | **1.9 × 10⁻⁴** |
   | RAI uptake: Yes vs No | **+0.035** | −0.440 to +0.511 | 0.88 |
   | driver class (fusion / pTERT / WT / other vs BRAF V600E) | all n.s. | — | 0.37–0.97 |

   The eight panel genes are thyrocyte-specific, so the score is dominated by thyrocyte
   content. Once that is held constant, the RAI-uptake association is essentially zero.
   This is the same confounder that reversed the direction of the Pu 2021 single-cell
   analysis (`results/reports/pu2021_rai_refractory_brief.md`).

Per-gene, the largest patient-level effects are TPO (d = +0.58, q = 0.51) and TG
(d = +0.58, q = 0.69); *SLC5A5*/NIS — the gene with the most direct mechanistic claim on
iodide transport — is flat (d = +0.06, P = 0.95). Leave-one-gene-out is stable
(d = +0.32 to +0.42), so no single gene carries or destroys the result.

## Power: what this cohort could never have shown

| Quantity | Value |
|---|---|
| Observed patient-level d | 0.373 |
| Power to detect that d at 16 vs 16 | **0.18** |
| d detectable at 80% power | **1.02** |
| Patients per group needed for the observed d at 80% power | **113** |
| **Total patients needed** | **≈ 226** |

GSE151179 could only ever have detected a very large effect (d ≥ 1.0). It is not evidence
against the hypothesis; it is evidence that this cohort cannot adjudicate it.

## Consequence for the manuscript

- The claim "the panel is associated with RAI uptake" is **not supported** by the only
  public lesion-level RAI-uptake label available. It must not appear.
- Any future RAI-response analysis **must report tumour purity / thyrocyte content as a
  covariate**, or a reviewer will correctly attribute the finding to composition.
- The design target is now quantitative: **≈ 226 patients with a lesion-level RAI uptake
  or response label and a purity estimate.** Mu et al. 2024 (NGDC HRA004166, n = 214,
  four-class uptake pattern) is the only identified cohort of approximately that size.

## Files

- `results/tables/gse151179_uptake_at_met_site_2026_08_06.tsv`
- `results/tables/gse151179_uptake_per_gene_2026_08_06.tsv`
- `results/tables/gse151179_uptake_logo_2026_08_06.tsv`
- `results/tables/gse151179_uptake_power_2026_08_06.tsv`
- `results/figures/figure_gse151179_uptake_at_met_site_2026_08_06.{png,pdf}`

Seed 20260806 · bootstrap 5,000 draws · pandas/scipy/statsmodels as in `requirements.txt`.
