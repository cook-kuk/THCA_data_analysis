# TCGA-THCA — 8-gene panel vs documented radioiodine treatment response

**Date** 2026-08-06 · **Script** `scripts/tcga_rai_best_response_2026_08_06.py`

## The dataset we did not know we had

Every prior statement in this project — including manuscript Limitation 6 — asserts that
no cohort links RAI outcome to per-sample DM1 calls. That is **false**. The GDC BCR Biotab
clinical supplement for TCGA-THCA contains a radiation file with per-course fields:

| Field | Content |
|---|---|
| `radiation_adjuvant_units` | **mCi 249** (millicuries = radioiodine), Gy 14, cGy 3 (external beam) |
| `radiation_total_dose` | cumulative activity; median **102.8 mCi** (IQR 98.0–150.0) |
| `treatment_best_response` | **CR 167 · PR 21 · SD 6 · Radiographic PD 7** |
| `radiation_therapy_type` | Systemic 263, RADIOISOTOPE 14, External 11 |

293 course records across 274 patients. **None of these fields are exposed by any
cBioPortal THCA study**, which is why they were missed for the life of the project.
Retrieval procedure and file path: memory `tcga-thca-rai-response-biotab-source`.

Joined to the R17 8-gene panel table: **235 RAI-treated patients with a panel score, 167
with an evaluable best response** (CR 143 / PR 15 / SD 4 / PD 5).

## Result: a clean null, not an underpowered one

Endpoint = failure to achieve complete response to RAI (non-CR vs CR). Overall survival is
not used: this cohort has **5 deaths among the 235 RAI-treated patients** (16/500 in the
full THCA set), and the ATA 2015 evidence base establishes an RAI survival benefit only for
T4 gross extrathyroidal extension and M1 disease, on observational data alone.

| Analysis | n (non-CR / CR) | Effect | 95% CI | P |
|---|---|---|---|---|
| **PRIMARY** panel z, non-CR vs CR | 24 / 143 | d = **−0.03** | −0.50 to +0.46 | 0.73 |
| RAI only, no external beam | 23 / 142 | d = −0.08 | −0.56 to +0.43 | 0.54 |
| DM1 vs DM2 × non-CR (Fisher) | 24 / 143 | OR = 1.31 | — | 0.64 |
| Ordinal gradient CR→PR→SD→PD | 167 | ρ = **−0.036** | — | 0.65 |
| Panel z vs cumulative mCi (Spearman) | 235 | ρ = −0.079 | — | 0.23 |
| PFI within RAI-treated (Cox) | 235, 36 events | HR = 0.76 | 0.44–1.33 | 0.34 |

Adjusted logistic model for non-CR (n = 167), covariates chosen as the standard confounders
of a differentiation score and of RAI prescribing:

| Term | OR | 95% CI | P |
|---|---|---|---|
| 8-gene panel z | 0.85 | 0.34–2.12 | 0.73 |
| Stage III/IV | 2.78 | 0.81–9.54 | 0.10 |
| BRAF V600E | 0.87 | 0.30–2.47 | 0.79 |
| Age | 0.97 | 0.94–1.02 | 0.23 |
| Leukocyte fraction (inverse purity) | 0.14 | 0.002–13.2 | 0.40 |
| Cumulative RAI activity (mCi) | 1.004 | 0.998–1.011 | 0.18 |

**Power:** with 24 vs 143 the study detects d ≥ 0.63 at 80% power. The observed effect is
d = −0.03 — indistinguishable from zero, and the confidence interval excludes any effect
larger than ±0.5. This is a genuine null over the range that would matter clinically, not
merely a failure to reach significance.

## Why the null is expected, and what it means

TCGA-THCA is overwhelmingly early-stage papillary carcinoma: **143/167 = 86% achieved
complete response** to a median 103 mCi. This is precisely the population in which
ESTIMABL2 (*NEJM* 2022;386:923-32) and IoN (*Lancet* 2025;406:52-62) showed that omitting
RAI is non-inferior. A differentiation-silencing signature has almost no discriminative
room when 86% of patients respond completely, and the ~14% who do not are dominated by
stage rather than by lineage state (stage III/IV OR = 2.78 versus panel OR = 0.85).

The published positive results for molecule-versus-RAI-response all come from
**metastatic or advanced** DTC, where the refractory fraction is 30–50%:
Boucai 2023 (*Clin Cancer Res* 29:1620-30; 24 metastatic, RECIST, median time to
progression 76 vs 11 months), Mu 2024 (*JCEM* 109:1231-40; 214 distant-metastatic, four
uptake classes, time-to-refractoriness), Siraj 2022 (*Cancers* 14:1584; 158 PTC, 66
refractory / 92 avid), Laschinsky 2023 (*J Nucl Med* 64:1865-8; TERT-mutant time-to-
refractoriness 0.7 vs 19.8 months).

**Conclusion:** the correct target population is advanced/metastatic DTC, not a
population-based primary-tumour cohort. TCGA-THCA can now be reported as an explicit,
adequately powered negative control for the RAI-response question — which is a stronger
position than silence, because it pre-empts the reviewer who would otherwise find these
fields and ask why we ignored them.

## Consequence for the manuscript

- Limitation 6 must be **rewritten**: the statement "no cohort in this study contains RAI
  outcome data" is factually wrong and would be embarrassing if a reviewer checked GDC.
  Replace with the tested null and its power.
- The panel must not be described as predicting or being associated with RAI response.
  "RAI-responsiveness panel" as a *name* is defensible only as a statement about the
  genes' canonical biology, and should be flagged as such at first use.
- The forward-looking claim becomes specific and testable: the hypothesis survives only in
  advanced/metastatic disease, and the cohorts that could adjudicate it are named above.

## Files

- `results/tables/tcga_rai_response_main_2026_08_06.tsv`
- `results/tables/tcga_rai_response_adjusted_2026_08_06.tsv`
- `results/tables/tcga_rai_response_dose_2026_08_06.tsv`
- `results/tables/tcga_rai_response_pfi_2026_08_06.tsv`
- `results/tables/tcga_rai_response_power_2026_08_06.tsv`
- `results/figures/figure_tcga_rai_best_response_2026_08_06.{png,pdf}`
- Source data: `/data/rai_atlas/raw/TCGA_THCA_biotab/nationwidechildrens.org_clinical_radiation_thca.txt`

Seed 20260806 · bootstrap 5,000 draws · lifelines Cox · statsmodels logit.
