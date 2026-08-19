# 00 — Executive decision

**Date** 2026-08-06 · **Branch** `paper9-perturbation-extension-20260506` · **No commit made.**
Working tree was already dirty (281 modified, 1,405 untracked) when this audit began, so no
branch was created and nothing was committed or overwritten; all audit output is new files.

Every verdict below is answered in the required format. Numbers are the audited values, not
the first-pass values, and where the two differ the difference is stated.

---

## 1. Is the eight-gene panel a predictor of initial radioiodine response?

**Verdict:** No, and this is now robust to the analytic choice that was previously implicit.

**Evidence:** TCGA-THCA, patients with an evaluable per-course response and a panel score.
Five index-course rules were tested because the first pass silently used "best response
across courses":

| Index rule | n (non-CR / CR) | Cohen's d | 95% CI | P |
|---|---|---|---|---|
| first course (earliest start day) | 26 / 141 | −0.000 | −0.445 to +0.457 | 0.91 |
| last course | 24 / 143 | −0.030 | −0.505 to +0.459 | 0.73 |
| highest-dose course | 24 / 143 | −0.030 | −0.505 to +0.459 | 0.73 |
| best response across courses (first pass) | 24 / 143 | −0.030 | −0.505 to +0.459 | 0.73 |
| worst response across courses | 26 / 141 | −0.000 | −0.445 to +0.457 | 0.91 |

**Critical caveat:** the first pass wrote that the cohort "detects d ≥ 0.63 at 80% power, so
clinically meaningful effects are excluded". That conflates a minimum-detectable-effect
calculation with an equivalence test, and no equivalence margin was pre-specified. **That
sentence is withdrawn.** The defensible statement is that the interval excludes effects
larger than |d| ≈ 0.46 and that no association was detected.

**Recommended manuscript wording:** "The thyroid-lineage differentiation score showed no
detectable association with the recorded response to the initial radioiodine course in this
predominantly complete-response cohort (Cohen's d = 0.00, 95% CI −0.45 to +0.46, P = 0.91;
141 of 167 evaluable patients achieved complete response)."

**Confidence:** high.

---

## 2. Is it a radioiodine-specific predictive biomarker?

**Verdict:** No, and this cannot be determined from any currently available data.

**Evidence:** among 265 patients with no recorded mCi course there are 2 new-tumour events
and 1 persistent-disease event. The treatment × score interaction is not estimable
(β = −2.23, P = 0.077 — a number that should not be quoted as a trend).

**Critical caveat:** a purely prognostic marker produces exactly the observed pattern. No
claim of differential benefit from radioiodine is permissible.

**Recommended manuscript wording:** "Because only two events occurred among patients without
a recorded radioiodine course, the treatment-by-score interaction could not be estimated;
the association is therefore reported as prognostic within a radioiodine-treated population
rather than as evidence of differential treatment benefit."

**Confidence:** high.

---

## 3. Is it associated with post-treatment structural outcomes?

**Verdict:** Yes, within radioiodine-treated patients, with the magnitude revised.

**Evidence:** Firth-penalised logistic regression, one covariate added at a time:

| Endpoint | Model | EPV | OR | 95% CI | P |
|---|---|---|---|---|---|
| New tumour event (9 events) | score only | 9.0 | **0.221** | 0.062–0.784 | 0.019 |
| | + stage | 4.5 | 0.239 | 0.069–0.829 | 0.024 |
| | + purity | 4.5 | 0.202 | 0.054–0.750 | 0.017 |
| | + stage + purity *(exploratory)* | 3.0 | 0.219 | 0.060–0.798 | 0.021 |
| Persistent disease (13 events) | score only | 13.0 | **0.224** | 0.061–0.826 | 0.025 |
| | + stage + purity *(exploratory)* | 4.3 | 0.237 | 0.064–0.874 | 0.031 |

**Difference from the first pass:** the first pass used ordinary maximum-likelihood logistic
with up to six covariates on 9 and 13 events (≈1.5–2 events per parameter) and reported
OR ≈ 0.19. Small-sample ML logistic is biased away from the null. Firth penalisation moves
the estimates to 0.221 and 0.224 — the direction and significance survive, the point
estimates were mildly overstated.

**Critical caveat:** `persistent disease within 3 months of surgery` overlaps the usual 4–12
week interval to radioiodine and is not a pure post-radioiodine endpoint.
`tumour present at final follow-up` (34 events) is null and mixes pre- and post-treatment
disease; it should be reported as an ambiguous endpoint, not buried.

**Confidence:** moderate.

---

## 4. Does the association survive rare-event correction?

**Verdict:** Yes.

**Evidence:** leave-one-event-out refits give OR 0.121–0.282 (new tumour event) and
0.156–0.274 (persistent disease); every refit stays below 1. Patient-level bootstrap
(2,000 refits) gives median OR 0.208 and 0.218, with 98.8% and 99.9% of refits below 1.

**Critical caveat:** stability is not power. Nine events remain nine events, and the
intervals are wide.

**Confidence:** moderate.

---

## 5. Does it survive purity adjustment?

**Verdict:** Yes for TCGA; no for GSE151179.

**Evidence:** TCGA — panel z correlates with leukocyte fraction (Spearman ρ = −0.383,
P = 1.3 × 10⁻⁹), and the Firth odds ratio is essentially unchanged with purity in the model
(0.221 → 0.202). GSE151179 — the same adjustment collapses the association
(β = +0.035, P = 0.88) while purity itself dominates (β = −0.88, P = 1.9 × 10⁻⁴).

**Critical caveat:** the first pass wrote that purity is "definitionally not a confounder
because it is not associated with the outcome". **That sentence is withdrawn.** P > 0.05 is
not evidence of no association. The only defensible statement is that adjustment for the
available purity estimate did not attenuate the association. The purity variable used is a
CIBERSORT leukocyte fraction, not ABSOLUTE tumour purity.

**Confidence:** moderate.

---

## 5b. Does GSE138042 support the association?

**Verdict:** No. It neither supports nor refutes, and one number previously reported from it
was answering a different question.

**Evidence:** the "composition-adjusted OR 0.288, P = 0.020" was fitted on all 95 libraries
with the outcome "is this library named `RAIR-*`", a comparison whose control arm contains
17 benign follicular adenomas and 3 medullary carcinomas. Applying the same covariates to
the paper's actual contrast (13 refractory versus 10 radiosensitive) gives OR 0.442
(0.091–2.140) P = 0.31 unadjusted, 0.204 (0.024–1.699) P = 0.14 with immune content, and
0.205 (0.023–1.838) P = 0.16 with immune and stromal content. **None is significant.**

**Critical caveat:** group membership is perfectly confounded with annotation batch. Every
refractory sample carries a histology string no other sample uses — "Papillary **thyroid**
cancer" (7), "Follicular **thyroid** cancer" (5), "Poorly differentiated thyroid cancer" (1)
— against "Papillary cancer" (35) and "Follicular cancer" (17) elsewhere. Library naming
(`RAIR-*` versus `TC-*`) points the same way. The two arms were annotated separately and
plausibly processed separately.

**Recommended manuscript wording:** "In an independent cohort with a documented radioiodine
outcome, refractory carcinomas scored lower than radiosensitive carcinomas but the
difference was not significant (Cohen's d = −0.42, 95% CI −1.44 to +0.42, P = 0.34), and
group membership was confounded with annotation batch."

**Confidence:** high (that it is uninformative).

**Detail:** `04c_GSE138042_ESTIMAND.md`

---

## 5c. Is GSE173248 usable?

**Verdict:** No for panel testing; yes as an inventory example.

**Evidence:** serum small-RNA sequencing — 2,580 mature miRNAs, 7,485 small-RNA loci, and no
mRNA of any kind. The 22 GEO samples are **11 patients × 2 timepoints** (basal versus 72 h
post-rhTSH), not 22 patients.

**Critical caveat:** three mutually inconsistent cohort sizes exist (GEO overall design says
37 patients, the paper says 30, the deposit contains 11) and the submission has **transposed
field names** — the field labelled `rhtsh stimulation:` holds the disease-status group and
the field labelled `cancer recurrence:` holds the TSH timepoint. Radioiodine treatment is a
background statement in the paper, not a variable in the deposit, so it must not be called a
radioiodine cohort.

**Confidence:** high.

**Detail:** `16_SCREENED_EXCLUSIONS.tsv` row EXC-12

---

## 6. Can we say the panel is superior to driver class?

**Verdict:** No. No head-to-head comparison exists.

**Evidence:** pooled random-effects analysis across Siraj 2022 (158), Zhang 2026 (113) and
Boucai 2023 (23) — BRAF/RAS-negative versus rest, OR 0.75, 95% CI 0.30–1.87, P = 0.53,
I² = 63%.

**Critical caveat:** three studies only, heterogeneous endpoint definitions (seven-criterion
refractoriness / scan avidity / RECIST), a wide interval that does not exclude moderate
associations in either direction, and the smallest cohort (Boucai, n = 23) is the only
nominally significant one and points the **opposite** way. **The claim "the driver axis is
dead" and the claim "all three cohorts were consistent" are both withdrawn.** No analysis
compared driver class against the eight-gene score in the same patients.

**Recommended manuscript wording:** "Across three clinically heterogeneous cohorts, driver
class did not consistently discriminate radioiodine-refractory disease, with substantial
between-study heterogeneity (I² = 63%)."

**Confidence:** moderate.

---

## 7. Does Zhang 2026 validate our panel?

**Verdict:** No. It is orthogonal conceptual support, not validation.

**Evidence:** within 113 advanced-DTC patients, proteomic consensus subtype associates with
refractoriness at χ² P = 1.4 × 10⁻⁷ (CC1 21% → CC2 57% → CC3 87%; CC1 vs CC3 OR 24.6,
P = 2.0 × 10⁻⁸) while driver class does not (P = 0.109) and histology does not (P = 0.248).

**Critical caveat:** the subtypes are whole-proteome consensus clusters, not the eight-gene
transcriptomic score, and the individual protein values are behind iProX. Using one cohort
and one endpoint removes cohort-level confounding but not patient-level confounding —
disease burden, metastatic site, stage, prior treatment and cellularity are not adjusted.
**The claim "cannot be explained by confounding" is withdrawn.**

**Recommended manuscript wording:** "Within the same advanced-disease cohort and clinical
endpoint, proteomic subtype showed a substantially stronger unadjusted association with
refractory disease than driver class."

**Confidence:** moderate.

---

## 8. Can the current manuscript title stand?

**Verdict:** No, if it contains "predicts radioiodine-refractoriness".

**Evidence:** the initial-response test is null under all five index rules, and no
treatment-by-score interaction is estimable. "Predicts" commits to a temporal predictive
design that the data do not support.

**Critical caveat:** the structural-outcome association is real but prognostic within a
treated population, which supports a risk-marker title, not a response-predictor title.

**Recommended manuscript wording:** see `07_MANUSCRIPT_CHANGE_PLAN.md` for title candidates.
The title should not be changed automatically; it is a decision for the author.

**Confidence:** high.

---

## 9. How much radioiodine analysis belongs in the current manuscript?

**Verdict:** A bounded amount — the corrected null, the structural association with its
rare-event caveats, and the external driver evidence in Discussion. Not the dataset atlas.

**Evidence:** the DM1 manuscript is about BRAF/RAS-negative papillary carcinoma; the external
cohorts are advanced/metastatic disease. Scope mismatch is the main risk of full integration.

**Confidence:** moderate.

---

## 10. Does a separate radioiodine paper stand up?

**Verdict:** Yes, on today's data alone.

**Evidence:** four assets that do not depend on new data — the recovered TCGA treatment layer
(235 patients, invisible in cBioPortal), a well-characterised null with its ceiling-effect
explanation, two demonstrations of how public radioiodine associations fail (composition
confounding in GSE151179, comparator contamination in GSE138042), and a design prescription.

**Critical caveat:** it is a resource-and-reappraisal paper, not a biomarker paper. It must
not be written as if it validates the eight-gene score.

**Confidence:** moderate.

---

## 11. What additional data is most valuable now?

**Verdict:** The Zhang individual proteomics matrix (iProX IPX0011848000) plus the matched
GSA sequencing.

**Evidence:** 113 advanced-DTC patients, 51% refractory, with a documented radioiodine
endpoint — the only identified cohort in the right disease setting and near the required
size. One contact, Xiao Shi, is both the paper's lead contact and the GSA data-access
committee contact.

**Critical caveat:** even 113 patients falls short of the ≈226 the observed effect sizes
require under equal allocation.

**Confidence:** high.

---

## 12. The five things to do in the next seven days

1. Read the P0 rows of `06_CLAIMS_REGISTRY.csv` and rewrite Limitation 6, which is factually
   wrong and checkable by any reviewer with GDC access.
2. Decide Strategy A versus B in `07_MANUSCRIPT_CHANGE_PLAN.md` before writing any new prose.
3. Open iProX IPX0011848000 in a browser with a free account and capture the file listing.
4. Pull the Chirra 2026 *Clin Cancer Res* full text through the SNU library and read only the
   Methods sentence on whether radioiodine-refractory status is a coded CODEai variable.
5. Send nothing yet. The four request drafts in
   `rai-response-genomics-atlas/docs/data_requests_drafts_2026_08_06.md` need the corrected
   effect sizes substituted before they go out — they currently quote OR 0.19.

---

## Claim tiers

**Tier 1 — directly supported by today's audited analyses**
- TCGA-THCA contains per-course radioiodine dose and response records, absent from cBioPortal.
- No detectable association between the score and initial recorded response in the available
  predominantly complete-response cohort.
- Lower differentiation score is associated with selected post-treatment structural endpoints
  within radioiodine-treated TCGA patients, subject to rare-event limitations.
- Driver class does not consistently classify refractory disease across three heterogeneous
  external cohorts.
- Proteomic subtype associates strongly with refractory disease in Zhang 2026.

**Tier 2 — conceptually supported, not directly validated**
- Molecular differentiation state may be more informative than driver class.
- Advanced/metastatic disease is the appropriate validation setting.
- Tissue composition is a major source of false-positive association in this literature.

**Tier 3 — not currently supported, must not appear**
- The eight-gene score predicts treatment benefit.
- The score is a validated radioiodine-specific predictive biomarker.
- Zhang 2026 directly validates the eight-gene panel.
- Driver mutations are biologically irrelevant to radioiodine handling.
- Methylation × radioiodine response data exist in Rodríguez-Lloveras 2025 — **confirmed
  false**, see `16_SCREENED_EXCLUSIONS.tsv` row EXC-01.
