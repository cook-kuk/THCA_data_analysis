# AUDIT 19 — RAI-anchored reanalysis, in response to Prof. Kang Minsu

**Date** 2026-08-06 · triggered by an external reviewer objection, not by our own audit plan.

## The objection

> "Progression free survival로 과연 RAI 반응성을 평가할 수 있는가? → 논리적 비약이 있는 것
> 같습니다. RAI 안 하고 갑상선 수술 후 재발한 사람들도 포함되어 있을 것이기 때문입니다."
> — 강민수, 2026-08-06

Two defects are named: **cohort contamination** (non-radioiodine-treated patients contribute
recurrence events) and, implicitly, **time-origin mismatch** (a diagnosis-anchored clock
cannot isolate what happened after radioiodine).

## 1 · The contamination is real, and we measured it

TCGA-THCA panel patients (n = 500) stratified by radiation record:

| Group | n | PFI events | Event rate |
|---|---|---|---|
| Radioiodine (mCi) | 235 | 36 | 15.3% |
| External beam only (Gy/cGy) | 11 | 5 | 45.5% |
| **No radiation record** | 229 | **9** | 3.9% |
| Radiation, unit not recorded | 25 | 1 | 4.0% |

**15 of 51 total PFI events (29.4%) occur in patients who never received radioiodine.**
The objection is quantitatively correct for any analysis run on the full cohort.

**Our analyses were not contaminated** — `tcga_rai_best_response_2026_08_06.py:111` performs
an inner join onto patients with at least one mCi-dosed course, so the cohort has been
restricted to n = 235 throughout. **But this was never stated in the manuscript**, which
means a reader has no way to distinguish our cohort from the contaminated one. That is a
reporting defect and must be fixed.

## 2 · The time-origin defect was real and we had it

Restricting the cohort does not fix the clock. A diagnosis-anchored interval bundles surgery,
the pre-radioiodine waiting period, and post-radioiodine follow-up into one number.

The GDC BCR Biotab radiation file records `radiation_therapy_started_days_to` for
**all 235 patients**. Re-anchoring:

    t = new_tumor_event_dx_days_to − radiation_therapy_started_days_to

| Quantity | Value |
|---|---|
| Analysis set | 231 |
| Recurrence events after radioiodine | 35 |
| Events occurring **before** the first radioiodine course | **0** |
| Median follow-up after radioiodine | 25.3 months |

## 3 · Result: null, and the one nominal signal does not survive

| Analysis | n | events | Result |
|---|---|---|---|
| **Cox, score continuous, RAI-anchored** | 231 | 35 | **HR 0.710 [0.398, 1.266], P = 0.246** |
| Cox + stage | 231 | 35 | HR 0.733, P = 0.289 |
| Cox + age | 231 | 35 | HR 0.722, P = 0.265 |
| Cox + cumulative dose | 231 | 35 | HR 0.698, P = 0.231 (dose P = 0.050) |
| Sensitivity: drop external-beam co-treated | 225 | 34 | HR 0.705, P = 0.260 |
| (contrast) Cox, diagnosis-anchored PFI | 235 | 36 | HR 0.761, P = 0.337 |

### The median split is an artefact of the cut-point — RETRACTED

An earlier run of this analysis reported a log-rank **P = 0.043** for a median split. That
number must not be used.

| | |
|---|---|
| Median split, as first reported | χ² = 4.083, **P = 0.043** |
| **Maximally selected log-rank over all cut-points** | χ²max = 6.437 at cut −0.2607 |
| **Permutation-corrected P (5,000 permutations of the score)** | **P = 0.142** |

The median is one of ~180 candidate cut-points. Selecting it is a free parameter, and once
that freedom is paid for, the signal disappears. **The correct statement is that the
RAI-anchored analysis is null.**

Supporting: tertile trend P = 0.099 (event rate 19.5% → 15.6% → 10.4%, monotone but not
significant); quartile trend P = 0.144 and **not** monotone (17.2%, 22.4%, 8.8%, 12.1%).

### Proportional hazards is violated

Formal Schoenfeld test (`lifelines.statistics.proportional_hazard_test`, rank time-transform):
**test statistic 7.481, P = 0.0062.** An in-repo Spearman check of score residuals against
event time agreed (ρ = 0.456, P = 0.006).

The hazard ratio is therefore not a stable summary over follow-up, and any Cox HR from this
cohort — including the null ones above — should be reported with that caveat rather than as a
single constant effect.

## 4 · A second TCGA source, and a conflict inside TCGA

The full BCR Biotab clinical supplement (downloaded 2026-08-06, 9 files, 883 KB) carries a
field the radiation file does not:

    i_131_radiation_first_tx_method ∈ {Thyroxine withdrawal, rhTSH,
                                       "Patient did not receive I-131 treatment"}

| | radiation file has mCi course | no mCi course |
|---|---|---|
| **explicit YES I-131** (withdrawal 59 / rhTSH 44) | 101 | 2 |
| **explicit NO I-131** | **18 ← conflict** | 54 |
| not recorded | 116 | 209 |

**18 patients (7.7% of our 235-patient cohort) are explicitly recorded as not having
received I-131 while carrying a millicurie-dosed radiation course.** TCGA contradicts itself.
This has to be disclosed, not silently resolved.

Cohort-rule sensitivity — the conclusion does not depend on how the conflict is resolved:

| Rule | n | events | HR | P |
|---|---|---|---|---|
| Permissive (any mCi course) | 231 | 35 | 0.710 | 0.246 |
| Strict (mCi course **and** not explicitly denied) | 214 | 32 | 0.644 | 0.153 |
| Explicit-only (explicit YES **and** mCi course) | 100 | 11 | 0.659 | 0.442 |

## 5 · The treatment × score interaction is estimable for the first time — and is not significant

The explicit field supplies a genuine untreated comparator rather than one inferred from a
missing record. This is what turns a prognostic question into a predictive one.

| Stratum | n | PFI events |
|---|---|---|
| Explicit I-131 | 103 | 11 |
| Explicit NO I-131 | 72 | 8 |

| | Estimate |
|---|---|
| **Interaction, score × treatment** | **HR 0.360 [0.087, 1.484], P = 0.157** |
| Within explicit I-131 | HR 0.689, P = 0.478 |
| Within explicit NO I-131 | **HR 1.930**, P = 0.189 |

The stratum-specific estimates point in **opposite directions**, which is the shape a
treatment-specific effect would produce. But with 19 events total, neither stratum is
significant and the interaction is not significant. **This is exploratory and must be
labelled as such.** It is a reason to seek an untreated comparator with more events — not a
predictive claim.

## 6 · TSH preparation method — an efficacy covariate nobody models

Thyroxine withdrawal 59 versus rhTSH 44. Preparation method affects iodine uptake, so it is a
plausible modifier of radioiodine efficacy and is absent from every analysis we have seen.

Cox on the RAI-anchored clock: score HR 0.636 (P = 0.411); rhTSH versus withdrawal
HR 2.201 (P = 0.223). Underpowered (11 events), directionally worth recording, not a finding.

## What changes in the manuscript

1. **State the cohort restriction explicitly.** "Restricted to the 235 patients with a
   millicurie-dosed radiation course" must appear wherever a radioiodine-related result does.
2. **Do not cite the median-split P = 0.043 anywhere.** It is retracted here.
3. **Report the RAI-anchored analysis as null**, with the diagnosis-anchored version shown
   only as a contrast.
4. **Disclose the 18-patient TCGA conflict** and show the three cohort rules.
5. **Keep "prognostic within treated patients"; do not write "predictive".** The interaction
   is now estimable and is not significant.
6. **Add the proportional-hazards violation** to limitations.

## Files

- `scripts/audit/08_time_from_rai_to_recurrence_2026_08_06.py`
- `scripts/audit/09_rai_anchored_robustness_2026_08_06.py`
- `scripts/audit/10_explicit_i131_label_2026_08_06.py`
- `results/tables/audit08_rai_anchored_recurrence_2026_08_06.tsv`
- `results/tables/audit08_rai_anchored_cohort_ledger_2026_08_06.tsv`
- `results/tables/audit09_rai_anchored_robustness_2026_08_06.tsv`
- `results/tables/audit10_explicit_i131_concordance_2026_08_06.tsv`
- `results/tables/audit10_explicit_i131_results_2026_08_06.tsv`
- `external_data/tcga_thca_biotab_full/` (9 files, SHA-256 manifest)

## Verification performed

**Permutation stability.** The corrected P was recomputed under five independent seeds
(20260806, 1, 42, 777, 20250101): **0.1420, 0.1390, 0.1404, 0.1490, 0.1484** — range
0.139–0.149, every value far above 0.05. The retraction does not depend on the seed.

**Cox and log-rank implementations cross-checked against `lifelines` 0.30.3**, which was not
used to produce any result above:

| Model | in-repo | lifelines |
|---|---|---|
| score only, RAI-anchored | HR 0.7103, P 0.2463 | HR 0.7103, P 0.2463 |
| score + stage | HR 0.7326, P 0.2887 | HR 0.7326, P 0.2887 |
| score + age | HR 0.7225, P 0.2651 | HR 0.7225, P 0.2651 |
| interaction term | HR 0.3599, P 0.1574 | HR 0.3600, P 0.1574 |
| log-rank, median split | χ² 4.0825, P 0.0433 | χ² 4.0825, P 0.0433 |

**Cohort counts recomputed from the raw files with no shared helper code**: total PFI events
51, of which 15 (29.4%) in non-radioiodine-treated patients; explicit NO I-131 = 72;
explicit YES = 103; conflict = 18 (7.7% of the 235-patient mCi cohort); missed = 2. All match.

seed 20260806 · 5,000 permutations · Cox by Newton-Raphson with Breslow ties, implemented
in-repo and validated against lifelines · no manuscript file modified · no commit made.
