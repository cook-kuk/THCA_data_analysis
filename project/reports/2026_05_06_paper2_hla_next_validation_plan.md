# Paper 2 HLA — next validation plan
**Date:** 2026-05-06
**Goal:** Define the minimum validation work that lets us graduate from "exploratory candidate-prioritization map" (current state) to a publishable Korean-PTC HLA-association statement.

---

## 0. Premise

Current best-allowed claim:
> "The six-allele HLA forest provides an exploratory candidate-prioritization map suggesting Korean PTC may have immune-genetic background signals requiring matched-control validation."

After strict metric harmonization (S3 direct allele-vs-allele Fisher, BH-q within 6):
- **0/6** alleles support an enrichment claim (A*02:07, B*46:01, DPB1*05:01, DRB1*07:01 fail q < 0.05 under strict comparison).
- **2/6** alleles show a depletion direction surviving q < 0.05: **C*01:02** (q = 7e−4) and **DQB1*02:01** (q = 5e−8, zero-cell caveat).

The validation plan must therefore (a) replicate the depletion signals on independent cohorts with matched controls, and (b) clarify whether any of the four S0-only enrichment signals are real but suppressed by HWE-violation/baseline-source-mismatch.

---

## 1. Required validation work (gate to graduate)

### 1.1 Matched-control NGS HLA typing (highest priority)
- **Population:** Independent Korean healthy adults (not cancer-patient adjacent normal — those are paired tumor controls, not population controls).
- **Typing pipeline:** Same as PTC pool (arcasHLA from RNA-seq OR equivalent NGS-based 4-digit typing). Resolution must match.
- **Minimum n:** ≥ 1,000 individuals (so each cohort has > 100 expected carriers for the 5%–35% AF range we care about).
- **Output schema required:** per-individual 2-allele genotype (not aggregated frequency). This lets us compute carrier frequency *and* allele frequency directly on both sides — eliminating metric-mismatch entirely.

### 1.2 Allele-level dosage harmonization (analytical)
- Both PTC and control sides must be expressed as *allele copy count over 2n_callable* in addition to carrier counts.
- The primary Fisher comparison must be allele-level (S3-equivalent) on both sides.
- Carrier-level comparison is reported as a sensitivity check, not as the primary metric.

### 1.3 Multi-cohort replication
- Minimum **2 independent Korean PTC cohorts** with the same NGS HLA typing pipeline.
- Each cohort analyzed independently; pooled meta-analysis only after per-cohort direction stability is confirmed.
- Heterogeneity test (Cochran's Q or I²) reported; subcohort effects flagged if I² > 50%.

### 1.4 Pre-registration
- Hypothesis must be specified per allele *and per direction* (e.g., "C*01:02 PTC depletion"; "DQB1*02:01 PTC depletion") before counting.
- Multiple-test correction rule (BH within candidate set or full-genome) decided in advance.
- Effect size threshold (minimal important difference, MID) declared in advance.

---

## 2. Graduation criteria (per allele)

An allele can graduate from "exploratory candidate" to "supported Korean-PTC HLA association" only when **all** of the following hold:

1. **Direction stability** — same direction across S0, S1, S2, S3 metric scenarios.
2. **Strict-metric q < 0.05** — S3 direct allele-vs-allele Fisher BH-q < 0.05 (within the pre-registered candidate set, or genome-wide-corrected if all-allele tested).
3. **Multi-cohort replication** — direction reproduces in ≥ 2 independent Korean PTC cohorts at p < 0.05.
4. **Effect size > MID** — pre-specified MID; 95% CI of OR excludes 1.0 in the matched-control comparison.
5. **No source-sensitivity collapse** — if the baseline is a meta-pool (like AFND), at least one alternative individual-cohort baseline must reproduce direction.
6. **No zero-cell collapse** — if the signal is driven by a zero-observation cell, the zero must reproduce in matched-control NGS (same pipeline).

Only after all six criteria pass for an allele can the manuscript / web page upgrade language from "exploratory candidate" to "supported Korean PTC HLA association."

---

## 3. Action items (next 4 weeks, by priority)

| # | Action | Owner | Eta | Dependency |
|---|---|---|---|---|
| 1 | Identify a public or in-house Korean healthy-population NGS HLA dataset matching PTC pipeline | Cook + Yu lab | 1 wk | none |
| 2 | If no public match, scope the cost / consent / typing pipeline for a small (n=200) pilot matched-control set | Cook | 2 wk | (1) negative |
| 3 | Re-run S3 direct Fisher on Kim 2014 reference panel as a *partial* matched control (already has F17_kim2014_reference_panel_carrier_validation_forest.png; needs S3 allele-level extension) | Cook | 1 wk | none |
| 4 | Run S3 direct Fisher on Baek 2021 PLOS NGS baseline (already has F19_allele_vs_baek2021_ngs_forest.png at carrier level; extend to allele level) | Cook | 1 wk | none |
| 5 | Pre-register hypothesis + analysis plan for matched-control validation (focus: C*01:02 depletion, DQB1*02:01 depletion) | Cook + Yu | 2 wk | (3) and (4) |
| 6 | Document per-pipeline 4-digit typing call rate per locus on PTC pool (callable n already varies: A=874, DPB1=779, DQB1=631) | Cook | 1 wk | none |
| 7 | Audit DQB1*02:01 zero-call: verify it isn't an arcasHLA imputation specific failure by spot-checking high-coverage RNA-seq samples | Cook | 2 wk | none |
| 8 | Decide whether genome-wide HLA testing (all 4-digit alleles passing call-rate threshold) is feasible from existing PTC pool data | Cook + Yu | 2 wk | (6) |

Actions 3, 4, 6, 7 can begin immediately with locally available data.

---

## 4. What we will *not* do during this phase

- We will **not** publish or pre-print the 6-allele forest as a stand-alone case-control association paper. It is a candidate-prioritization map only.
- We will **not** combine the HLA forest with the 8-gene RAI-lineage readout into a "joint immune-genetic risk model." Conditional independence has not been tested; no integrated risk model is supported.
- We will **not** expand to GWAS / non-HLA SNP analysis until the HLA-only validation is settled.
- We will **not** perform unblinded analyses on any new control cohort before the analysis plan is pre-registered.
- We will **not** touch Paper 1 (driver-orthogonal differentiation axis) or Paper 9 (perturbation / synthetic lethality) deliverables in service of Paper 2 HLA.

---

## 5. Connection to Paper 9 (future-work only)

If C*01:02 depletion is replicated in matched-control validation, the resulting hypothesis — that some Korean PTC tumors arise in HLA-C-restricted immune contexts — can motivate testing whether HLA-low samples co-occur with synthetic-lethality candidates from Paper 9 (e.g., GLS, metabolic dependencies). This is a *future-work* link only; no quantitative analysis is performed in this packet.

---

## 6. Status snapshot

| Item | Status |
|---|---|
| 6-allele exploratory forest (S0) | ✅ done, exploratory only |
| Metric-sensitivity (S0/S1/S2/S3) | ✅ done; see `paper2_hla_metric_sensitivity_strengthened.tsv` |
| BH-q within 6 per scenario | ✅ done; see `paper2_hla_bh_fdr_by_scenario.tsv` |
| Subcohort heterogeneity | ✅ done; F17_subcohort_frequency_heatmap.png |
| Direct allele-vs-allele Fisher (S3) | ✅ done; `paper2_hla_allele_vs_allele_direct.tsv` |
| Claim-boundary matrix | ✅ done; F18_claim_boundary_matrix.png |
| Validation roadmap diagram | ✅ done; F19_validation_roadmap.png |
| Matched-control NGS dataset acquisition | ❌ pending (Action 1) |
| Pre-registered validation plan | ❌ pending (Action 5) |
| Multi-cohort replication | ❌ pending |
| Genome-wide HLA testing | ❌ pending (Action 8) |

---

**End of validation plan. No association claim is implied or attached to any of the actions above.**
