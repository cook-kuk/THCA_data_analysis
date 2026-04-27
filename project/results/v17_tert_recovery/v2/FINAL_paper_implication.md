# v17 TERT recovery — paper implication

_Scenario A — 36 TCGA-THCA TERT promoter mutations recovered_

## Bottom line

The 4-group narrative (BRAF / RAS / +TERT / triple-neg) is now defensible
on TCGA-THCA. Survival separation is highly significant (logrank
p = 4.9e-6 TERT+ vs WT, p = 3.8e-5 four-group), driven by 6/36 events in
TERT+ vs 10/468 in WT (~8x event rate).

## Concrete next actions for v17 paper integration

1. **Sample master update**
   - Replace `tert_promoter` / `tert_status` columns in
     `project/results/v17/tables/sample_master_v17_tert.tsv` from
     `v2/sample_master_v17_tert_v2.tsv` (column `tert_promoter_integrated`).
   - 36 rows flip from "wildtype" / "" to "mutated" with provenance
     `thca_tcga_pub` (cBioPortal Sanger-validated MAF).

2. **Main-text figure (proposed)**
   - Panel A: 4-group KM curves (BRAF only / RAS only / TERT+ / triple-neg),
     with shaded 95% CI bands and at-risk table.
   - Panel B: TERT+ rate by quad-group bar chart (8.9% BRAF, 11.1% RAS, 2.8%
     triple-neg) — note absolute n's because absolutes drive interpretability.
   - Panel C: TDS16 differentiation score box plot, TERT+ vs WT (Δ = -0.64).

3. **Supplementary figure**
   - DM1 vs DM2 enrichment for TERT+ (4 vs 1 in DM subcohort, Fisher
     p = 0.65). Document as exploratory only because it's underpowered.
   - Stage distribution stacked bar (73% III/IV in TERT+).

4. **Survival CIs — present honestly**
   - 6 events in 36 TERT+ → wide HR CI. Use bootstrap or Cox with Firth
     correction. State n and event count next to every HR.

5. **External validation paragraph**
   - 1,361 MSK-IMPACT, 342 MSK-CHORD, 81 MSK PDTC/ATC TERT promoter records
     available. These cannot be patient-matched to TCGA but support a cohort-
     level prevalence sanity check (TERT+ is ~10-30% in advanced thyroid;
     our 7% in TCGA-THCA is consistent with well-differentiated TCGA
     enrichment).

## Honest limitations to add

- TERT calls come from the original 2014 *Cell* paper Sanger validation,
  not new sequencing. We did not re-call from BAMs.
- Liu 2017 *JCO* supplement (potentially 100+ additional patients) was not
  recovered after 8 URL attempts. Future versions should retrieve manually.
- Yoo 2019 Korean cohort and Pozdeyev 2024 Nature Cancer paper not
  recovered automatically; both would meaningfully extend external
  validation if obtained later.

## Probability deltas (informal)

- npj submission: ~70% → ~78% (4-group narrative now defensible, survival
  signal strong, but only a single TCGA cohort).
- Genome Med: ~40% → ~50% (TERT recovery + survival is meaningful but not
  uniquely positioning).
- Reviewer attack vector "you ignored TERT" is closed by the audit trail in
  `FINAL_recovery_audit.md` regardless of which venue.

## Outstanding follow-ups

1. Retry Liu 2017 supplement via manual library / interlibrary loan request.
2. Reach out to Pozdeyev / Fagin labs at MSK for Pozdeyev 2024 supplement.
3. Korean cohort (Yoo 2019): direct request to Yoo SK or KOBIC.
4. dbGaP DAR for TCGA-THCA WGS BAMs (470 cases listed; see
   `S7_gdc_controlled_application_guide.md`) — long timeline, not for
   current sprint.
