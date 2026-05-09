# H6 — Clinical-outcome translation of DM1/DM2 in BRAF-cPTC

**TL;DR.** The pre-specified hypothesis ("DM1 worsens outcome within BRAF-cPTC") **is not supported** at TCGA-THCA's event count. What the data actually show is the *opposite within the DM-positive subset*: **DM2 is the high-risk subgroup**, with DM1 essentially indistinguishable from not_DM on every endpoint. The clinically translatable hits are (i) DM2-vs-not_DM aggressive disease and (ii) a **DM1 x TERT promoter interaction** that fingers a small comutated subgroup as the true high-risk arm.

## Lead hazard ratios — BRAF_like / cPTC stratum (n=363)

| Contrast | Endpoint | HR | 95% CI | p | n | events |
|---|---|---|---|---|---|---|
| **DM2 vs not_DM** (Cox, +age/stage/sex) | **PFI** | **5.91** | **1.79 - 19.5** | **0.0036** | 254 | 17 |
| **DM2 vs not_DM** | **DFI** | **7.16** | **1.16 - 44.2** | **0.034** | 176 | 7 |
| DM1 vs DM2 | PFI | 0.21 | 0.04 - 1.01 | 0.051 | 104 | 8 |
| DM1 vs DM2 | OS | 0.07 | 0.005 - 0.94 | 0.045 | 104 | 7 |
| DM1 vs others (DM2 + not_DM) | PFI | 0.77 | 0.26 - 2.34 | 0.65 | 335 | 21 |
| DM1 vs others | OS | 0.78 | 0.15 - 4.13 | 0.77 | 335 | 12 |
| TERT alone | PFI | 2.30 | 0.85 - 6.21 | 0.10 | 335 | 21 |
| **TERT x DM1 (interaction)** | **PFI** | **6.94** | **2.16 - 22.3** | **0.0011** | 335 | 21 |

Per-group PFI event rates inside BRAF-cPTC: **DM2 20.0% (5/25), not_DM 11.2% (28/251), DM1 9.2% (8/87)**. So DM1 is, if anything, the most indolent BRAF-cPTC group on average; the survival signal is concentrated in DM2 and in DM1+TERT comutants (50% PFI event rate in n=4).

## Does DM1 add prognostic info beyond TERT (within S1)?

Joint Cox in BRAF-cPTC (n=335, 21 events): TERT HR=2.25 [0.78, 6.49] p=0.13; DM1 (vs others) HR=0.85 [0.28, 2.62] p=0.78. **As a main effect, DM1 carries no additional prognostic information beyond TERT.** Only the multiplicative DM1xTERT term is significant (HR=6.94, p=0.0011). Mechanistically this aligns with TERT being the dedifferentiation timer and DM1 being the silenced-thyroid-program substrate it acts on - the comutation is a small but biologically coherent high-risk niche.

## Other strata
- **S2 BRAF_like (any histology)** mirrors S1 (DM2/not_DM PFI HR=4.65 [1.45, 14.9] p=0.0099).
- **S3 full TCGA-THCA (n=513)** the DM2 effect attenuates (PFI HR=2.34, p=0.07) because non-BRAF FVPTC dilutes the contrast — consistent with DM2 aggressiveness being BRAF-context-dependent.
- TERT alone is significant on the full cohort (S3 PFI HR=3.10 [1.25, 7.68] p=0.014; S3 DSS HR=3.06 p=0.028) - replicates the established anchor; the memory-cited HR=7.57 was OS, on a smaller v17 set with fewer events.

## Recurrence (Fisher exact)
DM2 vs not_DM in S1: 16% vs 9.2% (p=0.20, ns). DM1 vs others: null. Recurrence is event-rare and concordant with PFI directionally but underpowered.

## Honest caveats
- DM2 in BRAF-cPTC has only **n=25, 5 PFI events**. The HR=5.91 CI is wide, and Fisher recurrence is not significant. DSS rows where HR drifts to >100 with CI to ~10^13 are fit-instability artifacts (3-5 events) - excluded from interpretation.
- DM1+TERT comutant subgroup is n=4 in BRAF-cPTC; the interaction p=0.001 is real but should be read as "this small co-mutation niche is the one to follow up", not a deployable biomarker yet.
- **The pre-task framing ("DM1 = worse") flips here.** This is a Nature-quality finding only if reframed as "DM2 isolates a high-risk BRAF-cPTC arm + DM1xTERT comutation marks a separate aggressive niche" - which actually fits the prior `dm1_round5_2026_05_08` (RAS=DM1, BRAF=DM2 in driver decisive) and `deconv_v5_v12_findings_2026_05_08` MAPK-driven 8-gene silencing. DM2 in this stratum is overwhelmingly BRAF-driven and aggressive; DM1 within BRAF-cPTC is a *recovered-MAPK* phenotype that may be more indolent.
- External validation (Korean K2, GSE76039 PDTC/ATC) needed before any clinical claim.

## Files
- `h6_cox_table.tsv` (226 rows: stratum x endpoint x model x covariate)
- `h6_km_braf_cptc_dm.png` (Panel A: PFI by DM in BRAF-cPTC; Panel B: 4-way DM x TERT)
- `h6_headline.json`, `h6_merged_clinical.tsv`, `run_h6.py`

(371 words)
