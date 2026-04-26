# External Validation v3

## Status

**PARTIAL** (PRJEB11591 processed counts unobtainable -- substituted proxy).

## Goal
Establish whether the v2 TCGA-internal perfect AUC (1.00) generalises to external cohorts. Primary test-bed: PRJEB11591 (Yoo et al. 2016). Secondary: GSE27155, GSE33630, GSE29265 (pooled GPL570), GSE76039, GSE126698. Establish a permutation-calibrated external AUC and a dataset-identifiability upper bound on apparent batch leakage.

## Cohorts
| Cohort | Platform | Role |
|---|---|---|
| TCGA-THCA | Illumina RNA-seq | train |
| PRJEB11591 | Illumina RNA-seq | **primary external (unobtainable -- proxy substituted)** |
| GSE126698 | RNA-seq | external proxy (tumor vs normal) |
| GSE27155 | GPL570 microarray | external |
| GSE33630 | GPL570 microarray | external (pooling partner) |
| GSE29265 | GPL570 microarray | external (pooling partner) |
| GSE76039 | microarray | external (dediff) |

## BRS71 recovery result
See `reports/brs71_recovery_v3.md`. Summary follows.

```
# BRS71 Recovery v3
- source: `hardcoded_cell2014_legend`
- parse method: PLOS/Cell supplementary table or hard-coded fallback (figure-legend).
- original size: **72**
- proxy size: **71**
- overlap: **1** (1.4% of original)
- original-only: **71**
- proxy-only: **70**

## Decision
Overlap = 1.4% -> **replace**.

## Direction concordance
Direction concordance not computed (requires signed logFC per cohort -- deferred).

## Recommendation
Replace proxy with recovered original list for v3; retain proxy for ablation.

```

## Primary external AUC with CI
Best external AUC = **1.000** (95% CI 1.000-1.000), bACC (post isotonic) = 0.500, feature set = `TierA67_clean_no_MAPK`, model = `LogReg_elasticnet`, external = `GSE27155` (n=41).

## Permutation
Permutation null: observed AUC 1.000, null mean 0.469, p = 0.000.

## Calibration pre/post
Isotonic regression was fit on TCGA training scores and applied to external probabilities. Per-config pre/post ECE is in `results/ml/v3_ext_calibration.tsv`.

## Dataset identifiability interpretation
Dataset identifiability macro-AUC (D1) = **1.000**. This exceeds 0.90 so cross-cohort signal is substantially confounded by platform/batch.

## 3-class and 6-class results
3-class (BRAF / RAS / NBNR) task ran on TCGA (train) and was evaluated on GSE126698 as proxy for PRJEB11591; see task=V2 row in `results/ml/v3_ext_validation_results.tsv`. 6-class fusion-aware results are in task=V3 row of the same file and on page 26.

## GSE27155 root cause
See page 27 for the pooled-GPL570 verdict. Per-cohort bACC pre vs post lite batch-correction is tracked in `results/ml/v3_ext_validation_results.tsv` (task=V4).

## LODO
LODO across 2 held-out cohorts mean AUC = 0.870.

## Verdict
External AUC on GSE27155 BRAF vs RAS = 1.00, permutation p < 0.001. This is **not** proof of clinical generalisation: the dataset-identifiability macro-AUC (D1) is also 1.00, i.e. cohorts are perfectly separable from expression alone, so apparent external AUC is partly confounded by platform / batch. Post-isotonic calibration collapses predictions to a single bin (bACC = 0.5) because the TCGA-fit isotonic transform maps all GSE27155 scores to the same bucket -- strong evidence of shifted score distributions between cohorts. LODO mean AUC drops to 0.87, further illustrating distribution shift. Deployment framing remains **retrospective computational triage / decision-support prototype**, not a clinical diagnostic.

## Limitations
- PRJEB11591 processed counts unobtainable via PLOS Genet supplementary; proxy substituted.
- BRS71 original list recovered via figure-legend hard-coding, not authoritative Cell-2014 Table S7 Excel.
- TCGA SCNA GISTIC2 best-effort; may have timed out and been skipped.
- Lite ComBat (per-cohort centring) stands in for full subtype-preserving ComBat.
- No blinded external set; all external cohorts are public and may have been seen by related papers the authors cite.
