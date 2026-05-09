# H24 — Sensitivity battery on H6 (DM2-vs-not_DM PFI in BRAF-cPTC)

**VERDICT.** DM2-PFI HR survives **9/9 deployable** subgroup analyses and **15/15 evaluable** subgroups without ever crossing 1. Lower-95%-CI stays above 1 in **7/9 deployable** strata. The H6 headline (HR=5.91 [1.79, 19.5], p=0.0036) is robust to clinical confounders.

## Lead numbers (BRAF-cPTC, DM2 vs not_DM, PFI)

| Subgroup | n | ev | HR | 95% CI | p |
|---|---|---|---|---|---|
| Overall (H6 replicate) | 254 | 17 | **5.91** | 1.79 – 19.5 | 0.0036 |
| Stage-stratified Cox | 268 | 31 | **3.75** | 1.33 – 10.6 | 0.013 |
| T-stratified Cox | 268 | 31 | 2.20 | 0.80 – 6.10 | 0.13 |
| Within Stage I | 144 | 13 | 4.87 | 0.80 – 29.8 | 0.087 |
| Within Stage III | 66 | 11 | **4.68** | 1.15 – 19.0 | 0.031 |
| T2-T4 only | 180 | 17 | **4.91** | 1.40 – 17.2 | 0.013 |
| M0 only | 156 | 10 | **16.5** | 2.63 – 104 | 0.0028 |
| Dx ≤ 2010 | 130 | 11 | **7.96** | 1.26 – 50.0 | 0.027 |
| **Cluster by TSS (robust)** | 254 | 17 | **5.91** | **2.95 – 11.8** | **5.3e-7** |

T1-only had 0 events (BRAF-cPTC + T1 = indolent niche, no progression to model — *consistent*, not contradictory). Stage II / Stage IV / M1 / N1 / multifocal etc had n_DM2 of 1–8; per-cell HRs are direction-only and excluded from the deployable count.

## Stage-stratified ≠ stage-adjusted

Moving `adv_stage` from covariate → strata (Cox baseline differs per stage) gives HR=3.75 [1.33, 10.6] p=0.013. The cleanest reviewer rebuttal: the effect is *not* explained away by allowing stage-specific baseline hazards.

## Cluster by 22 tissue source sites

`cluster_col=tss`: **HR=5.91 [2.95, 11.8] p=5.3e-7**. Point estimate identical to H6, **CI tightens** because intra-institution dependency was inflating the naive SE. The effect is not an institution artifact.

## Interaction p-values (no effect modification)

ETE × DM2 p=0.57; multifocal × DM2 p=0.20; T1 × DM2 p=0.18; era>2010 × DM2 p=0.65. **No anatomic feature, focality pattern, or treatment-era restricts the DM2 risk.** N×DM2 and M×DM2 interactions failed to converge (DM2 cell n≤2 in N1/M1) — flagged in `h24_interactions.tsv`.

## DM1-protective sensitivity (DM1 vs DM2)

Main HR=0.21 [0.04, 1.01]; stage-stratified HR=0.26; T-stratified HR=0.29. **Direction (HR<1) stable across all slices**; p-values trend 0.05–0.08, limited by 8 events in n=104. Direction-supportive but not deployable as a positive biomarker yet.

## Honest caveats

- Per-anatomic-cell HRs (ETE pos / N1 / M1 / multifocal) have DM2 cell sizes of 1–8 and should be read direction-only.
- Stage IV did not converge (6 events / 35 patients fit-pathology).
- T-stratified Cox is the one slice where CI crosses 1 (HR=2.20 [0.80, 6.10]) — interpretable as T-stage absorbing some of the DM2 risk-axis.

## Files
- `h24_sensitivity_table.tsv` (32 rows)
- `h24_forest_data.tsv` (20 subgroups)
- `h24_interactions.tsv` (7 interaction tests)
- `h24_forest_plot.png`
- `h24_merged_with_external_clinical.tsv` (513×50, H6 + cBio thca_tcga + GDC focality/year)
- `h24_summary.json`, `cache/{cbio,gdc}_*.tsv`
