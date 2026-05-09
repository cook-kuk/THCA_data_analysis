# H9 — MAPK-inhibitor Sensitivity in DM1 Thyroid Cell Lines

**Headline.** **Dabrafenib (BRAFi) and encorafenib (BRAFi) preferentially hit DM1-high thyroid lines (Cohen's d ≈ −0.6, GDSC2 + PRISM, signs concordant). Trametinib (MEKi) is universally effective across thyroid (median AUC 0.72) regardless of DM1 — consistent with deconv-v12 two-axis convergence.** Within-thyroid n=12–13 is below FDR power.

**Cohort.** DepMap 24Q2, n=22 thyroid lines, panel 25/25 resolved. DM1 = z(HT-13).sum − z(FA-12).sum. DM1-high: HTCC3, BCPAP, MB1, CCLFTHYR-006/-008. DM1-low: ASH3, TT, CAL62, BHT101, IHH4, FTC133/238, 8505C, 8305C. CRISPR n=5/6, PRISM n=6/7, GDSC2 n=6/6.

**MAPK class table (PRISM 24Q2 + GDSC2 + CTRP).**

| Drug | Class | Source | n | Cohen's d | Spearman r |
|---|---|---|---|---|---|
| dabrafenib | BRAFi | GDSC2 | 12 | **−0.57** | −0.21 |
| encorafenib | BRAFi | PRISM | 13 | **−0.60** | −0.03 |
| regorafenib | multi-RAFi | CTRP | 9 | **−1.31** (p=0.09) | −0.68 |
| RAF709 | pan-RAFi | PRISM | 13 | −0.10 | −0.18 |
| AZD-0364 | MEKi | PRISM | 13 | −0.09 | −0.15 |
| trametinib | MEKi | GDSC2 | 12 | +0.29 | −0.01 |
| LY3214996 | ERKi | PRISM | 13 | +0.23 | +0.09 |
| SCH772984 | ERKi | GDSC2 | 12 | +0.24 | +0.24 |

(Cohen's d < 0 → DM1-high more sensitive. Negative LFC/AUC = killed.)

**CRISPR (n=11, no gene FDR<0.25 by n).** MAPK/lineage hits direction-consistent with pancan: **BRAF d=−0.58, MAP2K1 d=−0.46, TSHR d=−1.09, FOXE1 d=−1.08** (all DM1-high more dependent). **MYC/KRAS/NAMPT flip sign inside thyroid (d≈+0.5)** — opposite to pancan Phase D — because thyroid lines are uniformly MAPK-driven so cross-lineage secondary genes don't stratify within-lineage.

**Cross-reference with `dm1_round6` pancan PRISM top-15.** Partial, attenuated. AZD-0364, RAF709, AZ-628, ML786, encorafenib retain correct sign but lose significance. Pancan signal was cross-lineage variance: AZD-0364 median LFC −1.875 in thyroid vs −2.082 in 886 non-thyroid lines. **Within thyroid, MAPK inhibition is uniform; only BRAFi (V600E-targeted) sub-class keeps DM1 selectivity.**

**Translation for Paper 1+2.** Within-lineage analysis splits cleanly: (1) thyroid-vs-other → MAPK pathway universally addictive (trametinib effective regardless of DM1); (2) within thyroid → BRAF V600E-targeted agents (dabrafenib, encorafenib) preferentially favor the HT-rich/FA-low DM1 stratum, matching deconv v12 sub-A (MAPK-active, BRAF V600E-enriched) as principal DM1-high driver. Clinical-grade combos require GDSC2 combo file or RAI-refractory PDX (out of scope).

**Files.** `h9_thyroid_dm1_score.tsv`, `h9_depmap_essentiality_top.tsv`, `h9_prism_drug_ranking.tsv`, `h9_mapk_drug_focus.tsv`, `h9_gdsc_ctrp_clinical_mapk.tsv`, `h9_summary.json`, `run_h9.py`.
