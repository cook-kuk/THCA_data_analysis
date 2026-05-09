# Paper 2 BRAF stratum H2 — fatty-acid + cytotoxic + TLS axes

**Headline — TWO-axis convergence confirmed.** Within BRAF_like/cPTC TCGA-THCA (n=41, DM1=26 / DM2=15), the **fatty-acid panel (12 genes) reaches pooled OOF AUC=0.969 [0.918, 1.000] as a solo classifier** while being statistically independent of the immune axis (HT↔FA r=−0.156; FA↔cyto r=−0.226; FA↔TLS r=−0.158; FA↔MAPK r=+0.007). Combining the two orthogonal axes (HT+FA, 25 features) gives **AUC=0.974 [0.926, 1.000], +0.048 over HT alone**, and HT+FA+cyto+TLS pushes the ceiling to **0.982 [0.944, 1.000]**. **DM1 vs DM2 within BRAF_like is therefore separable along at least two biologically distinct axes, not just immune** — supporting "two-axis convergence" framing for Paper 2.

**Independence verdict (mean-panel Pearson r).** The three "immune" panels collapse into a single axis (HT↔TLS r=+0.92, HT↔cyto r=+0.90, cyto↔TLS r=+0.89) and should be treated as one immune compartment with three readouts. **FA is genuinely orthogonal to all three** (|r| ≤ 0.23) AND to MAPK output (r=+0.007). MAPK↔immune is also low (r≈+0.18–0.23), giving three candidate axes: immune, fatty-acid metabolism, MAPK-output.

**Effect sizes (Cohen d, DM1 − DM2 on mean panel z).** ht +1.58, mapk +1.81, cyto +1.23, tls +1.07, **fa −0.64** (lipid genes silenced in DM1 — direction matches the thyroid-de-differentiation prediction).

**MAPK-output caveat.** mapk solo AUC=0.974 (d=+1.81, DM1>DM2) was UNEXPECTED. Within BRAF_like both groups are BRAF-driven, so this should be near-chance. Two non-exclusive explanations: (a) the 9 mapk-output genes (DUSP4/5/6, SPRY2/4, ETV4/5, PHLDA1, CCND1) co-vary with thyroid-differentiation programs (TPO/PAX8 axis modulates RTK→MAPK feedback), so they may be a downstream readout of the same DM1 program rather than driver-independent, OR (b) DM1 within BRAF_like has higher BRAF-V600E pathway output than DM2. This is consistent with `deconv_v5_v12_findings_2026_05_08`: MAPK output × 8-gene panel inverse correlation. We treat MAPK as a *secondary* axis here and lead with HT+FA for the cleanest two-axis claim.

**Combo ladder.** fa 0.969 / mapk 0.974 / ht 0.926 / tls 0.918 / cyto 0.903 (singletons). HT+FA 0.974 / ht+cyto 0.936 / ht+tls 0.923 / fa+cyto+tls 0.969 / ht+fa+cyto+tls **0.982**.

**Caveats.** N=41 with ~8 slides/fold; CI is the headline uncertainty (no point AUC differs by >0.05 from its CI midpoint among the top combos, so paper claim should be the floor of CI for HT+FA = 0.926). RNA-z is pan-TCGA-cohort-relative; cross-cohort transfer needs within-cohort z-rebuild. All 35 H2-target genes present in the TCGA-THCA RNA-z file (no drops).

**Files.** `h2_results.tsv` (combo AUC) · `h2_panel_correlations.tsv` (pairwise Pearson) · `h2_panel_effect_sizes.tsv` (Cohen d) · `h2_per_slide_preds.tsv` (OOF probs) · `run_h2_fa_axis.py` (script).
