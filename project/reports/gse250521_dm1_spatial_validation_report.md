# GSE250521 spatial validation of 8-gene RAI / DM1-like axis — advisor brief
**Run: 2026-05-03** · Dataset: GSE250521 (Lu et al., thyroid Visium ST, 16 slides) · Owner: Seungho Cook

## TL;DR — verdict: **MODERATE / supplementary-only**
- 8 RAI genes (TPO, DIO1, TSHR, PAX8, TG, FOXE1, NKX2-1, SLC5A5) detected in **8/8** in **16/16** slides ✓
- All-spot stage trend (PT→PTC→LPTC→ATC) for sample-mean RAI_8 / DM1-like / TDS-like: **no signal** (ρ ≈ 0, p > 0.6)
- Epithelial-enriched top-50% subset: **directionally consistent** (RAI ↓, DM1 ↑, TDS ↓) at borderline (sample-mean Spearman ρ ≈ ±0.44, p ≈ 0.09)
- After residualizing log(counts)+log(n_genes) within each sample: signal **weakens to ρ ≈ ±0.35, p = 0.18** → substantial sequencing-depth confound, but direction preserved
- **Proliferation gradient survives depth correction** (epithelial subset, ρ = +0.58, p = 0.018) → stage labels biologically coherent ✓
- DM1_like ~ TDS_like all-spot ρ = **−0.89** (internal consistency check passed)

## 1. Data download / parse
| step | status |
|---|---|
| download GSE250521_RAW.tar (1.27 GB) | ✓ `project/data/raw/GSE250521/GSE250521_RAW.tar` |
| extract | ✓ 123 files; 16 Visium + 9 paired scRNA-seq (`_sc_` files **excluded** from this analysis) |
| parse → per-sample `.raw.h5ad` | ✓ 16/16 |
| QC (in_tissue + n_genes ≥ 200) | ✓ 55,873 spots retained (97.5% of in-tissue spots) |
| 20 duplicate gene symbols disambiguated with `-N` suffix | logged per sample |

Sample metadata: `project/data/processed/GSE250521/sample_metadata.tsv`

## 2. Sample / stage mapping
Stage inferred from filename label (`N→PT`, `PTC`, `LPTC`, `ATC`); `metadata_stage_raw` and `inferred_stage` both stored. **No conflict in any of 16 samples.**

| stage | n_slides | n_spots_post_QC |
|---|---|---|
| PT (N-1..4) | 4 | 14,740 |
| PTC | 4 | 16,223 |
| LPTC | 4 | 13,634 |
| ATC | 4 | 11,276 |

## 3. 8-gene availability
All 8 RAI genes present in all 16 slides. NKX2-1 / SLC5A5 alias resolution succeeded directly (no need for fallback). Per-set availability identical across samples:

| set | found / total |
|---|---|
| RAI_8 | 8 / 8 |
| TDS_like | 6 / 6 |
| CAF_ECM | 10 / 10 |
| EMT | 10 / 10 |
| Hypoxia | 9 / 9 |
| Proliferation | 7 / 7 |
| Epithelial | 5 / 5 |

## 4. Core stage trends (sample-mean Spearman, n = 16 slides)
### Raw (no depth correction)

| score | all spots ρ | p | epi top 50% ρ | p |
|---|---|---|---|---|
| RAI_8 | 0.000 | 1.00 | **−0.437** | **0.091** |
| DM1_like | 0.000 | 1.00 | **+0.437** | **0.091** |
| TDS_like | −0.109 | 0.69 | **−0.461** | **0.072** |
| CAF_ECM | −0.218 | 0.42 | −0.206 | 0.44 |
| EMT | +0.291 | 0.27 | +0.279 | 0.30 |
| Hypoxia | +0.085 | 0.75 | +0.243 | 0.37 |
| Proliferation | −0.243 | 0.37 | **+0.800** | **2.0e-4** |

### Depth-corrected (residualize log_counts + log_ngenes within sample)

| score | all spots ρ | p | epi top 50% ρ | p |
|---|---|---|---|---|
| RAI_8 | +0.170 | 0.53 | −0.352 | 0.18 |
| DM1_like | −0.170 | 0.53 | +0.352 | 0.18 |
| TDS_like | +0.243 | 0.37 | −0.327 | 0.22 |
| CAF_ECM | −0.109 | 0.69 | **−0.521** | **0.038** |
| Proliferation | −0.024 | 0.93 | **+0.582** | **0.018** |

**Read:** the raw borderline epithelial RAI/DM1 trend (ρ ≈ ±0.44, p ≈ 0.09) loses ~20% magnitude after depth correction (ρ ≈ ±0.35, p = 0.18). **Direction is preserved**, but at n = 16 the effect is underpowered. Proliferation up-trend is the only RAI-relevant biology that fully survives correction. CAF_ECM going **down** with stage in epithelial spots after depth correction is a sign that the depth control is partially over-correcting genuine tumor stromal infiltration — interpret with care.

## 5. Correlations (technical + biology checks)
- DM1_like ~ TDS_like, all spots: ρ = −0.89, p ≈ 0; sample-mean ρ = −0.87 → **internal consistency confirmed**
- DM1_like ~ log(total_counts), all spots: ρ = −0.45 → **substantial depth confound**
- DM1_like ~ Hypoxia, all spots ρ = −0.26 (negative) but **sample-mean ρ = +0.71, p = 0.002** → Simpson's-paradox-like sign flip; aggregate-level alignment with hypoxia is real, spot-level overwhelmed by within-sample technical variance
- DM1_like ~ CAF_ECM, all spots ρ = +0.24 (weak positive); sample-mean +0.10 (NS)
- DM1_like ~ EMT, weak; DM1_like ~ Proliferation, sample-mean +0.00 (NS at sample level despite strong stage Proliferation trend) — Proliferation ≠ DM1 axis here

## 6. Representative figures
- `project/results/figures_for_advisor/fig_2x4_RAI_DM1_per_stage.png` — 2 rows (RAI_8 / DM1_like) × 4 cols (PT/PTC/LPTC/ATC, one rep slide per stage), spatial heatmap
- `project/results/figures_for_advisor/fig_sample_means_RAI_DM1_TDS.png` — 16-slide sample-mean strip plot per stage
- `project/results/02_stage_trend/violin_primary_scores.png` — all-spot violin
- `project/results/02_stage_trend/violin_primary_scores_epi.png` — epithelial top-50% violin
- `project/results/01_spatial_score/per_sample/*.png` — per-slide RAI / DM1 / TDS heatmaps (16 slides)
- `project/results/02_stage_trend/stage_trend_summary.csv` — full numeric table
- `project/results/02_stage_trend/correlation_summary.csv` — DM1 vs (TDS / counts / CAF / EMT / Hypoxia / Proliferation)
- `project/results/02_stage_trend/depth_corrected_trend.csv` — depth-residualized re-test

## 7. Technical caveats
1. **Sequencing depth confound is large** — DM1_like ~ log(counts) ρ = −0.45. Any DM1 claim must be reported alongside depth-corrected version.
2. **n = 16 slides → underpowered for moderate effects.** A true ρ ≈ −0.4 needs ~30 slides for p < 0.05.
3. **Mixed-effects model misspecified here**: stage is sample-level, so `(1|sample_id)` random intercept absorbs all stage variance → all p values from `mixedlm` ≈ 1. Stage trend should be done at sample-mean level (as primary test), not spot-level mixed model. Reported in summary table for transparency, not used for inference.
4. **Epithelial subset uses simple top-50% Epithelial_score within sample** — a crude proxy. A real tumor-area mask would need pathology annotation. TG/PAX8 alone NOT used as filter (per spec — DM1 spots can have low thyroid lineage).
5. **Within-sample z-score** is the correct normalization (cross-sample batch is large), but loses cross-stage absolute comparison. The sample-mean Spearman test was chosen specifically for stage trend over absolute scoring.
6. **Hires image used as H&E** — Visium hires PNG is downsampled from full-res scan. Tile resolution is ~110 µm at default scalefactor; sufficient for tissue-architecture POC, **insufficient for nuclear morphometry**.

## 8. Go / no-go
**Verdict: MODERATE — supplementary external validation only.**
- The 8-gene RAI / DM1-like axis is **directionally consistent** with PT→ATC dedifferentiation in the epithelial-enriched subset (RAI ↓, DM1 ↑, TDS ↓), but **not statistically significant after depth correction at n = 16** (p = 0.18).
- This dataset alone **cannot be a main-figure claim** for Paper 1.
- However, it **does not contradict** the 8-gene framework — direction matches, internal consistency strong (DM1 ~ TDS ρ = −0.89), proliferation gradient confirms stage labels.

## 9. Next actions (gated)
- **If aiming Paper 1 reach venue (Cell Rep Med / JCI Insight)**: add **HRA003537** (Chinese Visium thyroid, ~30 slides) and **GSE230424** (additional thyroid Visium) to the meta-trend test. Combined n ≈ 60 slides should bring the epithelial DM1 trend to p < 0.05 if the effect is real.
- **Pathology POC** is **GO**: tile metadata ready (3,200 tiles balanced 800/stage at `project/results/03_pathology_poc/tile_metadata.tsv.gz`); spot-level DM1 ~ TDS correlation strong (ρ = −0.89) → ResNet50 baseline should be feasible. Skeleton at `project/src/05_pathology_poc/README.md`. **Do not run on this single-dataset POC alone** — wait for HRA003537/GSE230424 ingestion to give honest LOSO across **multiple cohorts**.
- **Drop scope** if: external dataset ingestion is blocked AND Paper 1 reach is downgraded to Sci Rep base. Then this analysis lives as a supplementary figure showing "directionally consistent but underpowered single-cohort spatial signal".

## 10. Paper isolation guard (this run)
- Stayed within Paper 1 scope (8-gene driver-excluded sub-stratifier). No HT/TLS/BCR/AICDA (Paper 2). No GD/HLA/Graves' (Paper 3). Immune scores intentionally NOT computed in this run to avoid Paper 2 contamination per spec.

## 11. Reproducibility
```bash
bash project/src/00_download/download_gse250521.sh
python3 project/src/01_parse_visium/parse_gse250521.py
python3 project/src/02_qc_score/score_8gene_dm1.py
python3 project/src/03_spatial_plots/plot_spatial_scores.py
python3 project/src/04_stage_trend/stage_trend_analysis.py
python3 project/src/04_stage_trend/depth_corrected_trend.py
python3 project/src/05_pathology_poc/extract_tiles.py --max-per-sample 200
```
Total wall time end-to-end: ~5 minutes on this machine (no GPU).
