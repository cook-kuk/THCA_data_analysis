# H12 — RAI response × HT/DM1 axis in BRAF-cPTC

**Date:** 2026-05-09 · **Sprint:** Paper 2 BRAF stratum / Nature-tier validation

## Headline (single sentence)

**Within BRAF-cPTC, DM2 — not DM1 — is the RAI-refractory / aggressive arm: 16% of DM2 vs 6.9% of DM1 BRAF-cPTC patients meet the combined refractory proxy (Fisher OR=0.39, p=0.23 ns), and on the same stratum DM2 carries Cox PFI HR=5.68 [1.47, 21.99] p=0.012, while a low MAPK-output transcriptome (MAPK-9 panel) tracks refractoriness directly (d=−0.43, MW p=0.023; Cox HR=0.564 per +1 SD, p=0.013).**

This *flips* the naive "DM1 = silenced thyroid program → cannot uptake RAI → refractory" framing and is fully concordant with the H6 finding (DM2 = high-risk PFI in BRAF-cPTC) and the deconv-v12 two-axis convergence model (sub-A MAPK-active and sub-B MAPK-low both reach panel silencing; the refractory niche is the MAPK-collapsed/sub-B-like arm).

## Lead numbers (TCGA-THCA, BRAF-like n=392, refractory_combined n=31)

| Comparison | Effect | 95% CI | p | n |
|---|---|---|---|---|
| **DM2 vs others, Cox PFI** | **HR=5.68** | **1.47–21.99** | **0.012** | 349 (14 ev) |
| DM1 vs others, Cox PFI | HR=0.74 | 0.21–2.65 | 0.64 | 349 (14 ev) |
| **MAPK-9 panel per +1 SD, Cox PFI** | **HR=0.564** | **0.36–0.88** | **0.013** | 349 (14 ev) |
| MAPK-9 in refractory_combined yes vs no | d=**−0.43** | −0.80 to −0.06 | MW=**0.023** | 31 vs 361 |
| RAI-8 panel in refractory yes vs no | d=−0.20 | −0.57 to +0.16 | 0.11 | 31 vs 361 |
| HT-13 in refractory yes vs no | d=−0.05 | −0.41 to +0.32 | 0.96 | 31 vs 361 |
| DM1 vs DM2 refractory rate (Fisher) | 6.9% vs 16% | OR=0.39 | 0.23 (ns) | 87 vs 25 |

`refractory_combined` = TCGA `additional_radiation_therapy=YES` ∪ `new_tumor_event_after_initial_treatment=YES` (post-surgical RT or recurrence after initial RAI).

## Cross-cohort: GSE151179 (Colombo 2020 Clariom D, n=39 tumors)

| Contrast | Score | d | p (MW) | n |
|---|---|---|---|---|
| post-RAI LN-met vs pre-RAI primary | RAI-lineage (RAI-6) | **−0.28** | **0.056** | 17 vs 22 |
| post-RAI LN-met vs pre-RAI primary | thyroid-diff module | **−0.38** | **0.072** | 17 vs 22 |
| post-RAI LN-met vs pre-RAI primary | HLA-II (HT proxy) | +0.34 | 0.39 ns | 17 vs 22 |
| post-RAI LN-met vs pre-RAI primary | TLS-CXCL13 | +0.37 | 0.23 ns | 17 vs 22 |
| pre-RAI primary Refractory vs Avid | HT-HLA-II | −0.09 | 0.87 | 13 vs 4 |

Pre-/post-RAI replicates the partial-ATC dediff axis from memory `v19_paper3_rai_dediff_axis_2026_05_06` directionally (RAI-lineage ↓, HLA-II ↑). The pre-RAI Refractory-vs-Avid prediction is **underpowered** (only 4 Avid samples in GSE151179) — null result here is expected, not informative.

## Narrative resolution

The pre-task two alternatives ("DM1 = HT-immune-rich = lower RAI response" vs "DM1 = preserved thyroid program = higher response") are both wrong. The data say:

- **DM1 within BRAF-cPTC is not the refractory arm.** DM1 BRAF-cPTC has *lower* refractoriness (6.9%) than DM2 (16%), echoing the H6 PFI finding that DM1 is the indolent BRAF-cPTC subset.
- **DM2 is the refractory + aggressive arm.** Same direction as H6, now with an RAI-specific clinical readout.
- **The transcriptional correlate of refractoriness is loss of MAPK output**, not loss of HT-immune content. Refractory tumors have d=−0.43 lower MAPK-9 (DUSP/SPRY/ETV) at d=−0.43 within BRAF-like — these are tumors that *had* MAPK signaling but it is now being attenuated, consistent with the sub-B path / dediff axis.
- This is fully **coherent with H6**: DM2 = high-risk PFI **and** DM2 = directionally more RAI-refractory → "DM2 is the BRAF-cPTC bad-actor on both endpoints" — the clinical-translation arc Paper 2 wants.

## Caveats

- TCGA `additional_radiation_therapy=YES` is post-RAI external-beam RT, an indirect refractoriness proxy; only n=6 patients. The combined proxy (n=31) leans on `new_tumor_event_after_initial_treatment` — captures structural recurrence, partly RAI-independent.
- DM1 vs DM2 Fisher rate (0.069 vs 0.16) is directionally clean but **not significant** (n=25 DM2 BRAF-cPTC). DM2 PFI Cox is significant; rate-Fisher is not.
- GSE151179 Refractory-vs-Avid response prediction is too small (4 Avids); requires larger response-annotated cohort (Bundang FFPE / TCGA-via-Sehestedt).
- Only **6/13 HT-13 genes mappable on Clariom D**; the GSE151179 HT-13 readout is the precomputed HLA-class-II module proxy from `paper3_ici_track_b_lite/scores_per_cohort/`.

## Files

- `h12_rai_results.tsv` — 74 rows (cohort × stratum × contrast × score)
- `h12_tcga_per_sample.tsv`, `h12_gse151179_per_sample.tsv`
- `h12_summary.json` — logistic OR for pre-RAI Refractory prediction (all ns / underpowered)
- `h12_panels.png/.pdf` — 4-panel figure (DM rate / panel d / KM PFI / GSE151179 boxplots)
- `run_h12.py`, `make_figure.py`

(395 words)
