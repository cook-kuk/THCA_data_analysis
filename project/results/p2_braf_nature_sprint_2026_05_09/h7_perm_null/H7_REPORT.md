# H7 image rescue via attention map analysis (2026-05-09)

**Cohort.** N=41 BRAF_like/cPTC TCGA-THCA (DM1=26, DM2=15). Same 5-fold splits as phase2 v2 CLAM run; image-only pooled AUC=0.592 is the failure we are diagnosing.

## (1) Attention map vs HT-13 RNA
Per slide, the 200x1024 ImageNet ViT-L tiles were forwarded through the held-out CLAM checkpoint to extract the gated-attention softmax over tiles. Across 8 concentration metrics (entropy, top-k mass, max, effective N), **none correlated significantly with HT-13 RNA**. Strongest: `attn_top1_mass` Spearman rho=-0.262 p=0.097. DM1 vs DM2 entropy ratio: no metric p<0.05 (best p=0.30). Interpretation: the CLAM gated attention is essentially decoupled from the lymphoid axis separating DM1 from DM2.

## (2) Tile-level linear probe to HT-13
Per fold we fit a Ridge probe (alpha=10) on ~6,500 training tiles (1,024-d) with slide-level HT-13 z broadcast as per-tile target, then scored held-out tiles. Per-slide statistics (max, p90, attention-weighted mean) entered as new image features.

| Setting | pooled AUC | 95% CI | delta vs HT13 |
|---|---:|---|---:|
| ht13_only | 0.887 | [0.764, 0.982] | +0.000 |
| ht13_plus_tile_probe_max | 0.877 | [0.744, 0.974] | -0.010 |
| ht13_plus_tile_probe_attn_w | 0.877 | [0.754, 0.974] | -0.010 |
| ht13_plus_img_prob | 0.869 | [0.749, 0.967] | -0.018 |
| ht13_plus_tile_probe_p90 | 0.867 | [0.744, 0.969] | -0.021 |
| ht13_plus_tile_probe_all | 0.808 | [0.656, 0.931] | -0.079 |
| img_prob_only | 0.564 | [0.377, 0.738] | -0.323 |
| tile_probe_max_only | 0.531 | [0.349, 0.721] | -0.356 |

Univariate tile-probe stats AUCs sit at 0.42–0.55 (chance). Probe-stat vs HT-13 Spearman is at noise (|rho|<0.06, p>0.7).

## Verdict
**RESCUE-FAIL (neutral).** Every multimodal combo lands inside the HT-13-only bootstrap CI; the best image add-on is delta=-0.010. Univariate H&E features (CLAM prob_DM1, tile-probe max/p90/attn-w) are at chance. With N=41 and generic ImageNet ViT-L tiles, no linear aggregation recovers lymphoid signal beyond what HT-13 RNA already captures. Honest claim becomes "molecular DM1 stratification cannot be derived from histology alone in BRAF-cPTC; multimodal RNA+H&E shows no uplift over RNA alone."

## Caveats
- N=41, ~8/fold; per-fold AUCs are noisy. Headline = pooled OOF + bootstrap CI.
- Falsifies "rescue with current features" only. A thyroid- / TLS-pretrained foundation model (UNI, Virchow, Phikon) could still rescue and remains an open branch.
- Probe target is intentionally noisy (slide-level label broadcast to 200 tiles); max/p90 were chosen to surface rare lymphoid tiles. Negative result bounds what a generic linear probe can extract from this feature space.

## Files
- `h7_attention_metrics.tsv`, `h7_attention_corr_HT13.tsv`, `h7_attention_DM1vsDM2.tsv`
- `h7_tile_probe_diagnostics.tsv`, `h7_tile_probe_results.tsv`
- `fig_h7_attention_rescue.png/.pdf` 3-panel diagnostic
