# Paper 2 — BRAF_like/cPTC multimodal rescue (2026-05-09)

**Cohort.** N=41 BRAF_like / cPTC TCGA-THCA primary tumors with all three modalities available (CLAM image embedding + bulk-RNA z-score + HM450 8-gene beta). Class balance: 26 DM1 / 15 DM2.

**Baseline.** Image-only CLAM (ImageNet ViT-L tile features → gated-attention MIL) pooled OOF AUC = 0.592 on this stratum, vs the overall N=59 mean AUC 0.83 ± 0.14 reported in `phase2_tcga_clam/PHASE2_REPORT.md`. The RAS_like/FVPTC stratum (n=16) had AUC=1.00; BRAF_like/cPTC is the failure mode we are rescuing.

**Label-leak audit (DECISIVE for honest reporting).** DM1/DM2 was DEFINED upstream from the HM450 8-gene panel (DIO1, FOXE1, NKX2-1, PAX8, SLC5A5, TG, TPO, TSHR; see `audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv`). On the cohort: DM1 mean_8g_β=0.384 vs DM2=0.258 — near-perfect by construction. Therefore: (a) HM450 betas of the 8-gene panel **trivially** reproduce the label, (b) RNA z-scores of the SAME 8 genes are ~80% redundant with the label via the well-documented methylation↔expression inverse correlation. We split the RNA panel into a NON-LEAK block (9 MAPK-output + 13 HT/B-cell = 22 genes) and a LEAK block (8 thyroid-diff genes), and likewise mark methylation as TAUTOLOGY. Headline conclusions use NON-LEAK combos only.

**Setup.** Same 5 folds as the image-only run (loaded from `clam_per_slide_predictions.tsv`). Per fold we (a) re-load the matching CLAM checkpoint, (b) forward the 200×1024 tile features through gated-attention pooling to obtain the 512-d slide bag vector, (c) concatenate with the chosen modality blocks, (d) standardize each block on the training fold, (e) fit L2-regularized LogReg (C=0.5). Held-out fold predictions are pooled across all 5 folds for one OOF AUC. 95% CI is a 1,000-resample stratified bootstrap.

## Headline (NON-LEAK combos)

| Combo | n_feat | pooled AUC | 95% CI | per-fold mean ± std |
|---|---:|---:|---|---|
| img_prob_plus_rna_nonleak | 23 | 0.995 | [0.977, 1.000] | 1.000 ± 0.000 |
| rna_nonleak_only | 22 | 0.992 | [0.969, 1.000] | 1.000 ± 0.000 |
| img_prob_plus_rna_mapk | 10 | 0.977 | [0.931, 1.000] | 0.960 ± 0.080 |
| rna_mapk_only | 9 | 0.974 | [0.926, 1.000] | 0.960 ± 0.080 |
| rna_ht_only | 13 | 0.926 | [0.790, 1.000] | 0.956 ± 0.089 |
| img_prob_plus_rna_ht | 14 | 0.918 | [0.790, 0.997] | 0.936 ± 0.088 |
| image_emb_plus_rna_nonleak | 534 | 0.738 | [0.569, 0.877] | 0.960 ± 0.049 |
| image_only_loaded | 1 | 0.564 | [0.377, 0.738] | 0.743 ± 0.216 |
| image_emb_only | 512 | 0.451 | [0.274, 0.644] | 0.436 ± 0.207 |

**Best non-leak.** `img_prob_plus_rna_nonleak` AUC=0.995 (95% CI [0.977, 1.000]); Δ vs image-only = +0.403.

## Reference (label-leaking; included for transparency only)

| Combo | n_feat | pooled AUC | 95% CI | per-fold mean ± std |
|---|---:|---:|---|---|
| rna_leak_plus_meth_TAUTOLOGY | 17 | 0.987 | [0.956, 1.000] | 1.000 ± 0.000 |
| rna_leak_only_TAUTOLOGY | 8 | 0.959 | [0.882, 1.000] | 0.960 ± 0.049 |
| meth_only_TAUTOLOGY | 9 | 0.918 | [0.803, 0.992] | 0.937 ± 0.064 |
| all_three_LEAK | 551 | 0.905 | [0.797, 0.982] | 1.000 ± 0.000 |
| image_emb_plus_meth_LEAK | 521 | 0.621 | [0.444, 0.782] | 0.647 ± 0.238 |

These rows are tautologies (or near-tautologies) because the modality being tested IS, or is tightly coupled to, the modality used to define DM1/DM2. They are NOT a generalization claim — they are a sanity check that the data join and pipeline work end-to-end.

## Honesty caveats

- N=41 with held-out folds of ~8 slides — single-fold AUCs are very noisy; treat the bootstrap-CI on the pooled OOF AUC as the headline uncertainty.
- The CLAM checkpoint for fold k was trained without that fold's slides, so per-slide image embeddings and `img_prob_DM1_loaded` are valid OOF (matches the original v2 0.830 ± 0.139 mean AUC).
- `image_only_loaded` runs the held-out scalar prob_DM1 through the same per-fold StandardScaler + LogReg as the other combos — this can monotonically remap and thus give a slightly different pooled AUC from the raw 0.592. This is BY DESIGN (apples-to-apples) but means the reported `image_only_loaded` AUC is NOT the same number as the published 0.592; see the script header for the rationale.
- RNA z-scores are pan-TCGA-cohort-relative. Within-fold leakage is impossible (z-scores are fixed at the cohort level), but a Korean K2 validation will need to rebuild z-scores within the K2 cohort before transferring this classifier.
- All combos use identical CV / standardization / hyper-parameters. No per-combo tuning, no model selection on test folds.
- 512-d image embedding + 22 RNA = 534 features for N=41 is high-dim/low-N; the bootstrap CI is the relevant uncertainty, not the point AUC.
- The 9 MAPK-output genes (DUSP4/5/6, SPRY2/4, ETV4/5, PHLDA1, CCND1) are a downstream readout of BRAF V600E activation. In the BRAF_like/cPTC stratum these should NOT separate DM1 vs DM2 by driver alone — both groups are BRAF-like. The HT/B-cell genes (HLA-DRA/DRB1/DPA1/DPB1/DQA1/DQB1, CD79A/B, MS4A1, AICDA, CXCL13, CCR6, IFNG) capture the Hashimoto-overlap immune axis (per v17_D5P6_BCR_clonal_TLS) which IS expected to enrich in DM1.

## Files
- `braf_multimodal_results.tsv` — full table
- `braf_multimodal_per_slide_preds.tsv` — OOF prob per slide per combo
- `fig_braf_multimodal_uplift.png` / `.pdf` — bar+CI figure (gray = leaking)