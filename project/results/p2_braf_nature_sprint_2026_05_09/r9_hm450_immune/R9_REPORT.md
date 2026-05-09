# R9 — Pan-immune methylation deep-dive (BRAF-cPTC)

**Headline.** DM1 BRAF-cPTC does **NOT** show measurable immune-promoter de-methylation beyond the 8-gene thyroid-differentiation panel. The 8-gene panel **IS** capturing the dominant methylation layer in this stratum; the immune RNA gradient is collinear with 8-gene mean β and vanishes when 8-gene β is regressed out.

## Data audit
- Full TCGA-THCA HM450 matrix: **NOT on disk**. GPL13534 manifest: **NOT on disk**.
- 8-gene HM450 panel (n=503): on disk (`r5_2_sample_methylation_8gene.tsv`).
- GSE97466: only top-5000 variable probes; no immune-promoter mapping possible offline.
- Pivot per task §6: RNA-z proxy for promoter activity (caveat: bulk RNA cannot resolve promoter vs enhancer). 8-gene HM450 used as positive control.

## Stratum
BRAF-cPTC, TCGA-THCA: total n=363; DM1+DM2 contrast n=112 (DM1=87, DM2=25); all 112 have both 8-gene HM450 and RNA.

## Findings

**(A) 8-gene HM450 (positive control).** `mean_8g_beta` d=+2.18, p=1.3e-11. TPO d=+2.99, DIO1 +1.22, TG +1.18, TSHR +1.00, PAX8 +0.85, FOXE1 +0.74, NKX2-1 +0.57; SLC5A5 d=−0.22 ns. DM1 hyper-methylated.

**(B) Immune RNA proxy (32 genes, 7 modules).** 31/32 FDR<0.05; **all 32 d>0 (UP in DM1)**. Top hits: HLA-DRA d=+2.06, HLA-DRB1 +2.04, CIITA +1.85, HLA-B +1.84, HLA-DPB1 +1.83, HAVCR2 +1.82, HLA-C +1.65, CD274 +1.50, HLA-A +1.46, CCL5 +1.38. Module composites: ag_pres_II +2.02, ag_pres_I +1.63, treg/supp +1.17 (FOXP3 *up*, not down), checkpoint +1.21, effector +1.00, chemokine +1.09, b_cell +0.77.

**(C–D) Independence test — decisive.** `immune_demeth_proxy_up` raw d=+1.37 (p=2.6e-08) → **after residualizing 8-gene mean β: d=−0.017, p=0.80**. All 7 module proxies collapse to |d|<0.5, p>0.08 (b_cell raw 0.01 → ns FDR). The 8-gene HM450 axis explains the immune-RNA gradient.

**(E) RNA↔HM450 concordance.** All 8 8-gene Spearman r<0 (median −0.38; TPO −0.68). RNA-z direction is the right proxy semantic.

**(F) Causal layering (rai_score).** R² M0 (MAPK only)=0.956 → +8-gene meth Δ=2e-4 → +immune-proxy Δ=3e-3. HT axis is MAPK-determined; methylation layers add <0.5% R².

## Implication
The 8-gene panel is not an arbitrary subset — it captures the methylation layer that separates DM1 from DM2 within BRAF-cPTC. With disk data alone we cannot demonstrate an *additional* immune-promoter de-methylation layer. Fully closing the question requires the full HM450 + GPL13534 manifest (deferred — flag for re-execution when network/manifest is available).

Files: `r9_per_gene_methylation_d.tsv`, `r9_8gene_HM450_DM1_vs_DM2.tsv`, `r9_immune_promoter_composite.tsv`, `r9_methylation_rna_concordance.tsv`, `r9_independence_from_8gene_meth.tsv`, `r9_composite_DM1_vs_DM2.tsv`, `r9_causal_layering.tsv`, `r9_summary.json`, `run_r9.py`.
