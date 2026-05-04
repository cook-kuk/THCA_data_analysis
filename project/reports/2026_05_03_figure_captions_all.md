# Figure + Supplementary Figure + Supplementary Table Captions (all)

**Date:** 2026-05-03 (marathon scaffolding/infra)
**Format:** Cell Press caption format — paste-ready for v8 manuscript figure legends section
**Voice-protected:** Caption verb tone borderline (user may polish "establishes" / "reveals" / "identifies" verb consistency). All numerical content fact-only.

> **Scope status (2026-05-04 update).** This caption file was authored 2026-05-03 for the v8 5-pillar mega-paper.
> Captions for **F1 (HLA forest)** = Paper 4; **F2 (GSE286332)**, **F5 (autoimmune-PTC 4-panel)**, **F6 (TCGA Hashimoto)**,
> **F7 (sub-A vs sub-B)** = Paper 2. Only **F3 (driver neutrality)** and **F4 (pan-genome ARI)** of the legacy F1–F7 set
> belong to Paper 1.
> Paper 1 v8 has been re-figured per `2026_05_04_8gene_curated_vs_denovo_final_strategy.md` §9 and
> `2026_05_03_manuscript_v8_OUTLINE.md` §Figures (updated 2026-05-04). The new Paper 1 figures
> are F1 (driver orthogonality) / F2 (axis robustness) / F3 (axis anchor) / F4 (mechanism + TROP2) /
> F5 (outcome) / F6-Supp (spatial). Captions for these new Paper 1 figures will be drafted in a separate
> caption pass; they are NOT in this file yet. Existing F1–F7 captions below are retained for Paper 2 / Paper 4
> future use; do not delete.

---

## Main Figures (F1-F7) — `[LEGACY 5-pillar v8; carved per Paper 1 / 2 / 4 split, see banner above]`

### Figure 1 — Pan-Asian HLA forest meta (Pillar 1)

```
Figure 1. Pan-Asian autoimmune-thyroid HLA susceptibility allele continuum 
across Korean papillary thyroid carcinoma (PTC) and Han Chinese Graves' 
disease (GD).

(A) Three-cohort allele frequency comparison for 5 focus HLA alleles 
(DPB1*05:01, B*46:01, C*01:02, A*02:07, DRB1*07:01). Korean PTC pool 
(n=874; K2 PRJEB11591 n=235 + Lee 2024 GSE213647 n=630 + GSE286332-PTC 
n=9; red bars) compared with Han Chinese controls (n=1,490; gray) and 
Han Chinese GD cases (n=1,468; navy) from Chu et al. 2018 J Med Genet 
55(10):685-692. Error bars show 95% Wilson confidence intervals.

(B) Random-effects DerSimonian-Laird forest plot pooling Korean-PTC-vs-
Chu-control + Chu-GD-vs-control estimates per allele. Diamonds represent 
pooled odds ratios (OR) with 95% CI. Vertical dashed line indicates 
OR=1. Heterogeneity I² (range 0-91%) shown to the right of each row.

DPB1*05:01 is the most stable risk allele across cohorts: Korean PTC 
53.2% (95% CI 49.9-56.5%) > Han Chinese GD 44.0% > Han Chinese controls 
31.3%; pooled OR 2.16 (95% CI 1.65-2.83), I²=85%, with within-Korean 
sub-cohort homogeneity I²=0% (K2 56%, Lee 52%, GSE286332-PTC 56%).

Source data: results/p2_pillar1_forest/{forest_meta_results.tsv, 
random_effects_pooled.tsv, korean_subcohort_heterogeneity.tsv}.
```

### Figure 2 — GSE286332 PTC vs PTC+HT molecular dissection (Pillar 2)

```
Figure 2. Differential expression and pathway enrichment in PTC versus 
PTC+Hashimoto's thyroiditis (GSE286332, n=9 vs n=9, Korean cohort).

(A) Volcano plot of 29,672 genes tested by PyDESeq2 Wald test with 
Benjamini-Hochberg FDR correction. Red points (significant up): 
6,004 DEGs (padj<0.05, log2FC>0); navy points (significant down): 
4,376 DEGs. Top 8 up-regulated and top 5 down-regulated genes labeled. 
Total significant DEGs at padj<0.05: 10,380.

(B) MSigDB Hallmark gene set enrichment analysis (GSEA pre-ranked, 
1,000 permutations). Top 8 positively enriched (red bars) and top 4 
negatively enriched (navy bars) gene sets sorted by NES; FDR q-values 
shown. Hallmark Allograft Rejection NES=+2.12, Interferon Gamma 
Response NES=+1.80 (FDR=2e-4), Inflammatory Response (FDR=1.3e-4).

(C) Per-gene boxplots of 8-gene RAI panel members (SLC5A5, TPO, TG, 
TSHR, PAX8, NKX2-1, FOXE1, DIO1) in PTC (navy, n=9) vs PTC+HT (red, 
n=9). Pooled 8-gene Cohen's d = -1.60, p=0.008. PAX8 d=-2.32, NKX2-1 
d=-1.92, FOXE1 d=-1.75; SLC5A5 preserved (d=+0.25), indicating 
TF-driven dedifferentiation rather than transporter loss.

(D) HLA-II module heatmap (10 genes including HLA-DRA/B1, DPA1/B1, 
DQA1/B1, DMA/B, CIITA, HLA-DOB) per sample (z-score across cohort). 
Group separator (vertical line) at sample 9-10. HLA-II Cohen's 
d=+3.65, MW p=4e-4 (★ exceptional in this cohort).

Source data: results/p3_gse286332/{deg_ptcht_vs_ptc.tsv, 
gsea_MSigDB_Hallmark_2020.tsv, 8gene_panel_per_sample.tsv, 
hla_module_scores.tsv}.
```

### Figure 3 — Driver mRNA neutrality (Pillar 3)

```
Figure 3. BRAF/RAS/TERT driver mRNA expression is independent of 
mutation status; driver_anchor genes alone cannot define DM clusters.

(A) BRAF mRNA expression (log2 expression) by V600E mutation status 
in TCGA-THCA (V600E carriers n=273, wild-type n=182). Cohen's d = 
-0.044 (negligible), Mann-Whitney p = 0.567. Driver transcript ≠ 
driver mutation.

(B) Driver mRNA single-feature AUC for DM1 vs DM2 classification. 
BRAF 0.602, TERT 0.578, KRAS 0.525, NRAS 0.521, HRAS 0.500. 
Horizontal lines: AUC=0.5 (chance) and AUC=0.7 (threshold). All driver 
AUCs <0.7 confirm that driver transcripts cannot define DM clusters.

(C) TIERA67 67-gene Cohen's d ranking (top 24). 8-gene panel members 
(red, ranks 4-50; TPO #4, DIO1 #5, TG #18, PAX8 #19) cluster at top 
ranks; Driver_anchor genes (navy, ranks 15-67; CDKN2B #15, RET #22, 
CDKN2A #24, BRAF #52, TERT #56, KRAS #63, NRAS #66, HRAS #67) cluster 
at bottom ranks.

(D) Adjusted Rand Index (ARI) of cluster panels vs DM1/DM2 labels. 
Driver_anchor 12 alone yields ARI=-0.007 (random); 8-gene panel ARI=
0.489 (modest, clinically interpretable trade-off); TIERA67 ARI=0.903 
≈ pan-genome top-5000 MAD ARI=0.918 (cluster definition robust).

Source data: results/p1_driver_mrna_audit/{driver_mrna_mutation_audit.tsv, 
driver_mrna_dm_auc.tsv, top20_by_d_full_tiera67.tsv}.
```

### Figure 4 — Pan-genome cluster robustness (Pillar 4)

```
Figure 4. The DM1/DM2 axis is reproducible by pan-genome unrestricted 
gene selection.

(A) ARI vs DM1/DM2 cluster across 6 panel choices. Driver_anchor 12 
alone (-0.007) confirms drivers cannot define DM. 8-gene panel alone 
(0.489) modest. TIERA67 (0.903) ≈ pan-genome top-200 MAD (0.864) ≈ 
top-1000 (0.902) ≈ top-5000 (0.918) — cluster definition independent 
of candidate-pool restriction.

(B) Top-N coverage ladder showing percentage of TIERA67 (green) and 
8-gene panel (red) members in pan-genome top-N univariate Cohen's d 
ranks. Hypergeometric enrichment of TIERA67 in pan-genome top 100 = 
3 × 10⁻⁴ (★), justifying biological prior.

Source data: results/p4_pangenome_vs_tiera67/{ari_comparison.tsv, 
topN_coverage_ladder.tsv}.
```

### Figure 5 ★ — Autoimmune-PTC mechanism layer (Pillar 5)

```
Figure 5. Mechanistic dissection of the autoimmune-PTC overlap 
sub-population.

(A) GSE286332 P(DM1) probability spectrum across 18 samples (PTC n=9 
navy, PTC+HT n=9 red). Linear mixed-method mediation analysis: HLA-II 
mediates 140% (over-mediation/full-pathway, bootstrap 95% CI 
[-0.31, -0.03], p_emp=0.023) of the PTC+HT → P(DM1) ↓ relationship; 
8-gene RAI mediates 63% (p_emp=0.002).

(B) TCGA-THCA Hashimoto-like signature × DM cluster forest plot. 
Five thresholds (GMM, Otsu, top-20%, top-30%, residualized Otsu) all 
show DM2 enrichment with OR ranging 0.20-0.29 and Fisher p=6.4×10⁻¹⁰ 
to 8e-9. Hashimoto-like prevalence: 18-30% across thresholds. Vertical 
dashed line: OR=1 (no enrichment).

(C) Per-sample heatmap of BCR repertoire (IGHV total, IGHV clonality), 
TLS 12-gene signature, HLA-II/I modules, and 8-gene RAI score across 
18 GSE286332 samples (PTC n=9, PTC+HT n=9). Group separator (horizontal 
line) at sample 9-10. PTC+HT shows elevated IGHV clonality (d>0.5), 
TLS (Cabrita 2020 12-gene; Cohen d = +1.96), and AICDA (somatic 
hypermutation enzyme) consistent with antigen-driven clonal B-cell 
response.

(D) DM1 sub-A vs sub-B mutation status stacked bar (n=84 vs n=56). 
sub-A: 51 RAS+ (69% of mutation-tested; 61% of total n=84), 1 BRAF+; 
sub-B: 1 BRAF+, 2 RAS+, 53 mutation-negative (53/56, 94.6%). 
Hashimoto-like overlap (red bar): sub-B 12.5% vs sub-A 3.6% (4× higher).

Source data: results/d3p5_pdm1_gradient/{D3P5_summary.json, 
mediation_results.json}; results/d4p2_tcga_hashimoto_signature/
{D4P2_summary.json}; results/d5p6_bcr_repertoire/{per_sample_diversity.tsv, 
tls_score_per_sample.tsv}; results/d6p7_dm1_subcluster/
{dm1_subcluster_labels.tsv}.
```

### Figure 6 — TCGA Hashimoto-like distribution + cross-cohort

```
Figure 6. Cross-cohort generalization of the autoimmune-PTC signature.

(A) Density distribution of GSE286332 PTC+HT-derived signature score 
applied to TCGA-THCA (n=500, navy) and Korean GSE213647 (n=632, red) 
samples. Bimodality coefficients: TCGA 0.552, Korean 0.471 (boundary 
bimodal).

(B) Hashimoto-like signature prevalence by cohort. TCGA-THCA: GMM 18%, 
Otsu 19.6%; Korean GSE213647: GMM 22.8%, Otsu 28.2%; GSE286332 
PTC+HT only: 100%. The 22-28% Korean rate cohort-replicates the 
18-20% TCGA rate, confirming axis generalization.

Source data: results/d4p2_tcga_hashimoto_signature/, 
results/d8b_korean_replication/.
```

### Figure 7 — DM1 sub-A vs sub-B + K2 NBNR clinical phenotype

```
Figure 7. The DM1 sub-B sub-cluster represents the Korean BRAF-/RAS- 
NBNR (no-BRAF-no-RAS) population.

(A) Mutation breakdown stacked bar: sub-A (n=84) 51 RAS+ (69% of 
mutation-tested; 61% of total n=84) classical FVPTC; sub-B (n=56) 
53 mutation-negative (53/56, 94.6%).

(B) Hashimoto-like positive rate: sub-A 3.6% vs sub-B 12.5% (4× higher 
in sub-B, Fisher OR=0.26, p=0.09).

(C) Korean GSE213647 (n=632) sub-B-like signature transfer: GMM 47.2% 
/ Otsu 52.5% — approximately half of Korean PTC samples carry the 
sub-B signature, replicating Yu professor's K2 NBNR phenotype 
observation.

Source data: results/d6p7_dm1_subcluster/, 
results/d8c_dm1_subB_x_K2_NBNR/.
```

---

## Supplementary Figures (SF1-SF6)

### Suppl Figure 1 — TIERA67 pan-genome rank ladder

```
Suppl Figure 1. TIERA67 enrichment in pan-genome univariate Cohen's d 
top-N rankings. TIERA67 (green) and 8-gene panel (red) coverage as 
percentage of total TIERA67/8-gene members in top-N. Hypergeometric 
p (TIERA67 in top 100) = 3 × 10⁻⁴, justifying biological prior. Source: 
results/p4_pangenome_vs_tiera67/topN_coverage_ladder.tsv.
```

### Suppl Figure 2 — Mediation Baron-Kenny detail

```
Suppl Figure 2. Mediation analysis Baron-Kenny three-equation framework 
with 5,000-iteration bootstrap. Indirect effect (a × b) bootstrap 95% 
confidence intervals shown. HLA-II 140% mediation (★ p_emp=0.023) is 
the dominant pathway; 8-gene RAI 63% (p_emp=0.002) is parallel partial 
mediator; HLA-I 87% (p_emp=0.042); immune-proxy 88% (p_emp=0.120 NS). 
Source: results/d3p5_pdm1_gradient/mediation_results.json.
```

### Suppl Figure 3 — Korean GSE213647 replication detail

```
Suppl Figure 3. Korean GSE213647 (Lee 2024, n=632) Hashimoto-like 
signature cross-cohort replication.

(A) Score density distribution overlay (TCGA-THCA n=500 navy vs Korean 
GSE213647 n=632 red). Right-tail elevation in Korean cohort consistent 
with elevated Korean Hashimoto-like background.

(B) Threshold-stratified prevalence: GMM 18.0% TCGA vs 22.8% Korean; 
Otsu 19.6% TCGA vs 28.2% Korean. Korean rate elevated 1.3-1.4× over 
TCGA, consistent with population-level Hashimoto background but axis 
direction preserved.

Source: results/d8b_korean_replication/D8B_summary.json.
```

### Suppl Figure 4 — GSE286332 BCR repertoire detail

```
Suppl Figure 4. GSE286332 (n=18) per-sample BCR repertoire metrics.

(A) IGHV total expression (log10) by group. PTC n=9 vs PTC+HT n=9. 
PTC+HT ~5× higher (Cohen d > 1.5).

(B) IGHV clonality (1 - normalized Shannon entropy) by group. PTC+HT 
shows elevated clonality (Cohen d > 0.5), consistent with antigen-
driven clonal expansion rather than polyclonal infiltration.

Source: results/d5p6_bcr_repertoire/per_sample_diversity.tsv.
```

### Suppl Figure 5 — DM1 sub-cluster phenotype detail

```
Suppl Figure 5. DM1 sub-A (n=84) vs sub-B (n=56) phenotype detail.

(A) Mutation breakdown stacked bar: sub-A BRAF+ 1 + RAS+ 51 + neg 32; 
sub-B BRAF+ 1 + RAS+ 2 + neg 53.

(B) Hashimoto-like positive rate (Otsu threshold): sub-A 3.6% vs 
sub-B 12.5% (4× higher).

(C) Age at diagnosis distribution by sub-cluster. Sub-A vs sub-B age 
distributions overlap; minor age difference (Cohen d=0.34, MW p=0.12).

Source: results/d6p7_dm1_subcluster/, results/d4p2_tcga_hashimoto_signature/.
```

### Suppl Figure 6 — Pan-Asian sensitivity 4-scenario forest

```
Suppl Figure 6. Pan-Asian forest sensitivity analysis across 4 Korean 
PTC pool scenarios for 5 focus HLA alleles.

Korean PTC pool scenarios: (i) full pool n=874; (ii) excluding 
GSE286332-PTC n=865; (iii) Lee 2024 only n=630; (iv) K2 only n=235.

DPB1*05:01 carrier frequency: 53.2% / 53.2% / 52.1% / 56.2% (range 
±2%) — extremely robust. B*46:01 stability 10.0-10.3%. DRB1*07:01 
(protective) 9.8-11.8%. DQB1*02:01: 0% across all scenarios. 

Source: results/p2_pillar1_forest/sensitivity_4scenarios.tsv.
```

---

## Supplementary Tables (S1-S8)

### S1 — Cohort assembly

```
Suppl Table 1. Cohort assembly. Five primary cohorts and one published 
summary statistics dataset spanning 2 Korean (K2 PRJEB11591 n=235 
valid, Lee 2024 GSE213647 n=630, Lim 2025 GSE286332 n=18) + mixed-
ancestry (TCGA-THCA n=500) + Han Chinese (Chu et al. 2018 n=2,958 GD/
ctrl). Total cohort n=1,432 East Asian PTC + 2,958 Han Chinese GD/ctrl. 
DM cluster availability per cohort + HLA imputation modality + use 
case per pillar.
```

### S2 — TIERA67 candidate gene pool

```
Suppl Table 2. TIERA67 67-gene candidate pool with 7-category 
classification (TDS_core 16, MAPK_output_ERK 10, Driver_anchor 12, 
Aggressive_marker 10, Dediff_invasion 10, Immune_stromal_light 5, 
Thyroid_lineage_extra 4). 8-gene RAI panel (SLC5A5/TPO/TG/TSHR/PAX8/
NKX2-1/FOXE1/DIO1) is a subset of TDS_core indicated by in_8gene_panel 
flag.
```

### S3 — GSE286332 top 300 DEGs

```
Suppl Table 3. Top 300 differentially expressed genes (PyDESeq2 Wald 
test, BH-FDR padj<0.05) in GSE286332 PTC+HT (n=9) vs PTC (n=9). 
Gene name, log2 fold change, Wald statistic, p-value, padj, and 
predicted category (Ig V/J/C, HLA-II, B-cell, T-cell, Other immune, 
Other) provided. Full DEG table (29,672 genes) at GitHub repository.
```

### S4 — Pan-Asian HLA per-allele forest

```
Suppl Table 4. Per-allele forest meta-analysis. For 8 alleles (HLA-A 
A*02:07, B*46:01, C*01:02, DPA1*02:02, DPB1*05:01, DQA1*02:01, 
DQB1*02:01, DRB1*07:01): Chu 2018 published GD case n + freq + ctrl 
n + freq + OR + 95% CI + p-value; Korean PTC pool (n=874) carrier 
count + freq + 95% Wilson CI; Korean PTC vs Chu ctrl OR + p; Korean 
PTC vs Chu GD OR + p; random-effects DerSimonian-Laird pooled OR + 
95% CI + I² + τ²; Korean K2/Lee/GSE286332-PTC sub-cohort breakdown 
+ Cochran's Q.
```

### S5/S5b — Mediation analysis + OLS decomposition

```
Suppl Table 5. Baron-Kenny mediation analysis (5,000-iteration 
bootstrap) for 4 candidate mediators of PTC+HT → P(DM1) ↓ pathway 
in GSE286332 (n=18): a (treat→mediator), b (mediator→outcome|treat), 
c_total, c_direct, indirect (a×b), %mediated, bootstrap 95% CI, 
empirical p_emp.

Suppl Table 5b. OLS regression decomposition. P_DM1 ~ HLA_II + 
8-gene + immune (z-standardized predictors): per-predictor β, SE, t, 
p; R²=0.756, R²_adj=0.704, F-statistic=14.50, F p=1.4e-4.
```

### S6 — TCGA Hashimoto-like × DM cluster cross-tabulation

```
Suppl Table 6. TCGA-THCA Hashimoto-like signature carrier × DM cluster 
cross-tabulation. Six threshold methods (GMM 2-component, Otsu, 
top-10%, top-20%, top-30%, residualized Otsu after Stromal+immune 
residualization) with per-method n_hashi positive, DM1/DM2 hashi 
positive percentage, Fisher exact OR, and Fisher p-value. Top-30% 
threshold yields strongest DM2 enrichment (OR=0.20, Fisher p=6.4×10⁻¹⁰).
```

### S7 — Cross-cohort generalization

```
Suppl Table 7. Cross-cohort Hashimoto-like signature prevalence. 
TCGA-THCA n=500 (GMM 18.0% / Otsu 19.6%); Korean GSE213647 n=632 
(GMM 22.8% / Otsu 28.2%); GSE286332 PTC+HT n=9 (GMM 100% / Otsu 100%). 
Sub-B-like signature (DM1 sub-B vs sub-A signature transfer) Korean 
GSE213647: GMM 47.2% / Otsu 52.5%.
```

### S8a/b/c — DM1 sub-cluster A vs B characterization

```
Suppl Table 8a. DM1 sub-cluster score profile. For 7 score modules 
(g8_RAI, TDS16, HLA_I, HLA_II, TLS, B_cell, Dediff_EMT): sub-A median, 
sub-B median, Cohen's d (sub-B - sub-A), MW p-value.

Suppl Table 8b. DM1 sub-cluster clinical phenotype. For 3 metrics 
(age_at_diagnosis, tumor_size_mm, os_days): sub-A mean, sub-B mean, 
Cohen's d, MW p-value.

Suppl Table 8c. DM1 sub-cluster mutation × Hashimoto. sub-A and sub-B 
counts of: BRAF+, RAS+, mutation-negative, Hashimoto+ (Otsu). 
Per-cluster mutation-negative percentage and Hashimoto-positive 
percentage.
```

---

## Quality check ✅

- [x] All 7 main figures captioned with source data path
- [x] All 6 suppl figures captioned
- [x] All 8 suppl tables captioned (S1-S8 with sub-tables)
- [x] All numerical values cross-verified with manuscript outline + summary.json
- [x] Cell Press caption format compliant
- [x] No defensive over-claiming
- [x] Voice-protected verb tone marked borderline (user polish)
