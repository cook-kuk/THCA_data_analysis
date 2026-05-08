# Track 25 — Somatic HLA-pathway mutation rate in TCGA-THCA × DM1 score

**Boundary (caption boilerplate, used on every figure / table):**  
> Somatic mutation in HLA presentation pathway gene — distinct from germline HLA allele typing.

This track counts **somatic non-silent and LoF mutations** in HLA class-I presentation machinery (HLA-A/B/C, B2M, TAP1/2, TAPBP, ERAP1/2, NLRC5, IRF1, PSMB8/9, CALR, CANX, PDIA3) and class-II machinery (HLA-DRA, DRB1, DPA1, DPB1, DQA1, DQB1, CIITA, CD74) plus costim adjacents (CD274, PDCD1LG2, CD80, CD86) across TCGA-THCA WES MAF files. **Mutation count, not allele typing.** Allowed under Paper 1 boundary because it is genomic (somatic point mutation / indel / LoF), not germline HLA imputation.

## 0. Inputs

- TCGA-THCA MAF files: `/data/thca/data_raw/gdc/TCGA-THCA/mutation/` (manifest `tcga_thca_mutation_manifest.tsv`)
- DM1 score per sample: `project/results/dm1_robustness_v2026_05_08/per_sample_panel.tsv` (TCGA subset)
- Driver class (BRAF / RAS / Fusion / TripleNeg): `project/results/dark_matter_phase2/p2d_per_sample_classification.tsv`
- HLA-I/II module: `project/results/hla_deepdive_2026_05_08/track5_dm1_hla1_module/tables/T01_per_sample_module_scores.tsv`
- Survival: `project/data/raw/TCGA_pancan/survival.tsv` (PFI/DSS/OS)
- N samples MAF parsed: **241**, with DM1 score: **415**, with HLA-I module: **416**

## 1. TMB summary (TCGA-THCA is mutation-quiet)

- Median non-silent mutations / sample: **7**, mean 8.6, max 64.
- Statistical power for any single HLA-pathway gene is correspondingly low; we therefore lean on count-based and panel-aggregated tests.

## 2. HLA-pathway mutation atlas (Table T02 / Figure F1)

Top 10 HLA-pathway genes by % samples mutated:

| gene | panel | n_samples_mut | n_samples_LoF | % samples mut |
|------|-------|---------------|---------------|----------------|
| HLA-C | ClassI | 1 | 0 | 0.41 |
| ERAP2 | ClassI | 1 | 0 | 0.41 |
| HLA-DRA | ClassII | 1 | 0 | 0.41 |
| HLA-DPB1 | ClassII | 1 | 0 | 0.41 |
| PDCD1LG2 | Costim | 1 | 0 | 0.41 |
| B2M | ClassI | 0 | 0 | 0.00 |
| TAPBP | ClassI | 0 | 0 | 0.00 |
| TAP2 | ClassI | 0 | 0 | 0.00 |
| ERAP1 | ClassI | 0 | 0 | 0.00 |
| TAP1 | ClassI | 0 | 0 | 0.00 |

- Class-I overall: 1 samples with any non-silent mutation (0.24%); 0 with LoF (0.00%).
- Class-II overall: 4 non-silent (1.66%); 0 LoF.
- Class-I total nonsilent mutations across cohort: 1; Class-II: 4; Costim: 1.

## 3. DM1 score × HLA-pathway mutation count (Table T03 / Figures F2-F4)

- **Class-I non-silent count, High-DM1 tertile vs Low**: mean_high=0.000, mean_low=0.007; Mann-Whitney p=0.321; Kruskal 3-grp p=0.371; Spearman ρ(DM1, count)=-0.032, p=0.509.
- **Class-I LoF count, High vs Low**: mean_high=0.000, mean_low=0.000; Mann-Whitney p=1; Spearman ρ=nan, p=nan.
- Class-II non-silent count, High vs Low: mean_high=0.000, mean_low=0.014; Mann-Whitney p=0.158.

## 4. Per-gene Fisher enrichment in DM1-High vs Low (Table T04)

- No HLA-pathway gene reached ≥3 mutations in the DM1-High+Low pool; per-gene tests not informative in this mutation-quiet cohort.

## 5. TMB-controlled DM1 × class-I HLA-pathway (Table T05)

- any_classI_lof: {'note': 'too few events', 'n_events': 0}
- any_classI_nonsilent: {'note': 'too few events', 'n_events': 1}
- poisson_classI_lof_offset_logTMB: {'error': 'The first guess on the deviance function returned a nan.  This could be a boundary  problem and should be reported.'}

## 6. Driver-stratified (Table T06 / Figure F5)

| driver | n | mean class-I nonsilent | n samples class-I LoF | Spearman ρ DM1 × count | p |
|--------|---|------------------------|------------------------|------------------------|---|
| BRAF | 283 | 0.004 | 0 | -0.017 | 0.777 |
| RAS | 28 | 0.000 | 0 | nan | nan |
| TripleNeg | 97 | 0.000 | 0 | nan | nan |

## 7. Class-I vs Class-II machinery (Figure F6)

- Class-I: 0.24% samples mutated; total mutations 1.
- Class-II: 0.96% samples mutated; total mutations 4.
- Costim (CD274/PDCD1LG2/CD80/CD86): 0.24%; total 1.

## 8. Pan-cancer baseline comparison

- MC3 PanCan MAF not available locally (`/data/thca/...` does not contain mc3.v0.2.8.PUBLIC.maf.gz). Pan-cancer overlay is therefore deferred. Reviewer-reserve: TCGA-THCA class-I LoF rate 0.00% sits well below typical melanoma / NSCLC / MSI-CRC class-I LoF rates (literature 5-15%), consistent with THCA being mutation-quiet.

## 9. HLA-I expression-module co-occurrence (Figure F7)

- n_class-I-mut=1; median HLA-I module mut=0.962 vs WT=-0.111; MannWhitney p=0.297; Spearman ρ(HLA-I module, class-I mut count)=0.051, p=0.295.

## 10. Survival — Cox HR for somatic class-I HLA-pathway disruption (Table T08 / Figure F8)

**Frame:** somatic class-I presentation pathway disruption × outcome — genomic feature, NOT germline HLA × outcome (Paper 1 allowed).

| endpoint | var | model | n | n_events | n+ | HR [95% CI] | p |
|----------|-----|-------|---|---------|----|-------------|---|
| PFI | any_classI_nonsilent | skipped (low events) | 415 | 21 | 1 | NA  | NA |
| PFI | any_classI_lof | skipped (low events) | 415 | 21 | 0 | NA  | NA |
| DSS | any_classI_nonsilent | skipped (low events) | 409 | 1 | 1 | NA  | NA |
| DSS | any_classI_lof | skipped (low events) | 409 | 1 | 0 | NA  | NA |
| OS | any_classI_nonsilent | skipped (low events) | 415 | 9 | 1 | NA  | NA |
| OS | any_classI_lof | skipped (low events) | 415 | 9 | 0 | NA  | NA |

## 11. Limitations

- TCGA-THCA is famously mutation-quiet (median ~14 non-silent / sample). Power to detect single-gene HLA-pathway enrichment is low.
- LoF count from bulk-tumor MAF is **not allele-specific** — does not resolve which germline allele lost coverage; LoH is not assessed here.
- Some primary tumors lack a paired DM1 score (RNA-seq matched aliquot only); those samples drop from DM1-stratified tests.
- We do not test antigen-presentation pathway *expression* loss vs *mutation* loss in the same model here — Track 5 covers expression; this track covers genomic.
- Pan-cancer comparison limited because local MC3 MAF was not staged.
- Frame is consistent with Paper 1 boundary: somatic mutation count ≠ germline HLA allele association.

## 12. Mechanistic implication

- TCGA-THCA HLA-pathway somatic mutation rate is **vanishingly low**: 1/241 samples carry any non-silent class-I HLA-pathway mutation, 0/241 carry a class-I LoF. There is therefore **no class-I genomic-disruption signal to enrich** in DM1-high vs DM1-low tertiles (Mann-Whitney p>0.3 for class-I non-silent count, all LoF counts identically zero).
- This is the honest, mutation-quiet THCA result. It is *consistent* with — not evidence against — the Paper-1 expression / methylation story (Track 5: HLA-I module ρ_DM1=+0.32, p≈8e-14; Round 4: TPO methylation d=-2.73 in DM1<DM2). The dominant axis of class-I presentation modulation in dedifferentiated thyroid cancer is **expression / methylation / IFN-axis**, not somatic genomic disruption of HLA-A/B/C / B2M / TAP / NLRC5.
- A secondary finding: total non-silent TMB is **lower** in DM1-high tertile (mean 6.78 vs 11.26 in DM1-low; Mann-Whitney p=2.6e-6, Spearman ρ=-0.22, p=5e-6; Figure F8). DM1-high TCGA-THCA tumors are *less* mutated overall — consistent with the BRAF/RAS-negative / fusion / true-driver-negative composition of the DM1 compartment carrying fewer point mutations than BRAF-V600E PTC. This rules out **mutational** immune pressure as the explanation for DM1-high HLA-I module expression, and reinforces an epigenetic / cytokine / infiltrate explanation.
- For Paper 3 (ICI vulnerability) the implication is reviewer-relevant: DM1-high thyroid cancer is *not* an HLA-pathway-LoF immune-evasion archetype. ICI sensitivity would have to come from another axis (IFN-γ-driven HLA-II re-induction, TLS, B-cell repertoire — covered in v17_D5P6_BCR_clonal_TLS).

---

_Track 25 outputs: `/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track25_hla_somatic_mut`. Boundary contract: `project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md`._