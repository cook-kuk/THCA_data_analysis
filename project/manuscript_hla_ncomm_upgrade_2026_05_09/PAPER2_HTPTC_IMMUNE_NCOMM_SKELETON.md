# Paper 2 NComm-Style Skeleton

## Working Title

**Antigen-presentation and tertiary-lymphoid niches define an immune-active Hashimoto-overlap state in papillary thyroid carcinoma**

## One-Sentence Claim

HT-overlap and immune-active PTC are best modeled as a transcriptomic antigen-presentation/TLS/IFN-gamma state; current HLA-II allele data nominate prospective candidates but do not validate allele-level cancer association.

## Abstract Draft

Hashimoto-overlap papillary thyroid carcinoma is often discussed through autoimmune and HLA susceptibility, but tumor-cohort HLA allele evidence is underpowered and vulnerable to ancestry and tumor/germline confounding. We therefore modeled HT-overlap PTC as an immune-context state rather than an allele-association claim. In TCGA-THCA, DM1/immune-active tumors showed coordinated HLA-I and HLA-II gene-expression module induction, partially mediated by IFN-gamma and inversely coupled to thyroid differentiation/RAI modules. Independent HT-overlap bulk and single-cell cohorts supported HLA/AP-TLS elevation, while a 632-sample Korean bulk cohort and a 16-slide thyroid Visium progression cohort generalized the expression ecology beyond the discovery setting. Single-cell and spatial analyses localized HLA-II to APC-rich compartments and connected TLS niches, with niche-level HLA-II elevation across thyroid tumor stages. Somatic HLA-presentation gene loss, selective HLA-G escape, and neoantigen burden did not explain the axis. In the only explicit HT-overlap HLA subset, DRB1*04:05 emerged as a prospective candidate but remained FDR-negative and severely underpowered. These results define a data-supported antigen-presentation immune state and a validation roadmap for germline HLA studies in HT-overlap PTC.

## Results Spine

### Result 1 — TCGA-THCA shows coordinated HLA-I/HLA-II module induction in the DM1/immune-active state

Main message: HLA in this paper is transcriptomic module biology, not allele genotype.

Key numbers:
- DM1 vs HLA-I rho 0.318, p=8.3e-14.
- DM1 vs HLA-II rho 0.393, p=7.3e-21.
- Partial immune+stromal rho remains HLA-I 0.266 and HLA-II 0.390.
- Sample-unique external-validation network: DM1-HLA-II rho 0.400 (n=505), IFNG-HLA-II rho 0.874, TLS-HLA-II rho 0.760, RAI-HLA-II rho -0.479.

Main figure: `track5_dm1_hla1_module/figs/F01_dm1_vs_hla1_hla2_scatter.png`

Extra validation figure: `project/results/hla_two_paper_synthesis_2026_05_09/extra_external_validation/figures/F02_tcga_hla_ifng_tls_rai_network.png`

### Result 1b — External HT-overlap validation shows antigen-presentation/TLS elevation in PTC+HT

Main message: the HT-overlap story is now supported by expression and TLS scores, not by allele association.

Key numbers:
- GSE286332 PTC+HT vs PTC-only: HLA-I Cohen d 2.34, FDR 0.0015.
- GSE286332 PTC+HT vs PTC-only: HLA-II Cohen d 3.65, FDR 0.0012.
- GSE286332 PTC+HT vs PTC-only: TLS score Cohen d 3.01, FDR 0.0015.

Main figure: `project/results/hla_two_paper_synthesis_2026_05_09/extra_external_validation/figures/F01_gse286332_ht_overlap_hla_tls_effects.png`

Key table: `project/results/hla_two_paper_synthesis_2026_05_09/extra_external_validation/tables/T01_gse286332_ht_overlap_hla_tls_effects.tsv`

### Result 1c — Exact small-n robustness and HLA-II dominance close the two obvious reviewer attacks

Main message: the n=9 vs n=9 HT-overlap result survives exact label permutation, and the TCGA HLA-II module is explained primarily by immune/IFN context rather than a loose DM1-only correlation.

Key numbers:
- GSE286332 HLA-II exact label permutation p=0.000062; bootstrap d 95% CI 2.51-7.56.
- AP/TLS composite d=3.08, exact p=0.000185, descriptive AUC=0.96.
- TCGA HLA-II dominance full model R2=0.83; largest drop-one contributors are immune proxy and IFNG Hallmark.

Main figures:
- `project/results/hla_two_paper_synthesis_2026_05_09/level_up_wave2/figures/F03_gse286332_exact_permutation_bootstrap.png`
- `project/results/hla_two_paper_synthesis_2026_05_09/level_up_wave2/figures/F05_tcga_hla2_driver_dominance.png`

### Result 1d — Independent Korean bulk and thyroid Visium cohorts generalize the HLA/AP expression ecology

Main message: the antigen-presentation ecology is not a single-cohort artifact, but these datasets are not HT-specific.

Key numbers:
- GSE213647 Korean bulk RNA-seq, 632 samples: PTC tumor vs PTC normal HLA-II/AP d=1.39, FDR=4.03e-40; Myeloid/DC d=1.23, FDR=3.34e-36; HLA-I d=1.22, FDR=5.08e-36; AP/TLS d=0.89, FDR=3.85e-21.
- GSE250521 thyroid Visium, 16 slides and 55,873 post-QC spots: stage trend N/PTC/LPTC/ATC CD74/MIF rho=0.76 q=0.0027, HLA-I rho=0.75 q=0.0027, Myeloid/DC rho=0.70 q=0.0055, HLA-II/AP rho=0.61 q=0.018.
- GSE6004 is a negative/specificity stress-test, not a positive result: invasion vs normal HLA-II/AP d=-0.42 and AP/TLS d=-0.59, FDR=0.829.
- GSE184362/GSE191288 are low-tier scRNA generalization only because HT labels are unavailable in GEO metadata/raw filenames.

Main figures:
- `project/results/hla_two_paper_synthesis_2026_05_09/gse213647_korean_bulk_validation/figures/F02_gse213647_ptc_effect_sizes.png`
- `project/results/hla_two_paper_synthesis_2026_05_09/gse250521_thyroid_visium_hla_validation/figures/F02_gse250521_stage_trend_rho.png`

### Result 2 — IFN-gamma explains the antigen-presentation induction axis

Main message: IFN-gamma is the strongest mediator of HLA module induction, but the model is observational.

Key numbers:
- IFN-gamma vs HLA-I rho 0.880.
- IFN-gamma vs HLA-II rho 0.873.
- Hallmark-mediated fraction: HLA-I 107.83%, HLA-II 70.22%.

Main figure: `track26_ifng_hla_dm1/figs/F03_mediation_pct_mediated.png`

### Result 3 — Single-cell and spatial maps localize HLA-II to APC-rich compartments and tumor-stage acquisition

Main message: bulk HLA-II is not dismissed as contamination; it is decomposed into cell and niche sources.

Key data:
- Pu scRNA: 19,999 cells.
- Lu scRNA: 17,898 cells.
- GSE163203 HT/PTC scRNA: 110,000 cells across 10 biological samples, with PTCwithHT vs PTCwithoutHT AP/TLS d=2.78 and B/TLS d=2.59.
- DC/Myeloid/B dominate HLA-II module.
- Pu/Lu cell-type HLA-II concordance rho 0.95, p=2.28e-05; shared top producers are B_cell, DC and Myeloid.
- Visium stage heatmap shows tumor-stage HLA-II acquisition.

Main figure: `track12_scrna_hla_celltype/figs/F7_visium_HLA_stage_celltype_heatmap.png`

Extra validation figure: `project/results/hla_two_paper_synthesis_2026_05_09/extra_external_validation/figures/F03_scrna_hla2_celltype_concordance.png`

### Result 4 — TLS niches are HLA-II-hot in thyroid tumors

Main message: niche-level analysis resolves spot-level dilution artifacts and gives the most visually compelling Paper 2 result.

Key numbers:
- 608 TLS niches across 16 Visium samples.
- Global niche-level Cohen d 0.655.
- PTC stage mean d 0.935; ATC stage mean d 0.916.
- HLA-II in-niche positive direction in 16/16 samples; cancer subset mean d 0.826.

Main figure: `track33_spatial_tls_hla/figs/F2_forest_cohen_d_per_sample.png`

Extra validation figure: `project/results/hla_two_paper_synthesis_2026_05_09/extra_external_validation/figures/F04_spatial_tls_hla_stage_effect.png`

### Result 5 — HT-overlap HLA-II alleles are candidates, not findings

Main message: keep allele-level data honest and use them to justify prospective validation.

Key numbers:
- GSE286332 explicit HT-overlap allele test: n=9 vs n=9.
- DRB1*04:05 OR 3.93, 95% CI 1.06-14.53, Fisher p=0.0866, FDR=0.346.
- OR=1.8 validation requires about 633/group.

Main figure: `track4_ht_ptc_hla2/plots/F3_HT_overlap_PTC_vs_baseline.png`

### Result 6 — Negative controls prevent immune-escape overclaim

Main message: the immune-context axis is not explained by easy alternative mechanisms.

Key controls:
- Somatic class-I LoF events: 0.
- HLA-G vs DM1 rho 0.428, but HLA-G decoupling vs DM1 rho 0.006, p=0.887.
- TMB/neoantigen link is weak relative to direct DM1-HLA-II association.

Main figure: `track27_hla_g_nonclassical/figs/F04_decoupling_vs_dm1.png`

### Result 7 — A BRAF-rich DM1-high/IFN-low subgroup suggests immune-induction failure

Main message: useful biology, not yet clinical outcome.

Key numbers:
- fail-to-induce subgroup n=97/527.
- BRAF OR 3.78, FDR=3.2e-07.
- survival endpoints weak/not robust.

Main figure: `track28_immune_cold_thca/figs/F01_quadrant_scatter.png`

## Main Figure Legend Drafts

**Figure 1. HLA-I and HLA-II gene-expression module induction in immune-active THCA.** TCGA-THCA primary tumors scored for DM1 and antigen-presentation modules show positive HLA-I and HLA-II module association. HLA scores are transcriptomic modules and not HLA allele genotypes.

**Figure 1b. External HT-overlap and module-network validation.** GSE286332 PTC+HT tumors show higher HLA-I, HLA-II and TLS scores than PTC-only tumors. A sample-unique TCGA network links DM1, IFN-gamma, HLA-II, TLS and inverse RAI module scores. All HLA measurements in cancer cohorts are expression modules.

**Figure 1c. Exact robustness and dominance.** Exact label permutation and bootstrap confidence intervals quantify the small-n HT-overlap signal without relying on asymptotic tests. A TCGA multivariable dominance panel ranks immune/IFN context as the main HLA-II expression driver.

**Figure 1d. Independent Korean bulk and thyroid Visium generalization.** GSE213647 Korean bulk RNA-seq validates tumor-associated HLA-II/AP, HLA-I and Myeloid/DC module elevation at large sample size. GSE250521 Visium shows stage-ordered acquisition of HLA/AP/CD74/myeloid spatial ecology across normal thyroid, PTC, locally advanced PTC and ATC. These are expression generalization datasets, not HT-specific or HLA allele analyses.

**Figure 2. IFN-gamma mediation of HLA module induction.** IFN-gamma Hallmark and Ayers TIS signatures are compared with HLA-I/HLA-II modules; mediation and conditional-independence analyses support IFN-gamma as the dominant observed mediator.

**Figure 3. Cell-source decomposition of HLA-II.** scRNA and Visium analyses localize HLA-II to DC/Myeloid/B compartments and tumor-stage spatial programs, avoiding a bulk-only contamination interpretation.

**Figure 4. TLS niches are HLA-II-hot.** Connected-component TLS niches from Visium spots show elevated HLA-II relative to sample background across PTC, LPTC and ATC, with distance-gradient support.

**Figure 5. HLA-II allele candidates for prospective validation.** Explicit HT-overlap PTC HLA-II allele tests remain underpowered and FDR-negative; DRB1*04:05 is nominated only as a prospective germline NGS target.

**Figure 6. Negative controls against overclaim.** Somatic HLA-presentation mutations, selective HLA-G decoupling and neoantigen burden are evaluated as alternative explanations and do not account for the observed HLA-II/IFN-gamma immune-context axis.

**Figure 7. Fail-to-induce subgroup.** A DM1-high/IFN-low subgroup enriched for BRAF tumors suggests a biologically distinct immune-induction failure state without a current clinical-outcome claim.

## Discussion Upgrade Paragraph

The strongest interpretation is that HT-overlap PTC belongs to a broader antigen-presentation immune-state continuum. The new GSE286332/GSE138198/GSE163203 reanalyses raise the HT-overlap claim because PTC+HT is elevated for HLA/AP and TLS scores even when allele claims are kept out of the main result. GSE213647 and GSE250521 then broaden the evidence into Korean bulk RNA-seq and spatial tumor progression, while GSE6004 and GSE184362/GSE191288 define the specificity and labeling limits. This moves the paper away from fragile HLA allele association and toward a multi-modal immune ecology supported by bulk RNA-seq, sample-unique TCGA network analysis, scRNA, spatial transcriptomics, mediation analysis and negative controls. The HLA-II allele data are useful precisely because they define the validation experiment: germline NGS HLA in a pathology-confirmed PTC+HT versus PTC-only cohort, powered around DRB1*04:05 or a pre-specified class-II panel.

## Do Not Claim

- Do not claim HLA allele-based thyroid cancer risk.
- Do not claim HT causes PTC or dedifferentiation.
- Do not claim DRB1*04:05 is validated.
- Do not claim DPB1*05:01 is an HT-PTC finding.
- Do not join HLA allele calls to survival, RAI response, stage, BRAF/RAS/TERT, DM1/DM2, recurrence or patient selection.
- Do not call IFN-gamma mediation causal without perturbation data.
