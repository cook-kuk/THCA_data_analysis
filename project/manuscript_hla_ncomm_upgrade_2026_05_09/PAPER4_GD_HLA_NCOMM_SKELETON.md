# Paper 4 NComm-Style Skeleton

## Working Title

**Pan-Asian HLA architecture of autoimmune thyroid disease reveals an AITD-broad DPB1*05:01 axis and a Graves-specific C*01:02 signal**

## One-Sentence Claim

Pan-Asian Graves/AITD risk is organized by a DPB1*05:01-centered class-II axis and a partly independent class-I block, with C*01:02 showing the strongest Graves-specificity and DPB1*05:01 behaving as an AITD-broad rather than Graves-exclusive allele.

## Abstract Draft

Autoimmune thyroid disease has a strong HLA component, but East-Asian allele effects, Korean baseline frequencies, disease specificity, and tissue mechanism remain fragmented across historical case-control studies and reference panels. We curated Korean and Pan-Asian HLA evidence for Graves disease and autoimmune thyroid disease, harmonized carrier and allele-frequency metrics, built a healthy Korean/Pan-Asian HLA atlas, compared candidate alleles across autoimmune diseases and public trans-ancestry MHC summary statistics, and reanalyzed independent thyroid tissue transcriptomic/spatial datasets for HLA/AP mechanism. DPB1*05:01, B*46:01, A*02:07 and C*01:02 formed the Pan-Asian Graves/AITD core. DPB1*05:01 showed a robust Pan-Asian random-effects signal but broader AITD distribution, whereas C*01:02 had the highest Graves-specificity score in the cross-autoimmune panel. Healthy Korean HLA references and trans-ancestry MHC signals supported population-aware interpretation while highlighting portability limits. AITD spatial and array datasets showed concordant HLA-II/AP, CD74/MIF and TLS/B-cell tissue activation. These results define a Korean-anchored HLA architecture for autoimmune thyroid disease and establish the adult Korean germline NGS cohort required for definitive validation.

## Results Spine

### Result 1 — A harmonized Pan-Asian Graves HLA meta-analysis identifies the core alleles

Main message: DPB1*05:01, B*46:01, A*02:07 and C*01:02 are the core Pan-Asian GD/AITD alleles.

Key numbers:
- DPB1*05:01 v3 pooled OR 2.10, 95% CI 1.70-2.60, k=4.
- B*46:01 pooled OR 2.17, 95% CI 1.39-3.36, but I2=85.09%, so report as heterogeneous.
- A*02:07 pooled OR 2.12, 95% CI 1.73-2.61, upgraded to full Pan-Asian meta-ready.
- C*01:02 pooled OR 1.88, 95% CI 1.58-2.24.

Main figure: `track10_korean_lit/plots/forest_panasian_GD_v3.png`

### Result 2 — Korean literature extraction upgrades the evidence but exposes the adult-GD NGS gap

Main message: Korean pediatric/adult literature strengthens the Pan-Asian architecture, but a de novo Korean adult GD NGS cohort is still missing.

Key framing:
- A*02:07 grade changed from partial to full Pan-Asian meta-ready.
- DRB1*07:01 protective signal now has k=2 support.
- Korean adult NGS is the A+ blocker, not another meta-analysis pass.

Main figure: `track10_korean_lit/plots/forest_korean_only_v2.png`

### Result 3 — A healthy Korean/Pan-Asian HLA atlas anchors ancestry-aware interpretation

Main message: population baseline is not a nuisance; it is part of the central result.

Key numbers:
- Healthy atlas: 9 populations, 145 alleles.
- Korean platform concordance r≈0.86-0.93 across Baek2021 NGS, K2 normal RNA-seq, and GSE213647 normal RNA-seq.

Main figure: `track3_pan_asian_atlas/F9_cross_korean_concordance.png`

### Result 4 — Cross-autoimmune specificity separates GD-specific and AITD-broad alleles

Main message: DPB1*05:01 is not Graves-specific; C*01:02 is the best GD-specific candidate.

Key numbers:
- C*01:02 GD specificity score 14.97.
- DPB1*05:01 GD specificity score 2.89; AITD breadth ratio 3.31.

Main figure: `track31_cross_autoimmune/figures/F01_pleiotropy_heatmap.png`

### Result 5 — Trans-ancestry MHC signals corroborate class-II thyroid autoimmunity but do not transfer Korean effect sizes

Main message: FinnGen/Pan-UKBB support MHC class-II thyroid autoimmunity but also demonstrate why ancestry-specific interpretation is mandatory.

Key numbers:
- FinnGen GD strict rs9277534 beta +0.129, p=3.6e-07.
- FinnGen HT-AI rs9277534 beta -0.107, p=7.1e-42.
- Extra endpoint-sensitivity layer: strict Graves DPB1*05:01 proxy has 1 positive signal at p<1e-6 and 3 at p<1e-4; HT/hypothyroid has 5 genome-wide negative signals.
- Broad hyperthyroid/self-report endpoints show a negative genome-wide proxy signal, so they should be shown as phenotype-label sensitivity rather than mixed into strict Graves.

Main figure: `track11_finngen_panukbb/figures/F08_dpb1_0501_transancestry_forest.png`

Extra validation figure: `project/results/hla_two_paper_synthesis_2026_05_09/extra_external_validation/figures/F05_paper4_gwas_mhc_share.png`

Extra validation tables:
- `project/results/hla_two_paper_synthesis_2026_05_09/extra_external_validation/tables/T07c_track1_tag_snp_branch_threshold_summary.tsv`
- `project/results/hla_two_paper_synthesis_2026_05_09/extra_external_validation/tables/T08_afnd_korean_baseline_context.tsv`

### Result 5b — GWAS Catalog and AFND provide external calibration, not allele-level replication

Main message: public external resources now support the MHC-centered framing and population calibration, while keeping replication language precise.

Key numbers:
- GWAS Catalog GD MHC GWS share 0.128 (36/282).
- GWAS Catalog hypothyroidism MHC GWS share 0.305 (981/3213).
- AFND DPB1*05:01 South Korea weighted frequency 0.367; East-Asia pool 0.375.

Main figure: `project/results/hla_two_paper_synthesis_2026_05_09/extra_external_validation/figures/F05_paper4_gwas_mhc_share.png`

### Result 5c — Allele triangulation turns the synthesis into a reviewer-facing prioritization map

Main message: this is the level-up figure that prevents the paper from reading like a list of historical HLA associations.

Key numbers:
- C*01:02 is the specificity-led GD candidate in the triangulation map.
- DPB1*05:01 remains the AITD-broad class-II anchor, not the most GD-specific allele.
- Endpoint robustness heatmap separates strict Graves, autoimmune hyperthyroid, broad hyperthyroid, HT/hypothyroid and non-thyroid autoimmunity.

Main figures:
- `project/results/hla_two_paper_synthesis_2026_05_09/level_up_wave2/figures/F01_paper4_allele_triangulation_map.png`
- `project/results/hla_two_paper_synthesis_2026_05_09/level_up_wave2/figures/F02_paper4_endpoint_robustness_heatmap.png`

### Result 5d — Independent AITD tissue datasets validate the HLA-II/AP mechanism layer

Main message: the HLA genetics are now linked to tissue antigen-presentation biology, without calling expression data genotype replication.

Key numbers:
- GSE248205 AITD Visium, 8 samples and 16,985 in-tissue spots: AITD vs control CD74/MIF d=3.62, B/TLS d=2.57, AP/TLS d=2.26; HT spot burden above control p90: HLA-II/AP 0.964, B/TLS 0.997, AP/TLS 0.993.
- GSE29315 thyroid array, 71 samples: HT vs thyroid hyperplasia HLA-II/AP d=5.45, exact p=0.0007, FDR=0.0023; HLA-I d=4.40, AP/TLS d=3.52.

Main figures:
- `project/results/hla_two_paper_synthesis_2026_05_09/gse248205_aitd_spatial_validation/figures/F03_gse248205_spot_burden_heatmap.png`
- `project/results/hla_two_paper_synthesis_2026_05_09/gse29315_aitd_array_validation/figures/F02_gse29315_ht_effect_sizes.png`

### Result 6 — Tolerance and peptide-binding analyses define a functional follow-up model

Main message: mechanism is plausible but not yet functionally proven.

Use as discussion-facing or final main figure depending on venue:
- AIRE/thymus signal.
- thyroid antigen panel.
- peptide-binding hypothesis.

Main figure: `track19_aire_tolerance/plots/F8_DPB1_05_01_hypothesis_schematic.png`

## Main Figure Legend Drafts

**Figure 1. Pan-Asian Graves disease HLA random-effects meta-analysis.** Random-effects forest plot of curated Pan-Asian Graves/AITD HLA association rows after carrier/allele metric harmonization. DPB1*05:01, B*46:01, A*02:07 and C*01:02 form the core risk architecture. B*46:01 is reported with heterogeneity caveat.

**Figure 2. Korean evidence upgrade and replication readiness.** Korean-only forest and readiness scoreboard showing which alleles are publishable as synthesis-level evidence and which require Korean adult NGS replication.

**Figure 3. Healthy Korean/Pan-Asian HLA atlas.** Cross-platform concordance among Korean NGS and RNA-seq-imputed normal HLA references, used to anchor population-frequency interpretation.

**Figure 4. Cross-autoimmune specificity.** Heatmap and specificity summaries across autoimmune diseases showing C*01:02 as the most GD-specific focus allele and DPB1*05:01 as AITD-broad.

**Figure 5. Trans-ancestry MHC sensitivity.** FinnGen and Pan-UKBB MHC summary-statistic evidence at Track-1 tag SNPs, supporting class-II thyroid autoimmunity while highlighting ancestry and phenotype-label limitations. Strict Graves, autoimmune hyperthyroid, broad hyperthyroid and HT/hypothyroid endpoints are separated to avoid phenotype mixing.

**Figure 5b. External MHC and AFND calibration.** GWAS Catalog MHC share and AFND baseline-frequency panels show that thyroid autoimmunity is MHC-concentrated and that Korean/East-Asian allele frequencies must be explicitly calibrated. This is external context, not de novo Korean replication.

**Figure 5c. Allele triangulation and endpoint robustness.** Candidate alleles are mapped by Pan-Asian meta-analysis effect, GD-specificity, heterogeneity, AFND context and endpoint-specific tag-SNP behavior. The figure separates specificity-led candidates from broad AITD anchors and makes phenotype sensitivity explicit.

**Figure 5d. Tissue mechanism validation.** Independent AITD spatial and array transcriptomic datasets show HLA-II antigen-presentation, CD74/MIF and TLS/B-cell module expansion in autoimmune thyroid tissue. These data validate tissue mechanism, not allele-level replication.

**Figure 6. Mechanistic model.** Tolerance and thyroid-antigen presentation hypothesis linking AIRE/thymus expression, thyroid antigen availability, and HLA peptide-binding predictions. This figure is explicitly hypothesis-generating.

## Discussion Upgrade Paragraph

The central advance is not the rediscovery of isolated HLA alleles, but a population-resolved map that separates five often-conflated signals: Pan-Asian AITD susceptibility, strict Graves enrichment, broad thyroid endpoint heterogeneity, autoimmune pleiotropy, and tissue antigen-presentation mechanism. This distinction prevents DPB1*05:01 from being overinterpreted as Graves-specific and instead identifies C*01:02 as a sharper disease-specific candidate. The external validation layer strengthens the MHC-centered architecture while showing why East-Asian HLA work cannot be safely projected from European GWAS tags and why a Korean adult GD germline NGS case-control cohort is the decisive next step.

## Do Not Claim

- Do not claim Korean adult GD discovery without new Korean adult case-control NGS.
- Do not transfer GD/AITD HLA effects to PTC risk.
- Do not call DPB1*05:01 Graves-specific.
- Do not report B*46:01 as a single universal effect without heterogeneity caveat.
- Do not treat trans-ancestry tag-SNP effects as allele-level Korean replication.
- Do not call GSE248205 or GSE29315 genotype replication; they are tissue-expression mechanism support.
