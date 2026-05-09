# HLA 두 논문 NComm 업그레이드 패키지

**Generated:** 2026-05-09  
**목표:** 현재 데이터에서 논문 수준을 최대로 올리되, HLA allele-cancer boundary를 깨지 않는다.

## 한 줄 전략

- **Paper 4는 genetic architecture 논문으로 키운다.** 메인 강점은 Pan-Asian GD/AITD HLA meta + Korean literature extraction + healthy Korean HLA atlas + cross-autoimmune specificity + trans-ancestry MHC sensitivity다.
- **Paper 2는 immune ecology 논문으로 키운다.** 메인 강점은 HLA-I/II expression module, IFN-gamma mediation, scRNA/spatial HLA-II cell source, TLS-HLA-II niche, negative controls다. Allele association은 validation roadmap으로 낮춘다.

## 추천 제목

| Paper | Rank | Title |
|---|---:|---|
| Paper 4 | 1 | Pan-Asian HLA architecture of autoimmune thyroid disease reveals an AITD-broad DPB1*05:01 axis and a Graves-specific C*01:02 signal |
| Paper 4 | 2 | Population-resolved HLA architecture of Graves disease across East Asia |
| Paper 4 | 3 | A Korean-anchored synthesis of Pan-Asian HLA susceptibility in autoimmune thyroid disease |
| Paper 2 | 1 | Antigen-presentation and tertiary-lymphoid niches define an immune-active Hashimoto-overlap state in papillary thyroid carcinoma |
| Paper 2 | 2 | IFN-gamma-coupled HLA-II and TLS spatial niches in immune-active thyroid carcinoma |
| Paper 2 | 3 | Transcriptomic HLA induction and TLS niche architecture in Hashimoto-overlap papillary thyroid carcinoma |

## Paper 4 Main Figure Plan

- Fig. 1: **Pan-Asian Graves HLA random-effects architecture** — DPB1*05:01, B*46:01, A*02:07 and C*01:02 define the publishable Pan-Asian GD/AITD HLA core.  
  Asset: `project/results/hla_deepdive_2026_05_08/track10_korean_lit/plots/forest_panasian_GD_v3.png` (yes)
- Fig. 2: **Korean literature extraction and readiness upgrade** — Korean pediatric/adult literature upgrades A*02:07 and DRB1*07:01 while exposing the adult-GD NGS gap.  
  Asset: `project/results/hla_deepdive_2026_05_08/track10_korean_lit/plots/forest_korean_only_v2.png` (yes)
- Fig. 3: **Healthy Korean/Pan-Asian HLA atlas** — Healthy HLA atlas anchors population structure and shows Korean platform concordance for population-level frequency work.  
  Asset: `project/results/hla_deepdive_2026_05_08/track3_pan_asian_atlas/F9_cross_korean_concordance.png` (yes)
- Fig. 4: **Cross-autoimmune specificity map** — C*01:02 is most GD-specific; DPB1*05:01 is AITD-broad, not GD-exclusive.  
  Asset: `project/results/hla_deepdive_2026_05_08/track31_cross_autoimmune/figures/F01_pleiotropy_heatmap.png` (yes)
- Fig. 5: **Trans-ancestry MHC sensitivity** — FinnGen/Pan-UKBB corroborate MHC class-II thyroid-autoimmunity signals while showing ancestry-specific direction and portability limits.  
  Asset: `project/results/hla_deepdive_2026_05_08/track11_finngen_panukbb/figures/F08_dpb1_0501_transancestry_forest.png` (yes)
- Fig. EV1: **External MHC and endpoint-sensitivity validation** — GWAS Catalog MHC concentration, FinnGen/PanUKBB strict-GD versus broad-hyperthyroid sensitivity, and AFND baseline frequencies raise the external-validation layer without overcalling Korean replication.  
  Asset: `project/results/hla_two_paper_synthesis_2026_05_09/extra_external_validation/figures/F05_paper4_gwas_mhc_share.png` (yes)
- Fig. EV2: **Allele triangulation and endpoint robustness** — Meta-analysis, specificity, heterogeneity, AFND context and endpoint-specific tag-SNP behavior are integrated into a reviewer-facing allele-prioritization map.  
  Asset: `project/results/hla_two_paper_synthesis_2026_05_09/level_up_wave2/figures/F01_paper4_allele_triangulation_map.png` (yes)
- Fig. 6: **Mechanistic hypothesis: tolerance and thyroid-antigen presentation** — AIRE/thymus and peptide-binding analyses motivate functional follow-up without claiming in-vivo causality.  
  Asset: `project/results/hla_deepdive_2026_05_08/track19_aire_tolerance/plots/F8_DPB1_05_01_hypothesis_schematic.png` (yes)

## Paper 2 Main Figure Plan

- Fig. 1: **DM1/immune-active state and HLA-I/II module induction** — TCGA-THCA shows a reproducible HLA-I/HLA-II expression-module axis, not allele genotype.  
  Asset: `project/results/hla_deepdive_2026_05_08/track5_dm1_hla1_module/figs/F01_dm1_vs_hla1_hla2_scatter.png` (yes)
- Fig. EV1: **External HT-overlap HLA/TLS validation and module network** — GSE286332 validates higher HLA-I, HLA-II and TLS scores in PTC+HT; sample-unique TCGA network connects DM1, IFNG, HLA-II, TLS and inverse RAI.  
  Asset: `project/results/hla_two_paper_synthesis_2026_05_09/extra_external_validation/figures/F01_gse286332_ht_overlap_hla_tls_effects.png` (yes)
- Fig. EV2: **Exact small-n robustness and HLA-II dominance** — GSE286332 HT-overlap effects survive exact label permutation, and TCGA multivariable dominance ranks the immune/IFNG contributors to HLA-II expression.  
  Asset: `project/results/hla_two_paper_synthesis_2026_05_09/level_up_wave2/figures/F03_gse286332_exact_permutation_bootstrap.png` (yes)
- Fig. 2: **IFN-gamma mediation of antigen-presentation induction** — IFN-gamma explains most HLA-I induction and a major component of HLA-II induction.  
  Asset: `project/results/hla_deepdive_2026_05_08/track26_ifng_hla_dm1/figs/F03_mediation_pct_mediated.png` (yes)
- Fig. 3: **Single-cell and spatial cell-source decomposition** — DC/Myeloid/B compartments dominate HLA-II, with tumor-stage acquisition in malignant/spatial compartments.  
  Asset: `project/results/hla_deepdive_2026_05_08/track12_scrna_hla_celltype/figs/F7_visium_HLA_stage_celltype_heatmap.png` (yes)
- Fig. 4: **TLS-HLA-II spatial niche coupling** — Connected TLS niches are HLA-II-hot across thyroid tumor stages, resolving spot-level dilution artifacts.  
  Asset: `project/results/hla_deepdive_2026_05_08/track33_spatial_tls_hla/figs/F2_forest_cohen_d_per_sample.png` (yes)
- Fig. 5: **HT-overlap allele candidates as prospective validation targets** — DRB1*04:05 is the only defensible candidate, but current allele data are FDR-negative and underpowered.  
  Asset: `project/results/hla_deepdive_2026_05_08/track4_ht_ptc_hla2/plots/F3_HT_overlap_PTC_vs_baseline.png` (yes)
- Fig. 6: **Negative controls against overclaim** — The immune-context axis is not explained by somatic HLA LoF, selective HLA-G escape, or neoantigen load alone.  
  Asset: `project/results/hla_deepdive_2026_05_08/track27_hla_g_nonclassical/figs/F04_decoupling_vs_dm1.png` (yes)
- Fig. 7: **BRAF-rich fail-to-induce subgroup** — A DM1-high/IFN-low subgroup is BRAF-rich and biologically distinct, but not yet a clinical-outcome claim.  
  Asset: `project/results/hla_deepdive_2026_05_08/track28_immune_cold_thca/figs/F01_quadrant_scatter.png` (yes)

## 2026-05-09 추가 외부 검증 레이어

- **Paper 4:** strict Graves DPB1*05:01 proxy는 FinnGen/PanUKBB에서 p<1e-6 positive signal 1개, p<1e-4 positive signal 3개를 보이고, HT/hypothyroid branch는 genome-wide negative signal 5개다. Broad hyperthyroid endpoint의 negative signal은 phenotype heterogeneity panel로 분리한다.
- **Paper 4:** GWAS Catalog에서 GD MHC GWS share=0.128 (36/282), hypothyroidism MHC GWS share=0.305 (981/3213). AFND에서 DPB1*05:01 South Korea weighted frequency=0.367, East-Asia pool=0.375로 population baseline을 고정한다.
- **Paper 2:** GSE286332 PTC+HT vs PTC-only 재검정에서 HLA-I d=2.34, HLA-II d=3.65, TLS d=3.01이며 모두 FDR<0.002다.
- **Paper 2:** sample-unique TCGA network에서 DM1-HLA-II rho=0.400 (n=505), IFNG-HLA-II rho=0.874, TLS-HLA-II rho=0.760, RAI-HLA-II rho=-0.479다.
- **Paper 2:** scRNA Pu/Lu HLA-II cell-type concordance rho=0.95이고, spatial TLS는 HLA-II in-niche positive 16/16 samples다.
- Output: `project/reports/2026_05_09_HLA_EXTRA_EXTERNAL_VALIDATION_KR.md`

## Wave2 reviewer-proof 레이어

- **Paper 4:** allele triangulation map으로 meta OR, GD-specificity, heterogeneity, external tag-SNP direction, AFND context를 한 판에 묶었다. C*01:02는 specificity-led candidate, DPB1*05:01은 AITD-broad class-II anchor로 분리한다.
- **Paper 4:** endpoint robustness heatmap으로 strict Graves, autoimmune hyperthyroid, broad hyperthyroid, HT/hypothyroid를 분리해 FinnGen/PanUKBB direction conflict를 phenotype sensitivity로 처리한다.
- **Paper 2:** GSE286332 exact label permutation에서 HLA-II d=3.65, exact p=0.000062; AP/TLS composite d=3.08, exact p=0.000185, AUC=0.96이다.
- **Paper 2:** TCGA HLA-II dominance model은 full R2=0.83이며 immune proxy와 IFNG가 가장 큰 설명력을 가진다. 이 결과는 HLA-II가 tumor-cell-only signal이 아니라 immune/APC/TLS-rich context임을 강화한다.
- Output: `project/reports/2026_05_09_HLA_LEVEL_UP_WAVE2_KR.md`

## NComm Readiness Score

| Criterion | Paper 4 /10 | Paper 2 /10 | Upgrade Action |
|---|---:|---:|---|
| Big biological question | 8 | 8 | Frame Paper 4 as Pan-Asian AITD architecture; frame Paper 2 as immune-state ecology, not allele association. |
| Novelty over existing literature | 7 | 8 | Paper 4 needs Korean adult NGS or a stronger atlas/specificity angle; Paper 2 novelty is spatial TLS-HLA-II + IFNG mediation. |
| Independent validation | 8.5 | 8.5 | Added explicit external-validation and Wave2 robustness: Paper 4 uses GWAS Catalog, FinnGen/PanUKBB endpoint sensitivity, AFND and allele triangulation; Paper 2 uses GSE286332 exact permutation, sample-unique TCGA network, scRNA concordance, spatial TLS and ICI context. |
| Mechanistic depth | 5 | 8 | Paper 4 mechanism remains peptide/tolerance hypothesis; Paper 2 has IFNG mediation, scRNA, spatial TLS, and negative controls. |
| Reviewer safety | 8.5 | 8 | Wave2 adds endpoint-sensitivity, exact permutation/bootstrap and multivariable dominance layers. Keep HT-PTC allele claims in roadmap language. |
| NComm readiness today | 7.5 | 7.5 | External validation plus Wave2 robustness raises both to credible borderline NComm/resource-synthesis if claims stay disciplined; a decisive prospective cohort is still the clean path to 8.5+. |
| NComm readiness after one decisive dataset | 8.5 | 8.5 | Paper 4: Korean adult GD NGS. Paper 2: HT-PTC germline NGS plus HT-PTC spatial validation. |


## Reviewer Attack Matrix

| Paper | Attack | Response |
|---|---|---|
| Paper 4 | This is only a meta-analysis of old HLA papers. | Make the novelty four-layer: Korean extraction, healthy HLA atlas, cross-autoimmune specificity, trans-ancestry MHC sensitivity. Admit that de novo Korean adult GD NGS is the top-tier upgrade. |
| Paper 4 | B*46:01 is too heterogeneous to call one effect. | Do not call one universal effect. Report pooled significance with I2=85% and prediction interval crossing 1; use it as population-heterogeneous architecture. |
| Paper 4 | DPB1*05:01 is not Graves-specific. | Agree. Reframe as AITD-broad, with C*01:02 as the most GD-specific signal. |
| Paper 4 | FinnGen/PanUKBB DPB1*05:01 direction differs across thyroid endpoints. | Separate strict Graves, autoimmune hyperthyroid, broad self-reported hyperthyroid, and HT/hypothyroid endpoints. Report strict-GD positive proxy support and broad-endpoint sensitivity as phenotype heterogeneity, not as allele-level Korean replication. |
| Both | RNA-seq HLA imputation is not germline-grade. | Use RNA-seq HLA only for population-level/resource or exploratory ranking. Final validation requires germline NGS/SBT. |
| Paper 2 | HT-PTC allele data are underpowered. | Agree and move allele results to prospective candidate/roadmap. The main Paper 2 story is expression-module and spatial immune biology. |
| Paper 2 | HLA-II expression is just immune-cell contamination. | Use purity partial correlations, sample-unique TCGA network, GSE286332 HT-overlap validation, scRNA cell-source decomposition, and TLS niche-level spatial analysis. Do not claim tumor-cell-only origin. |
| Paper 2 | Mediation is not causality. | Call it observational mediation. The causal test is IFNG perturbation or cytokine stimulation in thyroid cancer/organoid models. |
| Paper 2 | HLA-G suggests immune escape. | Do not overclaim. HLA-G rises with DM1 but decoupling rho≈0; it tracks broad HLA-I induction rather than selective escape. |
| Paper 2 | Clinical outcome claims are weak. | Do not build main claims on survival/response. Keep outcome signals contextual unless externally validated. |

## Abstract급 스토리라인

### Paper 4

Autoimmune thyroid disease shows strong HLA genetic architecture, but East-Asian allele effects, disease specificity, and Korean baseline context remain fragmented across older case-control studies and reference panels. We curated Korean and Pan-Asian HLA evidence, harmonized carrier and allele-frequency metrics, built a healthy Korean/Pan-Asian HLA atlas, and compared AITD alleles across autoimmune diseases and public trans-ancestry MHC summary statistics. DPB1*05:01, B*46:01, A*02:07 and C*01:02 form the Pan-Asian Graves/AITD core. C*01:02 is the most GD-specific focus allele, whereas DPB1*05:01 is AITD-broad. The work establishes a Korean-anchored HLA architecture and defines the adult Korean GD NGS cohort needed for definitive validation.

### Paper 2

Hashimoto-overlap papillary thyroid carcinoma should not be reduced to HLA allele association from underpowered tumor cohorts. We instead model it as an antigen-presentation immune state. In TCGA-THCA, DM1/immune-active tumors show coordinated HLA-I/HLA-II module induction, partly mediated by IFN-gamma and inversely coupled to thyroid differentiation/RAI modules. Single-cell and spatial analyses localize HLA-II to APC-rich compartments and HLA-II-hot TLS niches across thyroid tumor stages. Somatic HLA-presentation mutations, selective HLA-G escape, and neoantigen burden do not explain the axis. DRB1*04:05 emerges only as a prospective HT-overlap allele candidate requiring germline NGS validation.

## 결정적 추가 데이터

1. **Paper 4:** Korean adult GD germline NGS HLA case-control, 최소 n>=250 cases / n>=500 controls.
2. **Paper 2:** PTC+HT vs PTC-only germline NGS HLA with explicit HT pathology/TPOAb/TgAb labels. DRB1*04:05 OR=1.8 검정이면 n≈633/group 필요.
3. **Paper 2:** HT-PTC Visium/Xenium/CosMx with TLS markers, HLA-II module, B/T/myeloid labels.
4. **Both:** arcasHLA vs NGS/SBT orthogonal typing subset to close reviewer attack.

## 금지선

- Paper 2에서 HLA allele을 thyroid cancer risk, prognosis, RAI response, DM1/DM2, BRAF/RAS/TERT, stage, LN, recurrence와 연결하지 않는다.
- Paper 4의 GD/AITD HLA effect size를 Korean PTC나 HT-PTC cancer risk로 옮기지 않는다.
- IFN-gamma mediation은 observational mediation으로만 쓴다.
- HLA-G는 selective escape가 아니라 broad HLA-I induction과 같이 움직인다고 쓴다.

## 생성 파일

- `project/results/hla_two_paper_synthesis_2026_05_09/ncomm_upgrade_pack/main_figure_asset_registry.tsv`
- `project/results/hla_two_paper_synthesis_2026_05_09/ncomm_upgrade_pack/ncomm_readiness_scorecard.tsv`
- `project/results/hla_two_paper_synthesis_2026_05_09/ncomm_upgrade_pack/reviewer_attack_response_matrix.tsv`
- `project/results/hla_two_paper_synthesis_2026_05_09/ncomm_upgrade_pack/title_options.tsv`
