# Paper Repo Map v2 (2026-05-10)

Local-only build. No /var/www deployment was performed.

- Branches: 11
- Hub pages scanned: 179
- Results files scanned: 5317
- Software/method contracts: 8

## Branches

| priority | branch | status | role | data | software | boundary |
|---:|---|---|---|---|---|---|
| 1 | paper1-main | SUBMIT_NOW | trunk | TCGA-THCA RNA/HM450/clinical/SV, Korean K2, Lee 2024, GSE76039, Lu 2023, Pu 2021, DepMap/PRISM reserve. | Python, pandas, scipy, statsmodels, scikit-learn, NuSVR, matplotlib, static HTML dossier builder. | No direct mechanism-of-silencing proof; no clinical diagnostic claim without prospective validation. |
| 2 | paper11-pancancer | NEXT_MAJOR_AFTER_P1 | release branch | TCGA pan-cancer, DepMap, PRISM, hallmark scores, survival meta, Korean multi-cohort support. | Pan-cancer module scoring, survival meta-analysis, lineage mitigation, PRISM/MAPK enrichment. | Must present lineage-specificity and residual confounding up front. |
| 3 | paper1b-braf-axis | MERGE_OR_STANDALONE_DECISION | feature branch | TCGA BRAF/RAS/TERT, Lee 2024, HM450, spatial BRAF inference, single-cell support, ICI panels. | Python, survival models, Fisher tests, deconvolution, spatial scoring, HTML synthesis. | Do not let BRAF-axis dilute the main DM1 novelty unless it resolves a reviewer attack. |
| 4 | paper2-wsi-image | VIABLE_PILOT | application branch | TCGA WSI embedded slides, GSE230424, GSE250521, GSE248205, pathology tile atlases and spatial overlays. | DINOv2/UNI embeddings, CLAM, scikit-learn, spatial module scoring, RunPod GPU jobs. | No clinical deployment, no universal pathology biomarker claim, no leakage-prone slide-level inflation. |
| 5 | paper9-ic50-perturbation | WETLAB_REQUIRED | therapeutic branch | PRISM drug response, DepMap CRISPR, thyroid cell line annotations, spatial drug target overlays. | PRISM differential ranking, DepMap dependency tests, target-state overlap matrices, spatial overlay scoring. | No synthetic-lethality or patient-treatment claim without wet-lab IC50/viability validation. |
| 6 | neoantigen-clean-neo | SEPARATE_UNIVERSE | separate repository | TESLA, CEDAR, dbPepNeo2, NeoDB, McPAS, IEDB, MD/TCR/pMHC audit outputs. | CLEAN-Neo pipeline, QK/rule-gated models, ESM2 features, external predictors, PMTnet, OpenMM/MD audit. | No THCA vaccine-readiness claim; public predictors are comparator-only when licenses/data leakage are uncertain. |
| 7 | paper2b-ht-immune | BIOLOGY_STRONG_ALLELE_WEAK | sidecar branch | TCGA HT labels, GSE286332, GSE248205, GSE250521, HLA deep-dive tracks. | arcasHLA outputs, Fisher/meta-analysis, spatial burden summaries, module scoring. | No Korean/allele-level HLA claim until adult Korean GD/PTC NGS exists. |
| 8 | paper3-ici-readiness | FROZEN_DESIGN | reserve branch | Public ICI cohorts, ATC/THCA immune modules, HLA/IFNG tracks, phase2/phase4 ICI summaries. | Cohort pooling, module scoring, response association tests, forest summaries. | No thyroid ICI response model without matched thyroid ICI outcomes. |
| 9 | paper4-gd-hla | UPGRADE_IF_DATA_ARRIVES | data-gated branch | Korean PTC pool, Chen 2018 Han Chinese GD, GSE286332 arcasHLA, HLA deep-dive tracks. | arcasHLA, allele-frequency harmonization, Fisher tests, meta-analysis forest plots. | No Korean GD association claim from proxy PTC/HT data. |
| 10 | paper10-12-network-atlas | ABSORB_INTO_P11 | atlas branch | paper12 network edges/modules, cancer gene matrix, HLA atlas, source atlas, figure/table atlases. | Correlation networks, module summaries, static atlas pages, visualization tables. | Do not split into low-novelty standalone paper unless it gains external validation. |
| 11 | paper13-spatial-drug | VISUAL_STRONG_DATA_NEEDS_LOCK | visual application branch | GSE250521/GSE230424 spatial overlays, target gene expression, drug-state scorecards, CV2 validation ROI packs. | Spatial module scoring, target-state overlap, ROI gallery generation, static visual dossiers. | Spatial co-localization is not drug efficacy; match/mismatch boards must stay explicit. |

## History / Dependency Edges

| parent | child | reason |
|---|---|---|
| paper1-main | paper1b-braf-axis | shared TCGA/DM axis; BRAF branch tests whether driver subtype explains or complements DM1 |
| paper1-main | paper9-ic50-perturbation | DM1/MAPK state seeds perturbation and PRISM/DepMap prioritization |
| paper1-main | paper11-pancancer | DM1 trunk becomes pan-cancer portability hypothesis |
| paper1-main | paper2-wsi-image | molecular DM1 labels seed image/pathology predictors |
| paper2-wsi-image | paper2b-ht-immune | spatial/pathology immune context branches into HT/HLA questions |
| paper2b-ht-immune | paper4-gd-hla | HT/HLA context defines Korean GD data gate |
| paper9-ic50-perturbation | paper13-spatial-drug | drug-response proxies are localized with target-state spatial overlays |
| paper9-ic50-perturbation | paper11-pancancer | PRISM/DepMap pan-cancer signal informs Paper 11 therapeutic appendix |
| paper11-pancancer | paper10-12-network-atlas | network/atlas layer supports portability interpretation |

## Software Interpretation Ledger

| method | software | inputs | outputs | interpretation | failure mode |
|---|---|---|---|---|---|
| 8-gene DM1/module scoring | Python, pandas, scipy, statsmodels | RNA expression matrices, cohort metadata, driver annotations | DM1 score, cohort rank, contrast tables, survival/driver overlays | Relative molecular-state axis; strongest when replicated across cohorts and split logic. | Batch, lineage and label leakage can make clean biology look stronger than it is. |
| Methylation-expression support | pandas, scipy correlation, statsmodels | TCGA HM450 probes, gene expression, DM calls | gene-level methylation-expression correlation and group contrasts | Supports epigenetic compatibility, not causal silencing. | Probe mapping and tumor purity can confound direction. |
| Deconvolution / cell-state attribution | NuSVR-style deconvolution, marker panels, custom Python | bulk RNA, marker signatures, spatial/scRNA references | cell-state fractions, marker module overlays, MAPK/TDS panel relationships | Explains which compartment might carry a signal. | Reference mismatch and collinearity; fractions are estimates, not sorted-cell truth. |
| Pathology foundation modeling | DINOv2/UNI, CLAM, scikit-learn, RunPod GPU | H&E whole-slide tiles, slide labels, spatial labels for bridge analyses | slide/tile embeddings, AUC, attention maps, pathology-spatial overlays | Image feasibility and localization signal. | Slide/site leakage, small N, inflated AUC without leave-one-site/sample tests. |
| Spatial transcriptomics overlay | scanpy-like matrices, custom scoring, matplotlib | Visium/spatial spot matrices, histology coordinates, module genes | spot maps, hotspot burden, target-state overlap, stage trend tables | Shows tissue localization and co-occurrence. | Spot resolution mixes cells; spatial co-occurrence is not causal interaction. |
| PRISM / DepMap IC50 proxy mapping | Python, pandas, scipy, multiple-testing correction | PRISM drug response, DepMap CRISPR, model lineage metadata, DM1 scores | drug rankings, MAPK enrichment, dependency tables, spatial target overlays | Hypothesis-prioritization map for wet-lab IC50/viability testing. | Cell-line N and lineage confounding; proxy response is not clinical efficacy. |
| Neoantigen leakage-aware benchmarking | CLEAN-Neo, rule-gated/QK models, external predictors, ESM2, OpenMM/MD audit | TESLA/CEDAR/dbPepNeo2/NeoDB/McPAS/IEDB, TCR/pMHC structures, public predictors | locked split metrics, abstention queues, candidate scorecards, MD evidence cards | Methods and candidate-prioritization framework. | Source overlap, allele leakage, license constraints and incomplete immunogenicity labels. |
| Static paper hub / dossier generation | Python HTML builders, CSS, TSV/JSON inventories | local result files, figure assets, paper metadata | paper1-style HTML dossiers, page inventories, figure/table/claim ledgers | Review and navigation layer, not primary analysis. | Can overstate if not tied to source paths and claim boundaries. |

## Current Paper Count Decision

Broad count: 11 active/near-active branches plus the Neoantigen separate repository branch.
Immediate manuscript lane: Paper 1 trunk. Next high-value lane: Paper 11. Wet-lab priority lane: Paper 9 IC50/PRISM/DepMap. Visual-pilot lane: Paper 2/13. Data-gated lane: Paper 4 GD-HLA. Separate methods lane: CLEAN-Neo/CROSS-Neo.
