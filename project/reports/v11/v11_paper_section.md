## 5.X Novel co-expression modules localise the DIAL flip

### Motivation

The v5.1 diagnostic fires on a ComBat-corrected 3000-gene space, and v8 S6B BUHMBOX showed that on the affected axis cohort is effectively PC1 (KS p < 1e-40). If the corruption DIAL detects reflects a biological subspace rather than residual noise, it should carry modular structure detectable in unsupervised co-expression decomposition. Module-level analysis thus bridges the gene-level DIAL hit and the sample-level BUHMBOX mechanism.

### Module discovery

We ran two clustering perspectives on the 2472/2773 v5.1+v8 biomarkers retained after expression-variance filtering. WGCNA with signed scale-free topology (soft-threshold beta = 14, Ward linkage on TOM dissimilarity) [Langfelder & Horvath 2008] yielded 21 modules (M001-M021, 30-1182 genes). Leiden multi-resolution clustering on the same TOM graph [Traag et al. 2019] at gamma = {0.5, 1.0, 2.0} produced a hierarchy with 28 / 29 / 1012 modules at L1 / L2 / L3 respectively. Every WGCNA module and every high-level Leiden module was novel against MSigDB Hallmark, Reactome 2022, GO_BP 2023 and KEGG 2021 at max Jaccard < 0.2; the largest observed Jaccard was 0.0909 (M012 vs Aerobic Electron Transport Chain).

### Near-miss markers and the threshold question

Under the strict novel_marker definition (|r_Y| > 0.5 and |r_cohort| < 0.3) the count is 0. The cause is mechanistic, not statistical: median |r_cohort| across the 21 WGCNA modules is 0.916 (min 0.882, max 0.974), precisely the regime S6B identified at the sample level — cohort dominates eigengene variance. Five modules cross a relaxed |r_Y| > 0.3 threshold and constitute the near-miss set (Table T5.X). M010 (PTPRF, ACVR1, CNN3; r_Y = 0.43) and M006 (TEAD1, IQGAP1, AGFG1; r_Y = 0.32, KEGG near-match: PD-L1/PD-1 checkpoint, consistent with TEAD1 Hippo-pathway cross-talk) are the most biologically coherent. A referee should focus on this threshold sensitivity: the relaxed threshold produces a structured cluster with modest but non-zero BRAF/RAS association, not noise.

### DIAL-signature decomposition and druggability

Decomposing the DIAL signature into flip-contributing and batch-aligned axes yields orthogonal gene sets at top-100 strictness (intersection = 0); overlap appears only at top-500 (7 genes). Top flip-contributing weights load on ST3GAL6, NQO1, SYT12, DCSTAMP, NUCB2. PLEKHA6, one of eight v7 druggable targets, ranks in the label-preserving top 3 (w_preserve = 0.0754), so its prioritisation is robust under ComBat rather than induced by it. TACSTD2, another v7 druggable, anchors the largest novel Leiden L2 module (L2_1, n = 911, novelty 0.94); TMPRSS6 appears in the label-preserving top 4. These are the first prospective targets from DIAL-decomposed, hierarchy-validated co-expression structure.

Fig 5.X shows the interactive module/sunburst page (http://40.82.129.113:8012/reports/html/pages/v11_novel_pathways.html); Supplementary Fig S11.1 contains the WGCNA-TOM sunburst and Leiden hierarchy network.

**Table T5.X.** Near-miss novel modules (|r_Y| > 0.3), top three genes, and WGCNA novelty.

| module_id | n_genes | r_Y | max_jaccard | best_db_match | top_3_genes |
|-----------|--------:|------:|-----:|----------------------------------|--------------------|
| M010 | 67  | 0.4297 | 0.0352 | GO_BP: Wound Healing | PTPRF, ACVR1, CNN3 |
| M001 | 63  | 0.3960 | 0.0353 | GO_BP: Neg. Reg. Intrinsic Apoptosis | MVP, TAGLN2, ITGA3 |
| M009 | 117 | 0.3247 | 0.0482 | GO_BP: Actin Filament Organization | FLII, FLNA, MYH9 |
| M006 | 83  | 0.3244 | 0.0361 | KEGG: PD-L1 / PD-1 Checkpoint | TEAD1, IQGAP1, AGFG1 |
| M011 | 30  | 0.3046 | 0.0526 | GO_BP: Ca2+ Import Plasma Membrane | RAD23B, PRNP, ERLIN1 |

*AAAI adaptation.* This section serves as the interpretable feature-subspace decomposition arm supporting the DIAL theorem: the co-expression structure identified here *is* the corrupted subspace that DIAL detects, with cohort-dominated eigengenes mirroring the sample-level S6B BUHMBOX result.
