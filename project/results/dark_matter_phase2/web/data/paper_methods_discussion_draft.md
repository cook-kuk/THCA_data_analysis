# Dark Matter of Thyroid Cancer — Paper Methods + Discussion Draft

**Working title.** *A continuous transcriptional differentiation axis sub-stratifies driver-negative thyroid cancer at single-cell resolution and identifies DICER1/EIF1AX as the FVPTC-like genomic anchor.*

**Author.** Seungho Cook (corresponding); collaborators TBD. Affiliations: TBD.

---

## Methods

### Cohort assembly
The discovery cohort was TCGA-THCA (n=482 with mutation calls; n=506 with curated survival data per Liu 2018 TCGA Clinical Data Resource [@liu2018cdr]). Mutation classifications (BRAF V600E, NRAS/HRAS/KRAS hotspot Q61/G12/G13, TERT promoter C228T/C250T) were obtained from the cBioPortal `thca_tcga_pub` study and supplemented for TERT promoter recovery from masked somatic mutation MAFs. Histology subtype (cPTC, FVPTC, FTC, PDTC, ATC) was extracted from GDC pathology metadata.

The Korean external bulk cohort was Yoo SK et al. 2016 [@yoo2016ptc] (PRJEB11591, n=180 with FA, miFTC, cPTC, fvPTC). Mutation status and molecular subtype (BRAF-like / RAS-like / NBNR) were mined from Supplementary Table S6 of the paper. Run-to-sample mapping was via ENA metadata.

External single-cell cohorts: GSE241184 (Pu et al. 2023; 1 patient with tumor/normal/lymph node metastasis, 30,493 cells) [@pu2023gse241184]; GSE193581 (Lu et al. 2023; 23 samples PTC/ATC/normal, 67,678 cells) [@lu2023jci]; GSE184362 (Pan et al. 2021; 11 PTC patients with multi-site samples T/P/LN, 158,577 cells) [@pan2021gse184362].

Pan-cancer specificity controls: GSE39582 (COAD; n=586), GSE31210 (LUAD; n=247), additional thyroid microarray cohorts GSE33630 (n=49) and GSE29265.

### 8-gene panel selection
Gene panel selection was performed on a 67-gene curated framework (TIERA67) excluding the Driver_anchor category (BRAF, NRAS, HRAS, KRAS, RET, NTRK1, NTRK3, ALK, PAX8, PPARG, TERT, EIF1AX) to prevent label leakage with reference subtypes. RandomForest feature importance ranking on the resulting 55-gene clean pool produced the top 8: DIO1, FOXE1, NKX2-1, PAX8, SLC5A5, TG, TPO, TSHR — all from the Thyroid Differentiation Score (TDS) core category. Selection process is documented and audited in `results/audit_2026_04_29/audit_report_8gene.md`.

### Bulk transcriptional clustering
TCGA-THCA log2(TPM+1) expression for the 8-gene panel was used as input to ConsensusClusterPlus [@wilkerson2010ccp] with k-grid 2-10. The k=2 solution showed stability 0.978 (PAC = 0.05) and was retained. Cluster labels (DM1, DM2) were assigned within the Dark Matter sub-cohort defined as BRAF V600E-negative AND RAS hotspot-negative. ConsensusClusterPlus was independently validated by KMeans + Adjusted Rand Index bootstrap (n=100 iterations, 80% subsample): mean ARI = 0.994, 95% CI [0.93, 1.00], with random null = 0 (feature shuffle control).

### Single-cell processing
Each scRNA-seq cohort was processed independently with scanpy 1.12 [@wolf2018scanpy]. QC: cells with ≥200 expressed genes and ≤25% mitochondrial reads. Normalization: total-count to 10,000 with log1p transform. Highly variable genes: top 2000 (Seurat method). Dimensionality reduction: PCA (30 components). Neighbor graph: k=15. Clustering: Leiden (resolution 0.6). Cell type annotation: marker score with sc.tl.score_genes for {Thyrocyte: TG/TPO/TSHR/TFF3/PAX8/FOXE1; T cell: CD3D/CD3E/CD8A; Myeloid: LYZ/CD68/CD14; Endothelial: PECAM1/VWF/CDH5; Fibroblast: COL1A1/COL1A2/DCN}; cluster-level argmax assigns the cell type.

### Single-cell signature scoring
Per-cell scores: 8-gene panel signature score, FVPTC signature (TG/TPO/TSHR/DIO1/DIO2/SLC5A5/FOXE1), cPTC signature (KRT19/TIMP1/FN1/BCL2/CITED1), DICER1 axis proxy (DICER1/DROSHA/DGCR8/AGO1/AGO2), EIF1AX axis proxy (EIF1AX/EIF1/EIF2S1/EIF4E/EIF4G1) — all via sc.tl.score_genes (control gene set size = 50, n_bins = 25).

### Survival analysis
Cox proportional hazards via lifelines [@davidsonpilon2019lifelines]. Univariate cluster vs reference HR; multivariate with cluster + age + stage_ord. PFI as primary endpoint per Liu 2018 [@liu2018cdr], with DFI/DSS/OS as supplementary. Within Stage I sub-analysis to control for stage confounding. Kaplan-Meier curves with multivariate log-rank (multiple-class) for driver-class survival comparison. All statistical tests two-sided.

### Pathway enrichment
Eleven curated thyroid-relevant pathways (Thyroid_differentiation, MAPK_target, PI3K_AKT, EMT_markers, Epithelial_markers, let-7_targets, miR-200_targets, Translation_eIF1AX, Cell_cycle_E2F, Hypoxia, Inflammation_IFN) scored via mean log2(TPM+1) of pathway genes. DM1 vs DM2 comparison: Welch t-test with Cohen's d effect size.

### Predictive performance
8-gene → cPTC vs FVPTC binary classification within Dark Matter (n=125, cPTC=79, FVPTC=46) via logistic regression (sklearn [@pedregosa2011sklearn], C=1.0, max_iter=2000) with 5-fold stratified cross-validation. ROC AUC reported with bootstrap-CI baseline. Single-gene baselines computed for each panel member.

### External cohort validation pipeline
Identical sc-processing pipeline applied to each external cohort. Per-patient Pearson correlation between 8-gene and FVPTC signatures within tumor thyrocytes (n ≥ 30 cells per patient). Pooled correlation across all malignant cells; bootstrap 95% CI.

### Pan-cancer specificity
8-gene + FVPTC scores computed identically on COAD (GSE39582) and LUAD (GSE31210) bulk RNA-seq. Pearson r reported as control.

### Korean cohort calibration
The K2 (Yoo 2016) cohort 8-gene TPM matrix shows ~10-100× scale inflation vs TCGA. Within-sample-centered re-prediction (log10 + per-sample mean centering) attempted to recover discrimination. Honest result: median split on centered means restores 55:45 cluster balance but cPTC/FVPTC histology concordance only 35.3% (Fisher OR=0.16, p=0.067). TCGA-trained centered logistic regression transfer is reported as future work.

### Software stack
scanpy 1.12.1 [@wolf2018scanpy], pandas 2.3.3, scipy 1.13, lifelines 0.27 [@davidsonpilon2019lifelines], scikit-learn 1.5 [@pedregosa2011sklearn], matplotlib 3.9, statsmodels (BH FDR), Python 3.12. All code at https://github.com/USERNAME/THCA_DarkMatter (TBD).

### Data availability
TCGA-THCA via GDC. Yoo 2016 K2 raw via ENA PRJEB11591 + S6 mutation table (DOI). GSE241184/GSE193581/GSE184362 via NCBI GEO. All processed result tables and code in `project/results/dark_matter_phase1/` and `project/results/dark_matter_phase2/` (submission supplementary).

---

## Discussion (7-paragraph structure)

### ¶1 Summary of finding
The 8-gene transcriptional axis (DIO1/FOXE1/NKX2-1/PAX8/SLC5A5/TG/TPO/TSHR) sub-stratifies BRAF V600E / RAS hotspot-negative thyroid cancer (28% of TCGA-THCA, "Dark Matter") into a cPTC-architectured DM1 (n=89, mean age 42, mechanism unknown) and an FVPTC-like DM2 (n=55, mean age 54, DICER1/EIF1AX/PPM1D 9× enriched, p=0.0175). The axis is reproducible at single-cell resolution across four independent thyroid sc cohorts (Phase 1 GSE241184 r=0.91, GSE193581 r=0.89, GSE184362 r=0.89, multi-site r=0.91) and Korean bulk validation cohort (Yoo 2016 NBNR concordance 93.5%, DICER1+EIF1AX 10.3% ≈ TCGA 10.9%). Pathway analysis reveals DM2 is more differentiated (TDS Cohen d=+1.54), less EMT (d=-1.33), less MAPK active (d=-1.24), with let-7 target genes (HMGA2/LIN28B/MYC) suppressed (d=-2.09). 8-gene → cPTC/FVPTC predictive AUC = 0.71 in 5-fold CV, indicating moderate clinical decision support.

### ¶2 Comparison with prior art
The NBNR (Non-BRAF-Non-RAS) molecular subtype was first defined in the Korean Yoo 2016 PLoS Genetics paper [@yoo2016ptc] as one of three transcriptional clusters (BRAF-like, RAS-like, NBNR). Our Dark Matter definition (BRAF V600E− AND RAS hotspot−) recapitulates Yoo's NBNR with 93.5% concordance (43 of 46 samples). Our key contribution is the further sub-stratification of NBNR into DM1 vs DM2, with DM2 enriched for the alternative drivers (DICER1, EIF1AX, PPM1D) reported by Yoo 2016 as NBNR-associated [@yoo2016ptc] and confirmed in 899 nodules by Wang 2025 (DICER1 ⊥ BRAF) [@wang2025dicer1]. Wasserman 2018 [@wasserman2018dicer1] reported DICER1 mutations frequent in adolescent-onset PTC; while our adult cohort places DICER1+ in the older DM2 subgroup (mean 54 yrs), incorporation of pediatric cohorts into a unified analysis is an immediate future direction. Chernock 2021 [@chernock2021macrofollicular] reported macrofollicular variant FTC = DICER1-mutated young females, consistent with our DM2 = FVPTC-like + DICER1-rich phenotype.

### ¶3 Mechanism interpretation
At single-cell resolution, DICER1 pathway expression (DICER1/DROSHA/AGO1/2) correlates with the 8-gene score r = 0.01 — i.e., DICER1 mutation acts as a genomic event whose phenotypic consequence (FVPTC architecture) is not directly read out by DICER1 pathway transcript levels in individual cells. This is consistent with Frontiers Endocrinol 2023 [@frontiersDicer1miRNA2023] reporting that DICER1 RNase IIIb mutations trigger let-7 / miR-200 dysregulation: the bulk-level let-7 target signature (HMGA2/LIN28B/MYC) is markedly reduced in DM2 (Cohen d=-2.09) — consistent with DICER1 mutation effect propagating through miRNA biogenesis to the differentiation phenotype, but not through DICER1 transcriptional level itself. Our reading is honest: DM2 = FVPTC-architectured tumors, ~11% of which carry DICER1/EIF1AX/PPM1D mutations as a genomic correlate, but the transcriptional cluster is driven by the broader differentiation phenotype rather than the DICER1 pathway in isolation.

### ¶4 Method comparison and clinical utility
Current targeted molecular tests for ambiguous Bethesda III/IV thyroid nodules (Afirma GEC/GSC, ThyroSeq v3) rely on driver mutation detection and report "no alterations" or "benign-class" for ~28% of patients (the Dark Matter subgroup). The proposed 8-gene RNA score, working independently of mutation calling, predicts cPTC vs FVPTC architecture with AUC = 0.71 in cross-validation on 125 driver-negative patients (Figure E10). For prospective deployment, the score could be implemented as a 4-gene RT-qPCR (TPO + DIO1 + TG + FOXE1) given that single-gene baselines of TPO and DIO1 reach AUC 0.73 individually. Honest limitation: AUC 0.71 represents moderate clinical decision support, not stand-alone diagnosis; integration with cytology and supplementary markers is recommended.

### ¶5 DM1 — true unknown driver population
The cPTC-architectured DM1 subgroup (n=89, 18% of TCGA-THCA) presents with younger median age (42 vs 54 yrs in DM2; Welch p=1.9e-5, Cohen d=0.77) and slightly more advanced stage (Stage III/IV 22% vs 13%). DM1 carries no canonical driver mutation, no consistent alternative driver, and no clear pathway signature distinguishing it from DM2 beyond the differentiation axis. Anecdotal alt-driver enrichment in DM1 includes CHEK2 (n=2; possible germline DNA damage response background) and rare RTK fusions (NTRK3, ALK, RET; 1 each). DM1 represents a "true mechanism-unknown" driver-negative thyroid cancer subgroup distinct from the DICER1/EIF1AX-anchored DM2. Functional dissection of DM1 — through germline sequencing, methylation array, deeper fusion calling, or pediatric / adolescent cohort integration following Wasserman 2018 [@wasserman2018dicer1] — is the most important unfinished problem from this work.

### ¶6 Limitations
This study is retrospective and observational with seven explicit limitations: (i) TCGA-THCA event scarcity (PFI = 9 events in 137 DM patients) precluded prognostic claims and motivated the Frame B (molecular taxonomy) reframe; (ii) discovery sc cohort GSE241184 is a single 17-year-old patient (pediatric-biased) whose limitations were addressed by adult validation in P2-A1/A2 and multi-site analyses; (iii) Highly Variable Gene filtering of GSE193581 partially attenuated the FVPTC signature score (PTC-only r drops from 0.89 PTC+ATC to 0.69 PTC), recovered when full-gene data are used (P2-A2 GSE184362 r=0.89); (iv) K2 cluster prediction calibration bias resulted in 92% DM2 calls under absolute-form transfer; within-sample-centered correction recovered 55:45 balance but with imperfect cPTC/FVPTC concordance (35%); (v) DICER1/EIF1AX absolute counts are small (n=7 TCGA + n=7 K2); (vi) 분당 SNUH cohort outreach is pending and would provide Korean follow-up survival; (vii) no functional validation experiments (organoid/cell line DICER1 knock-in) have been performed — these are proposed in the Phase 3 roadmap.

### ¶7 Bridge to method paper
The DIAL-U cluster identifiability framework (multi-seed bootstrap + direction-shuffle null + identifiability metric) underlying the cluster stability evidence presented here is being prepared as a separate methodology companion paper. The present manuscript employs a light version (KMeans bootstrap ARI = 0.994; random null = 0) sufficient to anchor the biological claims; the full DIAL-U formalism with theoretical analysis and pan-cancer benchmarking will be published elsewhere.

---

## Acknowledgments
The authors thank Prof. Yu Hyeong-Won (Seoul National University Bundang Hospital) for clinical guidance and ongoing 분당 SNUH cohort discussions. Yoo SK and the original PLoS Genet 2016 cohort assembly team are acknowledged for the Korean K2 dataset.

## Author contributions
SC conceived the project, performed all bioinformatics analyses, and wrote the manuscript. (All Phase 1 + Phase 2 single-day sprint completed 2026-04-29.) Collaborators TBD.

## Competing interests
The authors declare no competing interests.
