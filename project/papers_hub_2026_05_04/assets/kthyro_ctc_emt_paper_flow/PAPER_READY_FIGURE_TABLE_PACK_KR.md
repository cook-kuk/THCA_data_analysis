# Paper-Ready Figure/Table Pack

Purpose: 대표 figure/table을 실제 manuscript/proposal에 바로 배치할 수 있게 legend, results scaffold, section map, data used, where used, claim boundary를 완성한다.

Core boundary: 현재 public-data output은 tissue-state/marker/NGS-anchor prior framework다. Yu/hospital serial CTC/cfDNA/tissue NGS data 없이 postoperative CTC-EMT transition discovery라고 쓰지 않는다.

## Figure Legends
### Grand Overview. Data-to-Claim Grand Overview
- Assets: PF15_grand_overview_data_to_claim_map.png
- Legend: Grand overview of the K-Thyro CTC-EMT-NGS proposal logic. Public tissue, spatial, single-cell, proteomic, and published CTC evidence are separated from the future prospective serial blood experiment. The figure maps how public priors support representative figures and control tables, while Yu/hospital serial CTC, cfDNA, and matched tissue-normal NGS data are required to test thyroidectomy-induced CTC-EMT state transition.
- Results scaffold: The overview fixes the current manuscript as a public-data prior/framework and the Samsung proposal as a prospective perturbation test.
- Data used: All local evidence tables; public-pilot output summaries; Yu 2024 published anchor; future Samsung/Yu cohort design.
- Where used: 교수님 첫 설명, proposal opening schematic, manuscript overview figure.
- Boundary: Schematic is not a new statistical result and does not claim postoperative CTC transition has already been discovered.

### Main Fig. 1. Study Logic And Evidence Provenance
- Assets: PF01_manuscript_flow.png;PF02_data_provenance_map.png
- Legend: Study logic and evidence provenance. Yu 2024 is used as a published feasibility anchor for serial CTC-EMT measurement in PTC, whereas TCGA, spatial, scRNA, and proteomic resources are used as public marker and tissue-state priors.
- Results scaffold: Evidence sources were assigned distinct functions so that public omics are not treated as substitutes for serial CTC data.
- Data used: Yu 2024 published CTC feasibility; TCGA/spatial/scRNA/proteomics public pilot outputs; proposal design tables.
- Where used: Introduction end, proposal Slide 3-4, reviewer response to 'what data do you have?'
- Boundary: Yu 2024 is not local raw data; public data is not used as CTC proof.

### Main Fig. 2. Driver Mutation Alone Is Not Enough
- Assets: F03_tcga_braf_state_split.png;F04_tcga_driver_variance_boundary.png
- Legend: Driver mutation alone does not define a single biological state in PTC. TCGA-THCA public-pilot outputs show that BRAF-mutant tumors distribute across multiple state/vulnerability labels, with the largest BRAF label accounting for 33.2% of the BRAF-mutant subset.
- Results scaffold: This supports phenotype-state plus genotype interpretation rather than driver-only postoperative blood-signal interpretation.
- Data used: TCGA-THCA bulk RNA/mutation/clinical public-pilot tables; BRAF-only vulnerability label distribution; axis variance models.
- Where used: Results 1, Samsung Slide 5, rationale for phenotype + genotype integration.
- Boundary: Does not prove CTC shedding, postoperative transition, or recurrence prediction.

### Main Fig. 3. Spatial Tissue-State Organization
- Assets: F05_spatial_tissue_state_organization.png;spatial_coherence_summary.png
- Legend: Spatial tissue-state organization in public thyroid spatial transcriptomics. GSE250521-derived coherence and spatial summary outputs support non-random organization of tissue states across 16 slides and 57,144 spots.
- Results scaffold: Spatial organization supports marker/ROI logic for prospective CTC-EMT studies, not direct CTC-origin proof.
- Data used: GSE250521 spatial transcriptomics; 16 slides; 57,144 spots; coherence/adjacency/hotspot/interface summaries.
- Where used: Results 2, tissue-context rationale, reviewer response to marker-origin logic.
- Boundary: Spatial organization is not direct proof of CTC origin or shedding source.

### Main Fig. 4. Marker Module Cross-Layer Support
- Assets: F06_scrna_marker_context.png;F07_proteomics_dediff_direction.png;PF10_marker_module_support_heatmap.png
- Legend: Cross-layer support for candidate CTC-EMT marker modules. Single-cell marker context, proteomic direction, and marker-module matrices support epithelial, hybrid/mesenchymal, thyroid-lineage, stress/survival, and immune/stress modules.
- Results scaffold: The marker panel is biologically organized and suitable for prospective testing but remains hypothesis-generating.
- Data used: PTC scRNA marker-context summary; bulk proteomics proxy; marker module cross-layer support matrix.
- Where used: Results 3, CTC marker-panel design, Inocras/vendor assay discussion.
- Boundary: Marker support is hypothesis-generating; not a validated clinical CTC assay.

### Main Fig. 5. Phenotype-To-Genotype NGS Bridge
- Assets: PF03_marker_to_ngs_map.png;F08_ctc_emt_ngs_prior_panel.png
- Legend: Phenotype-to-genotype bridge for CTC-EMT-NGS interpretation. CTC epithelial, hybrid E/M, and mesenchymal states are mapped to matched tissue-normal NGS, cfDNA, and CTC-enriched sequencing logic.
- Results scaffold: The proposed blood readout requires genetic anchoring because phenotype alone can be nonspecific in early PTC.
- Data used: CTC-EMT marker prior table; marker-to-NGS map; Samsung proposal NGS logic.
- Where used: Results 4, proposal Aim 2, Inocras questionnaire, assay workflow figure.
- Boundary: Does not assume cfDNA/CTC-enriched NGS will work in every early PTC patient.

### Main Fig. 6. Missing Dataset And Future Data Unlock
- Assets: PF04_dataset_gap_bridge.png;PF12_future_data_unlock_map.png;PF05_claim_boundary_matrix.png
- Legend: Missing dataset and prospective data unlock map. No public dataset jointly provides serial pre/post-thyroidectomy CTC-EMT phenotype, matched tumor-normal NGS, cfDNA or CTC-enriched NGS, and postoperative follow-up.
- Results scaffold: The data gap justifies the prospective Samsung/Yu cohort and defines required CRF and assay fields.
- Data used: serial CTC-NGS dataset gap map; future data unlock table; claim boundary matrix; hospital CRF.
- Where used: Discussion/future prospective study, Samsung final ask, Prof. Yu data request.
- Boundary: Dataset absence justifies the prospective study; it is not itself a biological discovery.

### Defense/Supp. Reviewer Defense And Execution Control
- Assets: PF07_reviewer_to_figure_map.png;PF14_analysis_control_board.png;F09_publishability_decision.png
- Legend: Reviewer defense and execution-control figure set. Expected reviewer attacks, analysis-control tables, publishability boundaries, and stage-gated kill criteria are mapped to supporting figures and tables.
- Results scaffold: The defense layer keeps the proposal ambitious but bounded and falsifiable.
- Data used: reviewer risk heatmap; statistical analysis plan; stage-gate kill criteria; safe sentence bank.
- Where used: Supplement, internal review, Samsung Q&A rehearsal.
- Boundary: Defense artifacts are not new biological results.

## Table Captions
- Data Definition Master (data_dictionary_master.tsv): Data Definition Master. 각 데이터셋이 무엇이고 무엇에 쓰였는지 정의. This table is used for Methods/Data section; reviewer data provenance. Boundary: forbidden_claim 컬럼을 반드시 같이 본다.
- Caption And Results Scaffold (figure_caption_and_results_scaffold.tsv): Caption And Results Scaffold. 각 figure의 caption, results paragraph, methods note, reviewer sentence. This table is used for Figure legend drafting; manuscript Results scaffold. Boundary: 문장 scaffold이지 final claim expansion이 아니다.
- Representative Figure Overview (representative_figure_overview.tsv): Representative Figure Overview. 대표 figure spine과 data/claim/boundary. This table is used for 교수님 briefing; manuscript figure plan. Boundary: 이 파일 자체는 guide다.
- Data-To-Figure-Claim Map (data_to_figure_claim_map.tsv): Data-To-Figure-Claim Map. 데이터별 어디에 썼고 어떤 claim을 허용하는지. This table is used for Methods, rebuttal, grant review defense. Boundary: future data rows는 current result가 아니다.
- Claim Boundary Matrix (claim_boundary_matrix.tsv): Claim Boundary Matrix. Allowed/forbidden statement를 명시. This table is used for Reviewer defense and final editing. Boundary: forbidden 문장은 본문에 쓰지 않는다.
- Hospital Data Request CRF (hospital_data_request_crf.tsv): Hospital Data Request CRF. Yu/hospital에서 받아야 할 필드. This table is used for 7-day execution, IRB/CRF. Boundary: 요청 필드이지 이미 확보한 데이터가 아니다.
- Statistical Analysis Plan (statistical_analysis_plan.tsv): Statistical Analysis Plan. 현재 public-data analysis와 future prospective test 분리. This table is used for Methods and Samsung review. Boundary: future SAP는 현재 결과가 아니다.
- Stage-Gate Kill Criteria (stage_gate_kill_criteria.tsv): Stage-Gate Kill Criteria. assay/NGS/follow-up 실패 시 downgrade rule. This table is used for Budget/risk management. Boundary: kill criteria를 숨기지 않는다.
- Safe Sentence Bank (safe_sentence_bank.tsv): Safe Sentence Bank. 본문에 바로 쓸 수 있는 안전 문장. This table is used for Writing and reviewer answer. Boundary: safe sentence 외 확장 시 claim boundary 재확인.

## Section Map
- Title / Abstract: A public-data prior and prospective Science proposal for serial PTC CTC-EMT-NGS mapping after thyroidectomy. Assets: PF15_grand_overview_data_to_claim_map.png;PF11_manuscript_blueprint_board.png Write now: We define a staged framework for testing whether thyroidectomy perturbs CTC-EMT state and whether persistent postoperative genetic signals can be anchored to the resected tumor clone. Do not write: Do not say the postoperative CTC-EMT transition has already been discovered in our cohort.
- Introduction: PTC is usually curable, but postoperative residual-risk biology is not directly read in real time. Assets: PF01_manuscript_flow.png;slide01_problem.png;slide02_hypothesis.png Write now: Thyroidectomy provides a biologically interpretable perturbation around which serial blood and matched tissue genotyping can be organized. Do not write: Do not frame this as a broad AI/spatial/drug platform.
- Results 1: Driver mutation alone is insufficient to define PTC biological state. Assets: F03_tcga_braf_state_split.png;F04_tcga_driver_variance_boundary.png Write now: BRAF-mutant PTC tumors distributed across multiple state labels, supporting the need for phenotype-state interpretation beyond driver genotype. Do not write: Do not infer CTC shedding or recurrence from TCGA tissue-state heterogeneity.
- Results 2: Public spatial data support organized tissue-state context. Assets: F05_spatial_tissue_state_organization.png;spatial_coherence_summary.png Write now: Spatial tissue-state organization supports marker and ROI logic for prospective CTC-EMT studies. Do not write: Do not claim spatial hotspots are direct CTC shedding sources.
- Results 3: Marker modules have cross-layer tissue/cell/protein support. Assets: F06_scrna_marker_context.png;F07_proteomics_dediff_direction.png;PF10_marker_module_support_heatmap.png Write now: Epithelial, hybrid/mesenchymal, thyroid-lineage, stress/survival, and immune/stress modules are interpretable across public evidence layers. Do not write: Do not present the marker panel as a validated clinical CTC assay.
- Results 4: Phenotype must be genetically anchored by matched tumor-normal NGS and cfDNA/CTC-enriched sequencing. Assets: PF03_marker_to_ngs_map.png;F08_ctc_emt_ngs_prior_panel.png Write now: The proposed CTC-EMT readout is only interpretable as residual-risk biology when anchored to tissue-derived genetic context. Do not write: Do not assume low-input CTC sequencing succeeds in every early PTC case.
- Results 5 / Discussion: The exact public dataset needed for serial postoperative CTC-EMT genetic mapping is missing. Assets: PF04_dataset_gap_bridge.png;PF12_future_data_unlock_map.png;PF05_claim_boundary_matrix.png Write now: The missing joint dataset defines the prospective Samsung/Yu cohort requirement and its stage-gated endpoints. Do not write: Do not use dataset absence as a substitute for prospective biological evidence.
- Supplement / Reviewer Defense: Reviewer attack, claim boundary, and kill criteria are explicitly pre-specified. Assets: PF07_reviewer_to_figure_map.png;PF14_analysis_control_board.png;F09_publishability_decision.png Write now: The package is designed to be ambitious but falsifiable, with downgrade criteria for assay or data failure. Do not write: Do not hide feasibility risks or overstate current public-data claims.