# Representative Figure/Table/Data Overview

Purpose: 대표 figure, 대표 table, 전체 overview 그림, figure 흐름 설명, 데이터 사용 위치, claim boundary를 한 문서에 고정한다.

Core rule: public data는 prior/framework이고, thyroidectomy-induced CTC-EMT transition은 Yu/hospital serial blood + matched NGS로 prospective test한다.

## Representative Figures
### Grand Overview. Data-to-Claim Grand Overview
- Asset: PF15_grand_overview_data_to_claim_map.png
- Role: 전체 흐름 그림. clinical perturbation, public-data priors, representative figures, control tables, current/future claims를 한 장에 연결한다.
- Data used: All local evidence tables; public-pilot output summaries; Yu 2024 published anchor; future Samsung/Yu cohort design.
- Where used: 교수님 첫 설명, proposal opening schematic, manuscript overview figure.
- Claim: This is a staged Science case: public data builds priors; prospective serial CTC/cfDNA/tissue NGS tests the thyroidectomy perturbation question.
- Boundary: Schematic is not a new statistical result and does not claim postoperative CTC transition has already been discovered.

### Main Fig. 1. Study Logic And Evidence Provenance
- Asset: PF01_manuscript_flow.png;PF02_data_provenance_map.png
- Role: Yu 2024, public-data priors, and future hospital cohort를 구분한다.
- Data used: Yu 2024 published CTC feasibility; TCGA/spatial/scRNA/proteomics public pilot outputs; proposal design tables.
- Where used: Introduction end, proposal Slide 3-4, reviewer response to 'what data do you have?'
- Claim: Evidence layers have distinct roles: public priors, published feasibility, future prospective test.
- Boundary: Yu 2024 is not local raw data; public data is not used as CTC proof.

### Main Fig. 2. Driver Mutation Alone Is Not Enough
- Asset: F03_tcga_braf_state_split.png;F04_tcga_driver_variance_boundary.png
- Role: BRAF-mutant PTC도 하나의 biological state로 접히지 않음을 보여준다.
- Data used: TCGA-THCA bulk RNA/mutation/clinical public-pilot tables; BRAF-only vulnerability label distribution; axis variance models.
- Where used: Results 1, Samsung Slide 5, rationale for phenotype + genotype integration.
- Claim: Public PTC tumor states split beyond mutation group; genotype alone is insufficient as a state descriptor.
- Boundary: Does not prove CTC shedding, postoperative transition, or recurrence prediction.

### Main Fig. 3. Spatial Tissue-State Organization
- Asset: F05_spatial_tissue_state_organization.png;spatial_coherence_summary.png
- Role: 조직 state가 random이 아니라 spatially organized되어 marker/ROI logic을 지지한다.
- Data used: GSE250521 spatial transcriptomics; 16 slides; 57,144 spots; coherence/adjacency/hotspot/interface summaries.
- Where used: Results 2, tissue-context rationale, reviewer response to marker-origin logic.
- Claim: Spatial data support organized tissue-state context for marker selection.
- Boundary: Spatial organization is not direct proof of CTC origin or shedding source.

### Main Fig. 4. Marker Module Cross-Layer Support
- Asset: F06_scrna_marker_context.png;F07_proteomics_dediff_direction.png
- Role: Epithelial, EM/M, thyroid-lineage, survival/stress marker modules를 tissue/cell/protein context에서 정리한다. PF10은 별도 supplementary planning matrix.
- Data used: PTC scRNA marker-context summary; bulk proteomics proxy. (PF10 supplementary planning matrix은 별도 슬롯.)
- Where used: Results 3, CTC marker-panel design, Inocras/vendor assay discussion.
- Claim: Candidate CTC-EMT marker modules are biologically interpretable across public layers.
- Boundary: Marker support is hypothesis-generating; not a validated clinical CTC assay.

### Supplementary Planning Matrix. PF10 - Marker Module Planning Support Scores (hand-coded; not data-derived)
- Asset: PF10_marker_module_support_heatmap.png
- Role: Marker module별 prior 구성 유용성을 0/1/2로 손코딩한 manuscript-planning matrix. 분석 산출물이 아닌 기획 점수표다.
- Data used: marker_module_cross_layer_support_matrix.tsv (hand-coded). TCGA-THCA tissue (no CTC); GSE250521 spatial; scRNA module summary; bulk proteomics proxy; Yu 2024 published anchor summary (no unpublished/raw Yu CTC data).
- Where used: Supplement only. Reviewer-visible planning artifact.
- Claim: 수기로 매긴 manuscript-planning judgement이다.
- Boundary: 측정값/validation metric/cross-layer 분석 결과가 아니다. TCGA는 CTC 데이터를 포함하지 않으며 Yu 2024 컬럼은 published summary만 인용한 것이다.

### Main Fig. 5. Phenotype-To-Genotype NGS Bridge
- Asset: PF03_marker_to_ngs_map.png;F08_ctc_emt_ngs_prior_panel.png
- Role: CTC phenotype만으로 부족하므로 matched tumor-normal NGS/cfDNA anchoring이 필요함을 보여준다.
- Data used: CTC-EMT marker prior table; marker-to-NGS map; Samsung proposal NGS logic.
- Where used: Results 4, proposal Aim 2, Inocras questionnaire, assay workflow figure.
- Claim: CTC-EMT phenotype should be interpreted with genetic anchoring.
- Boundary: Does not assume cfDNA/CTC-enriched NGS will work in every early PTC patient.

### Main Fig. 6. Missing Dataset And Future Data Unlock
- Asset: PF04_dataset_gap_bridge.png;PF12_future_data_unlock_map.png;PF05_claim_boundary_matrix.png
- Role: 왜 병원/Yu serial CTC/cfDNA/tissue NGS 데이터가 필요한지 figure-driven으로 설명한다.
- Data used: serial CTC-NGS dataset gap map; future data unlock table; claim boundary matrix; hospital CRF.
- Where used: Discussion/future prospective study, Samsung final ask, Prof. Yu data request.
- Claim: No public dataset currently tests the full serial postoperative CTC-EMT genetic-map question.
- Boundary: Dataset absence justifies the prospective study; it is not itself a biological discovery.

### Defense/Supp. Reviewer Defense And Execution Control
- Asset: PF07_reviewer_to_figure_map.png;PF14_analysis_control_board.png;F09_publishability_decision.png
- Role: 예상 공격, kill criteria, safe claim sentence를 한 번에 방어한다.
- Data used: reviewer risk heatmap; statistical analysis plan; stage-gate kill criteria; safe sentence bank.
- Where used: Supplement, internal review, Samsung Q&A rehearsal.
- Claim: The proposal is ambitious but bounded and stage-gated.
- Boundary: Defense artifacts are not new biological results.

## Representative Tables
- Data Definition Master (data_dictionary_master.tsv): 각 데이터셋이 무엇이고 무엇에 쓰였는지 정의 / Where used: Methods/Data section; reviewer data provenance / Boundary: forbidden_claim 컬럼을 반드시 같이 본다.
- Caption And Results Scaffold (figure_caption_and_results_scaffold.tsv): 각 figure의 caption, results paragraph, methods note, reviewer sentence / Where used: Figure legend drafting; manuscript Results scaffold / Boundary: 문장 scaffold이지 final claim expansion이 아니다.
- Representative Figure Overview (representative_figure_overview.tsv): 대표 figure spine과 data/claim/boundary / Where used: 교수님 briefing; manuscript figure plan / Boundary: 이 파일 자체는 guide다.
- Data-To-Figure-Claim Map (data_to_figure_claim_map.tsv): 데이터별 어디에 썼고 어떤 claim을 허용하는지 / Where used: Methods, rebuttal, grant review defense / Boundary: future data rows는 current result가 아니다.
- Claim Boundary Matrix (claim_boundary_matrix.tsv): Allowed/forbidden statement를 명시 / Where used: Reviewer defense and final editing / Boundary: forbidden 문장은 본문에 쓰지 않는다.
- Hospital Data Request CRF (hospital_data_request_crf.tsv): Yu/hospital에서 받아야 할 필드 / Where used: 7-day execution, IRB/CRF / Boundary: 요청 필드이지 이미 확보한 데이터가 아니다.
- Statistical Analysis Plan (statistical_analysis_plan.tsv): 현재 public-data analysis와 future prospective test 분리 / Where used: Methods and Samsung review / Boundary: future SAP는 현재 결과가 아니다.
- Stage-Gate Kill Criteria (stage_gate_kill_criteria.tsv): assay/NGS/follow-up 실패 시 downgrade rule / Where used: Budget/risk management / Boundary: kill criteria를 숨기지 않는다.
- Safe Sentence Bank (safe_sentence_bank.tsv): 본문에 바로 쓸 수 있는 안전 문장 / Where used: Writing and reviewer answer / Boundary: safe sentence 외 확장 시 claim boundary 재확인.

## Data Used / Where Used
- TCGA-THCA public pilot: used in Main Fig. 2; Evidence Figures; Slide 5; allowed: PTC tissue-state heterogeneity beyond driver mutation; BRAF-mutant state split; genotype-alone insufficiency rationale.; not allowed: CTC shedding, postoperative CTC transition, recurrence prediction, clinical diagnostic claim.
- GSE250521 spatial transcriptomics: used in Main Fig. 3; Evidence Figures; spatial coherence panels; allowed: Spatial tissue-state organization; marker/ROI context rationale.; not allowed: CTC origin proof or direct tumor-cell shedding source.
- PTC scRNA marker context: used in Main Fig. 4 (F06); Supplementary Planning Matrix (PF10 - hand-coded planning score, not data-derived); allowed: Cell/tissue context for epithelial, mesenchymal, thyroid-lineage, stress/survival modules.; not allowed: Definition of circulating CTC state without blood data.
- Bulk proteomics proxy: used in Main Fig. 4; F07; proteomics supplementary figures; allowed: Protein-direction support for dedifferentiation/marker-prior logic.; not allowed: Direct CTC protein phenotype or clinical assay validation.
- Yu 2024 published CTC study: used in Main Fig. 1; Breakthrough Case; Slide 3; proposal text; allowed: Feasibility anchor: serial PTC CTC-EMT measurement around thyroidectomy has been reported.; not allowed: Local raw reanalysis; genetic-map evidence; our cohort result.
- Future Yu/hospital serial CTC dataset: used in Future Main Fig. serial transition; Samsung Aim 1; Future Data Unlock Map; allowed: Required to test postoperative CTC-EMT state transition.; not allowed: Do not write as existing result until data are obtained and QC passes.
- Matched tumor-normal NGS: used in Main Fig. 5; Samsung Aim 2; NGS genetic map; allowed: Required to anchor CTC/cfDNA signals to resected tumor clone.; not allowed: Do not assume all low-input CTC NGS succeeds.
- Postoperative follow-up: used in Future Aim 3; data unlock map; reviewer defense; allowed: Secondary association with residual-risk features after adequate follow-up.; not allowed: Immediate recurrence prediction claim before validation.