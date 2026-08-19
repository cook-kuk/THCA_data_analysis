# Figure Caption And Results Scaffold

## Fig. 1. Manuscript logic linking postoperative residual-risk uncertainty, published PTC CTC feasibility, public multi-omics prior construction, and the proposed prospective perturbation test.

**Image:** `PF01_manuscript_flow.png`

**Long caption:** The paper should open by separating what is already known from what remains untested. Yu 2024 supports the feasibility of serial CTC-EMT measurement around thyroidectomy, while public TCGA, spatial, single-cell, and proteomics resources support a tissue-state and marker-prior framework. The actual thyroidectomy-induced CTC-EMT transition remains a prospective test requiring hospital serial blood and matched genetic data.

**Results scaffold:** We first defined the study as a staged translational framework rather than a completed CTC-transition discovery analysis. The logic proceeds from the clinical gap after thyroidectomy to a published CTC feasibility anchor, then to public-data-derived tissue-state priors, and finally to a prospective serial blood and matched NGS experiment.

**Methods/provenance note:** This figure is a synthesis schematic. It does not present a new statistical test; it maps evidence layers to manuscript sections and claim boundaries.

**Reviewer defense sentence:** This is a Science-track biological perturbation proposal because thyroidectomy is the perturbation and CTC-EMT state transition is the primary biological readout.

## Fig. 2. Data provenance map defining how each public, published, planned, and vendor layer is used and bounded.

**Image:** `PF02_data_provenance_map.png`

**Long caption:** TCGA-THCA, spatial transcriptomics, scRNA marker context, bulk proteomics proxy, Yu 2024, the future Samsung/Yu cohort, and Inocras/vendor workflows have distinct roles. Public data provide tissue-state, spatial, marker, and protein-direction priors; Yu 2024 provides published CTC feasibility; the future cohort is required to test serial CTC-EMT and genetic mapping.

**Results scaffold:** We organized all evidence sources into explicit functional categories. No public dataset is treated as a surrogate for serial CTC measurement. Instead, each public layer contributes a defined prior that informs the design of a future liquid-biopsy perturbation experiment.

**Methods/provenance note:** Dataset roles were assigned from locally generated pilot outputs and proposal reports. The map is intentionally conservative and includes each layer's non-supported claim.

**Reviewer defense sentence:** The public pilot is not being used as CTC proof; it is being used to define a rational marker and genotype interpretation framework.

## Fig. 3. BRAF-mutant PTC tumors distribute across multiple public tissue-state labels rather than forming a single biological state.

**Image:** `F03_tcga_braf_state_split.png`

**Long caption:** In the TCGA-THCA BRAF-mutant subset, tumors were distributed across seven tissue-state/vulnerability labels, and the largest label accounted for only approximately one third of BRAF-mutant cases. This supports the central premise that driver mutation alone is insufficient to define the biological state relevant to postoperative blood-signal interpretation.

**Results scaffold:** Among BRAF-mutant PTC tumors, mutation status did not collapse tumors into a single state. The largest state represented only 33.2% of the BRAF-mutant subset, leaving substantial state heterogeneity within the same driver-defined group.

**Methods/provenance note:** Input table: tcga_braf_only_label_distribution.tsv generated from the public pilot extra analyses.

**Reviewer defense sentence:** This figure does not claim CTC biology; it justifies why tissue genotype must be interpreted together with state markers.

## Fig. 4. Driver-only models explain only part of public tissue-state axis variation.

**Image:** `F04_tcga_driver_variance_boundary.png`

**Long caption:** Driver grouping provides useful information but does not fully account for RAI, immune, HLA, CD8, or related state axes in the public TCGA analysis. This supports a phenotype-plus-genotype strategy rather than a mutation-only interpretation of postoperative blood signals.

**Results scaffold:** Driver-only models captured only a fraction of tissue-state variation. These residual state axes motivate the combined CTC-EMT phenotype and matched NGS design, where mutation identifies clone origin and phenotype captures transition state.

**Methods/provenance note:** Input table: tcga_axis_variance_models.tsv. Models are descriptive and not causal.

**Reviewer defense sentence:** The driver analysis is a boundary argument: mutation helps, but it is not enough for the state-transition question.

## Fig. 5. Spatial thyroid tissue states show non-random organization across public spatial transcriptomic slides.

**Image:** `F05_spatial_tissue_state_organization.png`

**Long caption:** GSE250521 spatial transcriptomics supports the organization of tissue-state labels in thyroid tissue. Same-niche coherence exceeded z>2 in all 16 slides analyzed, across 57,144 spots. This supports ROI and tissue-context logic for marker selection but does not establish the physical origin of circulating tumor cells.

**Results scaffold:** Public spatial data showed that tissue-state labels are spatially structured rather than random spot-level noise. This provides tissue-context support for the marker-prior framework that will later be tested in serial blood.

**Methods/provenance note:** Input table: spatial_coherence_statistics.tsv; source figure: spatial_coherence_summary.png.

**Reviewer defense sentence:** Spatial coherence is used as tissue-context evidence only, not as evidence of CTC shedding.

## Fig. 6. Single-cell marker context supports the biological plausibility of epithelial, mesenchymal, thyroid-lineage, and stress/survival modules.

**Image:** `F06_scrna_marker_context.png`

**Long caption:** Single-cell module attribution provides cell-context rationale for selecting CTC-EMT and thyroid-lineage marker modules. The goal is not to validate a blood assay from tissue scRNA data, but to avoid an arbitrary marker panel and to ensure that each module has an interpretable biological context.

**Results scaffold:** We used single-cell marker attribution to organize candidate CTC modules into epithelial, hybrid/mesenchymal, thyroid-lineage, and survival/stress compartments. This creates a structured marker prior for the prospective assay.

**Methods/provenance note:** Input table: scrna_module_celltype_attribution.tsv; source figure: scrna_module_celltype_attribution_heatmap.png.

**Reviewer defense sentence:** The scRNA layer supports marker selection, not direct circulating-cell validation.

## Fig. 7. Bulk proteomics proxy provides orthogonal support for dedifferentiation-related directionality.

**Image:** `F07_proteomics_dediff_direction.png`

**Long caption:** Bulk proteomics contrasts support loss of thyroid/RAI-related programs and gain of stress, myeloid, or TGF-beta-associated programs along dedifferentiation-related comparisons. This provides protein-level directionality for the tissue-state prior, while remaining distinct from spatial proteomics or CTC proteomics.

**Results scaffold:** Orthogonal proteomic summaries were consistent with the direction of dedifferentiation-prior modules. This supports including thyroid-lineage loss and survival/stress gain in the CTC-EMT panel design.

**Methods/provenance note:** Input tables: bulk_proteomics_kthyro_module_contrasts.tsv and bulk_proteomics_kthyro_dediff_trends.tsv.

**Reviewer defense sentence:** The proteomics layer strengthens biological plausibility but does not substitute for serial blood evidence.

## Fig. 8. Phenotypic CTC state modules are connected to tumor-informed genetic anchoring.

**Image:** `PF03_marker_to_ngs_map.png`

**Long caption:** The proposed interpretation framework separates phenotypic CTC state from clone identity. Epithelial, hybrid E/M, mesenchymal, thyroid-lineage, and survival/stress markers define candidate cell-state fractions; matched tissue-normal NGS and cfDNA/CTC-enriched NGS are required to determine whether persistent postoperative blood signals are genetically anchored to the resected tumor.

**Results scaffold:** We therefore defined the CTC-EMT panel as a two-layer system: phenotype first, genetic anchoring second. This prevents overinterpretation of nonspecific postoperative blood cells as tumor-derived residual-risk biology.

**Methods/provenance note:** Input sources: ctc_emt_marker_prior_table.tsv, specific aims, Samsung NGS map schematic.

**Reviewer defense sentence:** Matched NGS is the guardrail against overcalling CTC phenotype as tumor-derived biology.

## Fig. 9. CTC-EMT-NGS prior panel for post-thyroidectomy monitoring in PTC.

**Image:** `F08_ctc_emt_ngs_prior_panel.png`

**Long caption:** The candidate panel includes epithelial, hybrid E/M, mesenchymal, thyroid identity, survival/stress, immune/stress, and genetic-anchor modules. It is intended as a research-grade prior for prospective testing, not as a locked clinical diagnostic assay.

**Results scaffold:** The final prior panel integrates tissue-state evidence and thyroid-cancer genetic context into a candidate serial blood assay. The central readout is not CTC count alone, but the relationship among CTC state, thyroid-lineage signal, and tumor-informed genetic persistence.

**Methods/provenance note:** Input table: ctc_emt_marker_prior_table.tsv.

**Reviewer defense sentence:** This panel is a hypothesis-testing framework, not a validated commercial test.

## Fig. 10. The missing integrated serial CTC-EMT-NGS dataset defines the prospective study need.

**Image:** `PF04_dataset_gap_bridge.png`

**Long caption:** No local or public dataset currently combines serial pre/post-thyroidectomy blood, CTC EMT phenotype, matched tissue-normal NGS, cfDNA or CTC-enriched NGS, and postoperative follow-up. This gap is not a weakness of the public analysis; it is the direct justification for the Samsung Science proposal.

**Results scaffold:** The public-data analysis identifies the exact dataset that must be generated next. The missing elements are serial blood, CTC EMT subtyping, tissue-normal clone anchoring, postoperative cfDNA or CTC-enriched sequencing, and longitudinal clinical follow-up.

**Methods/provenance note:** Input table: serial_ctc_ngs_dataset_gap_map.tsv and proposal study design.

**Reviewer defense sentence:** The proposal is necessary precisely because the integrated dataset does not exist.

## Fig. 11. Allowed and forbidden claims for the current evidence package.

**Image:** `PF05_claim_boundary_matrix.png`

**Long caption:** The claim-boundary matrix separates publishable current statements from statements that require hospital serial CTC data or future validation. This is a core reviewer defense element because it shows that the proposal is ambitious without treating preliminary public data as clinical proof.

**Results scaffold:** We explicitly separated current evidence from future validation. Allowed claims include public tissue-state heterogeneity and published CTC feasibility; forbidden claims include new postoperative CTC transition, recurrence prediction, universal early-PTC NGS success, and clinical diagnostic readiness.

**Methods/provenance note:** Input table: claim_boundary_matrix.tsv and claim_boundaries_and_kill_criteria.md.

**Reviewer defense sentence:** The most credible version of this project is the version that states exactly what it cannot yet claim.

## Fig. 12. Reviewer attack map linking each expected criticism to figure/table support and a bounded answer.

**Image:** `PF07_reviewer_to_figure_map.png`

**Long caption:** Expected reviewer criticisms include category confusion, absence of public CTC data, mutation-only sufficiency, spatial-data relevance, CTC nonspecificity, and excessive breadth. Each criticism is paired with a specific figure/table and a constrained answer.

**Results scaffold:** We pre-specified reviewer attacks and matched them to evidence. This table should be used to keep oral and written responses narrow, especially when reviewers challenge the absence of local serial CTC data.

**Methods/provenance note:** Input sources: reviewer defense table, figure_table_manifest.tsv, claim_boundary_matrix.tsv.

**Reviewer defense sentence:** When attacked for lacking CTC data, answer directly: correct, and that is why public data are used only as prior evidence.

## Fig. 13. CV2 visual atlas summarizing the full figure-first manuscript package.

**Image:** `PF08_cv2_visual_atlas.png`

**Long caption:** The CV2 atlas is designed for fast human review: it places logic, provenance, marker-to-NGS interpretation, dataset gap, claim boundaries, reviewer defense, and key public-data figures into one large visual board. It is a navigation figure, not a substitute for the underlying tables.

**Results scaffold:** For internal and professor-facing communication, we generated a CV2 visual atlas to make the evidence structure inspectable at a glance. This supports rapid review while preserving links to source tables and individual figure explanations.

**Methods/provenance note:** Generated with cv2 from the deployed figure assets.

**Reviewer defense sentence:** Use this as an overview board, not as a standalone quantitative result.

## Fig. 14. Decision matrix distinguishing what can be published now from what requires Yu/hospital raw data.

**Image:** `F09_publishability_decision.png`

**Long caption:** A public-data framework or methods-rationale paper is feasible now. A postoperative CTC-EMT transition discovery paper or residual-risk prediction paper is not feasible without patient-level serial CTC, matched NGS, and follow-up data.

**Results scaffold:** The immediate manuscript should be framed as a public-data prior-map paper. The higher-impact biological discovery manuscript becomes feasible only after the prospective serial dataset is generated or obtained.

**Methods/provenance note:** Input table: publishability_decision_table.tsv.

**Reviewer defense sentence:** The current manuscript is intentionally scoped as a prior-map paper; the discovery paper is the next data-dependent stage.

## Fig. 15. Manuscript claim scorecard translating in silico outputs into writeable Results blocks.

**Image:** `PF09_manuscript_claim_scorecard.png`

**Long caption:** The scorecard separates data-backed public in silico claims from conceptual or future-validation claims. It identifies which results can be written now, which figures support each result, and why specific claims remain bounded.

**Results scaffold:** The current evidence package yields six writeable in silico result blocks: mutation/state heterogeneity, spatial tissue-state organization, marker-module rationale, phenotype-to-genotype anchoring, integrated dataset gap definition, and publishability scope.

**Methods/provenance note:** Input table: manuscript_claim_evidence_scorecard.tsv.

**Reviewer defense sentence:** This scorecard prevents moderate public-data evidence from being overstated as prospective CTC discovery.

## Supplementary planning matrix. Supplementary planning matrix of hand-coded support scores for CTC-EMT marker modules (manuscript-planning judgement, not data-derived).

**Image:** `PF10_marker_module_support_heatmap.png`

**Long caption:** Marker modules were hand-scored (0/1/2) across TCGA tissue-state context, spatial context, single-cell context, proteomics directionality, and the Yu 2024 published anchor. The scores reflect manuscript-planning judgement about how useful each evidence layer is for prior construction. They are not derived from raw data, are not measured statistics, and are not validation metrics. TCGA-THCA contains tissue (RNA/mutation/methylation/clinical) only and has no CTC, blood, or liquid biopsy data. The Yu 2024 column scores the published feasibility summary; no unpublished or raw Yu CTC data was used.

**Results scaffold:** As a supplementary planning matrix, we documented our judgement of how useful each public/published evidence layer was for constructing each marker-module prior. Thyroid-lineage modules carry the strongest public prior judgement, while EMT and survival/stress modules remain hypotheses requiring serial blood validation and matched genetic anchoring.

**Methods/provenance note:** Input table: marker_module_cross_layer_support_matrix.tsv. Scores are hand-coded manuscript-planning judgements (0/1/2), not measurements: 0=no current support, 1=indirect/context support, 2=directly useful for prior construction. TCGA column = TCGA-THCA tissue-state context (no CTC data). Yu 2024 column = scoring of the published Yu 2024 paper's reported feasibility (62 PTC, 87% CTC detection); no unpublished/raw Yu CTC data was accessed.

**Reviewer defense sentence:** This supplementary matrix is a manuscript-planning score, not a validated assay-performance matrix; readers should not interpret it as a cross-layer measurement.

## Fig. 17. Manuscript blueprint for the current no-Yu-data in silico prior/framework paper.

**Image:** `PF11_manuscript_blueprint_board.png`

**Long caption:** The blueprint board fixes the safe manuscript identity, recommended title, core thesis, allowed Results blocks, and boundaries. It separates a feasible translational framework paper from a future high-impact CTC-transition discovery paper that requires serial hospital data.

**Results scaffold:** We finalized the manuscript as a public-data prior/framework paper rather than a discovery paper. The writeable Results sequence is mutation/state heterogeneity, spatial tissue-state organization, marker-context rationale, marker support grading, phenotype-to-genotype anchoring, and integrated dataset-gap definition.

**Methods/provenance note:** Input tables: manuscript_blueprint.tsv and submission_strategy_matrix.tsv.

**Reviewer defense sentence:** This board is a scope-control device: it prevents the manuscript from drifting into unsupported CTC transition or recurrence-prediction claims.

## Fig. 18. Future data unlock map separating current public-data claims from future serial CTC-EMT discovery claims.

**Image:** `PF12_future_data_unlock_map.png`

**Long caption:** The data-unlock map specifies which hospital/Yu data elements are required to move beyond the current public-data prior manuscript. Serial CTC phenotype, CTC marker intensity, matched tissue-normal NGS, serial cfDNA, CTC-enriched NGS, pathology fields, and postoperative follow-up each unlock specific future figures and claims.

**Results scaffold:** We explicitly mapped the data elements required for the next-stage discovery manuscript. Serial CTC phenotype data would unlock transition plots, matched tissue-normal NGS would unlock clone anchoring, cfDNA would unlock molecular persistence kinetics, and postoperative follow-up would unlock risk-association analyses.

**Methods/provenance note:** Input table: future_data_unlock_map.tsv.

**Reviewer defense sentence:** Future data requirements are shown to prevent unsupported future claims from being imported into the current public-data manuscript.

## Fig. 19. Execution board for assembling the in silico manuscript package.

**Image:** `PF13_execution_board.png`

**Long caption:** The execution board combines the proposed main/supplementary figure layout, reviewer-risk severity, and a 72-hour execution sequence. It is designed to convert the current dossier into a manuscript draft and professor-facing data request without expanding unsupported claims.

**Results scaffold:** We converted the in silico dossier into an execution-ready package, defining main and supplementary figures, methods reproducibility checks, reviewer risks, safe sentences, hospital CRF requests, and a 72-hour action plan.

**Methods/provenance note:** Input tables: main_supp_figure_layout.tsv, reviewer_risk_heatmap.tsv, next_72h_manuscript_execution.tsv.

**Reviewer defense sentence:** Execution planning is separated from biological evidence to keep manuscript operations from becoming unsupported claims.

## Fig. 20. Analysis control board locking key numbers, statistical plan, and stage-gate criteria.

**Image:** `PF14_analysis_control_board.png`

**Long caption:** The analysis control board consolidates the key numerical ledger, statistical analysis plan, and stage-gate kill criteria. It is designed to prevent numeric drift in manuscript writing and to keep current public-data analyses separate from future prospective validation analyses.

**Results scaffold:** We created a final analysis-control layer that records key numbers, current versus future statistical questions, and stage-gate criteria. This converts the dossier from a figure collection into a controlled manuscript-preparation system.

**Methods/provenance note:** Input tables: key_number_ledger.tsv, statistical_analysis_plan.tsv, stage_gate_kill_criteria.tsv.

**Reviewer defense sentence:** The analysis control board makes clear which analyses were performed now and which require future serial hospital data.
