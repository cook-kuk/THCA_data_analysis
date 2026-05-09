# Figure and Table Guide for Samsung PPT

## Recommended Figure Order

1. Public Pilot GO decision
2. Data sources and analysis flow
3. TCGA vulnerability landscape
4. Mutation-only insufficiency
5. Spatial niche coherence
6. Representative spatial niche maps
7. Integrated public validation matrix
8. External validation expansion
9. Drug/perturbation class hypotheses
10. Experimental validation design

---

## Core Figures

### 1. Data Flow

File:

- `figures/figure1_public_pilot_data_flow_dark.png`

Use:

- Show input data layers and output concepts.

Caption:

> Public TCGA, spatial transcriptomics, external thyroid cohorts, scRNA/GeoMx references, and drug resources were integrated into a K-Thyro therapeutic vulnerability scoring framework.

Claim boundary:

> Multi-layer evidence supports feasibility and validation design, not clinical deployment.

---

### 2. TCGA Vulnerability Heatmap

File:

- `figures/tcga_vulnerability_axes_heatmap_dark.png`

Use:

- Show patient-level heterogeneity across RAI, HLA/APM, immune, stromal/barrier, drug-delivery proxy, and aggressive axes.

Caption:

> TCGA-THCA 505 primary tumors separate into multiple therapeutic vulnerability states based on expression-derived K-Thyro axes.

Allowed:

- Patient-level axis separability.
- Hypothesis generation.

Forbidden:

- Spatial heterogeneity proof.
- Treatment response prediction.

---

### 3. TCGA Label Distribution

File:

- `figures/tcga_label_distribution_dark.png`

Use:

- Show subtype counts.

Key numbers:

- RAI-readable differentiated: 127
- Mixed/Other: 125
- Drug-delivery barrier-high: 105
- HLA-visible inflamed: 44
- CD8-excluded myeloid/CAF-high: 42
- HLA-low immune-invisible: 32
- RAI-low dedifferentiated: 30

Caption:

> K-Thyro scoring nominates distinct patient-level vulnerability labels rather than a single average-risk continuum.

---

### 4. Mutation Not Enough: Driver Label Fraction

File:

- `figures/tcga_mutation_not_enough_driver_label_fraction.png`

Use:

- Strong reviewer-facing evidence that mutation alone is insufficient.

Key numbers:

- BRAF group: 274 patients, 7 labels.
- BRAF largest label: 33.2%.
- BRAF ambiguity: 66.8%.

Caption:

> BRAF/RAS driver grouping is associated with vulnerability state but does not determine it; BRAF tumors split across all seven K-Thyro vulnerability labels.

Forbidden:

- Do not say mutation is irrelevant.

Correct wording:

> Mutation is necessary but insufficient for therapeutic vulnerability stratification.

---

### 5. Mutation Not Enough: Axis Variance

File:

- `figures/tcga_mutation_not_enough_axis_variance.png`

Use:

- Show driver group explains only part of axis variance.

Caption:

> Driver group explains a limited fraction of K-Thyro axis variance, especially for immune visibility and CD8 exclusion axes.

---

### 6. Spatial Coherence Barplot

File:

- `figures/spatial_coherence_barplot_dark.png`

Use:

- Main spatial evidence.

Key numbers:

- 16 slides.
- 57,144 spots.
- 16/16 slides z-score > 2.
- z-score min/median/max: 6.09 / 22.00 / 43.16.

Caption:

> Therapeutic vulnerability niches show non-random spatial coherence across all 16 GSE250521 slides.

Forbidden:

- Do not treat spots as independent patient replicates.
- Do not claim clinical response prediction.

---

### 7. Representative Spatial Maps

File:

- `figures/spatial_representative_maps_dark.png`

Use:

- Visual proof-of-concept.

Caption:

> Representative spatial maps show RAI differentiation, HLA/APM visibility, CAF/myeloid barrier, and composite niche labels occupying distinct tissue territories.

Claim boundary:

> Spatial transcriptomics maps mRNA expression territories; protein-level and functional validation are required.

---

### 8. Integrated Evidence Matrix

File:

- `figures/integrated_public_pilot_evidence_matrix_dark.png`

Use:

- Tie all axes together.

Caption:

> Integrated public-data evidence supports multiple separable therapeutic vulnerability axes, each linked to a specific experimental validation requirement.

---

### 9. Public Validation v2 Support Counts

File:

- `figures/public_validation_expansion_v2_support_counts.png`

Use:

- Show expanded public validation beyond TCGA and GSE250521.

Key numbers:

- 301 tests.
- strong 144.
- moderate 30.
- opposite 2.

Caption:

> Additional public cohorts and resource layers provide broad but heterogeneous support for K-Thyro axes; weak layers define validation priorities.

---

### 10. Exact External K-Thyro Heatmap / Boxplots

Files:

- `figures/exact_kthyro_external_bulk_axis_heatmap.png`
- `figures/exact_kthyro_external_bulk_boxplots.png`

Use:

- External reproducibility slide or appendix.

Key numbers:

- 206 external samples.
- 168 contrast tests.
- strong 65, moderate 12.
- prespecified direction 42 match / 2 opposite.

Caption:

> Independent external thyroid cohorts reproduce major K-Thyro axis directions, especially reduced RAI differentiation and increased aggressive/barrier-associated scores in advanced histologies.

---

### 11. Axis Separability Heatmaps

File:

- `figures/axis_separability_correlation_heatmaps.png`

Use:

- Answer “are these just the same score?”.

Caption:

> K-Thyro axes are correlated where biology overlaps, but they do not collapse into a single score, particularly within spatial slides.

---

### 12. Drug/Perturbation Candidate Pilot

File:

- `figures/figure6_drug_perturbation_candidate_pilot.png`

Use:

- Put late in the deck, after spatial/validation evidence.

Caption:

> Public drug resources nominate perturbation classes for validation, not final therapies.

---

### 13. Experimental Validation Design

File:

- `figures/figure7_samsung_experimental_validation_design.png`

Use:

- Close the scientific loop.

Caption:

> K-Thyro converts public-data therapeutic vulnerability niches into FFPE/mIHC, GeoMx/ROI, fresh tissue perturbation, iodide uptake, HLA/APM rescue, and drug-delivery imaging validation.

---

## Core Tables

### TCGA patient vulnerability scores

File:

- `tables/tcga_thca_patient_vulnerability_scores.tsv`

Use:

- Source for all TCGA patient-level scores and labels.

Important columns:

- `primary_vulnerability_label`
- `rai_differentiation_score`
- `hla_i_apm_score`
- `immune_visibility_score`
- `cd8_exclusion_proxy`
- `myeloid_caf_barrier_score`
- `drug_delivery_failure_proxy`
- `aggressive_dedifferentiation_score`
- `driver_anchor`

---

### Spatial slide niche summary

File:

- `tables/spatial_slide_niche_summary.tsv`

Use:

- Slide-level spatial fractions and coherence stats.

Important columns:

- `sample_id`
- `condition`
- `n_spots`
- `same_niche_z`
- `same_niche_empirical_p`
- niche fractions

---

### Public validation matrix v2

File:

- `tables/external_public_validation_signal_matrix_v2.tsv`

Use:

- Master evidence matrix.

Important columns:

- `dataset`
- `validation_layer`
- `axis`
- `endpoint`
- `effect_size_cohens_d_or_reported`
- `support_strength`
- `caveat`

---

### Exact K-Thyro external contrasts

File:

- `tables/exact_kthyro_external_bulk_contrasts.tsv`

Use:

- External reproducibility and directionality.

Important columns:

- `dataset`
- `contrast`
- `axis_label`
- `effect_size_cohens_d`
- `fdr_q_value`
- `expected_direction`
- `direction_match`
- `support_strength`

---

### TCGA driver diversity

File:

- `tables/tcga_driver_within_group_vulnerability_diversity.tsv`

Use:

- Mutation-only insufficiency argument.

Important columns:

- `driver_group`
- `n_patients`
- `n_nonzero_vulnerability_labels`
- `largest_label_fraction`
- `mutation_only_ambiguity_fraction`

---

### TCGA driver axis variance

File:

- `tables/tcga_driver_axis_variance_explained.tsv`

Use:

- Quantify how little driver explains some axes.

Important columns:

- `axis_label`
- `eta_squared_driver_with_no_call_group`
- `unexplained_fraction_1_minus_eta_with_no_call_group`
