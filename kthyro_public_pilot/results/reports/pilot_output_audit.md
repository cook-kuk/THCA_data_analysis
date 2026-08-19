# Pilot Output Audit

## Expected Files

- tcga: EXISTS `results/tables/tcga_thca_patient_vulnerability_scores.tsv`
- spatial_spots: EXISTS `results/tables/spatial_spot_vulnerability_scores.tsv`
- spatial_slides: EXISTS `results/tables/spatial_slide_niche_summary.tsv`
- spatial_coherence: EXISTS `results/tables/spatial_coherence_statistics.tsv`
- drug_candidates: EXISTS `results/tables/drug_pilot_candidate_rankings.tsv`
- evidence_matrix: EXISTS `results/tables/integrated_public_pilot_evidence_matrix.tsv`
- proposal_figures: EXISTS `results/figures/proposal`
- public_report: EXISTS `results/reports/public_pilot_report.md`
- go_decision: EXISTS `results/reports/GO_NO_GO_DECISION.md`
- claim_boundaries: EXISTS `results/reports/claim_boundaries.md`
- validation_plan: EXISTS `results/reports/experimental_validation_plan.md`

## Row/Column Counts

- TCGA patient vulnerability: 505 rows x 74 columns
- Spatial spot vulnerability: 57144 rows x 33 columns
- Spatial slide summary: 16 rows x 22 columns
- Spatial coherence: 208 rows x 8 columns
- Drug candidates: 140 rows x 17 columns
- Integrated evidence matrix: 8 rows x 12 columns

## Identity And Duplication Checks

- TCGA canonical sample column detected: `sample_id_x`
- TCGA row count: 505
- TCGA unique patients: 505
- TCGA duplicated sample records: 0
- TCGA missing labels: 0
- Suspicious but non-blocking: TCGA table contains `sample_id_x/sample_id_y` from clinical merge. Use `sample_id_x` or `patient_id` for proposal summaries.
- Spatial spot rows: 57144
- Spatial slide count from spots: 16
- Spatial duplicated sample_id+spot_id records: 0
- Spatial spots are nested within slide via `sample_id`; group claims must be slide-level.
- Spatial slide summary rows: 16
- Spatial conditions: {'normal': 4, 'PTC': 4, 'locally_advanced_PTC': 4, 'ATC': 4}
- KNN same-niche coherence rows: 16
- Slides with same-niche z > 2: 16/16
- Minimum empirical p for same-niche test: 0.0019960079840319
- Drug candidate duplicated `drug_label` rows: 63
- Drug table is usable only after compound/class de-duplication; target-level duplicates are not proposal-ready.

## Missing Fields And Usability

- TCGA: missing fields = none; proposal usable = True
- Spatial spots: missing fields = none; proposal usable = True
- Spatial slides: missing fields = none; proposal usable = True
- Drug: missing fields = none; proposal usable = True
