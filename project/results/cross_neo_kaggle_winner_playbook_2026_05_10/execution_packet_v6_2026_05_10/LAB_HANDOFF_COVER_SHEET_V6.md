# CROSS-Neo v6 lab handoff cover sheet

Generated: 2026-05-10T13:51:53

This packet contains the files needed to run the v5 assay plan as a traceable execution package.

Core files:

- `peptide_hla_reagent_order_manifest_v6.tsv`: peptide/HLA/control/provenance order and audit lines.
- `plate_96_execution_map_v6.tsv`: 96-well map with blinded well IDs and result-entry fields.
- `candidate_result_entry_v6.tsv`: candidate-level PASS/FAIL sheet consumed by the interpreter.
- `well_result_entry_v6.tsv`: well-level raw/normalized signal entry sheet.
- `decision_worksheet_v6.tsv`: endpoint-level manual decision worksheet.
- `STATISTICAL_ANALYSIS_PLAN_V6.md`: frozen analysis rules.

Do not edit endpoint thresholds after assay results are visible. If a candidate is excluded, keep the row and mark `EXCLUDE` with an explicit reason.
