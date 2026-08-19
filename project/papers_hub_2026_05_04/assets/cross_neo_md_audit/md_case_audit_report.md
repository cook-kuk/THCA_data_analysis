# CROSS-Neo + MD Integration Audit

Generated from `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_md_audit_2026_05_10` and `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09`.

## Boundary

MD evidence is used as structural audit evidence only. It is not an immunogenicity label, clinical validation, or SOTA proof.

## Candidate-Level MD Evidence

| candidate              | md_representative_run_id   | MD_evidence_label   |   MD_evidence_score |   pMHC_stability_score |   TCR_recognition_score |   runtime_fraction |   live_time_ps | md_interpretation                                  |
|:-----------------------|:---------------------------|:--------------------|--------------------:|-----------------------:|------------------------:|-------------------:|---------------:|:---------------------------------------------------|
| GADGVGKSAL/HLA-C*08:02 | prod_10ns_6UON_2fs300K     | MD_MODERATE         |            0.556608 |               0.520263 |                  0.4404 |                  1 |          10000 | usable_structural_support_not_immunogenicity_proof |
| HMTEVVRHC/HLA-A*02:01  | prod_10ns_6VRN_1fs300K     | MD_STRONG           |            0.743344 |               0.666694 |                  1      |                  1 |          10000 | usable_structural_support_not_immunogenicity_proof |

## Joined Prediction Rows

- prediction rows with MD candidate key: 484
- MD-supported model-high rows: 182
- model-high / MD-low-or-insufficient rows: 0
- model-low / MD-high rows: 48

## Interpretation

- `GADGVGKSAL / HLA-C*08:02` has completed 10 ns explicit-solvent evidence and a 1 ns replicate screen, so it can be discussed as structural plausibility/wetlab prioritization evidence if the contact and QC plots remain stable.
- `HMTEVVRHC / HLA-A*02:01` is still running in the primary 10 ns job; state traces are useful for live QC, but trajectory-derived contact claims must wait for DCD sync and rerun of the analysis.
- No WT or scrambled-decoy trajectory is present yet, so mutant-specific TCR recognition and cross-reactivity claims remain forbidden.

## Outputs

- `integrated/cross_neo_md_joined.tsv`
- `integrated/md_supported_top_candidates.tsv`
- `integrated/model_high_md_low_cases.tsv`
- `integrated/model_low_md_high_cases.tsv`
