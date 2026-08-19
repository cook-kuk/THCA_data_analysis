# CROSS-Neo-TCR Manuscript Track Scaffold

Working title: **CROSS-Neo-TCR: TCR-aware diagnostic prioritization of cancer neoantigens using external recognition experts**

Status: scaffold / companion-track candidate. This is not final manuscript prose.

## Core Thesis

Presentation is not recognition. CROSS-Neo-TCR keeps the main CROSS-Neo pMHC neoantigen ranker intact and adds an optional recognition-aware branch for rows with TCR evidence. The branch uses TCR-pMHC registry linkage, external sequence recognition experts, structure-readiness auditing, and wetlab prioritization diagnostics.

## 1. Motivation: Presentation Is Not Recognition

Evidence hooks:

- Main pMHC prediction is necessary but not sufficient for TCR recognition.
- TCR data are sparse in neoantigen benchmarks, so TCR-aware modeling cannot replace the main ranker.
- Figure hook: `fig_tcr1_recognition_gap`.

Safe claim:

- CROSS-Neo-TCR adds a diagnostic recognition layer for TCR-available cases.

Unsafe claim:

- Do not claim TCR recognition from peptide-HLA alone.

## 2. Dataset Registry And Missingness Reality

Key facts:

- Local TCR registry rows: 537,618.
- Paired alpha/beta rows: 64,037.
- Any TCR sequence rows: 123,035.
- Peptide-HLA rows: 406,977.
- Peptide-HLA-TCR labeled rows: 98,718 with any TCR sequence; 56,329 with paired alpha/beta.

Primary outputs:

- `tcr_registry.parquet`
- `tcr_registry_source_counts.tsv`
- `tcr_registry_missingness.tsv`
- `tcr_registry_report.md`

Figure hook:

- `fig_tcr1_recognition_gap`

## 3. Leakage And Source Shift Audit

Rules:

- No identical peptide-HLA leakage.
- No identical TCR leakage.
- No peptide-only label transfer from pathogen epitopes to cancer neoantigens.
- Source-heldout results must be interpreted separately from pooled results.

Source-heldout ensemble results:

- VDJdb heldout: AUPRC 0.6523, AUROC 0.6040.
- McPAS-TCR heldout: AUPRC 0.8025, AUROC 0.8232.
- VDJdb 10x heldout: AUPRC 0.8408, AUROC 0.8340.

Interpretation:

- Signal persists in some heldout sources but is source-dependent.
- This supports diagnostic prioritization, not universal TCR-aware prediction.

## 4. TCR-Available Subset Construction

Linkage facts:

- CROSS-Neo rows: 2,715.
- Any TCR-resource overlap: 2,016.
- Exact paired TCR-pMHC rows: 101.
- Cancer-context evidence rows: 1,876.
- No TCR match: 699.

Primary outputs:

- `tcr_neo_linked_registry.parquet`
- `tcr_neo_linkage_summary.tsv`
- `tcr_neo_linkage_report.md`

## 5. TCR-pMHC Structure Prediction Pipeline

Current status:

- Structure job manifest: 638 jobs.
- P0 jobs: 507.
- Main blocker: missing full TCR and MHC sequences for most rows.
- Structure outputs are QC/parser validation and missingness diagnostics, not primary recognition evidence.

Supported or audited tools:

- TCRdock
- TCRmodel2
- AlphaFold-Multimer / ColabFold
- AlphaFold3
- Boltz / Boltz-2
- Chai-1
- tFold-TCR

Figure hooks:

- `fig_tcr3_structure_pipeline`
- `fig_tcr4_mutant_wt_interface_delta`

## 6. Mutant-WT Counterfactual Interface Modeling

Current status:

- Not yet supported as a result claim.
- Requires mutant and WT TCR-pMHC structures or high-confidence modeled complexes.

Planned features:

- peptide-TCR contact delta
- mutation residue contact flag
- mutant-WT interface confidence delta
- WT cross-reactivity risk

Safe claim:

- This is a planned diagnostic branch and wetlab prioritization feature.

Unsafe claim:

- Do not claim mutant-WT interface delta improves ranking until paired structure outputs exist.

## 7. Source-Robust And TCR-Heldout Evaluation

Fast external expert benchmark:

- pMTnet alone: AUPRC 0.7479, AUROC 0.7255.
- TEPCAM alone: AUPRC 0.6878, AUROC 0.6970.
- Mean pMTnet/TEPCAM ensemble: AUPRC 0.7841, AUROC 0.7652.
- Peptide-HLA GroupKFold logistic ensemble: AUPRC 0.7803, AUROC 0.7683.

Interpretation:

- pMTnet and TEPCAM are complementary enough to justify a two-expert diagnostic score.
- Synthetic decoy benchmark supports model selection only.

Figure hook:

- `fig_tcr5_tcr_available_subset_performance`

## 8. Case Studies

Wetlab candidate external expert rescoring:

- HMTEVVRHC / HLA-A*02:01: 4 exact TCR rows scored; external mean 0.7121; max 0.8029; top external-adjusted candidate.
- GADGVGKSAL / HLA-C*08:02: 5 exact TCR rows scored; external mean 0.5703; max 0.6790; second external-supported candidate.

Primary outputs:

- `wetlab_candidates_external_tcr_expert_ranked.tsv`
- `wetlab_external_tcr_expert_pair_scores.tsv`
- `05_wetlab_external_tcr_expert_scores.md`

Figure hook:

- `fig_tcr6_case_studies`

## 9. Claim Boundary

Allowed if supported:

- TCR-aware features improve or reprioritize cases on TCR-available diagnostic subsets.
- External TCR recognition experts provide useful wetlab prioritization evidence.
- Missing paired TCR data is a major limitation in current neoantigen benchmarks.

Forbidden:

- Universal TCR-aware neoantigen prediction.
- Clinical utility.
- AlphaFold or other structure prediction is always correct.
- External validation unless strict heldout datasets support it.
- TCR recognition inferred from peptide-HLA alone.

Figure hook:

- `fig_tcr7_claim_boundary`

## 10. Wetlab Validation Plan

Priority candidates:

1. HMTEVVRHC / HLA-A*02:01.
2. GADGVGKSAL / HLA-C*08:02.

Suggested assays:

- peptide-HLA multimer staining with candidate TCRs where available.
- TCR transduction or reporter assay for exact TCR-pMHC pairs.
- mutant versus WT peptide specificity assay when WT peptide exists.
- dose-response activation readout.
- negative decoy peptide controls.

Decision:

- Current evidence supports **PROMOTE_TO_DIAGNOSTIC_CASE_STUDY** and **PROMOTE_TO_WETLAB_PRIORITIZATION_TOOL**.
- Hold **PROMOTE_TO_TCR_AWARE_METHOD** until strict external paired-TCR benchmarks and, ideally, wetlab validation are complete.
