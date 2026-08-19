# CROSS-Neo-TCR Reviewer Packet

Compact packet for reviewing the TCR-aware diagnostic extension.

Decision: PROMOTE_TO_DIAGNOSTIC_CASE_STUDY and PROMOTE_TO_WETLAB_PRIORITIZATION_TOOL. Hold main-method promotion until strict external paired-TCR benchmarks and structure/wetlab validation are complete.

Key headline:

- pMTnet + TEPCAM mean ensemble: AUPRC 0.7841, AUROC 0.7652 on six fast challenge panels.
- Top external-supported wetlab candidate: HMTEVVRHC / HLA-A*02:01.
- Second external-supported candidate: GADGVGKSAL / HLA-C*08:02.
- MD escalation queue: P0 TCR-pMHC MD for HMTEVVRHC/HLA-A*02:01 and GADGVGKSAL/HLA-C*08:02; existing PDB pilots are available for both.
- OpenMM repair/minimization sanity check succeeded for representative P0 complexes 6VRN and 6UON.
- Local OpenMM smoke dynamics succeeded for both P0 representatives; production still requires CUDA explicit-solvent pilot MD.

Contents:

- `reports/`: decision report, literature/tool audit, registry/linkage reports, model audit, benchmark reports, manuscript-track scaffold.
- `tables/`: source counts, linkage summary, external model metrics, source-heldout metrics, wetlab ranked candidates, MD escalation queue.
- `figures/`: seven TCR extension figures as PNG and PDF.

Claim boundary:

- This packet supports an optional TCR-aware diagnostic/wetlab-prioritization branch.
- It does not support universal TCR-aware neoantigen prediction or clinical utility.
- Synthetic decoy benchmarks are model-selection diagnostics, not final external SOTA evidence.

File manifest:

| path | bytes |
|---|---:|
| `figures/fig_tcr1_recognition_gap.pdf` | 23186 |
| `figures/fig_tcr1_recognition_gap.png` | 240126 |
| `figures/fig_tcr2_architecture.pdf` | 22572 |
| `figures/fig_tcr2_architecture.png` | 342446 |
| `figures/fig_tcr3_structure_pipeline.pdf` | 26914 |
| `figures/fig_tcr3_structure_pipeline.png` | 221806 |
| `figures/fig_tcr4_mutant_wt_interface_delta.pdf` | 21758 |
| `figures/fig_tcr4_mutant_wt_interface_delta.png` | 151068 |
| `figures/fig_tcr5_tcr_available_subset_performance.pdf` | 27610 |
| `figures/fig_tcr5_tcr_available_subset_performance.png` | 201560 |
| `figures/fig_tcr6_case_studies.pdf` | 26530 |
| `figures/fig_tcr6_case_studies.png` | 239289 |
| `figures/fig_tcr7_claim_boundary.pdf` | 25417 |
| `figures/fig_tcr7_claim_boundary.png` | 256267 |
| `packages/openmm_p0_10ns_pilot_package_2026_05_09.tar.gz` | 5785734 |
| `reports/00_tcr_lit_tool_audit.md` | 18585 |
| `reports/01_tcr_model_discovery.md` | 8419 |
| `reports/02_tcr_model_runtime_audit.md` | 5201 |
| `reports/04_external_tcr_expert_benchmark.md` | 1993 |
| `reports/05_wetlab_external_tcr_expert_scores.md` | 1714 |
| `reports/CROSS_Neo_TCR_extension_decision_report.md` | 7301 |
| `reports/TCR_EXTENSION_REVIEWER_QA.md` | 4848 |
| `reports/manuscript_tcr_aware_track.md` | 5799 |
| `reports/md_escalation_report.md` | 2366 |
| `reports/md_job_stub_report.md` | 568 |
| `reports/openmm_p0_minimization_report.md` | 1037 |
| `reports/openmm_p0_pilot_package_README.md` | 567 |
| `reports/p0_openmm_smoke_report.md` | 606 |
| `reports/p0_structure_pilot_report.md` | 1857 |
| `reports/tcr_extension_figure_report.md` | 2420 |
| `reports/tcr_feature_report.md` | 1333 |
| `reports/tcr_neo_linkage_report.md` | 2474 |
| `reports/tcr_registry_report.md` | 3454 |
| `reports/tcr_structure_pipeline_report.md` | 1230 |
| `scripts/analyze_openmm_pilot.py` | 4387 |
| `scripts/run_openmm_pilot.py` | 6893 |
| `scripts/runpod_setup_openmm.sh` | 286 |
| `scripts/submit_p0_10ns_pilots.sh` | 577 |
| `tables/external_tcr_expert_benchmark_metrics.tsv` | 4414 |
| `tables/external_tcr_expert_source_heldout_metrics.tsv` | 477 |
| `tables/external_tcr_model_pilot_comparison.tsv` | 764 |
| `tables/md_engine_availability.tsv` | 158 |
| `tables/md_escalation_queue.tsv` | 11012 |
| `tables/md_escalation_queue_top20.tsv` | 8792 |
| `tables/md_job_manifest.tsv` | 6313 |
| `tables/md_simulation_tiers.tsv` | 746 |
| `tables/openmm_p0_minimization_qc.tsv` | 1343 |
| `tables/openmm_p0_pilot_manifest.tsv` | 766 |
| `tables/p0_exact_pdb_registry_rows.tsv` | 1496 |
| `tables/p0_md_pilot_ready_complexes.tsv` | 2966 |
| `tables/p0_openmm_smoke_qc.tsv` | 1399 |
| `tables/p0_pdb_chain_qc.tsv` | 2859 |
| `tables/submit_md_jobs.template.sh` | 11645 |
| `tables/tcr_model_discovery_matrix.tsv` | 7138 |
| `tables/tcr_neo_linkage_summary.tsv` | 423 |
| `tables/tcr_registry_label_balance.tsv` | 227 |
| `tables/tcr_registry_missingness.tsv` | 1116 |
| `tables/tcr_registry_source_counts.tsv` | 590 |
| `tables/wetlab_candidates_external_tcr_expert_ranked.tsv` | 17658 |
| `tables/wetlab_external_tcr_expert_pair_scores.tsv` | 2697 |
