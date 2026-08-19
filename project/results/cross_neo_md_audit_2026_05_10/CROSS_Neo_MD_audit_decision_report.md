# CROSS-Neo MD Audit Decision Report

## 1. Executive Verdict

- `GADGVGKSAL / HLA-C*08:02`: **PROMOTE_TO_WETLAB_PRIORITIZATION_TOOL**, with MD evidence label `MD_MODERATE`. The completed 10 ns run supports structural plausibility, not immunogenicity proof.
- `HMTEVVRHC / HLA-A*02:01`: **HOLD_AS_RUNNING_P0_MD_AUDIT** until the 10 ns DCD is synced and contact/QC analysis is rerun. Live state is stable so far, but trajectory-derived claims are not ready.
- Overall: keep MD as a **structure-aware diagnostic audit layer**. Do not move it into the main CROSS-Neo claim without WT/decoy controls, replicates, and external validation.

## 2. Which Simulations Completed

| run_id                                 | candidate              | condition                      |   time_ps |   target_ns |   temperature_k |   speed_ns_per_day |
|:---------------------------------------|:-----------------------|:-------------------------------|----------:|------------:|----------------:|-------------------:|
| prod_10ns_6VRN_1fs300K                 | HMTEVVRHC/HLA-A*02:01  | explicit_cuda_10ns             |   10000   |        10   |        299.891  |             21     |
| 6VRM_HMTEVVRHC_HLA-A0201_chainP_0p5ns  | HMTEVVRHC/HLA-A*02:01  | explicit_cuda_replicate_screen |     500   |         0.5 |        300.355  |             22.4   |
| prod_10ns_6VRN_1fs300K                 | HMTEVVRHC/HLA-A*02:01  | explicit_cuda_10ns             |   10000   |        10   |        299.891  |             21     |
| prod_10ns_6UON_2fs300K                 | GADGVGKSAL/HLA-C*08:02 | explicit_cuda_10ns             |   10000   |        10   |        299.997  |             41.9   |
| 6UON_GADGVGKSAL_HLA-C0802_chainC_1p0ns | GADGVGKSAL/HLA-C*08:02 | explicit_cuda_replicate_screen |    1000   |         1   |        301.096  |             43.1   |
| 6UON_GADGVGKSAL                        | GADGVGKSAL/HLA-C*08:02 | vacuum_smoke                   |       0.1 |        10   |         28.0727 |              0.417 |
| 6VRN_HMTEVVRHC                         | HMTEVVRHC/HLA-A*02:01  | vacuum_smoke                   |       0.1 |        10   |         30.9585 |              0.267 |

Completed primary evidence: `GADGVGKSAL / HLA-C*08:02` reached 10,000 ps. `HMTEVVRHC / HLA-A*02:01` was still partial at the latest sync, at 10000.0 ps.

## 3. Which Simulations Are Still Partial

- HMTEVVRHC primary 10 ns OpenMM CUDA run: latest synced time 10000.0 ps, temperature 299.89 K, speed 21.0 ns/day.
- HMTEVVRHC alternate replicate watcher should start after the primary exits; do not score replicate consistency until those DCDs are present.

## 4. QC Status

| run_id                                 | candidate              | condition                      |   total_time_ps |   peptide_rmsd_final_nm |   peptide_com_drift_final_nm | qc_status   |
|:---------------------------------------|:-----------------------|:-------------------------------|----------------:|------------------------:|-----------------------------:|:------------|
| prod_10ns_6VRN_1fs300K                 | HMTEVVRHC/HLA-A*02:01  | explicit_cuda_10ns             |         10000   |              0.0651099  |                   0.404641   | ok          |
| 6VRM_HMTEVVRHC_HLA-A0201_chainP_0p5ns  | HMTEVVRHC/HLA-A*02:01  | explicit_cuda_replicate_screen |           500   |              0.0378221  |                   0.043586   | ok          |
| prod_10ns_6VRN_1fs300K                 | HMTEVVRHC/HLA-A*02:01  | explicit_cuda_10ns             |         10000   |              0.0651099  |                   0.404641   | ok          |
| prod_10ns_6UON_2fs300K                 | GADGVGKSAL/HLA-C*08:02 | explicit_cuda_10ns             |         10000   |              0.138828   |                   0.659816   | ok          |
| 6UON_GADGVGKSAL_HLA-C0802_chainC_1p0ns | GADGVGKSAL/HLA-C*08:02 | explicit_cuda_replicate_screen |          1000   |              0.111929   |                   0.243361   | ok          |
| 6UON_GADGVGKSAL                        | GADGVGKSAL/HLA-C*08:02 | vacuum_smoke                   |             0.1 |              0.0100448  |                   0.00478767 | ok          |
| 6VRN_HMTEVVRHC                         | HMTEVVRHC/HLA-A*02:01  | vacuum_smoke                   |             0.1 |              0.00476196 |                   0.00222067 | ok          |

GADGVGKSAL primary peptide RMSD final: 0.139 nm. Peptide COM drift final: 0.660 nm. This is compatible with a stable pilot, not a definitive binding claim.

## 5. pMHC Stability Evidence

| run_id                                 | candidate              |   max_contact_occupancy |   sum_contact_occupancy |
|:---------------------------------------|:-----------------------|------------------------:|------------------------:|
| prod_10ns_6VRN_1fs300K                 | HMTEVVRHC/HLA-A*02:01  |                   0.97  |                   6.735 |
| prod_10ns_6VRN_1fs300K                 | HMTEVVRHC/HLA-A*02:01  |                   1     |                   6.23  |
| 6VRM_HMTEVVRHC_HLA-A0201_chainP_0p5ns  | HMTEVVRHC/HLA-A*02:01  |                   1     |                   7.5   |
| 6VRM_HMTEVVRHC_HLA-A0201_chainP_0p5ns  | HMTEVVRHC/HLA-A*02:01  |                   1     |                   7.5   |
| prod_10ns_6VRN_1fs300K                 | HMTEVVRHC/HLA-A*02:01  |                   0.97  |                   6.735 |
| prod_10ns_6VRN_1fs300K                 | HMTEVVRHC/HLA-A*02:01  |                   1     |                   6.23  |
| prod_10ns_6UON_2fs300K                 | GADGVGKSAL/HLA-C*08:02 |                   0.995 |                   6.525 |
| prod_10ns_6UON_2fs300K                 | GADGVGKSAL/HLA-C*08:02 |                   1     |                   7.515 |
| 6UON_GADGVGKSAL_HLA-C0802_chainC_1p0ns | GADGVGKSAL/HLA-C*08:02 |                   1     |                   6.7   |
| 6UON_GADGVGKSAL_HLA-C0802_chainC_1p0ns | GADGVGKSAL/HLA-C*08:02 |                   1     |                   8.3   |
| 6UON_GADGVGKSAL                        | GADGVGKSAL/HLA-C*08:02 |                   1     |                   7     |
| 6UON_GADGVGKSAL                        | GADGVGKSAL/HLA-C*08:02 |                   1     |                   7     |
| 6VRN_HMTEVVRHC                         | HMTEVVRHC/HLA-A*02:01  |                   1     |                   7.2   |
| 6VRN_HMTEVVRHC                         | HMTEVVRHC/HLA-A*02:01  |                   1     |                   6.1   |

GADGVGKSAL anchor contacts remained high in the completed primary and 1 ns replicate screens. HMTEVVRHC anchor evidence currently comes only from smoke/partial status until full trajectory analysis is available.

## 6. TCR-Peptide Recognition Evidence

| run_id                                 | candidate              | condition                      |   TCR_recognition_score |   tcr_peptide_contacts_tail | MD_evidence_label       |
|:---------------------------------------|:-----------------------|:-------------------------------|------------------------:|----------------------------:|:------------------------|
| prod_10ns_6VRN_1fs300K                 | HMTEVVRHC/HLA-A*02:01  | explicit_cuda_10ns             |                  1      |                       63.81 | MD_VERY_STRONG          |
| 6VRM_HMTEVVRHC_HLA-A0201_chainP_0p5ns  | HMTEVVRHC/HLA-A*02:01  | explicit_cuda_replicate_screen |                  1      |                       70.5  | MD_VERY_STRONG          |
| prod_10ns_6VRN_1fs300K                 | HMTEVVRHC/HLA-A*02:01  | explicit_cuda_10ns             |                  1      |                       63.81 | MD_VERY_STRONG          |
| prod_10ns_6UON_2fs300K                 | GADGVGKSAL/HLA-C*08:02 | explicit_cuda_10ns             |                  0.4404 |                       22.02 | MD_MODERATE             |
| 6UON_GADGVGKSAL_HLA-C0802_chainC_1p0ns | GADGVGKSAL/HLA-C*08:02 | explicit_cuda_replicate_screen |                  0.504  |                       25.2  | MD_MODERATE             |
| 6UON_GADGVGKSAL                        | GADGVGKSAL/HLA-C*08:02 | vacuum_smoke                   |                  0.54   |                       27    | MD_INSUFFICIENT_RUNTIME |
| 6VRN_HMTEVVRHC                         | HMTEVVRHC/HLA-A*02:01  | vacuum_smoke                   |                  1      |                       73.5  | MD_INSUFFICIENT_RUNTIME |

The GADGVGKSAL primary run has persistent TCR-peptide contact signal in the audit tables, but this remains structural compatibility evidence. It does not prove T-cell activation.

## 7. Mutant-WT Or Decoy Counterfactual Evidence

No WT, scrambled decoy, or same-HLA positive-control trajectory has been discovered yet. Mutant-specific recognition, specificity, and WT cross-reactivity claims are therefore blocked.

## 8. Replicate Consistency

GADGVGKSAL has one completed 1 ns alternate screen and one 10 ns primary run. The replicate is useful as an early instability filter only. A serious decision needs 3x10 ns for mutant, WT, and scrambled controls before any 50-100 ns escalation.

## 9. Integration With CROSS-Neo Predictions

- Joined prediction rows with MD candidate keys: 484
- MD-supported model-high rows: 182
- model-high but MD-low-or-insufficient rows: 0
- model-low but MD-high rows: 48

Top MD-supported rows:

| row_id     | source_dataset   | peptide   | hla_4digit   |   label |    score | model_name                | split_name           | MD_evidence_label   |   MD_evidence_score |   wetlab_priority_score_evidence_adjusted |
|:-----------|:-----------------|:----------|:-------------|--------:|---------:|:--------------------------|:---------------------|:--------------------|--------------------:|------------------------------------------:|
| CNV0_01149 | NEPdb            | HMTEVVRHC | HLA-A*02:01  |       1 | 0.802273 | v2_cf_plm_rf_hla_ranknorm | source_heldout_NEPdb | MD_STRONG           |            0.743344 |                                 nan       |
| CNV0_01150 | NEPdb            | HMTEVVRHC | HLA-A*02:01  |       1 | 0.802273 | v2_cf_plm_rf_hla_ranknorm | source_heldout_NEPdb | MD_STRONG           |            0.743344 |                                 nan       |
| CNV0_01151 | NEPdb            | HMTEVVRHC | HLA-A*02:01  |       1 | 0.802273 | v2_cf_plm_rf_hla_ranknorm | source_heldout_NEPdb | MD_STRONG           |            0.743344 |                                   1.82549 |
| CNV0_01152 | NEPdb            | HMTEVVRHC | HLA-A*02:01  |       1 | 0.802273 | v2_cf_plm_rf_hla_ranknorm | source_heldout_NEPdb | MD_STRONG           |            0.743344 |                                 nan       |
| CNV0_01172 | NEPdb            | HMTEVVRHC | HLA-A*02:01  |       1 | 0.802273 | v2_cf_plm_rf_hla_ranknorm | source_heldout_NEPdb | MD_STRONG           |            0.743344 |                                 nan       |
| CNV0_01176 | NEPdb            | HMTEVVRHC | HLA-A*02:01  |       1 | 0.802273 | v2_cf_plm_rf_hla_ranknorm | source_heldout_NEPdb | MD_STRONG           |            0.743344 |                                 nan       |
| CNV0_01177 | NEPdb            | HMTEVVRHC | HLA-A*02:01  |       1 | 0.802273 | v2_cf_plm_rf_hla_ranknorm | source_heldout_NEPdb | MD_STRONG           |            0.743344 |                                 nan       |
| CNV0_01183 | NEPdb            | HMTEVVRHC | HLA-A*02:01  |       1 | 0.802273 | v2_cf_plm_rf_hla_ranknorm | source_heldout_NEPdb | MD_STRONG           |            0.743344 |                                 nan       |

## 10. Allowed Claims

- Explicit-solvent MD was used as a structural audit layer for selected CROSS-Neo/TCR-aware candidates.
- GADGVGKSAL / HLA-C*08:02 showed moderate structural audit support in a completed 10 ns pilot, including peptide-MHC stability and persistent TCR-peptide contact signal.
- MD can help prioritize wetlab candidates and explain model/model-expert disagreements.

## 11. Forbidden Claims

- Do not claim MD proves immunogenicity, clinical efficacy, or SOTA.
- Do not claim one short trajectory proves binding.
- Do not claim mutant-specific recognition without WT/decoy comparison.
- Do not treat AlphaFold/TCRdock/template structures as ground truth.
- Do not use peptide-only TCR evidence as definitive neoantigen ground truth.

## 12. Wetlab Candidate Recommendation

- Immediate wetlab-priority candidate: `GADGVGKSAL / HLA-C*08:02`, as structural audit support now agrees with several CROSS-Neo/TCR-aware model rows. Use as prioritization evidence only.
- Keep `HMTEVVRHC / HLA-A*02:01` as the highest-priority pending MD candidate because it has strong TCR evidence and an active 10 ns run, but wait for completed trajectory contacts before moving it above GADGVGKSAL on MD evidence.
- Do not advance candidates lacking paired TCR/template support to TCR-pMHC MD claims; use pMHC-only stability screens for those.

## 13. Next Simulation Batch Recommendation

|   priority_rank | row_id     | peptide         | hla_4digit   | md_tier          | recommended_next_batch                                                                           | reason                                                                                                |
|----------------:|:-----------|:----------------|:-------------|:-----------------|:-------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------|
|               1 | CNV0_01151 | HMTEVVRHC       | HLA-A*02:01  | P0_MD_TCR_pMHC   | P0: let current 10ns finish, sync DCD, rerun full analysis, then launch WT/decoy 3x10ns          | highest TCR evidence and live state is stable, but trajectory-derived contacts are pending            |
|               2 | CNV0_01132 | GADGVGKSAL      | HLA-C*08:02  | P0_MD_TCR_pMHC   | P0: add WT, scrambled decoy, same-HLA positive control; run 3x10ns before any 50-100ns promotion | 10ns mutant TCR-pMHC completed with moderate structural support, but no counterfactual control exists |
|               3 | CNV0_01388 | MAWSLGVLVALPFPL | HLA-B*40:01  | P1_MD_pMHC_bulge | P1: pMHC-only 3x10ns bulge/stability screen; no TCR-recognition claim                            | long class-I peptide or conformation uncertainty without paired TCR                                   |
|               4 | CNV0_01384 | MAWSLGVLVALPFPL | HLA-A*02:01  | P1_MD_pMHC_bulge | P1: pMHC-only 3x10ns bulge/stability screen; no TCR-recognition claim                            | long class-I peptide or conformation uncertainty without paired TCR                                   |
|               5 | CNV0_01146 | FSLVFLVYSVFKNNV | HLA-A*11:01  | P1_MD_pMHC_bulge | P1: pMHC-only 3x10ns bulge/stability screen; no TCR-recognition claim                            | long class-I peptide or conformation uncertainty without paired TCR                                   |
|               6 | CNV0_01163 | MAWSLGVLVALPFPL | HLA-B*18:01  | P1_MD_pMHC_bulge | P1: pMHC-only 3x10ns bulge/stability screen; no TCR-recognition claim                            | long class-I peptide or conformation uncertainty without paired TCR                                   |
|               7 | CNV0_01164 | MAWSLGVLVALPFPL | HLA-A*25:01  | P1_MD_pMHC_bulge | P1: pMHC-only 3x10ns bulge/stability screen; no TCR-recognition claim                            | long class-I peptide or conformation uncertainty without paired TCR                                   |
|               8 | CNV0_01276 | MSSLAATTFHWKKCR | HLA-B*07:02  | P1_MD_pMHC_bulge | P1: pMHC-only 3x10ns bulge/stability screen; no TCR-recognition claim                            | long class-I peptide or conformation uncertainty without paired TCR                                   |
|               9 | CNV0_00930 | MSSLAATTFHWKKCR | HLA-A*68:01  | P1_MD_pMHC_bulge | P1: pMHC-only 3x10ns bulge/stability screen; no TCR-recognition claim                            | long class-I peptide or conformation uncertainty without paired TCR                                   |
|              10 | CNV0_00931 | MSSLAATTFHWKKCR | HLA-B*35:03  | P1_MD_pMHC_bulge | P1: pMHC-only 3x10ns bulge/stability screen; no TCR-recognition claim                            | long class-I peptide or conformation uncertainty without paired TCR                                   |
|              11 | CNV0_01391 | FSLVFLVYSVFKNNV | HLA-B*38:01  | P1_MD_pMHC_bulge | P1: pMHC-only 3x10ns bulge/stability screen; no TCR-recognition claim                            | long class-I peptide or conformation uncertainty without paired TCR                                   |
|              12 | CNV0_01144 | FSLVFLVYSVFKNNV | HLA-B*52:01  | P1_MD_pMHC_bulge | P1: pMHC-only 3x10ns bulge/stability screen; no TCR-recognition claim                            | long class-I peptide or conformation uncertainty without paired TCR                                   |

Minimum serious batch per candidate: mutant pMHC, WT pMHC, mutant TCR-pMHC if paired TCR/template exists, WT TCR-pMHC, scrambled peptide negative control, same/similar-HLA positive control, 3 replicates each, 10 ns each. Promote only the strongest and most consistent candidates to 50-100 ns.

## Exact Rerun Commands

```bash
# Update live status and regenerate all currently possible MD outputs.
python project/scripts/cross_neo_md/00_parse_openmm_status.py
python project/scripts/cross_neo_md/01_annotate_complex.py
python project/scripts/cross_neo_md/02_md_qc.py
python project/scripts/cross_neo_md/03_analyze_pmhc_contacts.py
python project/scripts/cross_neo_md/04_analyze_tcr_contacts.py
python project/scripts/cross_neo_md/05_counterfactual_md_analysis.py
python project/scripts/cross_neo_md/06_replicate_consistency.py
python project/scripts/cross_neo_md/07_md_evidence_score.py
python project/scripts/cross_neo_md/08_integrate_md_with_cross_neo.py
python project/scripts/cross_neo_md/09_generate_md_figures.py
python project/scripts/cross_neo_md/11_generate_md_extra_visuals.py
python project/scripts/cross_neo_md/10_build_md_visual_dossier.py
python project/scripts/cross_neo_md/12_write_md_decision_report.py

# While 6VRN is running, sync state only.
OUT=project/results/cross_neo_md_audit_2026_05_10/remote_sync/pod1_thca_neo_bayesian_aux
mkdir -p "$OUT"
ssh -i ~/.runpod/ssh/RunPod-Key-Go -p 20878 root@135.84.176.142 \
  'cd /workspace/openmm_pilot_10ns_package && tar --exclude=trajectory.dcd --exclude=final.chk -czf - prod_10ns_6VRN_1fs300K' \
  | tar -xzf - -C "$OUT"

# After 6VRN finishes, sync full trajectory and rerun the analysis stack above.
rsync -avP -e 'ssh -i ~/.runpod/ssh/RunPod-Key-Go -p 20878' \
  root@135.84.176.142:/workspace/openmm_pilot_10ns_package/prod_10ns_6VRN_1fs300K/ \
  project/results/cross_neo_md_audit_2026_05_10/remote_sync/pod1_thca_neo_bayesian_aux/prod_10ns_6VRN_1fs300K/
```
