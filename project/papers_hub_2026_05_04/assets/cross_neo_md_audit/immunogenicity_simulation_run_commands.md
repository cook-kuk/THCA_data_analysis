# Immunogenicity Simulation Escalation Commands

These commands are launch scaffolds. They do not prove immunogenicity and should be run only when the current queue/resources are intentionally allocated.

## Refresh status and reports

```bash
python project/scripts/cross_neo_md/13_watch_runpod_md_and_autoreport.py --once
python project/scripts/cross_neo_md/30_design_immunogenicity_simulation_escalation.py
python project/scripts/cross_neo_md/29_build_high_impact_decision_web.py
```

## Existing queued/ready OpenMM diverse simulation pack

```bash
cd project/results/cross_neo_md_audit_2026_05_10/diverse_simulation_plan
bash run_immediate_ready_batch.sh
bash run_after_hmtevv_sync.sh
```

## Structure-prep refresh before new WT/decoy/TCR jobs

```bash
python project/scripts/cross_neo_v2/13_prepare_tcr_structure_jobs.py
python project/scripts/cross_neo_md/14_prepare_md_counterfactual_batch.py
python project/scripts/cross_neo_md/17_prepare_counterfactual_structure_prep_jobs.py
python project/scripts/cross_neo_md/25_prepare_diverse_md_launch_pack.py
```

## Endpoint energy / MMGBSA placeholder after trajectories finish

```bash
# If AmberTools/MMPBSA.py is installed, run per completed trajectory.
# Otherwise use OpenMM interaction-energy decomposition from parsed frames.
python project/scripts/cross_neo_md/03_analyze_pmhc_contacts.py
python project/scripts/cross_neo_md/04_analyze_tcr_contacts.py
python project/scripts/cross_neo_md/05_counterfactual_md_analysis.py
```

## Expensive hold-only modules

- `alchemical_fep_mutant_to_wt_delta_delta_g`: run only for GADGVGKSAL/HMTEVVRHC after WT mapping and stable 10 ns controls.
- `umbrella_or_steered_md_tcr_unbinding_pmf`: run only after stable TCR-pMHC 10 ns and if multimer/TCR-binding assay decision depends on it.
