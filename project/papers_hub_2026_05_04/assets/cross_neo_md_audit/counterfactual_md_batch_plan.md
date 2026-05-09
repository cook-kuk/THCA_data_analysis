# Counterfactual MD Batch Design

## Boundary

This is a simulation-design manifest. It does not establish mutant-specific TCR recognition. WT rows with inferred sequences are marked for manual confirmation before any run.

## Summary

- candidate/control jobs: 56
- ready or near-ready rows: 14
- blocked rows: 26
- candidate sequence rows: 12

## Immediate P0 Plan

- `HMTEVVRHC / HLA-A*02:01`: finish current 10 ns, sync DCD, rerun full analysis, then add WT/decoy controls once WT is confirmed.
- `GADGVGKSAL / HLA-C*08:02`: add WT confirmation, anchor-preserved decoy, and same-HLA positive-control selection before long MD.

## RunPod Submission Template

Use this only after structure preparation has produced PDB inputs for each manifest row.

```bash
python run_openmm_pilot.py --input INPUT.pdb --outdir OUTDIR --mode explicit --platform CUDA --ns 10 --report-steps 50000 --minimize-iterations 1000 --timestep-fs 2.0 --temperature-k 300
```

## Outputs

- `counterfactual_md_batch_manifest.tsv`
- `candidate_control_sequences.tsv`
- `candidate_decoy_controls.fasta`
