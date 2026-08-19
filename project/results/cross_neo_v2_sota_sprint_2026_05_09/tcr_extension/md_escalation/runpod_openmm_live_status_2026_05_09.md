# RunPod OpenMM Live Status

Timestamp: 2026-05-09 22:12 KST.

## Primary 10 ns pilots

| candidate | pod | pid | condition | latest progress | latest temperature | speed | status |
|---|---|---:|---|---:|---:|---:|---|
| HMTEVVRHC / HLA-A*02:01 / 6VRN | thca-neo-bayesian-aux | 1956 | explicit solvent, CUDA, 1 fs, 300 K | 50 ps | 299.995 K | pending after first report | running |
| GADGVGKSAL / HLA-C*08:02 / 6UON | thca-img-dm1-a6000 | 1047 | explicit solvent, CUDA, 2 fs, 300 K | 300 ps | 300.245 K | 42.2 ns/day | running |

## Queued follow-up screens

| queue | pod | watcher pid | starts after | payload |
|---|---|---:|---|---|
| HMTEVVRHC replicate screens | thca-neo-bayesian-aux | 2213 | primary 6VRN pid 1956 exits | 7 alternate solved TCR-pMHC complexes, 0.5 ns each |
| GADGVGKSAL replicate screen | thca-img-dm1-a6000 | 1428 | primary 6UON pid 1047 exits | alternate 6UON complex, 1 ns |

## Interpretation

Both primary simulations have passed CUDA explicit-solvent startup. The GADGVGKSAL run already has multiple stable reports at 300 K. The HMTEVVRHC run was restarted with the corrected 1 fs step count and has reached its first 50 ps report at 300 K.

This is still MD pilot evidence, not immunogenicity proof. Final interpretation requires completed trajectories plus peptide RMSD, peptide-HLA contacts, and TCR-peptide contact analysis.
