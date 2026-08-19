# P0 OpenMM Pilot Smoke Test

Local CPU smoke tests used dry/vacuum dynamics only to verify that the packaged OpenMM runner can create trajectories, state logs, and downstream MDTraj contact summaries. These are not production MD results.

- Smoke runs found: 2
- Successful runs: 2

| smoke id | status | steps | final peptide RMSD nm | final pMHC contacts | final TCR contacts |
|---|---|---:|---:|---:|---:|
| 6UON_GADGVGKSAL | ok | 200 | 0.0144 | 245 | 58 |
| 6VRN_HMTEVVRHC | ok | 200 | 0.0091 | 383 | 27 |

Next step: run the same package in explicit-solvent CUDA mode on a GPU pod for 10 ns pilots.
