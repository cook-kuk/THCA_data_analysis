# CROSS-Neo-TCR Structure Job Manifest Report

Candidate paired TCR-pMHC jobs emitted: 638

## Tool Availability

| Tool | Local status |
|---|---|
| TCRdock | not_found |
| TCRmodel2 | not_found |
| AlphaFold-Multimer_colabfold_batch | /home/seungho/personal/THCA_data_analysis/.venv/bin/colabfold_batch |
| AlphaFold_python_module | available |
| AlphaFold3 | not_found |
| Boltz | not_found |
| Chai-1 | not_found |
| tFold-TCR | not_found |
| ImmuneBuilder | not_found |

## Priority Summary

| Priority | Tool | Claim status | n jobs |
|---|---|---|---:|
| P0 | blocked_missing_full_tcr_sequence | diagnostic_manifest_only_requires_full_tcr | 507 |
| P1 | blocked_missing_full_tcr_sequence | diagnostic_manifest_only_requires_full_tcr | 131 |

## Claim Boundary

- Jobs are retained even when full TCR or MHC sequences are missing.
- CDR3-only rows are useful for prioritization but are not sufficient for TCR-pMHC structure prediction with most tools.
- Local runnable path currently favors ColabFold/AlphaFold-Multimer only when full chain sequences are available.

## Parse Status

- Manifest jobs checked: 638
- Parsed structure-like outputs: 2
- Interface features are not claimable until actual TCR-pMHC outputs exist.
