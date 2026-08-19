# OpenMM P0 Minimization Sanity Check

Representative P0 TCR-pMHC PDB complexes were repaired with PDBFixer and dry-minimized in OpenMM.

- Tested complexes: 2
- Successful minimizations: 2

| peptide | HLA | PDB | status | atoms | initial kJ/mol | minimized kJ/mol | output |
|---|---|---|---|---:|---:|---:|---|
| GADGVGKSAL | HLA-C*08:02 | 6UON | ok | 12,778 | 124424.8 | -42550.5 | `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/md_escalation/p0_structures/openmm_minimized/6UON_GADGVGKSAL_HLA-C0802_chainF.openmm_minimized.pdb` |
| HMTEVVRHC | HLA-A*02:01 | 6VRN | ok | 12,775 | 987921.7 | 16663.6 | `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/md_escalation/p0_structures/openmm_minimized/6VRN_HMTEVVRHC_HLA-A0201_chainP.openmm_minimized.pdb` |

Boundary: this is not production MD. It confirms the representative structures can be repaired and passed through an OpenMM force-field/minimization workflow.
