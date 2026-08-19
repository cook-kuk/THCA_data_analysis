# P0 TCR-pMHC Structure Pilot

Existing PDB structures were fetched for P0 MD candidates and screened for exact peptide chains plus nearby MHC, beta-2 microglobulin, and TCR chains.

- P0 registry PDB rows: 9
- Extracted pilot complexes: 11
- Ready TCR-pMHC PDB pilots: 10

## Best Ready Pilot Complexes

| peptide | HLA | PDB | peptide chain | selected chains | TCR-peptide contacts | extracted PDB |
|---|---|---|---|---|---:|---|
| HMTEVVRHC | HLA-A*02:01 | 6VRN | P | A|B|D|E|P | 254 | `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/md_escalation/p0_structures/pilot_complexes/6VRN_HMTEVVRHC_HLA-A0201_chainP.pdb` |
| HMTEVVRHC | HLA-A*02:01 | 6VRM | P | A|B|D|E|P | 214 | `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/md_escalation/p0_structures/pilot_complexes/6VRM_HMTEVVRHC_HLA-A0201_chainP.pdb` |
| HMTEVVRHC | HLA-A*02:01 | 6VQO | P | A|B|D|E|P | 200 | `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/md_escalation/p0_structures/pilot_complexes/6VQO_HMTEVVRHC_HLA-A0201_chainP.pdb` |
| GADGVGKSAL | HLA-C*08:02 | 6UON | F | D|E|F|G|H | 90 | `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/md_escalation/p0_structures/pilot_complexes/6UON_GADGVGKSAL_HLA-C0802_chainF.pdb` |
| GADGVGKSAL | HLA-C*08:02 | 6UON | C | A|B|C|I|J | 73 | `/home/seungho/personal/THCA_data_analysis/project/results/cross_neo_v2_sota_sprint_2026_05_09/tcr_extension/md_escalation/p0_structures/pilot_complexes/6UON_GADGVGKSAL_HLA-C0802_chainC.pdb` |

## Boundary

These are experimentally solved complexes, so they are the right first MD pilots. They still need protonation/repair/force-field preparation before production MD.
