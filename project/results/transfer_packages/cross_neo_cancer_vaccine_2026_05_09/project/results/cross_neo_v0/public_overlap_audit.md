# Public / Local Overlap Audit

This audit is intentionally conservative.

Strict rows audited: 89
Local exact peptide-HLA hits: 0
Local exact peptide hits: 1
Local near hits similarity >=0.75: 1

Public predictor training overlap remains unresolved unless each public training corpus is downloaded and row-audited.

| source     | status                       | note                                                                          |
|:-----------|:-----------------------------|:------------------------------------------------------------------------------|
| MHCflurry  | unresolved_public_pretrained | IEDB/MS-ligand/affinity training corpus not row-audited locally               |
| NetMHCpan  | unresolved_public_pretrained | BA/EL public training corpus not row-audited locally                          |
| BigMHC     | unresolved_public_pretrained | presentation/immunogenicity release exists but not downloaded into this audit |
| PRIME      | unresolved_public_pretrained | public immunogenicity training set not row-audited locally                    |
| MixMHCpred | unresolved_public_pretrained | public ligand training corpus not row-audited locally                         |
| CEDAR      | local_train_pool_present     | exact/near overlap against local train pool auditable                         |
| TESLA      | local_train_pool_present     | exact/near overlap against local train pool auditable                         |
| NEPdb      | local_train_pool_present     | exact/near overlap against local train pool auditable                         |
