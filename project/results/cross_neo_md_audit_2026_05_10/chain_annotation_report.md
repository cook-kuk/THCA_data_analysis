# Chain Annotation Report

Chains are assigned by sequence and length heuristics: exact expected short peptide sequence, longest HLA-like chain, beta2m-like 80-120 residue chain, and 150-260 residue TCR candidate chains.

| run_id | candidate | expected_peptide | peptide_sequence_from_structure | peptide_match | peptide_chains | mhc_heavy_chain | beta2m_chain | tcr_candidate_chains | manual_review_flag |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| prod_10ns_6VRN_1fs300K | HMTEVVRHC/HLA-A*02:01 | HMTEVVRHC | HMTEVVRHC | True | E | A | B | C|D |  |
| prod_10ns_6UON_2fs300K | GADGVGKSAL/HLA-C*08:02 | GADGVGKSAL | GADGVGKSAL | True | E | C | D | A|B |  |
| 6UON_GADGVGKSAL_HLA-C0802_chainC_1p0ns | GADGVGKSAL/HLA-C*08:02 | GADGVGKSAL | GADGVGKSAL | True | E | C | D | A|B |  |
| 6UON_GADGVGKSAL | GADGVGKSAL/HLA-C*08:02 | GADGVGKSAL | GADGVGKSAL | True | E | C | D | A|B |  |
| 6VRN_HMTEVVRHC | HMTEVVRHC/HLA-A*02:01 | HMTEVVRHC | HMTEVVRHC | True | E | A | B | C|D |  |

Mutation residue positions are marked unknown unless a validated mutant/WT mapping is supplied. Do not claim mutant-site contact specificity without that mapping.
