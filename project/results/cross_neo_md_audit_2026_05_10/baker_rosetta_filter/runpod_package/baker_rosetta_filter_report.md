# Baker/Rosetta Cheap Structural Filter

## Purpose

This is the structural layer after the DL-first screen and before expensive MD. It prepares Rosetta/ProteinMPNN/RF-style jobs and runs a static contact fallback when Rosetta is not installed.

## Tool Availability

| environment   | host                  | tool_type     | tool      | path_or_version                                             |
|:--------------|:----------------------|:--------------|:----------|:------------------------------------------------------------|
| local         | YG1-cook-vm           | command       | python    | /home/seungho/personal/THCA_data_analysis/.venv/bin/python  |
| local         | YG1-cook-vm           | command       | python3   | /home/seungho/personal/THCA_data_analysis/.venv/bin/python3 |
| local         | YG1-cook-vm           | python_module | Bio       | 1.85                                                        |
| local         | YG1-cook-vm           | python_module | openmm    | 8.5.1                                                       |
| local         | YG1-cook-vm           | python_module | pennylane | 0.44.1                                                      |
| local         | YG1-cook-vm           | python_module | torch     | 2.11.0+cpu                                                  |
| runpod        | thca-neo-bayesian-aux | command       | python    | /usr/bin/python                                             |
| runpod        | thca-neo-bayesian-aux | command       | python3   | /usr/bin/python3                                            |
| runpod        | thca-neo-bayesian-aux | python_module | openmm    | 8.5.1                                                       |
| runpod        | thca-neo-bayesian-aux | python_module | torch     | 2.4.1+cu124                                                 |

## RunPod Status

| pod_id         | pod_name               | desired_status   | ssh_ip         |   ssh_port | ssh_available   |
|:---------------|:-----------------------|:-----------------|:---------------|-----------:|:----------------|
| r8v801csxxscw0 | thca-neo-bayesian-aux  | RUNNING          | 135.84.176.142 |      20878 | True            |
| uvp9i2r9s6l85y | thca-spark-dm-a6000-v4 | EXITED           |                |            | False           |
| ytnsyachant0a4 | thca-img-dm1-a6000     | RUNNING          | 194.68.245.88  |      22105 | True            |

## Manifest Summary

| candidate_key          | control_type                         | complex_kind   |   n |
|:-----------------------|:-------------------------------------|:---------------|----:|
| GADGVGKSAL|HLA-C*08:02 | mutant                               | TCR-pMHC       |   1 |
| GADGVGKSAL|HLA-C*08:02 | mutant                               | pMHC           |   1 |
| GADGVGKSAL|HLA-C*08:02 | same_or_similar_hla_positive_control | pMHC           |   1 |
| HMTEVVRHC|HLA-A*02:01  | mutant                               | TCR-pMHC       |   1 |
| HMTEVVRHC|HLA-A*02:01  | mutant                               | pMHC           |   1 |
| HMTEVVRHC|HLA-A*02:01  | same_or_similar_hla_positive_control | TCR-pMHC       |   1 |
| HMTEVVRHC|HLA-A*02:01  | same_or_similar_hla_positive_control | pMHC           |   1 |

## Fallback Static Interface Scores

| filter_id        | sequence   | candidate_key          | control_type                         | complex_kind   |   fallback_structural_score |   pmhc_residue_pair_contacts_4A |   tcr_peptide_residue_pair_contacts_4A |   cross_chain_clashes_2A |
|:-----------------|:-----------|:-----------------------|:-------------------------------------|:---------------|----------------------------:|--------------------------------:|---------------------------------------:|-------------------------:|
| BAKERFILTER_0001 | HMTEVVRHC  | HMTEVVRHC|HLA-A*02:01  | mutant                               | TCR-pMHC       |                        1    |                              40 |                                     15 |                        0 |
| BAKERFILTER_0003 | ELAGIGILTV | HMTEVVRHC|HLA-A*02:01  | same_or_similar_hla_positive_control | TCR-pMHC       |                        1    |                              43 |                                     16 |                        0 |
| BAKERFILTER_0000 | HMTEVVRHC  | HMTEVVRHC|HLA-A*02:01  | mutant                               | pMHC           |                        0.95 |                              40 |                                      0 |                        0 |
| BAKERFILTER_0002 | ELAGIGILTV | HMTEVVRHC|HLA-A*02:01  | same_or_similar_hla_positive_control | pMHC           |                        0.95 |                              43 |                                      0 |                        0 |
| BAKERFILTER_0004 | GADGVGKSAL | GADGVGKSAL|HLA-C*08:02 | mutant                               | pMHC           |                        0.95 |                              35 |                                      0 |                        0 |
| BAKERFILTER_0006 | GADGVGKSA  | GADGVGKSAL|HLA-C*08:02 | same_or_similar_hla_positive_control | pMHC           |                        0.95 |                              35 |                                      0 |                        0 |
| BAKERFILTER_0005 | GADGVGKSAL | GADGVGKSAL|HLA-C*08:02 | mutant                               | TCR-pMHC       |                        0.86 |                              35 |                                      6 |                        0 |

## Interpretation Boundary

- Rosetta/RF/ProteinMPNN scores are structural triage signals, not immunogenicity proof.
- Static fallback contacts are weaker than Rosetta energies and much weaker than replicated MD.
- TCR-discordant main-model positives should be treated as false-positive audit/control candidates.
- MD escalation should remain limited to candidates that survived DL uncertainty and TCR concordance.
