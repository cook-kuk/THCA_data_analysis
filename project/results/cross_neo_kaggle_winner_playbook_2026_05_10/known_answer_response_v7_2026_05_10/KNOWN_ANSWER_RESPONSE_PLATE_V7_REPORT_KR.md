# CROSS-Neo known-answer response visualization v7

Generated: 2026-05-10T20:34:53

## Core result

Using `stress_guarded_discovery_score` to rank candidates, the top-96 plate contains 91/96 known positive labels (94.8%). Top-10 is 10/10 and top-24 is 22/24.

## Boundary

This is not a prospective wetlab result. It is a known-answer visualization built from existing labels, useful for demo figures, positive-control board design, and failure-well inspection.

## Score board

| score_column                        |   top10_hits |   top10_precision |   top24_hits |   top24_precision |   top48_hits |   top48_precision |   top96_hits |   top96_precision | known_answer_boundary                                               |
|:------------------------------------|-------------:|------------------:|-------------:|------------------:|-------------:|------------------:|-------------:|------------------:|:--------------------------------------------------------------------|
| stress_guarded_discovery_score      |           10 |             1.000 |           22 |             0.917 |           46 |             0.958 |           91 |             0.948 | retrospective known-label metric; not prospective wetlab validation |
| bma_v2_discovery_score              |           10 |             1.000 |           22 |             0.917 |           45 |             0.938 |           89 |             0.927 | retrospective known-label metric; not prospective wetlab validation |
| finetuned_experiment_priority_score |           10 |             1.000 |           22 |             0.917 |           44 |             0.917 |           84 |             0.875 | retrospective known-label metric; not prospective wetlab validation |
| finetuned_score_booster_prob        |            9 |             0.900 |           20 |             0.833 |           42 |             0.875 |           83 |             0.865 | retrospective known-label metric; not prospective wetlab validation |
| impact_portfolio_score              |            5 |             0.500 |           15 |             0.625 |           26 |             0.542 |           67 |             0.698 | retrospective known-label metric; not prospective wetlab validation |
| v4_claim_unlock_score               |            5 |             0.500 |           17 |             0.708 |           17 |             0.354 |           58 |             0.604 | retrospective known-label metric; not prospective wetlab validation |

## Top false-positive wells

| well   |   plate_rank | candidate_id   | peptide     | hla_allele_4digit   | source_name   | leakage_risk_level   |   stress_guarded_discovery_score | known_answer_response_call   |
|:-------|-------------:|:---------------|:------------|:--------------------|:--------------|:---------------------|---------------------------------:|:-----------------------------|
| H2     |           16 | CNV0_00027     | GNNDVKEDP   | HLA-C*07:02         | CEDAR         | high                 |                            0.721 | KNOWN_NONRESPONDER           |
| E3     |           21 | CNV0_00006     | VKEDPKWEF   | HLA-C*07:02         | CEDAR         | high                 |                            0.716 | KNOWN_NONRESPONDER           |
| E7     |           53 | CNV0_00003     | LYGLLLEML   | HLA-A*02:01         | CEDAR         | high                 |                            0.694 | KNOWN_NONRESPONDER           |
| H7     |           56 | CNV0_00711     | RREPPHLARNF | HLA-C*07:02         | CEDAR         | high                 |                            0.693 | KNOWN_NONRESPONDER           |
| E9     |           69 | CNV0_00064     | LLAGLVSLL   | HLA-A*02:01         | CEDAR         | high                 |                            0.686 | KNOWN_NONRESPONDER           |

## Top responder wells

| well   |   plate_rank | candidate_id   | peptide         | hla_allele_4digit   | source_name   |   stress_guarded_discovery_score | known_answer_response_call   |
|:-------|-------------:|:---------------|:----------------|:--------------------|:--------------|---------------------------------:|:-----------------------------|
| A1     |            1 | CNV0_00169     | AALEDTLAETEAR   | HLA-B*37:01         | CEDAR         |                            0.763 | KNOWN_RESPONDER              |
| B1     |            2 | CNV0_00804     | SRSYTSGPGSRISSS | HLA-B*27:09         | CEDAR         |                            0.761 | KNOWN_RESPONDER              |
| C1     |            3 | CNV0_00688     | YYYGIKDLATVFF   | HLA-A*24:02         | CEDAR         |                            0.760 | KNOWN_RESPONDER              |
| D1     |            4 | CNV0_00837     | EQFLDGDGWTSR    | HLA-A*34:02         | CEDAR         |                            0.751 | KNOWN_RESPONDER              |
| E1     |            5 | CNV0_00874     | FPSEYLSSHLEA    | HLA-B*54:01         | CEDAR         |                            0.751 | KNOWN_RESPONDER              |
| F1     |            6 | CNV0_00694     | AEIEGLKGQRASLE  | HLA-B*35:02         | CEDAR         |                            0.750 | KNOWN_RESPONDER              |
| G1     |            7 | CNV0_00494     | QEMASVEAAVSL    | HLA-A*02:01         | CEDAR         |                            0.749 | KNOWN_RESPONDER              |
| H1     |            8 | CNV0_00771     | YYGTGETFLYTF    | HLA-A*24:02         | CEDAR         |                            0.740 | KNOWN_RESPONDER              |
| A2     |            9 | CNV0_00088     | YAYDNFGVLGL     | HLA-C*03:03         | CEDAR         |                            0.738 | KNOWN_RESPONDER              |
| B2     |           10 | CNV0_00120     | AEEEASAVSTAA    | HLA-B*45:01         | CEDAR         |                            0.737 | KNOWN_RESPONDER              |
| C2     |           11 | CNV0_00784     | ETVNALISDQKL    | HLA-A*68:02         | CEDAR         |                            0.737 | KNOWN_RESPONDER              |
| D2     |           12 | CNV0_00346     | PASTGLGWGSGL    | HLA-B*58:01         | CEDAR         |                            0.725 | KNOWN_RESPONDER              |
| E2     |           13 | CNV0_00245     | ATSPHLESLLK     | HLA-A*11:01         | CEDAR         |                            0.724 | KNOWN_RESPONDER              |
| F2     |           14 | CNV0_00306     | TGKEKVTSGSTTTTR | HLA-A*33:01         | CEDAR         |                            0.722 | KNOWN_RESPONDER              |
| G2     |           15 | CNV0_00733     | KPTSSKSSSEATL   | HLA-B*07:02         | CEDAR         |                            0.722 | KNOWN_RESPONDER              |
| H2     |           16 | CNV0_00027     | GNNDVKEDP       | HLA-C*07:02         | CEDAR         |                            0.721 | KNOWN_NONRESPONDER           |
| A3     |           17 | CNV0_00455     | NAPMTLEEF       | HLA-C*01:02         | CEDAR         |                            0.721 | KNOWN_RESPONDER              |
| B3     |           18 | CNV0_00487     | QELNELSAISL     | HLA-B*40:01         | CEDAR         |                            0.718 | KNOWN_RESPONDER              |
| C3     |           19 | CNV0_00822     | RSTIITLYNGAFDSS | HLA-A*01:01         | CEDAR         |                            0.716 | KNOWN_RESPONDER              |
| D3     |           20 | CNV0_02473     | KEFEDDIINW      | HLA-B*44:03         | ITSNdb_main   |                            0.716 | KNOWN_RESPONDER              |
