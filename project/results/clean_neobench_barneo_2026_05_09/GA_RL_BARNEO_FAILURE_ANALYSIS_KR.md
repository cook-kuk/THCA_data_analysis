# GA/RL -> BAR-Neo-NG Failure Analysis KR

## 결론

GA/RL 고득점 후보는 많지만, 실제 clean handoff는 T1 3개뿐이다. 나머지는 exact/near overlap, high leakage, source reuse, 또는 MHC-II separate-benchmark 사유로 막힌다.

## Failure analysis table

| analysis_group        | candidate_id   | ga_rl_track              | peptide         | hla_allele_4digit   | source_name   |   ga_rl_score |   barneo_ng_score | ga_rl_barneo_decision                | analysis_reason               | analysis_overlap_status   | leakage_risk_level   | audit_priority                              |
|:----------------------|:---------------|:-------------------------|:----------------|:--------------------|:--------------|--------------:|------------------:|:-------------------------------------|:------------------------------|:--------------------------|:---------------------|:--------------------------------------------|
| unique_t1_handoff     | CNV0_02407     | BioDarwin_GA_RL_academic | KLMNIQQKL       | HLA-A*02:01         | ITSNdb_main   |      0.740475 |         0.72      | GA_RL_hit_confirmed_T1_by_BAR_Neo_NG | BAR-Neo-NG T1 pass            | low_leakage_or_pass       | low                  | nan                                         |
| unique_t1_handoff     | CNV0_02504     | BioDarwin_GA_RL_academic | LLVDLAEEL       | HLA-A*02:01         | ITSNdb_main   |      0.685061 |         0.72      | GA_RL_hit_confirmed_T1_by_BAR_Neo_NG | BAR-Neo-NG T1 pass            | low_leakage_or_pass       | low                  | nan                                         |
| unique_t1_handoff     | CNV0_02410     | BioDarwin_GA_RL_academic | MLGEQLFPL       | HLA-A*02:01         | ITSNdb_main   |      0.638571 |         0.72      | GA_RL_hit_confirmed_T1_by_BAR_Neo_NG | BAR-Neo-NG T1 pass            | low_leakage_or_pass       | low                  | nan                                         |
| high_ga_blocked_audit | CNV0_02448     | BioDarwin_GA_RL_academic | ILDKVLVHL       | HLA-A*02:01         | ITSNdb_main   |      0.996333 |         0.499141  | GA_RL_hit_watchlist_manual_review    | manual_review_or_source_reuse | non_overlap_but_high_GA   | low                  | B_watchlist_high_GA_manual_review           |
| high_ga_blocked_audit | CNV0_02708     | KG_GA_evolved            | ALDPLLLRI       | HLA-A*02:01         | ITSNdb_Val    |      0.709598 |         0.521087  | GA_RL_hit_watchlist_manual_review    | manual_review_or_source_reuse | non_overlap_but_high_GA   | low                  | B_watchlist_high_GA_manual_review           |
| high_ga_blocked_audit | CNV0_02405     | KG_GA_evolved            | KMIGNHLWV       | HLA-A*02:01         | ITSNdb_main   |      0.460144 |         0.465655  | GA_RL_hit_watchlist_manual_review    | manual_review_or_source_reuse | non_overlap_but_high_GA   | low                  | B_watchlist_high_GA_manual_review           |
| high_ga_blocked_audit | CNV0_00688     | BioDarwin_GA_RL_academic | YYYGIKDLATVFF   | HLA-A*24:02         | CEDAR         |      1        |         0.122851  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_01132     | KG_GA_evolved            | GADGVGKSAL      | HLA-C*08:02         | NEPdb         |      1        |         0.0997747 | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00784     | BioDarwin_GA_RL_academic | ETVNALISDQKL    | HLA-A*68:02         | CEDAR         |      0.999559 |         0.111528  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00169     | BioDarwin_GA_RL_academic | AALEDTLAETEAR   | HLA-B*37:01         | CEDAR         |      0.9513   |         0.112326  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00804     | BioDarwin_GA_RL_academic | SRSYTSGPGSRISSS | HLA-B*27:09         | CEDAR         |      0.946788 |         0.11132   | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00874     | BioDarwin_GA_RL_academic | FPSEYLSSHLEA    | HLA-B*54:01         | CEDAR         |      0.939424 |         0.111181  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00791     | BioDarwin_GA_RL_academic | EVIDKNSGGWWYV   | HLA-A*68:02         | CEDAR         |      0.91656  |         0.106931  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00821     | BioDarwin_GA_RL_academic | NADPASHEIW      | HLA-B*53:01         | CEDAR         |      0.911398 |         0.119311  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00771     | BioDarwin_GA_RL_academic | YYGTGETFLYTF    | HLA-A*24:02         | CEDAR         |      0.910634 |         0.12171   | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00245     | BioDarwin_GA_RL_academic | ATSPHLESLLK     | HLA-A*11:01         | CEDAR         |      0.908619 |         0.109555  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00120     | BioDarwin_GA_RL_academic | AEEEASAVSTAA    | HLA-B*45:01         | CEDAR         |      0.906665 |         0.111753  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00704     | BioDarwin_GA_RL_academic | MTEVISSLENANY   | HLA-A*01:01         | CEDAR         |      0.905715 |         0.120643  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00779     | BioDarwin_GA_RL_academic | RTLMENQHW       | HLA-B*58:01         | CEDAR         |      0.90527  |         0.118176  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00812     | BioDarwin_GA_RL_academic | AEIDVIFKDFVNKY  | HLA-B*44:03         | CEDAR         |      0.904865 |         0.108313  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00231     | BioDarwin_GA_RL_academic | AMFGKLMTI       | HLA-A*02:01         | CEDAR         |      0.904159 |         0.113677  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00723     | BioDarwin_GA_RL_academic | EAMNYEGSPIKV    | HLA-A*68:02         | CEDAR         |      0.898748 |         0.108026  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00228     | BioDarwin_GA_RL_academic | ALQDIGKNIYTI    | HLA-A*02:01         | CEDAR         |      0.896869 |         0.124448  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00152     | BioDarwin_GA_RL_academic | EIFDSRGNPTVEV   | HLA-A*68:02         | CEDAR         |      0.89408  |         0.10816   | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00156     | BioDarwin_GA_RL_academic | ALPEVLAVIQV     | HLA-A*02:01         | CEDAR         |      0.892308 |         0.113476  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00820     | BioDarwin_GA_RL_academic | FMMPRIVNV       | HLA-A*02:01         | CEDAR         |      0.891626 |         0.140876  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00138     | BioDarwin_GA_RL_academic | YHEKGRAFL       | HLA-B*15:10         | CEDAR         |      0.890759 |         0.120472  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00852     | BioDarwin_GA_RL_academic | RRIILSGPSGTGKTY | HLA-B*27:05         | CEDAR         |      0.890244 |         0.113505  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00733     | BioDarwin_GA_RL_academic | KPTSSKSSSEATL   | HLA-B*07:02         | CEDAR         |      0.889592 |         0.118144  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00370     | BioDarwin_GA_RL_academic | NMASFVRQL       | HLA-A*02:02         | CEDAR         |      0.88265  |         0.120601  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00187     | BioDarwin_GA_RL_academic | VTENTTGKAATY    | HLA-A*01:01         | CEDAR         |      0.8823   |         0.120911  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00753     | BioDarwin_GA_RL_academic | KLISSDGHEFIVK   | HLA-A*03:01         | CEDAR         |      0.881309 |         0.119157  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00751     | BioDarwin_GA_RL_academic | GSFPENLRHLK     | HLA-A*11:01         | CEDAR         |      0.879285 |         0.108796  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00306     | BioDarwin_GA_RL_academic | TGKEKVTSGSTTTTR | HLA-A*33:01         | CEDAR         |      0.876719 |         0.109701  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00144     | BioDarwin_GA_RL_academic | TPIENMILRY      | HLA-B*35:01         | CEDAR         |      0.87623  |         0.121914  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_02473     | BioDarwin_GA_RL_academic | KEFEDDIINW      | HLA-B*44:03         | ITSNdb_main   |      0.876092 |         0.120539  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00487     | BioDarwin_GA_RL_academic | QELNELSAISL     | HLA-B*40:01         | CEDAR         |      0.875806 |         0.107394  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00906     | BioDarwin_GA_RL_academic | FLLRRIPTL       | HLA-A*02:01         | CEDAR         |      0.873567 |         0.139237  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00317     | BioDarwin_GA_RL_academic | GVADALLYR       | HLA-A*11:01         | CEDAR         |      0.87214  |         0.134054  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00074     | BioDarwin_GA_RL_academic | EASPLSSNKLILR   | HLA-A*33:03         | CEDAR         |      0.871032 |         0.117869  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00835     | BioDarwin_GA_RL_academic | RQMQNINPLTM     | HLA-B*27:09         | CEDAR         |      0.869598 |         0.109677  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00837     | BioDarwin_GA_RL_academic | EQFLDGDGWTSR    | HLA-A*34:02         | CEDAR         |      0.86932  |         0.110001  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00310     | BioDarwin_GA_RL_academic | KEYQETIGQIEL    | HLA-B*40:01         | CEDAR         |      0.868772 |         0.115617  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00498     | BioDarwin_GA_RL_academic | NATTIANHW       | HLA-B*53:01         | CEDAR         |      0.868279 |         0.119559  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00440     | BioDarwin_GA_RL_academic | MTHKLLSRY       | HLA-A*30:02         | CEDAR         |      0.867649 |         0.117511  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00128     | BioDarwin_GA_RL_academic | NAHMDSLQW       | HLA-B*53:01         | CEDAR         |      0.867307 |         0.117689  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00181     | BioDarwin_GA_RL_academic | AEFQDPLGY       | HLA-B*44:03         | CEDAR         |      0.867015 |         0.119745  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00494     | BioDarwin_GA_RL_academic | QEMASVEAAVSL    | HLA-A*02:01         | CEDAR         |      0.865846 |         0.127902  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00088     | BioDarwin_GA_RL_academic | YAYDNFGVLGL     | HLA-C*03:03         | CEDAR         |      0.865313 |         0.125601  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00243     | BioDarwin_GA_RL_academic | RSTHEAVARW      | HLA-B*58:01         | CEDAR         |      0.864846 |         0.119272  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00353     | BioDarwin_GA_RL_academic | QLSWLINRL       | HLA-A*02:02         | CEDAR         |      0.864352 |         0.119039  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00211     | BioDarwin_GA_RL_academic | STKSLQELF       | HLA-B*57:03         | CEDAR         |      0.862873 |         0.119218  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00816     | BioDarwin_GA_RL_academic | FLSEVWNTHTL     | HLA-A*02:01         | CEDAR         |      0.861913 |         0.139722  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00663     | BioDarwin_GA_RL_academic | AEVEGLGKGVA     | HLA-B*40:06         | CEDAR         |      0.861635 |         0.118196  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00411     | BioDarwin_GA_RL_academic | MEVDPIGHVYIF    | HLA-B*44:03         | CEDAR         |      0.861553 |         0.109567  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00722     | BioDarwin_GA_RL_academic | KEDSAVFYEL      | HLA-B*40:01         | CEDAR         |      0.86085  |         0.13258   | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00421     | BioDarwin_GA_RL_academic | MQWVLPKI        | HLA-B*52:01         | CEDAR         |      0.860765 |         0.119278  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00337     | BioDarwin_GA_RL_academic | REMFEVTGL       | HLA-B*40:01         | CEDAR         |      0.859095 |         0.132938  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00742     | BioDarwin_GA_RL_academic | YWNEYGGGLLW     | HLA-A*24:02         | CEDAR         |      0.855774 |         0.13375   | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |
| high_ga_blocked_audit | CNV0_00789     | BioDarwin_GA_RL_academic | LRHELISTL       | HLA-B*27:09         | CEDAR         |      0.854747 |         0.117469  | GA_RL_hit_blocked_from_clean_claim   | high_leakage                  | non_overlap_but_high_GA   | high                 | A_blocked_high_GA_review_leakage_or_overlap |

## Source summary

| analysis_group        | source_name   |   n |
|:----------------------|:--------------|----:|
| high_ga_blocked_audit | CEDAR         |  64 |
| high_ga_blocked_audit | NEPdb         |   4 |
| high_ga_blocked_audit | ITSNdb_main   |   3 |
| high_ga_blocked_audit | ITSNdb_Val    |   1 |
| unique_t1_handoff     | ITSNdb_main   |   3 |

## HLA summary

| analysis_group        | hla_allele_4digit   |   n |
|:----------------------|:--------------------|----:|
| high_ga_blocked_audit | HLA-A*02:01         |  12 |
| high_ga_blocked_audit | HLA-B*40:01         |   7 |
| high_ga_blocked_audit | HLA-A*24:02         |   4 |
| high_ga_blocked_audit | HLA-A*68:02         |   4 |
| high_ga_blocked_audit | HLA-B*27:09         |   4 |
| high_ga_blocked_audit | HLA-B*44:03         |   4 |
| high_ga_blocked_audit | HLA-A*11:01         |   3 |
| high_ga_blocked_audit | HLA-B*53:01         |   3 |
| high_ga_blocked_audit | HLA-A*01:01         |   2 |
| high_ga_blocked_audit | HLA-A*02:02         |   2 |
| high_ga_blocked_audit | HLA-A*30:02         |   2 |
| high_ga_blocked_audit | HLA-B*15:10         |   2 |
| high_ga_blocked_audit | HLA-B*35:02         |   2 |
| high_ga_blocked_audit | HLA-B*58:01         |   2 |
| high_ga_blocked_audit | HLA-A*03:01         |   1 |
| high_ga_blocked_audit | HLA-A*29:02         |   1 |
| high_ga_blocked_audit | HLA-A*33:01         |   1 |
| high_ga_blocked_audit | HLA-A*33:03         |   1 |
| high_ga_blocked_audit | HLA-A*34:02         |   1 |
| high_ga_blocked_audit | HLA-B*07:02         |   1 |
| high_ga_blocked_audit | HLA-B*27:05         |   1 |
| high_ga_blocked_audit | HLA-B*35:01         |   1 |
| high_ga_blocked_audit | HLA-B*37:01         |   1 |
| high_ga_blocked_audit | HLA-B*39:01         |   1 |
| high_ga_blocked_audit | HLA-B*40:06         |   1 |
| high_ga_blocked_audit | HLA-B*45:01         |   1 |
| high_ga_blocked_audit | HLA-B*49:01         |   1 |
| high_ga_blocked_audit | HLA-B*52:01         |   1 |
| high_ga_blocked_audit | HLA-B*54:01         |   1 |
| high_ga_blocked_audit | HLA-B*57:03         |   1 |
| high_ga_blocked_audit | HLA-C*03:03         |   1 |
| high_ga_blocked_audit | HLA-C*08:02         |   1 |
| high_ga_blocked_audit | HLA-C*14:02         |   1 |
| unique_t1_handoff     | HLA-A*02:01         |   3 |

## Interpretation

- `exact_overlap` and `near_overlap` are not clean claim territory.
- `high_leakage` rows are audit/failure cases even when GA/RL score is high.
- `unique_t1_handoff` rows are the only immediate assay-design queue.
