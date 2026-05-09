
# BAR-Neo Abstention Report

## Abstention Counts

- Candidates scored: 2715
- Abstain true: 2703
- Abstain false: 12

## Top Abstention Reasons

| abstention_reason                                                          |   count |
|:---------------------------------------------------------------------------|--------:|
| High leakage risk invalidates clean benchmark claim                        |    2304 |
| High disagreement across predictors                                        |    1235 |
| Low-prevalence source with poor historical top-k performance               |     915 |
| High source-shift risk: similar heldout source had poor top-k survival     |     915 |
| HLA allele underrepresented in benchmark                                   |     274 |
| Low confidence due to sparse method coverage                               |     183 |
| High disagreement between internal anchor and public pretrained predictors |      15 |
| no_abstention_reason                                                       |      12 |

## High Score / Low Confidence Cases

| candidate_id   |   barneo_score |   patient_gated_score |   confidence_score | confidence_bin   | abstain   | abstention_reason_primary                           |
|:---------------|---------------:|----------------------:|-------------------:|:-----------------|:----------|:----------------------------------------------------|
| CNV0_00120     |       0.985514 |              0.363851 |           0.922305 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00410     |       0.984203 |              0.386932 |           0.933276 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00694     |       0.983584 |              0.347737 |           0.713599 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00128     |       0.983515 |              0.401335 |           0.740182 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00835     |       0.983467 |              0.37891  |           0.929032 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00453     |       0.98345  |              0.405666 |           0.74231  | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00733     |       0.983248 |              0.348827 |           0.914014 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00465     |       0.98297  |              0.365801 |           0.922342 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00279     |       0.982962 |              0.396363 |           0.73752  | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00745     |       0.982803 |              0.417064 |           0.94775  | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00430     |       0.982801 |              0.414327 |           0.946389 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00545     |       0.982216 |              0.402827 |           0.940473 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00372     |       0.982183 |              0.414877 |           0.746452 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00724     |       0.982082 |              0.401636 |           0.939835 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00305     |       0.982009 |              0.42173  |           0.7498   | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00821     |       0.982004 |              0.411961 |           0.744941 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00810     |       0.981974 |              0.404203 |           0.941074 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00877     |       0.981969 |              0.397603 |           0.73779  | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00397     |       0.981762 |              0.41661  |           0.747171 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00837     |       0.981724 |              0.361784 |           0.719889 | high             | True      | High leakage risk invalidates clean benchmark claim |

## Public-Strong / Internal-Weak Caveat Examples

| candidate_id   |   barneo_score |   confidence_score | abstention_reason_primary                           | abstention_reason_all                                                                                                                                                                                          |
|:---------------|---------------:|-------------------:|:----------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| CNV0_02422     |       0.918055 |           0.917663 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                           |
| CNV0_01433     |       0.938424 |           0.917061 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                           |
| CNV0_01058     |       0.890885 |           0.893315 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                           |
| CNV0_02444     |       0.908405 |           0.897742 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                           |
| CNV0_01063     |       0.85944  |           0.878358 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                           |
| CNV0_01122     |       0.865553 |           0.673044 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; HLA allele underrepresented in benchmark; High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors |
| CNV0_01124     |       0.828199 |           0.644953 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; HLA allele underrepresented in benchmark; High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors |
| CNV0_01101     |       0.806193 |           0.819462 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                           |
| CNV0_02442     |       0.60044  |           0.723227 | High disagreement across predictors                 | High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                                                                                |
| CNV0_02712     |       0.569068 |           0.715743 | High disagreement across predictors                 | High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                                                                                |
| CNV0_01135     |       0.524853 |           0.694767 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                           |
| CNV0_02555     |       0.248894 |           0.548151 | High disagreement across predictors                 | High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                                                                                |
| CNV0_02520     |       0.164053 |           0.503532 | High disagreement across predictors                 | High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                                                                                |
| CNV0_01349     |       0.175675 |           0.49415  | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                           |
| CNV0_02564     |       0.11501  |           0.469333 | High disagreement across predictors                 | High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                                                                                |

## Source-Shift Abstention Examples

| candidate_id   |   barneo_score |   confidence_score | abstention_reason_primary                           | abstention_reason_all                                                                                                                                                                                                          |
|:---------------|---------------:|-------------------:|:----------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| CNV0_01491     |       0.197423 |           0.519767 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01497     |       0.190755 |           0.511219 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01967     |       0.181233 |           0.509685 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01560     |       0.169124 |           0.500217 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01495     |       0.157922 |           0.496578 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01552     |       0.15837  |           0.496123 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01865     |       0.156117 |           0.496587 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01572     |       0.155166 |           0.494237 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01776     |       0.153598 |           0.489819 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01771     |       0.149989 |           0.492461 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01606     |       0.150751 |           0.491656 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; High disagreement across predictors; Low-prevalence source with poor historical top-k performance |
| CNV0_01937     |       0.140956 |           0.491318 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01922     |       0.145354 |           0.482618 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; High disagreement across predictors; Low-prevalence source with poor historical top-k performance |
| CNV0_01699     |       0.137246 |           0.48603  | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01785     |       0.134595 |           0.487971 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01808     |       0.138973 |           0.480384 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01563     |       0.129886 |           0.483635 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01493     |       0.131478 |           0.481641 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01685     |       0.13304  |           0.479433 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01822     |       0.130615 |           0.478662 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
