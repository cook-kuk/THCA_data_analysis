
# BAR-Neo Abstention Report

## Abstention Counts

- Candidates scored: 2715
- Abstain true: 2702
- Abstain false: 13

## Top Abstention Reasons

| abstention_reason                                                          |   count |
|:---------------------------------------------------------------------------|--------:|
| High leakage risk invalidates clean benchmark claim                        |    2304 |
| High disagreement across predictors                                        |    1114 |
| Low-prevalence source with poor historical top-k performance               |     915 |
| High source-shift risk: similar heldout source had poor top-k survival     |     915 |
| HLA allele underrepresented in benchmark                                   |     274 |
| Low confidence due to sparse method coverage                               |     183 |
| High disagreement between internal anchor and public pretrained predictors |      15 |
| no_abstention_reason                                                       |      13 |

## High Score / Low Confidence Cases

| candidate_id   |   barneo_score |   patient_gated_score |   confidence_score | confidence_bin   | abstain   | abstention_reason_primary                           |
|:---------------|---------------:|----------------------:|-------------------:|:-----------------|:----------|:----------------------------------------------------|
| CNV0_00252     |       0.984759 |              0.409963 |           0.744893 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00305     |       0.983925 |              0.420643 |           0.749907 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00562     |       0.983826 |              0.419779 |           0.949444 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00460     |       0.983689 |              0.402047 |           0.740596 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00372     |       0.983412 |              0.413343 |           0.746109 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00293     |       0.983387 |              0.385407 |           0.932229 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00561     |       0.983327 |              0.413649 |           0.746231 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00874     |       0.982798 |              0.351429 |           0.915139 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00078     |       0.982555 |              0.412257 |           0.945276 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00178     |       0.982522 |              0.392043 |           0.93522  | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00553     |       0.982102 |              0.416759 |           0.94736  | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00848     |       0.981886 |              0.36444  |           0.921269 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00471     |       0.981737 |              0.407277 |           0.742521 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00453     |       0.981612 |              0.402675 |           0.740188 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00761     |       0.981581 |              0.416773 |           0.947191 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00545     |       0.981536 |              0.402314 |           0.939982 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00898     |       0.981317 |              0.395843 |           0.936687 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00211     |       0.981291 |              0.406142 |           0.741802 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00166     |       0.981277 |              0.416377 |           0.746891 | high             | True      | High leakage risk invalidates clean benchmark claim |
| CNV0_00473     |       0.981091 |              0.404439 |           0.740886 | high             | True      | High leakage risk invalidates clean benchmark claim |

## Public-Strong / Internal-Weak Caveat Examples

| candidate_id   |   barneo_score |   confidence_score | abstention_reason_primary                           | abstention_reason_all                                                                                                                                                                                          |
|:---------------|---------------:|-------------------:|:----------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| CNV0_02422     |       0.926148 |           0.922114 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                           |
| CNV0_02444     |       0.927551 |           0.908273 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                           |
| CNV0_01433     |       0.90912  |           0.901055 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                           |
| CNV0_01058     |       0.895602 |           0.89569  | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                           |
| CNV0_01063     |       0.827538 |           0.860996 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                           |
| CNV0_01122     |       0.839553 |           0.658713 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; HLA allele underrepresented in benchmark; High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors |
| CNV0_01124     |       0.829336 |           0.646767 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; HLA allele underrepresented in benchmark; High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors |
| CNV0_01101     |       0.836461 |           0.835475 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                           |
| CNV0_01135     |       0.574978 |           0.722161 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                           |
| CNV0_02442     |       0.610714 |           0.728175 | High disagreement across predictors                 | High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                                                                                |
| CNV0_02712     |       0.471438 |           0.661245 | High disagreement across predictors                 | High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                                                                                |
| CNV0_02555     |       0.242358 |           0.543857 | High disagreement across predictors                 | High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                                                                                |
| CNV0_02520     |       0.151862 |           0.496636 | High disagreement across predictors                 | High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                                                                                |
| CNV0_01349     |       0.163477 |           0.487717 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                           |
| CNV0_02564     |       0.116013 |           0.4699   | High disagreement across predictors                 | High disagreement across predictors; High disagreement between internal anchor and public pretrained predictors                                                                                                |

## Source-Shift Abstention Examples

| candidate_id   |   barneo_score |   confidence_score | abstention_reason_primary                           | abstention_reason_all                                                                                                                                                                                                          |
|:---------------|---------------:|-------------------:|:----------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| CNV0_01830     |       0.181295 |           0.507213 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; High disagreement across predictors; Low-prevalence source with poor historical top-k performance |
| CNV0_01657     |       0.185382 |           0.502406 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; High disagreement across predictors; Low-prevalence source with poor historical top-k performance |
| CNV0_01606     |       0.164067 |           0.500693 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01491     |       0.159879 |           0.501243 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01495     |       0.154681 |           0.495482 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01497     |       0.151404 |           0.491418 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01967     |       0.14614  |           0.492237 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01616     |       0.163265 |           0.479133 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; High disagreement across predictors; Low-prevalence source with poor historical top-k performance |
| CNV0_01521     |       0.149735 |           0.482848 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; High disagreement across predictors; Low-prevalence source with poor historical top-k performance |
| CNV0_01552     |       0.142008 |           0.489166 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01737     |       0.145049 |           0.486324 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01731     |       0.164157 |           0.469378 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; High disagreement across predictors; Low-prevalence source with poor historical top-k performance |
| CNV0_01808     |       0.142199 |           0.483696 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01699     |       0.137265 |           0.485813 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01560     |       0.135763 |           0.484976 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01937     |       0.133831 |           0.486683 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01943     |       0.130673 |           0.478397 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01820     |       0.124842 |           0.482654 | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01865     |       0.126671 |           0.47997  | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
| CNV0_01785     |       0.121764 |           0.48406  | High leakage risk invalidates clean benchmark claim | High leakage risk invalidates clean benchmark claim; High source-shift risk: similar heldout source had poor top-k survival; Low-prevalence source with poor historical top-k performance                                      |
