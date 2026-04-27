#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from v17p2_common import RPT, TAB


def j(name: str) -> dict:
    p = TAB / name
    return json.loads(p.read_text()) if p.exists() else {}


def t(name: str) -> pd.DataFrame:
    p = TAB / name
    return pd.read_csv(p, sep="\t") if p.exists() else pd.DataFrame()


def write(path: Path, text: str) -> None:
    path.write_text(text.strip() + "\n", encoding="utf-8")


def main() -> None:
    robust = j("robustness_summary.json")
    clinical = j("clinical_summary.json")
    biology = j("biology_summary.json")
    external = j("external_summary.json")
    dial = j("dial_method_summary.json")
    auc = t("external_classifier_comparison.tsv")
    gsea = t("gsea_dm1_vs_dm2.tsv")
    cox = t("cox_regression_results.tsv")

    top_auc = auc.sort_values("DIA_AUC", ascending=False).head(4) if not auc.empty else pd.DataFrame()
    top_gsea = gsea.sort_values("FDR q-val").head(8) if not gsea.empty else pd.DataFrame()
    cox_txt = cox.to_markdown(index=False) if not cox.empty else "_OS Cox unavailable_"

    outline = f"""
# v17p2 Paper Outline

## Title 후보
1. Orthogonal validation of DM1/DM2 thyroid cancer subtypes beyond the BRAF/RAS dichotomy
2. Reproducible Dark Matter taxonomy in papillary thyroid carcinoma across bulk and single-cell transcriptomes
3. Five-layer validation of driver-negative thyroid cancer subtypes with DIAL-guided external audit

## Abstract 골격
- Background: v14는 BRAF-like/RAS-like 축을 정리했지만 driver-negative 178명은 residual bucket으로 남아 있었다.
- Methods: Phase 2는 robustness, clinical, biology, external transfer, DIAL audit의 5-layer 검증을 수행했다.
- Results: bootstrap persistence는 {robust.get('median_persistence'):.3f}, 50% subsample ARI는 {robust.get('subsample_50_mean_ari'):.3f}, age difference p는 {clinical.get('age_p_welch'):.2e}, external robust cohort는 {external.get('cohorts_dia_auc_gt_0_85')}개였다.
- Results: pathway proxy에서 FDR<0.05 hallmark는 {biology.get('hallmark_fdr_lt_0_05')}개였고, scRNA {external.get('scrna_cells_scored')} cells에서 DM signature가 projection되었다.
- Interpretation: DM1/DM2는 purely unsupervised artifact가 아니라 age-stratified, biologically distinct, externally transferable subtype axis다.

## Figure 1-7
1. Dark Matter cohort and Phase 2 study design
2. Layer 1 robustness: multi-method ARI, bootstrap consensus, permutation null
3. Layer 2 clinical separation: age, stage, OS Cox
4. Layer 3 biology: pathway proxy, MAPK, thyroid differentiation
5. Layer 4 external validation: cohort transfer and scRNA projection
6. Layer 5 DIAL method audit: ComBat threshold and failure modes
7. Integrated model and reviewer-facing limitations

## Results 뼈대
1. DM1/DM2는 bootstrap과 subsampling에서도 유지된다.
2. DM2는 DM1보다 연령이 높고 differentiation state가 다르다.
3. 임상 endpoint power는 제한적이지만 age-independent transcriptional split은 유지된다.
4. Pathway proxy는 EMT/TGF-beta/inflammatory vs oxidative phosphorylation 축을 제시한다.
5. External transfer는 GSE27155와 GSE76039에서 강하고, scRNA projection은 cell-level heterogeneity를 보여준다.
6. DIAL audit는 ComBat over-correction을 정량화하며 signal destruction threshold를 제시한다.

## Discussion 포인트
- Strength: unsupervised cluster를 5개의 독립 axis로 검증했다.
- Weakness: robust external bulk cohort는 현재 2개, TERT/fusion raw recovery는 미완.
- Method angle: DIAL은 batch correction sensitivity audit로 보조 value가 있다.
- Next step: Korean cohort와 wet validation이 venue ceiling을 결정한다.

## Target venue
- 1순위: npj Precision Oncology
- 2순위: Genome Medicine
- 3순위: Bioinformatics
"""
    write(RPT / "v17p2_paper_outline.md", outline)

    response = f"""
# v17p2 Response To Reviewers

## Reviewer #1 expected concerns
1. Cluster artifact 주장: Layer 1 `bootstrap_consensus.html`, `robustness_method_consistency.html`로 대응.
2. Age confounding 주장: `age_adjusted_cluster_effect.tsv`에서 residual TDS p={t('age_adjusted_cluster_effect.tsv').iloc[-1]['pvalue'] if not t('age_adjusted_cluster_effect.tsv').empty else 'NA'}`.
3. Clinical irrelevance 주장: OS event는 적지만 age/stage/histology stratification을 제시.
4. Mechanism 부족 주장: `gsea_top_pathways_bar.html`, `mapk_activity_violin.html`, `thyroid_diff_score_violin.html`.
5. TCGA-only 주장: `external_4cohort_forest.html`, `scrna_dm_signature_umap.html`.

## Reviewer #2 expected concerns
1. ComBat artifact 주장: `dial_combat_threshold.html`에서 identifiability collapse를 정량화.
2. Alternative metric이면 충분하다는 주장: `dial_vs_alternatives_radar.html`.
3. Sample size 우려: `sample_subsample_stability.tsv`, bootstrap persistence.
4. scRNA projection이 bulk signature recycling이라는 주장: per-patient heterogeneity plot로 mixture를 제시.
5. Pan-cancer generalization 부족: 현재 supplementary limitation으로 명시.

## 예상 공격 10개와 답변
1. Random cluster? → bootstrap persistence {robust.get('median_persistence'):.3f}
2. Permutation null? → p={robust.get('permutation_p'):.4f}
3. Age only? → residual TDS after age p={t('age_adjusted_cluster_effect.tsv').iloc[-1]['pvalue'] if not t('age_adjusted_cluster_effect.tsv').empty else 'NA'}
4. Survival signal weak? → honest limitation, OS events={clinical.get('os_events')}
5. No mechanism? → hallmark proxy FDR<0.05 {biology.get('hallmark_fdr_lt_0_05')}
6. No external? → robust bulk cohorts={external.get('cohorts_dia_auc_gt_0_85')}
7. No single-cell? → scRNA cells scored={external.get('scrna_cells_scored')}
8. Batch effect? → lambda threshold={dial.get('combat_threshold_lambda_ident_lt_0_5')}
9. No pan-cancer method proof? → current limitation section으로 후퇴
10. No wet validation? → limitation에서 정직하게 인정
"""
    write(RPT / "v17p2_response_to_reviewers.md", response)

    dump = f"""
# v17p2 External Review Dump V2

## 1. Layer 1
- mean off-diagonal ARI: {robust.get('mean_ari_offdiag'):.3f}
- bootstrap persistence median: {robust.get('median_persistence'):.3f}
- permutation p-value: {robust.get('permutation_p'):.4f}
- 50% subsample ARI: {robust.get('subsample_50_mean_ari'):.3f}

## 2. Layer 2
- DM1/DM2 size: {clinical.get('dm1_n')} / {clinical.get('dm2_n')}
- age Welch p: {clinical.get('age_p_welch'):.2e}
- stage high Fisher p: {clinical.get('stage_high_p'):.3f}
- OS events: {clinical.get('os_events')}
- Cox table:
{cox_txt}

## 3. Layer 3
- hallmark/pathway proxy FDR<0.05 count: {biology.get('hallmark_fdr_lt_0_05')}
- MAPK activity DM2-DM1: {biology.get('mapk_dm2_minus_dm1'):.3f}
- top pathways:
{top_gsea.to_markdown(index=False) if not top_gsea.empty else '_NA_'}

## 4. Layer 4
- bulk external tested cohorts: {external.get('cohorts_tested')}
- robust cohorts (DIA-AUC>0.85): {external.get('cohorts_dia_auc_gt_0_85')}
- scRNA cells scored: {external.get('scrna_cells_scored')}
- best transfer:
{top_auc.to_markdown(index=False) if not top_auc.empty else '_NA_'}

## 5. Layer 5
- ComBat threshold (identifiability < 0.5): {dial.get('combat_threshold_lambda_ident_lt_0_5')}
- pancancer tested: {dial.get('pancancer_tested')}
- max pan-cancer DIA-AUC: {dial.get('max_pancancer_dia_auc')}

## 6. Self-assessment
| Item | Score | Note |
|---|---:|---|
| Statistical robustness | 4/5 | persistence strong, ARI moderate |
| Clinical stratification | 2/5 | age strong, survival weak |
| Biological interpretability | 4/5 | pathway proxy clear |
| External validation | 3/5 | 2 bulk cohorts + scRNA |
| Method audit | 3/5 | THCA audit strong, pan-cancer partial |
| Wet validation | 0/5 | 없음 |

## 7. Venue
- npj Precision Oncology: 55%
- Genome Medicine: 25%
- Bioinformatics: 75%
- Nature Communications: 8%

## 8. 다음 결정 3개
1. GSE213647 gene harmonization을 더 깊게 밀어 bulk robust cohort를 3개 이상으로 회복할지.
2. TERT/raw fusion recovery를 별도 1주 sprint로 분리할지.
3. 지금 바로 manuscript draft로 들어갈지, 아니면 OS/clinical metadata를 더 보강할지.
"""
    write(RPT / "v17p2_external_review_dump_v2.md", dump)


if __name__ == "__main__":
    main()
