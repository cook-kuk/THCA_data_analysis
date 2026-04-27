# v17p2 Paper Outline

## Title 후보
1. Orthogonal validation of DM1/DM2 thyroid cancer subtypes beyond the BRAF/RAS dichotomy
2. Reproducible Dark Matter taxonomy in papillary thyroid carcinoma across bulk and single-cell transcriptomes
3. Five-layer validation of driver-negative thyroid cancer subtypes with DIAL-guided external audit

## Abstract 골격
- Background: v14는 BRAF-like/RAS-like 축을 정리했지만 driver-negative 178명은 residual bucket으로 남아 있었다.
- Methods: Phase 2는 robustness, clinical, biology, external transfer, DIAL audit의 5-layer 검증을 수행했다.
- Results: bootstrap persistence는 0.903, 50% subsample ARI는 0.757, age difference p는 6.19e-07, external robust cohort는 2개였다.
- Results: pathway proxy에서 FDR<0.05 hallmark는 10개였고, scRNA 66015 cells에서 DM signature가 projection되었다.
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
