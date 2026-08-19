# K-Thyro Extra Public Analyses Round 2

생성 시점: 2026-05-09  
위치: `/home/seungho/personal/THCA_data_analysis/kthyro_public_pilot/results/extra_analyses`

목적: 첫 추가 분석 이후, 심사위원 반박 가능성이 높은 세 질문을 더 정량화.

1. Spatial interface/hotspot: 서로 다른 치료취약성 state가 실제로 가까이/멀리 배치되는가?
2. TCGA multivariable/confounding: driver, stage, epithelial score로 축이 다 설명되는가?
3. External meta-consistency: 외부 thyroid cohort에서 방향성이 cohort 간 일관적인가?

---

## 1. Spatial interface and hotspot analysis

### 질문

Spatial niche가 응집한다는 것에서 한 단계 더 나아가, 특정 치료취약성 state끼리 **공간적으로 가까운 interface** 또는 **co-localized hotspot**을 만드는가?

### 방법

- Input: `spatial_spot_vulnerability_scores.tsv`
- Dataset: GSE250521, 16 slides
- 각 slide에서 top/bottom quartile mask 생성:
  - tumor_high
  - CD8 high/low
  - HLA-low tumor
  - RAI-low tumor
  - barrier high
  - delivery failure high
  - aggressive high
  - vascular high
  - hypoxia high
- Interface distance:
  - source mask에서 target mask까지 nearest distance 계산
  - target label randomization 75 permutations
  - negative z = target이 random보다 더 가까움
- Hotspot co-localization:
  - 두 mask의 overlap enrichment 계산
  - hypergeometric p는 spot-level proxy로만 사용

### 주요 결과: interface distance

Median distance z by metric:

| Metric | Median distance delta | Median z |
|---|---:|---:|
| barrier_to_cd8_high | -1.19 | -2.52 |
| rai_low_tumor_to_aggressive | -14.09 | -0.93 |
| rai_low_tumor_to_delivery_failure | 9.12 | 0.83 |
| delivery_failure_to_cd8_low | 1.87 | 1.53 |
| tumor_to_cd8_high | 0.00 | 1.63 |
| hla_low_tumor_to_cd8_high | 30.43 | 6.98 |
| hla_low_tumor_to_barrier_high | 90.77 | 10.23 |
| delivery_failure_to_vascular_high | 41.83 | 13.96 |

### 해석

쓸 수 있는 메시지:

- Barrier-high territory와 CD8-high territory는 여러 slide에서 random보다 가까운 interface를 보이는 경향이 있다. 이는 CD8 exclusion을 단순 absence가 아니라 **interface biology**로 검증할 수 있음을 시사한다.
- RAI-low tumor와 aggressive-high territory도 일부 slide에서 가까운 경향을 보인다.
- HLA-low tumor가 CD8-high 또는 barrier-high와 항상 붙는 것은 아니며, HLA-low는 별도의 immune-invisible territory일 수 있다.
- Delivery-failure high가 vascular-high와 가깝지 않고 오히려 멀게 나오는 경향은 “delivery failure proxy” 가설과 방향성이 맞다.

주의:

- nearest-distance 분석은 spot coordinate proxy다.
- negative/positive distance z를 mechanistic exclusion proof로 해석하지 않는다.
- 실제 검증은 IF/mIHC에서 tumor-HLA, CD8, CAF/myeloid marker의 interface distance로 해야 한다.

### 주요 결과: hotspot co-localization

Median hotspot enrichment:

| Hotspot | Median log2 overlap enrichment | Median Jaccard |
|---|---:|---:|
| delivery-failure + hypoxia-high | 1.094 | 0.364 |
| HLA-high + IFN-high | 0.489 | 0.218 |
| RAI-low tumor + aggressive-high | 0.455 | 0.024 |
| HLA-high + CD8-high | 0.400 | 0.250 |
| HLA-low tumor + CD8-low | 0.225 | 0.043 |
| RAI-low tumor + delivery-failure | 0.105 | 0.018 |
| barrier-high + CD8-low | -0.112 | 0.205 |
| HLA-low tumor + barrier-high | -0.874 | 0.017 |

### 해석

가장 좋은 문장:

> Spatial hotspot analysis showed co-localization of delivery-failure proxy with hypoxia-high territories and of HLA-high with IFN/CD8-high territories, supporting biologically interpretable spatial vulnerability states.

더 조심할 문장:

> HLA-low tumor and barrier-high did not robustly co-localize in this public dataset, suggesting that immune invisibility and stromal barrier may be separable spatial programs.

### Outputs

- `tables/spatial_interface_distance_tests.tsv`
- `tables/spatial_hotspot_colocalization_tests.tsv`
- `figures/spatial_interface_distance_z_by_condition.png`
- `figures/spatial_hotspot_colocalization_heatmap.png`

---

## 2. TCGA multivariable / confounding checks

### 질문

K-Thyro axes가 driver mutation, stage, epithelial/tumor-lineage score, age, sex로 다 설명되는가?

### 방법

각 axis에 대해 descriptive linear model R² 계산:

1. driver_only
2. driver_stage
3. driver_stage_epithelial
4. driver_stage_epithelial_age_sex

또한 가장 큰 모델의 residual이 vulnerability label 사이에서 여전히 다른지 Kruskal test.

### R² 결과

| Axis | Driver only | Driver+stage | Driver+stage+epithelial | Full model |
|---|---:|---:|---:|---:|
| RAI differentiation | 0.309 | 0.329 | 0.448 | 0.459 |
| Drug-delivery failure | 0.332 | 0.359 | 0.401 | 0.405 |
| Aggressive dediff. | 0.273 | 0.280 | 0.342 | 0.396 |
| Myeloid/CAF barrier | 0.182 | 0.190 | 0.239 | 0.249 |
| CD8 exclusion | 0.162 | 0.179 | 0.224 | 0.224 |
| HLA-I/APM | 0.156 | 0.162 | 0.182 | 0.203 |
| Immune visibility | 0.125 | 0.126 | 0.147 | 0.169 |
| Proliferation | 0.107 | 0.108 | 0.146 | 0.171 |

### 핵심 해석

Full model로도 상당 부분이 설명되지 않음:

- RAI differentiation unexplained fraction 약 54.1%
- Drug-delivery failure unexplained fraction 약 59.5%
- HLA-I/APM unexplained fraction 약 79.7%
- Immune visibility unexplained fraction 약 83.1%
- CD8 exclusion unexplained fraction 약 77.6%

### Residual vulnerability label separation

Driver/stage/epithelial/age/sex를 보정한 residual도 vulnerability label 사이에서 강하게 다름:

| Axis | Residual label Kruskal p | FDR q |
|---|---:|---:|
| Drug-delivery failure | 6.36e-38 | 5.09e-37 |
| Myeloid/CAF barrier | 3.98e-30 | 1.59e-29 |
| RAI differentiation | 8.63e-30 | 2.30e-29 |
| CD8 exclusion | 7.05e-29 | 1.41e-28 |
| HLA-I/APM | 7.17e-26 | 1.15e-25 |
| Immune visibility | 9.68e-26 | 1.29e-25 |
| Aggressive dediff. | 1.76e-24 | 2.01e-24 |
| Proliferation | 2.67e-08 | 2.67e-08 |

### 제안서용 문장

> TCGA multivariable checks showed that driver group, stage, epithelial-lineage score, age, and sex explain only part of K-Thyro therapeutic axes. Even after these covariates are regressed out, vulnerability labels remain strongly separated, supporting the need for a dedicated therapeutic vulnerability atlas rather than driver-only stratification.

### Claim boundary

- Descriptive regression, not causal inference.
- Missing variables include fusion, CNV, methylation, tumor purity, treatment history.
- Use as robustness check only.

### Outputs

- `tables/tcga_axis_variance_models.tsv`
- `tables/tcga_residual_vulnerability_label_tests.tsv`
- `tables/tcga_axis_covariate_correlations.tsv`
- `figures/tcga_axis_variance_models.png`
- `figures/tcga_axis_covariate_correlation_heatmap.png`

---

## 3. External exact K-Thyro meta-consistency

### 질문

External GPL570 cohorts에서 K-Thyro axis 방향성이 cohort마다 재현되는가?

### 방법

Input:

- `exact_kthyro_external_bulk_contrasts.tsv`

대상:

- ATC vs normal
- ATC vs PTC
- PTC vs normal
- PDTC vs PTC where available

예상 방향:

- ATC/advanced: RAI differentiation down
- ATC/advanced: RAI-restorable down
- ATC/advanced: aggressive/proliferation/hypoxia/drug-delivery-failure/myeloid-CAF up

Signed Stouffer p-value and direction-match fraction 계산.

### 결과: consistent strong

| Contrast | Axis | Direction match | Median effect | FDR q |
|---|---|---:|---:|---:|
| ATC vs normal | RAI differentiation down | 3/3 | -3.10 | 1.49e-12 |
| ATC vs normal | Drug-delivery failure up | 3/3 | 2.61 | 1.59e-12 |
| ATC vs normal | RAI-restorable down | 3/3 | -2.39 | 2.97e-12 |
| ATC vs normal | Aggressive dediff. up | 3/3 | 2.43 | 3.74e-12 |
| ATC vs normal | Myeloid/CAF barrier up | 3/3 | 4.59 | 3.74e-12 |
| PTC vs normal | RAI differentiation down | 3/3 | -1.92 | 5.40e-11 |
| ATC vs normal | Proliferation up | 3/3 | 2.48 | 5.40e-11 |
| ATC vs normal | Hypoxia up | 3/3 | 1.83 | 3.23e-10 |
| ATC vs PTC | Drug-delivery failure up | 2/2 | 2.36 | 4.96e-08 |
| ATC vs PTC | RAI differentiation down | 2/2 | -2.95 | 9.69e-08 |
| ATC vs PTC | Myeloid/CAF barrier up | 2/2 | 1.85 | 7.82e-07 |
| ATC vs PTC | Proliferation up | 2/2 | 2.13 | 1.74e-06 |
| ATC vs PTC | Aggressive dediff. up | 2/2 | 1.80 | 2.57e-05 |
| ATC vs PTC | RAI-restorable down | 2/2 | -1.90 | 2.71e-05 |
| ATC vs PTC | Hypoxia up | 2/2 | 1.49 | 4.38e-05 |

Weak/inconsistent:

- PDTC vs PTC has only one dataset and is underpowered/inconsistent for aggressive dediff.

### 제안서용 문장

> Exact K-Thyro rescoring of independent external thyroid cohorts showed highly consistent ATC-associated reduction of RAI differentiation and increase of aggressive, myeloid/CAF, hypoxia, and drug-delivery-failure proxy axes, supporting cross-cohort reproducibility of the therapeutic vulnerability framework.

### Claim boundary

- External cohorts are bulk microarray and histology contrasts.
- This validates expression-axis directionality, not therapeutic response or spatial niche biology.

### Outputs

- `tables/external_exact_kthyro_meta_consistency.tsv`
- `figures/external_exact_kthyro_meta_direction_consistency.png`

---

## 4. Round 2 strongest takeaways

1. **External reproducibility is now very strong.**  
   ATC vs normal/PTC contrasts reproduce RAI-low, aggressive-high, drug-delivery-failure-high, hypoxia-high, myeloid/CAF-high directions across independent cohorts.

2. **K-Thyro axes are not just driver/stage/epithelial score.**  
   Full covariate models explain only part of axis variance, and residuals still separate vulnerability labels.

3. **Spatial hotspot biology became more concrete.**  
   Delivery-failure + hypoxia-high and HLA-high + IFN/CD8-high co-localization patterns are biologically interpretable.

4. **Some axes are separable, not always co-localized.**  
   HLA-low tumor + barrier-high was not robustly co-localized, which supports the idea that immune invisibility and stromal barrier are separable programs.

---

## 5. Round 2 claim boundaries

Do not overclaim:

- Spatial interface distance is not proof of functional exclusion.
- Hypergeometric hotspot p-values are spot-level proxies under spatial autocorrelation.
- TCGA regression is descriptive, not causal.
- External meta-consistency validates axis direction, not clinical response.

Strong but honest wording:

> Public data consistently support separable, spatially organized, and externally reproducible thyroid cancer therapeutic vulnerability axes. These axes require FFPE/mIHC, GeoMx/ROI, fresh tissue perturbation, iodide uptake, HLA/APM rescue, and drug/nanoparticle penetration imaging for functional validation.

