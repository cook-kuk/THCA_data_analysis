# K-Thyro Extra Public Analyses Summary

생성 시점: 2026-05-09  
위치: `/home/seungho/personal/THCA_data_analysis/kthyro_public_pilot/results/extra_analyses`

목적: 기존 public-data pilot 이후, 삼성 리뷰어가 물을 가능성이 높은 추가 질문을 공개데이터로 더 검증.

추가 분석 6개:

1. Spatial niche adjacency / exclusion analysis
2. Spatial condition trend analysis
3. TCGA BRAF-only therapeutic heterogeneity analysis
4. TCGA clinical-risk exploratory association
5. GSE151179 direct RAI label reanalysis
6. scRNA module cell-type attribution

---

## 1. 핵심 결론

추가 분석 후에도 **GO 유지**입니다.

기존 강점:

- TCGA 505명에서 vulnerability label 분리.
- GSE250521 spatial 16 slide에서 16/16 non-random niche coherence.
- Public validation matrix v2 301 tests.
- External exact K-Thyro validation 206 samples.

추가로 강해진 점:

- Spatial niche는 단순히 응집될 뿐 아니라 **특정 niche끼리 이웃하는 adjacency pattern**을 보임.
- BRAF-mutant TCGA tumors 내부에서도 vulnerability axes가 매우 강하게 갈라짐.
- Stage III/IV TCGA tumors에서 **RAI-low, drug-delivery-failure-high, CD8-exclusion-high** exploratory association이 관찰됨.
- scRNA attribution은 module들이 예상 cell type과 대체로 맞는 방향임을 지지함.

약한/방어해야 할 점:

- GSE151179 direct RAI avid/refractory label은 강한 validation이 아님.
- RAI response는 public label로 overclaim하지 말고 prospective iodide uptake assay 필요성으로 전환해야 함.
- Spatial adjacency 통계는 spot-level proxy이며 biological replicate p-value로 주장하면 안 됨.

---

## 2. Spatial niche adjacency analysis

### 질문

기존 분석은 “같은 niche가 공간적으로 응집되는가?”를 보았다.  
추가 분석은 “서로 다른 therapeutic niche가 어떤 공간 이웃 관계를 만드는가?”를 물었다.

### 방법

- Input: `spatial_spot_vulnerability_scores.tsv`
- Dataset: GSE250521
- Slides: 16
- Spots: 57,144
- 각 slide에서 spot coordinate 기반 KNN, k=6
- Source niche -> neighbor niche adjacency 계산
- Expected fraction은 slide-level label frequency 기반
- Output은 log2 enrichment proxy와 binomial z proxy

### 결과

가장 강한 median adjacency enrichment:

| Source niche | Neighbor niche | Median log2 enrichment proxy |
|---|---|---:|
| RAI-restorable niche | RAI-restorable niche | 1.350 |
| APC-rich niche | APC-rich niche | 1.286 |
| HLA-low invisible tumor niche | HLA-low invisible tumor niche | 1.194 |
| CD8-excluded / myeloid-CAF niche | CD8-excluded / myeloid-CAF niche | 0.564 |
| HLA-visible inflamed niche | HLA-visible inflamed niche | 0.540 |
| Drug-delivery failure proxy niche | Drug-delivery failure proxy niche | 0.432 |
| RAI-low dedifferentiated niche | Drug-delivery failure proxy niche | 0.244 |
| Drug-delivery failure proxy niche | RAI-low dedifferentiated niche | 0.244 |

### 해석

Observed data:

- 각 niche는 자기 자신과 가장 강하게 이웃하는 spatial coherence를 보임.
- RAI-low dedifferentiated niche와 drug-delivery failure proxy niche가 서로 이웃하는 패턴이 관찰됨.

Exploratory hypothesis:

- dedifferentiated/RAI-low territory가 barrier 또는 delivery-failure territory와 공간적으로 맞물릴 수 있음.
- 이 패턴은 fresh tissue slice에서 iodide uptake와 drug penetration imaging을 같이 봐야 검증 가능.

Claim boundary:

- 이 분석은 spot-neighborhood proxy이다.
- Spots are not independent biological replicates.
- p-value 기반 biological claim이 아니라 ROI nomination evidence로 사용해야 한다.

### Outputs

- `tables/spatial_niche_adjacency_enrichment.tsv`
- `tables/spatial_functional_adjacency_metrics.tsv`
- `figures/spatial_niche_adjacency_enrichment_heatmap.png`
- `figures/spatial_condition_and_adjacency_trends.png`

---

## 3. Spatial condition trend analysis

### 질문

Normal -> PTC -> locally advanced PTC -> ATC 순서로 slide-level niche fraction이나 adjacency metric이 변하는가?

### 방법

- Slide-level summary만 사용.
- n=4 slides per condition.
- Spearman trend across condition order and Kruskal test.
- Exploratory only.

### 주요 결과

| Metric | Spearman rho | p | FDR q | Normal mean | PTC mean | LPTC mean | ATC mean |
|---|---:|---:|---:|---:|---:|---:|---:|
| fraction_CD8-excluded / myeloid-CAF niche | -0.776 | 0.000408 | 0.0118 | 0.152 | 0.0806 | 0.0784 | 0.0422 |
| delivery_high_neighbor_cd8_low_fraction | -0.715 | 0.00183 | 0.0266 | 0.829 | 0.712 | 0.727 | 0.364 |
| hla_low_tumor_neighbor_cd8_high_fraction | -0.616 | 0.0110 | 0.0926 | 0.802 | 0.392 | 0.583 | 0.137 |
| same_niche_z | 0.606 | 0.0128 | 0.0926 | 11.49 | 25.46 | 30.16 | 26.29 |
| fraction_HLA-low invisible tumor niche | 0.509 | 0.0439 | 0.141 | 0.0083 | 0.0263 | 0.0167 | 0.0362 |
| fraction_RAI-low dedifferentiated niche | 0.292 | 0.273 | 0.439 | 0.0012 | 0.0013 | 0.0006 | 0.0234 |

### 해석

쓸 수 있는 메시지:

- Spatial niche coherence는 normal보다 cancer/advanced slides에서 더 커지는 경향이 있다.
- ATC에서 HLA-low invisible tumor niche와 RAI-low dedifferentiated niche fraction이 증가하는 방향성이 보인다.

주의할 메시지:

- CD8-excluded/myeloid-CAF niche fraction은 condition order에서 감소 방향으로 나왔다. 이건 단순 stage progression claim으로 쓰면 안 된다.
- 조건별 slide 수가 4개라 group comparison은 exploratory이다.
- Slide annotation, tissue composition, ROI/tissue capture 면적 차이가 영향을 줄 수 있다.

### Best use

PPT에서는 condition trend를 메인 claim으로 앞세우지 말고, **spatial heterogeneity and ROI nomination** 보조근거로 사용.

---

## 4. TCGA BRAF-only therapeutic heterogeneity

### 질문

BRAF mutation만 알면 치료취약성 상태가 결정되는가?

### 결과

BRAF-mutant 274명 내부 label 분포:

| Vulnerability label | n | Fraction of BRAF |
|---|---:|---:|
| Drug-delivery barrier-high | 91 | 0.332 |
| Mixed/Other | 85 | 0.310 |
| CD8-excluded myeloid/CAF-high | 34 | 0.124 |
| HLA-visible inflamed | 19 | 0.069 |
| HLA-low immune-invisible | 18 | 0.066 |
| RAI-low dedifferentiated | 18 | 0.066 |
| RAI-readable differentiated | 9 | 0.033 |

모든 K-Thyro axis가 BRAF label 사이에서 유의하게 다름:

| Axis | Kruskal p | FDR q |
|---|---:|---:|
| Drug-delivery failure proxy | 1.48e-38 | 1.18e-37 |
| CD8 exclusion | 1.22e-24 | 3.54e-24 |
| Myeloid/CAF barrier | 1.33e-24 | 3.54e-24 |
| RAI differentiation | 2.92e-17 | 5.84e-17 |
| Immune visibility | 1.72e-15 | 2.75e-15 |
| HLA-I/APM | 2.38e-15 | 3.18e-15 |
| Aggressive dedifferentiation | 2.85e-15 | 3.26e-15 |
| Proliferation | 6.74e-07 | 6.74e-07 |

### 해석

이 결과는 매우 강하게 쓸 수 있음:

> BRAF mutation은 중요한 driver이지만, BRAF-mutant 갑상선암 내부에서도 RAI, HLA/APM, CD8 exclusion, myeloid/CAF barrier, drug-delivery failure proxy가 서로 다른 치료취약성 상태로 갈라진다.

Claim boundary:

- TCGA mutation/expression association.
- 치료반응 예측이나 BRAF-targeted therapy response를 증명하지 않음.

### Outputs

- `tables/tcga_braf_only_label_distribution.tsv`
- `tables/tcga_braf_only_axis_heterogeneity_tests.tsv`
- `figures/tcga_braf_only_vulnerability_axis_heatmap.png`

---

## 5. TCGA clinical-risk exploratory association

### 질문

K-Thyro axes가 TCGA clinical-risk feature와 방향성 있게 연관되는가?

### Stage III/IV vs I/II 결과

| Axis | n high stage | n low stage | Cohen's d | MW p | FDR q |
|---|---:|---:|---:|---:|---:|
| Drug-delivery failure proxy | 166 | 336 | 0.471 | 3.0e-06 | 2.4e-05 |
| RAI differentiation | 166 | 336 | -0.419 | 7.0e-06 | 2.7e-05 |
| CD8 exclusion | 166 | 336 | 0.381 | 1.98e-04 | 5.28e-04 |
| Aggressive dedifferentiation | 166 | 336 | 0.291 | 0.00778 | 0.0156 |
| Myeloid/CAF barrier | 166 | 336 | 0.241 | 0.0153 | 0.0245 |
| HLA-I/APM | 166 | 336 | 0.164 | 0.0252 | 0.0337 |

### Survival exploratory

PFI:

- Proliferation: Cox HR per SD 1.51, p=0.00273, FDR q=0.0436.
- RAI differentiation: median split logrank p=0.00708 but Cox FDR not significant.
- Aggressive dedifferentiation: logrank p=0.0339 but Cox FDR not significant.
- Drug-delivery failure proxy: Cox p=0.0936, not significant after FDR.

### 해석

쓸 수 있는 메시지:

- Higher stage tumors show lower RAI differentiation and higher drug-delivery failure/CD8-exclusion/aggressive axes.
- Proliferation has exploratory PFI association.

주의:

- This is TCGA exploratory association, not a validated prognostic model.
- Survival event counts are limited in THCA.
- Do not claim treatment response prediction.

### Outputs

- `tables/tcga_clinical_axis_associations_expanded.tsv`
- `tables/tcga_survival_exploratory_axis_associations.tsv`
- `figures/tcga_clinical_stage_axis_associations.png`

---

## 6. GSE151179 direct RAI label reanalysis

### 질문

Public RAI avid/refractory label이 K-Thyro RAI score를 직접 검증해 주는가?

### 결과

강한 검증은 아님.

Best directional signals:

| Context | Comparison | Score | n A | n B | Cohen's d | p | FDR q | Support |
|---|---|---|---:|---:|---:|---:|---:|---|
| all_tumor | no_uptake_vs_uptake | TPO | 20 | 19 | -0.589 | 0.119 | 0.989 | weak_expected |
| primary_pre_rai_tumor | no_uptake_vs_uptake | TG | 6 | 11 | -0.577 | 0.591 | 0.989 | weak_expected |
| metastatic_or_nonprimary_tumor | no_uptake_vs_uptake | DIO1 | 14 | 8 | -0.569 | 0.402 | 0.989 | weak_expected |
| primary_pre_rai_tumor | no_uptake_vs_uptake | TPO | 6 | 11 | -0.533 | 0.301 | 0.989 | weak_expected |

Before/after paired:

- Only **2 paired patients**, not enough for inference.
- RAI8, TG, TPO, TSHR, PAX8 all after-minus-before negative in these 2 patients, but p-value not testable.

### 해석

이 결과는 앞세우면 안 됨.

가장 좋은 사용법:

> Public RAI response labels are noisy and underpowered; therefore K-Thyro must prospectively validate RAI-restorable states using iodide uptake assays and perturbation experiments.

Claim boundary:

- RAI-restorable score is hypothesis only.
- True validation requires iodide uptake assay after MAPK/RET/NTRK/BRAF/MEK-axis perturbation.

### Outputs

- `tables/gse151179_rai_direct_label_reanalysis.tsv`
- `tables/gse151179_rai_before_after_paired_reanalysis.tsv`
- `figures/gse151179_rai_direct_label_reanalysis.png`

---

## 7. scRNA module cell-type attribution

### 질문

Spatial module score가 어떤 cell-type signal과 관련되는가?

### 결과

| Module | Top cell type | Top mean | Second cell type | Delta |
|---|---|---:|---|---:|
| RAI differentiation | Epithelial cell | 3.987 | Malignant cell | 3.269 |
| HLA-II/APC | Myeloid cell | 1.439 | B cell | 0.653 |
| Cytotoxic T | NK cell | 1.625 | T cell | 1.322 |
| Myeloid/TAM | Myeloid cell | 1.296 | Fibroblast | 1.375 |
| CAF/ECM/TGFB | Fibroblast | 2.089 | Malignant cell | 1.412 |
| Vascular delivery proxy | Endothelial cell | 4.087 | Fibroblast | 3.832 |
| Tumor epithelial | Epithelial cell | 2.922 | Malignant cell | 1.946 |

### 해석

쓸 수 있는 메시지:

- RAI/tumor-epithelial module은 epithelial/malignant lineage와 맞음.
- HLA-II/APC module은 myeloid/B cell contribution 가능성이 큼.
- CAF/ECM module은 fibroblast에서 높고 vascular proxy는 endothelial에서 높음.
- Spatial score는 cell composition + activation state mixture로 해석해야 함.

주의:

- scRNA는 Visium deconvolution을 직접 수행한 것이 아니다.
- HLA-I/APM coverage가 약했던 processed scRNA object이므로 HLA-I/APM validation으로 쓰면 안 됨.

### Outputs

- `tables/scrna_module_celltype_attribution.tsv`
- `figures/scrna_module_celltype_attribution_heatmap.png`

---

## 8. 삼성 제안서에 추가할 수 있는 새 문장

### 강한 문장

> 추가 spatial adjacency 분석에서 RAI-restorable, APC-rich, HLA-low invisible, CD8-excluded/myeloid-CAF, drug-delivery failure niche는 각각 자기 자신과 강하게 이웃하는 spatial self-enrichment를 보였고, RAI-low dedifferentiated niche와 drug-delivery failure proxy niche 사이의 상호 adjacency도 관찰되었다. 이는 치료취약성 niche가 단순 spot-level noise가 아니라 조직 내 치료생태계 구조로 배열될 가능성을 지지한다.

### Mutation-only 방어 문장

> BRAF-mutant TCGA tumors 274명 내부에서도 7개 vulnerability label이 관찰되었고, drug-delivery failure, CD8 exclusion, myeloid/CAF barrier, RAI differentiation, HLA/APM 축이 BRAF 내부 label 사이에서 모두 유의하게 달랐다. 따라서 driver mutation은 중요하지만, 치료취약성 공간상태를 대체하지 못한다.

### Clinical cautious 문장

> TCGA exploratory clinical analysis에서 stage III/IV tumors는 stage I/II tumors 대비 drug-delivery failure proxy, CD8 exclusion, aggressive dedifferentiation이 높고 RAI differentiation이 낮은 방향성을 보였다. 이는 clinical deployment 근거가 아니라, high-risk tissue에서 공간 치료취약성 validation을 우선 수행해야 할 근거이다.

### RAI boundary 문장

> Public RAI avid/refractory labels alone did not provide strong direct validation of RAI-restorable scores, underscoring the need for prospective iodide uptake assays after MAPK-axis perturbation.

---

## 9. 추가 분석 figure/table paths

### Tables

- `tables/spatial_niche_adjacency_enrichment.tsv`
- `tables/spatial_functional_adjacency_metrics.tsv`
- `tables/spatial_condition_trend_summary.tsv`
- `tables/tcga_braf_only_label_distribution.tsv`
- `tables/tcga_braf_only_axis_heterogeneity_tests.tsv`
- `tables/tcga_clinical_axis_associations_expanded.tsv`
- `tables/tcga_survival_exploratory_axis_associations.tsv`
- `tables/gse151179_rai_direct_label_reanalysis.tsv`
- `tables/gse151179_rai_before_after_paired_reanalysis.tsv`
- `tables/scrna_module_celltype_attribution.tsv`

### Figures

- `figures/spatial_niche_adjacency_enrichment_heatmap.png`
- `figures/spatial_condition_and_adjacency_trends.png`
- `figures/tcga_braf_only_vulnerability_axis_heatmap.png`
- `figures/tcga_clinical_stage_axis_associations.png`
- `figures/gse151179_rai_direct_label_reanalysis.png`
- `figures/scrna_module_celltype_attribution_heatmap.png`

---

## 10. 최종 판단

이번 추가 분석 후 판단:

**GO 유지.**

가장 강한 보강:

1. Spatial niche adjacency pattern.
2. BRAF-only vulnerability heterogeneity.
3. Stage III/IV exploratory association with RAI-low/barrier-high axes.
4. scRNA cell-type attribution supporting marker logic.

가장 중요한 boundary:

1. RAI public direct labels are weak.
2. Spatial adjacency is ROI hypothesis, not independent-replicate inference.
3. Clinical associations are exploratory.
4. Drug/perturbation remains class-level validation hypothesis.

