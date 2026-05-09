# K-Thyro Public Pilot Analysis Flow

## 1. Overall Logic

질문:

> 갑상선암 치료취약성이 평균 driver mutation만으로 설명되는가, 아니면 RAI 분화상태, HLA/APM 면역가시성, CD8 exclusion, myeloid/CAF barrier, drug-delivery failure proxy가 별도의 공간 niche로 존재하는가?

분석 구조:

1. TCGA bulk에서 patient-level therapeutic vulnerability axes 계산.
2. GSE250521 spatial transcriptomics에서 spot-level axes와 niche label 계산.
3. Slide-level spatial coherence로 niche가 random spot noise인지 검정.
4. External thyroid cohorts에서 동일 K-Thyro gene set으로 axis reproducibility 검증.
5. TCGA driver group 내부 heterogeneity로 mutation-only insufficiency 정량화.
6. scRNA/GeoMx로 cell-type/ROI-level sanity check.
7. DepMap/PRISM로 final drug가 아닌 perturbation class hypothesis 도출.
8. 모든 결과를 experimental validation plan으로 연결.

---

## 2. Stage A: TCGA-THCA Bulk Scoring

### Input

- Public TCGA-THCA primary tumor expression.
- Clinical/driver metadata where available.
- 최종 분석 샘플: **505 primary tumors**.

### Scores

Gene set mean-z scoring으로 다음 module score 계산:

- RAI differentiation
- MAPK/dedifferentiation stress
- HLA-I/APM
- HLA-II/APC niche
- CD8/cytotoxic T cell
- IFN activation
- Checkpoint/suppression
- Myeloid/TAM
- CAF/ECM/TGF-beta barrier
- Vascular/drug-delivery proxy
- Hypoxia/metabolism
- Proliferation
- Tumor epithelial/thyroid lineage

Derived axes:

- immune_visibility_score
- hla_low_invisible_score
- cd8_exclusion_proxy
- rai_restorable_score
- drug_delivery_failure_proxy
- aggressive_dedifferentiation_score
- myeloid_caf_barrier_score

### Outputs

- `tables/tcga_thca_patient_vulnerability_scores.tsv`
- `figures/proposal/tcga_vulnerability_axes_heatmap_dark.png`
- `figures/proposal/tcga_label_distribution_dark.png`
- `figures/proposal/tcga_therapeutic_quadrant_dark.png`

### Main Observed Result

TCGA 505명은 다음 vulnerability label로 분리됨:

- RAI-readable differentiated: 127
- Mixed/Other: 125
- Drug-delivery barrier-high: 105
- HLA-visible inflamed: 44
- CD8-excluded myeloid/CAF-high: 42
- HLA-low immune-invisible: 32
- RAI-low dedifferentiated: 30

### Interpretation

Observed data:

- Bulk expression에서 치료취약성 관련 module score가 환자별로 다르게 나타남.

Inferred score:

- RAI, immune visibility, CD8 exclusion, CAF/myeloid barrier, drug-delivery failure proxy 등.

Exploratory hypothesis:

- 같은 thyroid cancer라도 치료취약성 생태계가 다르므로, 평균 driver mutation만으로 치료 전략을 결정하기 어렵다.

Validation needed:

- FFPE/mIHC, GeoMx/ROI, iodide uptake, HLA/APM rescue, drug penetration imaging.

---

## 3. Stage B: Spatial Transcriptomics Niche Mapping

### Input

- Public GSE250521 spatial transcriptomics.
- Slides: **16**
- Spots: **57,144**
- Conditions:
  - normal 4
  - PTC 4
  - locally advanced PTC 4
  - ATC 4

### Processing

- Spot QC.
- Library-size normalization/log transform.
- Same K-Thyro gene set scoring per spot.
- Quantile threshold 기반 niche label:
  - RAI-restorable niche
  - RAI-low dedifferentiated niche
  - HLA-visible inflamed niche
  - HLA-low invisible tumor niche
  - CD8-excluded / myeloid-CAF niche
  - APC-rich niche
  - Drug-delivery failure proxy niche
  - Mixed/Other

### Spatial Coherence Test

- KNN same-niche enrichment.
- 각 slide에서 random label permutation과 비교.
- Spots는 독립 biological replicate로 취급하지 않음.

### Outputs

- `tables/spatial_spot_vulnerability_scores.tsv`
- `tables/spatial_slide_niche_summary.tsv`
- `tables/spatial_coherence_summary_cleaned.tsv`
- `figures/proposal/spatial_representative_maps_dark.png`
- `figures/proposal/spatial_coherence_barplot_dark.png`

### Main Observed Result

- **16/16 slide에서 same-niche coherence z-score > 2**
- z-score min/median/max: **6.09 / 22.00 / 43.16**

### Interpretation

Observed data:

- Spot-level expression states and niche labels show spatial organization.

Inferred score:

- Spatial RAI, HLA/APM, CD8/cytotoxic, CAF/myeloid, drug-delivery proxy.

Exploratory hypothesis:

- 치료취약성은 bulk average가 아니라 조직 내 territory/niche 단위로 존재할 가능성이 있다.

Validation needed:

- Protein-level spatial validation and functional assays.

---

## 4. Stage C: External Public Validation Expansion

### Input

Exact K-Thyro rescoring on external GPL570 thyroid cohorts:

- GSE29265
- GSE33630
- GSE53157
- GSE65144

Total: **206 samples**

### Processing

- 기존 local gene-level expression matrix 사용.
- 동일 `config/gene_sets.yaml` 기반 mean-z scoring.
- Contrast tests:
  - ATC vs PTC
  - ATC vs normal
  - PTC vs normal
  - PDTC vs PTC
  - PDTC vs normal
  - 기타 가능한 histology contrast

### Outputs

- `tables/exact_kthyro_external_bulk_scores.tsv`
- `tables/exact_kthyro_external_bulk_contrasts.tsv`
- `figures/exact_kthyro_external_bulk_axis_heatmap.png`
- `figures/exact_kthyro_external_bulk_boxplots.png`

### Main Result

- 168 contrast tests.
- strong: 65
- moderate: 12
- weak: 57
- weak_or_none: 32
- opposite: 2
- Prespecified direction: 42 match / 2 opposite.

### Interpretation

External GPL570 cohorts support reproducibility of axis directionality, especially:

- ATC shows lower RAI differentiation than PTC/normal.
- ATC shows higher aggressive/dedifferentiation and barrier/drug-delivery-failure proxy in strong cohorts.

Claim boundary:

- External bulk validates axis directionality/separability, not treatment response or spatial niche biology.

---

## 5. Stage D: Public Validation Matrix v2

### Included Evidence Layers

- TCGA bulk scores.
- GSE250521 spatial coherence.
- External GPL570 exact K-Thyro rescoring.
- Previous external bulk expression tests.
- PDTC/ATC bulk layer.
- GeoMx ROI spatial layer.
- Korean scRNA pseudobulk layer.
- RAI before/after and direct RAI avid/refractory public labels.
- Proteomics dedifferentiation trend.
- DepMap/PRISM drug/resource layer.

### Result

Total validation tests: **301**

- strong: 144
- moderate: 30
- weak: 73
- weak_or_none: 52
- opposite: 2

### Interpretation

This is not a single definitive validation. It is a multi-layer feasibility matrix showing which therapeutic axes are supported, weak, or still need experimental validation.

---

## 6. Stage E: Axis Separability

### Question

Are the axes truly separable, or are they just the same immune/stromal signal renamed?

### Analysis

- Pairwise Spearman correlations among main axes.
- TCGA patient-level.
- External GPL570 sample-level.
- GSE250521 spot-level within each slide.

### Result

- Axis-pair tests: **330**
- TCGA median abs rho: **0.609**
- External GPL570 median abs rho: **0.637**
- Spatial slide-level spot correlations were lower; no slide had abs rho >= 0.75 among main axis pairs.

### Interpretation

Axes are biologically related but not collapsed into one score. Spatial data show stronger local separability, supporting niche-level logic.

---

## 7. Stage F: Mutation-Only Insufficiency

### Question

Does BRAF/RAS mutation alone explain therapeutic vulnerability state?

### Result

- BRAF: 274 patients, 7 vulnerability labels.
- BRAF largest label: Drug-delivery barrier-high, 91/274 = 33.2%.
- BRAF mutation-only ambiguity: 66.8%.
- Driver-label Cramer’s V: 0.464, p = 8.05e-40.

Driver group variance explained:

- Drug-delivery failure proxy eta² 0.332.
- RAI differentiation eta² 0.309.
- HLA-I/APM eta² 0.156.
- Immune visibility eta² 0.125.

### Interpretation

Driver mutation is associated with vulnerability state, but does not determine it. This directly supports the need for AI-spatial therapeutic vulnerability atlas.

---

## 8. Stage G: scRNA and GeoMx Sanity Check

### scRNA

- GSE193581 h5ad.
- 67,678 cells.
- 8 cell types.

Use:

- Confirm which modules are epithelial, fibroblast, endothelial, myeloid, T/B cell enriched.

Caveat:

- HLA-I/APM coverage in this processed scRNA object was weak/zero for some genes, so scRNA is supportive for cell-type context, not definitive HLA-I/APM validation.

### GeoMx

Use:

- ROI-level spatial sanity check.
- Helps distinguish epithelial/stromal/APC-rich territories.

Caveat:

- ROI labels and selection matter; not equivalent to Visium spot-level map.

---

## 9. Stage H: Drug/Perturbation Pilot

### Input

- Public DepMap/CCLE/PRISM-derived local resources.
- Thyroid cancer cell lines where available.

### Output

- Candidate table cleaned and deduplicated.
- Perturbation classes, not final drugs.

### Correct Interpretation

Use as:

- MAPK-axis modulation for RAI redifferentiation hypothesis.
- Epigenetic/APM/IFN-related perturbation for HLA/APM rescue hypothesis.
- Cell-cycle stress for proliferative/aggressive niche.
- Myeloid/CAF barrier modulation as literature/proxy-guided class.
- Drug-delivery validation with fluorescent drug/nanoparticle imaging.

Do not claim:

- Any compound is a proven clinical therapy from this pilot.
- JAK inhibitors restore HLA/APM.
- Mycophenolic acid is immune-restorative.
- Cell-line PRISM validates CAF/myeloid barrier or drug delivery.
