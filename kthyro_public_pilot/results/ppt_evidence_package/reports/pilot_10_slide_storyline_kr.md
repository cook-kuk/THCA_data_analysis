# Pilot 10-Slide Storyline (KR)

## 1. Public-data pilot GO decision
- One-line message: 공개데이터만으로도 K-Thyro 핵심 가설은 GO이다.
- Figure: `results/reports/GO_NO_GO_DECISION.md` 또는 `integrated_public_pilot_evidence_matrix_dark.png`
- Bullets: TCGA 505명; spatial 16 slide; 5/5 main axes; 16/16 slide coherent; claim은 hypothesis-generation.
- Speaker note: “임상 적용을 주장하는 것이 아니라, 삼성 과제로 검증할 충분한 신호가 있음을 보였습니다.”
- Claim boundary: public pilot, not deployment.

## 2. Data sources and analysis flow
- One-line message: bulk, spatial, scRNA, drug resource가 하나의 검증 흐름으로 연결된다.
- Figure: `figure1_public_pilot_data_flow_dark.png`
- Bullets: TCGA bulk; GSE250521 spatial; scRNA sanity check; DepMap/PRISM-derived perturbation classes; output tables.
- Speaker note: “기존 공개데이터를 사용해 신규 병원 코호트 투입 전 위험을 낮췄습니다.”
- Claim boundary: public data have no paired perturbation assay.

## 3. TCGA 505-patient vulnerability landscape
- One-line message: 치료취약성은 평균 bulk 하나가 아니라 다축 landscape이다.
- Figure: `tcga_vulnerability_axes_heatmap_dark.png`
- Bullets: 505 primary tumors; five main axes; axes are correlated but not identical; patient-level labels generated.
- Speaker note: “driver mutation 중심 분류를 넘어 치료 intervention이 가능한 상태축을 정의했습니다.”
- Claim boundary: bulk cannot prove spatial heterogeneity.

## 4. Patient-level subtype counts
- One-line message: vulnerability labels are large enough to motivate prospective validation.
- Figure: `tcga_label_distribution_dark.png`
- Bullets: label counts {'RAI-readable differentiated': 127, 'Mixed/Other': 125, 'Drug-delivery barrier-high': 105, 'HLA-visible inflamed': 44, 'CD8-excluded myeloid/CAF-high': 42, 'HLA-low immune-invisible': 32, 'RAI-low dedifferentiated': 30}; Mixed/Other remains intermediate; label is exploratory inference; clinical association is secondary.
- Speaker note: “제안서에는 label count를 preliminary evidence로, clinical prediction은 목표로 둡니다.”
- Claim boundary: labels are not clinical subtypes yet.

## 5. Spatial ST 16-slide niche coherence
- One-line message: therapeutic niches are spatially organized, not random spot noise.
- Figure: `spatial_coherence_barplot_dark.png`
- Bullets: GSE250521 16 slides; KNN same-niche permutation; 16/16 z > 2; spot nested within slide.
- Speaker note: “공간분석을 해야 하는 이유가 여기서 나옵니다.”
- Claim boundary: no clinical-response prediction.

## 6. Representative spatial niche maps
- One-line message: RAI, HLA/APM, CAF/myeloid, niche labels show distinct territories.
- Figure: `spatial_representative_maps_dark.png`
- Bullets: representative PTC/LPTC/ATC; four-panel map; tumor/immune/stromal axes separate; ROI validation targets.
- Speaker note: “병리 ROI와 기능실험 위치를 지정할 수 있는 지도가 됩니다.”
- Claim boundary: expression states, not protein or peptide presentation.

## 7. Integrated evidence matrix
- One-line message: each axis has different evidence strength and validation requirement.
- Figure: `integrated_public_pilot_evidence_matrix_dark.png`
- Bullets: RAI and HLA axes moderate/strong; delivery proxy weak; scRNA supports cell-type interpretation; validation assays specified.
- Speaker note: “강한 주장과 약한 주장을 분리해 reviewer risk를 낮춥니다.”
- Claim boundary: matrix is proposal-readiness, not final proof.

## 8. Drug/perturbation hypotheses
- One-line message: drug names are secondary; perturbation classes are the proposal logic.
- Figure: `figure6_drug_perturbation_candidate_pilot.png` plus `drug_pilot_candidate_rankings_cleaned.tsv`
- Bullets: MAPK-axis; epigenetic/APM modulation; proliferation stress; barrier modulation; delivery imaging.
- Speaker note: “약효 주장 대신 functional rescue assay로 닫겠습니다.”
- Claim boundary: DepMap/PRISM nominate hypotheses only.

## 9. Experimental validation plan
- One-line message: AI-spatial inference is closed by FFPE, GeoMx, and fresh tissue functional assays.
- Figure: `figure7_samsung_experimental_validation_design.png`
- Bullets: core marker panel; extended panel; iodide uptake; HLA/APM rescue; fluorescent delivery imaging.
- Speaker note: “수술 검체를 바로 기능실험으로 연결하는 것이 과제의 차별점입니다.”
- Claim boundary: validation required before clinical use.

## 10. Samsung 30억/3년 execution logic
- One-line message: 3년 과제는 public signal을 clinical-spatial-functional proof로 전환한다.
- Figure: `figure8_pilot_conclusion_slide.png`
- Bullets: Year 1 FFPE/mIHC; Year 2 GeoMx/fresh tissue assay; Year 3 integrated AI-spatial theranostic model; 50억/5년 확장 가능.
- Speaker note: “30억은 platform proof, 50억은 atlas-scale clinical deployment pre-stage입니다.”
- Claim boundary: deployment is beyond pilot.
