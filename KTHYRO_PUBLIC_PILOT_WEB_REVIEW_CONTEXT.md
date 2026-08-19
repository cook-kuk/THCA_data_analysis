# K-Thyro Public Pilot 상황 판단용 문서

이 문서는 웹 기반 LLM, 외부 리뷰어, 삼성 제안서 검토자에게 현재 상황을 빠르게 판단시키기 위한 요약 파일입니다.

작성 위치:

`/home/seungho/personal/THCA_data_analysis/KTHYRO_PUBLIC_PILOT_WEB_REVIEW_CONTEXT.md`

분석 루트:

`/home/seungho/personal/THCA_data_analysis/kthyro_public_pilot`

---

## 1. 과제 개념

제안 과제명:

**K-Thyro Surgical-Spatial Theranostic Perturbation Atlas**

핵심 아이디어:

갑상선암 치료반응은 단일 driver mutation 또는 평균 bulk expression만으로 설명되지 않는다. 수술 전후 액체생검, 공간오믹스, 병리 AI, fresh tissue perturbation, drug-delivery/response imaging을 통합해 치료취약성 niche를 지도화한다.

관심 치료축:

1. RAI 분화상태 / RAI redifferentiation
2. HLA/APM immune visibility
3. HLA-low immune-invisible tumor state
4. CD8 exclusion
5. Myeloid/CAF barrier
6. Drug-delivery failure proxy
7. Aggressive dedifferentiation / proliferation

---

## 2. Public-data Pilot 목적

즉시 목표는 제안서를 쓰는 것이 아니라, 공개데이터만으로 핵심 가설이 말이 되는지 검증하는 것이었다.

검증 질문:

1. 갑상선암 치료취약성이 평균 bulk mutation만으로 설명되지 않는가?
2. RAI, HLA/APM, CD8 exclusion, myeloid/CAF, drug-delivery proxy가 분리 가능한 축인가?
3. Spatial transcriptomics에서 치료 관련 niche가 공간적으로 응집되어 보이는가?
4. Public drug/perturbation resource가 실험 후보 class를 제안할 수 있는가?
5. 이 결과가 삼성 제안서 preliminary evidence와 실험 설계로 연결되는가?

---

## 3. 완료된 작업

아래 디렉토리에 public-data pilot pipeline이 생성되고 실행 완료되었다.

`/home/seungho/personal/THCA_data_analysis/kthyro_public_pilot`

실행 스크립트:

`/home/seungho/personal/THCA_data_analysis/kthyro_public_pilot/scripts/run_all.sh`

추가 보강/감사 스크립트:

`/home/seungho/personal/THCA_data_analysis/kthyro_public_pilot/scripts/11_samsung_audit_strengthen.py`

주요 산출물:

- reproducible code pipeline
- public data manifest
- TCGA patient-level vulnerability table
- spatial spot-level vulnerability table
- slide-level niche summary
- cleaned drug/perturbation candidate table
- Samsung proposal-ready dark figures
- Korean preliminary result text
- Korean 10-slide storyline
- final executive summary

---

## 4. 사용 데이터

### TCGA-THCA bulk

- Public/local TCGA pancancer expression에서 THCA primary tumor 추출
- 최종 분석 샘플: **505 unique primary tumor samples**
- unique patient: **505**
- duplicated sample record: **0**

주요 파일:

`kthyro_public_pilot/results/tables/tcga_thca_patient_vulnerability_scores.tsv`

감사 결과:

`kthyro_public_pilot/results/reports/pilot_output_audit.md`

### Spatial transcriptomics

- GSE250521 thyroid spatial transcriptomics
- 최종 분석: **16 slides**
- spot 수: **57,144 spots**
- condition distribution:
  - normal: 4
  - PTC: 4
  - locally advanced PTC: 4
  - ATC: 4

중요한 해석 원칙:

**Spots are nested within slides. Spots must not be treated as independent biological replicates.**

주요 파일:

`kthyro_public_pilot/results/tables/spatial_spot_vulnerability_scores.tsv`

`kthyro_public_pilot/results/tables/spatial_slide_niche_summary.tsv`

`kthyro_public_pilot/results/tables/spatial_coherence_summary_cleaned.tsv`

### scRNA reference

- Local GSE193581 h5ad
- 67,678 cells
- 8 annotated cell types summarized

목적:

Spatial module이 어떤 cell type에서 주로 나오는지 sanity check.

중요한 해석:

HLA-II/CD74 signal은 tumor cell이 아니라 APC/myeloid/B cell에서 올 수 있다.

### DepMap/PRISM-derived drug resource

- Thyroid cell lines: 22
- Raw candidate rows: 140
- Cleaned candidate rows: 74
- compound + perturbation class duplicate: 0

중요한 해석:

Drug 이름을 최종 치료제로 주장하면 안 된다.  
**Perturbation class hypothesis**로만 사용해야 한다.

주요 파일:

`kthyro_public_pilot/results/tables/drug_pilot_candidate_rankings_cleaned.tsv`

---

## 5. TCGA 핵심 결과

TCGA-THCA 505명 primary tumor에서 vulnerability label이 아래처럼 분리되었다.

| Label | Count |
|---|---:|
| RAI-readable differentiated | 127 |
| Mixed/Other | 125 |
| Drug-delivery barrier-high | 105 |
| HLA-visible inflamed | 44 |
| CD8-excluded myeloid/CAF-high | 42 |
| HLA-low immune-invisible | 32 |
| RAI-low dedifferentiated | 30 |

해석:

갑상선암은 단일 driver mutation 또는 평균 bulk expression만으로 설명되지 않고, RAI, immune visibility, stromal exclusion, delivery proxy, aggressive dedifferentiation이 조합된 치료취약성 상태로 나뉠 가능성이 있다.

주의:

TCGA bulk는 spatial heterogeneity를 증명하지 못한다.  
TCGA 결과는 patient-level expression state evidence이다.

제안서용 안전 문장:

> TCGA-THCA 505명 공개 bulk transcriptome 분석에서 갑상선암은 단일 driver mutation 또는 평균 발현값으로 설명되지 않고, RAI 분화상태, HLA/APM 면역가시성, CD8 exclusion, stromal/drug-delivery barrier, proliferative dedifferentiation이 분리된 치료취약성 상태로 나뉘었다.

---

## 6. Spatial 핵심 결과

GSE250521 spatial transcriptomics 16 slides에서 치료취약성 niche label을 계산하고 KNN same-niche permutation을 수행했다.

결과:

**16/16 slides showed same-niche coherence z-score > 2.**

상위 coherence slide 예시:

| sample_id | condition | n_spots | same_niche_z | empirical_p |
|---|---|---:|---:|---:|
| GSM7980870_LPTC-3 | locally_advanced_PTC | 4355 | 43.1644 | 0.001996 |
| GSM7980864_PTC-1 | PTC | 2243 | 38.1139 | 0.001996 |
| GSM7980873_ATC-2 | ATC | 2656 | 37.2604 | 0.001996 |
| GSM7980869_LPTC-2 | locally_advanced_PTC | 4129 | 33.2437 | 0.001996 |
| GSM7980871_LPTC-4 | locally_advanced_PTC | 4688 | 28.8097 | 0.001996 |
| GSM7980866_PTC-3 | PTC | 4475 | 25.8125 | 0.001996 |
| GSM7980872_ATC-1 | ATC | 4672 | 23.9135 | 0.001996 |
| GSM7980875_ATC-4 | ATC | 2566 | 22.3509 | 0.001996 |

해석:

치료취약성 niche는 무작위 spot noise가 아니라 조직 내에서 공간적으로 응집된 expression state로 보인다.

주의:

이 결과는 임상반응 예측을 증명하지 않는다.  
Spatial transcriptomics는 peptide presentation이나 실제 drug delivery를 직접 증명하지 않는다.

제안서용 안전 문장:

> 공간전사체 pilot에서 치료취약성 niche는 무작위 spot noise가 아니라 조직 내에서 공간적으로 응집된 구조를 보였으며, 이는 FFPE/mIHC, GeoMx ROI, fresh tissue 기능실험으로 이어지는 수술-공간 theranostic platform의 필요성을 뒷받침한다.

---

## 7. Drug / Perturbation 해석

Drug 후보를 약 이름 중심으로 쓰면 위험하다.

제안서에서는 아래 class 중심으로만 사용해야 한다.

| Perturbation class | 목적 |
|---|---|
| RAI redifferentiation / MAPK-axis modulation | RAI-low/MAPK-high niche에서 iodide uptake 회복 검증 |
| HLA/APM restoration / immune visibility modulation | HLA-low tumor territory에서 HLA-I/B2M/TAP1 회복 검증 |
| Proliferation stress vulnerability | aggressive/proliferative niche 취약성 검증 |
| Myeloid/CAF barrier modulation | CD8 exclusion / stromal barrier 해소 검증 |
| Drug-delivery or nanoparticle distribution validation strategy | 실제 drug/liposome/nanoparticle penetration imaging |

주의할 후보:

- JAK inhibitor 계열: IFN-driven HLA/APM induction을 억제할 수도 있다. HLA restoration drug로 주장 금지.
- Mycophenolic acid: immunosuppressive 성격. immune restoration therapy로 주장 금지.
- AZD7762 / Alisertib: general proliferation toxicity일 수 있다.
- Selumetinib combination: RAI redifferentiation hypothesis일 뿐, iodide uptake assay 전까지는 효과 주장 금지.

안전 문장:

> Public drug resource는 최종 약물을 제시한 것이 아니라, RAI redifferentiation, HLA/APM rescue, proliferative stress, myeloid/CAF barrier modulation, delivery imaging으로 이어지는 perturbation class hypothesis를 제안한다.

---

## 8. GO / NO-GO 판정

최종 판정:

**GO**

근거:

1. TCGA-THCA 505명에서 5/5 main therapeutic vulnerability axis가 계산 가능하고 label이 분리됨.
2. GSE250521 spatial 16 slides에서 16/16 slides가 non-random same-niche coherence를 보임.
3. 분석 결과가 바로 FFPE/mIHC, GeoMx ROI, fresh tissue perturbation, iodide uptake, HLA/APM rescue, drug-delivery imaging 실험 가설로 이어짐.
4. 삼성 제안서용 preliminary figure와 report가 생성됨.

중요한 boundary:

Public-data pilot은 clinical deployment evidence가 아니라 hypothesis-generation evidence이다.

---

## 9. Proposal-ready Figure 파일

아래 파일들은 VSCode 또는 이미지 뷰어에서 바로 확인 가능하다.

```text
kthyro_public_pilot/results/figures/proposal/tcga_label_distribution_dark.png
kthyro_public_pilot/results/figures/proposal/tcga_vulnerability_axes_heatmap_dark.png
kthyro_public_pilot/results/figures/proposal/tcga_therapeutic_quadrant_dark.png
kthyro_public_pilot/results/figures/proposal/spatial_coherence_barplot_dark.png
kthyro_public_pilot/results/figures/proposal/spatial_niche_fraction_by_condition_dark.png
kthyro_public_pilot/results/figures/proposal/spatial_representative_maps_dark.png
kthyro_public_pilot/results/figures/proposal/integrated_public_pilot_evidence_matrix_dark.png
```

기존 figure 1-8도 아래에 있다.

```text
kthyro_public_pilot/results/figures/proposal/
```

---

## 10. Proposal-ready Report 파일

```text
kthyro_public_pilot/results/reports/pilot_output_audit.md
kthyro_public_pilot/results/reports/GO_NO_GO_DECISION.md
kthyro_public_pilot/results/reports/public_pilot_report.md
kthyro_public_pilot/results/reports/public_pilot_executive_summary_final_kr.md
kthyro_public_pilot/results/reports/samsung_preliminary_results_section_kr.md
kthyro_public_pilot/results/reports/pilot_10_slide_storyline_kr.md
kthyro_public_pilot/results/reports/tcga_preliminary_result_text_kr.md
kthyro_public_pilot/results/reports/spatial_preliminary_result_text_kr.md
kthyro_public_pilot/results/reports/claim_boundaries.md
kthyro_public_pilot/results/reports/experimental_validation_plan.md
```

---

## 11. Validation Panel

### Core panel

| Axis | Markers |
|---|---|
| Tumor/thyroid | PAX8, TG, TPO, NIS/SLC5A5 |
| HLA/APM | HLA-I/HLA-ABC, B2M, TAP1 |
| T cell | CD3, CD8, GZMB |
| Myeloid | CD68, CD163 |
| CAF/ECM | ACTA2/alphaSMA, FAP, COL1A1 |
| Checkpoint | PD-L1 |
| Hypoxia/delivery | CA9, VEGFA |

### Extended panel

| Axis | Markers |
|---|---|
| Thyroid differentiation | TSHR, FOXE1, NKX2-1 |
| HLA/APM regulation | NLRC5, PSMB8, PSMB9 |
| T cell exhaustion | PD-1 |
| Myeloid state | MRC1, SPP1 |
| CAF/ECM | POSTN |
| Hypoxia/metabolism | GLUT1/SLC2A1 |
| Proliferation | Ki-67 |

---

## 12. Experimental Validation Plan

1. FFPE cohort 구성
2. mIHC/IF panel validation
3. GeoMx ROI validation
4. Fresh tissue organoid/slice culture
5. Iodide uptake assay
6. HLA/APM rescue assay
7. PBMC/TIL co-culture if feasible
8. Fluorescent drug/liposome/nanoparticle distribution imaging

핵심 연결 구조:

**AI/spatial inference -> pathology validation -> GeoMx ROI -> fresh tissue perturbation -> functional assay**

---

## 13. 삼성 제안서 핵심 문장

### One-paragraph conclusion

공개데이터 pilot은 갑상선암 치료반응이 단일 driver mutation이나 평균 bulk expression으로 설명되지 않으며, RAI 분화상태, HLA/APM 면역가시성, CD8 exclusion, myeloid/CAF barrier, drug-delivery failure proxy가 환자 및 조직 공간 수준에서 분리된 치료취약성 축으로 나타남을 보여준다. 특히 GSE250521 공간전사체 16개 slide에서 치료취약성 niche가 비무작위적으로 응집되어, 본 과제의 핵심 가설인 “갑상선암 공간 치료취약성 생태계”가 공개데이터 수준에서 지지되었다.

### Final one-line message

본 과제는 갑상선암을 평균 유전체 질환이 아니라, 수술 perturbation과 공간 치료취약성 niche의 조합으로 치료반응이 결정되는 동적 생태계로 재정의한다.

---

## 14. Recommended PPT Insertion Order

1. Public-data pilot GO decision
2. Data sources and analysis flow
3. TCGA 505-patient vulnerability landscape
4. Patient-level subtype counts
5. Spatial ST 16-slide niche coherence
6. Representative spatial niche maps
7. Integrated evidence matrix
8. Drug/perturbation hypotheses
9. Experimental validation plan
10. Samsung 30억/3년 execution logic

---

## 15. 외부 리뷰어에게 물어볼 질문

아래 질문을 웹 LLM 또는 외부 reviewer에게 던지면 된다.

1. 이 public-data pilot 결과가 삼성육성과제 preliminary evidence로 충분히 강한가?
2. TCGA 505명과 spatial 16 slide 결과에서 과도하게 주장하는 부분이 있는가?
3. Spatial same-niche coherence 16/16 slides 결과를 제안서에서 어떻게 가장 안전하게 표현해야 하는가?
4. Drug/perturbation 후보를 약 이름이 아니라 class 중심으로 제시하는 전략이 적절한가?
5. 30억/3년 과제 범위에서 어떤 실험 validation이 가장 reviewer 설득력이 높은가?
6. 50억/5년 확장 논리로 Nature급 proof를 주장하려면 어떤 paired data가 반드시 필요한가?

---

## 16. Reviewer용 판단 요청 Prompt

아래 문단을 웹 LLM 또는 외부 reviewer에게 그대로 붙여넣어도 된다.

```text
You are a ruthless Samsung Future Technology Program reviewer and computational oncology validator.

Please judge whether the following public-data pilot is strong enough as preliminary evidence for a Samsung Future Technology Program proposal.

The proposed project is “K-Thyro Surgical-Spatial Theranostic Perturbation Atlas,” aiming to integrate perioperative liquid biopsy, spatial omics, pathology AI, and fresh tissue perturbation/drug-delivery imaging to map thyroid cancer therapeutic vulnerability niches.

Completed public-data pilot:
- TCGA-THCA bulk RNA-seq: 505 unique primary tumor samples.
- Patient-level therapeutic vulnerability labels:
  RAI-readable differentiated: 127
  Mixed/Other: 125
  Drug-delivery barrier-high: 105
  HLA-visible inflamed: 44
  CD8-excluded myeloid/CAF-high: 42
  HLA-low immune-invisible: 32
  RAI-low dedifferentiated: 30
- GSE250521 spatial transcriptomics: 16 slides, 57,144 spots.
- Condition distribution: normal 4, PTC 4, locally advanced PTC 4, ATC 4.
- KNN same-niche permutation: 16/16 slides showed same-niche coherence z-score > 2.
- DepMap/PRISM-derived public resource was cleaned into perturbation classes, not final drug claims.

Claim boundaries:
- Public data support hypothesis generation, not clinical deployment.
- TCGA bulk cannot prove spatial heterogeneity.
- Spatial transcriptomics supports spatial expression states, not direct peptide presentation.
- Drug-delivery failure is RNA proxy until fluorescent drug/liposome/nanoparticle imaging validates it.
- RAI-restorable score is hypothesis until iodide uptake assay validates it.
- Drug candidates must be framed as perturbation classes, not final therapies.

Question:
Is this sufficient for a GO decision as Samsung proposal preliminary evidence?
What are the strongest claims, weakest risks, and safest proposal wording?
```

---

## 17. 최종 판단

현재 상태는 **GO**이다.

다만 제안서에서 앞에 세울 것은 drug 후보가 아니다.

앞에 세울 메시지:

> TCGA 505명 + GSE250521 spatial 16 slide에서 치료취약성 niche가 분리되고, 그 niche가 공간적으로 응집된다.

Drug는 뒤에서:

> 이 niche를 재프로그래밍하기 위한 perturbation hypothesis

로 배치해야 한다.

