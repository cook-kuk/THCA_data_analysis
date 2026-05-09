# 2026-05-09 HLA 두 논문 Nature Communications급 합성 브리프

## 결론

두 논문은 같이 키우되, **절대 같은 주장으로 섞으면 안 된다.**

1. **Paper 4: Graves/AITD HLA**
   - 지금 가장 논문-ready에 가깝다.
   - 핵심은 `DPB1*05:01`, `B*46:01`, `A*02:07`, `C*01:02` 중심의 Pan-Asian Graves/AITD HLA 구조.
   - 현재 데이터만으로도 meta/resource/synthesis 논문은 강하다.
   - Nature Communications급으로 올리는 병목은 **Korean adult GD germline NGS HLA case-control 신규 코호트 부재**다.

2. **Paper 2: HT-overlap PTC HLA/immune-context**
   - allele association 논문으로 쓰면 약하다.
   - 대신 **HT-overlap/immune-active PTC의 antigen-presentation, IFN-gamma, TLS, HLA-II expression-module biology**로 쓰면 강하다.
   - allele-level 후보는 `DRB1*04:05` 하나가 제일 낫지만, 현재 n=9 vs 9라서 “검증된 발견”이 아니라 prospective validation target이다.

## Paper 4 핵심 가설

| ID | 가설 | 현재 검증 |
|---|---|---|
| GD-H1 | `DPB1*05:01`은 Pan-Asian Graves/AITD의 가장 강한 class-II risk allele이다. | v3 pooled OR 2.10 (1.70-2.60), k=4, I2=55.16%; manuscript-grade caveat 가능. |
| GD-H2 | `B*46:01`은 진짜 신호지만 population heterogeneity가 크다. | OR 2.17 (1.39-3.36), I2=85.09%, PI 0.52-8.94. single-effect 과장 금지. |
| GD-H3 | `A*02:07`은 Korean literature extraction 후 full Pan-Asian meta-ready로 승격됐다. | OR 2.12 (1.73-2.61), grade B -> A. |
| GD-H4 | `C*01:02`은 cross-autoimmune panel에서 가장 GD-specific하다. | specificity score 14.97, AITD breadth ratio 14.97. |
| GD-H5 | `DPB1*05:01`은 GD-specific이 아니라 AITD-broad allele이다. | GD specificity 2.89, AITD breadth ratio 3.31, outside-GD strong hit 2개. |
| GD-H6 | `DPB1*05:01`은 `B*46:01/A*02:07` haplotype passenger만은 아니다. | K2 LD: DPB1*05:01 vs B*46:01 r2≈0.077; B*46:01 vs A*02:07 r2≈0.337. |
| GD-H7 | FinnGen/Pan-UKBB는 MHC class-II thyroid-autoimmunity signal을 trans-ancestry로 뒷받침한다. | FinnGen GD rs9277534 beta +0.129, p=3.6e-07; HT-AI beta -0.107, p=7.1e-42. effect-size transfer는 금지. |
| GD-H8 | healthy Korean HLA atlas가 arcasHLA population-frequency 사용을 보강한다. | 9 populations, 145 alleles; Korean platform concordance r≈0.86-0.93. |
| GD-H9 | `DRB1*07:01/DQB1*02:01` protective axis는 보이지만 allele별 증거 강도가 다르다. | DRB1*07:01 OR 0.43, k=2; DQB1*02:01은 아직 single-source. |
| GD-H10 | Korean adult GD NGS cohort가 NComm급 업그레이드의 핵심이다. | 현재 meta는 강하지만 de novo Korean discovery claim은 아직 불가. |
| GD-H11 | AITD 조직에서는 HLA-II/AP, CD74/MIF, TLS/B-cell 축이 공간적으로 확장된다. | GSE248205 AITD vs control: CD74/MIF d=3.62, B/TLS d=2.57, AP/TLS d=2.26; HT spot burden above control p90 HLA-II/AP 0.964. |
| GD-H12 | 독립 older-array HT 조직에서도 HLA/AP 축이 매우 강하다. | GSE29315 HT vs hyperplasia: HLA-II/AP d=5.45, exact p=0.0007, FDR=0.0023; HLA-I d=4.40. |

## Paper 2 핵심 가설

| ID | 가설 | 현재 검증 |
|---|---|---|
| PTC-H1 | HT-overlap PTC allele 후보는 `DRB1*04:05`가 1순위다. | OR 3.93 (1.06-14.53), p=0.0866, FDR=0.346; n/group for OR=1.8 ≈ 633. |
| PTC-H2 | `DPB1*05:01`은 현재 HT-PTC enrichment allele로 지지되지 않는다. | OR 1.08 (0.40-2.90), p=1.0; carrier OR 0.65. |
| PTC-H3 | GSE286332 within-cohort HT vs non-HT allele comparison은 유의하지 않다. | all FDR=1.0. |
| PTC-H4 | TCGA-THCA에서 DM1/immune-active state는 HLA-I/II module과 양의 상관이다. | HLA-I rho 0.318, p=8.3e-14; HLA-II rho 0.393, p=7.3e-21. |
| PTC-H5 | DM1-HLA 관계는 purity/stromal artifact가 아니다. | immune+stromal partial rho: HLA-I 0.266, HLA-II 0.390. |
| PTC-H6 | IFN-gamma가 DM1 -> HLA induction의 주요 mediator다. | IFNG vs HLA-II rho 0.873; Hallmark-mediated fraction HLA-II 70.22%. |
| PTC-H7 | HLA/IFN activation은 RAI differentiation module과 반대로 움직인다. | DM1 vs RAI rho -0.886; HLA-II vs RAI rho -0.475. |
| PTC-H8 | pan-cancer에서도 DM1-HLA co-induction은 lineage-portable하다. | 10,996 samples, 32 lineages; sign coherence HLA-I 100%, HLA-II 96.88%. |
| PTC-H9 | scRNA/spatial은 HLA-II가 APC compartment와 tumor-stage acquisition에 걸려 있음을 보여준다. | Pu 19,999 cells, Lu 17,898 cells; DC/Myeloid/B가 HLA-II dominant. |
| PTC-H10 | TLS niche는 HLA-II-hot spatial niche다. | 608 niches; global Cohen d 0.655; PTC d 0.935, ATC d 0.916. |
| PTC-H11 | somatic HLA-presentation LoF는 DM1/HLA state의 주 원인이 아니다. | class-I LoF 0건, nonsilent 1건. |
| PTC-H12 | HLA-G는 올라가지만 selective HLA-G escape 모델은 아니다. | HLA-G vs DM1 rho 0.428; HLA-G decoupling vs DM1 rho 0.006, p=0.887. |
| PTC-H13 | BRAF-rich DM1-high/IFN-low fail-to-induce subgroup이 있다. | 97/527; BRAF OR 3.78, FDR=3.2e-07; survival은 약함. |
| PTC-H14 | neoantigen/TMB는 보조 설명이고 중심축은 아니다. | TMB vs DM1 rho 0.252; DM1 -> HLA-II direct beta 0.410, p=6.0e-22. |
| PTC-H15 | AIRE/thymic tolerance와 peptide-binding은 functional follow-up 가설이다. | HPA thymus AIRE rank 1, GTEx thyroid donors 653, peptidomics present. |
| PTC-H16 | ICI response에서 HLA-I expression은 modest predictor, HLA-II는 약하다. | HLA-I meta OR 1.35, p=0.0143; HLA-II OR 1.02, p=0.839. |
| PTC-H17 | prospective HT-PTC germline NGS cohort가 allele-level 검증 병목이다. | 현재 explicit HT-overlap allele test는 n=9 vs 9. |
| PTC-H18 | 독립 bulk HT/PTC와 scRNA HT/PTC에서 HLA/AP-TLS elevation이 재현된다. | GSE138198 PTCwithHT vs PTCwithoutHT: T/IFNG d=1.23, HLA-II/AP d=1.05; GSE163203 PTCwithHT vs PTCwithoutHT: AP/TLS d=2.78, B/TLS d=2.59. |
| PTC-H19 | Korean bulk RNA-seq는 HLA/AP expression ecology가 대형 독립 코호트에서도 강함을 보인다. | GSE213647 PTC tumor vs PTC normal: HLA-II/AP d=1.39 FDR=4.03e-40, HLA-I d=1.22 FDR=5.08e-36, Myeloid/DC d=1.23 FDR=3.34e-36. |
| PTC-H20 | thyroid Visium progression은 HLA/AP/CD74/myeloid spatial acquisition을 보인다. | GSE250521 N/PTC/LPTC/ATC stage trend: CD74/MIF rho=0.76 q=0.0027, HLA-I rho=0.75 q=0.0027, HLA-II/AP rho=0.61 q=0.018. |
| PTC-H21 | 모든 dataset이 positive일 필요는 없고, negative stress-test가 specificity를 만든다. | GSE6004 invasion vs normal은 HLA/AP/TLS 증가가 없음: HLA-II/AP d=-0.42, AP/TLS d=-0.59, FDR=0.829. GSE184362/GSE191288은 HT label 부재로 low-tier supplement만 가능. |

## 지금 쓰면 되는 강한 문장

**Paper 4:**  
Pan-Asian Graves/AITD HLA architecture is centered on DPB1*05:01 with a partly independent B*46:01/A*02:07/C*01:02 class-I axis; C*01:02 appears most GD-specific, whereas DPB1*05:01 is AITD-broad rather than GD-exclusive, and independent AITD tissue datasets support a class-II antigen-presentation/CD74/TLS mechanism layer.

**Paper 2:**  
HT-overlap and immune-active PTC are best framed as a transcriptomic antigen-presentation/TLS/IFN-gamma state; HLA-II alleles are prospective candidates only, while HLA-I/II expression modules are supported across explicit HT-overlap bulk/scRNA cohorts, Korean bulk RNA-seq, spatial progression, and TLS-HLA-II niches.

## 절대 금지

- HLA allele이 thyroid cancer risk/outcome을 예측한다고 쓰지 말 것.
- TCGA/K2/Lee cancer cohort HLA allele table을 survival, RAI, BRAF/RAS/TERT, DM1/DM2, stage, LN, recurrence와 join하지 말 것.
- `DPB1*05:01`을 HT-PTC allele finding이라고 부르지 말 것.
- `DRB1*04:05`를 validated finding이라고 부르지 말 것.
- HLA-G를 selective immune escape라고 과장하지 말 것.

## 생성 산출물

- `project/results/hla_two_paper_synthesis_2026_05_09/HLA_TWO_PAPER_NCOMM_SYNTHESIS.md`
- `project/results/hla_two_paper_synthesis_2026_05_09/hypothesis_validation_matrix.tsv`
- `project/results/hla_two_paper_synthesis_2026_05_09/paper_priority_table.tsv`
- `project/reports/2026_05_09_HLA_TWO_PAPER_NCOMM_SYNTHESIS.md`
- `project/reports/2026_05_09_HLA_TWO_PAPER_NCOMM_SYNTHESIS_KR.md`
- `project/reports/2026_05_09_HLA_DATA_HUNT_WAVE3_KR.md`
