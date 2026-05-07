# Paper 2 HLA — 교수님 보고용 정밀 패킷 (KR)
**날짜:** 2026-05-06
**작성자:** Seungho Cook
**분류:** Paper 2 Pillar I HLA · 탐색적(exploratory) 후보 우선순위 지도
**원칙:** 과대주장(overclaim) 금지. 본 자료는 *causal HLA susceptibility / clinical risk prediction / patient selection*의 근거가 **아닙니다**.

---

## 0. 한 줄 요약

한국인 PTC 풀(n = 874명, 6개 HLA 후보 대립유전자)의 1차 forest는
**"PTC carrier 빈도 vs 공개 baseline allele 빈도(2n)"**라는 metric 불일치 비교였으며,
**엄격한 metric 정합화(직접 allele 대 allele 비교)** 후에는 6개 중 4개의 enrichment 신호가 사라집니다.
나머지 2개(C*01:02 depletion, DQB1*02:01 depletion)만 BH-q < 0.05로 살아남습니다.
따라서 본 결과는 **matched-control NGS validation이 필요한 후보 우선순위 지도**이지 case-control association이 아닙니다.

---

## 1. 핵심 caveat (반드시 먼저 명시)

| 축 | PTC측 | Baseline측 |
|---|---|---|
| 단위 | individual 단위 carrier 수 | 출판된 allele frequency × 2n |
| 분모 | n_individuals = 874 | 2n_alleles |
| 사건 정의 | 적어도 한 copy를 가진 사람 | allele 수 |

**1차 분석은 이 두 단위를 직접 비교**했으므로 OR/Fisher/forest 결과는 **탐색적(exploratory) hypothesis-generating** 수준에서만 의미가 있습니다.
HWE 가정 하에서 carrier 빈도와 allele 빈도는 환산 가능하지만, 실제 분포가 HWE를 따른다는 보장은 없습니다.

| 환산식 | 의미 |
|---|---|
| `carrier_freq = 1 − (1 − AF)²` | baseline AF → 기대 carrier 빈도 |
| `AF = 1 − √(1 − carrier_freq)` | PTC carrier → 함의된 AF |

본 패킷에서는 위 환산을 모두 적용한 4개 시나리오(S0-S3)와 PTC 실제 allele dosage count로부터 직접 계산한 S3(직접 allele vs allele)를 동시에 보고합니다.

---

## 2. 4-시나리오 metric sensitivity 결과

각 시나리오별 OR(Haldane 보정) + 2-sided Fisher exact p, n=6 후보 내 BH-q.

### 2.1 시나리오 정의

- **S0 primary** — PTC carrier 수(개인 단위) vs baseline allele 수(2n). 단위 불일치, 1차 forest의 기준.
- **S1** — baseline AF를 HWE 가정 하 기대 carrier 빈도로 환산 후 비교 (단위: carrier).
- **S2** — PTC carrier 빈도를 HWE 가정 하 함의된 AF로 환산 후 비교 (단위: AF).
- **S3 direct** — PTC 실제 allele dosage count(2n)와 baseline allele count(2n)을 직접 비교. **HWE 가정 불필요**, 가장 robust한 sanity check.

### 2.2 4-시나리오 forest (S0 → S3)

`F16_metric_sensitivity_forest.png` — 6 alleles × 4 시나리오 horizontal forest.

| Allele | S0 OR | S1 OR | S2 OR | S3 direct OR | 방향 안정성 |
|---|---|---|---|---|---|
| A*02:07 | 2.43★ | 1.20 | 1.20 | **1.26** (q=0.43) | enrichment_stable |
| B*46:01 | 2.16★ | 1.06 | 1.04 | **1.05** (q=0.80) | enrichment_stable |
| C*01:02 | 1.47★ | **0.66**◆ | **0.69**◆ | **0.68** (q=7×10⁻⁴) | **direction_flips** |
| DPB1*05:01 | 1.96★ | **0.76**◆ | **0.80**◆ | **0.98** (q=0.80) | **direction_flips** |
| DQB1*02:01 | 0.026 (zero) | 0.013 | 0.013 | **0.018** (q=5×10⁻⁸) | depletion_stable |
| DRB1*07:01 | 1.73★ | **0.84** | **0.84** | **0.84** (q=0.43) | **direction_flips** |

★ = S0에서 enrichment 방향, ◆ = harmonization 후 depletion 방향으로 flip.

### 2.3 BH-q (n=6 candidates) per scenario

| Allele | S0 q | S3 direct q |
|---|---|---|
| A*02:07 | **1.1×10⁻⁵** | 0.43 |
| B*46:01 | **1.1×10⁻⁵** | 0.80 |
| C*01:02 | **4.3×10⁻⁴** | **7.4×10⁻⁴** ✓ |
| DPB1*05:01 | **9.9×10⁻¹⁴** | 0.80 |
| DQB1*02:01 | **2.5×10⁻⁶** | **5.3×10⁻⁸** ✓ |
| DRB1*07:01 | **4.3×10⁻⁴** | 0.43 |

**S0(primary)에서 6/6 q<0.05 → S3(직접 allele vs allele)에서 2/6 q<0.05만 유지.**

---

## 3. Subcohort heterogeneity

`F17_subcohort_frequency_heatmap.png` — 6 alleles × 3 subcohorts (K2, Lee2024, GSE286332_PTC).

| Allele | min cohort / freq | max cohort / freq | χ² p (3 cohorts) | Flag |
|---|---|---|---|---|
| A*02:07 | Lee2024 / 7.5% | GSE286332_PTC / 33.3% | **0.016** ★ | small_subcohort + heterog. |
| B*46:01 | Lee2024 / 10.0% | GSE286332_PTC / 33.3% | 0.073 | small_subcohort |
| C*01:02 | K2 / 18.3% | GSE286332_PTC / 55.6% | **0.006** ★ | dominated by Lee2024; heterog. |
| DPB1*05:01 | Lee2024 / 52.1% | K2 / 56.2% | 0.554 | (균질) |
| DQB1*02:01 | K2 / 0% | K2 / 0% | n/a | zero cells |
| DRB1*07:01 | K2 / 9.8% | GSE286332_PTC / 33.3% | 0.084 | small_subcohort |

GSE286332_PTC는 n = 18로 매우 작습니다. extreme value(33%, 55%)는 N에 의한 분산일 수 있어 **본 분석의 주력 신호로 사용하지 않습니다.**

---

## 4. DPB1*05:01 source-sensitivity 노트

DPB1*05:01의 baseline은 다른 5개와 달리 In JW 2015 (Ann Lab Med, n=613) 가 아닌 **AFND South Korea pool (n=680)** 에서 가져왔습니다.

- AFND는 여러 개별 출판물의 메타-pool로 typing protocol/해상도가 출처마다 상이.
- Jung 2023 Korean 데이터에서는 DPB1*05:01:01이 35.1%로 보고됨 (현재 baseline 36.7%와 매우 근접).
- 따라서 DPB1*05:01의 1차 forest 신호는 baseline 출처 변경에 따라 변동 가능성이 큽니다.
- **S3 직접 allele 비교에서 OR=0.98 (q=0.80) — 사실상 PTC와 baseline이 동일** → 1차 forest의 OR=1.96 enrichment는 metric 불일치 + 제한적 baseline pool의 합동 artifact.

---

## 5. DQB1*02:01 zero-cell stress 노트

PTC 풀(n=874 중 callable 631)에서 DQB1*02:01 carrier가 0명 관측 → S0/S3 모두 강한 depletion 방향.

- Haldane-Anscombe (+0.5) 보정 OR = 0.018 (95% CI 0.001-0.295), Fisher p = 8.9×10⁻⁹.
- **그러나** zero observation은:
  1. 실제 생물학적 결손일 수 있음(antigen presentation/AITD context와 일관),
  2. arcasHLA imputation pipeline의 specific 호출 실패일 수도 있음 (technical artifact),
  3. typing 해상도 차이(2-digit vs 4-digit)로 인한 분류 충돌일 수 있음.
- **단일 cohort + zero observation은 stand-alone 결론을 지지하지 않습니다.** 매칭된 healthy Korean control NGS에서 동일 0% 또는 매우 낮은 carrier 빈도가 재현되어야 신뢰 가능.

---

## 6. Claim grade per allele — v2 (F18_claim_boundary_matrix.png)

> **v2 core statement.** S0 carrier-vs-allele forest는 6/6 significant 신호를 보였지만, strict S3 allele-vs-allele harmonization 후에는 **C*01:02, DQB1*02:01** 두 개 depletion 후보만 남습니다. **어떤 enrichment 후보도 strict metric harmonization을 통과하지 못합니다.**

| Allele | S0 primary q | S3 direct q | S0 direction | S3 direction | Claim grade (v2) | 다음 단계 우선순위 |
|---|---|---|---|---|---|---|
| A*02:07 | 1.1e-5 | 0.43 (n.s.) | enrichment | enrichment | **S0_only_exploratory** | medium |
| B*46:01 | 1.1e-5 | 0.80 (n.s.) | enrichment | enrichment | **S0_only_exploratory** | medium |
| C*01:02 | 4.3e-4 | **7.4e-4 ✓** | enrichment | **depletion** | **survives_strict_metric_but_direction_flips** | **high** |
| DPB1*05:01 | 9.9e-14 | 0.80 (n.s.) | enrichment | depletion | **source_sensitive_artifact_likely** | medium |
| DQB1*02:01 | 2.5e-6 | **5.3e-8 ✓** | depletion | depletion (zero-cell) | **survives_strict_metric_but_zero_cell_caution** | **high** |
| DRB1*07:01 | 4.3e-4 | 0.43 (n.s.) | enrichment | depletion | **S0_only_exploratory** | medium |

**DQB1*02:01 는 2026-05-06 zero-cell QC 감사 결과 typing-resolution artifact 가능성이 높아 일시 suspend 상태**입니다 (`2026_05_06_paper2_dqb1_zero_cell_qc_report.md` 참조; DQB1 callability 72.2%, *02:01 mass가 *02:02 + 40개 rare *02:xxx subtype으로 분산). orthogonal NGS-typer validation 후 재평가.

Reference TSV: `project/results/p2_pillar1_forest_v2/paper2_hla_claim_grade_v2.tsv`.

---

## 7. 가장 강하게 말할 수 있는 결론 (허용)

> "한국인 PTC 874명 풀에서 carrier-vs-allele 메트릭 불일치로 6-allele HLA exploratory forest가 만들어졌고, 엄격한 metric 정합화 후 다음과 같은 후보 우선순위 지도가 남는다:
> (1) C*01:02와 DQB1*02:01의 PTC 측 **depletion** 신호가 4개 시나리오 + 직접 allele 비교에서 모두 q < 0.05로 살아남는다.
> (2) 1차 forest에서 가장 두드러졌던 enrichment 후보(A*02:07, B*46:01, DPB1*05:01, DRB1*07:01)는 strict harmonization에서 살아남지 않는다.
> (3) 따라서 6-allele forest는 면역-유전학적 background 신호의 **탐색적 후보 우선순위 지도**이며, **matched-control NGS validation 없이는 case-control association으로 보고할 수 없다.**"

이 한 문장이 본 자료에서 허용된 최강 표현입니다.

---

## 8. 절대 사용 금지 표현 (forbidden claims)

- "Korean PTC HLA association이 확인되었다" / "case-control association 결과"
- "DPB1*05:01 / A*02:07 / B*46:01이 한국인 PTC 위험 allele이다"
- "HLA에 의한 인과(causal) 감수성"
- "임상 위험 예측" / "환자 선별(patient selection)"
- "carrier 빈도와 allele 빈도를 동등하게 비교 가능"
- "8-gene index와 통합 위험 모델로 사용 가능"
- "AITD/HT/Graves' 메커니즘과 직접 연결됨" *(이건 다른 paper의 영역; 여기서 주장하지 않음)*

---

## 9. 다음 검증 단계 (validation roadmap)

`F19_validation_roadmap.png` 참조.

### 9.1 필수 (Required)
1. **Matched-control NGS HLA typing** — 독립 한국인 healthy population, PTC와 동일한 typing pipeline (arcasHLA 또는 동등) + 동일 해상도(4-digit) + 동일 callable threshold.
2. **Allele-level dosage harmonization** — 양측 모두 allele copy 수 / 2n_callable 단위로 통일.
3. **Multi-cohort replication** — 최소 2개 독립 한국인 PTC cohort에서 동일 typing pipeline으로 검증.
4. **Pre-registration** — 가설(어느 allele이 어느 방향인지), 분석 plan, multiple-test 보정 규칙을 데이터 카운팅 전에 등록.

### 9.2 졸업(graduation) 기준
- S0/S1/S2/S3 4시나리오 모두에서 동일 방향 + matched-control 비교에서 q<0.05.
- 최소 2개 독립 cohort에서 재현.
- 효과크기 > 사전 정의된 minimal important difference (MID).
- 95% CI가 1.0을 포함하지 않음.

졸업 후에야 "Korean PTC HLA association" 표현 가능. 그 전까지는 본 자료 § 7 표현 이상으로 강하게 말하지 않습니다.

### 9.3 보조 분석 (optional, 강한 후보를 위한)
- Kim 2014 Korean reference panel을 baseline으로 한 carrier-matched validation forest (이미 partial하게 진행 — F17_kim2014_*).
- Baek 2021 PLOS NGS-baseline에 대한 allele-level forest (F19_allele_vs_baek2021_ngs_forest 참조).
- 8-gene readout과의 conditional independence 검정 — 본 패킷에서는 *수행하지 않았으며* claim하지 않음.

---

## 10. Paper 9 perturbation 연결 (future work only)

C*01:02 depletion이 matched-control validation에서 살아남는다면, 해당 allele가 결손된 PTC sample의 immune evasion / antigen presentation 결손 phenotype이 Paper 9의 합성치사(synthetic lethality) 후보 (예: GLS, metabolic dependencies)와 교차할 가능성이 있습니다. **단, 본 패킷에서는 가설-수준의 향후 작업으로만 언급하며 어떤 quantitative 연결도 주장하지 않습니다.**

---

## 11. 첨부 자료

- 그림: F13, F16, F17, F18, F19_validation_roadmap.png (모두 `project/papers_hub_2026_05_04/assets/paper2_hla/`)
- 표: `project/results/p2_pillar1_forest_v2/`
  - `paper2_6allele_or_fisher_primary.tsv` — S0 1차
  - `paper2_hla_metric_sensitivity_strengthened.tsv` — S0+S1+S2+S3 종합 + 안정성
  - `paper2_hla_allele_vs_allele_direct.tsv` — S3 직접 비교
  - `paper2_hla_bh_fdr_by_scenario.tsv` — 시나리오별 BH-q
  - `paper2_hla_subcohort_heterogeneity.tsv` — subcohort χ²
  - `paper2_hla_claim_grade.tsv` — per-allele 등급
  - `paper2_matched_control_validation_gate.tsv` — gate 기준
- 웹 페이지: http://40.82.129.113/paper2-hla/

---

**End of professor packet. Claim level: exploratory candidate-prioritization map only. No final association reported.**
