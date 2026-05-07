# Paper 2 HLA — 1-page summary (KR)
**2026-05-06 · Seungho Cook · 탐색적(exploratory) 후보 지도, association 결과 아님**

## 데이터
한국인 PTC 풀 n = 874명, 6개 후보 HLA allele(2/4-digit collapse).
3개 subcohort (K2, Lee2024, GSE286332_PTC=18명).
Baseline: In JW 2015 Ann Lab Med (n=613, 5 alleles) + AFND South Korea pool (n=680, DPB1).

## 핵심 caveat
1차 forest는 **PTC carrier(개인) vs baseline allele(2n)** 메트릭 불일치 비교.
HWE 환산: `carrier = 1−(1−AF)²`, `AF = 1−√(1−carrier)`.
4 시나리오(S0/S1/S2/S3) + 직접 allele dosage 비교 모두 보고.

## 한 문장 결론
> 6-allele HLA forest는 한국인 PTC의 면역-유전학적 background 신호에 대한 **탐색적 후보 우선순위 지도**이며, matched-control NGS validation 없이는 case-control association으로 보고할 수 없다.

## strict 결과 (S3 direct allele vs allele, BH-q within 6)

| Allele | 1차 OR (S0) | strict OR (S3) | strict q | 방향 안정성 | 등급 |
|---|---|---|---|---|---|
| A*02:07 | 2.43 | **1.26** | 0.43 | 안정 (작아짐) | exploratory_only |
| B*46:01 | 2.16 | **1.05** | 0.80 | 안정 (~null) | exploratory_only |
| C*01:02 | 1.47 | **0.68** | **7×10⁻⁴** | **flips → depletion** | **strict 통과 (depletion)** |
| DPB1*05:01 | 1.96 | **0.98** | 0.80 | flips, source-sensitive | source 의심 artifact |
| DQB1*02:01 | 0.026 (zero) | **0.018** | **5×10⁻⁸** | depletion 안정 (zero-cell) | **strict 통과 + zero-cell caveat** |
| DRB1*07:01 | 1.73 | **0.84** | 0.43 | flips | exploratory_only |

**6/6 → S3 strict에서 2/6만 q<0.05를 통과 (모두 depletion).**

## 주요 figure
- **F16** — 4-시나리오 forest (3/6 alleles direction-flip)
- **F17** — subcohort heatmap (Lee2024 dominance, GSE286332_PTC n=18 caution)
- **F18** — claim-boundary matrix
- **F19** — validation roadmap (matched-control NGS → 졸업 기준)

## 가장 강하게 허용 (max claim)
"Exploratory candidate-prioritization map suggesting Korean PTC may have immune-genetic background signals requiring matched-control validation."

## 절대 금지
case-control association · causal HLA susceptibility · clinical risk prediction · patient selection · carrier–allele 단위 동등시 · DPB1*05:01을 위험 allele로 보고

## 다음 단계
1. Matched-control NGS HLA typing (독립 한국인 healthy population, 동일 pipeline)
2. Allele-level dosage harmonization (양측 모두 2n)
3. Multi-cohort replication (≥ 2 independent Korean PTC)
4. Pre-registered hypothesis/analysis plan

## 웹
http://40.82.129.113/paper2-hla/
