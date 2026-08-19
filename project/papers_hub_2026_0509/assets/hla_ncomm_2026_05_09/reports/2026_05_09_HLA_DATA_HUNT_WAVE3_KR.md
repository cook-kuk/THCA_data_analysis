# HLA data hunt Wave3: 추가 외부데이터 검증

## 결론

2개 paper 전략은 그대로 가는 게 맞다. 다만 Wave3에서 Paper 2와 Paper 4의 약점이 많이 줄었다.

Paper 2는 이제 `GSE286332 -> GSE138198 -> GSE163203 -> GSE213647 -> GSE250521 -> TCGA-THCA -> single-cell/spatial/TLS`로 이어지는 expression validation ladder가 생겼다. 암 코호트에서 HLA allele을 주장하지 않고, `HLA-II antigen-presentation / TLS expression ecology`로 밀면 훨씬 단단하다.

Paper 4는 `Korean/non-Korean HLA genetics -> BBJ/KoGES population genetics -> Korean HLA reference -> GSE248205 AITD spatial tissue mechanism` 구조로 올라간다. GSE248205는 allele replication은 아니지만, HLA class-II/AP 유전학이 조직에서 어떤 antigen-presenting thyroid ecosystem으로 연결되는지 보여주는 강한 기전 보강이다.

## 새로 실제 분석 완료

| dataset | paper | n | 핵심 결과 | 논문 역할 |
|---|---|---:|---|---|
| GSE138198 bulk HT/PTC | Paper 2 | 36 samples | PTCwithHT vs PTCwithoutHT: T_IFNG d=1.23 p=0.037, HLA-II/AP d=1.05 p=0.082, AP/TLS d=1.00 p=0.090. HT vs TN: HLA-II/AP d=3.36 p=0.0036. | HT-background PTC bulk validation |
| GSE163203 scRNA PTC+HT | Paper 2 | 110,000 cells, 10 biological samples | PTCwithHT vs PTCwithoutHT: AP/TLS d=2.78 p=0.035, B/TLS d=2.59 p=0.035, HLA-II/AP d=2.57 p=0.053. | independent single-cell validation |
| GSE213647 Korean bulk RNA-seq | Paper 2 | 632 samples | PTC tumor vs PTC normal: HLA-II/AP d=1.39 FDR=4.03e-40, Myeloid/DC d=1.23 FDR=3.34e-36, HLA-I d=1.22 FDR=5.08e-36, AP/TLS d=0.89 FDR=3.85e-21. | large Korean expression generalization |
| GSE250521 thyroid Visium | Paper 2 | 16 Visium slides, 55,873 spots | Stage trend N/PTC/LPTC/ATC: CD74/MIF rho=0.76 q=0.0027, HLA-I rho=0.75 q=0.0027, Myeloid/DC rho=0.70 q=0.0055, HLA-II/AP rho=0.61 q=0.018. | spatial progression generalization |
| GSE248205 AITD spatial | Paper 4 mechanism, Paper 2 context only | 8 Visium samples, 16,985 spots | AITD vs control: CD74/MIF d=3.62, B/TLS d=2.57, AP/TLS d=2.26. HT spots above control p90: HLA-II/AP 0.964, B/TLS 0.997, AP/TLS 0.993. | AITD spatial antigen-presentation mechanism |
| GSE29315 AITD array | Paper 4 mechanism | 71 samples; HT 6 vs hyperplasia 8 primary | HT vs hyperplasia: HLA-II/AP d=5.45 p=0.0007, HLA-I d=4.40 p=0.0007, AP/TLS d=3.52 p=0.0007. | independent older-array corroboration |
| GSE6004 PTC invasion array | Paper 2 stress-test | 18 samples | Invasion vs normal did not support HLA/AP/TLS increase: HLA-II/AP d=-0.42, AP/TLS d=-0.59, all FDR=0.829. | negative/specificity stress-test only |
| GSE184362 + GSE191288 scRNA | Paper 2 low-tier supplement | 23 + 7 samples | HT label unavailable in GEO metadata/raw filenames. GSE191288 qualitative tumor vs one non-tumor direction: Myeloid/DC +1.36, HLA-I +1.02, HLA-II/AP +1.01. | expression-only generalization; not HT evidence |

## 새 산출물

| output | path |
|---|---|
| GSE138198 report | `project/reports/2026_05_09_GSE138198_HTPTC_BULK_VALIDATION_KR.md` |
| GSE138198 figures/tables | `project/results/hla_two_paper_synthesis_2026_05_09/gse138198_htptc_bulk_validation/` |
| GSE163203 report | `project/reports/2026_05_09_GSE163203_HTPTC_SCRNA_VALIDATION_KR.md` |
| GSE163203 figures/tables | `project/results/hla_two_paper_synthesis_2026_05_09/gse163203_htptc_scrna_validation/` |
| GSE213647 report | `project/reports/2026_05_09_GSE213647_KOREAN_BULK_VALIDATION_KR.md` |
| GSE213647 figures/tables | `project/results/hla_two_paper_synthesis_2026_05_09/gse213647_korean_bulk_validation/` |
| GSE250521 report | `project/reports/2026_05_09_GSE250521_THYROID_VISIUM_HLA_VALIDATION_KR.md` |
| GSE250521 figures/tables | `project/results/hla_two_paper_synthesis_2026_05_09/gse250521_thyroid_visium_hla_validation/` |
| GSE248205 report | `project/reports/2026_05_09_GSE248205_AITD_SPATIAL_VALIDATION_KR.md` |
| GSE248205 figures/tables | `project/results/hla_two_paper_synthesis_2026_05_09/gse248205_aitd_spatial_validation/` |
| GSE29315 report | `project/reports/2026_05_09_GSE29315_AITD_ARRAY_VALIDATION_KR.md` |
| GSE29315 figures/tables | `project/results/hla_two_paper_synthesis_2026_05_09/gse29315_aitd_array_validation/` |
| GSE6004 report | `project/reports/2026_05_09_GSE6004_PTC_INVASION_ARRAY_VALIDATION_KR.md` |
| GSE6004 figures/tables | `project/results/hla_two_paper_synthesis_2026_05_09/gse6004_ptc_invasion_array_validation/` |
| GSE184362/GSE191288 report | `project/reports/2026_05_09_GSE184362_GSE191288_SCRNA_GENERALIZATION_KR.md` |
| GSE184362/GSE191288 figures/tables | `project/results/hla_two_paper_synthesis_2026_05_09/gse184362_gse191288_scrna_generalization/` |

## 추가로 찾은 데이터 후보

| dataset/source | 즉시성 | paper | 판단 |
|---|---|---|---|
| GSE184362 + Frontiers Endocrinology 2024 HT/non-HT reanalysis | 낮음 | Paper 2 | Raw와 GEO metadata에는 HT label이 없어 직접 HT/non-HT 재분석 불가. scRNA expression-only stress-test로만 보관. |
| GSE213647 Korean thyroid bulk RNA-seq | 완료/높음 | Paper 2 | South Korea 632 fresh frozen tissues. HT-specific은 아니지만 Korean thyroid cancer HLA/AP expression generalization으로 강함. |
| GSE250521 thyroid cancer progression spatial/scRNA | 완료/높음 | Paper 2 | normal/PTC/LPTC/ATC Visium progression에서 HLA/AP/CD74/T-cell/myeloid stage trend가 양성. |
| GSE6004 PTC invasion array | 완료/낮음 | Paper 2 | invasive front vs normal에서 HLA/AP/TLS 증가가 없어 positive evidence가 아니라 negative stress-test. |
| GSE191288 PTC scRNA | 완료/낮음 | Paper 2 | non-tumor n=1이라 inferential use 불가. qualitative expression-only supplement. |

## 논문 급 판단

Paper 2는 이제 “HLA allele cancer association”을 버리고, `HT-overlap PTC HLA-II/AP-TLS immune ecology`로 쓰면 훨씬 강하다. 독립 bulk 2개, scRNA 1개, TCGA network, spatial/TLS가 묶여서 전문 oncology/immunology 저널 상위권 또는 Nature Communications 도전권으로 올릴 수 있다. 단, HT label과 clinical outcome을 더 정밀하게 붙이면 더 올라간다.

Paper 4는 genetics paper로 더 강하다. Korean AITD HLA, KoGES MHC/PheWAS, BBJ/East Asian GD genetics, AFND/Korean HLA reference, GSE248205 tissue mechanism까지 들어가면 HLA/AITD 분야에서는 상위권 narrative가 된다. Nature Communications급을 노리려면 자체 Korean adult GD case-control HLA-NGS 또는 KoGES/BBJ summary-level 재분석이 결정타다.

## 절대 금지선

1. Paper 2에서 `HLA allele is associated with thyroid cancer risk/prognosis/RAI/BRAF`라고 쓰지 않는다.
2. Paper 2의 HLA는 `HLA/AP expression module`이다.
3. GSE248205는 AITD spatial expression 기전자료이지 HLA genotype replication이 아니다.
4. Paper 4에서만 allele/genotype을 말한다.

## Sources

- GSE138198: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138198
- GSE163203: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163203
- GSE248205 / Nature Communications 2024: https://www.nature.com/articles/s41467-024-50192-5
- GSE184362 HT/non-HT reanalysis: https://www.frontiersin.org/journals/endocrinology/articles/10.3389/fendo.2024.1339473/full
- GSE213647: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213647
- GSE250521: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250521
- GSE29315: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29315
- GSE6004: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6004
