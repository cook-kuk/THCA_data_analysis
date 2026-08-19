# HLA Two-Paper Extra External Validation Pack (2026-05-09)

## Executive Upgrade

이번 추가 분석은 두 논문을 같은 HLA라는 단어로 억지 연결하지 않고, 서로 다른 증거축으로 수준을 올리는 구조입니다.
Paper 4는 비암 GD/AITD germline HLA allele 논문이고, Paper 2는 암 코호트의 HLA gene-expression module/immune-context 논문입니다.

## Paper 2: HT-PTC Immune/HLA Module

- GSE286332 HT-overlap 재검정: HLA-I d=2.34, p=0.0015; HLA-II d=3.65, p=0.0004, FDR=0.0012; TLS d=3.01, p=0.0011.
- TCGA-THCA module network: DM1-HLA-II rho=0.400, p=8.09e-21, n=505; IFNG-HLA-II rho=0.874, p=2.18e-159, n=505; TLS-HLA-II rho=0.760, p=1.26e-94, n=497; RAI-HLA-II rho=-0.479, p=2.91e-30, n=505.
- scRNA 외부 검증: Pu/Lu 두 scRNA 코호트의 HLA-II cell-type rank concordance rho=0.95, p=2.28e-05; top producers는 B_cell,DC,Myeloid.
- Spatial validation: GSE250521 TLS niche 안쪽 HLA-II가 바깥보다 높은 샘플 16/16, binomial p=1.53e-05; cancer subset mean d=0.83.
- ICI pan-cancer context: HLA-I response OR=1.35, p=0.0143; HLA-II OR=1.02, p=0.8394. 이 결과는 HLA-I translational context는 살리고, HLA-II는 thyroid/TLS biology 중심으로 제한하는 데 쓰는 것이 안전합니다.

## Paper 4: GD/AITD HLA Allele Genetics

- GWAS Catalog 외부 검증: GD genome-wide significant signals 중 MHC share=0.128 (36/282); hypothyroidism MHC share=0.305 (981/3213).
- FinnGen/PanUKBB tag-SNP 검증: DPB1*05:01 proxy는 strict Graves branch에서 positive signal p<1e-6 1개, p<1e-4 3개이고, HT/hypothyroid branch에서 genome-wide negative signal 5개입니다. Broad hyperthyroid endpoint에서는 negative genome-wide signal 1개가 있어 endpoint sensitivity panel로 분리해야 합니다. 직접 HLA imputation이 아니라 tag-SNP proxy로 표기해야 합니다.
- AFND population context: DPB1*05:01 South Korea weighted frequency=0.367, East-Asia pool=0.375; disease association이 아니라 background-frequency calibration으로만 사용합니다.

## Manuscript-Level Actions

- Paper 4 main Figure: Pan-Asian forest + GWAS MHC concentration + FinnGen/PanUKBB DPB1 branch divergence + AFND baseline context.
- Paper 2 main Figure: HT-overlap HLA/TLS effect + TCGA IFNG-HLA network + scRNA APC source + spatial TLS-HLA niche + ICI context panel.
- Paper 2 allele 후보(DRB1*04:05 등)는 n=9 vs 9에서 FDR-negative이므로 main claim이 아니라 power/future cohort box로 내려야 합니다.

## Output Index

- Tables: `project/results/hla_two_paper_synthesis_2026_05_09/extra_external_validation/tables`
- Figures: `project/results/hla_two_paper_synthesis_2026_05_09/extra_external_validation/figures`
- Validation ladder: `project/results/hla_two_paper_synthesis_2026_05_09/extra_external_validation/tables/T10_external_validation_ladder.tsv`

## Boundary

Cancer-cohort HLA는 gene-expression module입니다. HLA allele association은 Paper 4의 비암 GD/AITD 유전학에서만 주장합니다.
