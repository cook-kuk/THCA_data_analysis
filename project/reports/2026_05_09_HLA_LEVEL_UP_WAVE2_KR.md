# HLA Two-Paper Level-Up Wave 2

## 결론

두 논문은 계속 2개로 유지하는 것이 맞습니다. 이번 wave는 새 논문을 쪼개기보다, 각 논문 안에 reviewer-proof 정량 방어층을 추가했습니다.

## Paper 4 강화

- Allele triangulation: C*01:02는 GD-specificity lead로 triangulation score 46.3; DPB1*05:01은 AITD-broad class-II anchor로 score 28.5.
- Endpoint robustness: strict Graves, autoimmune hyperthyroid, broad hyperthyroid, HT/hypothyroid endpoint를 분리해 tag-SNP direction conflict를 phenotype sensitivity로 처리합니다.
- 이 레이어가 Paper 4를 단순 meta-analysis에서 population-resolved HLA architecture/resource 논문으로 올립니다.

## Paper 2 강화

- GSE286332 exact permutation: HLA-II d=3.65, exact p=0.0001; AP/TLS composite d=3.08, exact p=0.0002, AUC=0.96.
- TCGA HLA-II dominance: full model R2=0.83; top drop-one contributor는 Immune_proxy (R2 loss=0.050).
- 이 레이어가 Paper 2를 underpowered allele story가 아니라 HT-overlap antigen-presentation/TLS immune-state 논문으로 고정합니다.

## 새 산출물

- Tables: `project/results/hla_two_paper_synthesis_2026_05_09/level_up_wave2/tables`
- Figures: `project/results/hla_two_paper_synthesis_2026_05_09/level_up_wave2/figures`
- `T01_paper4_allele_triangulation_score.tsv`
- `T03_gse286332_exact_permutation_bootstrap.tsv`
- `T04_tcga_hla2_driver_dominance.tsv`
- `T05_wave2_level_up_ladder.tsv`

## Boundary

Paper 2의 HLA는 cancer cohort에서 expression module입니다. Allele association은 prospective validation roadmap에만 둡니다.
