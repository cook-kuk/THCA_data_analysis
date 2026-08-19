# No-Yu-Data In Silico Paper Feasibility

## 결론

**가능하다.** 다만 지금 쓸 수 있는 논문은 `CTC-EMT transition discovery paper`가 아니라, **public-data 기반 CTC-EMT-NGS prior/framework paper**다.

가장 안전한 제목:

> Tissue-State and Genetic Prior Mapping for Post-Thyroidectomy CTC-EMT Monitoring in Papillary Thyroid Cancer

## 지금 주장 가능한 것

- Yu 2024는 PTC에서 serial CTC-EMT 측정 가능성을 이미 보여준 published anchor다.
- TCGA-THCA public pilot에서 BRAF-mutant 274례도 7개 tissue-state label로 갈라진다. 가장 큰 label도 33.2%라 mutation-only ambiguity가 66.8%다.
- Driver-only model은 RAI/HLA/immune/CD8 등 tissue-state axis를 일부만 설명한다.
- GSE250521 spatial data는 16 slides, 57144 spots에서 tissue-state organization을 지지한다. niche label same-neighbor coherence는 16/16 slides에서 z>2다.
- scRNA module attribution은 CTC marker panel의 cell-context rationale을 지원한다.
- thyroid bulk proteomics proxy는 dedifferentiation 방향에서 RAI protein loss와 myeloid/TGFB protein gain을 지지한다.

## 지금 절대 주장하면 안 되는 것

- 우리 데이터에서 thyroidectomy 후 CTC-EMT transition을 새로 증명했다.
- public pilot이 CTC shedding을 증명한다.
- CTC/cfDNA NGS가 early PTC에서 universal하게 성공한다.
- postoperative recurrence/residual-risk prediction을 이미 검증했다.
- 이 결과가 clinical diagnostic이다.

## 논문 급 판단

**Preprint / methods-framework / translational rationale paper는 지금 가능.**
병원 raw CTC table 없이 high-impact mechanism paper는 어렵다. Yu raw data가 들어오면 같은 framework가 바로 `serial CTC-EMT transition + tissue genetic anchor` 논문으로 업그레이드된다.

## Generated outputs

- `tables/in_silico_evidence_matrix.tsv`
- `tables/ctc_emt_marker_prior_table.tsv`
- `tables/publishability_decision_table.tsv`
- `tables/serial_ctc_ngs_dataset_gap_map.tsv`
- `figures/F01_no_yu_claim_ladder.png` ... `figures/F09_publishability_decision.png`
