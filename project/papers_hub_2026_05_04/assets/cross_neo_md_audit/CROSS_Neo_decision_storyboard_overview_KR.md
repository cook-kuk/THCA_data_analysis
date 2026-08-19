# CROSS-Neo Decision Storyboard Package

## 한 줄 결론

이 패키지는 CROSS-Neo 후보를 cheap DL, uncertainty, TCR branch, Baker/static structure, MD audit 순서로 줄이는 전체 흐름을 발표/논문용으로 설명한다.

## 대표 숫자

- 시작 후보: 649개 (positive 161, negative 488).
- HMTEVVRHC: 10 ns explicit-solvent MD 완료, current MD label `MD_VERY_STRONG` (score 0.822).
- GADGVGKSAL: 10 ns explicit-solvent MD 완료, current MD label `MD_MODERATE` (score 0.557), 추가 replicate/control batch 실행 중.
- 현재 default MD escalation: 2개 (TP 2, FP 0).
- 현재 default wetlab shortlist: 1개 (TP 1, FP 0).
- No-FP optimizer preset: 13 TP, 0 FP, precision 1.000, recall 0.081.

## 추천 Figure 흐름

1. `fig_md24_story_overview_flow`: 전체 후보 축소 흐름과 핵심 숫자.
2. `fig_md18_dl_first_candidate_funnel`: DL-first gate별 후보 감소.
3. `fig_md22_threshold_precision_recall_frontier`: threshold preset 탐색 결과.
4. `fig_md26_representative_candidate_board`: GADGVGKSAL/HMTEVVRHC 대표 후보 근거와 최신 MD evidence.
5. `fig_md25_data_usage_matrix`: 어떤 데이터가 어느 단계에 쓰였는지.

## 데이터 사용처

| artifact                                |   rows | used_for                                                            | pipeline_stage          | claim_boundary                                                       |
|:----------------------------------------|-------:|:--------------------------------------------------------------------|:------------------------|:---------------------------------------------------------------------|
| TCR wetlab candidate table              |    649 | candidate universe, labels, pMHC/TCR branch scores                  | input_registry          | labels are dataset labels, not direct clinical immunogenicity proof  |
| CROSS-Neo prediction ensemble           |  81081 | main DL ensemble aggregation and rank evidence                      | DL_first_screen         | model score is prioritization evidence, not a calibrated probability |
| Bayesian/dropout uncertainty tables     |    649 | remove unstable or low-posterior candidates before expensive layers | uncertainty_gate        | uncertainty summaries are internal triage metrics                    |
| DL-first funnel with labels             |     13 | show how positive/negative composition changes across gates         | dashboard_and_reporting | stage enrichment is retrospective on available labels                |
| Threshold optimizer presets             |      6 | label-aware operating point selection for slider presets            | threshold_optimization  | optimized on current labels; external validation required            |
| TCR/structure/MD escalation candidates  |      2 | identify candidates to send into structure/MD or wetlab controls    | escalation              | escalation means higher priority, not confirmed positive             |
| Baker/Rosetta fallback interface scores |      7 | cheap structure sanity check before Rosetta/MD                      | cheap_structural_filter | static contacts are not Rosetta energies or immunogenicity proof     |
| OpenMM MD status and evidence           |      7 | late structural audit of peptide-MHC/TCR-pMHC stability             | MD_audit                | MD supports structural plausibility, not immune recognition proof    |

## Figure/Table 설명

| id                | title                                        | main_message                                                                                                                              | claim_boundary                                                           |
|:------------------|:---------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------|:-------------------------------------------------------------------------|
| Fig. A / fig_md24 | End-to-end DL-first recognition-aware funnel | The pipeline turns a broad labeled candidate set into a tiny, auditable wetlab/MD shortlist.                                              | Enrichment and prioritization, not proof of immunogenicity.              |
| Fig. B / fig_md25 | Data usage and evidence lineage matrix       | Every claim has a traceable source table and a stated boundary.                                                                           | Lineage map does not validate model performance by itself.               |
| Fig. C / fig_md26 | Representative candidate evidence board      | HMTEVVRHC now has the strongest completed MD signal, while GADGVGKSAL remains the more balanced DL+TCR+MD candidate with controls queued. | Candidate evidence does not establish clinical utility.                  |
| Table 1           | Stage-by-stage label composition             | Positive fraction rises from 24.8% in the full set to 100% in MD escalation/wetlab shortlist under current settings.                      | Current-label enrichment only; external holdout needed for final claims. |
| Table 2           | Threshold operating presets                  | A no-false-positive preset found 13 labeled positives with 0 false positives in the current table.                                        | Preset is optimized on current labels and must be validated externally.  |
| Table 3           | Data usage manifest                          | The pipeline is auditable from input label table to final UI.                                                                             | Documentation table, not performance evidence.                           |

## Claim Boundary

- 이 결과는 후보 우선순위/실험 설계 근거다.
- MD와 구조 점수는 면역원성 증명이 아니다.
- threshold preset은 현재 label table에서 최적화된 값이므로 외부 validation 전까지 clinical threshold로 주장하면 안 된다.
