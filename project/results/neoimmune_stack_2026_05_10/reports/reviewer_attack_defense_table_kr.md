# Reviewer Attack Defense Table KR

| Expected attack | Answer | Figure/table support | Claim boundary |
|---|---|---|---|
| 이건 binding predictor 조합 아닌가? | clean track과 production track을 분리했다. clean track은 public predictor score를 feature로 쓰지 않는다. | `model_registry.tsv`, `leaderboard_clean_track_strict_no_overlap.tsv` | production stack은 실용 레이어, clean novelty는 local-only track에서만 주장 |
| leakage가 심한 public benchmark 아닌가? | 맞다. 그래서 exact peptide, pHLA, study, patient, HLA, public-tool training risk를 별도 audit하고 strict no-existing-overlap view를 만들었다. | `leakage_audit.md`, `leakage_summary.tsv` | random split 숫자는 smoke test로만 사용 |
| patient-level ranking이라면서 patient ID가 없는데? | 현재 public scaffold는 patient ID가 부족하므로 patient-level endpoint를 claim하지 않는다. 병원 데이터 수집이 다음 gate다. | `patient_topN_report.md`, `hospital_data_request_sheet.tsv` | primary endpoint는 병원 per-patient table 확보 후 |
| clinical vaccine efficacy를 예측하나? | 아니다. 후보 우선순위화 시스템이다. efficacy claim은 금지한다. | `claim_ladder.tsv`, `claim_ladder.png` | no clinical efficacy claim |
| presentation을 증명했나? | 아니다. NetMHC/MHCflurry/BigMHC는 predictor다. MS 없으면 presentation confirmed라고 말하지 않는다. | `algorithm_integration_matrix.tsv` | MS 필요 |
| immunogenicity를 증명했나? | T-cell assay label이 있는 경우에만 benchmark label로 쓴다. 새로운 후보 immunogenicity는 wet-lab 전까지 hypothesis다. | `failure_case_audit.md` | T-cell assay 필요 |
| quantum kernel은 진짜 살아남나? | 현 단계에서는 promising branch다. source-heldout/no-overlap에서 독립적으로 살아남아야 한다. | `leaderboard_clean_track_strict_no_overlap.tsv` | overclaim 금지 |
| Wave8/TCR이 왜 중요한가? | vaccine bottleneck은 binding만이 아니라 TCR-visible immunogenicity다. strict no-overlap clean view에서 TCR/self-sim branch가 상위권이다. | strict clean leaderboard | 다른 candidate set comparator와 직접 비교 금지 |
| LLM 쓰면 환각/black box 아닌가? | LLM은 점수 feature가 아니다. structured evidence를 읽어 contradiction/rationale/report만 생성한다. | `latest_llm_copilot_blueprint.md` | no LLM-as-predictor claim |
