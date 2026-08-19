# 기존 아이디어 x Kaggle 우승 비법 통합안

Date: 2026-05-10

## 한 줄 결론

기존 아이디어의 방향이 맞다. `BAR-Neo-BMA/X + stress gate + hard decoy + TCR expert + MD/structure audit + patient gate`를 코어로 두고, Kaggle 우승작의 비법은 새 모델을 마구 추가하는 용도가 아니라 다음 네 가지를 잠그는 용도로 쓰는 게 가장 좋다.

1. source/HLA/near-peptide에서 안 깨지는 validation contract
2. 모델 다양성 기반 BMA weight
3. 생물학적 validity DAG cap
4. label/public-overlap/metadata 리스크를 claim-safe abstention으로 보내는 장치

즉, 목표는 "neoantigen SOTA predictor"가 아니라 **reviewer-safe cancer vaccine prioritization engine**이다.

## 통합 아키텍처

### Layer 0. Intake / patient-context gate

기존 아이디어:

- PAAD/THCA patient-gated triage
- HLA LOH, HLA expression, antigen processing, expression, VAF, clonality, IFNG/TLS/cytolytic context
- missing context면 clinical-use claim 차단

Kaggle 보강:

- Open Problems single-cell perturbation식 context embedding
- one-hot + target encoding + missingness bit
- MoA식 non-scored meta-feature stacking

실행:

- patient metadata는 점수에 더하는 feature가 아니라 **상한 cap**으로 먼저 작동해야 한다.
- `missing_required_patient_context`는 0점이 아니라 `research_triage_only` abstention reason이다.
- 질병별 branch:
  - THCA: BRAF/RET/NTRK/TERT, TLS/IFNG, HLA presentation, RAI/redifferentiation context
  - PAAD: KRAS/MRD/adjuvant, low TMB, stromal/immune-cold penalty, vaccine+ICI/chemo sequencing caveat

### Layer 1. Clean pMHC / antigen evidence

기존 아이디어:

- Structure_LR honest anchor
- W7A/W7B/stacked internal candidates
- public comparator는 overlap audit 전까지 caveated comparator
- QK는 bounded fallback/fusion

Kaggle 보강:

- BELKA식 small transformer + task-relevant SSL pretraining
- CAFA식 hierarchy-aware postprocessing
- Ribonanza/OpenVaccine식 1D sequence + 2D pairwise structural bias

실행:

- P0는 새 거대 모델이 아니라 현재 score table에 `validity_dag_cap`을 붙이는 것.
- P1로 peptide/HLA SSL 모델:
  - masked peptide/HLA pseudo-sequence
  - descriptor prediction: anchor, length, charge, hydrophobicity, instability, cleavage/TAP-like descriptors
  - benchmark labels/public predictor score distillation 금지
- pMHC pairwise feature:
  - peptide position x HLA pocket contact
  - anchor-preserved decoy contact delta
  - MD contact occupancy summary

### Layer 2. Hard decoy / low-prevalence FP control

기존 아이디어:

- hard-decoy repair가 top-k triage를 살렸지만 standalone motif predictor는 위험
- no-FP preset이 high-impact decision package에서 강함

Kaggle 보강:

- PANDA식 label-noise audit
- online hard negative mining
- Ribonanza/OpenVaccine식 distribution-shift-aware private validation

실행:

- hard decoy v2는 5종으로 분리한다.
  - `anchor_preserved_decoy`: anchor는 보존, non-anchor scramble
  - `near_positive_decoy`: positive와 유사하지만 label/source/HLA가 다른 row
  - `hla_confuser_decoy`: 같은 peptide, 관련 supertype HLA
  - `tcr_shuffle_decoy`: 같은 peptide-HLA, shuffled TCR
  - `source_shift_decoy`: CEDAR-like positive vs NEPdb/TESLA-like negative
- no-FP threshold는 "논문 claim"이 아니라 "wetlab plate design preset"으로 둔다.
- low prevalence TESLA slice에서는 AUPRC보다 `top5/top10 precision`, FP count, enrichment를 우선 본다.

### Layer 3. TCR recognition diagnostic branch

기존 아이디어:

- pMTnet/TEPCAM external expert ensemble
- paired alpha/beta registry
- exact paired TCR evidence가 있는 후보만 TCR-supported로 올림
- TCR branch는 supplement/diagnostic, main pMHC claim 대체 금지

Kaggle 보강:

- Ribonanza two-track idea: sequence representation과 pairwise representation이 서로 업데이트
- MoA ensemble diversity: external experts 간 correlation/discordance를 weight에 반영
- CAFA conditional logic: presentation gate 아래에서만 TCR evidence가 작동

실행:

- TCR score는 pMHC presentation cap을 넘을 수 없다.
- pMTnet/TEPCAM 평균만 쓰지 말고 `expert_disagreement`, `exact_tcr_count`, `pathogen_only_dependency`, `cancer_context_evidence`를 함께 둔다.
- TCR-high / pMHC-low row는 promotion이 아니라 "diagnostic disagreement queue"로 간다.

### Layer 4. MD / Baker / Rosetta / structure audit

기존 아이디어:

- GADGVGKSAL/HLA-C*08:02: MD_MODERATE
- HMTEVVRHC/HLA-A*02:01: high priority but trajectory/contact completion boundary
- MD는 structural audit layer, immunogenicity proof 아님

Kaggle 보강:

- Ribonanza pairwise attention bias: MD contact occupancy를 model feature가 아니라 interpretability/diagnostic bias로 사용
- PANDA external validation discipline: reproduced/controlled validation 전에는 clinical claim 금지

실행:

- MD evidence label ladder:
  - `MD_INSUFFICIENT_RUNTIME`
  - `MD_PLAUSIBLE_PMHCI`
  - `MD_MODERATE_TCR_CONTACT`
  - `MD_STRONG_WITH_REPLICATE`
  - `MD_SPECIFIC_WITH_WT_DECOY`
- stronger claim에 필요한 최소 세트:
  - mutant, WT, scrambled peptide, same-HLA positive control
  - 3x10 ns each
  - top 후보만 50-100 ns escalation
- MD는 BMA score를 직접 boost하기보다 wetlab plate ordering과 reviewer explanation에 쓴다.

### Layer 5. BAR-Neo-BMA/X decision controller

기존 아이디어:

- stress-guarded method downweight
- BAR-Neo-X explanation and claim-safe reranking
- public weight/fallback weight/candidate confidence

Kaggle 보강:

- MoA winner: ensemble member selection by low prediction correlation
- Ribonanza: top blend useful하지만 single distilled model은 claim boundary 필요
- CAFA: postprocess outputs to obey biology

실행:

- BMA weight formula에 추가할 항목:
  - `mean_AUPRC`
  - `stress_min_floor`
  - `source_heldout_delta`
  - `hla_stress_delta`
  - `korean_hla_delta`
  - `low_prevalence_top10`
  - `pairwise_oof_corr_penalty`
  - `incremental_topk_gain`
  - `public_overlap_penalty`
  - `fallback_role_cap`
- BAR-Neo-X explanation은 reviewer page에서 "왜 올렸는가"보다 "왜 안 올렸는가"가 더 중요하다.

## 가장 강한 merged story

`CROSS-Neo`는 하나의 classifier가 아니라 cancer vaccine 후보를 줄이는 **multi-layer claim-safe triage system**이다.

- pMHC layer가 후보를 잡는다.
- hard-decoy/source/HLA stress layer가 가짜 generalization을 깎는다.
- TCR expert layer가 recognition 가능성을 진단한다.
- MD/structure layer가 top 후보의 구조적 말이 되는지를 본다.
- patient gate가 실제 환자/질병 context 없이는 claim을 막는다.
- BAR-Neo-BMA/X가 이 모든 증거를 합쳐 "promote, review, abstain"을 낸다.

이건 Kaggle식으로 보면 leaderboard 모델이 아니라 **private-LB shake-up 방지형 operating system**이다.

## P0 merged sprint

1. `BMA diversity/stress weight v2`
   - 기존 stress-guarded ranker에 MoA식 prediction-correlation penalty 추가.
   - 목표: W7A/W7B처럼 버티는 expert는 살리고, 같은 slice에서 같이 깨지는 expert는 중복 weight 제거.

2. `Validity DAG cap`
   - CAFA hierarchy trick을 vaccine biology에 맞춘다.
   - `claim_safe_score <= min(antigen_gate, presentation_gate, tcr_gate_or_unknown_cap, immune_context_gate, patient_metadata_cap, overlap_clean_cap)`.

3. `Label-noise / source-conflict audit`
   - PANDA식 OOF disagreement queue.
   - output: trusted positive, trusted negative, ambiguous assay, source conflict, manual review.

4. `Hard-decoy v2`
   - anchor-preserved, near-positive, HLA-confuser, TCR-shuffle, source-shift decoys.
   - standalone claim 금지, bounded auxiliary only.

5. `Wetlab plate locked preset`
   - high-impact no-FP 13 candidates를 wetlab plate logic으로 재정렬.
   - GADGVGKSAL/HMTEVVRHC는 MD/TCR supported tier.
   - 나머지는 pMHC-led tier with required WT/decoy controls.

## P1 merged sprint

1. `Peptide-HLA SSL small transformer`
   - BELKA식 simple encoder + MLM/descriptor pretrain.
   - 목표는 SOTA가 아니라 source/HLA stress에서 feature robustness.

2. `Pairwise contact-bias table`
   - Ribonanza식 2D pair signal을 DL이 아니라 먼저 tabular summary로 넣는다.
   - pMHC pocket contact, TCR-peptide contact, MD occupancy.

3. `Patient-context embedding/cap`
   - Open Problems single-cell식 context features.
   - THCA/PAAD branch를 score boost가 아니라 cap + explanation으로 사용.

4. `Length-heldout contract`
   - 9-mer overfit 차단.
   - 8/9/10/11/12+ bucket별 top-k, AUPRC, FP count.

## 명확한 no-go

- public predictor score를 clean feature로 쓰지 않는다.
- unresolved public training corpus overlap이면 clean comparator로 올리지 않는다.
- pseudo-label은 discovery triage 전용이다.
- MD/TCR/synthetic-decoy 결과로 immunogenicity claim을 하지 않는다.
- QK는 계속 bounded fallback/fusion이다.
- no-FP result는 wetlab prioritization preset이지 external validation claim이 아니다.

## 제일 좋은 논문/사업 포지션

논문용:

- "We present a leakage-aware, stress-tested neoantigen prioritization framework that integrates pMHC, decoy robustness, TCR-recognition evidence, and structure-aware audit layers, with explicit abstention under unresolved patient or public-overlap context."

사업용:

- "Fast intake triage -> claim-safe review queue -> wetlab plate ordering -> structure/TCR-supported top candidates."

둘 다 같은 시스템에서 나오되, claim level이 다르다.

