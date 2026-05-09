# Cancer Vaccine Algorithm Decision Context - 2026-05-09

작성 목적: 지금 암 백신/neoantigen 알고리즘 개발 방향을 결정하기 위한 내부 판단 MD.

범위: local repo snapshot + 주요 임상 anchor의 빠른 공개 출처 확인. 이 문서는 임상 조언이 아니며, 환자 치료 선택 문서가 아니다. 현재 모델은 후보 prioritization, triage, abstention용 연구 도구로만 다룬다.

## 0. One-line decision

**현재 알고리즘 개발 verdict: HOLD for claim, CONTINUE for CLEAN-Neo v0 build.**

즉, "우리가 새 SOTA neoantigen predictor를 만들었다"는 claim은 보류한다. 대신 **contamination-controlled, leakage-aware, evidence-calibrated prioritization framework**로 재설계한다. 지금 가장 가치 있는 산출물은 더 높은 AUROC가 아니라 **어떤 후보를 신뢰하지 말아야 하는지 알려주는 abstention/uncertainty와 reviewer-safe benchmark contract**다.

## 1. 현재 전체 상황 요약

| 축 | 현재 상태 | 판단 |
|---|---|---|
| 임상 전략 | PAAD는 resected/MRD low-burden vaccine window가 가장 강함. THCA는 routine PTC vaccine-first가 약하고, ATC/PDTC/progressive RR-DTC patient-selection-first가 맞음. | 알고리즘은 질병별 patient gate를 먼저 통과시켜야 함. |
| Local class-I predictor | `Structure_LR`가 가장 깨끗한 local baseline. `MHCflurry`는 no-overlap 성능은 좋지만 public pretrained training-overlap 미해결. | `Structure_LR`를 internal anchor로 둔다. |
| Frozen ESM2 Bayesian | ITSNdb no-overlap AUROC 0.411로 below chance. Cross-source LOSO와 entropy OOD signal은 개선. | ranking 모델로 promote 금지. uncertainty branch로만 보존. |
| Wave 8 TCR/self similarity | no-overlap AUROC 0.735로 apparent leader지만 exact TCR-reference hit 제거 후 Structure_LR보다 약함. | retest candidate. 현재 claim 금지. |
| Quantum/QK branch | strict internal repeated CV에서 강한 파일럿이 있음. CROSS-Neo v1에서 fold-safe fusion 개선도 있음. 하지만 harm negatives가 많고 source-heldout collapse가 남음. | bounded fallback/fusion branch로만 사용. quantum advantage claim 금지. |
| CROSS-Neo v1 lockdown | primary locked split 4/4에서 anchor 대비 AUPRC/top-k 개선. 하지만 NEPdb/TESLA source-heldout top-k collapse와 public overlap 미해결. | `HOLD`. 내부 locked-split framework로만 표현. |
| Class II | 별도 comparable MHC-II peptide immunogenicity benchmark가 없음. | Class I과 pooling 금지. Class II는 별도 module로 설계. |
| Korean translation | HLA-A*24:02/A*11:01/A*02:01 등 allele-specific gap이 중요. ESM2 Bayesian은 A*11:01에서 소폭 개선(+0.020)만 있음. | Korean cohort용 allele coverage/abstention이 우선. |

## 2. 지금 선택해야 할 개발 방향

### Recommended path: CLEAN-Neo v0

`CLEAN-Neo`를 새 architecture 이름으로 잠정 고정한다.

Core design:

1. **Strict benchmark first**
   - Exact peptide-HLA holdout.
   - Near peptide similarity holdout.
   - Source-protein/window holdout.
   - Study/patient holdout.
   - HLA allele/supertype holdout.
   - Time/public-corpus overlap audit.

2. **Model branch는 간결하게 시작**
   - Mutant-WT delta sequence and biophys features.
   - HLA pseudo-sequence features.
   - Structure proxy features as fallback/uncertainty, not standalone proof.
   - Leakage-safe train-only retrieval evidence.
   - Calibration/abstention head.

3. **Public tool scores는 feature로 금지**
   - MHCflurry, NetMHCpan, BigMHC, PRIME, MixMHCpred, NetMHCstabpan 등은 comparator 또는 caveated upper-bound.
   - 이 도구들의 training corpus overlap audit 전에는 clean external baseline이라고 쓰지 않는다.

4. **Deep model은 v0 이후**
   - ESM2/ESM2-650M fine-tune, SaProt/ProSST, geometry encoder 등은 strict v0가 baseline을 이긴 뒤 붙인다.
   - 지금처럼 frozen ESM2 head만 바꾸는 반복은 효율이 낮다.

### What to stop now

| Stop item | 이유 |
|---|---|
| Frozen ESM2 MLP head 튜닝 반복 | ITSNdb no-overlap below chance를 못 넘었고, architecture head만 바꾸면 leakage/in-domain 개선으로 흐를 위험이 큼. |
| "quantum advantage" framing | QK가 rescue도 하지만 harmed negatives가 훨씬 많음. reviewer-safe하지 않음. |
| Class I/Class II pooling score | benchmark와 biology가 다르다. MHC-II는 separate task. |
| Public pretrained tool을 clean comparator로 부르는 것 | official training corpus row-level audit가 미완료. |
| Clinical vaccine selection claim | HLA LOH/B2M/presentation, expression, clonality, immune state, safety gate 없이는 임상 선택 주장 불가. |

## 3. Evidence snapshot - clinical/product strategy

### PAAD

| segment | best current clinical anchor | algorithm implication |
|---|---|---|
| Resected/MRD low-burden personalized | autogene cevumeran + atezolizumab + mFOLFIRINOX. Nature 2025 long follow-up reports vaccine-induced T-cell responders had median RFS not reached versus 13.4 months in non-responders, HR 0.14, P=0.007. | Algorithm should prioritize manufacturable class I/II neoantigen set, low-burden timing, and immune response monitoring. |
| MRD-positive KRAS-public | ELI-002 2P/7P mKRAS amphiphile vaccine. Nature Medicine 2025 final phase 1 reports above-threshold mKRAS T-cell response associated with RFS not reached versus 3.02 months, HR 0.12, P=0.0002. | Shared KRAS candidates need driver mutation, HLA match, ctDNA/MRD context, and immunologic pharmacodynamic monitor. |
| Resected KRAS broad | mKRAS long-peptide vaccine + nivolumab/ipilimumab phase I. | Useful immunogenicity anchor but toxicity and small-n prevent strong deployment claim. |

PAAD algorithm gate:

`resected/MRD low-burden` -> `HLA/presentation intact` -> `KRAS-public or private neoantigen manufacturable` -> `expression/clonality high` -> `class I + class II balance` -> `uncertainty not high`.

### THCA

| segment | current backbone | algorithm implication |
|---|---|---|
| BRAF V600E ATC | FDA-approved dabrafenib + trametinib; FDA page reports ORR 61% in evaluable ATC patients. | Vaccine is add-on/trial concept only. Fast targeted backbone first. |
| Progressive RR-DTC/PDTC | lenvatinib + pembrolizumab phase 2 signal; reported PR 65.5% and median PFS 26.8 months in lenvatinib-naive cohort. | TKI/VEGF + PD-1 backbone is more defensible than vaccine alone. |
| BRAF WT ATC | NEO-COMBAT XL type trial logic: multi-kinase/VEGF-axis inhibitor + PD-1. | Algorithm should find antigen-positive, presentation-intact, immune-active patients for trial-style add-on, not routine PTC. |
| HT/TLS immune-active high-risk PTC | local HLA/TLS/IFN biology supports plausibility only. | Low-risk PTC no-go. High-risk recurrent/progressive only, with autoimmune safety gate. |

THCA algorithm gate:

`ATC/PDTC/progressive RR-DTC/high-risk recurrence` -> `driver/backbone defined` -> `HLA/presentation intact` -> `immune-active or combination backbone available` -> `allele-specific Korean HLA coverage` -> `autoimmunity risk manageable`.

## 4. Evidence snapshot - local algorithm results

### 4.1 Core local benchmark

Primary source: `project/results/p_neo_bayesian_2026_05_09/curation_2026_05_09/CLASS1_CLASS2_ALL_METHOD_COMPARISON.md`

| method | role | strict/no-overlap result | interpretation |
|---|---|---:|---|
| `Wave8_TCR_SelfSim_full` | candidate feature model | AUROC 0.735, AUPRC 0.598, n=106/pos=33 | Apparent leader but exact-reference sensitive. |
| `MHCflurry` | public pretrained comparator | AUROC 0.668, AUPRC 0.553 | Useful upper-bound comparator; training overlap unresolved. |
| `Structure_LR` | primary local baseline | AUROC 0.653, AUPRC 0.489, n=103/pos=32 | Cleanest known-local-split comparator. |
| `BigMHC_IM` | public pretrained comparator | AUROC 0.626, AUPRC 0.380 | In-master inflated; public overlap unresolved. |
| `ESM2_Bayesian` | local deep Bayesian | AUROC 0.411, AUPRC 0.274 | Negative headline and leakage-failure baseline. |

Strict subset after removing TCR/self exact hits:

| method | AUROC | AUPRC | read |
|---|---:|---:|---|
| `Structure_LR` | 0.679 | 0.435 | strict-set leader |
| `MHCflurry` | 0.657 | 0.488 | caveated public comparator |
| ESMFold min peptide pLDDT | 0.644 | 0.313 | weak univariate structure proxy |
| `Wave8_TCR_SelfSim_full` | 0.632 | 0.403 | drops after exact-reference removal |
| ESMFold 3D LR CV | 0.569 | 0.279 | not enough as replacement |

Decision: `Structure_LR` remains the honest local anchor until CLEAN-Neo v0 beats it on strict external/top-k metrics.

### 4.2 Frozen ESM2 Bayesian PoC

Primary source: `project/results/p_neo_bayesian_2026_05_09/BAYESIAN_NEO_REPORT.md`

| benchmark | RF baseline | Bayesian ESM2 | delta | read |
|---|---:|---:|---:|---|
| ITSNdb no-overlap | 0.431 | 0.411 | -0.020 | no rescue; both below chance |
| ITSNdb combined | 0.734 | 0.784 | +0.050 | inflated by overlap; not headline |
| VenusVaccine test top10_mean | 0.779 | 0.738 | -0.041 | protein-level plateau not solved |
| Cross-source LOSO mean | 0.484 | 0.583 | +0.099 | real robustness gain |
| HLA-A*11:01 LOSO | 0.576 | 0.596 | +0.020 | small Korean-relevant improvement |
| entropy OOD-AUROC | n/r | 0.638 | n/a | useful abstention signal |

Decision: keep Bayesian/ESM2 as uncertainty/OOD branch, not as primary ranking branch.

### 4.3 CROSS-Neo v0/v1/v1-lockdown

Primary sources:

- `project/results/cross_neo_v0/CROSS_Neo_v0_decision_report.md`
- `project/results/cross_neo_v1/CROSS_Neo_v1_decision_report.md`
- `project/results/cross_neo_v1_lockdown/CROSS_Neo_v1_lockdown_decision_report.md`
- `project/results/cross_neo_v1_lockdown/CROSS_Neo_v1_reviewer_audit_pack.md`

v0 best eligible internal model:

| split | feature/model | n/pos | AUPRC | AUROC | top10 |
|---|---|---:|---:|---:|---:|
| exact peptide-HLA holdout | `C_counterfactual` / `rf_secondary` | 89/21 | 0.525 | 0.715 | 0.600 |

v1 lockdown primary delta versus anchor:

| split | selected method | AUPRC | anchor AUPRC | top10 | anchor top10 | read |
|---|---|---:|---:|---:|---:|---|
| exact peptide-HLA holdout | `prespecified_rf_qk_no_anchor_w0.5` | 0.587 | 0.525 | 0.700 | 0.600 | pass |
| HLA stratified group 5-fold | `nested_rf_qk_quantum_only_train_selected` | 0.555 | 0.505 | 0.600 | 0.600 | pass by AUPRC |
| HLA supertype heldout | `prespecified_lr_qk_quantum_only_w0.5` | 0.545 | 0.478 | 0.600 | 0.500 | pass |
| near peptide cluster holdout | `rule_gate_rf_qk_fallback_train_selected` | 0.519 | 0.442 | 0.700 | 0.500 | pass |

But source-heldout remains blocking:

| heldout source | n/pos | AUPRC | AUROC | top10 | root issue |
|---|---:|---:|---:|---:|---|
| CEDAR | 913/851 | 0.922 | 0.439 | 1.0 | high prevalence, source/HLA/miscalibration |
| NEPdb | 886/354 | 0.396 | 0.479 | 0.0 | source/label shift, missing WT/source window |
| TESLA_mmc4 | 610/37 | 0.055 | 0.465 | 0.0 | low prevalence, source shift |
| TESLA_mmc7 validation | 319/6 | 0.017 | 0.294 | 0.0 | very low prevalence, source shift |

QK boundary:

| metric | value |
|---|---:|
| Total audited rescued positives | 119 |
| Total audited harmed negatives | 419 |
| Allowed interpretation | bounded fallback/fusion component |
| Forbidden interpretation | standalone quantum advantage |

Decision: CROSS-Neo v1 is useful as a locked-split internal framework, but final disposition remains `HOLD`.

## 5. Claim boundary

### Allowed now

- "Internal locked-split neoantigen prioritization framework."
- "Mutant-WT counterfactual encoding plus fold-safe fallback/fusion improves internal/HLA/near-cluster ranking versus clean anchor on locked splits."
- "The system is source-shift limited and therefore uses uncertainty/abstention."
- "`Structure_LR` is the clean local baseline; public tools are caveated comparators."
- "PAAD vaccine algorithm should prioritize resected/MRD low-burden windows; THCA requires patient-selection-first."

### Not allowed now

- External validation.
- Clinical vaccine selection claim.
- Quantum advantage.
- SOTA predictor claim.
- Clean public comparator claim before official training-corpus overlap audit.
- THCA routine PTC vaccine-first claim.
- Class I/Class II pooled peptide predictor claim.
- Hiding NEPdb/TESLA source-heldout collapse.

## 6. Immediate development plan

### E0 - Freeze and label current evidence

Output:

- `CLEAN_NEO_E0_LOCKED_BASELINE.md`
- One table with `Structure_LR`, `MHCflurry`, `BigMHC`, `ESM2_Bayesian`, `Wave8`, `CROSS-Neo v1`.

Rule:

- Every row gets `clean`, `caveated_public`, `reference_sensitive`, `internal_only`, or `negative_baseline`.

### E1 - Public corpus overlap audit

Priority targets:

| source | action |
|---|---|
| MHCflurry | download/parse official training corpus or package data |
| NetMHCpan | document available BA/EL training references; row-level if possible |
| BigMHC | download release training data and row-audit peptide/HLA |
| PRIME/MixMHCpred/NetMHCstabpan | parse public training-like corpus if available |
| IEDB/CEDAR/NEPdb/TESLA | keep local row-level exact/near overlap flags |

Success criterion:

- Every benchmark row has exact peptide, exact peptide-HLA, near peptide cluster, source protein/window, study, and public-corpus overlap flags where metadata allow.

### E2 - CLEAN-Neo v0 strict class-I

Input:

- `project/results/p_neo_bayesian_2026_05_09/clean_neo_v0_2026_05_09/clean_neo_v0_classI_strict_seed.tsv`
- Current strict set: n=89, positives=21.

Model:

- Start with interpretable tabular model.
- Mutant-WT sequence deltas.
- HLA pseudo features.
- Structure proxy.
- Train-only retrieval evidence.
- Calibration and conformal/coverage abstention.

Go criterion:

- Beats `Structure_LR` on strict AUPRC/top-k, not only AUROC.
- Does not collapse on near peptide/HLA heldout.
- Does not worsen source-heldout top-k versus anchor.

### E3 - Patient-level triage layer

For PAAD/THCA use, peptide score alone is insufficient. Add patient-level gates:

| layer | required fields |
|---|---|
| disease timing | resected/MRD versus metastatic/progressive |
| mutation/antigen | KRAS/shared driver or private expressed neoantigens |
| HLA/presentation | HLA typing, HLA LOH, B2M, HLA-A/B/C, antigen-processing genes |
| tumor context | expression, clonality/VAF, RNA support, TMB, immune state |
| safety | ECOG, immunosuppression, autoimmune history, thyroid autoimmunity for THCA |
| model reliability | OOD distance, calibration, abstention flag |

Output:

- `patient_candidate_score = disease_gate * presentation_gate * antigen_gate * immune_context_gate * model_confidence`

No patient should receive a high rank if any hard gate fails.

### E4 - Class II separate module

Do not blend into class-I score. Build separate MHC-II benchmark:

- DR/DQ/DP allele normalization.
- 13-25mer handling.
- Weak-label or MIL framing.
- CD4/helper evidence and uncertainty reported separately.

## 7. Decision matrix for Seungho

| choice | upside | risk | decision |
|---|---|---|---|
| Keep tuning frozen ESM2 Bayesian | quick GPU iteration; uncertainty branch exists | no-overlap below chance; likely wasted cycles | **No** |
| Promote CROSS-Neo v1 | locked split gains; reviewer audit pack exists | source-heldout collapse and public-overlap block | **Hold** |
| Build CLEAN-Neo v0 strict tabular | reviewer-safe; solves core leakage problem | less glamorous than deep model | **Yes** |
| Add ESM2-650M/fine-tuning now | could improve representation | before strict data contract, may just inflate in-master | **Later** |
| Push QK as main novelty | distinct method story | harmed negatives and no source survival | **No** |
| Use QK as fallback/fusion | selected internal gains | needs harm audit | **Yes, bounded** |
| Build PAAD/THCA clinical triage score | matches translational product use | must avoid clinical overclaim | **Yes, non-clinical research triage** |
| Build Class II now | necessary for vaccines | no comparable benchmark yet | **Separate workstream** |

## 8. Files that currently matter

### Strategy and triage

- `project/results/p_cancer_vaccine_strategy_2026_05_09/SUMMARY.md`
- `project/results/p_cancer_vaccine_strategy_2026_05_09/evidence_matrix.tsv`
- `project/results/p_cancer_vaccine_strategy_2026_05_09/patient_triage_matrix.tsv`
- `project/results/p_cancer_vaccine_strategy_2026_05_09/combination_priority_matrix.tsv`
- `project/papers_hub_2026_05_04/cancer_vaccine_full_dossier.html`
- `project/papers_hub_2026_05_04/paad_thca_vaccine_combo_triage_2026_05_09.html`

### Algorithm blueprint and wave curation

- `project/results/p_neo_bayesian_2026_05_09/NEOANTIGEN_METHOD_BLUEPRINT_2026_05_09.md`
- `project/results/p_neo_bayesian_2026_05_09/BAYESIAN_NEO_REPORT.md`
- `project/results/p_neo_bayesian_2026_05_09/curation_2026_05_09/WAVE_CURATION_SUMMARY.md`
- `project/results/p_neo_bayesian_2026_05_09/curation_2026_05_09/CLASS1_CLASS2_ALL_METHOD_COMPARISON.md`
- `project/results/p_neo_bayesian_2026_05_09/clean_neo_v0_2026_05_09/README.md`
- `project/results/p_neo_bayesian_2026_05_09/clean_neo_v0_2026_05_09/FALLBACK_ALGORITHM_SWEEP_SUMMARY.md`

### CROSS-Neo

- `project/results/cross_neo_v0/CROSS_Neo_v0_decision_report.md`
- `project/results/cross_neo_v1/CROSS_Neo_v1_decision_report.md`
- `project/results/cross_neo_v1_lockdown/CROSS_Neo_v1_lockdown_decision_report.md`
- `project/results/cross_neo_v1_lockdown/CROSS_Neo_v1_reviewer_audit_pack.md`
- `project/results/cross_neo_v1_lockdown/CROSS_Neo_v1_method_card.md`
- `project/results/cross_neo_v1_lockdown/CROSS_Neo_v1_claim_boundary.md`

### Scripts

- `scripts/build_cross_neo_v0_master_table.py`
- `scripts/evaluate_cross_neo_v0_contract.py`
- `scripts/train_cross_neo_v0.py`
- `scripts/evaluate_cross_neo_v1_contract.py`
- `scripts/evaluate_cross_neo_v1_lockdown.py`
- `scripts/write_cross_neo_v1_lockdown_report.py`
- `project/results/p_neo_bayesian_2026_05_09/build_clean_neo_seed.py`
- `project/results/p_neo_bayesian_2026_05_09/build_method_training_audit.py`

## 9. External anchors checked for this memo

- Nature 2025 PDAC RNA neoantigen vaccine long follow-up: https://www.nature.com/articles/s41586-024-08508-4
- Nature Medicine 2025 AMPLIFY-201 ELI-002 2P final phase 1: https://www.nature.com/articles/s41591-025-03876-4
- Nature Communications 2026 mKRAS vaccine + nivolumab/ipilimumab phase I: https://www.nature.com/articles/s41467-026-68324-4
- ClinicalTrials.gov NCT05968326 autogene cevumeran phase II: https://clinicaltrials.gov/study/NCT05968326
- ClinicalTrials.gov NCT05726864 ELI-002 7P: https://clinicaltrials.gov/study/NCT05726864
- FDA dabrafenib + trametinib ATC approval: https://www.fda.gov/drugs/resources-information-approved-drugs/fda-approves-dabrafenib-plus-trametinib-anaplastic-thyroid-cancer-braf-v600e-mutation
- Lenvatinib + pembrolizumab RR-DTC phase 2: https://pmc.ncbi.nlm.nih.gov/articles/PMC11883846/
- NCI NEO-COMBAT XL trial page: https://www.cancer.gov/clinicaltrials/NCI-2025-01421

## 10. Final operating rule

For the next coding session:

1. **Do not chase a bigger model first.**
2. **Build the strict overlap-audited data contract first.**
3. **Make `Structure_LR` the local anchor.**
4. **Treat public tools as caveated comparators.**
5. **Use Bayesian/entropy and QK only as bounded support branches.**
6. **Promote only if strict AUPRC/top-k improves and source-heldout top-k does not collapse.**

Current disposition: **algorithm research continues, publication/clinical claim remains HOLD.**
