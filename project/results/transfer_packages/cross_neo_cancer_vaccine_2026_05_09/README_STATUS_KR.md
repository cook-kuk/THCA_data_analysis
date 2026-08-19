# CROSS-Neo v0 + Cancer Vaccine 상황 정리 - 2026-05-09

## 결론

- CROSS-Neo v0는 구현/실행 산출물이 `project/results/cross_neo_v0/`에 있다.
- 현재 판정은 **HOLD**다. 내부/locked split 성능은 확인됐지만 외부검증이나 quantum advantage는 주장하지 않는다.
- best clean model은 `C_counterfactual / rf_secondary`이며 HLA-stratified split에서 AUPRC 0.388, AUROC 0.628, top10 precision 0.400이다.
- fixed late-fusion에서는 prespecified equal-weight 후보가 일부 split에서 더 강하지만, exploratory 가중치와 구분해서 해석해야 한다.
- Cancer vaccine 정리는 `project/papers_hub_2026_05_04/cancer_vaccine_full_dossier.html` 및 `project/results/p_cancer_vaccine_strategy_2026_05_09/`에 있다.

## 포함된 주요 경로

- `configs/cross_neo_v0.yaml`: locked config.
- `scripts/*cross_neo_v0*` 및 `scripts/run_cross_neo_v0_all.sh`: CROSS-Neo v0 파이프라인.
- `project/results/cross_neo_v0/`: master table, retrieval audit, embeddings, structure geometry, QK fallback, OOF predictions, metrics, figures, decision report.
- `project/results/p_cancer_vaccine_strategy_2026_05_09/`: PAAD/THCA vaccine strategy matrix and summary.
- `project/papers_hub_2026_05_04/cancer_vaccine_full_dossier.html`: cancer vaccine full dossier.
- `project/papers_hub_2026_05_04/cancer_vaccine_agent.html`: cancer vaccine session/page artifact.
- `project/papers_hub_2026_05_04/lumenix_cancer_vaccine_index.html`: vaccine index page.
- `project/results/p_neo_bayesian_2026_05_09/NEOANTIGEN_METHOD_BLUEPRINT_2026_05_09.md`: local neoantigen method blueprint.
- `project/results/p_neo_bayesian_2026_05_09/WEB_SEARCH_PROMPT_BETTER_NEOANTIGEN_ALGORITHMS_2026_05_09.md`: algorithm search prompt/context.

## CROSS-Neo headline numbers

- Strict class-I set: n = 89, positives = 21.
- Best internal/locked model: `C_counterfactual / rf_secondary` on `hla_stratified_group_5fold`.
- Best internal metrics: AUPRC 0.388, AUROC 0.628, top10 precision 0.400, enrichment@10 1.695.
- Repeated stratified 5x5 best: `C_counterfactual / rf_secondary`, AUPRC 0.353, AUROC 0.610, top10 precision 0.500.
- HLA supertype heldout best: `C_counterfactual / rf_secondary`, AUPRC 0.368, AUROC 0.691, top10 precision 0.400.
- Quantum fallback best table row: `qk_quantum_only_gamma1` on near-peptide holdout, AUPRC 0.460, AUROC 0.679. This is a fallback branch, not a quantum-advantage claim.
- Prespecified late fusion example: `counterfactual_rf_plus_qk_no_anchor_gamma1_w0.5` on HLA-stratified split, AUPRC 0.473, AUROC 0.692, top10 precision 0.500.
- Decision report: `project/results/cross_neo_v0/CROSS_Neo_v0_decision_report.md`.

## Leakage / interpretation guardrails

- Public benchmark predictor scores such as MHCflurry, NetMHCpan, BigMHC, PRIME, MixMHCpred, NetMHCstabpan are not used as CROSS-Neo model features.
- Retrieval features are train-fold only.
- Scaling, fitting, calibration, and QK gamma are fold-safe/fixed.
- Source-window, patient, and time holdouts are unavailable in the current strict bundle because the metadata are missing or single-valued.
- Public corpus overlap audit is partial/unresolved for public predictor training corpora.
- Do not claim external validation.
- Do not claim quantum advantage.

## Cancer vaccine bottom line

- PAAD: vaccine-first is most defensible in resected/MRD low-burden settings.
- THCA: vaccine-first is not routine-PTC defensible; use patient-selection-first logic for ATC/PDTC/progressive RAI-refractory DTC with intact presentation and rational backbone therapy.
- Neoantigen patient selection should require WES/RNA/HLA typing plus HLA LOH/B2M/presentation checks.
- Local neoantigen model should be used for ranking/abstention, not definitive immunogenicity proof.

## Live/local dossier locations on server

- Local repo HTML: `project/papers_hub_2026_05_04/cancer_vaccine_full_dossier.html`
- Live mirror: `/var/www/papers/papers_hub_2026_05_04/cancer_vaccine_full_dossier.html`
- Local CROSS-Neo results: `project/results/cross_neo_v0/`

