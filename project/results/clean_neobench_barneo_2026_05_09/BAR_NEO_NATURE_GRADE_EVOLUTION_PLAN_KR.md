# BAR-Neo Nature-Grade Evolution Plan KR

## 한 줄 결론

지금 Nature급으로 진화시키는 길은 더 큰 모델이 아니라 **claim ladder + validation gate + T1 execution packet**이다. 현재는 Nature Methods/Nature Biomedical Engineering 스타일의 benchmark/resource core까지 왔고, T1 3개는 pHLA assay-design 후보로 올라갔다. Translational/Nature Cancer급 claim은 metadata와 assay/lockbox가 붙어야 한다.

## 현재 핵심 숫자

| 항목 | 값 |
|---|---:|
| candidates | NA |
| methods | NA |
| split metrics | NA |
| reviewer kill pass | 3 |
| T1 pHLA assay ready | 3 |
| T1 translational ready | 0 |
| public clean comparators allowed | NA |

## Claim ladder

| claim_level                | claim                                                        | current_status          | evidence_now                                                  | required_to_upgrade                                                                                      | forbidden_overclaim                                        |
|:---------------------------|:-------------------------------------------------------------|:------------------------|:--------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------|
| L1_current_safe_core       | Leakage-aware AI neoantigen predictor benchmarking framework | supported               | 0 candidates; 0 methods; 0 split metrics                      | freeze input manifests and add public training-corpus row-level overlap files                            | new SOTA predictor                                         |
| L2_current_safe_algorithm  | Benchmark-adaptive reliability ranking with abstention       | supported               | BAR-Neo/BMA/stress-guarded layers; reviewer kill-audit pass=3 | prospective lockbox or independently held-out source validation                                          | external validation proven                                 |
| L3_current_execution_ready | pHLA assay-design candidates for manual research review      | supported_with_boundary | T1 assay-ready=3; translational-ready=0                       | WT peptide, gene, mutation, expression, clonality, patient context, HLA-LOH/B2M                          | clinical vaccine selection                                 |
| L4_near_term_upgrade       | Source/HLA robust reliability framework                      | partially_supported     | source/HLA stress boards and win/loss reports exist           | pre-register source-heldout/HLA-heldout success thresholds and rerun after public overlap audit          | Class I and Class II unified predictor                     |
| L5_future_translational    | Translational neoantigen research triage system              | blocked                 | patient-gated demo exists but patient metadata are sparse     | real PAAD/THCA patient metadata and at least one source-backed antigen/presentation/immune context layer | clinical utility                                           |
| L6_future_experimental     | Experimental support for selected T1 candidates              | not_started             | assay design matrix created                                   | binding/presentation/immunogenicity readout with WT comparator and controls                              | validated vaccine target without assay and safety evidence |

## Validation gaps

| gap_id                      | blocks_claim                       | current_status                                            | minimum_fix                                                                                      | success_criterion                                                                  | impact_if_fixed                                                                 |
|:----------------------------|:-----------------------------------|:----------------------------------------------------------|:-------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------|:--------------------------------------------------------------------------------|
| GAP01_public_overlap_audit  | public tools as clean comparators  | 9 public methods unresolved; clean public comparators=0   | drop public training corpora into project/data/public_training_corpora and rerun row-level audit | clean_comparator_allowed_after_row_audit=true or method remains caveated           | turns comparator section from caveated to auditable                             |
| GAP02_t1_antigen_identity   | T1 antigen identity                | T1 translational-ready=0; current T1 rows=3               | fill WT peptide, gene, mutation_id, protein/source window for each T1 row                        | all T1 rows have source-backed antigen identity fields                             | upgrades T1 from pHLA assay-design candidate to antigen identity-supported lead |
| GAP03_t1_antigen_expression | tumor antigen evidence             | expression/VAF/clonality absent in current T1 master rows | link expression_tpm, mutant_expression, VAF and clonality or mark not available                  | source-backed expression/clonality evidence or explicit no-data caveat             | enables translational triage score rather than benchmark-only score             |
| GAP04_presentation_safety   | presentation and safety gate       | HLA-LOH/B2M/processing/WT comparator absent               | link HLA-LOH, B2M, antigen processing status and WT comparator review                            | no known presentation hard fail; WT comparator caveat resolved                     | enables stronger pHLA-to-patient triage framing                                 |
| GAP05_independent_lockbox   | external/source-heldout robustness | source/HLA stress exists but not prospective lockbox      | freeze candidate/method manifests and evaluate on held-out source or new audit drop              | predefined AUPRC/top-k/calibration thresholds met without modifying method weights | moves from framework demonstration toward stronger validation paper             |

## Figure blueprint

| figure       | title                                         | panel_plan                                                                      | source_outputs                                                                                                     | nature_grade_message                                                                    |
|:-------------|:----------------------------------------------|:--------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------|
| Fig1         | CLEAN-NeoBench benchmark contract             | candidate schema; leakage flags; split contracts; public-tool caveat logic      | clean_neobench_master.tsv; clean_neobench_overlap_flags.tsv; CLEAN_NEOBENCH_METHOD_CARD.md                         | the contribution is a leakage-aware benchmark contract, not another black-box predictor |
| Fig2         | Method zoo under strict split contracts       | overall leaderboard; clean internal board; source/HLA heldout; calibration      | clean_neobench_leaderboard.tsv; clean_neobench_split_metrics.tsv; clean_neobench_winloss_method_summary.tsv        | performance is shown by split and failure mode, not pooled AUROC                        |
| Fig3         | BAR-Neo reliability and abstention            | BMA weights; contextual weights; abstention funnel; failure-aware downweighting | barneo_bma_method_weights.tsv; barneo_contextual_bma_method_weights.tsv; barneo_failure_aware_candidate_scores.tsv | benchmark behavior becomes calibrated reliability and abstention                        |
| Fig4         | Reviewer kill-audit and top-pass evidence     | 12 high-impact leads; pass/manual/blocked disposition; T1 method support        | barneo_high_impact_reviewer_kill_audit.tsv; barneo_top_pass_reviewer_evidence.tsv                                  | the system kills attractive but leaky candidates instead of hiding them                 |
| Fig5         | T1 translational readiness and execution plan | assay-ready vs metadata-blocked; intake template; assay matrix; go/no-go gates  | barneo_t1_translational_readiness.tsv; barneo_t1_assay_design_matrix.tsv; barneo_t1_go_nogo_criteria.tsv           | T1 candidates are converted into executable research tasks with claim boundaries        |
| ExtendedData | Distribution error and public overlap audit   | source prevalence shift; low-prevalence errors; public corpus audit intake      | clean_neobench_method_distribution_vulnerability.tsv; clean_neobench_public_overlap_audit_intake.tsv               | failure analysis and comparator caveats are first-class evidence                        |

## Reviewer attack matrix

| reviewer_attack                           | defense_now                                                                              | remaining_gap                                         | next_artifact                                          |
|:------------------------------------------|:-----------------------------------------------------------------------------------------|:------------------------------------------------------|:-------------------------------------------------------|
| Your public baselines are contaminated    | public pretrained methods are caveated; clean allowed count remains 0                    | row-level public training corpora are not yet dropped | clean_neobench_public_training_row_overlap_summary.tsv |
| Your model is just learning source labels | source-heldout/study-heldout metrics and source/HLA stress boards are explicit           | prospective lockbox or new source drop still needed   | frozen source-heldout validation report                |
| Top candidates are leakage artifacts      | reviewer kill-audit blocks 8/12 high-impact rows and keeps only 3 T1                     | manual provenance confirmation for T1 source identity | T1 metadata intake template completion                 |
| High score does not mean usable antigen   | T1 translational readiness separates pHLA assay-design from patient/translational claims | WT/gene/mutation/expression/clonality missing         | T1 antigen identity and expression evidence table      |
| This is clinical vaccine selection        | all reports explicitly state research triage and assay-design boundary                   | avoid clinical language in manuscript and demos       | claim boundary checklist                               |

## 30-day sprint

|   priority | workstream              | task                                                                                           | expected_output                                  | upgrade_unlocked                                  |
|-----------:|:------------------------|:-----------------------------------------------------------------------------------------------|:-------------------------------------------------|:--------------------------------------------------|
|          1 | T1 antigen identity     | Fill WT peptide, gene, mutation_id, source_protein_window for CNV0_02407/CNV0_02504/CNV0_02410 | barneo_t1_metadata_intake_template.tsv completed | T1 antigen identity-supported leads               |
|          2 | Public overlap audit    | Collect/drop public training corpora for NetMHCpan/MHCflurry/BigMHC/PRIME/MixMHCpred/etc.      | public training row-overlap audit summary        | auditable comparator framing                      |
|          3 | T1 assay planning       | Prepare peptide synthesis/HLA-A*02:01 binding-stability planning table                         | assay design matrix with owner/status fields     | research assay execution packet                   |
|          4 | Lockbox validation      | Freeze runner outputs and evaluate on new held-out source or locked patient/study split        | source/HLA lockbox validation report             | stronger external robustness claim                |
|          5 | Manuscript architecture | Build figure panels from current TSVs without writing voice-protected prose                    | figure blueprint and factual result blocks       | Nature Methods/NBME-style resource paper scaffold |

## 결론

지금은 `Nature-grade methods/resource core`다. 다음 업그레이드는 T1 3개에 WT/gene/mutation/expression/clonality/patient/presentation metadata를 붙이고, public overlap audit과 source/HLA lockbox를 완료하는 것이다.

## Claim boundary

clinical vaccine selection, new SOTA predictor, external validation proven, quantum advantage claim은 금지한다.
