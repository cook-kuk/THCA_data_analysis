# THCA markdown organizer

Generated: 2026-05-14  
Scope: repo 안의 md 문서를 프로젝트별로 묶어 정리한 실행형 인덱스.  
Read coverage: `paper1`, `paper2 image/histology`, `paper2 HLA bridge`, `paper4 HLA germline`, `Samsung CTC-EMT`, `npj submission`, plus repo-level portfolio docs.  
Excluded from the pass: `.venv`, `node_modules`, `vendor`, package README noise.

## 0. 먼저 열 문서

| order | file | reason |
|---:|---|---|
| 1 | `project/HANDOFF_2026_05_09.md` | 최신 handoff. 지금 무엇이 끝났고 무엇이 남았는지 한 번에 잡힘 |
| 2 | `project/results/paper_portfolio_2026_05_09/PAPER_PORTFOLIO_ROADMAP_KR.md` | 전체 우선순위. 무엇을 지금 ship 해야 하는지 |
| 3 | `project/results/paper_repo_portfolio_2026_05_10/PAPER_REPO_MAP_V2_KR.md` | branch / boundary / software ledger |
| 4 | `project/REPRODUCIBILITY.md` | 재현성 기준, data/path 규칙 |
| 5 | `project/STATUS_PAPER1_2026_05_02.md` | Paper 1 현재 상태 요약 |
| 6 | `project/submission/npj/SUBMISSION_READY.md` | npj submission 패키지의 상태와 버전 |

## 1. 프로젝트별 정리

### A. Paper 1 / DM1

**한 줄:** 갑상선 cancer dark matter / 8-gene panel / fusion-driven silencing trunk.

**현재 상태:** `submit now` 급. Manuscript voice만 author가 마무리하면 되는 구간.

**문서 수:** `37 md`

**핵심 md**
- `project/manuscript_v8/README.md`
- `project/manuscript_v8/02_outline.md`
- `project/manuscript_v8/04_results.md`
- `project/manuscript_v8/05_figure_captions.md`
- `project/manuscript_v8/07_star_methods.md`
- `project/manuscript_v8/09_reviewer_qa.md`
- `project/manuscript_v8/10_full_manuscript_compiled.md`
- `project/manuscript_v8/13_supplementary_tables.md`
- `project/manuscript_v8/14_self_verification_report.md`
- `project/manuscript_v8/16_paper_quality_upgrade_pack.md`

**읽는 순서**
1. `README.md`로 구조 확인
2. `02_outline.md`로 figure / result spine 확인
3. `04_results.md`와 `05_figure_captions.md`로 서사 연결
4. `09_reviewer_qa.md`로 reviewer 방어 확인
5. `07_star_methods.md`와 `13_supplementary_tables.md`로 제출형식 확인

**관련 웹 페이지**
- `project/papers_hub_2026_05_04/paper1.html`
- `project/papers_hub_2026_05_04/paper1_reviewer_defense_dashboard.html`
- `project/papers_hub_2026_05_04/paper1_fig8_mechanism_dossier.html`
- `project/papers_hub_2026_05_04/paper1_methods_reproducibility_dossier.html`

---

### B. Paper 2 / histology / image-DM1

**한 줄:** tissue image branch. CLAM / foundation model / K2 validation / pathology-spatial bridge.

**현재 상태:** `PASS_LAUNCH_PAPER2`.

**문서 수:** `72 md`

**핵심 md**
- `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/SPRINT_RESULT.md`
- `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/PAPER2_OUTLINE_2026_05_08.md`
- `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/PAPER2_KOREAN_K2_VALIDATION_PLAN_2026_05_08.md`
- `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/GIGATIME_COMPARISON_2026_05_08.md`
- `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/manuscript_draft/02_methods_results.md`
- `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/manuscript_draft/08_cover_letter_scaffold.md`
- `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase1_gse250521/PHASE1_REPORT.md`
- `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam/PHASE2_REPORT.md`
- `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase3_integration/SPRINT_RESULT.md`
- `project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp/paper2_external_validation_ready_2026_05_10/SUMMARY.md`

**읽는 순서**
1. `SPRINT_RESULT.md`
2. `PAPER2_OUTLINE_2026_05_08.md`
3. `manuscript_draft/02_methods_results.md`
4. `analysis_supp/paper2_external_validation_ready_2026_05_10/SUMMARY.md`
5. `analysis_supp/paper2_impact_upgrade_2026_05_10/SUMMARY.md`

**관련 웹 페이지**
- `project/papers_hub_2026_05_04/paper2_image_dm1.html`
- `project/papers_hub_2026_05_04/paper2_cv2_visual_summary.html`
- `project/papers_hub_2026_05_04/paper2_external_validation_ready_pack.html`
- `project/papers_hub_2026_05_04/paper2_master_view.html`

---

### C. Paper 2 / HLA bridge / HT-overlap

**한 줄:** cancer-HLA bridge. allele claim은 quarantine, transcriptomic immune-context만 사용.

**현재 상태:** HT-overlap expression / TLS / antigen-presentation axis는 강함. allele-level cancer claim은 금지.

**문서 수:** `22 md`

**핵심 md**
- `project/manuscript_p2_brief/README.md`
- `project/manuscript_p2_brief/SITUATION_OVERVIEW.md`
- `project/manuscript_p2_brief/DATA_SOURCES_INDEX.md`
- `project/manuscript_p2_brief/COHORT_ACCESS_GUIDE.md`
- `project/manuscript_p2_brief/REFERENCES_BIB_AUDIT.md`
- `project/manuscript_p2_brief/PROMPT_DECISION_LOG.md`
- `project/manuscript_p2_brief/CHANGELOG.md`
- `project/manuscript_p2_brief/POST_COMMIT_STATUS.md`
- `project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md`

**읽는 순서**
1. `README.md`
2. `HLA_CANCER_SEPARATION_RULES.md`
3. `SITUATION_OVERVIEW.md`
4. `COHORT_ACCESS_GUIDE.md`
5. `DATA_SOURCES_INDEX.md`

**관련 웹 페이지**
- `project/papers_hub_2026_05_04/paper2_hla.html`
- `project/papers_hub_2026_05_04/paper2_submission_overview.html`
- `project/papers_hub_2026_05_04/paper2_high_impact_journal_strategy.html`

---

### D. Paper 4 / HLA germline / Graves-AITD

**한 줄:** Korean / Pan-Asian HLA architecture for autoimmune thyroid disease.

**현재 상태:** DPB1*05:01 중심 architecture는 강함. Korean adult GD NGS만 있으면 더 커짐.

**문서 수:** `9 md`  
(`project/results/p2_pillar1_forest`: 3 md, `project/results/p2_pillar1_forest_v2`: 6 md)

**핵심 md**
- `project/manuscript_hla_ncomm_upgrade_2026_05_09/PAPER4_GD_HLA_NCOMM_SKELETON.md`
- `project/manuscript_hla_ncomm_upgrade_2026_05_09/PAPER2_HTPTC_IMMUNE_NCOMM_SKELETON.md`
- `project/results/p2_pillar1_forest_v2/PILLAR1_FOREST_V2_SUMMARY.md`
- `project/results/p2_pillar1_forest_v2/paper2_hla_claim_boundary.md`
- `project/results/p2_pillar1_forest_v2/paper2_hla_metric_sensitivity_report.md`
- `project/results/p2_pillar1_forest_v2/paper2_hla_multiple_testing_report.md`
- `project/results/p2_pillar1_forest_v2/paper2_hla_subcohort_heterogeneity_report.md`
- `project/results/p2_pillar1_forest/PILLAR1_FOREST_SUMMARY.md` *(deprecated v1, keep for record only)*
- `project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md`

**읽는 순서**
1. `PAPER4_GD_HLA_NCOMM_SKELETON.md`
2. `PILLAR1_FOREST_V2_SUMMARY.md`
3. `HLA_CANCER_SEPARATION_RULES.md`
4. `paper2_hla_multiple_testing_report.md`
5. `paper2_hla_subcohort_heterogeneity_report.md`

**관련 웹 페이지**
- `project/papers_hub_2026_05_04/paper4_gd_hla.html`
- `project/papers_hub_2026_05_04/hla_ncomm_upgrade_2026_05_09.html`
- `project/papers_hub_2026_05_04/hla_execution_board_2026_05_09.html`

---

### E. Samsung CTC-EMT proposal

**한 줄:** thyroidectomy-induced CTC-EMT state transition proposal.

**현재 상태:** public pilot은 feasibility / marker-selection까지만. CTC biology proof는 아직 아님.

**문서 수:** `58 md`

**핵심 md**
- `kthyro_ctc_emt_samsung_proposal/README.md`
- `kthyro_ctc_emt_samsung_proposal/inputs/source_context_index.md`
- `kthyro_ctc_emt_samsung_proposal/outputs/reports/FINAL_IMPACT_BRIEF_KR.md`
- `kthyro_ctc_emt_samsung_proposal/outputs/reports/MANUSCRIPT_FLOW_MASTER_KR.md`
- `kthyro_ctc_emt_samsung_proposal/outputs/reports/MANUSCRIPT_BLUEPRINT_KR.md`
- `kthyro_ctc_emt_samsung_proposal/outputs/reports/MANUSCRIPT_EXECUTION_PACKAGE_KR.md`
- `kthyro_ctc_emt_samsung_proposal/outputs/reports/FINAL_SUBMISSION_WAR_ROOM_KR.md`
- `kthyro_ctc_emt_samsung_proposal/outputs/reports/CLAIM_BOUNDARY_NO_YU_DATA.md`
- `kthyro_ctc_emt_samsung_proposal/outputs/reports/FUTURE_DATA_UNLOCK_MAP_KR.md`
- `kthyro_ctc_emt_samsung_proposal/outputs/reports/samsung_ctc_emt_specific_aims_kr.md`
- `kthyro_ctc_emt_samsung_proposal/outputs/reports/samsung_ctc_emt_preliminary_results_3page_kr.md`
- `kthyro_ctc_emt_samsung_proposal/outputs/reports/samsung_ctc_emt_reviewer_defense_table_kr.md`

**읽는 순서**
1. `README.md`
2. `outputs/reports/FINAL_IMPACT_BRIEF_KR.md`
3. `outputs/reports/CLAIM_BOUNDARY_NO_YU_DATA.md`
4. `outputs/reports/MANUSCRIPT_FLOW_MASTER_KR.md`
5. `outputs/reports/samsung_ctc_emt_specific_aims_kr.md`

**관련 웹 페이지**
- `project/papers_hub_2026_05_04/kthyro_ctc_emt_samsung_science_dossier.html`
- `project/papers_hub_2026_05_04/kthyro_ctc_emt_final_impact_brief.html`
- `project/papers_hub_2026_05_04/kthyro_ctc_emt_paper_flow_master_dossier.html`
- `project/papers_hub_2026_05_04/kthyro_ctc_emt_submission_war_room.html`

---

### F. npj submission / submission bundle

**한 줄:** DM1 submission bundle. current manuscript v7 and defense assets.

**현재 상태:** submission-ready bundle. 버전 관리가 핵심.

**문서 수:** `43 md`

**핵심 md**
- `project/submission/npj/README.md`
- `project/submission/npj/SUBMISSION_READY.md`
- `project/submission/npj/SUBMISSION_CHECKLIST.md`
- `project/submission/npj/SUBMIT_INSTRUCTIONS.md`
- `project/submission/npj/SUBMIT_INSTRUCTIONS_v7.md`
- `project/submission/npj/manuscript_v7.md`
- `project/submission/npj/manuscript_v6_ULTIMATE.md`
- `project/submission/npj/cover_letter_v6_ULTIMATE.md`
- `project/submission/npj/reviewer_bundle.zip` 관련 설명 문서들
- `project/submission/npj/voice_polish_audit.md`
- `project/submission/npj/v17p35_REVIEWER_DEFENSE.md`

**읽는 순서**
1. `SUBMISSION_READY.md`
2. `SUBMISSION_CHECKLIST.md`
3. `SUBMIT_INSTRUCTIONS_v7.md`
4. `manuscript_v7.md`
5. `voice_polish_audit.md`

**관련 웹 페이지**
- `project/papers_hub_2026_05_04/paper1.html`
- `project/papers_hub_2026_05_04/paper1_reviewer_defense_dashboard.html`

---

### G. Separate universe / Neoantigen / CROSS-Neo

**한 줄:** neoantigen / vaccine / ICI platform. thyroid main track와 분리.

**현재 상태:** separate universe. cancer vaccine / combo-therapy page는 참고용.

**대표 md**
- `project/results/p_neo_ici_combo_2026_05_11/NEOICI_ALGORITHM_GAUNTLET_SUMMARY.md`
- `project/results/p_neo_ici_nature_cancer_package_2026_05_11/NEOICI_NATURE_CANCER_PACKAGE_SUMMARY.md`
- `project/results/p_neo_ici_nature_cancer_synthesis_2026_05_11/NEOICI_NATURE_CANCER_SYNTHESIS_SUMMARY.md`
- `project/results/paper11_pancancer/README.md`
- `project/results/transfer_packages/cross_neo_cancer_vaccine_2026_05_09/README_STATUS_KR.md`

**관련 웹 페이지**
- `project/papers_hub_2026_05_04/neoici_nature_cancer_package.html`
- `project/papers_hub_2026_05_04/cancer_vaccine_full_dossier.html`
- `project/papers_hub_2026_05_04/cross_neo_total_impact_v8.html`

---

## 2. Deprecated / archive / mirror only

- `project/results/p2_pillar1_forest/PILLAR1_FOREST_SUMMARY.md` — v1 deprecated, v2가 현재 기준.
- `project/manuscript_p2_brief/lit_enrich_2026_05_02/*` — literature enrichment background, 현재 중심 문서는 아님.
- `project/manuscript_v8/_archive_2026_05_03_sprint_draft.md` — archive.
- `project/results/transfer_packages/*` — mirror / handoff only. source of truth 아님.
- `project/submission/npj/manuscript_v1.md` ~ `manuscript_v6_ULTIMATE.md` — version diff 용도.
- `project/papers_hub_2026_0509/*` — older hub snapshot. 지금은 `project/papers_hub_2026_05_04/` 쪽이 live 중심.

## 3. 지금 회의용 우선순위

1. Samsung CTC-EMT
   - `FINAL_IMPACT_BRIEF_KR.md`
   - `CLAIM_BOUNDARY_NO_YU_DATA.md`
   - `MANUSCRIPT_FLOW_MASTER_KR.md`

2. Paper 1 / DM1
   - `README.md`
   - `02_outline.md`
   - `09_reviewer_qa.md`
   - `paper1_reviewer_defense_dashboard.html`

3. Paper 2 image / histology
   - `SPRINT_RESULT.md`
   - `PAPER2_OUTLINE_2026_05_08.md`
   - `paper2_external_validation_ready_pack.html`

4. HLA bridge / HLA germline
   - `HLA_CANCER_SEPARATION_RULES.md`
   - `PILLAR1_FOREST_V2_SUMMARY.md`
   - `PAPER4_GD_HLA_NCOMM_SKELETON.md`

5. npj submission
   - `SUBMISSION_READY.md`
   - `manuscript_v7.md`
   - `SUBMIT_INSTRUCTIONS_v7.md`

## 4. 메모

- 이 인덱스는 “모든 md를 하나씩 읽은 뒤”의 요약본이 아니라, **실제 작업에 필요한 md만 추려 정리한 실행형 지도**다.
- vendored dependency 문서와 `.venv` README는 의도적으로 제외했다.
- HLA는 `HLA_CANCER_SEPARATION_RULES.md`를 boundary contract로 둔다.
- Paper 1은 cancer, Paper 2 image는 histology, Paper 2 HLA bridge는 transcriptomic immune-context, Paper 4는 germline HLA architecture로 분리한다.
