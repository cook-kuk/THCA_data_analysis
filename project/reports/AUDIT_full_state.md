# AUDIT_full_state

## 0. Audit Metadata

| Item | Value |
|---|---|
| 생성 일시 UTC | `2026-04-26 23:21:16 UTC` |
| 생성 일시 KST | `2026-04-27 08:21:16 KST` |
| Git HEAD | `fe719967dad767d3fd207878c5700afda9a3244f` |
| Repo size | `project/ = 3.8G` |
| Python | `3.12.3` |
| pandas / numpy / scikit-learn / scanpy / plotly | `2.3.3 / 2.4.4 / 1.8.0 / 1.12.1 / 6.7.0` |

최근 수정 파일 top 10:

| Rank | Path | Modified |
|---|---|---|
| 1 | `project/logs/http_secure.log` | `2026-04-27 06:01` |
| 2 | `project/reports/v17p3/v17p3_response_to_reviewers_v2.md` | `2026-04-27 02:03` |
| 3 | `project/reports/v17p3/v17p3_paper_outline_FINAL.md` | `2026-04-27 02:03` |
| 4 | `project/reports/v17p3/v17p3_external_review_dump_v3.md` | `2026-04-27 02:03` |
| 5 | `project/reports/v17p3/index.html` | `2026-04-27 02:03` |
| 6 | `project/results/v17p3/tables/A3_tumor_predicted_response.tsv` | `2026-04-27 02:02` |
| 7 | `project/results/v17p3/tables/A3_top_drugs_dm2_selective.tsv` | `2026-04-27 02:02` |
| 8 | `project/results/v17p3/tables/A3_top_drugs_dm1_selective.tsv` | `2026-04-27 02:02` |
| 9 | `project/results/v17p3/tables/A3_summary.json` | `2026-04-27 02:02` |
| 10 | `project/results/v17p3/tables/A3_pathway_drug_coherence.tsv` | `2026-04-27 02:02` |

## 1. Repo Snapshot

### 1a. 결과물 규모

| Directory | Files |
|---|---:|
| `project/results/v14_strengthening` | 14 |
| `project/results/v15_neurips` | 14 |
| `project/results/v17` | 65 |
| `project/results/v17p2` | 48 |
| `project/results/v17p3` | 67 |
| `project/reports/v17` | 8 |
| `project/reports/v17p2` | 4 |
| `project/reports/v17p3` | 4 |
| `project/reports/html` | 709 |

### 1b. 핵심 리포트 경로

| Phase | Main report |
|---|---|
| v17 Phase 1 | `project/reports/v17/index.html` |
| v17 Phase 2 | `project/reports/v17p2/index.html` |
| v17 Phase 3 | `project/reports/v17p3/index.html` |
| Full audit | `project/reports/AUDIT_full_state.md` |

### 1c. 메타데이터 핵심 파일

| File | Shape / Size |
|---|---|
| `project/metadata/sample_master_v3.tsv` | `1509 x 27` |
| `project/results/v17/tables/sample_master_v17_full.tsv` | `1509 x 39` |
| `project/metadata/tierA67_genes.txt` | TierA panel |
| `project/metadata/tds16_genes.txt` | TDS16 panel |
| `project/metadata/brs71_genes.txt` | BRS71 panel |

## 2. Sample Master 진단

### 2a. `sample_master_v3` vs `sample_master_v17_full`

| Item | Value |
|---|---|
| `sample_master_v3` row 수 | `1509` |
| `sample_master_v17_full` row 수 | `1509` |
| 새 컬럼 수 | `12` |
| 새 컬럼 | `braf_v600e, driver_anchor_v17, fusion_classes, fusion_pair, mutation_genes, ras_hotspot, tcga12, tert_promoter, tert_status, v17_dark_cluster, v17_data_source, v3_anchor_6class` |

### 2b. TCGA tumor 기준 `driver_anchor_v17` 분포

| Group | n |
|---|---:|
| `BRAF` | 284 |
| `RAS` | 54 |
| `unknown` | 164 |
| `DICER1_EIF1AX_PPM1D` | 7 |
| `RET_fusion` | 1 |
| `NTRK_fusion` | 1 |
| `ALK_fusion` | 1 |
| `TP53` | 1 |

### 2c. Dark Matter 정의 일관성

| Definition | n |
|---|---:|
| Phase 1 Dark Matter cohort (`driver_anchor not in BRAF/RAS`) | `178` |
| Current `driver_anchor_v17 == unknown` only | `164` |

해석:
- `178`은 Phase 1의 원래 Dark Matter 분석 cohort 크기다.
- `164`는 Phase 3 시점 `driver_anchor_v17`로 extra driver와 fusion을 재배치한 뒤 남은 `true unknown`이다.
- 즉 현재 taxonomy에서는 `14`명이 Dark Matter bucket에서 다른 anchor로 재분류됐다.

## 3. v14 Baseline 상태

| Item | Status |
|---|---|
| TierA67_clean panel | 위치: `project/metadata/tierA67_genes.txt` |
| TDS16 panel | 위치: `project/metadata/tds16_genes.txt` |
| BRS71 panel | 위치: `project/metadata/brs71_genes.txt` |
| 5-cohort harmonization | `project/results/v14_strengthening/` |
| DIAL framework code | `project/notebooks_or_scripts/v17_idea4_dial_audit.py`, `v5_dial_*`, `v5p1_dial.py`, `v5p2_dial_proper.py` 계열 |

v14 핵심 요약:
- `TierA67_clean` 54-gene panel로 BRAF-like vs RAS-like 표현형 축을 복원했다.
- 5-cohort 통합과 ComBat/transfer audit를 통해 cross-cohort generalization을 점검했다.
- DIAL(direction-invariant audit)로 단순 AUC가 아닌 identifiable transfer를 정량화했다.

## 4. v15 상태

| Item | Status |
|---|---|
| 주제 | `BRAF × TROP2 × thyroid` NeurIPS draft |
| 데이터 | `TCGA, CCLE, PRISM, GTEx, GSE58545` |
| 산출 위치 | `project/reports/v15_neurips/` |
| 현재 단계 | `draft 단계` |
| v17과의 관계 | 같은 THCA 기반이지만 narrative는 분리 가능 |

분리도:
- v15는 `TROP2/ADC/translational drug angle` 중심이다.
- v17은 `DM1/DM2 transcriptional taxonomy` 중심이다.
- 중첩 데이터는 있으나 claim은 동일하지 않다.

## 5. v17 Phase 1

### 5a. 산출 inventory

| Item | Count |
|---|---:|
| `project/results/v17/figs/*.html` | `23+` |
| `project/results/v17/tables/*.tsv|json` | `18+` |
| `project/reports/v17/*` | `8` |

### 5b. 핵심 메트릭

| Metric | Value |
|---|---|
| Dark Matter n | `178` |
| Best K | `2` |
| Consensus stability | `0.9788` |
| DM1 n | `109` |
| DM2 n | `69` |
| DM1 age mean | `41.82` |
| DM2 age mean | `54.71` |
| DM1 TDS16 mean | `6.81` |
| DM2 TDS16 mean | `8.13` |
| DM1 RAI mean | `7.58` |
| DM2 RAI mean | `9.20` |

대표 marker:

| Cluster | Top markers |
|---|---|
| DM1-like (`cluster 0`) | `DUSP5(+2.79), DUSP6(+2.38), DUSP4(+2.16), MET(+2.13), FOXP3(+1.10)` |
| DM2-depleted in same contrast | `SLC5A8(-3.02), DIO1(-3.89), TPO(-3.91), DIO2(-2.12)` |

Phase 1 한계:
- `TERT` 실질적 비관측
- fusion은 `proxy` 재사용
- external robust transfer는 사실상 `2 cohort`
- TCGA survival power 부족

## 6. v17 Phase 2

### 6a. 산출 inventory

| Item | Count |
|---|---:|
| `project/results/v17p2/figs` | `21` |
| `project/results/v17p2/tables` | `21` |
| `project/reports/v17p2/*` | `4` |

### 6b. Layer별 핵심 메트릭

| Layer | Metric | Value |
|---|---|---|
| L1 Robustness | median persistence | `0.9035` |
| L1 Robustness | 5-method mean ARI | `0.6624` |
| L1 Robustness | permutation p | `0.0033` |
| L1 Robustness | 50% subsample ARI | `0.7571` |
| L2 Clinical | age Welch p | `6.19e-07` |
| L2 Clinical | stage high p | `0.0919` |
| L2 Clinical | OS events | `8` |
| L3 Biology | Hallmark FDR<0.05 | `10` |
| L3 Biology | MAPK DM2-DM1 | `-1.0801` |
| L4 External | cohorts tested | `2` |
| L4 External | DIA-AUC > 0.85 cohorts | `2` |
| L4 External | scRNA cells scored | `66015` |
| L5 DIAL | ComBat identifiability<0.5 threshold | `lambda=1.0` |
| L5 DIAL | pancancer tested | `3` |
| L5 DIAL | max DIA-AUC | `0.9912` |

Phase 2 결론:
- DM1/DM2가 random artifact는 아니라는 통계적 방어는 확보했다.
- 그러나 `external bulk 2-cohort`, `survival power 부족`, `biology proxy fallback`이 남았다.

## 7. v17 Phase 3

### 7a. 산출 inventory

| Item | Count |
|---|---:|
| `project/results/v17p3/tables` | `42` |
| `project/results/v17p3/figs` | `24` |
| `project/reports/v17p3/*` | `4` |

### 7b. Fix task 결과

| Task | Goal | Result | 판정 |
|---|---|---|---|
| F1 external recovery | `GSE213647 overlap recovery` | TierA overlap `55`, robust cohort 수 `2` | `부분 성공` |
| F2 proper GSEA | `MSigDB proxy 제거` | Hallmark FDR<0.01 `15`, Reactome FDR<0.05 `19` | `성공` |
| F3 clinical reframe | `OS 대신 대체 endpoint` | 유의 endpoint `4`, best=`rai_score_recalc` | `성공` |

F1 cohort별 핵심:

| Cohort | Best observed DIA-AUC | Note |
|---|---:|---|
| `GSE27155` | `0.8491` | SVM-RBF가 최고 |
| `GSE76039` | `1.0000` | LogReg 최고 |
| `GSE213647` | `0.9852` | gene-ID recovery 후 매우 강함 |

해석:
- `GSE213647 gene overlap 0` 문제는 해결됐다.
- 다만 robust cohort count summary는 여전히 `2`로 남아 있어, narrative상 `3+ cohort robust`라고 쓰면 과장이다.

### 7c. Amplify task 결과

| Task | 핵심 메트릭 | 판정 |
|---|---|---|
| A1 scRNA heterogeneity | `7 patients`, mixed `6`, dominance `0.757–0.996` | `성공` |
| A2 DM in BRAF/RAS+ | corr(BRAF, probDM2)=`-0.578`, `BRAF/DM2-like=2`, `RAS/DM1-like=15` | `성공` |
| A3 drug response | FDR<0.1 selective drug `0` | `실패/약함` |
| A4 pan-cancer | applicable cancers `5` | `성공` |
| A5 genomic/immune | TMB DM2-DM1 `0.092` | `보강은 됐으나 약함` |
| A6 TF/RAI model | PDTC low-RAI AUC `0.0118` | `실패` |

### 7d. Phase 3에서 실제로 강해진 claim

| Claim | 상태 |
|---|---|
| external gene-ID failure는 technical artifact였다 | `강화` |
| biology layer는 proper GSEA로 방어 가능하다 | `강화` |
| DM1/DM2는 age/histology/RAI-relevant axis다 | `강화` |
| DM axis는 BRAF/RAS mutation axis와 직교한다 | `강화` |
| single-cell 수준에서도 mixed-state dominance로 보인다 | `강화` |
| drug actionability가 강하다 | `강화 실패` |
| RAI prediction model이 강력하다 | `강화 실패` |

## 8. Cross-Phase 일관성 점검

| Item | Status |
|---|---|
| Phase 1 DM cohort size | `178` |
| Phase 2 DM cohort size | `178` |
| Phase 3 base DM cohort size | `178` |
| Current `true unknown` after v17 anchor update | `164` |
| TierA panel | 동일 panel 계열 유지 |
| Phase 1 vs 2 vs 3 narrative | `발견 → 검증 → 임팩트 증폭`으로 일관 |

주의:
- Phase 3에서 `driver_anchor_v17`가 업데이트되면서 original Dark Matter `178`과 current `unknown 164`를 혼용하면 안 된다.
- 논문 본문에서는 `discovery cohort (n=178)`와 `post-expansion true unknown (n=164)`를 구분해서 써야 한다.

## 9. Reviewer 공격 vs 현재 답변 가능도

| 공격 | 답변 근거 | 강도 |
|---|---|---|
| Cluster artifact | Phase 2 L1 persistence/permutation/subsample | `강` |
| Age confound | age Welch + multivariable reframing | `중강` |
| Mechanism 부재 | Phase 3 proper GSEA + TF/thyroid axis | `중강` |
| TCGA-only | `GSE76039`, `GSE213647` 강함, 나머지는 약함 | `중` |
| Method 검증 없음 | DIAL L5 + pancancer feasibility | `중강` |
| Sample size 작음 | bulk는 작고 scRNA 보강 존재 | `중` |
| Batch effect | DIAL/ComBat sensitivity 있음 | `중강` |
| Clinical relevance | age/histology/RAI score는 강함, OS는 약함 | `중` |
| Hypothesis-driven 아님 | de-differentiation / driver-orthogonal axis로 정리 가능 | `중강` |
| Reproducible | Phase 1-3 산출 정리, HTML/report 다수 | `강` |
| Drug actionability 부족 | A3 failed | `약` |
| RAI predictor validation 부족 | A6 failed | `약` |

## 10. Self-Assessment Score

| Item | Score (0-10) | Rationale |
|---|---:|---|
| Sample size per cluster | 7 | DM1/DM2 `109/69`는 최소 논문화 가능 |
| Cross-cohort validation | 6 | 강한 cohort는 2개, recovery는 됐지만 완전하지 않음 |
| Method novelty | 8 | DIAL + direction-invariant audit는 차별점 있음 |
| Biological discovery | 8 | MAPK/differentiation/driver-orthogonal axis가 명확 |
| Wet validation | 0 | 부재 |
| Clinical actionability | 4 | A3/A6가 약해 translational hook이 약함 |
| Korean cohort | 0 | 부재 |
| Statistical rigor | 8 | Phase 2 robustness와 Phase 3 proper GSEA로 보강 |
| Mechanism interpretability | 7 | Hallmark/Reactome/thyroid TF axis는 설명 가능 |
| Reproducibility | 8 | script/results/report가 전부 남아 있음 |

총점: `56 / 100`

## 11. Venue Probability

| Venue | Probability | 이유 |
|---|---:|---|
| Bioinformatics short | 85% | discovery + method + audit는 충분 |
| npj Precision Oncology | 65% | taxonomy + biology + external 일부 + translational hint |
| Genome Medicine | 35% | actionability / validation depth가 아직 부족 |
| Nature Communications | 10% | wet validation/Korean cohort/clinical endpoint 부족 |
| Nature Cancer | 3% | 현 상태로는 과도 |
| Cancer Cell | 1% | 불가에 가까움 |

## 12. 현재 시점의 솔직한 결론

- **지금 strongest narrative**: `BRAF/RAS mutation axis와 직교하는 DM1/DM2 transcriptional taxonomy`  
  `proper GSEA`, `age/histology/TDS/RAI axis`, `single-cell dominance`, `partial external recovery`로 방어 가능.
- **지금 weakest narrative**: `drug repositioning`, `RAI predictor`, `hard clinical outcome`.
- 따라서 main paper는 `biology + taxonomy + audit` 중심으로 가고, `A3/A6`는 supplementary 혹은 future-work로 내리는 편이 안전하다.

## 13. 다음 결정 5개

1. `v17` main paper target을 `npj Precision Oncology` 1차로 둘지 결정
2. `Genome Medicine` 도전용으로 actionability를 더 보강할지 결정
3. `A3` drug response를 더 고칠지, 현재는 supplementary로 둘지 결정
4. `A6` RAI model은 현재 결과가 역방향이므로 narrative에서 비중 축소할지 결정
5. Korean cohort 답 전 submit할지, external depth를 더 쌓을지 결정

## 14. Critical Next Action 3개

본인이 직접 해야 할 것:
1. `v17p3_external_review_dump_v3.md`를 기준으로 paper narrative 우선순위를 정하기
2. `A3/A6`를 main claim에서 뺄지 유지할지 결정하기
3. Korean cohort / 협업자 피드백 타이밍 결정하기

Claude Code에 위임 가능한 것:
1. `title + abstract + cover letter` 초안 작성
2. `Figure 1-7 main vs supplementary` 재배치
3. `Results / Discussion` 논문 문단 초안 작성

외부 멘토에게 던질 질문:
1. `npj`로 충분한가, 아니면 `Genome Medicine`을 위해 무엇이 더 필요한가
2. `A3/A6` 실패를 limitation으로 둘지 supplementary로 둘지
3. `164 true unknown`와 `178 discovery cohort`를 본문에서 어떻게 가장 깔끔하게 서술할지
