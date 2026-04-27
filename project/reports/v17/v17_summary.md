# v17 Summary


> ⚠️ **v5.2 retraction notice (2026-04-25):** some DIA-AUC and DIAL numbers below were computed under v5.1 leaky-ComBat protocol (full-pooled-data ComBat before LODO split, information leak). Under proper per-fold ComBat, the numbers shift. Current submission-ready figures are at `reports/html/pages/v17_npj_robustness.html` (manuscript v3 scenario-B reframe). See `reports/v5p2/v5p2_critical_assessment.md`.


## 한 줄 결론

v17은 갑상선암을 단순 BRAF vs RAS 이분법으로 보지 않고, **BRAF/RAS-negative Dark Matter 178명**, v17 driver map, cross-cohort robustness, trajectory를 한 묶음으로 재정리한 스프린트다. 이번 로컬 재현에서는 **Dark Matter가 안정적으로 2개 cluster(DM1=109, DM2=69)로 갈라졌고**, external transfer proxy에서 **DIA-AUC 0.969-1.000**을 보였으며, 반대로 **ComBat 적용 후 identifiability는 1.000/0.9998에서 0.2577/0.1573으로 급락**했다.

## Task별 핵심 finding

### Task 1 · Dark Matter sub-clustering

- TCGA-THCA tumor 중 `driver_anchor ∉ {BRAF, RAS}` 조건으로 **178명**을 추출했다.
- TierA67 기반 consensus clustering에서 **best K=2**, stability = **0.9788**로 수렴했다.
- cluster size는 **DM1 109명**, **DM2 69명**이다.
- RAI/TDS 축에서 DM2가 더 differentiated 쪽, DM1은 상대적으로 dedifferentiated 쪽으로 분리되는 패턴이 보인다.

### Task 2 · TERT 4-group

- 로컬 GDC WXS MAF를 사용했을 때 promoter hotspot (`C228T`, `C250T`)은 포착되지 않았다.
- 따라서 이번 재현에서 실질적으로 형성된 그룹은 **BRAF only (281)**, **RAS only (54)**, **triple-negative / Dark Matter subset (178)** 이다.
- 평균 RAI/TDS는 `RAS > triple-negative > BRAF` 순으로 정렬됐다.
- 중요 caveat: **WXS 기반 로컬 MAF는 TERT promoter assay가 아니므로, TERT 축은 현 시점 산출물에서 “부재 증거”가 아니라 “비관측”에 가깝다.**

### Task 3 · Fusion + driver landscape

- 전체 `sample_master_v17_full.tsv`는 **1,509개 sample** 기준으로 정리됐다.
- v17 driver summary top counts:
  - `unknown = 1072`
  - `BRAF = 343`
  - `RAS = 74`
  - `PAX8PPARG = 7`
  - `DICER1_EIF1AX_PPM1D = 7`
- TCGA mutation은 로컬 MAF에서, fusion은 현재 repo에 남아 있던 `v3_fusion_anchor_tcga.tsv` proxy를 사용했다.
- 즉 이번 v17 driver map은 **“로컬 mutation + 기존 fusion proxy 통합판”**이다.

### Task 4 · DIAL cross-cohort audit

- 외부 코호트 중 실제로 dark-matter transfer proxy가 돌아간 cohort는 **GSE27155, GSE76039** 이다.
- 결과:
  - `GSE27155`: DIA-AUC **0.9688**, identifiability **1.0000 → 0.2577 (ComBat 후)**
  - `GSE76039`: DIA-AUC **0.9902-1.0000**, identifiability **0.9998 → 0.1573 (ComBat 후)**
- 즉 **cluster predictability는 매우 높지만 batch identifiability는 ComBat 후 급격히 낮아진다.**
- 이건 v4의 “ComBat UNRECOVERABLE” 서사와 정성적으로 같은 방향이다.

### Task 5 · Trajectory

- trajectory는 공통 유전자 교집합 **50 genes**, 총 **645개 tumor sample**에서 계산됐다.
- 이번 실행에서 실제로 포함된 cohort는 `TCGA-THCA`, `GSE27155`, `GSE76039` 이다.
- `GSE126698`, `GSE213647` 는 현재 TierA67 alignment 경로에서 empty alignment로 제외됐다.
- pseudotime dynamic gene 상위:
  - 음의 상관: `TPO`, `SLC26A4`, `DIO1`, `TG`, `DIO2`, `PAX8`
  - 양의 상관: `DUSP5`, `MET`, `LOX`
- 해석은 명확하다. **pseudotime이 진행될수록 thyroid differentiation module은 무너지고, MAPK / invasion-like axis가 상대적으로 올라간다.**

## 산출물 위치

- Tables: `project/results/v17/tables/`
- Figures: `project/results/v17/figs/`
- Figure hub: `project/reports/v17/index.html`
- Methods: `project/reports/v17/v17_methods_detailed.md`
- THYRAI page: `project/reports/html/pages/v17_subtyping.html`

## 현재 한계

1. TERT promoter는 로컬 WXS MAF만으로는 제대로 assay되지 않아, canonical `+TERT` 그룹이 실질적으로 비어 있다.
2. fusion call은 이번 실행에서 공개 원천파일 fetch 대신 **기존 v3 fusion proxy**를 썼다.
3. trajectory는 현재 5-cohort full 성공이 아니라 **3-cohort usable subset** 결과다.
4. DIAL audit도 dark-matter sample과 marker overlap이 충분한 외부 코호트만 남아 **2개 cohort 결과**가 핵심이다.

## 다음 paper narrative 골격

1. **문제 제기**: BRAF/RAS 이분법은 TCGA-THCA의 driver-negative 178명을 설명하지 못한다.
2. **새 관찰**: Dark Matter는 단일 잔여군이 아니라, expression 상 **두 개의 안정적 subtype**로 갈라진다.
3. **강건성**: 이 분리는 external transfer에서 DIA-AUC ~ 0.97-1.00으로 재현된다.
4. **방법론 메시지**: 그러나 ComBat은 biology-preserving transfer와 cohort identifiability 사이에서 긴장을 만든다.
5. **진행축 메시지**: trajectory 상에서 differentiation genes는 pseudotime과 강한 음의 상관을 보이며, 일부 dark-matter cluster는 고위험 tail 쪽으로 사상된다.
