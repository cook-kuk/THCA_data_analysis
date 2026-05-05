---
title: "HLA work — current state + strengthening options (situational audit)"
date: 2026-05-04
purpose: HLA 쪽 연구 강화 의향에 대한 현 상황 정리. Marathon-compliant — spec/audit only, no execution.
url_referenced: http://40.82.129.113/papers_hub_2026_05_04 — fetch timeout (외부 환경 unreachable). 본인 local network 에서만 접근 가능 예상. 본인 직접 read 권장.
scope: cross-paper situational (Paper 2 HT-overlap + Paper 4 GD backlog + cross-paper boundary)
status: AUDIT ONLY — no analysis execution, no file edits, no manuscript prose
---

# HLA work — current state + strengthening options

---

## 1. 현재 HLA 작업 state (what's done)

### 1.1 Existing result directories

| Directory | Content | Status |
|---|---|---|
| `project/results/d4p1_panasian_meta/` | **v1 forest meta (DEPRECATED)** — Korean PTC pool n=908 vs Chu 2018 GD direct comparison | ❌ **DEPRECATED** (2026-05-04, conflation 1) |
| `project/results/p2_pillar1_forest_v2/` | **v2 forest meta (4-allele Korean PTC vs Korean baseline)** — Korean PTC pool n=874 (K2 235 + Lee 630 + GSE286332-PTC 9) | ✅ Already executed; 2 alleles (C*01:02, DQB1*02:01) deferred for Lee 2014 fill |
| `project/results/v17_hla/` | TCGA HLA per-sample + Korean GSE213647 HLA per-sample + lymph annotations | ✅ Available (Paper 2 substrate) |
| `project/results/v17_hla_autoimmune/` | HA1 carrier freq vs population, HA2 lymph per allele, HA3 DM crosstab, HA4 Asian vs White, HA5 Cox summary | ✅ Available (Paper 2 substrate, may need HT-isolation review) |
| `project/results/p5_8gene_vs_hla_autocorr/` | 8-gene vs HLA-II module autocorrelation (Paper 1 audit) | ✅ Used for Paper 1 §2.4 residualization |
| `project/results/d3p5_pdm1_gradient/` | Mediation analysis (HLA-II 140% mediation) | ✅ Available (Paper 2 mechanism) |

### 1.2 Cohort inventory (HLA-typeable)

| Cohort | n typeable | Imputation method | Disease substrate | Paper |
|---|---|---|---|---|
| K2 / PRJEB11591 (Yoo 2016) | 235 | arcasHLA (RNA-seq) | Korean PTC, ~20-30% HT comorbidity baseline | Paper 2 |
| Lee 2024 / GSE213647 | 630 | arcasHLA (RNA-seq) | Korean PTC, Hashimoto-like 22-28% (Otsu) | Paper 2 |
| GSE286332 PTC arm | 9 | arcasHLA (RNA-seq) | Korean PTC (NOT PTC+HT) | Paper 2 |
| GSE286332 PTC+HT arm | 9 | arcasHLA (RNA-seq) | Korean PTC+HT | Paper 2 (Pillar II/III mediation) |
| TCGA-THCA | n~500 | OptiType / arcasHLA from RNA-seq | Mixed ancestry primary PTC | Paper 1 + Paper 2 (cross-cohort) |
| AFND South Korea baseline | TBD pooled | Reference (KOTRY + KBP + cookHLA + Lee 2014) | Korean general population | Paper 2 baseline |
| Chu 2018 Han Chinese | 1,468 GD + 1,490 ctrl | Published summary stats only | Chinese GD | **Paper 4 reserved** |
| Bundang Graves' BTC cohort | 0 (outreach pending) | TBD | Korean GD | **Paper 4 backlog gated** |

### 1.3 Locked decisions (5/4 sweep)

- ✅ Pillar I v2 spec frozen (PAPER2_PILLAR1_V2_FOREST_META_SPEC_2026_05_04.md)
- ✅ G3 = YES Lee 2014 fill (6-allele forest path)
- ✅ Paper 2 = Hashimoto-overlap PTC ONLY (memory `v18_paper2_HT_isolated`)
- ✅ Paper 4 = GD backlog, gated 4/4 (Paper 1 bioRxiv + Paper 2 A/B/C + Bundang n>50 + Yu 합의)
- ✅ "autoimmune-PTC" → "Hashimoto-overlap PTC" terminology cleanup (8 Paper 1 manuscript files done)
- ✅ Forest v1 vs Chu GD direct comparison DEPRECATED

### 1.4 Open decision gates (Pillar I v2 spec §8)

| Gate | Status |
|---|---|
| G1 advisor (Yu professor) approval of v2 spec | ⏸ pending |
| G2 AFND South Korea baseline pool composition | ⏸ pending |
| G3 Lee 2014 fill | ✅ YES locked 5/4 |
| G4 Harbin Korean fallback exclusion (sensitivity) | ⏸ pending |
| G5 Figure 1 re-render with v2 forest | ⏸ pending |
| G6 manuscript paragraph framing | ⏸ pending |

---

## 2. Cross-paper boundary (HLA work split)

### Paper 2 (Hashimoto-overlap PTC) — HLA scope
- **Pillar I**: Korean PTC HLA susceptibility (DPB1*05:01 + 5 other alleles, vs Korean baseline)
- **Pillar II/III**: GSE286332 PTC+HT vs PTC mediation (HLA-II module 140% mediation per d3p5_pdm1_gradient)
- **HLA imputation**: arcasHLA from RNA-seq (n=874 Korean PTC pool)
- **Mechanism focus**: antigen-driven HLA-II–mediated dedifferentiation
- **Allowed terms**: HLA-II, antigen presentation, Korean PTC vs Korean baseline, AFND South Korea baseline

### Paper 4 (Korean GD HLA / Pan-Asian) — BACKLOG
- **Substrate**: Bundang Graves' BTC cohort (when arrives) + Chu 2018 Chinese GD published
- **Reserved cohorts**: Chu 2018 (1,468 GD + 1,490 ctrl)
- **Reserved alleles**: GD-specific (e.g., TSH receptor associated)
- **Gating**: 4/4 must clear (Paper 1 bioRxiv + Paper 2 A/B/C + Bundang n>50 + Yu 합의)
- **Forbidden in Paper 2 context**: Graves, GD, TSAb, TSI, hyperthyroidism, etc.

### Cross-paper Discussion-level coordination
- DPB1*05:01 = shared East-Asian thyroid autoimmunity allele (HT + GD overlap)
- Paper 2 Discussion: 1-line GD reference background only ("shared HLA background with GD; disease distinct — see Paper 4")
- Paper 4 Discussion (future): same 1-line reference back to Paper 2

---

## 3. HLA strengthening options (categorized + Marathon-scored)

### Category A — Paper 2 Pillar I (HT-only) execution path

| ID | Option | Marathon safety | Notes |
|---|---|---|---|
| **A1** | Execute Pillar I v2 6-allele forest (after G2 + advisor) | ★★★ | Already spec'd; awaits G1 + G2 only |
| **A2** | Lee 2014 baseline frequency extraction (C*01:02 + DQB1*02:01) | ★★★ | Single literature lookup, no analysis |
| **A3** | AFND South Korea baseline pool composition decision (G2) | ★★★ | Decision only |
| **A4** | Sensitivity analysis (Harbin Korean fallback exclusion) | ★★★ | Spec extension, no new data |
| **A5** | Per-sub-cohort forest (K2 / Lee / GSE286332 separate) | ★★★ | Spec extension within Pillar I v2 |

### Category B — Paper 2 mechanism deepening (HLA-II axis)

| ID | Option | Marathon safety | Notes |
|---|---|---|---|
| **B1** | Re-run mediation analysis with HT-isolated substrate (currently d3p5_pdm1_gradient may include autoimmune-PTC framing) | ★★ | Re-execution risk; spec-level review first |
| **B2** | HLA-II module activation pathway analysis (GSE286332 PTC+HT vs PTC) | ★★ | New analysis, but builds on existing P3 work |
| **B3** | Single-cell HLA-II expression in Pu 2021 / Lu 2023 thyrocyte clusters | ★★ | sc external data already loaded, spec extension |
| **B4** | BCR/TLS × HLA-II crosstab (D5-P6 BCR clonal + TLS work integration) | ★★ | Spec OK, execution = Paper 2 territory |
| **B5** | Antigen presentation gene set scoring (HLA-A/B/C + class II genes) | ★★ | Module score, marathon-safe spec |

### Category C — HLA imputation methodology strengthening

| ID | Option | Marathon safety | Notes |
|---|---|---|---|
| **C1** | cookHLA SNP integration (currently deferred per `v17_arcasHLA_korean_k2`) | ★ | Requires SNP data acquisition (currently absent for K2/Lee). Major data work. |
| **C2** | arcasHLA confidence cutoff sensitivity (n=235 typeable from K2 n=260; verify dropouts) | ★★★ | QC review, no new analysis |
| **C3** | Cross-validation arcasHLA vs HIBAG vs OptiType (TCGA only — has both) | ★★ | Spec only; execution = method paper material |
| **C4** | TCGA HLA imputation method consistency check | ★★★ | Existing data audit |

### Category D — Pathology image multimodal layer (user expertise)

| ID | Option | Marathon safety | Notes |
|---|---|---|---|
| **D1** | TCGA WSI HLA-II IHC pattern integration | ★ | New multi-modal analysis; requires image + HLA join. Major new direction. **Paper 1 strategy: H&E-DM1 dropped** — pathology image scope must be carefully scoped to Paper 2 HLA-II IHC, not DM1 inference. |
| **D2** | Pathology image multimodal scoping memo | ★★★ | Spec only, defines what image layer means for Paper 2 |

### Category E — Cross-paper coordination (Paper 2 + Paper 4 boundary)

| ID | Option | Marathon safety | Notes |
|---|---|---|---|
| **E1** | Discussion-level 1-line GD reference (Paper 2 §3.x) | ★★★ | Single sentence scaffolding, voice-protected by 본인 |
| **E2** | Paper 2 Discussion + Paper 4 Discussion mutual cross-link spec | ★★★ | Spec only, voice = 본인 |
| **E3** | DPB1*05:01 shared East-Asian thyroid autoimmunity allele framing | ★★★ | Single paragraph spec, no execution |

### Category F — Paper 4 backlog (BLOCKED)

| ID | Option | Marathon safety | Notes |
|---|---|---|---|
| F1 | Chu 2018 Chinese GD vs Chinese baseline forest | ❌ **BLOCKED** | Paper 4 territory, 4/4 gating not met |
| F2 | Bundang Graves' cohort outreach follow-up | ❌ **BLOCKED** | Outreach work, not analysis. Paper 4 territory |
| F3 | Okada 2015 Japanese / Taiwan CMUH GD HLA literature search | ❌ **BLOCKED** | Paper 4 territory |
| F4 | Cross-population GD allele OR meta-analysis | ❌ **BLOCKED** | Paper 4 territory |

**범례**:
- ★★★ = full marathon-safe (spec/decision/audit, no execution)
- ★★ = marathon-safe with discipline (spec OK, execution after specific approval)
- ★ = major new direction (data acquisition or new analysis required)
- ❌ = BLOCKED (gating not met)

---

## 4. 강화 우선순위 (Marathon-compliant)

본인 의향 "HLA 쪽 연구 더 강화" 에 대한 **marathon-safe + paper-blocking 우선순위**:

### Tier 1 — 즉시 실행 가능 (advisor 결정만 필요)

1. **A3 G2 AFND baseline pool composition decision** — Yu professor 1-question 답변
2. **A1 Pillar I v2 6-allele forest execution** — G1 + G2 받으면 즉시 실행 (1-2 hr 분석)
3. **A2 Lee 2014 baseline frequency lookup** — 30 min 문헌 lookup + Paperclip cite verify

### Tier 2 — 마라톤 spec 작성 (실행 보류)

4. **B5 Antigen presentation gene set scoring spec** — Paper 2 mechanism 보강
5. **B3 sc HLA-II thyrocyte cluster analysis spec** — Pu 2021 / Lu 2023 활용
6. **B4 BCR/TLS × HLA-II crosstab spec** — 기존 D5-P6 work 통합
7. **C4 TCGA HLA imputation consistency audit spec** — QC 강화

### Tier 3 — 큰 결정 필요 (advisor + 시간 + 데이터)

8. **C1 cookHLA SNP integration** — SNP data 확보 필요 (Korean cohort access)
9. **D1 Pathology image multimodal** — Paper 1 strategy 와 boundary 검토 필요 (H&E-DM1 dropped 이지만 HLA-II IHC 는 다른 layer)

### Tier 4 — Paper 4 backlog (BLOCKED)

10. F1-F4 — 4/4 gating 없으면 진행 X

---

## 5. 본인 specific 결정 요청

본인이 HLA 강화 의향에서 어느 방향 우선:

- **(α)** Tier 1 only — Pillar I v2 execution 마무리 (A1+A2+A3, advisor 1-question + 6-allele forest)
- **(β)** Tier 1 + Tier 2 spec — Pillar I 마무리 + Paper 2 mechanism 보강 spec (B5+B3+B4)
- **(γ)** Tier 2 spec only — execution 미루고 spec 위주 (Marathon-safest)
- **(δ)** D2 pathology image scoping memo — multimodal layer Paper 2 boundary 정의
- **(ε)** E3 cross-paper coordination spec — Paper 2/4 Discussion 1-line 시너지
- **(ζ)** 다른 조합 또는 (specify)

**본인 specific 선택 받기 전 까지 — execution 0. spec/plan 작성도 specific approval 후 시작.**

---

## 6. URL note (본인 직접 read 필요)

`http://40.82.129.113/papers_hub_2026_05_04` — 외부 환경 fetch timeout. 본인 local network (Azure burst Korea Central?) 에서만 접근 가능 예상. 본인이 직접 열어서 다음 확인:

1. Papers hub 의 HLA 관련 section 이 있는지
2. 위 audit 와 일치 / 불일치
3. 추가로 본 audit 에 반영할 내용 (강화 후보 / 우선순위 / 본인 brand 정렬)

본인 read 후 결과 공유하시면 본 audit md 에 통합 가능.

---

## 7. Marathon mode reaffirmation

- 5/4–6/13 active
- 새 분석 실행 = paper-blocking 만 (Tier 1 = paper-blocking, Tier 3 ★ 1개는 Paper 1 strategy boundary 검토 필요)
- 새 데이터 다운로드 금지 (C1 cookHLA SNP = data acquisition gate)
- voice-protected prose generation 금지 (E1 Discussion 1-line = 본인 voice)
- spec/audit only without specific approval

# END SITUATION AUDIT — execution 0. 본인 specific 선택 받기 전 까지 대기.
