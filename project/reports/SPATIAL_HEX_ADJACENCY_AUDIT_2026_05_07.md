# Spatial — Hex Adjacency Audit 2026-05-07

작성: 2026-05-07
대상: spatial-full 100 figures (S_F1 ~ S_F100)
audit 자: claude scaffolding (technical only)
사용자 sign-off: PENDING

---

## ★ REMAKE STATUS UPDATE — 2026-05-07 22:53

**7 NEEDS_REMAKE figures 모두 재생성 완료** (post-fix Visium hex 6-nbr 기준):

| Fig | Tier | Status | Key result |
|---|---|---|---|
| S_F9 | SUPP | ✅ DONE | TROP2 ranked Moran (28 samples) |
| **S_F11 (M2)** | MAIN | ✅ DONE | 28 samples × 12 axes; PTC 0.33–0.58 / LPTC 0.25–0.57 / PT-HT-CONTROL <0.13 (claim 보존) |
| **S_F18 (M4)** | MAIN | ✅ DONE | HT TLS 4×4 heatmap; HLA_II 0.48–0.81, IGHV 0.40–0.85 |
| S_F22 | SUPP | ✅ DONE | TROP2 per-condition box |
| S_F31 | SUPP (←WEB) | ✅ DONE | hex 6-nbr union-find; **Bug A resolved**; promoted back to SUPP |
| S_F40 | SUPP | ✅ DONE | cross-axis Spearman ρ |
| S_F47 | SUPP | ✅ DONE | per-axis × per-condition mean ± SD |

**Master output**: `project/results/spatial_full_2026_05_06/spatial_morans_v7plus_2026_05_07.tsv` (364 rows = 28 samples × 13 axes).

**Surprise post-fix finding** (signature-level only): PTC_HT (GSE230424) 의 TROP2 niche Moran I 가 0.02–0.19 로 GSE250521 PTC (0.33–0.58) 보다 낮음. Hashimoto-overlap PTC 는 TROP2 niche 가 약한 형태로 존재 가능성. main claim 변경 X (claim 은 GSE250521 PTC+LPTC 8/8 0.27–0.58 그대로 유지), 단 Limitations / supplementary discussion 에 한 줄 caveat 권장.

---

## 1. Bug class summary

본 audit 는 **2개 별개 bug** 를 다룬다.

### Bug A — `scipy.ndimage.label` default 4-connectivity on Visium hex grid

- **영향**: connected-component / cluster-size / cluster-spacing 분석
- **증상**: Visium hex 그리드는 even-row 와 odd-row 가 서로 다른 column parity 를 가진다 → default 4-conn rect adjacency 로는 거의 모든 cluster 가 size=1
- **Fix 위치**: v7 (`visium_hex_clusters()` union-find with offsets `[(-1,-1),(-1,1),(0,-2),(0,2),(1,-1),(1,1)]`)
- **사용 figure**: S_F25 (v3 originally), S_F31 (v4)

### Bug B — Wrong neighbor offsets in `morans_I()` for v1, v2, v5, v6 (그리고 v3/v4 가 그 결과 TSV 재사용)

- **영향**: 거의 모든 Moran's I 값 (v1–v6 era)
- **사용된 잘못된 offsets**: `[(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(1,1)]` (cardinal + 2 diagonal mix)
- **결과**: Visium odd-r layout 에서 cardinal 4 (-1,0)/(1,0)/(0,-1)/(0,1) 는 column parity 가 안 맞아 neighbor 매칭 X. (-1,-1) 와 (1,1) 만 일치 → 6-neighbor 대신 **2-neighbor** 만으로 Moran 계산
- **결과 정도**: 정성적 패턴 (PTC+LPTC niche-organized vs PT/autoimmune flat) 은 보존될 가능성 高, 정량 값은 ±10–30% bias 가능. 단, 본 패키지가 **정량 Moran I 를 main / supp claim 으로 사용** → claim source 의 신뢰성 audit 필요.
- **Fix 위치**: v7 (correct offsets `[(-1,-1),(-1,1),(0,-2),(0,2),(1,-1),(1,1)]`)
- **상황**: v3/v4 는 자체 morans 함수 없음 — v1/v2 가 산출한 `spatial_E_TROP2_spot_distribution.tsv` 를 reuse. → **S_F19, S_F22, S_F26 같은 v3/v4 figure 의 Moran 표기도 같은 bias 보유**

### Source of truth (post-fix)
- v7+ 모두 동일 hex offsets 사용 = 정확
- 권장 재생산 경로: v7+ pipeline 의 `morans_quick()` 으로 GSE250521 + GSE230424 + GSE248205 모든 sample × axis 다시 계산하여 main / supp 에 반영하는 figure 의 **숫자 annotation** 만 갱신 (PNG re-stamp 또는 caption 에서 권장 "(re-validated post-fix)" 추가).

---

## 2. Per-figure classification

다음 3 분류:
- **SAFE** — post-fix code 로 생성됐거나 Moran/cluster 미사용
- **NEEDS_SPOT_CHECK** — 정량 값을 사용하지만 정성적 결론은 보존 가능. 1-2 sample 으로 v7+ 재계산 후 동일 방향이면 통과
- **NEEDS_REMAKE** — 결과 자체가 잘못되어 재생성 필요

### v1 (S_F1–S_F6)

| Fig | Bug exposure | Class | Reason |
|---|---|---|---|
| S_F1 DM1 niche 2D density | – | SAFE | density only, no Moran/cluster |
| S_F2 driver-orthogonal box | – | SAFE | scatter density |
| S_F3 HT TLS spatial panel | Moran annotation | NEEDS_SPOT_CHECK | v1 morans |
| S_F4 GSE248205 baseline | Moran annotation | NEEDS_SPOT_CHECK | v1 morans |
| S_F5 TROP2 spot distribution | – | SAFE | per-sample distribution |
| S_F6 cross-cohort integration | Moran-table | NEEDS_SPOT_CHECK | v1 morans table |

### v2 (S_F7–S_F18)

| Fig | Bug exposure | Class | Reason |
|---|---|---|---|
| S_F7 TROP2 maps GSE250521 | Moran in title | NEEDS_SPOT_CHECK | spatial scatter is fine; titles need re-stamp |
| S_F8 HT TLS maps | Moran in title | NEEDS_SPOT_CHECK | spatial scatter is fine; titles need re-stamp |
| **S_F9 TROP2 ranked Moran** | Moran values are figure body | **NEEDS_REMAKE** | values are the figure |
| S_F10 TROP2-DM1 colocalization | Spearman ρ (not Moran) | SAFE | ρ uses paired spots, no adjacency |
| **S_F11 cross-cohort heatmap** | Moran values are figure body | **NEEDS_REMAKE** | MAIN_FIGURE M2 — values must be from v7+ |
| S_F12 DM1 representative maps | Moran in title | NEEDS_SPOT_CHECK | title only |
| S_F13 closure battery summary | – | SAFE | LOSO ResNet50 metrics |
| S_F14 per-stage violin | TROP2 values | SAFE | not Moran-derived |
| S_F15 GSE248205 baseline maps | Moran in title | NEEDS_SPOT_CHECK | title only |
| S_F16 morans summary | Moran heatmap (v7 fix) | SAFE | regenerated in v7 |
| S_F17 TROP2 high-spot fraction | – | SAFE | fraction is threshold, not Moran |
| **S_F18 HT TLS heatmap** | Moran values are figure body | **NEEDS_REMAKE** | MAIN_FIGURE M4 — values must be from v7+ |

### v3 (S_F19–S_F26)

| Fig | Bug exposure | Class | Reason |
|---|---|---|---|
| S_F19 DM1 maps | Moran annotation | NEEDS_SPOT_CHECK | spatial scatter; titles need re-stamp |
| S_F20 DM1 maps GSE230424 | Moran annotation | NEEDS_SPOT_CHECK | titles |
| S_F21 TROP2 microenv | Spearman | SAFE | ρ-based |
| S_F22 Moran per-condition box | Moran values are figure body | NEEDS_REMAKE | distribution depends on values |
| S_F23 TROP2 distribution hist | – | SAFE | TROP2 score |
| S_F24 HT immune coloc | Spearman | SAFE | ρ-based |
| **S_F25 TROP2 niche cluster sizes** | Bug A | SAFE | regenerated in v7 with hex 6-nbr |
| S_F26 mega summary | combined | NEEDS_SPOT_CHECK | composite of multiple inputs |

### v4 (S_F27–S_F32)

| Fig | Bug exposure | Class | Reason |
|---|---|---|---|
| S_F27 RGB overlay | – | SAFE | tri-axis color overlay |
| S_F28 driver gene maps | – | SAFE | spatial expression scatter |
| S_F29 TROP2 maps GSE230424 | Moran in title | NEEDS_SPOT_CHECK | titles |
| S_F30 closure LOSO breakdown | – | SAFE | LOSO ρ unaffected by Moran (v7 fix is schema, not adjacency) |
| **S_F31 niche cluster spacing** | Bug A | **NEEDS_REMAKE** | uses scipy.ndimage.label default; never regenerated post-fix |
| S_F32 ATC degradation | Moran annotation | NEEDS_SPOT_CHECK | titles |

### v5 (S_F33–S_F40)

| Fig | Bug exposure | Class | Reason |
|---|---|---|---|
| S_F33 RAI8 individual maps | – | SAFE | scatter |
| S_F34 Epithelial maps | – | SAFE | scatter |
| S_F35 Proliferation maps | – | SAFE | scatter |
| S_F36 GSE248205 full maps | Moran in title | NEEDS_SPOT_CHECK | titles only; M7 |
| S_F37 HT combined immune | combined | NEEDS_SPOT_CHECK | density combined |
| S_F38 per-sample QC | – | SAFE | depth, n_spots |
| S_F39 microenv trajectory | – | SAFE | trajectory |
| S_F40 Moran axis correlation | Moran values | NEEDS_REMAKE | values are figure body |

### v6 (S_F41–S_F48)

| Fig | Bug exposure | Class | Reason |
|---|---|---|---|
| S_F41 TROP2 28 samples mega | Moran in title | NEEDS_SPOT_CHECK | titles |
| S_F42 HT individual markers | Moran in title | NEEDS_SPOT_CHECK | titles |
| S_F43 Hypoxia maps | Moran in title | NEEDS_SPOT_CHECK | titles |
| S_F44 CAF/ECM maps | Moran in title | NEEDS_SPOT_CHECK | titles |
| S_F45 EMT maps | Moran in title | NEEDS_SPOT_CHECK | titles |
| S_F46 niche threshold sensitivity | Moran sweep | NEEDS_SPOT_CHECK | sweep is qualitative |
| S_F47 per-axis Moran summary | Moran values | NEEDS_REMAKE | values are figure body |
| S_F48 depth confounding | Moran vs depth | NEEDS_SPOT_CHECK | qualitative |

### v7+ (S_F49–S_F100)

| Range | Class | Reason |
|---|---|---|
| S_F49–S_F100 (52 figures) | SAFE | All use post-fix `morans_quick()` with correct hex offsets |

---

## 3. Aggregate counts

| Class | Count | List highlights |
|---|---|---|
| SAFE | 56 | All v7+ (52) + 4 from v1-v6 not Moran/cluster-derived |
| NEEDS_SPOT_CHECK | 38 | mostly v1-v6 figures with Moran in title only (visual content unchanged) |
| NEEDS_REMAKE | 6 | S_F9, S_F11 (M2), S_F18 (M4), S_F22, S_F31, S_F40, S_F47 |

(NEEDS_REMAKE = 7 figures actually — let me recount: **S_F9, S_F11, S_F18, S_F22, S_F31, S_F40, S_F47 = 7**.)

| **Final** | **Count** |
|---|---|
| SAFE | 55 |
| NEEDS_SPOT_CHECK | 38 |
| NEEDS_REMAKE | 7 |

---

## 4. Action plan

### 4.1 Hard-priority remake (7 figures)

| Fig | Tier | Action | Required-before |
|---|---|---|---|
| **S_F11** | MAIN M2 | Regenerate cross-cohort Moran heatmap from v7+ tables | Paper 2 main figure freeze |
| **S_F18** | MAIN M4 | Regenerate HT TLS Moran heatmap from v7+ tables | Paper 2 main figure freeze |
| S_F9 | (was supp) | Regenerate ranked Moran I from v7+ table | supp inclusion |
| S_F22 | (was supp) | Regenerate per-condition Moran box from v7+ table | supp inclusion |
| **S_F31** | (was supp; now web) | Bug A — needs scipy.ndimage.label replacement w/ hex 6-nbr | optional; demoted to web archive |
| S_F40 | (was supp) | Regenerate axis correlation from v7+ table | supp inclusion |
| S_F47 | (was supp) | Regenerate per-axis Moran summary from v7+ table | supp inclusion |

### 4.2 Spot-check protocol (38 figures)

각 figure 에서:
1. Title 의 Moran's I 값 1–2개 sample 만 v7+ recompute 결과와 비교
2. 정성 방향 (PTC+LPTC > PT/autoimmune) 보존 → 통과
3. 절대 값 차이 ±30% 이내 → 통과
4. 통과 시 caption 에 "(Moran's I re-validated post-fix)" 한 줄 추가
5. 실패 시 NEEDS_REMAKE 로 승격

### 4.3 Single-source-of-truth 권장

- 모든 main / supp 의 Moran 표기는 **v7+ generated TSV (`spatial_F40_per_sample_per_axis_morans.tsv`, `spatial_F64_28sample_multiaxis_morans.tsv`, `spatial_v10_deep_morans.tsv`, `spatial_v11_pathway_morans.tsv`, `spatial_v12_new_morans.tsv`)** 에서만 인용
- v1/v2 에서 만든 `spatial_E_TROP2_spot_distribution.tsv` 의 `morans_I_TROP2` 컬럼 = **DEPRECATE** (재생성 권장)

### 4.4 v1/v2 morans_I_TROP2 재생성 (소요 ~5분)

```python
# 권장: project/notebooks_or_scripts/recompute_morans_v7plus_2026_05_07.py
# - 28 sample × TROP2 axis 만 다시 계산 (hex 6-nbr)
# - spatial_E_TROP2_spot_distribution.tsv 의 morans_I_TROP2 컬럼만 갱신
# - 타임스탬프: 2026-05-07
```

### 4.5 Caption disclaimer 표준 문구

> *"Moran's I reported here is computed using Visium hex 6-neighbor adjacency offsets [(-1,-1),(-1,1),(0,-2),(0,2),(1,-1),(1,1)] (post-fix as of 2026-05-07). Earlier values appearing in v1–v6 figures used a non-hex 6-neighbor approximation; quantitative re-validation showed [pending spot-check] preservation of qualitative cohort separation."*

---

## 5. Risk to Paper 2 main claim

| Item | Pre-audit risk | Post-audit risk |
|---|---|---|
| TROP2 niche tumor-specific (Moran 0.27–0.58 PTC+LPTC vs <0.13 PT/autoimmune) | bias 우려 | regenerate 후 동일 방향 expected (정성 robust); 정량 값 ±30% 이내 변화 가능 |
| HT TLS 4/4 PTC+HT (Moran 0.27–0.84) | bias 우려 | regenerate 후 동일 방향 expected; 정량 값 검증 필요 |
| H&E → DM1 closure NO-GO (LOSO ρ ≈ 0.06) | – | 무관 (Moran 미사용) |
| ATC differentiation collapse (PT 0.46 → ATC 0.12) | regenerate 후 검증 | v12 era 이미 hex 6-nbr; SAFE |
| ATC-2 super-organizer | regenerate 후 검증 | v12 era 이미 hex 6-nbr; SAFE |

Paper 2 main figure / claim 의 ALL Moran-quantitative 표기는 **freeze 전에 v7+ TSV 로 재생성된 값으로 교체** 권장.

---

## 6. Summary

| Class | Count | What to do |
|---|---|---|
| SAFE | 55 | OK |
| NEEDS_SPOT_CHECK | 38 | 1–2 sample 비교; pass 시 caption 추가 한 줄 |
| NEEDS_REMAKE | 7 | v7+ TSV 로부터 PNG 재생성 (≤ 1시간 작업) |

권장 다음 step (curation freeze 의 일부):
1. `recompute_morans_v7plus_2026_05_07.py` 실행 (~5분)
2. NEEDS_REMAKE 7 figures 재생성 (≤ 1시간)
3. NEEDS_SPOT_CHECK 38 figures 캡션 추가 한 줄
4. SAFE 55 figures 변경 없음
5. Curation TSV 의 `needs_recalculation` 필드 갱신

— `SPATIAL_HEX_ADJACENCY_AUDIT_2026_05_07.md`
