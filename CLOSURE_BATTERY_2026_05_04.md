# Closure battery — H&E → DM1_like_score_resid (final decision: D 폐기)

**Run:** 2026-05-04 · **Owner:** Seungho Cook
**Wall time:** ~85 min CPU · No usable GPU (Azure quota=0, RunPod $30 funded but pod stuck pulling image)

---

## 한 줄 결론

> **D — 원본 image-DM1 idea 완전 폐기.** ResNet50 ImageNet 224 px tile에서 depth-residualized DM1_like axis 신호가 random panel과 구분 안 됨 (real 0.022 ≈ random max 0.031). 5가지 실패 원인 중 (3) "depth-residualization이 morphology 신호 제거"는 부분적 진실이지만 그 신호 자체가 stain↔depth artifact라 residualization은 정당. (4) "spot-level ST label과 H&E tile content가 fundamentally mismatched"가 dominant cause (stage ordinal r=−0.05 같이 죽음, ATC만 AUROC 0.69 — anaplastic gross histology이라 자명). Foundation model 업그레이드 (~+0.05~+0.15 gain 기대)도 GO threshold (+0.28) 못 건넘. Paper 2/IP에서 image-DM1 angle drop, manuscript marathon으로 pivot.

## 의사결정 매트릭스

| 옵션 | 통과 기준 | 결과 |
|---|---|---|
| **A. 원본 image-DM1 유지** | DM1_resid Spearman ≥ 0.20 OR AUROC ≥ 0.65 + 컨트롤 이김 + color/depth artifact 아님 | ✗ best 0.080 (simple RGB), 모든 setting 미달 |
| **B. raw thyroid-lineage morphology pivot** | raw Spearman ≥ 0.30 OR AUROC ≥ 0.70 + resid 약함 + 컨트롤 약함 | ✗ raw best 0.224/0.702인데 stain↔depth artifact, spec상 무효 |
| **C. gross histology / stage 모델 pivot** | stage 또는 ATC 강함 + molecular 약함 | △ ATC AUROC 0.689 (BORDERLINE) — 그러나 anaplastic 진단은 prior art 풍부 |
| **D. full 폐기** | 모든 molecular r < 0.15 AND AUROC < 0.60 + stage 약함 + 컨트롤 비슷 | ✓ molecular ≤0.080, stage ordinal r=−0.05, ATC만 trivial |
| **E. foundation model retry** | ResNet raw/stage 신호 OR DM1_resid r ≥ 0.10 | ✗ raw 0.06, resid 0.022 — spec "otherwise no" |

→ **D 채택. E는 향후 옵션으로만 (HF UNI access 신청 + RunPod 재가동 시)**.

## 가장 중요한 숫자 5개

1. DM1_like_score_resid Ridge 224: **r=0.022 / AUROC=0.511** (random max 0.031, indistinguishable)
2. DM1_like_score_raw + 9-dim RGB feature: **r=0.224 / AUROC=0.702** ← 2048-dim ResNet50 (0.061)을 압도. 신호가 H&E 색상 artifact임을 증명
3. RAI_8 raw real (50 random panel 비교): **0.067 ≈ random max 0.062** — raw도 ResNet50 space에선 random과 구분 안 됨
4. ATC vs non-ATC LogReg AUROC: **0.689** (BORDERLINE, 그러나 anaplastic 진단은 trivial prior art)
5. stage ordinal Ridge: **r=−0.047** (chance) — 단일 tile에서 stage gradient 안 보임

## 실험별 요약 (closure battery 6개 + 2 neg ctrl 모드)

| 섹션 | 실험 | best metric | 해석 |
|---|---|---|---|
| Phase A | DM1_resid 224 ResNet50 LOSO | r=0.022 | NO-GO baseline |
| **A** | raw DM1/RAI/TDS sensitivity | r=0.061 (ENet 음수 0.09) | raw도 ResNet50 space 약함 |
| **B** | resid variants (5종) | best B4 rank-norm r=0.060 | residualization formula 문제 아님 |
| **C** | tile size 448 | DM1_resid r=**0.017** (vs 224 0.022, delta −0.005) | hypothesis (2) **REJECTED** — 더 크게 만들어도 안 좋아짐 |
| **C** | tile size 672 | _cancelled_ (spec gate +0.05 미달) | n=2,386 tile만 extracted, embed 안 함 |
| **D** | stage / morphology sanity | ATC AUROC 0.689; stage r=−0.047 | 모델은 ATC 가능, gradient 불가능 |
| **E** | slide-level aggregation (n=16) | DM1_resid top25_mean r=+0.209 | n 작아 marginal, 신호 부재 |
| **F** | simple RGB 9-dim baseline | DM1_raw r=0.224 / AUROC 0.702 | ResNet50보다 잘 — color artifact 증명 |
| **G** | HF UNI/CONCH/Virchow2 access | 없음 (token X, cache X) | 향후 옵션 |
| neg_ctrl resid | 30 random + housekeeping | random mean 0.002, max 0.031 | real (0.022) ≤ max — 노이즈 |
| neg_ctrl raw | 50 random + housekeeping + RAI_real | random mean 0.019, max 0.062; RAI real 0.067 | real ≈ max — ResNet50 raw도 노이즈 |

## 5가지 실패 원인 — 증거 기반 ranking

| 원인 | 기여도 | 증거 |
|---|---|---|
| (1) ResNet50/ImageNet 약함 | 부분 | 9-dim RGB가 raw에서 압도 (0.224 vs 0.061) — ResNet50 ImageNet feature는 이 task에 mismatched |
| (2) 224 px tile 너무 작음 | **REJECTED** | 448 r=0.017 (224의 0.022보다 약간 악화). Tile context 키워도 신호 안 살아남 — 224 작아서가 dominant 아님 |
| (3) depth-residualization이 morphology 제거 | 부분 (정당함) | resid 0.022 ↔ raw 0.06: residualization이 신호의 65%를 제거. 그러나 그 신호는 stain↔depth artifact였음 (F section) → residualization은 옳은 design |
| (4) spot-level ST label이 H&E와 mismatched | **dominant** | (a) stage ordinal r=−0.05; (b) epi/prolif raw r≈0.03; (c) DM1_resid 모든 setting 약함; (d) ATC만 됨 — gross 만 |
| (5) hires image 해상도 / single cohort | 부분 | hires는 fullres의 0.17× downsample, nuclear morphometry 없음. n=16 cohort 1개 |

→ Foundation 모델은 (1)을 일부 풀 수 있으나 (3)+(4)는 못 풂. 따라서 retry value 낮음.

## Paper 2 / IP 영향

- ✅ **Pillar I/II/III 영향 없음** (Korean PTC vs Korean baseline, HT-overlap, Yu 2026-05-04)
- ❌ **image-DM1 / pathology-AI angle Paper 2에서 drop**
- 📝 **본 NO-GO를 Methods supplementary**로 — pre-registered failed hypothesis로 molecular-only 논리 강화
- 💰 **IP 무가치** (color-artifact는 patentability X, ATC classifier prior art 풍부)
- ⏱ **4-6주 절약** (multi-cohort + TCGA WSI download path 제거)

## 다음 액션

```bash
# 1. scp 결과 본인 PC로
scp seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/CLOSURE_BATTERY_2026_05_04.md .
scp seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/project/reports/pathology_dm1_closure_battery_2026_05_04.md .
scp -r seungho@40.82.129.113:/home/seungho/personal/THCA_data_analysis/project/results/03_pathology_poc .

# 2. Paper 2 brief에서 image-DM1 angle 삭제

# 3. 향후 옵션 (block 안 됨):
#    HF UNI access 신청: https://huggingface.co/MahmoodLab/UNI
#    승인 후 RunPod에 swap 가능 (현재 $30 잔액 살아있음)

# 4. Manuscript marathon 5/4-6/13 재개 (Pillar I forest)
```

## 산출물 (모든 파일)

```
project/results/03_pathology_poc/
  embeddings_resnet50_224.npz                    11.6 MB
  embeddings_resnet50_448.npz                    9.9 MB     (LOSO 진행중)
  embeddings_resnet50_672.npz                    예정       (CPU 진행중)
  loso_metrics_resnet50.tsv                      Phase A
  loso_predictions_resnet50.tsv.gz               per-tile preds
  closure_battery_metrics.tsv                    26 rows (A/B/D/E/F)
  closure_battery_predictions.tsv.gz             per-tile preds
  closure_battery_summary.png                    horizontal bar plot
  negative_controls_summary.tsv                  30 random + HK (resid)
  negative_controls_raw_summary.tsv              50 random + HK + RAI (raw)
  pred_vs_obs_resnet50.png                       224 scatter
  pred_vs_obs_resnet50_per_slide.png             16-slide bar
  simple_features.npy                            9-dim RGB stats
  tile_metadata*.tsv.gz                          8,521 tiles 224/448/672

project/reports/
  pathology_dm1_phaseA_cpu_verdict_2026_05_04.md       Phase A NO-GO 상세
  pathology_dm1_closure_battery_2026_05_04.md          closure 상세

CLOSURE_BATTERY_2026_05_04.md                          ← 이 파일 (top-level summary)
PHASE_A_NOGO_REPORT_2026_05_04.md                      Phase A summary
PATHOLOGY_DM1_FEASIBILITY_2026_05_03.md                초기 plan
```

---

*Total compute: ~85 min CPU on D8as_v5 (no GPU). All experiments under spec rules: no random tile split, no raw-success claim, no Phase B/C, no TCGA WSI download, no RET fusion claim. RunPod $30 funded but A4000 pods stuck at image pull (community + secure both tried).*
