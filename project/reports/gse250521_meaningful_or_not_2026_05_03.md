# GSE250521: 의미 있는 게 뭐냐 — 솔직한 판단
**2026-05-03 · single dataset n=16 slides, 55,873 spots post-QC**

판정: **3개 meaningful + 2개 negative + 1개 unclear**. Paper 1 main figure 로는 부족, **supplementary 한 줄 + Discussion 한 문장** 으로는 가치 있음. Multi-cohort meta 까지 확장하지 말고 marathon 진행 권장.

---

## ✅ Meaningful (3개) — Paper 1 supplementary 에 쓸 가치 있음

### 1. RAI_8 ↔ TDS_like 가 거의 완벽한 negative correlation
- All-spot ρ = −0.89, p ≈ 0
- Sample-mean ρ = −0.87, n=16, p ≈ 1e-5
- **의미**: 8-gene RAI score 가 spatial 환경에서도 thyroid-lineage axis 를 그대로 측정한다는 internal consistency 확인. bulk 에서 정의한 score 가 spatial 에서 깨지지 않음. Paper 1 reviewer Q "Why these 8 genes?" 에 spatial-level 1-line 답변 가능.

### 2. Proliferation gradient 가 stage 와 strong + survives depth correction
- Epithelial subset, sample-mean ρ = +0.80, p = 2e-4 (raw)
- Depth-corrected ρ = +0.58, p = 0.018
- **의미**: 이 dataset 의 stage label (PT→PTC→LPTC→ATC) 이 생물학적으로 일관 → cohort 자체 신뢰 가능. RAI/DM1 신호가 약한 게 stage label noise 때문이 아님 (=cohort 탓이 아님 → 8-gene 이 spatial 에서 약하다는 게 진짜 발견).

### 3. 8-gene 모두 16/16 slides 에서 detect (dropout 0)
- TPO, DIO1, TSHR, PAX8, TG, FOXE1, NKX2-1, SLC5A5 전부 발현
- **의미**: 8-gene panel 이 Visium platform 에서 robust → future ST validation 에서 gene dropout 으로 실패할 우려 없음. Methods Note 한 줄 가치.

---

## ❌ Negative (2개) — Paper 1 main figure 못 만듦

### 4. PT→ATC RAI/DM1 stage trend (sample-mean, all spots)
- RAI_8 ρ = 0.000, p = 1.0
- DM1_like ρ = 0.000, p = 1.0
- TDS_like ρ = −0.11, p = 0.69
- **의미**: 단일 dataset all-spot 차원에서는 8-gene RAI signal 이 stage progression 과 함께 떨어진다는 주장 **증명 불가**. n=16 underpowered 이지만 ρ=0.0 은 power 문제가 아니라 **신호 자체가 약함**.

### 5. Epithelial subset 에서 borderline 신호도 depth correction 후 약화
- Raw: RAI_8 ρ = −0.44, p = 0.091; DM1_like +0.44 p=0.091
- Depth-corrected: RAI_8 ρ = −0.35, p = 0.18; DM1_like +0.35 p=0.18
- **의미**: 살짝 나오던 epithelial trend 의 ~20% magnitude 가 sequencing-depth 인공물 이었음. 방향성은 유지되지만 (RAI ↓, DM1 ↑, TDS ↓) significance 에 한참 못 미침. 단독 claim 불가.

---

## ⚠️ Unclear (1개) — 해석 보류

### 6. DM1_like ~ Hypoxia Simpson flip
- Spot-level ρ = −0.26 (negative)
- Sample-mean ρ = +0.71, p = 0.002 (positive)
- **의미**: aggregate-level 에서는 DM1 high tumor 가 hypoxic 환경과 연관 (생물학적으로 그럴듯), 하지만 spot 내부에서는 반대. Within-sample 해석 vs between-sample 해석이 다름. **단독으로는 못 씀**, multi-cohort 으로 재현되면 의미 있는 mechanism claim 가능. 일단 보류.

---

## 🎯 Paper 1 활용법 (구체적으로)

### Supplementary figure 1장
- `fig_2x4_RAI_DM1_per_stage.png` (2행 × 4 stage representative slides) + violin epi panel
- Caption: "GSE250521 (n=16) Visium spatial validation. RAI_8 ↔ TDS_like ρ=−0.89 confirms 8-gene panel measures thyroid lineage at spatial resolution. Stage trend in epithelial subset is directionally consistent (RAI ↓, DM1 ↑) but underpowered at n=16 (Spearman ρ=−0.35 after depth correction, p=0.18). Multi-cohort meta-analysis required for inferential claim."

### Discussion 한 문장 (본인 voice 영역)
"Single-cohort spatial validation in GSE250521 (n=16 Visium slides) showed direction-consistent but underpowered DM1-axis trend (epithelial-enriched ρ=−0.35 after depth correction); independent confirmation will require multi-cohort spatial meta-analysis."

### Methods 한 줄
"8 RAI genes were detected in 16/16 GSE250521 slides with within-sample z-score; sample-mean Spearman against ordinal stage was used as the primary trend test."

---

## ❓ 확장 여부 판단 (HRA003537 + GSE230424)

**No, 지금 marathon mode 에서는 확장하지 마세요.**

이유:
1. 5/4 marathon Day 1, manuscript writing 6주가 paper-blocking
2. 추가 dataset 2개 ingest = 1주 분석 + 1주 통합 → 마라톤 절반 소진
3. 만약 합쳐서 p<0.05 나와도 Paper 1 main figure 가 아닌 supplementary 한 칸 추가에 그침 (cohort heterogeneity 때문에 reviewer 가 single-figure claim 인정 가능성 낮음)
4. 본 분석 결과 (8-gene 이 spatial 에서 직접 stage signal 못 보여줌) 자체가 reviewer 에게 **약한 결과** — 추가해도 odds 안 좋음

**대신**: 본 supplementary 만 넣고 marathon 진행. Paper 1 reach venue 정당화는 spatial 이 아니라 **driver-excluded sub-stratifier framing + DM1 sub-B NBNR bridge + TERT recovery v2** (이미 STATUS_PAPER1 에 적힌 방향) 으로.

**Pathology POC**: 3,200 tile + DM1~TDS ρ=−0.89 → 회귀 target 명확, ResNet50 baseline 1-2일 작업. 이건 marathon 끝난 후 (6/13 이후) 별도 sprint 로 가능. 지금은 skeleton 만 (`project/src/05_pathology_poc/README.md`) 그대로 두고 진행.

---

## 한 줄 요약

> **Internal consistency 강하고 stage label 신뢰 가능, 그런데 8-gene RAI 가 single-cohort spatial 에서 stage 와 직접 연관성 underpowered. Supplementary figure 1장 + Discussion 한 문장으로 정직하게 쓰고 marathon 진행. 확장 금지.**
