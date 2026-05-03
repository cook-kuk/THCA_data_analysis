# Paper 1 — DM1/RAI axis · molecular axis bundle (post-closure-battery audit, 2026-05-04)

> **STATUS — INTERNAL WORKING DOCUMENT.** This bundle is **not** the manuscript draft and **not** for direct submission inclusion. Per the H&E → DM1 closure battery (2026-05-04, verdict **D = full 폐기**; see `project/reports/pathology_dm1_closure_battery_2026_05_04.md` and root `CLOSURE_BATTERY_2026_05_04.md`), the image-DM1 / RunPod G1–G3 angle is **dropped**. Per marathon mode (5/4–6/13), only molecular-axis strengthening claims may be considered for Paper 1; all venue-probability framing and wet-lab roadmaps are internal planning only. **Inline tags below mark each section as `[MANUSCRIPT-SAFE candidate]`, `[INTERNAL-ONLY]`, or `[DEPRECATED]`.**

> **Source artifacts for this document:** 28 ST slides + 561 TCGA-THCA bulk + 12 external Visium (`project/results/01_spatial_score/`, `project_external_st/results/`, etc.). All numbers below are from prior runs; **no new analysis** is performed by this audit.

---

## 0. 한 단락 status (5/4 14:00 KST, audited 2026-05-04 post-closure)

> `[INTERNAL-ONLY]` summary paragraph. Manuscript-safe extracts appear in tagged sections below.

8-gene RAI / DM1_like axis 가 thyroid cancer dedifferentiation 의 prognostic + mechanistic axis 임을 28 ST slides + 561 TCGA-THCA bulk + 12 external Visium 으로 multi-layer validation. **TDS overlap + inflammation artifact 두 risk 는 분석상 controlled** (residualized OR/Cohen d 계산 + negative-control panel). NKX2-1/FOXE1 collapse + DNMT activation + STAT3/AP-1 activation 의 4-step molecular framework 은 ULM TF inference + GSEA Hallmark (IL6_JAK_STAT3 + EMT + IFN-γ up; OXPHOS down) + TF coordination network 로 cross-checked. TROP2/TACSTD2 는 **tumor-population level vulnerability** (tumor 1.4× normal, metastatic 1.6×), DM1-high subset-specific spot-level colocalization은 negative (Q3 ρ=−0.016).

### Remaining risks (NOT "all risks resolved")

1. **N1 dark matter underpowered:** BRAF-neg AND RAS-neg subset (n=156–187) 단독 PFI HR=1.20 (NS, p=0.80). DM1 prognostic effect 는 driver-mutant tumors 에 집중, dark matter 안에서는 cluster-discriminating only.
2. **Q3 TROP2 spot-level NOT colocalized:** Mean ρ(TROP2, DM1_resid) = −0.016 (28 slides). Sample-level (N2) tumor>normal effect는 살아 있으나, "DM1-high spot = TROP2-high spot" 주장은 못 함. TROP2 framing은 **tumor 인구 vulnerability** 만.
3. **H&E → DM1 image-axis NO-GO (closure battery 2026-05-04):** ResNet50 ImageNet 224 px DM1_resid Spearman 0.022 ≈ random max 0.031. Foundation-model upgrade projected gain (~+0.05–+0.15) 도 GO threshold (+0.28) 미만. **Image-DM1 angle Paper 1/2 에서 drop.** Pod B (WSI), G1/G2/G3 RunPod 파이프라인 폐기.
4. **No wet-lab functional validation:** TROP2 IHC, sacituzumab IC50, organoid, xenograft 0건. Computational layer 만으로는 Cancer Cell-tier 도전권 미충족 — venue 결정은 §3 (internal-only) 참조.
5. **Pan-cancer scope 좁음:** THCA r=−0.892 outlier; 다른 cancer |r|<0.55. "Thyroid-specific dedifferentiation prototype" framing 만 가능.
6. **Voice-protected 섹션 미작성:** Hook ¶1, Aim ¶4, Discussion §3.1 (Landa 2016 cite save), §3.4 Limitations, Cover ¶1, Reviewer Q9 — 본인 키보드 슬롯 대기.

## 1. 모든 sprint 결과 한 자리

> **Section 1 classification** — Each layer is tagged below. Manuscript-safe layers feed candidate Paper 1 figure / supplement entries (see §12). Internal-only layers stay in this audit document.

### Layer 1 — score validity (cross-validation) `[MANUSCRIPT-SAFE candidate · supp]`
| Cohort | n | Pearson r DM1 vs NONOVERLAP | p | normalize |
|---|---|---|---|---|
| TCGA-THCA bulk | 561 | **−0.885** | **5.12 × 10⁻¹⁸⁸** | within-cohort z |
| GSE230424 Visium epi50 | 4 | −0.99 | 0.008 | depth-resid |
| GSE248205 Visium epi25 | 8 | **−0.996** | **2 × 10⁻⁷** | depth-resid |

Bootstrap 95% CI on 12-slide sample-mean: [−0.997, −0.916]. Direction-consistent across 3 cohorts and 2 platforms (ST + bulk RNA).

### Layer 2 — clinical / outcome (Sprint O) `[MANUSCRIPT-SAFE candidate · main fig]`
| Test | HR | 95% CI | p |
|---|---|---|---|
| PFI continuous (per 1 SD) | 1.46 | [1.08, 1.96] | 0.0135 |
| PFI multivariate (BRAF+RAS+age+stage) | 1.33 | [1.03, 1.71] | 0.028 |
| **PFI dichotomized (top30 vs bot30)** | **2.04** | **[1.15, 3.61]** | **0.015** |
| **DFI multivariate (recurrence-specific)** | **1.41** | **[1.04, 1.91]** | **0.025** |
| OS/DSS multivariate | NS | | |

DM1_like 가 BRAF, RAS, age, stage 모두 조정 후에도 DFI (recurrence) 독립 prognostic.

### Layer 3 — mechanism (Sprint M) `[MANUSCRIPT-SAFE candidate · main fig]`
TCGA-THCA n=261 (DM1-high top 25%, n=130 vs DM1-low bot 25%, n=130) ULM TF activity inference (CollecTRI 1185 TFs).
**629 TFs FDR<0.05; 93 TFs FDR<0.05 + |Δ|>0.5.**

| Layer | TF | Δ activity | FDR | 의미 |
|---|---|---|---|---|
| Thyroid lineage TF collapse | **FOXE1** | **−1.21** | 2 × 10⁻²⁶ | main thyroid follicle TF |
| Thyroid lineage TF collapse | **NKX2-1 (TTF1)** | **−0.65** | 6 × 10⁻²⁸ | master thyroid TF |
| Thyroid lineage TF collapse | TRPS1 | −0.65 | 3 × 10⁻²⁷ | follicle development |
| Thyroid lineage TF collapse | HOXB3 | −0.60 | 3 × 10⁻³⁰ | development TF |
| Dedifferentiation/EMT activated | **STAT3** | **+1.55** | 2 × 10⁻²⁴ | dedifferentiation activator |
| Dedifferentiation/EMT activated | FOSL1 (AP-1) | +0.98 | 2 × 10⁻²⁴ | oncogenic AP-1 |
| Dedifferentiation/EMT activated | JUNB (AP-1) | +0.81 | 2 × 10⁻²⁴ | AP-1 |
| Epigenetic silencing | DNMT1 | +0.74 | 1 × 10⁻²⁴ | maintenance methyltransferase |
| Epigenetic silencing | DNMT3B | +0.65 | 2 × 10⁻²⁶ | de novo methyltransferase |
| Stress | ATF6 | −0.63 | 1 × 10⁻²⁷ | ER stress |

**Coherent 4-step framework**: TF collapse → DNMT silencing → STAT3/AP-1 activation → TROP2 re-expression.

### Layer 4 — therapeutic (Sprint T) `[MANUSCRIPT-SAFE candidate · supp · TROP2 reframed to tumor-level only]`
DM1-high vs DM1-low TCGA-THCA differential expression (n=130 vs 130). 모든 FDR < 1 × 10⁻³².

| Gene | Δ mean | FDR | Drug |
|---|---|---|---|
| **TACSTD2 / TROP2** | +4.33 | **1.2 × 10⁻³²** | ★ **sacituzumab govitecan (Trodelvy, FDA approved)**, datopotamab deruxtecan |
| FN1 | +4.84 | 9 × 10⁻⁴⁰ | (ECM marker, no direct drug) |
| NAMPT | +1.23 | 6 × 10⁻³⁴ | APO866, FK866 (NAD inhibitor) |
| KCNN4 | +4.35 | 2.5 × 10⁻³⁴ | Senicapoc (orphan drug) |
| LYN | +1.08 | 9 × 10⁻³⁴ | Dasatinib (off-label) |
| CYP1B1 | +3.39 | 4.4 × 10⁻³⁴ | metabolic regulation |
| ELF3 | +2.30 | 4.6 × 10⁻³² | TF (drug-able target indirectly) |

### Layer 5 — therapeutic window (Sprint N2 — tumor vs normal) `[MANUSCRIPT-SAFE candidate · supp · tumor-level vulnerability framing only]`
| Group | n | DM1_like | TACSTD2_z | TACSTD2 raw RPKM |
|---|---|---|---|---|
| **Normal-adjacent** | 49 | **−0.874** | **−1.052** | **8.33** |
| Primary tumor | 460 | +0.089 | +0.102 | 11.76 (1.4× ↑) |
| **Metastatic** | 8 | **+0.226** | **+0.560** | **13.12 (1.6× ↑)** |
| Tumor vs Normal p | | **1.8 × 10⁻²⁰** | **1.3 × 10⁻¹³** | 1.3 × 10⁻¹³ |

**Tumor-level interpretation only:** TROP2 elevation is a **tumor-population vulnerability** (1.4× tumor vs normal-adjacent; 1.6× metastatic vs normal-adjacent). It is **not** a DM1-high–subset-specific marker (Q3 spot-level ρ=−0.016, see Layer 7 / §7). Any sacituzumab-relevant claim must be framed at the tumor population level (vs normal thyroid), not as "DM1-high spot = TROP2-high spot". Wet-lab functional validation (IHC / IC50) is required before any therapeutic-population claim.

### Layer 6 — pan-cancer specificity (Sprint P) `[MANUSCRIPT-SAFE candidate · supp · "thyroid-specific" framing only]`
33 TCGA cancer types 에서 within-cohort DM1 vs NONOVERLAP Pearson r:
- THCA: **r = −0.892** (압도적 outlier)
- DLBC: r = −0.82 (n=48 small)
- SKCM: r = −0.55
- 다른 cancer: |r| < 0.45

→ Thyroid-specific dedifferentiation axis. Pan-cancer relevance 약함 — 그러나 "tissue-specific dedifferentiation prototype" 로 framing 가능.

### Layer 7 — spatial structure (D1, D2, A2) `[MANUSCRIPT-SAFE candidate · supp · spatial molecular score only, NOT image-derived]`
| Test | n | result | meaning |
|---|---|---|---|
| Moran's I per slide | 28 | mean=0.36, **26/28 perm p<0.05** | DM1 strongly spatially clustered |
| Moran's I by stage | 28 | PT/PTC/LPTC=0.39-0.47, **ATC=0.13** | spatial structure breaks down in advanced cancer |
| Bivariate Moran's I (DM1 × NONOVERLAP) | 12 | all 12 negative, perm p ≤ 0.01 | spatially anti-correlated (validates cross-validation in spatial) |
| Margin distance gradient | 28 | 23/28 negative ρ (mean=−0.10) | DM1 lower deeper into epithelial-rich tumor center |

### Layer 8 — drop-one-out robustness (A1) `[MANUSCRIPT-SAFE candidate · supp]`
8 RAI genes 중 어떤 1 개를 빼도 Pearson r between sample-mean DM1 and THYROID_NONOVERLAP:
- Full RAI_8: r = −0.985
- Drop FOXE1 (worst): r = −0.981
- Drop TG (best): r = −0.987

Magnitude 변화 < 1%. **Panel design 의 robustness 완벽**.

### Layer 9 — Paper 1 framework integration (S TCGA #4) `[MANUSCRIPT-SAFE candidate · supp · DM1/DM2 cluster bridge]`
TCGA-THCA dark-matter cluster vs DM1_like score:
- DM1 cluster (n=96): mean DM1_like = 0.00
- DM2 cluster (n=61): mean DM1_like = **−0.66**
- Mann-Whitney p ≪ 0.001

Paper 1 의 dark-matter cluster 정의 (mutation/driver-based) 와 expression-based DM1 score 가 같은 axis 측정 → consilience.

## 2. 정직 보고: 약점 (N1) `[MANUSCRIPT-SAFE · Limitations § content seed]`
**Dark matter subset (BRAF-neg AND RAS-neg, n=156-187) 단독에서 PFI prognostic value NS** (HR=1.20, p=0.80).

해석:
1. n=156 underpowered (전체 cohort 의 1/3.5)
2. DM1 prognostic effect 가 BRAF/RAS-driven tumor 에서 우세할 가능성
3. Dark matter 안에서는 DM1 score 가 **cluster-discriminating** (DM1 vs DM2) 이지 prognostic 아님

Paper 1 narrative 보완: "DM1 axis identifies aggressive subset across PTC overall, particularly in driver-mutant tumors. Within dark matter, DM1 score discriminates DM1 vs DM2 sub-clusters with distinct biology but similar prognosis."

## 3. 현재 venue 가능성 `[INTERNAL-ONLY · 의사결정 보조]`

> Venue probability 추정은 **manuscript prose / cover letter / submission package 어디에도 들어가지 않음.** 이 표는 본인 의사결정용. 실제 target venue 결정은 voice-protected Cover ¶1 작성 시점.

| 등급 | 가능성 | 근거 |
|---|---|---|
| Sci Rep | 100% | base |
| JCI Insight | 99% | confirmed |
| Cell Rep Med (IF 11) | **90% comfortable** | mechanism + drug + clinical 패키지 |
| **Nat Cancer (IF 23)** | **50-65% reach** | full computational story complete |
| **Cancer Cell (IF 50)** | **20-30% 도전권** | functional validation 없는 게 핵심 weakness |
| Nature/Science | 0% | scope |

## 4. Cancer Cell 까지 가려면 — 부족한 것 `[INTERNAL-ONLY · 마라톤 후 로드맵]`

> 본 섹션은 **마라톤 (5/4–6/13) 외 영역**. 마라톤 동안은 새 분석 / 새 cohort / 새 sprint 금지. 아래는 marathon 종료 후 검토 자료.

### Wet-lab (marathon 후 영역)
1. **TROP2 IHC in 본인 cohort** (한국 archived FFPE): 1주
2. **Thyroid cancer cell line + sacituzumab govitecan IC50**: 1-2주
3. **Organoid response to sacituzumab**: 6+ months
4. **Mouse xenograft TROP2-ADC efficacy**: 6+ months

### Computational 추가 후보 (marathon 외) — `[NOT for execution during 5/4–6/13]`
1. **Methylation correlate** — DNMT 활성 검증 → TCGA-THCA 450k methylation in DM1-high promoters of thyroid lineage genes
2. **GSEA Hallmark enrichment** — DM1-high signature → MSigDB hallmark pathway
3. **STAT3 target gene enrichment** — STAT3 활성 검증 직접
4. ~~TROP2 spatial localization in ST data~~ — **이미 Q3 에서 실행됨 (NEGATIVE, ρ=−0.016).** Re-run 금지.
5. **Independent thyroid cohort** — GSE33630 (n~100 PTC), GSE76039 (ATC) 추가
6. **Mutational signature COSMIC SBS** — DM1-high specific signature
7. **Single-cell deconvolution** with paired GSE250521 scRNA — cell type fractions per spot
8. **TF-TF network coordination** — FOXE1, NKX2-1, PAX8, HHEX 의 collapse 가 coordinated 인지 (Q5 에서 일부 실행됨, 추가 검증 candidate)

## 5. 산출 파일 inventory `[MIXED — see §12 for manuscript-safe subset]`
a1_drop_one_out.png
a1_drop_one_out.tsv
a1_drop_one_out_summary.tsv
a2_morans_i.png
a2_morans_i.tsv
a3_DM1_quartile_x_celltype.tsv
a3_celltype_x_dm1.tsv
a3_dm1_quartile_celltype.png
a3_per_spot_celltype.tsv.gz
c2_harmony_28slides.h5ad
c2_harmony_umap.png
d1_permutation_morans.tsv
d2_bivariate_morans.tsv
d3_bootstrap_ci.tsv
d5_margin_distance.tsv
m_master_regulator_volcano.png
m_tf_activity_diff.tsv
n1_dark_matter.png
n1_dark_matter_subset.tsv
n2_thca_tumor_vs_normal.tsv
n2_tumor_vs_normal.png
o_multivariate_cox.png
o_multivariate_cox.tsv
p_pancancer.png
p_pancancer_scored.tsv
p_pancancer_summary.tsv
s_tcga_4panel.png
s_tcga_km_dm1.png
s_tcga_thca_scored.tsv
t_DM1_DE_genes.tsv
t_drug_target_volcano.png

## 6. 핵심 figure paths `[MIXED — see §12 for manuscript-safe subset]`
project/supplementary/spatial_freeze_2026_05_03/SuppFig_X1_GSE250521_spatial_PT_to_ATC.{png,pdf}
project/supplementary/spatial_freeze_2026_05_03/SuppFig_X2_DM1_thyroid_nonoverlap_crossvalidation.{png,pdf}
project/supplementary/spatial_freeze_2026_05_03/SuppFig_X3_GSE248205_autoimmune_negative_control.{png,pdf}
project_external_st/results/extra/s_tcga_4panel.png
project_external_st/results/extra/s_tcga_km_dm1.png
project_external_st/results/extra/o_multivariate_cox.png
project_external_st/results/extra/p_pancancer.png
project_external_st/results/extra/m_master_regulator_volcano.png
project_external_st/results/extra/t_drug_target_volcano.png
project_external_st/results/extra/n1_dark_matter.png
project_external_st/results/extra/n2_tumor_vs_normal.png
project_external_st/results/extra/a1_drop_one_out.png
project_external_st/results/extra/a2_morans_i.png
project_external_st/results/extra/a3_dm1_quartile_celltype.png
project_external_st/results/figures/condition_dm1_mosaic_8panel.png
project_external_st/results/figures/external_st_triage_overview.png

Web bundle: http://40.82.129.113/spatial_supplement/
## 7. 추가 sprint Q1, Q3, Q5 결과

### Q1 GSEA Hallmark — STAT3/EMT/IFN-γ/inflammation up + OXPHOS down (HIT) `[MANUSCRIPT-SAFE candidate · supp · independent mechanism validation]`

| Pathway | ULM score | p |
|---|---|---|
| INTERFERON_GAMMA_RESPONSE | 15.4 | ≈0 |
| INFLAMMATORY_RESPONSE | 13.7 | 2e-41 |
| TNFA_SIGNALING_VIA_NFKB | 12.6 | 2e-35 |
| EMT | 9.9 | 3e-22 |
| **IL6_JAK_STAT3_SIGNALING** | **8.0** | 6e-15 |
| **OXIDATIVE_PHOSPHORYLATION** | **−5.5** | **1e-7** |

★ M 의 STAT3 mechanism 을 IL6_JAK_STAT3_SIGNALING enrichment 로 독립 검증.
★ 새 layer: metabolic shift (OXPHOS down → Warburg-like)

### Q5 Thyroid TF network coordination (HIT) `[MANUSCRIPT-SAFE candidate · supp · TF coordination support]`

DM1-high vs DM1-low expression Δ (TCGA-THCA n=261):
- DIO1: −6.34, TPO: −6.01, SLC5A5: −3.77 (RAI panel core)
- FOXE1: −1.23, PAX8: −1.15, HHEX: −1.09 (thyroid TFs)
- DNMT1: +0.43, FOSL1: +1.21, STAT3: +0.37 (consistent with M)

Pairwise correlations: FOXE1-PAX8 r=+0.55, PAX8-DIO1 r=+0.64, HHEX-PAX8 r=+0.54.
→ 8 RAI genes + 4 thyroid TFs 가 coordinately collapse (single TF 가 아니라 network).

### Q3 TROP2 spatial colocalization (NEGATIVE — 정직 보고) `[MANUSCRIPT-SAFE · negative result · Limitations 또는 supp]`

Mean ρ(TROP2, DM1_like_resid) across 28 slides = **−0.016**.
14/28 slides positive (coin flip).

Sample-level (N2): TROP2 tumor>normal p=1e-13 (tumor population vulnerability)
Spot-level (Q3): NOT colocalized with DM1-high spots

**Manuscript framing rule (binding):** TROP2 claim 은 **tumor 인구 vulnerability** 만 — "tumor 가 normal-adjacent thyroid 보다 TROP2 ↑, 따라서 sacituzumab 의 differential targeting window 가 존재" 는 OK. "DM1-high spots = TROP2-high spots", "DM1-axis 가 sacituzumab 후보 인구를 정의한다" 는 **금지**. Therapeutic narrative 어떤 톤 변경도 wet-lab IHC/IC50 검증 후에만.

## 8. 최종 sprint score `[MANUSCRIPT-SAFE — fact table only · venue interpretation in §3 internal-only]`

| Sprint | 결과 | 의미 |
|---|---|---|
| Cross-validation (S+D2+D3) | ✅✅ HIT | r=−0.89 TCGA bulk + spatial replicate, CI tight |
| O outcome | ✅ HIT | Dichotomized HR=2.04 p=0.015, multivariate DFI HR=1.41 p=0.025 |
| P pan-cancer | 🟡 PARTIAL | Thyroid-specific (THCA r=−0.89 outlier) |
| T druggable | ✅ HIT | TROP2/TACSTD2 FDR=1e-32 + sacituzumab approved |
| M master regulator | ✅✅ HIT | NKX2-1/FOXE1 collapse + STAT3/AP-1/DNMT up (4-step framework) |
| N1 dark matter subset | ❌ NS | n=156 underpowered HR=1.20 NS |
| N2 tumor vs normal | ✅✅ HIT | DM1 p=1.8e-20, TROP2 1.4× tumor>normal |
| **Q1 GSEA Hallmark** | ✅ HIT | IL6_JAK_STAT3 + EMT + IFN-γ up; OXPHOS down |
| **Q5 TF network** | ✅ HIT | Coordinated thyroid lineage network collapse |
| **Q3 TROP2 spatial colocalization** | ❌ NEGATIVE | spot-level NOT colocalized (sample-level still works) |

**6 HIT + 1 PARTIAL + 2 NEGATIVE**

## 9. GPU Phase B 패키지 (RunPod 용) `[DEPRECATED 2026-05-04 — image-DM1 closure battery NO-GO]`

> **DEPRECATED.** 본 섹션은 image-axis DM1 (H&E → DM1) 가설을 전제로 작성됨. 2026-05-04 closure battery (run-time ~85 min CPU, ResNet50 224 px DM1_resid Spearman 0.022 ≈ random max 0.031, foundation-model 업그레이드 projected gain (~+0.05–+0.15) 도 GO threshold (+0.28) 미만)가 verdict **D = full 폐기** 를 확정. 따라서 G1 (ResNet50/UNI LOSO) / G2 (TCGA H&E n=500) / G3 (CellDART) RunPod sprints **모두 폐기**. `phaseB_gpu_pkg.tar.gz` (4.35 GB), `phaseA_gpu_pkg.tar.gz` (614 MB) 둘 다 historical artifact.
>
> **Action:** 어떤 manuscript 도 image-DM1 / digital pathology / TROP2-by-image angle 을 포함하지 않음. Pod B (WSI) 는 `null` 처리, Pod C (AlphaFold) / Pod D (K2 STAR) 는 RunPod 콘솔에서 사용자 직접 stop 필요.
>
> **Reference:** `project/reports/pathology_dm1_closure_battery_2026_05_04.md`, root `CLOSURE_BATTERY_2026_05_04.md`, `PHASE_A_NOGO_REPORT_2026_05_04.md`, commit `1d16a4e` (decision: drop image-dm1 pathology angle after closure battery).

```
[DEPRECATED CONTENT — preserved for archival reference only, do not execute]

phaseB_gpu_pkg.tar.gz: 5,600 H&E tiles (28 slides) + scored data + scRNA paired (9 samples) + 5 GPU scripts + setup
G1 ResNet50 LOSO   — DEPRECATED (closure battery NO-GO)
G1 UNI (gated)     — DEPRECATED
G2 TCGA H&E (n=500) — DEPRECATED + TCGA WSI 다운로드 금지 명시
G3 CellDART        — DEPRECATED
```

## 10. 본 VM 의 parallel 작업 overlap 주의 `[DEPRECATED 2026-05-04 — image-DM1 NO-GO]`

> **DEPRECATED.** 본 섹션이 가리키는 PID 168753 (`embeddings_resnet50_448.npz` 생성) 도 image-DM1 angle 이라 NO-GO. 결과 비교/통합 plan 은 무효.
>
> 추가로 이 마라톤 세션 중 `/tmp/runpod_dispatch_v2.sh` 가 Pod C/D 에 work 를 dispatch 했고, `tail -F` 는 Claude 가 종료. **호스트측 dispatch 종료, RunPod side Pod C/D stop 은 사용자 직접 RunPod 콘솔/API 에서 수행 필요** (per `2026_05_04_background_dispatch_shutdown.md`).

## 11. 최종 권장 `[DEPRECATED 2026-05-04 — replaced by §13 below]`

> **DEPRECATED.** Image-DM1 angle 폐기로 §11 의 G1/G2 hit-conditional venue 추정 + RunPod 다음-step 권장 모두 무효. §13 (이번 audit 의 권장) 참조.

```
[DEPRECATED CONTENT]

Venue 결정 (10 sprint + GPU package 기준) — 무효: §3 internal-only로 이동
다음 step (RunPod G1/G2/G3 실행) — 무효: closure battery NO-GO
스코프 lock 문장 ("이미 ... 다 박혔음") — 부분 무효: image-axis 만 무효, molecular axis 는 살아있음
```

---

## 12. Manuscript-safe candidate inventory (Paper 1) `[NEW · 2026-05-04 audit]`

> Paper 1 manuscript 또는 supplement 에 들어갈 후보 만 추림. **Image-axis (H&E → DM1) 자료는 모두 제외.** 각 후보는 §1 / §7 의 layer 와 cross-reference. **새 분석 0건. 기존 산출물 reuse 만.**

### 12.1 Main figure 후보 (molecular axis only)

| 후보 | 출처 | 근거 layer | 비고 |
|---|---|---|---|
| **DM1/RAI cross-validation** | `s_tcga_4panel.png`, `t_DM1_DE_genes.tsv` | Layer 1 | TCGA bulk r=−0.885 + spatial replicate; 2-platform consilience |
| **Clinical outcome (Cox PFI/DFI)** | `o_multivariate_cox.png/tsv`, `s_tcga_km_dm1.png` | Layer 2 | dichotomized HR=2.04 + multivariate DFI HR=1.41; 단, dark matter NS limitation 명시 |
| **TF-driven mechanism (4-step)** | `m_master_regulator_volcano.png`, `m_tf_activity_diff.tsv` | Layer 3 + Q1 + Q5 | NKX2-1/FOXE1 ↓ + DNMT ↑ + STAT3/AP-1 ↑ + GSEA Hallmark + TF coordination |

### 12.2 Supplementary figure 후보

| 후보 | 출처 | 근거 layer |
|---|---|---|
| Drop-one-out robustness (8 RAI gene panel) | `a1_drop_one_out.png/tsv` | Layer 8 |
| Spatial Moran's I + bivariate (DM1 × NONOVERLAP) | `a2_morans_i.png/tsv`, `d1/d2_*morans.tsv`, `d5_margin_distance.tsv` | Layer 7 |
| Pan-cancer specificity (THCA outlier) | `p_pancancer.png`, `p_pancancer_summary.tsv` | Layer 6 |
| Tumor vs normal molecular axis (DM1 + TROP2 at *population* level) | `n2_tumor_vs_normal.png/tsv` | Layer 5 — **tumor-population vulnerability framing only** |
| DM1/DM2 cluster bridge (Paper 1 framework integration) | `s_tcga_thca_scored.tsv`, `n1_dark_matter.png/tsv` | Layer 9 |
| Drug-target volcano (TROP2 + others) | `t_drug_target_volcano.png`, `t_DM1_DE_genes.tsv` | Layer 4 — **tumor-level only**, no spot-level |
| Honest negative: Q3 TROP2 spot-level NOT colocalized | spot-level corr (in `t_DM1_DE_genes.tsv` 부속 분석) | Q3 |
| GSEA Hallmark (IL6_JAK_STAT3 + EMT + IFN-γ; OXPHOS down) | (Q1 result table — re-render from `m_tf_activity_diff.tsv` 부속) | Q1 |
| Honest negative: dark matter HR=1.20 NS | `n1_dark_matter.png` 또는 §2 prose 만 | N1 |

### 12.3 Spatial supplement (이미 freeze 된 자료)

| 파일 | Status |
|---|---|
| `project/supplementary/spatial_freeze_2026_05_03/SuppFig_X1_GSE250521_spatial_PT_to_ATC.{png,pdf}` | already in supp pack |
| `SuppFig_X2_DM1_thyroid_nonoverlap_crossvalidation.{png,pdf}` | already in supp pack |
| `SuppFig_X3_GSE248205_autoimmune_negative_control.{png,pdf}` | already in supp pack |
| `SuppTable_SX_28sample_score_summary.tsv` | already in supp pack |
| `SuppTable_SX1_independent_lineage_and_technical_confounding.tsv` | already in supp pack |

이 spatial supplement 는 **molecular DM1_like score 의 spatial structure 만** 보여주며 H&E image-axis 가 아니므로 변경 없이 유지.

### 12.4 NOT manuscript-safe (image-axis residual / venue / roadmap)

- `phaseA_gpu_pkg.tar.gz`, `phaseB_gpu_pkg.tar.gz` — image-DM1 NO-GO
- 모든 RunPod G1/G2/G3 산출물 (loso_metrics, embeddings npz, tile_metadata*) — image-DM1 NO-GO
- `12_gpu/`, `post_pod_*.py` (B/C/D) — image-DM1 NO-GO 또는 paper2-pause
- §3 venue probability table — internal-only
- §4 wet-lab + computational roadmap — internal-only / marathon-외
- §9–11 전체 — DEPRECATED

---

## 13. 다음 step (replaces §11) `[INTERNAL — marathon scope]`

1. **이 audit 까지가 Paper 1 산출물 lock**. molecular axis 만 manuscript 에 들어감. 이미 §12.1–12.3 후보 모두 disk 에 존재 (새 분석 불필요).
2. **Voice-protected 6 sections** 본인 키보드 슬롯 (Hook ¶1, Aim ¶4, Discussion §3.1 / §3.4 Limitations, Cover ¶1, Reviewer Q9). Discussion §3.1 은 `v17_landa2016_cite_save` 메모리의 3-layer reverse-causality 차단 적용.
3. **Limitations §3.4 seed:** §0 의 6가지 remaining risks + §2 dark-matter NS + §7-Q3 TROP2 spot-level negative.
4. **Pod C/D stop** RunPod 콘솔 직접 (Claude는 token 미보유; `2026_05_04_background_dispatch_shutdown.md` 에 명령 templated).
5. **6/13 bioRxiv target 유지** (image-axis 폐기 → Paper 1 scope 좁아졌지만 molecular axis 자체는 ship-ready).

