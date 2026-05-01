---
title: "Paper 1 manuscript v8 — single-page outline (Mermaid + bullets) — v3"
date: 2026-04-30
revised: 2026-04-30 (v3 — web Claude 2차 review 후 Landa cite save + 5 추가 cites + dial-back framing + Fig 1 4 panels)
author: Seungho Cook
target_venue: Cell Reports Medicine (1순위) → JCI Insight + Nat Commun dual reach → npj Precision Oncology (fallback)
total_target_words: ~5,500-7,000 (excluding STAR Methods + figure captions)
status: v3 — Landa 2016 cite save + Yoo/TCGA/Pu/Bradley/Wirth 5 cites + Fig 1 4 panels + Hook Alt B refined + ATA dual cite
---

# Paper 1 (manuscript v8) — single-page outline (v3)

## 0a. v2 → v3 변경 7건 (web Claude 2차 review)

| # | 변경 | 근거 |
|---|---|---|
| 1 | Hook: Alt B refined ("intermediate-risk ~20%" + "~23% of cases" specific) | ATA 3.2 narrative arc thread 강화; 본인 voice "compass" 옵션 hybrid 보존 |
| 2 | Discussion 3.1 cite ★ Landa 2016 JCI 정정 (Krishnamoorthy 2025 Nat Comm misattribution catch — 4 fact 모두 틀림) | Landa I et al 2016 JCI 126(3):1052-1066 = Krishnamoorthy GP 9번째 공저자 + 8-gene 5/8 직접 1:1 overlap (TG/TSHR/TPO/PAX8/DIO1) |
| 3 | Discussion 3.1 dial-back framing ("identifies an upstream signature consistent with" 대신 "extends paired-cancer continuum upstream") + speculative tone ("may identify") | reviewer pushback 차단; Landa wording 에 가까운 톤 |
| 4 | Reverse-causality 3-layer 차단 strategy: Methods + Discussion 3.1 convergence + Cover Letter | 8-gene 5/8 = Landa ATC silenced gene list overlap → reviewer "panel-driven 의심" 사전 차단; **independent design 의 convergent validation 으로 strength 변환** |
| 5 | 5 추가 cite: Yoo 2016 PLOS Genet (panel origin, Methods + Intro 1.2 + Disc 3.1), TCGA 2014 Cell (Intro 1.2 canonical), Pu 2021 Nat Commun (Methods + Disc 3.3), Bradley 2010 (Disc 3.1 BRAF V600E HLA-I counter-intuitive), Wirth 2020 NEJM (Disc 3.2 LIBRETTO-001 specific) | web Claude Q7 catch — Yoo 2016 + TCGA 2014 빠지면 reviewer 1번에 잡힘; Bradley 2010 은 본 paper R3 finding 14 contradict cite; Wirth 2020 은 Disc 3.2 actionability specific |
| 6 | ATA 2025 (Ringel 2025 *Thyroid* 35(8):841-985, PMID 40844370) 보조 cite Discussion 3.2 | 2026 submission 시 ATA 2025 미반영 = reviewer Q risk; 1순위 ATA 2015 + 보조 ATA 2025 dual |
| 7 | **Fig 1: 5 panels → 4 panels** (Heatmap C → S1 supp 이동) | Cell Press first impression density 낮춤; Heatmap 은 standard, supp 가능 |

## 0b. v1 → v2 변경 5건 (web Claude 1차 review 적용 — 보존)

| # | 변경 | 근거 |
|---|---|---|
| 1 | Title: Hybrid Cand 1 (108 chars, "An 8-gene panel reveals fusion-driven, epigenetically silenced dark matter in BRAF/RAS-negative thyroid cancer") | venue ladder 전체 cover; Cand 1 약점 (no panel) + Cand 2 약점 (no mechanism) 동시 보완 |
| 2 | Results ordering 대안 A 채택: 2.1 → **2.3 → 2.4** → 2.2 → 2.5 | mechanism-first narrative arc; biomarker-first 인식 risk 회피; Cell Press editorial 선호 (mechanism story first, clinical second) |
| 3 | Section 2.4 두 sub-paragraph split: 2.4a (epigenetic R5-2) + 2.4b (immune-overlap teaser, R4-2 fusion- DM1) | Pillar 5 main story (HLA n=874 mediation, BCR/TLS, sub-B NBNR) 는 **Paper 2 reserve**. Paper 1 에서는 fusion- DM1 = older + immune-hot 까지만 teaser. Paper 2 differentiation 보호. |
| 4 | Fig 8 (D) HMA schematic → **supplementary 이동**. Main Fig 8 = 3 panels (β heatmap + mean β bar + fusion-independent). | Schematic 은 mechanism 그림이지 data 아님 → Cell Rep Med editor "data figure 에 schematic = weakness" 인식 risk |
| 5 | Discussion 3.1 + Krishnamoorthy 2025 PDTC/ATC Nat Comm cite (paired cancer continuum framing). 3.2 + ATA 2015 guideline alignment + LIBRETTO-001 specific. Limitations 6 paragraph → 4 main + 2 supp. Venue ladder language: "Nat Comm reach" → "Nat Comm stretch" + JCI Insight dual reach 명시. | web Claude Q5/Q7/Q8: ATA/Krishnamoorthy missing critical cite, Limitations 4 main + 2 supp 가 Cell Rep Med 표준, Nat Comm 단독 reach 는 30-45% (분당 X) — JCI Insight (75-85%) dual reach 가 정직 |

---

## 1. Narrative arc (Mermaid v2)

```mermaid
flowchart TD
  T[Title v2: Hybrid 'An 8-gene panel reveals<br/>fusion-driven, epigenetically silenced<br/>dark matter in BRAF/RAS-negative thyroid cancer']
  T --> AB[Abstract v2 — 153 words structured<br/>DM1 76.8% fusion+ · TPO d=2.30 · HR 2.53<br/>· mean β 0.385 vs 0.253 · pan-genome ARI 0.92]

  AB --> S1[Section 1 · Introduction · 600-900 words]

  S1 --> S1A[1.1 PTC clinical heterogeneity<br/>RAI decision unmet need · Bethesda III/IV]
  S1A --> S1B[1.2 Existing molecular framework<br/>TCGA 2014 BRAF-like/RAS-like · Yoo 2016 BRS]
  S1B --> S1C[1.3 Dark matter concept<br/>Xing 2014 · BRAF-/TERT- 23-37%]
  S1C --> S1D[1.4 Aim + preview<br/>8-gene panel · DM1/DM2 · 3-layer DM1]

  S1D --> S2[Section 2 · Results · 2,500-4,000 words<br/>★ 대안 A ordering ★]

  S2 --> S2A[2.1 Panel + DM1/DM2 cluster discovery<br/>Pillars 3+4 · Fig 1]
  S2A --> S2C[★ 2.3 DM1 fusion paradigm<br/>76.8% RET/NTRK/ALK/BRAF · OR 7.41 · Fig 7 A-C]
  S2C --> S2D1[★ 2.4a Epigenetic silencing<br/>TPO d=2.30 · fusion-independent · Figs 7D, 8]
  S2D1 --> S2D2[2.4b Fusion- DM1 = immune-overlap teaser<br/>Pillar 5 hint · Paper 2 reserve]
  S2D2 --> S2B[2.2 Clinical aggressiveness + Meta<br/>HR 2.53 [1.31, 4.89] · Figs 2, 6]
  S2B --> S2E[2.5 Cross-cohort + Reflex algorithm<br/>4,300+ EA · Figs 3, 5]

  S2E --> S3[Section 3 · Discussion · 1,000-1,500 words]

  S3 --> S3A[3.1 Three-layer DM1 + Krishnamoorthy 2025<br/>paired cancer continuum framing]
  S3A --> S3B[3.2 ATA 2015 alignment + LIBRETTO-001<br/>+ HMA + RAI re-induction rationale]
  S3B --> S3C[3.3 East-Asian generalizability]
  S3C --> S3D[3.4 Limitations · 4 main paragraphs]

  S3D --> SM[STAR Methods · unlimited<br/>cohorts · pipelines · stats · code release]
  S3D --> FG[Figures 1-8 · Cell Press main 8 + Suppl 8]
  S3D --> ST[Supplementary Tables S1-S10 + 2 limitation paragraphs]
```

---

## 2. Section-by-section bullet outline (v2)

### Section 1 · Introduction (600-900 words) — v3 cite 강화

| # | Sub | Length | Content | Citations (v3) |
|---|---|---|---|---|
| 1.1 | Clinical context | 150-200 w | PTC global incidence ↑; RAI ablation over-treatment concern; Bethesda III/IV indeterminate FNA unmet need; BRAF V600E status alone insufficient | Haugen 2016 *Thyroid* (ATA 2015); Ringel 2025 *Thyroid* (ATA 2025) |
| 1.2 | Existing molecular framework | 150-200 w | **TCGA 2014 BRAF-like/RAS-like spectrum (canonical origin cite)**; **Yoo 2016 BRS 273-gene + TDS-core (8-gene origin paper)**; 23% BRAF-/RAS- gap | ★ TCGA 2014 *Cell* 159(3):676-690 (필수); ★ Yoo 2016 *PLOS Genet* 12(8):e1006239 (필수) |
| 1.3 | Dark matter concept + paired-cancer continuum | 200-250 w | Xing 2014 NEJM; BRAF-/TERT- 23-37%; East Asian 37.8% (Korean) vs 28.4% (TCGA); existing panels insufficient; **Landa 2016 JCI advanced thyroid (PDTC+ATC) dedifferentiation trajectory + thyroid differentiation transcript suppression** | Xing 2014 *NEJM*; Liu 2017; ★ **Landa I et al 2016 *JCI* 126(3):1052-1066** (Krishnamoorthy 2025 Nat Comm 정정) |
| 1.4 | Aim + brief preview | 100-150 w | 8-gene RAI panel (independently selected from canonical RAI biology, Yoo 2016); DM1/DM2 axis orthogonal to BRAF/RAS framework; 3-layer DM1 pathology; multi-cohort validation + clinical actionability | self |

**★ Hook 1줄 (v3 본인 voice 적용 영역) — Alt B refined 권장:**

> "Despite a >98% 5-year overall survival in differentiated thyroid carcinoma, structural disease recurrence reaches ~20% in ATA 2015 intermediate-risk patients (Haugen et al., 2016), and BRAF/RAS-negative tumors — accounting for ~23% of cases — complicate radioiodine decisions in the absence of mechanistic sub-stratification."

**Hybrid 옵션 ("compass" 본인 voice 보존):**

> "Despite a >98% 5-year overall survival, structural disease recurrence reaches ~20% in ATA 2015 intermediate-risk papillary thyroid carcinoma (Haugen et al., 2016), and clinicians continue to make radioiodine decisions on heterogeneous BRAF/RAS-negative tumors — ~23% of cases — without a mechanistic compass."

본인 결정 영역. v2 의 "5-15%" → "intermediate-risk ~20%" specific number; "~23% of cases" 추가 → BRAF/RAS-negative dark matter scope 정량화 (Xing 2014 cite gateway).

### Section 2 · Results — ★ 대안 A ordering 채택 (mechanism-first → clinical-last)

5 sub-results × 500-800 w (총 2,500-4,000w)

| # | Sub-result | Pillars | Key claims | Figures | Source |
|---|---|---|---|---|---|
| **2.1** | 8-gene panel + DM1/DM2 cluster discovery | P3 + P4 | TIERA67 7-category curation; drivers excluded by design (NOT exclusion — implicit category restriction); 8-gene = TDS_core subset (Yoo 2016); ΔAUC 8 vs 16 = 0.013 NS; **pan-genome top-5000 ARI 0.92 ≈ TIERA67 0.90; 8-gene alone 0.49; Driver_anchor ARI=−0.007**; BRAF mRNA d=−0.04 mutation-neutral | Fig 1 (5 panels) | v4 §5.1; GRAND P1, P4 |
| **2.3** ★ | DM1 fusion paradigm | (R3-F4 paradigm) | DM1 76.8% fusion+ (63/82): RET 33, NTRK 10, ALK 4, BRAF 5; DM1 vs DM2 OR 7.41 (sensitivity 7.18-9.07); MSK SV cross-validation (RET 5, ALK 3, PAX8-PPARG 3); missingness MAR (chi² p=0.56) | Fig 7 A-C | v4 §5.0 R3-F4, R4-1, R5-1 |
| **2.4a** ★ | DM1 epigenetic silencing of differentiation machinery | (R5-2 NEW) | Promoter hypermethylation TPO d=2.30, DIO1 d=1.24, TSHR d=1.20, PAX8 d=0.97, TG d=0.86, FOXE1 d=0.84, NKX2-1 d=0.63 (SLC5A5 NS d=0.22); mean 8-gene β 0.385 (DM1) vs 0.253 (DM2); methylation **fusion-independent** (within-DM1 fusion+ vs - d=−0.36 NS); 3-layer model L1+L2+L3 | Figs 7D, 8 (3 panels) | v4 §5.0 R5-2 |
| **2.4b** | Fusion-negative DM1 = immune-overlap subtype (teaser) | P5 partial | Sub-A vs sub-B: age (37 vs 51), fusion% (84.7 vs 57.9), stage III/IV (15.3% vs 44.4%), CD8/IFN-γ low vs high (d=-0.5 to -0.6), Hashimoto-like 40% vs 67% (OR 0.34, p=0.064 trend); **Pillar 5 mediation/HLA n=874/BCR clonal/sub-B NBNR is Paper 2 backbone — Paper 1 stops at "older + immune-hot DM1 sub-B is mechanistically distinct"** | Fig 7D | v4 §5.0 R4-2 (Paper 1 scope only) |
| **2.2** | Clinical aggressiveness + Meta-analysis | (Tier 1) | Xing 73% rescue (180 BRAF-/TERT- → 131 sub-stratified); TCGA HR 2.30 [0.77, 6.88]; MSK HR 2.67 [1.17, 6.10]; **Pooled HR 2.53 [1.31, 4.89] I²=0%** (DerSimonian-Laird) | Figs 2, 6 | v4 §5.4; R2-N1 |
| **2.5** | Cross-cohort + Reflex testing algorithm | P1 + (translational) | sc external (Pu 2021 r=0.798-0.886; Lu 2023 thyrocyte-intrinsic); FFPE robust (KS p=0.44); Korean K2/Lee/GSE286332 n=874 + GSE213647 n=632 (28% Hashimoto-like); **DM1 captures 81.8% TCGA RET+ → reflex algorithm**; population estimate 48 selpercatinib-eligible per 1000 PTC | Figs 3, 5 | v4 §5.5, R4-3; D8-B |

**각 sub-result 끝 take-home one-sentence (Cell Press style, narrative arc):**
- 2.1: "...defining a transcriptional axis orthogonal to canonical BRAF/RAS classification."
- 2.3: "...elevating DM1 from biomarker to fusion-driven mechanism."
- 2.4a: "...adding an epigenetic mechanism layer independent of fusion driver presence."
- 2.4b: "...with fusion-negative DM1 representing a mechanistically distinct immune-overlap subtype warranting separate inquiry." → Paper 2 hook
- 2.2: "...quantifying clinical aggressiveness within Xing 2014 dark matter."
- 2.5: "...enabling immediate clinical use as a reflex fusion-testing trigger."

### Section 3 · Discussion (1,000-1,500 words) — v3 cite 강화

| # | Sub | Length | Content | Tone | New citations (v3) |
|---|---|---|---|---|---|
| 3.1 | Three-layer DM1 pathology + ★ Landa 2016 dial-back framing + Bradley 2010 contradict + Yoo 2016 convergence | 350-450 w | L1 genetic (76.8% fusion+) · L2 mechanism heterogeneity (sub-A vs sub-B) · L3 epigenetic (TPO d=2.30 fusion-independent); ★ **Landa et al 2016 JCI 126(3):1052-1066** characterized advanced thyroid cancer (84 PDTC + 33 ATC) genomic landscape, establishing PDTCs/ATCs arise from well-differentiated tumors with TERT promoter stepwise (9% PTC → 40% PDTC → 73% ATC) and **profound suppression of differentiation transcripts (TG, TSHR, TPO, PAX8, SLC26A4, DIO1, DUOX2)**; DM1 framework "**identifies an upstream signature consistent with this dedifferentiation trajectory**" — same machinery silenced at ATC end already epigenetically attenuated in DM1 PTCs. **Convergence framing**: 8-gene panel (Yoo 2016 PLOS Genet RAI biology) + Landa 2016 ATC silenced gene list — independent designs sharing 5/8 genes = convergent biological validation, NOT panel-driven circular inference. **Bradley 2010 contradict**: BRAF V600E HIGHER HLA-I (counter-intuitive, R3 finding 14) — alternate immune evasion mechanism in BRAF-driven PTC | mechanism-forward, dial-back | ★ Landa 2016 *JCI*; Yoo 2016 *PLOS Genet*; Bradley 2010 (citation verify); Krishnamoorthy 2025 *JEM* (broad context closing line only) |
| 3.2 | Clinical actionability — ★ ATA 2015+2025 dual + Wirth 2020 LIBRETTO specific + HMA rationale | 300-400 w | **ATA 2015 (Haugen 2016) — recently updated as ATA 2025 (Ringel 2025) — incorporates BRAF V600E as sole molecular risk modifier** (fusion drivers and epigenetic silencing not reflected in current risk tiers); BRAF/RAS-negative dark matter (~23%) currently defaults to clinico-pathological-only stratification; **DM1 reflex algorithm = orthogonal molecular axis** providing ATA intermediate-risk tier additional input. 57% DM1 actionable (RET 40%, NTRK 12%, ALK 5%, BRAF 6%). DM1 captures 81.8% TCGA RET+; 48 selpercatinib-eligible/1000 PTC. **Selpercatinib FDA accelerated approval May 2020 based on LIBRETTO-001 (Wirth et al 2020 NEJM, ORR 79% in RET-fusion+)** — DM1 RNA-score positivity provides upstream candidate pool for prospective evaluation. HMA + RAI re-induction rationale: TPO/DIO1/TSHR promoter hypermethylation → epigenetic-targeted RAI re-induction (decitabine retrospective ATC trials NCT00085293, NCT01065090); SLC5A5/NIS exception → combination strategies (HMA + lithium for NIS membrane trafficking). **Evaluation, not trial.** | translational | Haugen 2016 ★; ★ Ringel 2025 *Thyroid*; ★ Wirth 2020 *NEJM* |
| 3.3 | East-Asian generalizability + Pu 2021 cross-cohort | 200-300 w | 4,300+ East Asian PTC tumors covered; Korean DM 37.8% vs TCGA 28.4%; K2 NBNR mixed phenotype (vascular invasion); Pan-Asian Hashimoto-like replication (Korean GSE213647 22-28% ≈ TCGA 18-20%); **single-cell external validation in 6 PTC patients (Pu et al 2021 Nat Commun, GSE184362) confirmed thyrocyte-intrinsic DM1 signal** + Lu 2023 thyrocyte-intrinsic | population genomics | Wang 2024 *Endocr Connect*; Liu 2017; ★ Pu 2021 *Nat Commun* 12:6058 |
| 3.4 | Limitations — **4 main paragraphs** | 200-300 w | (i) Statistical: TCGA event rate 3.2% — mitigated by N1 meta. (ii) Cohort: MSK advanced enriched, Korean retrospective, Bundang prospective 0%. (iii) Methodological: HLA d 1.5+ partial autocorr (Δd=0.30); K2 mini-index calibration mismatch (R4-4) propagates GSE286332 reference. (iv) Translation gap: clinical cutoff TBD, FDA companion diagnostic TBD; **DM1 framework prospective integration with metastasis-specific molecular axes (Krishnamoorthy 2025 RBM10) merits future work**. **2 추가 limitation paragraphs는 supplementary 로 이동.** | honest disclosure | self + Krishnamoorthy 2025 *JEM* (1 line forward) |

**Supplementary 로 이동 (2 paragraph 분량):**
- SV missingness 2.7% MAR (이미 reviewer Q5 답변됨)
- GSE286332 calibration mismatch (Paper 2 critical issue, Paper 1 에선 Paper 2 deferred 만 명시)
- Lu 2023 r=0.97 autocorrelation (이미 N7 disclosure)
- Wang 2024 Scenario C
- R5-2 single-probe-per-gene methylation aggregate

**본인 voice 보호 영역:** 3.1 첫 paragraph (paper 의 진짜 mechanism story) + 3.4 Limitations (본인 정직 disclosure 정신).

### STAR Methods (unlimited)

5 sub-sections (drafted in Prompt 6):
1. Key Resources Table
2. Resource availability
3. Method details (8 sub-modules: cohort assembly, panel selection, DM cluster, sc validation, survival/meta, cBioPortal API, Korean processing, software)
4. Quantification + statistical analysis (Cohen's d, MW, Fisher, BH-FDR, bootstrap, mediation, random-effects meta)
5. Additional resources (GENCODE v44, IMGT/HLA v3.55, TIERA67 67-gene definition)

### Figures (Cell Press main 8 + Suppl 8) — v3: Fig 1 5→4 panels + Fig 8 (D) supp 이동

| Fig | Title | Panels (v3) | Source code |
|---|---|---|---|
| **Fig 1 (v3 = 4 panels)** | 8-gene panel definition + DM1/DM2 discovery | A-Sankey TIERA67 / B-UMAP / **C-Driver mRNA neutrality (BRAF V600E vs WT box, d=−0.04)** / **D-Pan-genome ARI bar (8-gene 0.49 vs TIERA67 0.90 vs top-5000 0.92 vs Driver-only −0.007)**. ★ **Heatmap → S1 supp 이동** (Cell Press first impression density 낮춤) | v17_8gene_figs_v2.py + p4_pangenome_vs_tiera67.py |
| **Fig 2** | Clinical aggressiveness within Xing rescue | A-Xing Sankey 73% / B-Korean DM 37.8% vs TCGA 28.4% / C-DM1 vs DM2 KM | v17_4way_figure.py + new |
| **Fig 3** | Single-cell validation | A-Pu 2021 6 patients heatmap / B-Lu 2023 thyrocyte-intrinsic UMAP / C-author independence GSE241184 vs GSE184362 | v17_GSE241184_scrna.py + existing sc figs |
| **Fig 4** | Mutation × TERT × outcome | A-8-cell decomposition stacked / B-BRAF_TERT+ vs OTHER_TERT- KM / C-4-patient OTHER_TERT+ caveat box | (existing) |
| **Fig 5** | Korean validation + FFPE | A-GSE213647 HLA d=0.95 / B-K2 NBNR mixed phenotype / C-FFPE vs FF KS p=0.44 | v17_KOREAN_K2_v260_figure.py + existing |
| **Fig 6** | META-analysis forest | A-TCGA HR 2.30 + MSK HR 2.67 + Pooled 2.53 / B-I²=0% / C-DM1 sub-A vs sub-B Cox | v17_D4P1_forest_meta.py |
| **Fig 7 ★** | DM1 mechanism (paradigm) | A-Sub-A vs Sub-B silhouette 0.584 / B-Fusion enrichment 76.8% vs 30.9% / C-Fusion partner stack (RET CCDC6/NCOA4 + NTRK + ALK + BRAF) / D-Fusion+ vs fusion- phenotype / E-Reflex algorithm 81.8% capture | v17_audit_F2_F3_F4.py + R4 + new |
| **Fig 8 ★ (v2 = 3 panels)** | Epigenetic layer (NEW R5-2) | A-Per-gene β heatmap DM1 vs DM2 (8-gene + DIO2/SLC26A4) / B-Mean 8-gene β bar (DM1 0.385 vs DM2 0.253 vs not_DM 0.356) / **C-Within-DM1 fusion+ vs - methylation NS (d=−0.36, p=0.31)** | v17_audit_R5_all.py + new |

**Suppl figures (S1-S9, v3 — Fig 1 Heatmap 추가 이동):**
- ★ **S1 8-gene heatmap sorted by P_DM1** (Fig 1 C 에서 이동) — TCGA + cohort comparison
- S2 sc gene panel similarity matrix
- S3 pan-genome ARI ladder (panel sizes 8/16/67/200/1000/5000)
- S4 HLA-II residualization (DM1 vs DM2)
- S5 arcasHLA Korean Pan-Asian forest (3-arm with **Chu X et al. 2018** — 본인 catch 5/3 정정, Chen 2018 NOT) — Pillar 1, but Paper 1 supp; main role in Paper 2
- S6 BCR clonal + TLS heatmap — Pillar 5 partial supp; Paper 2 main
- S7 DM1 sub-B × K2 NBNR signature transfer — Paper 2 main
- S8 ★ **HMA + RAI re-induction schematic + literature meta** (Fig 8 D 이동) — schematic + decitabine-RAI ongoing trials review (NCT00085293, NCT01065090)
- S9 K2 mini-index calibration FAIL diagnostic + alternate evidence chain

### Supplementary Tables (S1-S10)

| # | Title | Source |
|---|---|---|
| S1 | TIERA67 67-gene definition (7 categories) | metadata/tierA67_genes.txt |
| S2 | Cohort overview (n / RNA / WGS / methylation / SV / clinical) | v4 §9 |
| S3 | TCGA-THCA 8-gene per-sample scores + DM call + clinical | results/tables/ |
| S4 | TIERA67 univariate Cohen d ranking (full 67) | p1_driver_mrna_audit/top20_by_d_full_tiera67.tsv |
| S5 | Pan-genome ARI ladder (panel sizes 8/16/67/200/1000/5000) | p4_pangenome_vs_tiera67/ari_comparison.tsv |
| S6 | DM1 fusion details (per-sample, fusion class, partner) | round3/cbio_sv_thca.tsv + R4-3 reflex |
| S7 | DM1 vs DM2 per-gene methylation β (8-gene + DIO2/SLC26A4) | round5/r5_2_per_gene_methylation_DM.tsv |
| S8 | Meta-analysis raw inputs (TCGA + MSK Cox + DerSimonian-Laird) | results/round2/n1_meta.json |
| S9 | Korean Pan-Asian HLA per-sample arcasHLA genotypes | d4p1_panasian_meta/korean_PTC_pool_n908.tsv |
| S10 | Reviewer Q&A pre-empt 12 items | manuscript_v8/09_reviewer_qa.md |

---

## 3. Word-budget cumulative (v2)

| Section | Target | Cumulative |
|---|---|---|
| Title + Abstract | 153 | 153 |
| Introduction | 750 | 903 |
| Results | 3,250 | 4,153 |
| Discussion | 1,250 | 5,403 |
| Limitations (main 4 paragraph) | 250 | 5,653 |
| STAR Methods | 2,000-3,000 | 7,653-8,653 |

**Cell Reports Medicine main text constraint:** ~5,500-7,000 words (excluding STAR Methods). Current draft fits.

---

## 4. Venue ladder (v2 정직화)

| Venue | v1 평가 | v2 평가 (web Claude) | 통과 확률 (분당 X / 분당 O) |
|---|---|---|---|
| Cell Reports Medicine | 1순위 | **1순위 유지** | 65-75% (분당 X), 70-80% (분당 O) |
| **JCI Insight** | (언급 없음) | **★ Reach dual** | **75-85%** — 진짜 합리적 alternative |
| Nature Communications | reach | **stretch dual** | 30-45% (분당 X), 50-65% (분당 O) — Bundang 결정 |
| **Genome Medicine** | (언급 없음) | 예비 | 60-70% |
| npj Precision Oncology | fallback | **fallback 유지** | 95%+ 안전권 |
| ~~Nature Medicine~~ | ~~lottery~~ | **제외** | 5-10% — 본인 reality check 정렬 |

**v2 strategy**: Cell Rep Med 1순위 + **JCI Insight & Nat Comm dual stretch** (parallel pre-submission inquiry) + npj Prec Onco fallback.

**Bundang 0% 가 venue 결정 maximum risk** — 분당 prospective 추가 시 Nat Comm reach 진짜 가능 (50-65%).

---

## 5. Decision rule for next prompt (Prompt 2 — Introduction)

### v3 confirm 사항 (5/4-5/10 W1 작업)

- [x] Title Hybrid Cand 1 (verb: reveals 권장)
- [x] Abstract 153 words v2 — Conclusions (b) 채택
- [x] Results 대안 A ordering 2.1 → 2.3 → 2.4a → 2.4b → 2.2 → 2.5
- [x] Section 2.4 split (2.4a epigenetic R5-2 + 2.4b immune-overlap teaser)
- [x] Pillar 5 main story Paper 2 reserve
- [x] Fig 8 (D) → S8 supp
- [x] **★ v3: Landa 2016 JCI cite (Krishnamoorthy 2025 Nat Comm misattribution catch)**
- [x] **★ v3: Discussion 3.1 dial-back framing ("identifies an upstream signature consistent with...")**
- [x] **★ v3: Reverse-causality 3-layer 차단 (Methods Yoo 2016 first cite + Disc 3.1 convergence + Cover Letter)**
- [x] **★ v3: 5 추가 cite (Yoo 2016, TCGA 2014, Pu 2021, Bradley 2010, Wirth 2020 NEJM)**
- [x] **★ v3: ATA 2015 + ATA 2025 dual cite Discussion 3.2**
- [x] **★ v3: Fig 1 5 panels → 4 panels (Heatmap → S1 supp)**
- [x] Limitations 4 main + 2 supp
- [x] Venue ladder Cell Rep Med 1순위 + JCI Insight + Nat Comm dual stretch + npj fallback

### 본인 read 일정 (5/4 부터)

- [ ] Hook Alt B refined vs hybrid "compass" 보존 — 본인 voice 결정
- [ ] Title verb 결정 (reveals / unmasks / resolves / defines)
- [ ] Landa 2016 JCI 직접 read (1-1.5 hr) — Discussion 3.1 framing 정확도
- [ ] Yoo 2016 PLOS Genet review (30 min) — Methods + reverse-causality 차단
- [ ] ATA 2015 PDF + ATA 2025 PDF (1.5 hr 합산) — institutional access
- [ ] Bradley 2010 citation verify (15 min)

---

## 6. Bug / risk flags (v3 update)

| Risk | Mitigation (v3) |
|---|---|
| Pillar 5 Paper 1 vs Paper 2 boundary 모호 | 2.4b 마지막 sentence "mechanistically distinct immune-overlap subtype warranting separate inquiry" → Paper 2 hook 명시. Paper 2 backbone reserve: K2 + Pan-Asian HLA n=874 + GSE286332 mediation + BCR/TLS + sub-B NBNR |
| K2 mini-index calibration FAIL (D7-P3) | Limitations 3.4 (iii) 정직 disclosure; K2 alternate evidence (94.6% DM2 + HLA n=874) 로 cover; S9 supp 추가 |
| Bundang prospective 0% | Limitations 3.4 (ii) 정직 disclosure; East-Asian 4,300+ tumors 로 cover; venue ladder "stretch" 표현 |
| ★ Krishnamoorthy 2025 Nat Comm misattribution catch | v3 정정: 진짜 cite는 Landa 2016 JCI 126(3):1052-1066 (Krishnamoorthy GP 9번째 공저자). Discussion 3.1 dial-back framing 적용 |
| ★ 8-gene 5/8 = Landa ATC silenced gene list overlap (reverse-causality risk) | 3-layer 차단: Methods Yoo 2016 first cite (panel origin) + Discussion 3.1 convergence framing + Cover Letter explicit statement ("panel predates Landa access") |
| ATA 2025 fusion/epigenetic 반영 여부 미확인 | Working assumption: 미반영 가능성 높음. 본인 PDF 검증 30 min 후 framing 조정 (시나리오 1/2/3) |
| Bradley 2010 citation verification 필요 | 본인 PDF 또는 PubMed search 15 min |
| ATA 2015 alignment claim 강도 | Discussion 3.2 dual cite (ATA 2015 + ATA 2025) specific paragraph + Wirth 2020 NEJM LIBRETTO-001 quantitative reference |
| 28-point limitations | 4 main + 2 supp condensed (web Claude 권장 유지) |
| Figure 1 5-panel density (Cell Press first impression) | v3 적용: 4 panels (Heatmap → S1 supp); A-Sankey + B-UMAP + C-Driver neutrality + D-Pan-genome ARI |
| Figure 7 5-panel + Fig 8 3-panel | Cell Press 4-6 per figure 표준 cap 내 |

---

## 7. References (v3) — 7 hard cites + 2 보조

### ★★★ Hard cites (Introduction + Discussion 필수)

| # | Citation | Cite 위치 | Role |
|---|---|---|---|
| 1 | **Cancer Genome Atlas Research Network. Integrated genomic characterization of papillary thyroid carcinoma. *Cell* 2014;159(3):676-690. PMID 25417114** | Intro 1.2 | TCGA-PTC canonical BRAF-like/RAS-like spectrum origin |
| 2 | **Yoo SK et al. Comprehensive analysis of the transcriptional and mutational landscape of follicular and papillary thyroid cancers. *PLOS Genet* 2016;12(8):e1006239** | Methods (panel origin) + Intro 1.2 (BRS 273-gene) + Disc 3.1 (convergence framing reverse-causality 차단) | 8-gene origin paper, **reverse-causality 차단 enabler** |
| 3 | **Xing M et al. Association between BRAF V600E mutation and recurrence of papillary thyroid cancer. *NEJM* 2014;371(15):1456-1457** (또는 *JCO* 2015) | Intro 1.3 | Dark matter concept origin |
| 4 | **Haugen BR et al. 2015 American Thyroid Association Management Guidelines... *Thyroid* 2016;26(1):1-133. PMID 26462967** | Intro 1.1 + Disc 3.2 | ATA 2015 risk stratification |
| 5 | ★ **Ringel MD et al. 2025 American Thyroid Association Management Guidelines for Adult Patients with Differentiated Thyroid Cancer. *Thyroid* 2025;35(8):841-985. PMID 40844370** | Intro 1.1 + Disc 3.2 dual cite | ATA 2025 update |
| 6 | ★ **Landa I et al. Genomic and transcriptomic hallmarks of poorly differentiated and anaplastic thyroid cancers. *J Clin Invest* 2016;126(3):1052-1066. PMID 26878173** | Intro 1.3 + Disc 3.1 | **PDTC/ATC dedifferentiation trajectory + 8-gene 5/8 overlap** (Krishnamoorthy GP 9번째 공저자) |
| 7 | **Wirth LJ et al. Efficacy of selpercatinib in RET-altered thyroid cancers. *NEJM* 2020;383(9):825-835. PMID 32846061** | Disc 3.2 | LIBRETTO-001 specific (FDA approval, ORR 79%) |

### ★ 보조 cites

| # | Citation | Cite 위치 | Role |
|---|---|---|---|
| 8 | Pu W et al. Single-cell transcriptomic analysis of the tumor ecosystems underlying initiation and progression of papillary thyroid carcinoma. *Nat Commun* 2021;12:6058. PMID 34663816 | Methods + Disc 3.3 | sc external validation (GSE184362 6 PTC) |
| 9 | Bradley CA et al. (2010, citation verify needed) — BRAF V600E immune escape thyroid | Disc 3.1 | counter-intuitive BRAF V600E HIGHER HLA-I (R3 finding 14) |
| 10 | Krishnamoorthy GP et al. RBM10 loss promotes metastases by aberrant splicing... *J Exp Med* 2025;222(5):e20241029. PMID 39992626 | Disc 3.4 limitations | Forward implication only (broad context) |
| 11 | Liu Z et al. (2017 Asian thyroid cohort) | Intro 1.3 + Disc 3.3 | East-Asian generalizability |
| 12 | Wang YL et al. (Wang 2024 *Endocr Connect* 13(11):e240301, PMID 39235852) | Disc 3.3 | East-Asian Wang/Fudan (D3-P4 audit context) |
| 13 | (optional) Pozdeyev N et al. Genetic Analysis of 779 Advanced Differentiated and Anaplastic Thyroid Cancers. *Clin Cancer Res* 2018;24(13):3059-3068. PMID 29615459 | Disc 3.1 보조 | 더 광범위 advanced cohort context (Landa 2016 보강 가능) |
| 14 | (optional) Nat Commun 2025 PDTC/ATC proteogenomic (s41467-025-58910-3) | Disc 3.1 supporting | Recent PDTC/ATC subtype taxonomy 보조 cite (본인 직접 read 안 함) |

### 본인 read 일정 (필수)

- [ ] Landa 2016 JCI 직접 read (1-1.5 hr) — ★★★ Disc 3.1 framing 정확도
- [ ] Yoo 2016 PLOS Genet review (30 min) — Methods + reverse-causality 차단
- [ ] ATA 2015 PDF Tables 11-15 verbatim (1 hr, institutional access)
- [ ] ATA 2025 PDF fusion/epigenetic 반영 여부 (30 min, institutional access)
- [ ] Bradley 2010 citation verify (15 min, PubMed)

---

## 8. ★ 5/4-5/10 W1 schedule (marathon, NOT sprint)

### 5/3 EOD 까지 (today + 3 days)
- [x] Citation correction (v3 outline + prep files) — ✅ 완료
- [ ] 본인 read: Landa 2016 JCI (1.5 hr)
- [ ] 본인 read: Yoo 2016 PLOS Genet (30 min)
- [ ] Memory entry: Landa save + schedule + Pillar 1 STRONG

### 5/4 (Mon)
- [ ] ATA 2015 PDF Tables 11-15 verbatim verify (institutional access)
- [ ] ATA 2025 PDF fusion/epigenetic 반영 여부 (institutional access)
- [ ] Bradley 2010 citation verify (15 min PubMed)
- [ ] Yu professor 미팅 약속 잡기 (분당 outreach 합의용)

### 5/5-5/10 W1 (manuscript writing 시작)
- [ ] Yu professor 미팅 → 분당 outreach 발송 합의
- [ ] v2.5 미세 조정 (Hook Alt B vs hybrid compass + Title verb 결정 + Fig 1 4 panels confirm)
- [ ] Manuscript v8 Prompt 1 (Title + Abstract + Outline) 본인 confirm
- [ ] **Prompt 2 (Introduction 600-900 words draft) 시작**

### D9, D10, D11 분석 sprint 없음

본인 4/30 약속 그대로 — marathon, sprint 안 함. Pillar 1 STRONG 확정 = 분석 끝.

---

## 9. v3 한 줄 정리

**v3 = web Claude 2차 review (7 권장 + 5 추가 cite + Landa misattribution 정정) 모두 적용. Hook Alt B refined + Landa 2016 JCI dial-back framing + Yoo 2016 reverse-causality 3-layer 차단 + 5 hard cites (Yoo/TCGA 2014/Pu 2021/Bradley 2010/Wirth 2020) + ATA dual + Fig 1 4 panels. 본인 read 1.5 hr (Landa) + ATA PDF 1.5 hr + Yoo 30 min = 3.5 hr prep → 5/5-5/10 W1 Prompt 2 (Introduction draft) 진입.**
