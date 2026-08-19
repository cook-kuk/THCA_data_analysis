# Paper 2 (HT-isolated) — outline skeleton — 2026-05-13

## Scope freeze (Yu 2026-05-04)
- HT-overlap PTC only. No GD/Graves. No TSAb/TSI/TED. No Chu 2018.
- GD/Chu forest → Paper 4/backlog (renumbered, 4/4 gating)
- Forbidden words (mirror v18_paper2_HT_isolated memory): Graves, GD, TSAb, TSI, TED, Chu 2018, thyroid-eye, TRAb-positive
- Audit gate: any draft text mentioning these → reroute to Paper 4 backlog

## Target venue ladder
- Reach: Nature Communications / Cell Reports Medicine
- Base: JCI Insight
- Fallback: npj Precision Oncology
- Decision criterion: image-DM1 Pillar III seed-sensitivity disclosure outcome (main vs supp) gates the reach-tier choice

## Title candidates (3)
1. "H&E pathology projection of a molecular dark-matter axis in Hashimoto-overlap papillary thyroid cancer"
2. "Three-pillar evidence (HLA-II, HT biology, image-DM1) for a Hashimoto-driven subclass of PTC"
3. "Foundation-model histology recovers a molecular DM1/DM2 axis in Hashimoto-overlap PTC and links it to HLA-II–mediated tumor immunity"

## Section structure

### 1. Introduction (400-600 words)
- 1.1 PTC heterogeneity + HT-overlap clinical question
  - HT-overlap PTC prevalence; clinical-management ambiguity; Korean overdiagnosis paradigm context (cite Paper 1 framing, do not duplicate)
  - **VOICE-PROTECTED Hook** — author keyboard only
- 1.2 Existing HT-PTC mechanistic limits
  - Bulk-only studies; missing pathology projection; HLA mediation untested at n>800
  - Reviewer Q hooks: why H&E now, why foundation features
- 1.3 Paper 1 → Paper 2 bridge: 8-gene panel + DM1/DM2 axis (Paper 1 cite only)
  - One-paragraph bridge; do not restate Paper 1 mechanism; cite Paper 1 preprint slot
  - Round 7 / Round 8 robustness numbers cited as "see Paper 1" only
- 1.4 Aim
  - Project DM1/DM2 onto H&E pathology + HT biology context + HLA-II mediation
  - Three-pillar architecture stated up front

### 2. Results (2,000-3,000 words)
- 2.1 Pillar I — HLA
  - 6-allele readiness matrix: HLA-A*02:01, A*24:02, B*44:02, B*44:03, DRB1*01:01, DQB1*02:01
  - kim2014 Korean reference panel forest; source-coverage heatmap; allele × source × n cells
  - n=874 HLA-II mediation analysis: HT biology → HLA-II → DM1 axis path
  - Honest boundary: K2 vs AFND ratios 0.93–1.22 all p>0.16 (mediation, not GWAS-style association)
- 2.2 Pillar II — HT biology
  - GSE286332 PTC vs PTC+HT n=18 main cohort
  - TLS d=+1.96; IGHV clonality d>0.5; AICDA up; IFN-γ FDR=2e-4
  - 4/4 sign consistency across TCGA / GSE286332 / Lu 2023 / Korean K2 projection
  - TCGA HT signature forest + Lu 2023 single-cell scatter as confirmatory layers
- 2.3 Pillar III — image-DM1
  - ViT-L + UNI foundation features + CLAM attention-MIL pipeline
  - N=59 (29 DM1 / 30 DM2) mean 5-fold CV AUC = 0.830 ± 0.139
  - Pooled cross-fold AUC = 0.746
  - Phase 1 GSE250521 spatial overlay |ρ|=0.610 max
  - GSE230424 pathology thyroid-axis scout (confirmatory)
  - **Residual seed-sensitivity disclosure** (AUC=1.000 audit artifact) — decision point §D2 below
- 2.4 Cross-pillar synthesis
  - HLA-II × HT-biology × image-DM1 convergence diagram
  - 84,549 labeled-record benchmark anchor (from manuscript_artifacts/draft.md)
  - Decision-support framing: H&E → DM1 score → HT-context interpretation

### 3. Discussion (800-1,200 words)
- 3.1 Three-pillar synthesis — **VOICE-PROTECTED** author keyboard only
- 3.2 Image-DM1 as decision-support layer (not standalone diagnosis)
  - Explicit non-claim: not a clinical diagnostic; not replacement for IHC/molecular panel
  - Decision-support framing aligned with Paper 1 Reviewer Q3 safety net
- 3.3 GSE248205 negative control caveat
  - AITD spatial → NO_GO_LOCAL_HE_CAVEAT_ONLY
  - Pooled rho=0.590 fails sample-centered; treat as honest negative
  - Reviewer-defense use only
- 3.4 Limitations (4 paragraphs) — **VOICE-PROTECTED** author keyboard only
  - Para A: Cohort size (N=59 image, n=18 GSE286332 main) — honesty paragraph
  - Para B: HLA mediation directionality (cross-sectional, not causal)
  - Para C: Pillar III seed-sensitivity; residual AUC audit
  - Para D: HT-isolation scope (GD/TSAb explicitly out; Paper 4 forward reference)

### 4. STAR Methods + Supplementary
- 4.1 HLA imputation (arcasHLA pipeline; reference: v17_arcasHLA_korean_k2 memory)
- 4.2 HT biology cohort assembly + DEG protocol (GSE286332 n=18 anchor)
- 4.3 image-DM1: ViT-L + UNI features + CLAM attention-MIL hyperparameters + seed schedule
- 4.4 Spatial overlay: GSE250521 / GSE230424 / GSE248205 protocols
- 4.5 84,549-record benchmark sheet (from manuscript_artifacts/draft.md)
- 4.6 Supplementary figure plan: extended HLA forest, GSE286332 TLS gallery, image-DM1 attention maps, GSE248205 negative-control board

## Voice-protected sections (mirror Paper 1 rule)
- Hook §1.1
- Discussion §3.1 + §3.4
- Cover letter para 1
- Reviewer Q&A author opinion section
- Aim sentence §1.4 (final author phrasing only)

## Decision points for Yu (3)
1. **Scope freeze 재확인** — HT-isolated 유지? (현재 v18_paper2_HT_isolated 메모리 기준 freeze 상태; Yu가 GD/Chu 재포함 원하시면 Paper 4 backlog 해제 + 4/4 gating 재검토 필요)
2. **Pillar III image-DM1 seed-sensitivity 본문 vs supp 배치** — residual AUC=1.000 audit를 §2.3 본문 단락으로 정직 노출할지, Supp Note로 빠질지. 본문 노출 = reviewer trust 상승 + reach venue 안전; supp 이동 = narrative cleaner. (권고: 본문 단락 1개 + supp 보강 figure)
3. **Cell Rep Med vs JCI Insight 1순위** — Cell Rep Med은 H&E foundation-model narrative 강함 (image-DM1이 주연); JCI Insight는 HT biology + HLA-II mediation 강함 (Pillar I+II가 주연). Pillar 가중치 결정 필요.

## Asset map
- HLA: project/papers_hub_2026_05_04/assets/paper2_hla/
- HT: project/papers_hub_2026_05_04/assets/paper2_ht/
- image-DM1: project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/
- path2space: project/papers_hub_2026_05_04/assets/paper2_path2space/
- he_overlay_viewer: project/papers_hub_2026_05_04/assets/paper2_he_overlay_viewer/ (2026-05-11)
- moonshot_editorial_pack: project/papers_hub_2026_05_04/assets/paper2_moonshot_editorial_pack/ (2026-05-10)
- external_validation_ready pack: project/papers_hub_2026_05_04/assets/paper2_external_validation_ready/
- manuscript draft: project/papers_hub_2026_05_04/manuscript_artifacts/draft.md (internal review, 84,549 labeled records)

## Cross-paper references
- Paper 1 (DM1 molecular dark matter) — cite only; do not duplicate mechanism
- Paper 3 (ICI vulnerability dark thyroid cancer) — Track A FROZEN; forward-ref only if Pillar II HT-immune axis intersects ICI vulnerability discussion (avoid scope creep)
- Paper 4 backlog (Korean GD HLA Pan-Asian) — explicit forward reference in §3.4 Para D

## Next steps after Yu meeting (2026-05-14)
- Lock decisions D1/D2/D3 → freeze §02_section_anchors.md
- Author-keyboard Hook + Aim + Limitations drafts
- Claude scaffolds §2.1/§2.2/§2.3 results prose from existing asset captions
- Internal review checkpoint before bioRxiv submission gate
