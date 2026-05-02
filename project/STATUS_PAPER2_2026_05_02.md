# Paper 2 Status — Hashimoto-overlap PTC mechanism
**Snapshot: 2026-05-02** · Owner: Seungho Cook

## One-line
**5/4 Yu 결정으로 scope re-isolated**: Paper 2 = Hashimoto-overlap PTC ONLY. GD/Graves' phenotype 는 Paper 3 영역. Pillar I v2 forest (Korean PTC vs **AFND South Korea baseline**) committed. Marathon 5/4 시작.

## 5/4 Scope reformulation (Yu advisor)
- **Before (v1, 5/3 STRONG)**: Korean PTC pool n=874 vs **Chu 2018 Han Chinese GD** random-effects forest. "Pan-Asian autoimmune-thyroid susceptibility continuum" framing.
- **After (v2, 5/4)**: Korean PTC pool n=874 vs **AFND South Korea baseline** (sample-weighted pool; China Harbin Korean fallback; Lee 2014 Tissue Antigens reference). HT-overlap mechanism only.
- Chu 2018 GD forest → **Paper 3 reserve** (`chinese_GD_vs_chinese_ctrl_forest_paper3.json`)
- Reason: HT (destructive infiltration → dedifferentiation) ≠ GD (TSAb hyperthyroidism). 같은 paper 에서 conflate 시 reviewer scope confusion.

## Forbidden in Paper 2 (등장 시 STOP + 재작성)
- GD-specific: `Graves`, `GD`, `TSAb`, `TSI`, `hyperthyroidism`, `thyrotoxicosis`, `exophthalmos`, `Graves' ophthalmopathy`, `TED`, `thyroid eye disease`
- Cancer-deep-dive: BRAF V600E mechanism, TERT promoter kinetics, DM1 76.8% fusion paradigm (Paper 1 main; Paper 2 supplementary only)
- Conflations: `autoimmune-PTC` → `Hashimoto-overlap PTC`; `autoimmune-thyroid axis` → `Hashimoto-thyroid overlap axis`; `Pan-Asian autoimmune-thyroid susceptibility` → `Pan-Asian thyroid HLA susceptibility (HT context)`

## Done
- **Pillar I v2 forest** committed: `p2_pillar1_forest_v2/PILLAR1_FOREST_V2_SUMMARY.md`, forest PDF/PNG/JSON
- **Terminology audit** (Task B): `terminology_correction_2026_05_04/before_after_diff.md`
- v1 script DEPRECATED banner added: `v17_paper2_pillar1_forest.py`
- v2 script: `v17_paper2_pillar1_forest_v2.py`
- GSE286332 (Korean PTC vs PTC+HT n=18): 10,380 DEGs, HLA-II d=+3.65, IFN-γ FDR=2e-4
- TCGA Hashimoto-like generalization: GSE286332 sig → TCGA DM2 enrichment OR up to 5×, p=6e-10
- BCR clonal + TLS (D5-P6): TLS d=+1.96, IGHV clonality d>0.5, AICDA up — antigen-driven B cell response
- DM1 sub-B = NBNR cluster (D6-P7): n=56, 96% mutation-neg, 4× Hashimoto-like rate

## Marathon Task A/B/C (5/4 → 6/13)
- Task A: Pillar I v2 forest — **DONE** (5/4 commit)
- Task B: terminology audit — **DONE** (5/4 commit)
- Task C: ?? (memory unclear, 진행 시 확인 필요)
- Voice-protected: Hook · Aim · Disc 3.1 · Limitations · Cover Para 1 · Q9 (본인 키보드)

## Cross-paper reciprocal Discussion (1 line)
"Korean PTC HLA susceptibility (this paper) shares background with Korean GD (Paper 3, Cook et al. in prep), but PTC ≠ GD as diseases — PTC = neoplasia, GD = autoimmune hyperthyroidism, HT = autoimmune destructive infiltration. Pan-Asian thyroid HLA axis is shared, disease mechanisms distinct."

## Venue
- 미상 (TBD Yu 미팅) — 가능 후보: JCI Insight, Cell Rep Med (HT-overlap mechanism + TLS + IGHV + TCGA generalization 패키지)

## Key files
- `project/manuscript_p2_brief/` — brief HTML/PDF, README, COHORT_ACCESS_GUIDE
- `project/results/p2_pillar1_forest_v2/` — Pillar I v2 anchor
- `project/results/terminology_correction_2026_05_04/` — Task B audit
- `project/notebooks_or_scripts/v17_paper2_pillar1_forest_v2.py` — reproducible
- Memory: `v18_paper2_HT_isolated`, `v17_paper2_pillar1_forest_strong` (deprecated framing), `v17_GSE286332_strong_go`, `v17_D4P2_tcga_hashimoto_generalization`, `v17_D5P6_BCR_clonal_TLS`, `v17_D6P7_dm1_subB_NBNR`
