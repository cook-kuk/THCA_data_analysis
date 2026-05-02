# Paper 3 Status — Korean Graves' disease HLA, Pan-Asian
**Snapshot: 2026-05-02** · Owner: Seungho Cook

## One-line
**DOCK** — 4/4 진입 조건 미충족. 분당 SNUH 응답 또는 6주 freeze (~6/15) 까지 대기. 본 session 진입 X, Paper 1/2 marathon 그대로.

## 진입 조건 (4/4 필수, 현재 0/4)
| # | 조건 | Status |
|---|---|---|
| 1 | Paper 1 bioRxiv 제출 완료 | ❌ target 6/13 |
| 2 | Paper 2 Task A/B/C 완료 | 🟡 A/B done (5/4), C 미상 |
| 3 | 분당 SNUH Graves' n>50 응답 | ❌ outreach 무응답 |
| 4 | Yu Paper 3 trajectory 합의 | ❌ |

## Disease 정의 (advisor 인용)
"갑상선을 자극하는 자가항체가 생성되어 갑상선 호르몬이 과다 분비되는 자가면역 질환"  
GD = TSAb-mediated hyperthyroidism. HT 의 destructive infiltration 와 mechanism **정반대**.

## Forbidden in Paper 3 (등장 시 STOP)
- HT 영역 (Paper 2): `Hashimoto`, `HT`, `lymphocytic thyroiditis`, `TLS`, `tertiary lymphoid structure`, `AICDA`, `BCR clonality`, `PTC+HT`, `Hashimoto-overlap`
- Cancer 영역 (Paper 1): BRAF V600E, RET fusion, DM1/DM2 cluster, 8-gene RAI panel mechanism, LIBRETTO-001, selpercatinib
- Conflations: `autoimmune-PTC`, `PTC+GD overlap`, `thyroid neoplasia` framing

## 5-Pillar (GD-only)
| Pillar | Source | Status |
|---|---|---|
| P1 Korean GD HLA cohort | Bundang prospective n>50 | 데이터 없음 (대기) |
| P2 Chu 2018 Han Chinese GD replication | n=1,468 | starting point ready (Paper 2 v1 forest 이전) |
| P3 Pan-Asian random-effects forest meta | DerSimonian-Laird | P1+P2 confirm 후 |
| P4 cookHLA SNP-based 4-digit imputation | Bundang SNP if available | 데이터 없음 |
| P5 GD-specific HLA haplotype | DRB1-DQA1-DQB1 trio + DPB1 LD | downstream |

## Tasks (진입 후 4-6주)
- T1: Bundang Graves' arcasHLA → `project/results/paper3_bundang_GD_HLA/` (1주)
- T2: Pan-Asian forest meta → `project/results/paper3_panasian_forest/` (1-2일)
- T3: cookHLA SNP cross-platform validation → `project/results/paper3_cookHLA_validation/` (1주, SNP 가능 시)
- T4: Outline (본인 voice 키보드, Claude generate X)

## Venue commit ladder (Bundang n 기반, Yu 미팅 직전)
- **n > 50 STRONG GO**: J Autoimmunity (IF 12-14) reach default
- **n = 30-50 BORDERLINE**: J Autoimmunity 시도 + Front Immunol (IF 5-7) parallel-ready, cover letter 두 버전
- **n < 30 또는 6주 freeze**: J Autoimmunity 포기. Front Immunol 또는 HLA (IF 3-4) 직행, Bundang supplementary cohort demote
- Commit timing: T1 + T2 완료 후 (n + effect size 확인), Yu 미팅 안건

## Cross-paper reciprocal Discussion (1 line, 본인 voice)
"The Pan-Asian HLA susceptibility axis identified here for Graves' disease shares background alleles with Korean PTC (Paper 2, Cook et al. in prep) and HT-overlap PTC, but mechanism diverges — GD = TSAb-driven hyperthyroidism, HT = destructive lymphocytic infiltration, PTC = thyroid neoplasia. Shared HLA risk reflects population-level susceptibility; disease specificity is determined downstream."

## Source files (진입 후 read-only)
- `project/results/p2_pillar1_forest/chu2018_allele_summary.tsv` (Paper 2 → Paper 3 이동)
- `project/results/p2_pillar1_forest/chinese_GD_vs_chinese_ctrl_forest_paper3.json` (Paper 2 v1 산출물, P3 starting point)
- Bundang Graves' arcasHLA results (응답 후 generation)
- cookHLA Nat Commun 본인 paper (Pan-Asian extension unique angle)

## NOT 진행
- Paper 1/2 file 직접 modify 금지
- HT mechanism (TLS, BCR, AICDA) 분석 금지
- PTC+GD overlap phenotype 분석 금지
- Manuscript section generate 금지 (본인 voice)

## Memory anchor
`v19_paper3_GD_gating.md` (gating + scope), `v18_paper2_HT_isolated.md` (Chu 2018 forest 이전 근거), `v17_K2_vs_bundang_distinction.md` (K2 ≠ Bundang)
