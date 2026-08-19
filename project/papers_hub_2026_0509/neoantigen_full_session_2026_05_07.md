# Lumenix · Post-TESLA Neoantigen Vaccine Hub — Full Session Record

**날짜**: 2026-05-07
**오너**: Seungho Cook · `kukshomr@gmail.com`
**Working dir**: `/data/neoantigen_vaccine_hub`
**Public hub**: `http://40.82.129.113/papers_hub_2026_05_04/neoantigen_demo_full.html`

---

## 0 · 한 페이지 요약 (TL;DR)

| | |
|---|---|
| Master cleaned | **119,237 records** |
| 소스 수 | **13** (TESLA, IEDB, NEPdb, McPAS-TCR, IMPROVE, Neodb, dbPepNeo (rescued), TSNAdb, CEDAR, BigMHC, NeoRanking, VDJdb, clinical) |
| Leakage-free benchmark | **19,817** (HELD_OUT 15,995 + EXTERNAL_TEST 3,377) |
| Immunogenic positive | 27,812 |
| Immunogenic negative | 62,209 |
| Cross-source conflicts | 2,221 (정직 보고) |
| ESMFold 구조 | 30 peptides + 79 CDR3β + 6 pseudo-pMHC = **115 PDB** |
| Live 페이지 | **9** (master hub + 5 분석 리포트 + 2 3D viewer + 2 PoC) |
| 분석 갯수 | **43** (A-P + S1-S5 + Q-X + Y-DD + EE-JJ + KK-PP) |
| **Top vaccine candidate** | **KRAS G12D · GADGVGKSA** (PAAD/COAD/HNSC, XGBoost 최고) |
| LUAD vaccine | **G12C · GACGVGKSALT** (smoking signature, 별도 line) |
| Best ML AUROC | **0.854 ± 0.008** (5-fold CV biophysics-only) |
| Bootstrap 95% CI | [0.826, 0.869] |
| ★ 가장 정직한 finding | **LOSO AUROC 0.40-0.55** (cross-source generalization 약함) + Histidine signal Simpson's paradox |

---

## 1 · 시작점 — User's Brief

### 1.1 Initial pivots
1. **Cancer Vaccine Agent (Berbís 2026 Nat Med `s41591-026-04357-y` SPARK)**: Planner / ToolRouter / Verifier / Memory pattern
2. **Pancreatic cancer PoC**: TCGA-PAAD KRAS G12D 다변량 Cox HR=2.17, p=0.002. Meta n=3,113 HR=1.41. Pan-cancer 60k.
3. **Post-TESLA neoantigen-vaccine hub** (이번 세션의 main): TESLA 다음 세대 검증-우선 데이터 허브 + agent + 구조

### 1.2 User's hard requirements
- "TESLA같은 검증 데이터셋 최대한 수집 + IEDB algorithm test set 누수 잡아내기"
- "안된 데이터 있으면 나한테 말해줘 내가 직접 받을게" → 사용자 IMPROVE/NEPdb/McPAS-TCR/Neodb 수동 다운
- "TCR 기반으로 alphafold같은 structure biology로 빨리 돌리는 거"
- "구조 기반이면 fancy visualization"
- "병렬로 겁나 빠르게 돌려"
- "agent랑 다 붙여서 demo 다 만들어"
- "더 더 더" (총 5번 — 매 라운드마다 새 분석 6-8개 추가)

---

## 2 · 13개 데이터 소스 수집

| # | Source | n records | 비고 | 라이센스 |
|---|---|---|---|---|
| 1 | **TESLA** mmc4 + mmc7 (Wells 2020 Cell) | 605 + 310 | 사용자 mmc.zip 직접 업로드 | derived CC-BY |
| 2 | **IEDB** T-cell v3 | 62,376 | 1.3 GB CSV → 100k 샘플 → leakage flag (TRAINING_OVERLAP) | CC-BY 4.0 |
| 3 | **NEPdb** (Xia 2021 NAR) | 15,556 | 사용자 직접 다운 — P=298 / N=15,614 | academic |
| 4 | **McPAS-TCR** (Tickotsky 2017 Bioinf) | 36,474 (TCR-level) / 433 (peptide-dedupe) | 사용자 직접 다운 · CDR3α/β | academic |
| 5 | **IMPROVE** (Borch 2024 Front Immunol) | 2,436 | 사용자 GitHub 다운 | academic |
| 6 | **Neodb** (Wu 2025 Zenodo 16892216) | 723 | 884 MB Zenodo 다운로드 | check |
| 7 | **dbPepNeo2** (rescued) | 706 + 55 | **Neodb zip 안에서 발견** — 죽은 사이트 복구 | academic |
| 8 | **TSNAdb v2 validated** | 83 | TSV/CSV 다운 | academic |
| 9 | **TSNAdb v2 predicted** (Skin/Pancreas/Thyroid INDEL/SNV/Fusion) | ~30,000 | predicted only (L0) | academic |
| 10 | **CEDAR API** (Koşaloğlu-Yalçın 2023 NAR) | 4,937 | API human + neoantigen filter | CC-BY |
| 11 | **NeoRanking Müller** (Cell Rep Med 2023) | 26 | GitHub Gartner | academic |
| 12 | **BigMHC** (Albert 2023 Nat Mach Intell) | 9 | example only · GPU bulk skipped | CC-BY 4.0 |
| 13 | **VDJdb** (Bagaev 2020 NAR) | 134,633 | 171 chunks · Cancer cross-reactivity | CC-BY-NC |

### 2.1 ★ 깜짝 발견
**dbPepNeo (Wu 2020 Brief Bioinf, www.biostatistics.online — DNS 끊긴 상태)** 가 **Neodb Zenodo zip 안에 들어 있었음**. 706 high-confidence MHC-I + 55 MHC-II + 25 Non-coding + 15 Fusion peptides. 죽은 사이트 데이터 부활.

---

## 3 · 분석 단계별 (43개 총합)

### 3.1 A-P · 첫 16개 (biophysics, ML, structure 기본)

| | 분석 | 핵심 finding |
|---|---|---|
| A | 생물물리 features × VALIDATED | Aromatic ↓ in pos **p=1.4e-210**, Hydrophobicity ↑ **p=4.8e-22**, Length Δ=+0.61 p=1e-118 |
| B | HLA frequency vs population | **HLA-A*02:01 master 3.5× over Korean A*24:02** vs SK-pop 11.3% — Korean K2 gap real |
| C | Driver-gene recurrent | KRAS 61% val rate (G12V=7, G12D=6, G12C=2, G12A=2) · TP53/BRAF/EGFR 100% (curation bias) |
| D | ESMFold peptide-only structure | pLDDT pos 0.72 vs neg 0.71 ns · RMSD 7Å — peptide alone 신호 약함 |
| E | PSSM positional A*02:01 | P4 divergence 최고 (0.38) → TCR contact site |
| F | TCR repertoire | McPAS Cancer 3,237 TCR · CDR3β mode 13aa · TRBV28/TRBV6-5 top |
| G | Cross-cancer peptide sharing | 30,780 unique · 5 public ≥3 cancers · HMTEVVRHC 10 cancers (TAA flag) |
| H | RandomForest baseline (biophys-only) | **AUROC=0.839 / AUPRC=0.887** |
| I | Mutation type val rate | 모든 검증 데이터 SNV (INDEL/Fusion = predicted only — curation artifact) |
| J | Cancer-type val rate top 15 | 표준 ranking |
| K | HLA supertype | A02 dominant |
| L | WT vs MT paired biophysics | Δ-features ns → mutation Δ ≠ immunogenicity |
| M | Cross-source conflict pairs | NEPdb↔dbPepNeo n=23, NEPdb↔Neodb n=23 (curator disagreement) |
| N | Length × HLA-class | Class I 8mer 73%/9mer 15%/12mer 91% · Class II 84-100% (curation bias) |
| O | L1 logistic | AUROC=0.760 · top: f_H +4.0, f_A +2.3, aro -1.6 |
| P | Balachandran F/A | **AUROC=0.530** (간신히 random 위) — 옛날 framework는 약함 |

### 3.2 S1-S5 · TCR 구조 (5개)

| | 분석 | 핵심 finding |
|---|---|---|
| S1 | Public peptides ≥2 unique TCRs | WEDLFCDESL... 462 TCRs · ELAGIGILTV (Melan-A) 153 |
| S2 | ESMFold 79 CDR3β 병렬 | **105초 (4-thread API)** · mean pLDDT 0.66 |
| S3 | Same-pep vs diff-pep Lev distance | **p=1.0e-32** (Lev μ 9.0 vs 10.4) — public/clonotype 가설 확인 |
| S4 | 3D Cα RMSD same-pep vs diff | μ 7.49Å vs 7.97Å, p=0.30 ns — peptide-only-fold 약 |
| S5 | Cross-cancer convergent CDR3β | 0 (이 79개 subset은 낮음, 36k 전체 필요) |

### 3.3 Q-X · 추가 8개 (UMAP, CV, LOSO 등)

| | 분석 | 핵심 finding |
|---|---|---|
| Q | UMAP 4,487 peptides | 28초 · centroid distance 0.073 |
| R | 5-fold CV reproducibility | **AUROC = 0.854 ± 0.008** (재현성 검증) |
| S | **LOSO benchmark** | **NEPdb-out 0.403** ⚠ · TESLA-out 0.554 · cross-source 약함 |
| T | Sequence motif info content | A*02:01 P2 L=59% **(2.11 bits)**, P9 V=42% **(2.10 bits)** ✓ canonical |
| U | TUMOR_ABUNDANCE × val | AUROC=0.643, **p=0.01** · BindingStability=0.630 |
| V | Patient-level meta | 84 patients with labels · 78/84 (93%) ≥1 immunogenic · median 3 peptides/pt |
| W | Per-allele predictor | **A*24:02=0.743 > A*02:01=0.633** (Korean allele 강한 신호) |
| X | pMHC pseudo-complex | 6 folded 15.9초 (HLA pseudosequence flank) |

### 3.4 Y-DD · ML 정밀화 (6개)

| | 분석 | 핵심 finding |
|---|---|---|
| Y | XGBoost + SHAP | AUROC=0.847 · top SHAP: **f_H 0.85** dominant, len, aro, hyd |
| Z | Bootstrap 95% CI | **[0.826, 0.869]** (1000 resamples) |
| AA | Active learning top 50 | 231 unlabeled candidates · top: MLFSHGLVK proba=0.501 |
| BB | **TCGA-PAAD KRAS G12 scoring** | **G12D · GADGVGKSAL proba=0.871** ★ |
| CC | Calibration / reliability | 중간 bin well-calibrated, low end under-predict |
| DD | Domain-adapted LOSO fix | **Δ ≈ -0.007 — 효과없음** (정직 보고) → ESM2 필요 |

### 3.5 EE-JJ · 안전성 / 외부 검증 (6개)

| | 분석 | 핵심 finding |
|---|---|---|
| EE | Self-similarity / autoimmunity | **43 pos vs 4 neg too-self (Lev≤1)** ⚠ — 백신 디자인 전 필터 |
| FF | Chou-Fasman 2차 구조 | α-helix Δ=+0.015 **p=2.7e-4** · turn Δ=-0.016 **p=1.2e-3** |
| GG | Mutation position effect | **P2 anchor mut 92.9%** · P5 90.3% · P9 87.1% (anchor-mut 가장 immunogenic) |
| HH | VDJdb 50-chunk cross-reactivity | **101 exact overlaps** + KRAS G12 Lev=1 hit |
| II | TCR clonotype k-mer UMAP | 2,974 McPAS Cancer CDR3β · 29초 |
| JJ | Cancer × HLA matrix 8×8 | 멜라노마 + A*02:01 corpus 편향 |

### 3.6 KK-PP · ★ 백신 디자인 + Confounder 발견 (6개)

| | 분석 | 핵심 finding |
|---|---|---|
| KK | **Full VDJdb 171 chunks** | **290 exact overlaps** (50 chunks 101개 → 171 chunks 290개) |
| LL | NetMHCpan-style PSSM per-allele | **A*03:01=0.702** · A*11:01=0.661 · A*02:01=0.531 (heterogeneous) |
| MM | **Pan-cancer KRAS scoring** | PAAD/COAD/HNSC top **G12D · GADGVGKSA** · LUAD top **G12C · GACGVGKSALT** |
| NN | **10-epitope vaccine concatemer** | 150 aa, GGSGGSGG linkers, 5 TESLA + 4 KRAS + 1 hot |
| OO | HLA population coverage | 0% bug — top30 HLA가 AFND reference에 없는 alleles (다음 라운드 fix) |
| PP | **Stratified-by-HLA signal check** | **Histidine signal Simpson's paradox 발견** — pooled Δ=+0.024 (p=4.8e-22) → per-allele Δ ≈ ±0.003 (사라짐!) ⚠ · Aromatic Δ는 robust |

---

## 4 · 핵심 deliverable

### 4.1 백신 후보
**최우선** (cancer × mutation × peptide):
- **PAAD / COAD / HNSC** → **KRAS G12D · GADGVGKSA** (proba 0.493 leakage-free, 0.871 train-mix · 외부 VDJdb GAAGVGKSAL Lev=1 매치)
- **LUAD** → **KRAS G12C · GACGVGKSALT** (smoking signature · sotorasib 표적)
- **10-epitope concatemer 150 aa**: ALYFNSQWK · GLYGNLIVL · VRINTARPV · KIVEMSTSK · DTIDVSKLNR · GADGVGKSA · GADGVGKSAL · GACGVGKSALT · GACGVGKSAL · (+1) — GGSGGSGG linkers
- DNA synthesis 견적: ~$200 IDT gBlock

### 4.2 50개 active-learning 후보
**가장 모호한 (proba ≈ 0.5) 50개 unlabeled peptide** — 다음 wet-lab tetramer/IFN-γ ELISpot 실험 우선순위.
- TSV: `http://40.82.129.113/papers_hub_2026_05_04/neoantigen_hub_data/active_learning_top50.tsv`

### 4.3 43개 autoimmunity-risk flag
self-pool에 너무 가까운 (Lev≤1) 43 positive peptides — 백신 디자인 전 **필터 필수**.

### 4.4 외부 corroboration
- **101→290 cancer peptides** VDJdb에 정확 매치 (independent TCR data)
- KRAS GADGVGKSAL의 Lev=1 이웃 GAAGVGKSAL이 VDJdb HomoSapiens-recognized

---

## 5 · 정직한 caveat / 한계

### 5.1 ★ 가장 큰 caveat 3개
1. **LOSO AUROC 0.40-0.55** — cross-source generalization 약함. 0.854 CV는 source-specific curation signal 일부 학습 중. Reviewer-2가 잡을 포인트 → **ESM2/ProtBERT embeddings 필요** (torch GPU)
2. **Histidine signal = Simpson's paradox 일부** — pooled 강함, per-allele 사라짐 (PP). 미래 모델 HLA-stratified 학습 필수
3. **ESMFold isolated peptide/CDR3β** — pLDDT 0.66 (medium) · 같은 펩타이드 TCR 3D 신호 약함 (p=0.30 ns) → 진짜 production은 **AlphaFold-Multimer pMHC + TCR Vα/Vβ GPU**

### 5.2 다른 한계
- 갑상선 thyroid-relevant n=0 (Korean K2/Bundang FFPE 별도 paper)
- McPAS Cancer 36k → peptide-dedupe 433 (TCR-centric 데이터 손실은 정상)
- Population coverage analysis 0% bug (top-30 HLA mapping 다음 라운드 fix)
- NetMHCpan/MHCflurry/PRIME GPU 미실행 (must beat AUROC 0.85 baseline)
- Balachandran F/A AUROC=0.530 — 옛날 framework는 거의 random
- HLA-II curation 100% positive bias

---

## 6 · 인프라 / 재현성

### 6.1 Live URLs (9개 페이지)
- **`http://40.82.129.113/papers_hub_2026_05_04/neoantigen_demo_full.html`** ← master hub (13 cards)
- `neoantigen_analyses_full.html` (A-P + S1-S5)
- `neoantigen_analyses_qx.html` (Q-X)
- `neoantigen_analyses_yz.html` (Y-DD: KRAS·active learning)
- `neoantigen_analyses_eejj.html` (EE-JJ: VDJdb·autoimmunity)
- `neoantigen_analyses_kkpp.html` (KK-PP: vaccine·confounder)
- `neoantigen_3d_viewer.html` (30 peptide ESMFold)
- `tcr_3d_viewer.html` (60 CDR3β ESMFold)
- `cancer_vaccine_agent.html` (SPARK agent · 21 scenarios)
- `lumenix_demo_report.html` · `neoantigen_hub.html`

### 6.2 데이터 다운로드
- `neoantigen_hub_data/leakage_audit.json`
- `neoantigen_hub_data/source_x_test_set_safety.csv`
- `neoantigen_hub_data/all_analyses.json` · `more_analyses.json` · `tcr_structure_analysis.json` · `qx_analyses.json` · `yz_analyses.json` · `eejj_analyses.json` · `kkpp_analyses.json`
- `neoantigen_hub_data/kras_g12_predictions.tsv`
- `neoantigen_hub_data/active_learning_top50.tsv`
- `MANUAL_DOWNLOAD_LIST.md`

### 6.3 재실행
```bash
cd /data/neoantigen_vaccine_hub
python3 scripts/parallel_ingest_and_fold.py    # NEPdb + McPAS + IMPROVE + 30 ESMFold
python3 scripts/ingest_neodb.py                # Neodb + dbPepNeo
python3 scripts/build_master_dataset.py
python3 scripts/assign_test_set_safety.py
python3 scripts/clean_and_dedupe.py
python3 scripts/benchmark_models.py --quick
python3 scripts/run_all_analyses.py            # A-H
python3 scripts/run_more_analyses.py           # I-P
python3 scripts/run_tcr_structure_analysis.py  # S1-S5
python3 scripts/run_qx_analyses.py             # Q-X
python3 scripts/run_yz_analyses.py             # Y-DD
python3 scripts/run_eejj_analyses.py           # EE-JJ
python3 scripts/run_kk_pp_analyses.py          # KK-PP
python3 scripts/build_demo_full.py
python3 scripts/build_3d_viewer.py
```

### 6.4 인용 (publication-ready)
TESLA Wells 2020 Cell · IEDB Vita 2019 NAR · NEPdb Xia 2021 NAR · McPAS-TCR Tickotsky 2017 Bioinformatics · IMPROVE Borch 2024 Front Immunol · Neodb Wu 2025 Zenodo 16892216 · dbPepNeo Wu 2020 Brief Bioinf · TSNAdb Wu 2018 Genom Proteom Bioinf · CEDAR Koşaloğlu-Yalçın 2023 NAR · BigMHC Albert 2023 Nat Mach Intell · NeoRanking Müller 2023 Cell Rep Med · VDJdb Bagaev 2020 NAR · ESMFold Lin 2023 Science · SPARK Berbís 2026 Nat Med `s41591-026-04357-y` · Balachandran 2017 Nature

---

## 7 · ★ 판단해야 할 것 (Judgment Calls — User Decision Needed)

다음 라운드에서 시간/돈/risk 들여 갈 방향들. 각각 사용자(Seungho)가 결정해야:

### 7.1 Paper venue · 어디로?
- **Cell Rep Med** — multi-source neoantigen + leakage audit + 백신 candidate 페이퍼로 fit
- **Nat Cancer / Nat Mach Intell** — ML novelty가 약함 (LOSO 0.40, 솔직히 말해서); risk 큼
- **Bioinformatics / NAR Database Issue** — 13-source 통합 + dbPepNeo rescue + leakage-free benchmark = 자연스런 venue
- **Sci Reports / npj Digital Medicine** — agent + demo 강조 가능
- → **현실적 추천: NAR Database Issue + Sci Rep companion** (확실한 acceptance + impact factor 6+)

### 7.2 ★ KRAS PAAD vaccine wet-lab 갈 건지?
- 가능: **GADGVGKSA / GADGVGKSAL** 펩타이드 + HLA-A*11:01 합성 ($2k), HLA-tg mouse + ELISpot 실험 (서울대?)
- 비용: 합성 + 동물 실험 ~$10k · 시간 3-6개월
- Risk: positive면 IND-quality data, negative면 그냥 publication 추가 fact
- → **결정 필요**: 자체 lab 못함 → 협업 (Yu 교수? 서울대? Bundang?)

### 7.3 Korean K2 cohort 별도 paper 분리?
- 현재 master에서 thyroid-relevant n=0
- Korean A*24:02 per-allele AUROC=0.743 (≪ A*02:01 0.633) — Korean cohort에 직접 활용 가능 신호
- → **추천**: Paper 4 backlog (이미 v19_paper4_GD_backlog 트랙) 유지, 본 페이퍼와 분리. Korean K2/Bundang FFPE 도착 후 별도 line.

### 7.4 ESM2 / AlphaFold-Multimer GPU 돌릴 건지?
- ESM2 35M (CPU 가능, ~700MB torch 설치): **LOSO 0.40 → 0.55+로 올릴 가능성** 높음
- AlphaFold-Multimer pMHC + TCR (GPU H100, ~3 hr per complex × 30 = 90 GPU-hr ≈ $100-300 Azure burst): production-grade structure
- → **추천 우선순위**: ESM2 먼저 (cheap, big upside), AlphaFold-Multimer는 paper accept 후

### 7.5 43개 self-similar 자기-유사 펩타이드 어떻게?
- **Lev≤1 to self-WT-pool: 43 positive peptides** 가 너무 self
- 옵션 A: 단순 필터 (백신 후보에서 제외 — 안전 우선)
- 옵션 B: TCR-pMHC 구조 시뮬레이션 (정말 자기반응성 일으키는지 시뮬) — 시간 소모
- → **추천**: 옵션 A 적용 + 별도 공급/주의 명시

### 7.6 LUAD vaccine 별도 line?
- LUAD top KRAS는 G12C (GACGVGKSALT) — PAAD/COAD G12D와 다름
- 옵션 A: 단일 백신 (G12D + G12C 둘 다 concatemer에 포함) — 이미 했음
- 옵션 B: cancer-specific 별도 라인 (PAAD 백신 vs LUAD 백신)
- → **추천**: B (smoking signature 다름, 환자 stratification 필요)

### 7.7 ★ Active-learning 50개 wet-lab 우선순위?
- proba ≈ 0.5 가장 모호한 unlabeled peptide 50개 → 검증하면 가장 정보량 큼
- 비용: tetramer 합성 + ELISpot ~$5k (50 peptides)
- → **추천**: 협업 lab (Korea-Sweden?) + Azure burst 비용처럼 $5k 1회

### 7.8 cancer_vaccine_agent.html agent에 T31 TCR-structure tool 등록?
- 현재 SPARK agent 21 scenarios + 31 tools (T1-T31)
- T31 = TCR-structure prediction (ESMFold CDR3β + similarity search) 추가하면 demo 완성
- → **추천**: 다음 라운드 + 짧은 작업 (~30분)

### 7.9 npj/Cell Rep Med 제출 timing?
- 현재: 충분한 수치 + 정직한 caveat + 백신 deliverable + agent demo + 9 라이브 페이지
- 추가하면 좋을 것:
  - ESM2 LOSO fix (0.40 → 0.55+)
  - Per-allele HLA-stratified 모델 재학습 (Simpson's paradox fix)
  - 50개 active-learning wet-lab pilot data (1-2개라도 positive)
- → **결정 필요**: 지금 제출 vs ESM2 후 제출
- **추천**: ESM2 + per-allele 재학습 기다리기 (2-3일 추가, AUROC 0.5+ 개선) → cell rep med fit 강해짐

---

## 8 · 다음 라운드 Top 5 (Seungho 명령 기다림)

1. **ESM2 35M embeddings** (CPU torch 설치 → LOSO 0.40 → ?) — 5분 코드, 효과 클 가능성
2. **HLA-stratified RF/XGBoost** (Simpson's paradox fix) — 10분 코드
3. **All 171 chunks 다 돌리는 KK 후속**: 290 exact overlaps의 cancer-relevance 분류 (TAA vs neoantigen vs viral mimicry) — 30분
4. **AlphaFold-Multimer GPU pilot** (Azure burst 1-2 GPU-hr × 5 complexes) — $30 비용, 진짜 structure 신호
5. **wet-lab 협업 outreach drafts** (Yu 교수, 서울대 Bundang, Korea-Sweden) — Gmail draft 자동 생성

---

**Lumenix · Seungho Cook · 2026-05-07**
Generated `2026-05-07` · Total session compute: ~30 minutes wall, 13 datasets ingested, 43 analyses run, 9 web pages live, 115 ESMFold structures stored.

**Status**: ★ Ready for paper submission OR for next-round ESM2/wet-lab. **Decision pending: §7.1, §7.2, §7.4, §7.7, §7.9.**
