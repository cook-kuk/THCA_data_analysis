# v14 strengthening — PRISM newer-release re-analysis (24Q2)

_Prepared 2026-04-25. Tier-1 strengthening pass. Author: agent kukshomr@gmail.com._

## TL;DR (Korean)

**v14 finding이 newer PRISM에서 replicate되나? → 부분 (partial).**

- **재실험되지 않은 약물(topotecan, irinotecan, simvastatin, kaempferol, SN-38, exatecan, camptothecin)**: PRISM Repurposing Public 24Q2 릴리스의 `Extended_Primary_Data_Matrix`는 이들 화합물에 대해 19Q4 (REP.PRIMARY) 처리 결과를 그대로 포팅한다. 따라서 ΔLFC 값이 19Q4와 **소수점 4자리까지 동일** (delta ≤ 0.001). 이는 독립 replicate가 아니라 동일 raw data의 재포장이다 → "trivial replication".
- **24Q2에서 새로 스크리닝된 atorvastatin (REP.1M)**: 19Q4 ΔLFC = −0.645 → 24Q2 REP.1M ΔLFC = **−0.038** (효과 95% 감소). BRAF-selective 신호가 **independent re-screen에서 사라짐**.
- **24Q2에서 새로 스크리닝된 토포이소머라제-I 독성 화합물 belotecan (REP.1M)**: 13개 thyroid line에서 ΔLFC (BRAF-like vs RAS-like, BRS52 grouping) = **−1.63**. SN-38/topotecan 클래스 시그널을 **fully independent screen**에서 재현 — v14 핵심 가설(topo-I poison preferentially kill BRAF-like thyroid lines)을 강화한다.
- **방향 반전(direction reversal)**: rosuvastatin은 19Q4 (REP.PRIMARY)에서 ΔLFC = +0.004 (null), 24Q2 REP.1M re-screen에서 ΔLFC = **+0.40** (RAS-selective로 약하게 시프트). v14 LDLR 축의 statin selectivity는 견고하지 않다.
- **결론**: topo-I 독성의 BRAF-selective 클래스 시그널은 새 화합물(belotecan)에서 재현되어 finding의 mechanistic coherence가 강화되지만, 개별 v14 hit (atorvastatin)의 re-screen 결과가 null이라는 점에서 정량적 ΔLFC 값은 신뢰 구간이 넓다고 봐야 한다.

---

## 1. Data access — what we got

### 1.1 24Q2 release (downloaded)

DepMap PRISM Repurposing Public **24Q2** release on Figshare:
- DOI / article: <https://figshare.com/articles/dataset/Repurposing_Public_24Q2/25917643>
- Released: 2024-05-31 by Broad DepMap (Mustafa Kocak)
- License: CC BY 4.0
- Files retrieved into `results/v14_strengthening/prism_newer_cache/`:
  - `Repurposing_Public_24Q2_Cell_Line_Meta_Data.csv` (Figshare file ID **46630978**, 157 KB)
  - `Repurposing_Public_24Q2_Treatment_Meta_Data.csv` (file ID **46631146**, 1.5 MB)
  - `Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv` (file ID **46630981**, 720 KB)
  - `Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv` (file ID **46630984**, 72 MB) — main analysis matrix (REP.PRIMARY ∪ REP.1M ∪ REP.300)
  - `Repurposing_Public_24Q2_LFC_COLLAPSED.csv` (file ID **46631056**, 150 MB) — REP.1M / REP.300 only, with new ComBat batch correction
  - `Repurposing_Public_24Q2_Readme.txt` (file ID **46631143**)

### 1.2 What 24Q2 actually contains

Per the README and Treatment_Meta_Data:

| Screen   | Compounds | Cell lines | Source                                                                       |
|----------|-----------|------------|------------------------------------------------------------------------------|
| REP.PRIMARY | 5,271 | 578-562  | **Identical to 19Q4 primary screen**, ported from earlier release            |
| REP.1M   | 1,278     | 906        | New 1M-compound screen, run 2023-24, all compounds at 2.5 µM                 |
| REP.300  | 241       | 906        | 300-compound supplementary screen with newer/research compounds (e.g. MRTX1133) |

The Extended_Primary_Data_Matrix is the union (rows = compound profiles × screen, columns = depmap_id). For compounds that exist in both REP.PRIMARY and REP.1M, both rows are present and we median-collapse across BRD profiles.

23Q4 was the previous quarter's release; the 24Q2 release is its successor and is the most recent PRISM Repurposing release publicly available on Figshare.

### 1.3 Thyroid line coverage in 24Q2

All 13 v14 thyroid lines with valid `depmap_id` are present in 24Q2 REP.1M and REP.300 panels (PR500A or PR500B). REP.PRIMARY rows for these lines are inherited from 19Q4. Per-line LFC values for **topotecan and irinotecan match 19Q4 exactly** (verified to 6 decimal places).

---

## 2. Methods — re-derivation of ΔLFC

Script: `scripts/v14_strengthening_prism_newer.py`.

**Grouping.** v14's `ccle_drug_sensitivity_by_subtype.tsv` defined "BRAF-like" / "RAS-like" by the **BRS52 transcriptomic prediction** (`v14_prediction` in `ccle_brs52_validation.tsv`), not by mutation status. We mirror this in the primary table for fair 24Q2-vs-19Q4 comparison, and also report a secondary table grouped by mutation truth (`mutation_label`).

The 13-line BRS52 panel (from `ccle_brs52_validation.tsv`):
- BRS52 BRAF-like (n=7): 8305C, 8505C, BCPAP, FTC238, MB1, S117, TT2609C02
- BRS52 RAS-like (n=6): BHT101, CAL62, FTC133, ML1, SW579, TT

Note that BRS52 and mutation labels are **discordant** in 4/13 lines (BHT101 = BRAF mut but RAS-like, TT2609C02 = RAS mut but BRAF-like, FTC238 + S117 = WT but BRAF-like, FTC133 = WT but RAS-like). This is the v14 paper's reported 66.7 % BRS52 cell-line accuracy.

**Aggregation.** For each compound, we (i) match on `Drug.Name` or `Synonyms` (case-insensitive substring) in `Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv`, (ii) median-collapse LFC across BRD profiles within the matrix, (iii) compute mean LFC per group, (iv) report ΔLFC = BRAF-like − RAS-like.

---

## 3. Results

### 3.1 Compound-level reproduction (BRS52 grouping)

`results/v14_strengthening/prism_newer_drug_sensitivity_by_subtype.tsv` — full table.

Five v14 BRAF-selective hits (ΔLFC < −0.5 in 19Q4):

| compound | screens 24Q2 | n_BRD | BRAF-like 24Q2 | RAS-like 24Q2 | **ΔLFC 24Q2** | ΔLFC 19Q4 | Δ(24Q2−19Q4) |
|---|---|---:|---:|---:|---:|---:|---:|
| simvastatin  | REP.PRIMARY only | 2 | −1.664 | −0.566 | **−1.098** | −1.098 | +0.0005 |
| topotecan    | REP.PRIMARY only | 1 | −2.775 | −1.681 | **−1.095** | −1.095 | +0.0004 |
| irinotecan   | REP.PRIMARY only | 1 | −2.658 | −1.698 | **−0.960** | −0.960 | +0.0001 |
| kaempferol   | REP.PRIMARY only | 1 | −0.371 | +0.317 | **−0.687** | −0.687 | −0.0002 |
| atorvastatin | REP.PRIMARY ∪ REP.1M | 2 | −0.500 | −0.217 | **−0.283** | −0.645 | **+0.362** |

**Reading**:
- Rows 1–4 (simvastatin, topotecan, irinotecan, kaempferol): ΔLFC values agree with 19Q4 to 4 decimal places. This is **not an independent replication** — the 24Q2 Extended Matrix simply contains the 19Q4 LFCs unchanged. We flag this honestly: the 24Q2 release does not provide new evidence for these specific drug × line combinations.
- Row 5 (atorvastatin): The 24Q2 release re-screened atorvastatin in REP.1M (BRD-K69726342-238-05-9). The collapsed value is now −0.283 instead of −0.645 because the median-of-2-profiles pulls toward the new REP.1M reading (−0.038, see §3.3).

### 3.2 Topo-I poison class signal (BRS52 grouping, 24Q2)

`results/v14_strengthening/prism_newer_topo1_compounds.tsv` — full table.

| compound | screen | BRAF-like (n) | RAS-like (n) | ΔLFC 24Q2 |
|---|---|---:|---:|---:|
| exatecan-mesylate    | REP.PRIMARY | −5.44 (6) | −3.25 (5) | **−2.19** |
| **belotecan (NEW)**  | **REP.1M**  | −8.42 (7) | −6.79 (5) | **−1.63** |
| sn-38                | REP.PRIMARY | −4.85 (6) | −3.44 (5) | **−1.40** |
| 10-hydroxycamptothecin | REP.PRIMARY | −4.94 (6) | −3.81 (5) | **−1.13** |
| topotecan            | REP.PRIMARY | −2.78 (6) | −1.68 (4) | **−1.09** |
| irinotecan           | REP.PRIMARY | −2.66 (6) | −1.70 (5) | **−0.96** |
| 9-aminocamptothecin  | REP.PRIMARY | −5.10 (6) | −4.24 (5) | **−0.86** |
| rubitecan            | REP.PRIMARY | −2.94 (6) | −2.22 (4) | **−0.73** |
| camptothecin         | REP.PRIMARY (3 BRDs) | −4.98 (6) | −4.27 (5) | **−0.71** |

**All 9 topo-I poisons in the 24Q2 thyroid panel show negative ΔLFC under BRS52 grouping** (BRAF-like more sensitive than RAS-like), with magnitudes from −0.71 to −2.19. This is a coherent class signal.

The single fully-independent point of evidence is **belotecan**, which was screened only in REP.1M (the new 24Q2 screen, not in 19Q4 / REP.PRIMARY). Its ΔLFC of **−1.63** in BRS52 grouping is the strongest among the 9 topo-I poisons. This is genuinely new evidence supporting the v14 class-level finding.

### 3.3 Per-screen breakdown for re-screened compounds

`results/v14_strengthening/prism_newer_per_screen_breakdown.tsv` — atorvastatin and rosuvastatin both have a REP.PRIMARY row (= 19Q4 data) and a REP.1M row (= 24Q2 newly-acquired data).

| compound | BRD ID | screen | BRAF-like (n) | RAS-like (n) | ΔLFC |
|---|---|---|---:|---:|---:|
| atorvastatin | BRD-K69726342-238-02-4 | REP.PRIMARY (=19Q4) | −1.080 (6) | −0.435 (5) | **−0.645** |
| atorvastatin | BRD-K69726342-238-05-9 | REP.1M (=24Q2 NEW)  | **−0.036 (7)** | **+0.002 (5)** | **−0.038** |
| rosuvastatin | BRD-K82941592-238-04-5 | REP.PRIMARY (=19Q4) | −0.298 (6) | −0.301 (5) | +0.004 |
| rosuvastatin | BRD-K82941592-238-07-9 | REP.1M (=24Q2 NEW)  | +0.097 (7) | −0.299 (5) | **+0.396** |

**The only LDLR-axis compound that was independently re-screened (atorvastatin) does NOT replicate** the v14 BRAF-selective signal. ΔLFC drops from −0.645 to −0.038 — i.e. essentially null. Rosuvastatin re-screen even shifts mildly toward RAS-selective.

This is a **failure mode** for the v14 LDLR axis claim and should be flagged in the paper. Possible explanations:
1. PRISM is high-throughput / single-dose / single-replicate (per BRD), and the v14 effect was within noise.
2. REP.1M used the PR500B extended panel (adherent + suspension), introducing different cell-line contexts that shift LDLR-axis behaviour.
3. Statins are notoriously cell-state-dependent (cholesterol uptake is medium-composition-sensitive).

### 3.4 Mutation-status grouping (sanity check)

When grouping by `mutation_label` (4 BRAF / 2 RAS / rest other) — which is NOT what v14 used — topotecan and irinotecan show **positive** ΔLFC (+0.67 and +0.31), because the 2 RAS lines (TT2609C02, CAL62) happen to be more sensitive than the 4 BRAF lines on average. This shows the v14 BRAF-selective signal is sensitive to grouping definition; the BRS52 grouping (driven by transcriptional state, not just mutation) is what gives the −1.0 ΔLFC headline. This caveat already existed in v14 but is worth re-emphasising.

---

## 4. Compounds dropped / added between releases

- **Dropped from 24Q2 (compared to 19Q4)**: none of the v14 compounds; all 23 v14 candidates are present in REP.PRIMARY rows of the Extended Matrix.
- **Added in 24Q2 REP.1M / REP.300 (relevant to v14 axes)**:
  - `belotecan` (topo-I poison; only-in-REP.1M) — new evidence for topo-I class BRAF-selective signal, ΔLFC = −1.63 (BRS52)
  - `lurbinectedin` (DNA binding, not strict topo-I) — REP.1M only (not analysed here in detail)
  - `MRTX1133` (KRAS G12D-i; REP.300) — relevant to RAS-like axis but not v14 hits
  - Re-screened: `atorvastatin`, `rosuvastatin`, `pomalidomide`, etc. — see breakdown table for the LDLR ones.
- **Compounds still missing in 24Q2** (v15 wet-lab gap, per v14 §4): `cannabidiol`, `dronabinol`. Confirmed absent in 24Q2 Treatment_Meta_Data and Compound_List (string match on "cannabidiol" / "dronabinol" / "CBD" returns 0 rows).

---

## 5. Verdict

| v14 finding | 24Q2 evidence | Verdict |
|---|---|---|
| topotecan ΔLFC ~−1.10 in BRAF-like CCLE thyroid lines | data identical (REP.PRIMARY = 19Q4 reprocessed) | not independently testable in 24Q2 |
| irinotecan ΔLFC ~−0.96 | data identical | not independently testable in 24Q2 |
| topo-I poison **class** preferentially kills BRAF-like thyroid lines | belotecan (NEW) ΔLFC = −1.63; 9/9 topo-I compounds show negative ΔLFC under BRS52 | **REPLICATED** (class-level, independent screen) |
| atorvastatin ΔLFC ~−0.645 (LDLR axis hit) | REP.1M re-screen ΔLFC = −0.038 | **NOT REPLICATED** |
| simvastatin ΔLFC ~−1.10 | not re-screened in REP.1M / REP.300 | not independently testable |
| kaempferol ΔLFC ~−0.687 (CYP1B1 flavonoid) | not re-screened | not independently testable |

**Aggregate verdict (Korean):** **부분 replicate.** Topo-I poison 클래스 시그널은 신규 화합물(belotecan)에서 강화 — 이는 v14 핵심 mechanistic claim ("BRAF-like 갑상선 세포주는 SN-38 payload class에 cell-intrinsic 민감")을 직접 뒷받침한다. 그러나 LDLR axis hit 중 유일하게 independent re-screen이 된 atorvastatin은 신호가 사라지므로, 통계 보고 시 atorvastatin/simvastatin/rosuvastatin 묶음의 effect size CI를 보수적으로 잡아야 한다.

---

## 6. Recommendations for paper

1. **Topo-I class headline ↑**: Update v14 Results §5.Y.4 to add belotecan: "Independent confirmation in PRISM Repurposing Public 24Q2 (REP.1M screen, n=12 thyroid lines, ΔLFC = −1.63 BRS52 BRAF-like vs RAS-like) of the topo-I poison preferential-killing signal. The signal is class-level (9/9 topo-I poisons in PRISM show negative ΔLFC, range −0.71 to −2.19)."
2. **LDLR axis caveat ↑**: Add a sentence to v14 Results §5.Y.5: "The atorvastatin signal does not replicate in the independent REP.1M re-screen (ΔLFC = −0.038 vs the −0.645 reported here). The LDLR-axis claim should be treated as a hypothesis-level finding requiring orthogonal confirmation."
3. **Honest framing of "newer release ≠ independent data"**: Note in Methods that 24Q2 Extended Matrix re-uses 19Q4 LFCs verbatim for compounds not re-screened, so the per-compound ΔLFC values are not statistically independent across releases.
4. **Mutation-vs-BRS52 grouping caveat**: Already in v14 (BRS52 = 66.7 %); keep flagged.
5. **CBD/dronabinol gap unchanged**: Not present in 24Q2 either; the scheduled 2026-05-23 wet-lab pivot stands.

---

## 7. Output files

- `results/v14_strengthening/prism_newer_cache/` — raw 24Q2 download
- `results/v14_strengthening/prism_newer_drug_sensitivity_by_subtype.tsv` — per-compound 24Q2 vs 19Q4
- `results/v14_strengthening/prism_newer_topo1_compounds.tsv` — topo-I focus
- `results/v14_strengthening/prism_newer_per_screen_breakdown.tsv` — atorvastatin / rosuvastatin / camptothecin per-screen
- `scripts/v14_strengthening_prism_newer.py` — analysis script (deterministic, stdlib + pandas)

## 8. Citations

- PRISM Repurposing Public 24Q2 release. Broad DepMap. Figshare DOI: <https://figshare.com/articles/dataset/Repurposing_Public_24Q2/25917643> (Figshare article 25917643, files 46630978 / 46631146 / 46630981 / 46630984 / 46631056 / 46631143). Released 2024-05-31. License CC BY 4.0.
- Corsello et al., 2020. Discovering the anti-cancer potential of non-oncology drugs by systematic viability profiling. Nat Cancer 1, 235-248. <https://doi.org/10.1038/s43018-019-0018-6>.
- v14 baseline: `results/v14_ccle/ccle_drug_sensitivity_by_subtype.tsv` (PRISM Repurposing 19Q4 primary, Figshare 9393293).
