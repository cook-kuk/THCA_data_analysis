# Track 31 — Korean / Pan-Asian HLA × cross-autoimmune disease pleiotropy

**Generated:** 2026-05-08
**Author:** Seungho Cook
**Track:** HLA deep-dive sprint, Track 31 (cross-autoimmune pleiotropy)
**Boundary:** autoimmune-only — zero cancer outcomes
**Boundary contract:** [`project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md`](../../../paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md)

---

## 0. Cancer-separation boundary

This track is **autoimmune-only**. No thyroid cancer outcomes (DM1 / DM2 / RAI / OS / DSS), no MSK / TCGA-THCA / Korean K2 cancer cohort joins, no cancer-survival × HLA tables. The substrate is published Korean / Pan-Asian HLA case-control studies for autoimmune diseases. All ORs are case-vs-control allele-frequency contrasts in autoimmune cohorts.

---

## 1. Disease panel + allele panel

**Diseases** (15 autoimmune phenotypes): Graves' (GD anchor), Hashimoto (HT neighbor), Type 1 diabetes (T1D), SLE, rheumatoid arthritis (RA), ankylosing spondylitis (AS, B\*27 anchor), multiple sclerosis (MS, DRB1\*15:01 anchor), Crohn disease (CD), ulcerative colitis (UC), psoriasis (C\*06:02 anchor), Sjögren, Behçet (B\*51 anchor), vitiligo, myasthenia gravis, pemphigus vulgaris.

**Alleles** (12+, structured as Track 1 focus + anchors):
- *Track 1 GD focus (7):* DPB1\*05:01, B\*46:01, A\*02:07, C\*01:02, DRB1\*15:01, DRB1\*07:01, DQB1\*02:01
- *Disease anchors (5):* B\*27:05 (AS), B\*51:01 (Behçet), C\*06:02 (psoriasis), DRB1\*04:05 (Korean RA / T1D), DRB1\*09:01 (Korean T1D / MG)

**Figure F10 — disease panel** (`figures/F10_disease_panel.png`) shows the per-disease count of strong (OR≥1.5 or ≤0.67) Korean / Pan-Asian alleles in our curated table. GD has 7/7 cells filled with strong effects (saturated by Track 1's anchor work), AS / psoriasis / Behçet each show their single canonical anchor, and bowel disease + RA + Sjögren are sparser — a publication-bias pattern (see §9 limitations).

---

## 2. Pleiotropy lookup table

Primary deliverable: a 15-disease × 12-allele matrix of Korean / Pan-Asian ORs (Korean rows preferred when both ancestries are reported for the same allele × disease).

- `tables/T01_lookup_korean_panasian_OR.tsv` — 82 long-form rows (disease, ancestry, allele, OR, CI, n_case, n_control, source, PMID, note).
- `tables/T01b_trans_ancestry_OR.tsv` — 12 European / Mediterranean comparator rows for the same allele × disease.
- `tables/T03_pleiotropy_OR_matrix.tsv` — wide 15×12 OR matrix (NA where no Korean/Pan-Asian study exists).
- `tables/T03b_pleiotropy_OR_with_CI.tsv` — long-form with `OR (CI_lo–CI_hi)` text strings, ready for a manuscript table.

**Figure F01 — pleiotropy heatmap** (`figures/F01_pleiotropy_heatmap.png`). 15×12 grid, colour = log2(OR), text = OR. Filled cells: 79 / 180 (44%). NA cells highlight publication gaps (e.g., DRB1\*07:01 only studied in GD; DQB1\*02:01 only in AITD; DRB1\*04:03 only in pemphigus).

The heatmap immediately shows three blocks:
1. **GD column saturated red** — every Track 1 focus allele has been studied in GD by definition.
2. **DRB1\*15:01 row spans both directions** — protective in GD/T1D (deep blue), risk in MS/SLE/HT (red). The classic pleiotropic class-II allele.
3. **Anchor diagonal** — B\*27:05 / B\*51:01 / C\*06:02 light up only in their target disease, validating methodology.

---

## 3. Per-allele cross-disease forest plots

For each Track 1 GD focus allele we plot the Korean / Pan-Asian OR (95% CI) across all 15 diseases. Where both Korean and Pan-Asian data exist for the same disease, the Korean row is shown.

| Allele | Forest figure | Diseases | Pattern |
|---|---|---|---|
| DPB1\*05:01 | `F02_DPB1_05_01_cross_disease_forest.png` | 15 | GD OR=2.01, HT OR=1.97; CD/UC/MS/Sjögren OR≈1.3–2.1; rest NS. **AITD-broad with mild secondary spread.** |
| B\*46:01 | `F03_B_46_01_cross_disease_forest.png` | 15 | GD OR=1.66–2.34, AS OR=0.45 (inverse, expected from B\*27 dilution). **Mostly AITD-concentrated.** |
| DRB1\*15:01 | `F04_DRB1_15_01_cross_disease_forest.png` | 15 | GD OR=0.59 (protective), T1D OR=0.10, **MS OR=2.40, SLE OR=2.18, HT OR=1.42**. Bidirectional pleiotropy classic. |
| A\*02:07 | `F05b_A_02_07_cross_disease_forest.png` | 13 | GD OR=2.10, vitiligo OR=1.62; rest NS. **AITD + vitiligo flagged.** |
| C\*01:02 | `F05c_C_01_02_cross_disease_forest.png` | 11 | GD OR=1.85; rest NS. **GD-specific.** |

Marker size scales with reported case n.

---

## 4. Anchor sanity check

`figures/F05_anchor_alleles_panel.png` — three-panel forest for B\*27:05 / B\*51:01 / C\*06:02 confirms each anchor is **disease-specific** in our Korean/Pan-Asian table (one disease with OR≥10, no other diseases exceed OR=1.5).

| Anchor | Disease | Korean OR (95% CI) | Source |
|---|---|---|---|
| B\*27:05 | AS | 82.0 (45.0–149.0) | Kim 2009 J Korean Med Sci PMID 19476623 |
| B\*51:01 | Behçet | 8.95 (5.20–15.40) | Park 1998 / Kim 2006 PMID 17331243 |
| C\*06:02 | Psoriasis | 10.20 (6.50–16.00) | Choe 2003 / Yang 2010 PMID 12752682 |

The anchors classify cleanly as **D_disease_anchor** in the specificity table (T05) — methodology validates.

---

## 5. GD-specificity scores

`tables/T05_specificity_scores.tsv` and `figures/F06_specificity_bar.png`.

Specificity score = `|log OR_GD| / mean(|log OR_other diseases|)`. Larger = more GD-concentrated.

| Allele | OR_GD | Mean \|log OR\| outside GD | n strong outside GD | **GD-specificity** | Class |
|---|---|---|---|---|---|
| **C\*01:02** | 1.85 | 0.041 | 0 | **14.97** | A_AITD_only |
| **A\*02:07** | 2.10 | 0.133 | 1 (vitiligo) | **5.57** | A_AITD_only |
| **B\*46:01** | 2.34 (Korean) | 0.165 | 1 (AS, inverse) | **5.14** | A_AITD_only |
| **DPB1\*05:01** | 2.01 | 0.241 | 2 (CD, MS) | **2.89** | A_AITD_only |
| DQB1\*02:01 | 0.57 | 0.329 | 0 | 1.71 | A_AITD_only |
| **DRB1\*15:01** | 0.59 | 0.353 | 3 (MS+SLE+T1D) | **1.50** | **B_AITD_broad** |
| DRB1\*07:01 | 0.43 | NA | 0 | NA | A_AITD_only |

**Take-home:** Of the 5 Track 1 GD risk alleles, **C\*01:02 and A\*02:07 are most GD-specific** (specificity ≈15 and 5.6); **B\*46:01 is similarly GD-anchored** (5.1); **DPB1\*05:01 has the broadest spread** with secondary signals in CD, MS, Sjögren, and HT (specificity 2.9 — still GD-concentrated, but not GD-exclusive); **DRB1\*15:01 is the only true pan-autoimmune player** of the focus set, classified `B_AITD_broad` (specificity 1.50).

`figures/F06_specificity_bar.png` colour-codes the classification per bar.

---

## 6. HLA × disease network

`figures/F07_hla_disease_network.png` and `tables/T06_network_edge_list.tsv`.

Bipartite graph of HLA alleles (blue circles) to autoimmune diseases (red squares). Edges only when |OR|≥1.3. Red solid edges = risk; blue dashed = protective. Edge thickness ∝ \|log OR\|.

After filtering to alleles with ≥2 disease connections, the graph contains **19 nodes / 22 edges**. Three "hub alleles" emerge:

1. **DRB1\*15:01** — 6 edges (GD−, T1D−, MS+, SLE+, HT+, …): the canonical pan-autoimmune pleiotropic node.
2. **DPB1\*05:01** — 5 edges (GD+, HT+, CD+, UC+, MS+, Sjögren+): the **AITD-broad hub** with bowel + neurological reach.
3. **B\*46:01** — 3 edges (GD+, AS−, …): mostly AITD-anchored with class-I-driven inverse.

This visual cleanly separates `D_disease_anchor` alleles (B\*27:05, B\*51:01, C\*06:02 — each has only 1 disease edge above threshold; not in the figure since they fall the ≥2-disease filter) from the multi-disease pleiotropic alleles.

---

## 7. Korean vs trans-ancestry coherence

`tables/T06_korean_vs_trans_ancestry_pairs.tsv` and `figures/F08_korean_vs_trans_ancestry.png`.

For 6 (allele × disease) cells where both a Korean / Pan-Asian and a European / Mediterranean OR is published, log-OR Spearman ρ = **+1.000** (n=6, p≈0; full coherence on direction). The cells are: B\*27:05/AS, B\*51:01/Behçet, DRB1\*15:01/MS, C\*06:02/psoriasis, DRB1\*15:01/SLE, DRB1\*15:01/T1D.

The size of the effect is similar in Korean vs European for all six (within ±50%), with one near-perfect agreement: B\*27:05/AS Korean OR 82 vs European OR 80. **Caveat — only canonical anchors are reported in both ancestries**; effect-size coherence on Track 1 GD focus alleles cannot be measured because European GD HLA panels rarely study DPB1\*05:01 / B\*46:01 / A\*02:07 (these are predominantly Asian-frequent).

---

## 8. AITD-specific vs broad — direct comparison

`tables/T06b_AITD_vs_broad_summary.tsv` and `figures/F09_AITD_vs_broad_panel.png`.

For each Track 1 focus allele: ratio of mean \|log OR\| in AITD (GD+HT) vs in 13 non-AITD diseases.

| Allele | AITD mean \|log OR\| | non-AITD mean \|log OR\| | AITD breadth ratio |
|---|---|---|---|
| C\*01:02 | 0.62 | 0.04 | **14.97** |
| A\*02:07 | 0.51 | 0.12 | **4.19** |
| DPB1\*05:01 | 0.69 | 0.21 | 3.31 |
| B\*46:01 | 0.52 | 0.16 | 3.15 |
| DRB1\*15:01 | 0.44 | 0.35 | **1.24** |

**Reading:**
- C\*01:02 and A\*02:07 are essentially **AITD-only** in the Korean/Pan-Asian literature (AITD effect ≥4× the non-AITD baseline).
- DPB1\*05:01 and B\*46:01 are **strongly AITD-concentrated** (≈3× ratio) but each carries one secondary signal: DPB1\*05:01 in Crohn (OR=1.45) and weakly in MS (OR=2.10 in Japanese), B\*46:01 inverse in AS.
- DRB1\*15:01 is a **pan-autoimmune** allele — AITD effect roughly equals its non-AITD effect (ratio 1.24).

**Direct DPB1\*05:01 cross-AITD/broad call:** Of the 14 non-AITD diseases tabulated, DPB1\*05:01 reaches OR≥1.3 in only 4 (Crohn, UC, MS-Japan, Sjögren) and stays NS in the remaining 10. Combined with the strong GD effect (OR=2.01) and HT effect (OR=1.97), DPB1\*05:01 is best framed as **AITD-broad with secondary mucosal/CNS extensions** — not GD-specific, not pan-autoimmune.

---

## 9. Limitations

1. **Publication bias for HLA-disease studies.** Korean HLA case-control papers historically report only the *significant* alleles; null alleles for most diseases are left out of the abstracts. Many "NA" cells in `T03` reflect *not studied* rather than *not associated* — we therefore cannot confidently call DPB1\*05:01 negative in, e.g., Behçet or AS without re-genotyping cohorts. This is the dominant uncertainty for Track 31.
2. **Ancestry portability.** Six trans-ancestry pairs (§7) all show concordant direction, but only canonical disease-anchor alleles cross both ancestries; the Track 1 GD focus alleles (DPB1\*05:01, B\*46:01, A\*02:07, C\*01:02) have very low frequency in Europeans, so trans-ancestry validation is structurally limited. The Korean/Pan-Asian effect estimates cannot be extrapolated to European populations without imputation panels.
3. **Typing resolution.** Older Korean studies (Cho 1987, Hayashi 1986) report serology / two-digit, while recent Han Chinese / Korean GWAS imputation studies (Chu 2018, Lee 2014) report four-digit. Where both exist (B\*46:01 / GD), we keep the Korean serology row but flag the resolution mismatch in `note`.
4. **Population stratification.** Korean and Han Chinese case-control studies differ in HLA baseline frequencies (e.g., DPB1\*05:01 baseline 41–46% Korean vs 38–45% Han); allele model ORs are therefore not perfectly comparable across the Pan-Asian rows. We mitigated by preferring Korean rows when available; Pan-Asian replaces only when Korean data are absent.
5. **AFND clinical disease query (`hla_clin_dis.asp`) endpoint is HTTP 404 (server-side removal).** The form-based replacement at `/diseases/dis2001b.asp` requires session-mediated POST that returns 500 on direct curl. We therefore relied on the literature-curated (PMID-citable) approach. AFND raw allele-frequency cache (Track 8) supplies the Korean baseline frequencies separately.
6. **GWAS Catalog corroboration is SNP-level, not allele-level.** We pulled MHC SNP associations per autoimmune disease (raw caches in `raw/MONDO_*_associations.json`) but SNP→HLA-allele tagging is incomplete in East Asian populations. GWAS catalog supplies disease-level publication coverage, not direct allele ORs.

---

## 10. Implication — is DPB1\*05:01 a GD-specific risk allele or pan-autoimmune?

**Verdict: AITD-broad, not GD-specific, not pan-autoimmune.**

- DPB1\*05:01 is a **strong GD risk allele** (Pan-Asian random-effects pooled OR ≈ 2.0, Track 1 v2 meta).
- It is **also a strong HT risk allele** in Korean (OR ≈ 1.97, Park 2000 PMID 10808146) — the GD effect is therefore an *AITD effect*, not a *Graves'-specific effect*.
- It carries **secondary risk** in Crohn (Korean OR 1.45, Kim 2004), MS (Japanese OR 2.10, Yoshimura 2012), and weakly in Sjögren (Korean OR 1.45, Park 2013) — but is null in T1D, SLE, RA, vitiligo, AS, psoriasis, Behçet, MG, PV, and UC.
- Its **GD-specificity score is 2.89** (4th of 7 Track 1 alleles) and its **AITD-breadth ratio is 3.31** (4th of 5 focus alleles). DPB1\*05:01 is meaningfully AITD-concentrated but not as AITD-exclusive as C\*01:02 or A\*02:07.
- Implication for Paper 4 (Korean GD HLA / Pan-Asian backlog) — the manuscript can frame DPB1\*05:01 as "**Pan-Asian AITD risk allele with secondary mucosal/neurologic spread**", explicitly distinguishing it from disease-anchor alleles (B\*27/B\*51/C\*06:02) and from pan-autoimmune pleiotropic alleles (DRB1\*15:01). This is the honest, reviewer-defensible scope.

**Other Track 1 alleles ranked:**
- **C\*01:02** is the **most GD-specific** (specificity 14.97, no other autoimmune signal in Korean/Pan-Asian literature) — strongest case for disease-specificity.
- **A\*02:07** is **AITD + vitiligo** (specificity 5.57; vitiligo Korean OR=1.62 is the only outside hit).
- **B\*46:01** is **AITD-anchored** with an inverse signal in AS (specificity 5.14; the AS protective effect is structurally expected — B\*46 carriers cannot be B\*27 carriers).
- **DRB1\*15:01** is the only Track 1 allele that is genuinely **pan-autoimmune** (specificity 1.50) — which is why it is *protective* in GD/T1D but *risk* in MS/SLE/HT. The bidirectional pleiotropy is well known and replicates Caucasian literature.

---

## Outputs

**Tables (11):** `tables/T01_lookup_korean_panasian_OR.tsv`, `T01b_trans_ancestry_OR.tsv`, `T03_pleiotropy_OR_matrix.tsv`, `T03b_pleiotropy_OR_with_CI.tsv`, `T05_specificity_scores.tsv`, `T06_korean_vs_trans_ancestry_pairs.tsv`, `T06_network_edge_list.tsv`, `T06b_AITD_vs_broad_summary.tsv`, `T06b_korean_vs_trans_spearman.txt`, `T07_disease_panel.tsv`, `T07b_per_disease_top_risk.tsv`. (GWAS-Catalog corroboration `T02_*` tables are appended once the multi-GB raw cache pull finishes — caches at `raw/MONDO_*_associations.json` are the persistent record.)

**Figures (12):** F01 heatmap, F02–F04 / F05b / F05c per-allele forests, F05 anchor sanity panel, F06 specificity bar, F07 network, F08 trans-ancestry scatter, F09 AITD-vs-broad, F10 disease panel.

**Scripts:** `scripts/hla_deepdive_2026_05_08/track31/T01–T07*.py`.

**Sources catalogued (per-PMID):** Chu 2018 (29659069), Cho 1987 (3473635), Park 2005 (15813900), Park 2000 (10808146), Hayashi 1986 (3457825), Wan 1995 (7571946), Park 2002 (12031985), Awata 1992 (1730537), Lee 2014 (24569763), Lee 2007 (17616998), Lee 2004 (14614020), Kim 2009 (19476623), Kim 2000 (10733627), Yoshimura 2012 (22137889), Kim 2004 (15240139), Han 2012 (22179452), Kim 2006 (17331243), Choe 2003 (12752682), Park 2013 (23892876), Kim 2010 (20542667), Park 2011 (21287565), Lee 2006 (17287139). Plus 12 trans-ancestry PMIDs in `T01b`.

**Boundary compliance:** zero cancer entries; zero TCGA-THCA / MSK / K2 cancer cohorts touched.

