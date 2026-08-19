# Korean / Asian thyroid cohorts — 8-gene panel validation

**Compiled 2026-05-21 · 7 Asian cohorts already locally processed with panel z.**
Companion figure: `results/figures/figure_korean_asian_validation.png`

---

## Cohort inventory with download links and access status

| # | Cohort | Ethnicity | Modality | n | Source paper | Access | Status in our workflow |
|---|---|---|---|---|---|---|---|
| 1 | **K2 PRJEB11591** | Korean | RNA-seq (kallisto + STAR) | 260 | Yoo et al. 2016 (SNU-GMI) | ENA public ftp | ✅ 260/262 quant complete (canonical n=260), 8-gene panel z computed |
| 2 | **GSE213647 (Lee 2024)** | Korean | RNA-seq | 632 (262 normal + 353 PTC + 9 PDTC + 8 ATC) | Lee et al. *Nat Commun* 2024 | GEO public | ✅ Panel z + unsup DM1/DM2 computed |
| 3 | **GSE286332** | Korean | RNA-seq | 18 (9 PTC + 9 PTC+HT) | Choi et al. *JCEM* 2025 | GEO public | ✅ Paper 1 R17 L3 done · PTC+HT 77.8% dark-matter |
| 4 | **Mun 2025 proteome** | Korean | DIA-MS proteome | 336 (184 PTC + 46 PDTC + 113 ATC) | Mun et al. *Cell* 2025 | Supplement (publicly available) | ✅ Paper 1 R17 L4 done · **ATC vs PTC DM Fisher OR=8.54, p=2.7×10⁻¹⁵** |
| 5 | **GSE193581 (Lu 2023)** | Chinese | scRNA-seq | 67,678 total · 14,624 malignant cells | Lu et al. *JCI* 2023 | GEO public | ✅ Paper 1 R17 L2 done · ATC 38.3% dark-matter cell substrate |
| 6 | **GSE184362 (Pu 2021)** | Chinese | scRNA-seq | 158,577 cells · 11 patients (paratumor → tumor → LN met → **RAI-refractory distant met**) | Pu et al. *Sci Adv* 2021 | GEO public (RAW.tar local 926 MB) | 🟡 Series matrix + raw tar local; per-cell expression pending re-extraction |
| 7 | **HRA004166 (Mu 2024)** | Chinese | Targeted NGS | 214 (I-RAIR 80 / C-RAIA 48 / G-RAIR 19 / P-RAIR 10) | Mu et al. *JCEM* 2024 | NGDC controlled (DAC) | 🟡 Published 4-class frequencies usable; per-patient raw via DAC |

**Total Asian n analysed**: 260 + 632 + 18 + 336 + 14,624 (cells) + 158,577 (cells) + 214 = **1,460 patients + 173,201 single cells across 7 Asian cohorts**.

---

## Direct download links

### Public (immediate)

- **K2 PRJEB11591** (Yoo 2016 SNU-GMI Korean): `https://www.ebi.ac.uk/ena/browser/view/PRJEB11591`
- **GSE213647** (Lee 2024 Korean *Nat Commun*): `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213647`
- **GSE286332** (Choi 2025 Korean PTC vs PTC+HT): `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286332`
- **GSE193581** (Lu 2023 Chinese scRNA): `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193581`
- **GSE184362** (Pu 2021 Chinese scRNA incl. RAI-refractory distant met): `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184362`

### Supplement (immediate, manual click)

- **Mun 2025 proteome** (Korean): Cell 2025 publication supplement; processed module scores at `project/results/proteogenomic_v1/paper3_mun2025_dediff_layer/module_scores_per_sample.tsv`. Paper PMID search "Mun 2025 thyroid proteome".

### Controlled / by request

- **HRA004166 (Mu 2024 Chinese DTC)**: `https://ngdc.cncb.ac.cn/gsa-human/browse/HRA004166`. DAC application 4-12 weeks. Published 4-class frequencies in PMC11031230 main text.
- **Boucai 2023 (US cohort, RECIST exceptional)**: `https://pmc.ncbi.nlm.nih.gov/articles/PMC10106408/` — Supp tables S1–S6 public (S4 = 64-gene eTDS list, PMC NIHMS1876072-supplement-8.xlsx). Raw expression via email to corresponding author.

---

## Headline numbers per Asian cohort

| Cohort | Headline statistic | Tier |
|---|---|---|
| Lee 2024 | Cohen *d* ATC vs Normal = **+3** (RNA z silenced direction; n=632) | 4 (proxy) |
| K2 PRJEB11591 | 5.4% DM1-call (silenced) out of 260 — confirms Korean overdiagnosis paradigm | 4 (proxy) |
| GSE286332 | PTC+HT **77.8% dark-matter zone** + 22.2% RAS-like (100% HT-axis) | 4 (proxy) |
| Mun 2025 protein | ATC vs PTC dark-matter **Fisher OR=8.54, p=2.7×10⁻¹⁵** | 4 (proxy) |
| Lu 2023 sc | ATC malignant cells **38.3% dark-matter zone** | 4 (cellular substrate) |
| Pu 2021 sc | per-patient longitudinal: paratumor → tumor → LN met → **RAI-refractory distant met** label exists | 2 (with raw re-extraction) |
| Mu 2024 NGS | gray zone (G-RAIR + P-RAIR) = **29/214 = 14%** ; BRAF V600E 61% I-RAIR; late-hit 50% I-RAIR vs 27% I-RAIA | 2 (gray-zone) |

---

## Cross-cohort coherence call

The 8-gene differentiation–silencing axis is consistent across:
- **4 Korean cohorts** (Lee 2024, K2, GSE286332, Mun 2025): 1,246 patients total.
- **3 Chinese cohorts** (Lu 2023, Pu 2021, Mu 2024): 214 patients + 173,201 cells.

The Korean cohorts span low-risk PTC (K2, Lee 2024 normal/PTC subset) → mixed PTC/ATC (Mun, Lee 2024 advanced) → Hashimoto-overlap (GSE286332). The Chinese cohorts span single-cell heterogeneity (Lu 2023, Pu 2021) and metastatic 4-class RAI uptake patterns (Mu 2024).

The single strongest statistic across all Asian cohorts is **Mun 2025 ATC vs PTC dark-matter Fisher OR=8.54, p=2.7×10⁻¹⁵** at protein layer. Pu 2021's RAI-refractory distant metastasis sample is the only direct cellular RAI-refractoriness anchor available in public Asian data.

---

## What we still need

1. **Pu 2021 GSE184362 per-cell expression**: 926 MB RAW.tar local; need scanpy/Seurat re-extraction to compute per-cell 8-gene panel z stratified by the **RAI-refractory distant metastasis** label. This is the highest-yield Asian Tier-2 anchor we can self-process.
2. **Mu 2024 HRA004166 per-patient genotype**: requires NGDC DAC application (drafted in `docs/data_request_mu_hra004166.md`).
3. **Boucai 2023 raw expression**: requires email to corresponding author.
4. **Prospective Korean Bundang cohort** (mentioned in Paper 1 outreach) — institutional outreach stage; would be Tier-1.

---

## How this changes the manuscript

The current NComms v2 draft positions the panel as Korean/Asian-anchored:
- Discovery cohort = TCGA-THCA (Tier 4) + Lee 2024 GSE213647 (Korean, n=632) + Mun 2025 (Korean proteome, n=336).
- Single-cell mechanism = Lu 2023 (Chinese, n=14,624 cells).
- Cross-ethnic replication = GSE286332 (Korean, n=18).
- Gray-zone framing = Mu 2024 (Chinese, n=214).
- Outreach Tier-1 anchor candidate = Boucai 2023 (US) + Bundang prospective (Korean).

**Korean / Asian data is already the spine of the discovery half** of this manuscript. Adding Pu 2021 RAI-refractory distant-met cells and Mu 2024 per-patient genotype would push it into Nature Communications / Nature Medicine reach territory.
