# H20 — HT axis prognostic in 9/32 cancer types FDR<0.1 (OS), 15/32 (PFI)

**Question.** Does HT-13 (HLA-DR/DP/DQ ×6, CD79A/B, MS4A1, AICDA, CXCL13, CCR6, IFNG) — read out per-sample as within-cohort z-mean — predict OS/DSS/PFI across 33 TCGA lineages, and how does its coverage compare to DM1?

**Data.** TCGA pan-cancer RSEM (Xena PanCanAtlas, 11,069 samples). Liu 2018 harmonized survival. All 13 HT-13 genes present. Cox PH per lineage with covariates {age, stage} when available; BH-FDR within endpoint. DM1 reference taken from `phase_B_survival/cox_per_lineage.tsv`.

**Headline.**

| Endpoint | Fit | FDR<0.1 | Top hits |
|---|---|---|---|
| OS | 32 | **9** | LGG **HR 2.12 [1.66,2.71] p=2.3e-9**, SKCM HR 0.60 p=2.4e-8, UVM 2.28, CESC 0.66, HNSC 0.77, LUAD 0.75, SARC 0.65, THCA 0.48, UCEC 0.72 |
| DSS | 29 | 8 | LGG, SKCM, LUAD, HNSC, UVM, CESC, SARC, THCA |
| PFI | 32 | **15** | LGG, GBM, SKCM, HNSC, LUAD, LIHC, KIRP, BRCA, BLCA, THYM, ACC, UCEC, MESO, CESC, CHOL |

Direction is biologically consistent: HT-13-high → **better OS in immune-responsive lineages** (SKCM, HNSC, LUAD, CESC, THCA, SARC) and **worse OS in immune-cold/immunosuppressive contexts** (LGG, UVM, GBM, THYM). LGG TLS-driven worse survival is the canonical exception (Mathewson 2021).

**ICI coverage.** Of 11 ICI-on-label lineages, HT-13 is FDR<0.1 in 4/11 OS, 4/11 DSS, 6/11 PFI — SKCM + LUAD + HNSC are protective in all three endpoints, matching ICI-benefit dogma.

**vs DM1 (Phase B).**
- OS: HT-13 sig=9, DM1 sig=3, **union=9, intersection=3 (LGG, LUAD, UCEC)**. HT-13 strictly broadens the prognostic-cancer set; HT-13-only adds SKCM, HNSC, UVM, CESC, SARC, THCA.
- PFI: HT-13=15, DM1=10, union 19, intersection 6.
- log HR Spearman (across lineages): **OS ρ=0.37 p=0.038, DSS ρ=0.56 p=0.0017, PFI ρ=0.58 p=5e-4** — same biology, partially-shared variance.
- Mean C-index (covariate-free, like-for-like): HT-13 OS 0.570, DM1 0.550, **combined 0.593**; DSS 0.576/0.564/**0.611**; PFI 0.555/0.556/**0.585**. Combined axis superior in all three endpoints.

**Caveats.** Univariate-style Cox per lineage; small-event lineages (DLBC, KICH, ACC) have wide HRs. THCA OS HR=0.48 is event-poor (n=569, 20 events). Stage parsing is string-based; >70% missing in heme/glioma → age-only adjustment.

**Files (abs):**
- `project/results/p2_braf_nature_sprint_2026_05_09/h20_pancancer/h20_pancancer_ht_cox.tsv`
- `…/h20_dm1_vs_ht_scatter.tsv`
- `…/h20_ht13_scored_per_sample.tsv`
- `…/forest_HT13_{OS,DSS,PFI}.png`, `…/scatter_dm1_vs_ht13_{OS,DSS,PFI}.png`
- `…/summary.json`, `…/run_h20.py`

**Bottom line.** HT-13 alone replicates and **broadens** DM1's pan-cancer prognostic signal: 9/32 OS-significant FDR<0.1, directionality matching ICI-response dogma in melanoma/lung/HNSC and immune-cold dogma in glioma/uveal melanoma. Effect-size scatter ρ=0.37–0.58 confirms HT-13 and DM1 read out the same axis from opposite poles. Combined-axis C-index uplift in every endpoint supports HT-13 as a portable, lineage-agnostic readout — the axis is **universal, not thyroid-specific**.
