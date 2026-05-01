# D8-D — 5-Pillar Paper Figure Inventory + Integration Plan

**Date:** 2026-05-02 (D8 wrap-up)
**Goal:** Map existing 100+ interactive figures + tables to the 5-Pillar paper structure; identify what's main / suppl / new figure needs.

---

## ★ Existing figure inventory (rough count)

| Category | Count | Location |
|---|---|---|
| HLA mega (deep + final + cross + KOR + pool) | 30+ | `reports/html/figs_interactive/v17/hla_mega/` |
| K2/Korean dashboards | 8+ | `reports/html/figs_interactive/v17/Korean_K2_*.html` |
| 8-gene panel figures | 6+ | `reports/html/figs_interactive/v17/v17_8gene_*.html` |
| Lu 2023 multi-cohort | 12+ | `reports/html/figs_interactive/v17/v17_lu2023_*.html` |
| TERT/BRAF/age | 8+ | `reports/html/figs_interactive/v17/v17_*tert*.html`, etc |
| ML / panel-size / forest | 10+ | `reports/html/figs_interactive/v17/v17_q*.html` |
| Wang 2023 population | 1 | `v17_wang2023_population.html` |
| Pozdeyev landscape | 1 | `v17_pozdeyev_landscape.html` |
| Multi-section HTML pages | 19 | `reports/html/pages/0X_*.html` |

**Total:** ~100+ interactive HTML figures + 19 long-form HTML pages.

---

## 1. Pillar 1 — Korean Pan-Asian HLA cohort (n=874)

### Main figure: Pan-Asian HLA forest meta
- **NEW (need to build):** Forest plot with 3 arms
  - Korean PTC pool n=874 (K2 235 + Lee 630 + GSE286332-PTC 9)
  - GSE286332 PTC+HT n=9 (exploratory)
  - Chen 2018 Han Chinese GD (case=1,468 ctrl=1,490)
  - 5 focus alleles: DPB1\*05:01, B\*46:01, DRB1\*15:01, DRB1\*03:01, DRB1\*04:01
  - Wilson 95% CI per allele, OR axis
  - Source data: `results/d4p1_panasian_meta/panasian_forest_meta.tsv`

### Supplementary figures (existing, ready)
- `hla_mega/FINAL_01_headline.html` — All-in-one HLA dashboard
- `hla_mega/FINAL_02_carrier_forest.html` — K2+Lee carrier forest
- `hla_mega/POOL_01_meta.html` — pooled meta
- `hla_mega/DEEP_01_pan_asian.html` — Pan-Asian replication (existing)
- `hla_mega/FINAL_05_allele_freq_grid.html` — Allele frequency grid
- `hla_mega/KOR_01_top_alleles.html` — Top Korean alleles

### Suppl tables
- `results/d4p1_panasian_meta/korean_PTC_pool_n908.tsv` (per-sample 4-digit alleles)
- `results/v17_korean/arcasHLA_GSE286332/{18 × genotype.json}`

---

## 2. Pillar 2 — GSE286332 PTC vs PTC+HT molecular dissection

### Main figure: Multi-panel autoimmune-PTC signature
- **Panel A:** Volcano plot (10,380 DEGs, top labels: IGHV3-66, BLK, HLA-DOB, etc)
- **Panel B:** GSEA bar plot (Hallmark IFN-γ + KEGG Type I diabetes + Reactome TCR)
- **Panel C:** 8-gene RAI score boxplot (PTC vs PTC+HT, Cohen d=−1.60)
- **Panel D:** HLA-II module heatmap (per-sample × HLA-II 10 genes)
- Source data: `results/p3_gse286332/{deg_*.tsv, gsea_*.tsv, 8gene_panel_per_sample.tsv, hla_module_scores.tsv}`

### Supplementary
- Per-gene 8-gene comparison: PAX8 d=−2.32, etc (table)
- DM1/DM2 prediction distribution (18/18 → DM2)
- AICDA + TLS module scores (D5-P6 results)

---

## 3. Pillar 3 — Driver mRNA neutrality

### Main figure (small): driver mRNA × mutation status
- **Panel A:** BRAF mRNA boxplot (V600E vs WT, Cohen d=−0.04, p=0.57)
- **Panel B:** RAS-family mRNA × RAS-mutation
- **Panel C:** Driver single-feature AUC bar (BRAF 0.602, TERT 0.578, KRAS 0.525, NRAS 0.521, HRAS 0.500)
- **Panel D:** TIERA67 univariate Cohen d ranking (8-gene labels at #4-50, drivers at #52-67)
- Source data: `results/p1_driver_mrna_audit/`

### Existing figures
- `v17_external_BRAF_TERT.html` — external BRAF/TERT validation
- `v17_quad_tert_stack.html` — TERT stratification
- `v17_8gene_cross_analysis.html` — 8-gene cross-validation

---

## 4. Pillar 4 — Pan-genome cluster robustness

### Main figure (sup): 8-panel ARI table
- Bar plot: 8-gene 0.49, TIERA67 0.90, top-5000 0.92, Driver-only 0.00
- Hypergeometric p (TIERA67 enriched in pan-genome top 100): 3e-4
- Source: `results/p4_pangenome_vs_tiera67/ari_comparison.tsv`

### Existing figures
- `v17_panel_size_sensitivity.html` — 8/10/12/16 panel comparison
- `v17_q21_AUC_matrix.html` — multi-panel AUC matrix

---

## 5. ★ Pillar 5 — Autoimmune-PTC mechanism layer (NEW, paper-changing)

### Main figure: 4-panel mechanism narrative
- **Panel A:** GSE286332 P_DM1 spectrum + mediation diagram (HLA-II 140% mediation)
  - Source: `results/d3p5_pdm1_gradient/`
- **Panel B:** TCGA Hashimoto-like × DM cluster forest plot
  - 4 thresholds (GMM, Otsu, top20, top30) × DM2 enrichment OR
  - Top30: OR=0.20, p=6e-10
  - Source: `results/d4p2_tcga_hashimoto_signature/`
- **Panel C:** BCR clonal + TLS heatmap
  - 18 sample × {IGHV clonality, TLS score, AICDA expression}
  - Source: `results/d5p6_bcr_repertoire/`
- **Panel D:** DM1 sub-A vs sub-B mutation × signature stacked bar
  - sub-A: 69% RAS+, sub-B: 96% mut-neg + 4× Hashimoto-like
  - Source: `results/d6p7_dm1_subcluster/`

### Supplementary
- TCGA Hashimoto-like × age, sex, stage cross-tab
- HLA-II d residualization table (full vs Hashimoto-excluded vs within-Hashimoto)
- Korean GSE213647 replication: 22.8% Hashimoto-like (Otsu 28.2%) ≈ TCGA 18-20%
- Korean sub-B-like rate 47-53% (D8-C signature transfer)

---

## 6. New figures (need to build)

| # | Figure | Status | Pillar | Tool |
|---|---|---|---|---|
| F1 | Pan-Asian HLA forest meta (3-arm) | Build | 1 | plotly forest |
| F2 | GSE286332 multi-panel (A/B/C/D) | Build | 2 | matplotlib subplots |
| F3 | Driver mRNA neutrality summary | Build | 3 | plotly + table |
| F4 | Pan-genome ARI bar | Build | 4 | plotly bar |
| F5 ★ | Autoimmune-PTC mechanism 4-panel | Build | 5 | matplotlib + ven |
| F6 | TCGA Hashimoto-like distribution | Build | 5 | violin + density |
| F7 | DM1 sub-A vs sub-B stacked phenotype | Build | 5 | stacked bar |
| F8 | Cross-cohort Hashimoto-like prevalence | Build | 5 | bar comparison |

---

## 7. Existing 19 HTML page sections — which to repurpose for paper

| Page | Content | Paper relevance |
|---|---|---|
| `00_all_sections.html` | Master overview | Discussion intro |
| `01_overview.html` | Project overview | Discussion |
| `02_datasets.html` | Cohort table | Methods (Suppl Table 1) |
| `03_sample_master.html` | Sample-level master | Suppl Table 2 |
| `04_gene_panels.html` | 8-gene + TIERA67 | Methods + Pillar 3 |
| `05_eda.html` | Exploratory | Suppl |
| `06_scores.html` | Score profiles | Pillar 1-5 various |
| `07_ml_baseline.html` | ML benchmarking | Pillar 4 |
| `08_panel_comparison.html` | 8/10/12/16 | Pillar 4 |
| `09_shap.html` | SHAP attribution | Pillar 3 (driver neutrality) |
| `10_gene_explorer.html` | Per-gene | Suppl |
| `11_cohort_compare.html` | Multi-cohort | Pillar 1, 2, 5 |
| `12_business.html` | Clinical applicability | Discussion |
| `13_reports.html` | Sub-reports | Suppl directory |
| `14_caveats.html` | Limitations | Discussion |
| `15_drug_discovery.html` | RAI / agent targets | Discussion (forward) |
| `16_quantum.html` | Quantum-inspired panel selection | Methods footnote |
| `17_biomarker_insights.html` | Biomarker validation | Pillar 1-5 |
| `18_pathway_immune_meth.html` | Pathway analysis | Pillar 5 |
| `19_cmap_network.html` | CMap drug network | Discussion forward |

---

## 8. Build order priority (D9–D11 if user wants)

Priority 1 (immediate, paper-shaping):
- **F5**: Autoimmune-PTC mechanism 4-panel ★
- **F2**: GSE286332 multi-panel
- **F1**: Pan-Asian HLA forest

Priority 2 (table-heavy supp):
- **F6**: TCGA Hashimoto-like distribution
- **F7**: DM1 sub-A/B stacked phenotype

Priority 3 (small supp):
- **F3**: Driver neutrality summary
- **F4**: Pan-genome ARI bar
- **F8**: Cross-cohort prevalence

---

## 9. Outputs

This file: `project/reports/2026_05_03_D8D_figure_inventory.md`
Existing inventory artifacts: ~100+ HTML in `reports/html/figs_interactive/v17/` + 19 sections in `reports/html/pages/`
