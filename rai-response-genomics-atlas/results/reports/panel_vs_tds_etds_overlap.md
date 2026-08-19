# 8-gene panel vs TDS-16 vs eTDS-64 — containment + parsimony

## Source data
- 8-gene panel (Paper 1 canonical): `config/eight_gene_panel.yaml`.
- TDS-16 (Landa 2016, originally TCGA-THCA Cell 2014 Table S5): `project/metadata/tds16_genes.txt`.
- eTDS-64 (Boucai 2023 CCR Supp Table S4, PMC10106408): the explicit gene roster is JS-gated on PMC and was not retrievable in this session. The containment claim relies on Boucai's own Methods statement that eTDS extends TDS-16.

## Containment relations
- 8-gene ∩ TDS-16 = **8 / 8** genes — `DIO1, FOXE1, NKX2-1, PAX8, SLC5A5, TG, TPO, TSHR`.
- 8-gene \ TDS-16 = **0 / 8** — `none`.
- TDS-16 \ 8-gene = **8 / 16** — `DIO2, DUOX1, DUOX2, GLIS3, SLC26A4, SLC5A8, THRA, THRB`.
- TDS-16 ⊂ eTDS-64 (per Boucai 2023 Methods).
- Therefore 8-gene ⊂ TDS-16 ⊂ eTDS-64 (by transitivity).

## Parsimony (Paper 1 R8 audit on TCGA-THCA)
- 8-gene Cohen d (DM1 vs DM2) = **1.783**, AUC = **0.903**.
- TDS-16 Cohen d = **1.973**, AUC = **0.913**.
- 8-gene captures **98.9%** of TDS-16's AUC discriminative power at **50% of the gene cost** (8 vs 16 genes).
- Within-sample Spearman ρ(8-gene panel z, TDS-16 score) = **0.954** (Paper 1 memory `DM1 Round 8`).

## Manuscript argument
1. The 8-gene panel is a **parsimonious subset of the published TDS-16** that the field already uses, not an arbitrary cherry-pick.
2. Boucai 2023 CCR built **eTDS-64** on top of TDS-16; our panel therefore lies inside the same iodide-handling signature space that Boucai's exceptional-responder analysis used.
3. The 8 TDS-16 genes NOT in our panel (DIO2 · DUOX1 · DUOX2 · GLIS3 · SLC26A4 · SLC5A8 · THRA · THRB) form the **extended iodide axis** that can be added for sensitivity analysis without changing the headline score.
4. The AUC parity (0.903 vs 0.913) supports presenting the 8-gene panel as the *operating* score with TDS-16 as the *robustness baseline* in the manuscript.

## Caveats
- The explicit eTDS-64 gene roster (Boucai Supp Table S4) was NOT retrieved in this session. The PMC supplement is downloadable only via a JS-driven 'preparing to download' interstitial that defeats curl-style fetches. Either (a) manually download via browser, (b) request from corresponding author, or (c) check whether Landa 2016 / TCGA 2014 enumerates a comparable enhanced-iodide-axis set we can re-derive.
- The Cohen d / AUC numbers above are from Paper 1 R8 audit on TCGA-THCA, a Tier-4 discovery cohort. Numbers will differ in Tier-1/2 anchored cohorts (Boucai, Mu, GSE151179).