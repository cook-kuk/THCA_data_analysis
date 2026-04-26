# Top 20 novel candidate biomarkers — PubMed deep dive

These genes are present in our 2,773-replicated biomarker set but absent from
all 10 reference gene lists (Chakravarty BRS52, Landa BRS, TCGA THCA 2014,
Yoo 2019, Costa 2015, Pu 2021 scRNA, Agrawal 2014, COSMIC CGC thyroid,
OncoKB thyroid, DisGeNET thyroid carcinoma).

Per-gene PubMed counts pulled live from NCBI Entrez. ` < 5 thyroid papers ` = potentially novel.

| # | Gene | log2FC | FDR | thyroid PMIDs | cancer PMIDs | BRAF PMIDs | classification |
|---|------|-------:|-----|--------------:|-------------:|-----------:|---------------|
| 1 | **DCSTAMP** | +5.35 | 1.7e-82 | 5 | 18 | 2 | emerging |
| 2 | **KCNN4** | +3.78 | 4.5e-59 | 5 | 63 | 0 | emerging |
| 3 | **PLEKHA6** | +2.35 | 3.0e-54 | 0 | 5 | 2 | novel |
| 4 | **CYP1B1** | +3.47 | 9.0e-53 | 16 | 1008 | 4 | emerging |
| 5 | **BNC1** | +1.65 | 5.0e-66 | 1 | 34 | 2 | novel |
| 6 | **LDLR** | +2.24 | 2.8e-53 | 44 | 357 | 3 | established |
| 7 | **GABRB2** | +3.70 | 2.9e-47 | 9 | 18 | 2 | emerging |
| 8 | **ST6GALNAC5** | +2.91 | 1.5e-61 | 0 | 24 | 0 | novel |
| 9 | **B3GNT3** | +3.78 | 8.0e-52 | 0 | 44 | 0 | novel |
| 10 | **PTPRE** | +2.38 | 3.8e-45 | 4 | 17 | 1 | novel |
| 11 | **KCNQ3** | +2.48 | 5.0e-51 | 3 | 13 | 0 | novel |
| 12 | **CREB5** | +2.09 | 3.6e-49 | 4 | 53 | 3 | novel |
| 13 | **LY6E** | +2.13 | 5.8e-43 | 0 | 49 | 0 | novel |
| 14 | **TAGLN2** | +1.22 | 2.3e-42 | 5 | 76 | 0 | emerging |
| 15 | **PDLIM4** | +3.44 | 3.1e-40 | 4 | 37 | 2 | novel |
| 16 | **ITGA3** | +1.52 | 5.0e-42 | 16 | 170 | 2 | emerging |
| 17 | **BID** | +1.35 | 1.7e-41 | 73 | 2598 | 44 | established |
| 18 | **SPOCK2** | +2.75 | 9.4e-49 | 0 | 28 | 1 | novel |
| 19 | **CST6** | +3.87 | 1.2e-41 | 3 | 56 | 1 | novel |
| 20 | **SYT12** | +4.88 | 1.1e-66 | 6 | 11 | 0 | emerging |

## Summary of top-20 novel-candidate landscape

- **11/20** genes have < 5 thyroid-cancer papers — these are the strongest novelty candidates
  surfaced by THYRAI's biomarker pipeline that the field has not yet flagged.
- **7/20** are emerging (5–24 thyroid papers).
- **2/20** are established in the broader thyroid literature but were missed by
  the 10 reference signature lists — suggesting these signature lists may themselves be
  incomplete and our pipeline is filling that gap.

## Claim

*Of the 20 highest-ranked novel candidates, 11 have fewer than five thyroid-cancer
papers indexed in PubMed. These represent genuinely under-investigated loci surfaced by
the THYRAI BRAF-vs-RAS differential expression pipeline — candidate additions to the
thyroid molecular subtype signature space whose biological role merits prospective study.*
