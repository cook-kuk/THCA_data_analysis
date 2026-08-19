# Liu 2024 — RR-PTC primary-tumour proteomics

**Added** 2026-08-06 · obtained by the author directly (the medRxiv host blocks scripted access)

Liu H, Wang J, Zhou Y, Kong D, Hu P, Li L, Yang D, Xu Z, Sun Y, Chen C.
*Molecular and Proteomic Profiles of Radioiodine Refractory Papillary Thyroid Cancer.*
medRxiv 2024.09.22.24314143 (v4, 2024-12-12). Renmin Hospital of Wuhan University +
Westlake University. Corresponding: Chuang Chen `chenc2469@whu.edu.cn`,
Yaoting Sun `sunyaoting@westlake.edu.cn`. **Still a preprint — no journal version as of 2026-08-06.**

## Local location

```
external_data/liu_rrptc_2024/
├── source/          liu_rrptc_original_bundle.zip            (17.4 MB, immutable)
│                    liu_rrptc_supplementary_methods_figures.pdf  (19.9 MB)
│                    liu_rrptc_supplementary_tables.xlsx      (35 KB, Tables S1–S7)
├── extracted/       liu_rrptc_calculation_data.xlsx          (14.7 MB, 26 sheets)
│                    liu_rrptc_additional_files.pdf           (2.8 MB)
├── processed/       liu_rrptc_primary_48_analysis_table.tsv
└── manifest.sha256
```

The two PDFs have different SHA-256 hashes and are different documents — both retained.
`source/` and `extracted/` are gitignored; only code and result tables are tracked.

## What the workbook contains

26 sheets. Full inventory at `results/liu_rrptc_2024/workbook_inventory.tsv`. The four that matter:

| Sheet | Unit | Size | Content |
|---|---|---|---|
| `6-1` | patient | 73 × 33 | clinical — group, RAI uptake, disease persistence, loss of uptake at first I-131, cumulative I-131 dose, stage, recurrence risk, follow-up |
| `8-1` | sample | 48 × 69 | per-sample computed scores — **TDS**, TumorPurity, ERK score, TIS, IIS, CYT, CIBERSORTx, ssGSEA |
| `8-2` | sample | 168 × 6 | ESTIMATE purity for all proteomic samples |
| `9-1` | gene | 141 × 6 | gene definitions for TDS, BRAF-RAS score, ERK score, cytolytic score, with a per-gene "identified in the current analysis" flag |

## Critical limitation: no per-protein abundance matrix

All 26 sheets were checked. The 7,394-row sheets (`2-1`, `3-1`, `4-1`, `5-3`) are
**differential-expression summary statistics** (logFC, AveExpr, t, P, adj.P), not sample-level
values. The only sheets with sample columns are GSVA matrices (`2-3`, `3-3`, `4-3`;
1,733 gene sets × samples) and the score sheets above.

**Consequence: our eight-gene or five-gene overlap score cannot be recomputed here.**
Only the authors' own TDS can be re-analysed. To compute our panel we would need the
normalised protein abundance matrix (≈7,394 proteins × 168 samples) plus sample annotation
from iProX.

## Relationship to our eight-gene panel

The authors' TDS is defined over 16 genes but **only 7 were detected in the proteome**:
DUOX1, DUOX2, FOXE1, NKX2-1, PAX8, TG, TPO.

| Our panel gene | In TDS definition | Detected in proteome |
|---|---|---|
| TPO, TG, PAX8, NKX2-1, FOXE1 | yes | **yes** (5 shared) |
| **SLC5A5 (NIS)** | yes | **no** |
| **TSHR** | yes | **no** |
| **DIO1** | yes | **no** |

The sodium-iodide symporter — the gene with the most direct mechanistic claim on iodide
transport — is not measurable at protein level here, the same limitation as the HTG EdgeSeq
panel used in E-MTAB-12837/12900.

## Cohort

| Stage | n | Detail |
|---|---|---|
| Full clinical series | 73 | RAIR 60 / Control 13 |
| Primary lesion unavailable | 25 | |
| **Included in the current analysis** | **48** | RAIR 37 / Control 11 |
| RAI uptake | 48 | absent 27 / present 21 |
| Four behaviour subgroups | 48 | Id 11 · ID 10 · iDF 9 · iDG 18 |

Subgroup rules reconstructed from the clinical fields and reproducing the published counts
exactly: **Id** = Control; **ID** = RAIR with uptake present; **iDF** = RAIR, uptake absent,
loss at first I-131 = Yes; **iDG** = RAIR, uptake absent, loss at first = No.

## Circularity found

`Disease persistence` is Yes 60 / No 13, **identical to** `Group` (RAIR 60 / Control 13).
It is a component of the group definition, not an independent endpoint, and must not be used
as one.

## Results (published TDS, n = 48)

| Test | n | Effect | P | Verdict |
|---|---|---|---|---|
| RAIR vs Control | 37 / 11 | d = −0.37 | 0.33 | null |
| **RAI uptake absent vs present** | 27 / 21 | **d = −0.51** | **0.088** | trend |
| Uptake absent vs present, within RAIR | 27 / 10 | d = −0.49 | 0.17 | null |
| Four subgroups (Kruskal-Wallis) | 48 | H = 2.97 | 0.40 | null |
| TDS × tumour purity | 48 | ρ = +0.15 | 0.30 | **no confounding** |
| TDS × cumulative I-131 dose | 48 | ρ = +0.16 | 0.27 | unrelated |
| **TDS × ERK output score** | 48 | **ρ = −0.42** | **0.0029** | significant |

Subgroup medians decline monotonically in the hypothesised order —
Id +0.183 → ID +0.119 → iDF −0.084 → iDG −0.280 — but the test is null.

## What this supports and what it does not

**Supports.** Independent proteomic corroboration that primary-tumour differentiation state
alone does not classify radioiodine behaviour well (RAIR vs Control P = 0.33), consistent with
the TCGA initial-response null. The uptake contrast is larger than the group contrast
(−0.51 vs −0.37), matching our TCGA finding that iodine-behaviour phenotypes are better
endpoints than binary treatment response. The TDS–ERK inverse correlation is mechanistically
coherent with the MAPK–differentiation axis seen in our own deconvolution work.

**Does not support.** This is **not** validation of our eight-gene panel — it is the authors'
seven-protein score, missing NIS, TSHR and DIO1, and the underlying protein values are not
available. Forbidden phrasings: "direct validation of our eight-gene panel",
"treatment-benefit predictor", "definitive RAI-specific biomarker", "no effect was proven".

Acceptable phrasing: *"The published seven-protein thyroid differentiation score showed no
significant difference between radioiodine-refractory and control patients in a small
primary-tumour proteomic cohort (n = 48, P = 0.33), with a non-significant trend by
radioiodine uptake status (d = −0.51, P = 0.088)."*

## Outstanding request

From iProX, the minimum useful files would be: the normalised gene-symbol protein abundance
matrix, sample annotation, technical-replicate/QC annotation, and the patient-to-sample map.
Raw DIA files are not needed.

## Files generated

- `results/liu_rrptc_2024/workbook_inventory.tsv`
- `results/liu_rrptc_2024/rai_subgroup_manifest.tsv`
- `results/liu_rrptc_2024/published_TDS_effects.tsv`
- `processed/liu_rrptc_2024/liu_rrptc_primary_48_analysis_table.tsv`
- `figures/liu_rrptc_2024/figure_liu_tds_2026_08_06.{png,pdf}`
