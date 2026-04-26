# Preprint scoop scan — Europe PMC (bioRxiv, medRxiv, Research Square, ChemRxiv)

_Date: 2026-04-25 · Author: Seungho Cook · Source: Europe PMC `SRC:PPR` filter, 10 queries._
_Companion to_ `novelty_assessment.md` _and_ `sn38_braf_priorart.md`.

## Why this scan

PubMed E-utilities (Task 1–2) covers indexed journal literature only. Preprint servers (bioRxiv, medRxiv, Research Square, ChemRxiv) are where a competing group would post first. This pass closes the **scoop-risk** gap before paper submission.

## Queries and hit counts

| Label | Query | Preprints |
|---|---|---:|
| P1 | `(BRAF AND TROP2 AND thyroid) AND SRC:PPR` | 0 |
| P2 | `(BRAF AND TACSTD2 AND thyroid) AND SRC:PPR` | 0 |
| P3 | `(V600E AND TROP2) AND SRC:PPR` | 1 |
| P4 | `(thyroid AND sacituzumab) AND SRC:PPR` | 0 |
| P5 | `(thyroid AND TROP2 AND (ADC OR antibody-drug conjugate)) AND SRC:PPR` | 0 |
| P6 | `(BRAF AND irinotecan AND thyroid) AND SRC:PPR` | 0 |
| P7 | `(BRAF AND topotecan AND thyroid) AND SRC:PPR` | 1 |
| P8 | `(PRISM AND BRAF AND thyroid) AND SRC:PPR` | 3 |
| P9 | `(thyroid AND TROP2) AND SRC:PPR` | 0 |
| P10 | `("papillary thyroid" AND TROP2) AND SRC:PPR` | 0 |
| **Total** | | **5** |

## Triage of the 5 preprint hits

| # | Year | Repo | Title | Verdict |
|---|---|---|---|---|
| 1 | 2025 | medRxiv | A clinically relevant morpho-molecular classification of lung neuroendocrine tumours | **No conflict** — lung NET, V600E + TROP2 incidentally co-mentioned in different sections. |
| 2 | 2024 | bioRxiv | Connectivity Map and perturbation-based sensitivity analysis identifies MEK inhibitors as senolytics in human lung fibroblasts | **No conflict** — lung fibroblast senescence, not thyroid cancer. |
| 3 | 2025 | bioRxiv | Repurposing passenger amplifications as Trojan horses identifies MPZL1 as a potent target for solid cancers | **No conflict** — pan-cancer MPZL1 target via amplification; uses PRISM but not BRAF×TROP2 axis. |
| 4 | 2025 | Research Square | Lung tumoroids as a testing platform for precision CAR T-cell therapy | **No conflict** — CAR T in lung tumoroids. |
| 5 | 2021 | medRxiv | Leveraging sequences missing from the human genome to diagnose cancer | **No conflict** — cfDNA neomer cancer detection. |

## Verdict

**No preprint scoop risk for the v14 BRAF × TROP2 × thyroid ADC stratification thesis as of 2026-04-25.**

The closest neighbour preprint (Hit #3, MPZL1 passenger amplification) uses PRISM data but in a pan-cancer amplification framework distinct from our BRAF-subtype × ADC-target axis. Hit #1 (lung NET) shares two of our keywords (V600E, TROP2) but in a different tumour and different mechanistic narrative.

## Caveats

- Europe PMC preprint coverage indexes the major repos (bioRxiv, medRxiv, Research Square, ChemRxiv, SSRN, Wellcome Open Research) but **lag of indexing** ranges 1–14 days.  Re-run within 48 h of submission.
- Conference abstracts (ASCO, ESMO, AACR, ATA) are **not** in this filter and have not been swept here. ATA 2025 / ASCO 2026 abstracts should be hand-checked before submission.
- Patent landscape: see `patent_scan_status.md`.

## Reproducibility

- Script: `notebooks_or_scripts/v14_preprint_patent_scan.py`
- Raw output: `results/v14_priorart/preprint_raw.json`
- Flat TSV: `results/v14_priorart/preprint_scan.tsv` (5 rows)
