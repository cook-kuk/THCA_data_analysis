# Manual download required · sources we cannot auto-fetch

Generated: 2026-05-07
Maintainer: Lumenix · Seungho Cook

This file lists data sources where automated download fails or is forbidden.
For each: exact URL, what to download, where to place it, format expected.

---

## 1. TESLA · Cell paper supplementary tables

**Status:** PARTIAL — GitHub portion is reachable (`https://github.com/ParkerICI/tesla`); the Cell paper supplementary tables (especially Table S4 / S7 with the full 608 pMHC list + immunogenicity labels) require manual download from Cell Press.

**Action:**
1. Visit https://www.cell.com/cell/fulltext/S0092-8674(20)31170-5 (or via institutional access).
2. Download all supplementary tables (Excel: typically `mmc1.xlsx` through `mmc7.xlsx`).
3. Save to `/data/neoantigen_vaccine_hub/data_raw/tesla/cell_supp_tables/`.
4. Run `python scripts/parse_tesla_supp.py` after dropping in.

**Citation:** Wells et al. 2020 Cell, PMID 32169175, DOI 10.1016/j.cell.2020.09.015.

---

## 2. dbPepNeo / dbPepNeo2.0

**Status:** BROKEN_URL — `www.biostatistics.online` DNS does not resolve as of 2026-05-07.

**Recovery options:**
1. Try alternate URL: https://web.archive.org/web/2024*/dbpepneo.org/
2. Email the original authors (corresponding contact in dbPepNeo paper).
3. Try the GitHub mirror if one exists: search "dbPepNeo" on GitHub.
4. Check if it is now hosted under the IEDB CEDAR umbrella (some Chinese databases have migrated).

**Files needed:**
- `dbPepNeo_LC.tsv` (LC = MS-only confidence)
- `dbPepNeo_MC.tsv` (MC = MS + WES/WGS)
- `dbPepNeo_HC.tsv` (HC = T-cell-validated)

Expected schema: peptide, hla, tumor_type, gene, mutation, validation_method, reference.

---

## 3. Neodb

**Status:** BROKEN_URL — `nat.zju.edu.cn/neodb` DNS does not resolve as of 2026-05-07.

**Recovery options:**
1. Search bioRxiv / GitHub for "Neodb" — paper authors often deposit a release on Zenodo.
2. Email the corresponding author of the Neodb paper.
3. Try mirror: https://pgx.zju.edu.cn/neodb (different subdomain at the same institution).

**Files needed:**
- `neodb_val_neo.csv` (Val-Neo = experimentally validated immunogenic)
- `neodb_driver_neo.csv` (Driver-Neo = recurrent driver-mutation-derived shared candidates)

---

## 4. NEPdb

**Status:** ATTEMPTED — site reachable (HTTP 200) but bulk download path unclear.

**Recovery procedure:**
1. Visit https://nep.whu.edu.cn/ (use Chrome/Firefox; site may render in Chinese).
2. Look for a "Download All" button or a link to a TSV/Excel under "Data" or "Bulk Download".
3. Save to `/data/neoantigen_vaccine_hub/data_raw/nepdb/`.
4. Required files: `nepdb_positive.tsv`, `nepdb_negative.tsv` (positive AND negative T-cell-validated neoepitopes).

If the download is gated behind a registration form, register and document credentials in a separate `.env` file (NEVER commit credentials).

---

## 5. TSNAdb v2.0

**Status:** ATTEMPTED — site reachable (HTTP 200), bulk-download path:

**Procedure:**
1. Visit https://pgx.zju.edu.cn/tsnadb/.
2. Navigate to "Download" tab. Look for SNV / INDEL / Fusion-derived neoantigen files.
3. Save raw files to `/data/neoantigen_vaccine_hub/data_raw/tsnadb/`.
4. Files we need (typically): `tsnadb_snv_predicted.tsv`, `tsnadb_snv_validated.tsv`, `tsnadb_indel.tsv`, `tsnadb_fusion.tsv`.

**Critical:** Predicted vs validated subsets MUST be kept separate. Our pipeline assigns LEVEL_0 vs LEVEL_4 respectively — do NOT merge.

---

## 6. McPAS-TCR

**Status:** ATTEMPTED — typically requires a form submission (email address) before download.

**Procedure:**
1. Visit https://friedmanlab.weizmann.ac.il/McPAS-TCR/
2. Click "Download data" → fill in name + email + research purpose.
3. Save the resulting CSV to `/data/neoantigen_vaccine_hub/data_raw/tcr/mcpas/`.

**Field warning:** Do NOT redistribute — academic use only per the Friedman lab terms.

---

## 7. PRIME / MixMHCpred binaries

**Status:** PARTIAL — the GitHub repo is cloneable (data files included), but the predictor binaries are gated by a license check.

**Procedure for license:**
1. Visit https://github.com/GfellerLab/PRIME → README → license link.
2. Submit institutional information for academic-use license.
3. Save the binary to `/data/neoantigen_vaccine_hub/scripts/external/prime/`.

For our benchmark, we use only the **published scores reported in PRIME paper** — we do not re-run PRIME on our own peptides until the license is obtained. Document the score source as cited.

**License flag:** BLOCKED_BY_LICENSE for commercial use.

---

## 8. IMPROVE supplementary

**Status:** PROBABLE — Frontiers in Immunology papers usually have CC-BY supplementary tables.

**Procedure:**
1. Visit the IMPROVE paper at frontiersin.org (Borch et al. 2024).
2. Download all supplementary tables (Excel) under "Supplementary Material".
3. Save to `/data/neoantigen_vaccine_hub/data_raw/improve/`.
4. Required: 17,500 tested neoepitope table + 467 T-cell-recognized subset.

---

## Summary table

| Source | Status | Manual? | Why |
|---|---|---|---|
| TESLA | PARTIAL | YES (Cell suppl) | Cell paper paywall on suppl |
| CEDAR | SUCCESS | NO | API works |
| IEDB | SUCCESS | NO | API works |
| BigMHC | SUCCESS | NO | Mendeley + GitHub |
| NeoRanking | SUCCESS | NO | GitHub clone |
| IMPROVE | PROBABLE | YES | suppl xlsx |
| NEPdb | ATTEMPTED | YES (form) | site OK but no bulk URL |
| dbPepNeo | BROKEN_URL | YES | DNS down |
| TSNAdb | ATTEMPTED | YES | download button on site |
| Neodb | BROKEN_URL | YES | DNS down |
| PRIME | BLOCKED_BY_LICENSE | YES | commercial license |
| VDJdb | SUCCESS | NO | GitHub clone |
| McPAS-TCR | ATTEMPTED | YES (form) | email registration |
| NeoTCR | PROBABLE | optional | GitHub fallback |
| TCR3d | PROBABLE | optional | bulk index download |
| clinical vaccine | curated list | YES | not raw data — paper citations |
