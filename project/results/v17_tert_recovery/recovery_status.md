# v17 TERT promoter recovery — status

_Run: 2026-04-27T00:24:24Z_  
_Wall: 3.5 s_

## Source-by-source result

### cBioPortal datahub
- URL(s): `https://github.com/cBioPortal/datahub/raw/master/public/thca_tcga_pan_can_atlas_2018/data_mutations.txt`
- Bytes: 13012773
- TERT rows: 4
- OK: **✅**

### Liu 2017 JCO suppl
- URL(s): `all candidate URLs failed`
- Bytes: 0
- TERT rows: 0
- OK: **❌**

### ICGC THCA-US
- URL(s): `https://dcc.icgc.org/api/v1/projects/THCA-US, https://dcc.icgc.org/api/v1/projects/THCA-US/mutations?size=10`
- Bytes: n/a
- TERT rows: 0
- OK: **✅**

## Cross-tab with v17 DM1/DM2
```json
{
  "DM1_TERT+": 0,
  "DM1_TERT-": 0,
  "DM2_TERT+": 0,
  "DM2_TERT-": 0,
  "other_TERT+": 0,
  "other_TERT-": 513,
  "n_total_cross_tabbed": 513,
  "n_unique_tert_promoter_cases": 0,
  "n_raw_tert_rows": 4,
  "n_promoter_filtered_rows": 0
}
```

## Verdict
**SUCCESS** — cBioPortal datahub returned a TERT mutation MAF. Cross-tab written to `tert_dm_crosstab.tsv`. Next step: integrate `tert_mutations_cbio.tsv` into `sample_master_v17_tert.tsv` (column `tert_promoter`) and refresh downstream survival / pathway analyses.