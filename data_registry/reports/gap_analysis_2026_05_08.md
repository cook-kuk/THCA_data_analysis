# Catalog gap analysis — 2026-05-08

**Master:** data_registry/manifests/master_dataset_catalog.csv (357 rows)
**Claim matrix:** data_registry/manifests/claim_matrix.csv (15 claims)

## 1. Zone-level verification rates

| Zone | Total | OK | OK-XCancer | Manual | Unverified | OK% |
|---|---:|---:|---:|---:|---:|---:|
| `BRIDGE_QUARANTINE_HT_PTC` | 27 | 6 | 0 | 5 | 16 | 22.2% |
| `PAPER1_CANCER_CORE` | 221 | 37 | 2 | 47 | 135 | 17.6% |
| `PAPER2_HLA_AUTOIMMUNE` | 40 | 1 | 0 | 17 | 22 | 2.5% |
| `REFERENCE_POPULATION_HLA` | 49 | 0 | 0 | 7 | 42 | 0.0% |
| `RESTRICTED_OR_REJECTED` | 20 | 0 | 0 | 0 | 20 | 0.0% |

## 2. Claim-by-claim support

Each claim's `allowed_datasets` list is enumerated against master verification status. **WEAK** = no VERIFIED-OK or VERIFIED-OK-XCANCER datasets back the claim.

| Claim | Paper | Strength | Total | OK | OK-XC | Manual | Unverified | Missing-from-master | Status |
|---|---|---|---:|---:|---:|---:|---:|---:|---|
| `P1_C1` | Paper 1 | STRONG | 8 | 5 | 0 | 0 | 3 | 0 | ✓ |
| `P1_C2` | Paper 1 | STRONG | 5 | 1 | 0 | 0 | 4 | 0 | ✓ |
| `P1_C3` | Paper 1 | STRONG | 9 | 4 | 0 | 0 | 2 | 3 | ✓ |
| `P1_C4` | Paper 1 | STRONG | 2 | 1 | 0 | 0 | 1 | 0 | ✓ |
| `P1_C5` | Paper 1 | MEDIUM | 4 | 3 | 0 | 0 | 0 | 1 | ✓ |
| `P1_C6` | Paper 1 | WEAK | 6 | 1 | 0 | 0 | 3 | 2 | ✓ |
| `P1_C7` | Paper 1 | MEDIUM | 4 | 2 | 0 | 0 | 1 | 1 | ✓ |
| `P2_C1` | Paper 2 | EXPLORATORY | 12 | 1 | 0 | 5 | 6 | 0 | ✓ |
| `P2_C2` | Paper 2 | STRONG | 10 | 0 | 0 | 3 | 7 | 0 | ⚠ WEAK |
| `P2_C3` | Paper 4 reserve | EXPLORATORY | 9 | 1 | 0 | 5 | 3 | 0 | ✓ |
| `P2_C4` | Paper 2 (sensitivity) | EXPLORATORY | 9 | 0 | 0 | 4 | 5 | 0 | ⚠ WEAK |
| `BRIDGE_C1` | Bridge | HYPOTHESIS-ONLY | 12 | 6 | 0 | 0 | 4 | 2 | ✓ |
| `BRIDGE_C2` | Bridge | HYPOTHESIS-ONLY | 3 | 1 | 0 | 0 | 2 | 0 | ✓ |
| `P1_C8` | Paper 1 | MEDIUM | 7 | 3 | 0 | 0 | 4 | 0 | ✓ |
| `P1_C9` | Paper 1 | MEDIUM | 5 | 1 | 0 | 0 | 2 | 2 | ✓ |

## 3. Per-claim dataset support detail

### `P1_C1` — Paper 1 DM1/DM2 resolves the BRAF/RAS-negative PTC compartment
_Paper 1 · strength=STRONG · validation=Korean K2+Lee n=865 already validates_

- ○ unverified `P1_TCGA_THCA`
- ✅ `P1_TCGA_PANCAN_2018`
- ✅ `P1_TCGA_THCA_PUB`
- ○ unverified `P1_KOREAN_K2_LEE`
- ✅ `P1_GSE286332`
- ✅ `P1_GSE33630`
- ✅ `P1_GSE76039`
- ○ unverified `P1_MSK_LANDA_2016`

### `P1_C2` — Paper 1 DM1 enriches fusion-positive tumors
_Paper 1 · strength=STRONG · validation=Cross-cohort fusion-status overlay required for full validation_

- ○ unverified `P1_TCGA_THCA`
- ✅ `P1_TCGA_PANCAN_2018`
- ○ unverified `P1_KOREAN_K2_LEE`
- ○ unverified `P1_GENIE_THYROID`
- ○ unverified `P1_MSK_LANDA_2016`

### `P1_C3` — Paper 1 DM1 has thyroid-differentiation / RAI-lineage silencing
_Paper 1 · strength=STRONG · validation=RAI-refractory cohorts independently confirm_

- ○ unverified `P1_TCGA_THCA`
- ✅ `P1_TCGA_PANCAN_2018`
- ○ unverified `P1_KOREAN_K2_LEE`
- ✅ `P1_GSE151179`
- ❌ missing-from-master `P1_GSE104260`
- ❌ missing-from-master `P1_GSE76096`
- ✅ `P1_GSE54958`
- ✅ `P1_GSE51090`
- ❌ missing-from-master `P1_GSE53051`

### `P1_C4` — Paper 1 DM1 sub-B (TCGA n=56) is NBNR-equivalent and partially Hashimoto-like
_Paper 1 · strength=STRONG · validation=Validated within Paper-1 zone via expression-only; HLA strictly excluded_

- ○ unverified `P1_TCGA_THCA`
- ✅ `P1_GSE286332`

### `P1_C5` — DM1 spatial niche organization replicates across cohorts
_Paper 1 · strength=MEDIUM · validation=Cross-platform spatial validation pending_

- ✅ `P1_GSE250521`
- ✅ `P1_GSE248205`
- ❌ missing-from-master `P1_GSE242530`
- ✅ `P1_GSE184362`

### `P1_C6` — RAI re-induction is feasible in DM1 sub-B (translational candidate)
_Paper 1 · strength=WEAK · validation=Prospective trial required_

- ○ unverified `P1_TCGA_THCA`
- ✅ `P1_GSE151179`
- ❌ missing-from-master `P1_GSE104260`
- ❌ missing-from-master `P1_GSE76096`
- ○ unverified `P1_DEPMAP_THYROID`
- ○ unverified `P1_CCLE_THYROID`

### `P1_C7` — DM1 promoter methylation is conserved across external cohorts
_Paper 1 · strength=MEDIUM · validation=Cross-platform methylation validation pending_

- ○ unverified `P1_TCGA_THCA`
- ✅ `P1_GSE54958`
- ✅ `P1_GSE51090`
- ❌ missing-from-master `P1_GSE53051`

### `P2_C1` — Certain HLA alleles are candidates for future validation in AITD/HT-overlap context
_Paper 2 · strength=EXPLORATORY · validation=Prospective Korean HT-overlap PTC HLA cohort required_

- ⚪ manual `P2_SHIN_2019`
- ⚪ manual `P2_CHO_2011`
- ⚪ manual `P2_PARK_2005`
- ✅ `P2_JANG_2011`
- ⚪ manual `P2_BAEK_2021`
- ○ unverified `P2_AFND_KOR`
- ○ unverified `P2_AFND_EAS`
- ○ unverified `P2_KMDP_HLA`
- ⚪ manual `P2_KOR_5802`
- ○ unverified `P2_GWAS_HASHIMOTO`
- ○ unverified `P2_GWAS_GRAVES`
- ○ unverified `P2_GWAS_TPO_AB`

### `P2_C2` — Korean population HLA baseline differs from other ancestry references
_Paper 2 · strength=STRONG · validation=Distinguish allele frequency vs carrier frequency_

⚠ **WEAK** — no VERIFIED-OK datasets currently back this claim. All supporting datasets are unverified or only manual-lookup-needed.

- ○ unverified `P2_AFND_KOR`
- ○ unverified `P2_AFND_EAS`
- ○ unverified `P2_AFND_GLOBAL`
- ⚪ manual `P2_KOR_5802`
- ⚪ manual `P2_PARK_2010_KORHLA`
- ⚪ manual `P2_BAEK_2021`
- ○ unverified `P2_KMDP_HLA`
- ○ unverified `P2_HLA_NET_KOR`
- ○ unverified `P2_IPD_IMGT_HLA`
- ○ unverified `P2_PANUKBB`

### `P2_C3` — GD-specific HLA signal is distinct from HT signal in Korean cohorts
_Paper 4 reserve · strength=EXPLORATORY · validation=Paper 4 reserve only_

- ⚪ manual `P2_PARK_2005`
- ✅ `P2_JANG_2011`
- ⚪ manual `P2_CHU_2018_GD`
- ⚪ manual `P2_HEWARD_UK_GD`
- ⚪ manual `P2_VAIDYA_UK_GD`
- ⚪ manual `P2_TOMER_GD_HLA`
- ○ unverified `P2_GWAS_GRAVES`
- ○ unverified `P2_FINNGEN_GRAVES`
- ○ unverified `P2_AFND_GRAVES_REF`

### `P2_C4` — East-Asian and non-Korean replication context strengthens AITD HLA prior
_Paper 2 (sensitivity) · strength=EXPLORATORY · validation=Ancestry caveat mandatory_

⚠ **WEAK** — no VERIFIED-OK datasets currently back this claim. All supporting datasets are unverified or only manual-lookup-needed.

- ○ unverified `P2_AFND_EAS`
- ⚪ manual `P2_HU_CHINESE_AITD`
- ⚪ manual `P2_CHU_2018_GD`
- ○ unverified `P2_BBJ_THYROID`
- ○ unverified `P2_FINNGEN_HASHIMOTO`
- ○ unverified `P2_FINNGEN_GRAVES`
- ○ unverified `P2_PANUKBB`
- ⚪ manual `P2_SIMMONDS_UK_AITD`
- ⚪ manual `P2_BEDNARCZUK_GD`

### `BRIDGE_C1` — HT+PTC expression datasets show immune-context overlap
_Bridge · strength=HYPOTHESIS-ONLY · validation=Used for hypothesis generation only_

- ✅ `BR_GSE138198`
- ○ unverified `BR_GSE163203`
- ○ unverified `BR_GSE230424`
- ✅ `BR_GSE29315`
- ✅ `BR_GSE6339`
- ✅ `BR_GSE71956`
- ✅ `BR_GSE136709`
- ○ unverified `BR_E_MEXP_2612`
- ○ unverified `BR_GSE54958_HTPTC`
- ❌ missing-from-master `BR_GSE112079`
- ❌ missing-from-master `BR_GSE36841`
- ✅ `P1_GSE286332`

### `BRIDGE_C2` — B-cell / TLS / AICDA enrichment in PTC+HT supports antigen-driven response
_Bridge · strength=HYPOTHESIS-ONLY · validation=No validation required at bridge level_

- ✅ `P1_GSE286332`
- ○ unverified `BR_GSE163203`
- ○ unverified `BR_GSE230424`

### `P1_C8` — DM1/DM2 epigenetic re-programming is mechanistically distinct from BRAF-driven oncogenesis
_Paper 1 · strength=MEDIUM · validation=DepMap perturbation evidence pending_

- ○ unverified `P1_TCGA_THCA`
- ○ unverified `P1_KOREAN_K2_LEE`
- ✅ `P1_GSE54958`
- ✅ `P1_GSE51090`
- ○ unverified `P1_DEPMAP_THYROID`
- ○ unverified `P1_CCLE_THYROID`
- ✅ `P1_GSE286332`

### `P1_C9` — DM1 prognostic stratification within BRAF-RAS-negative tumors
_Paper 1 · strength=MEDIUM · validation=External Korean OS/DSS pending_

- ○ unverified `P1_TCGA_THCA`
- ✅ `P1_TCGA_PANCAN_2018`
- ○ unverified `P1_KOREAN_K2_LEE`
- ❌ missing-from-master `P1_GSE104260`
- ❌ missing-from-master `P1_GSE76096`

## 4. ⚠ WEAKLY SUPPORTED CLAIMS (priority fill list)

These claims need at least one VERIFIED-OK supporting dataset before they should ride into Paper 1/2/3/4 prose.

- `P2_C2` (Paper 2, strength=STRONG): Korean population HLA baseline differs from other ancestry references
- `P2_C4` (Paper 2 (sensitivity), strength=EXPLORATORY): East-Asian and non-Korean replication context strengthens AITD HLA prior

## 5. Top MANUAL-NEEDED datasets by claim-citation count

These manual-needed entries underpin multiple claims; resolving their PMIDs / citations would unlock the most claim-support.

| Dataset | Cited by # claims |
|---|---:|
| `P2_PARK_2005` | 2 |
| `P2_BAEK_2021` | 2 |
| `P2_KOR_5802` | 2 |
| `P2_CHU_2018_GD` | 2 |
| `P2_SHIN_2019` | 1 |
| `P2_CHO_2011` | 1 |
| `P2_PARK_2010_KORHLA` | 1 |
| `P2_HEWARD_UK_GD` | 1 |
| `P2_VAIDYA_UK_GD` | 1 |
| `P2_TOMER_GD_HLA` | 1 |
| `P2_HU_CHINESE_AITD` | 1 |
| `P2_SIMMONDS_UK_AITD` | 1 |
| `P2_BEDNARCZUK_GD` | 1 |
