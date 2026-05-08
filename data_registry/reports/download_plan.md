# Download plan — THCA data registry

**Last updated:** 2026-05-08
**Source of truth:** `data_registry/manifests/master_dataset_catalog.csv`

This plan groups datasets by access mode and download cost. **Approval gate:** any download >5 GB or any access requiring credentials needs explicit Seungho approval before execution. The `Run?` column flags this.

---

## A. Already-present locally (verify, do not re-download)

| dataset_id | local_path | size_estimate | Run? | Notes |
|---|---|---|---|---|
| P1_TCGA_THCA | `project/data/external/tcga_thca/` | ~50 GB | No | Verify with `findmnt` and dir listing |
| P1_KOREAN_K2_LEE | `/data/thca/v17_korean/` + `/data/thca/PRJEB11591_fastq/` | ~33 GB raw + processed | No | n=865 already integrated |
| P1_GSE250521 | `project/data/external/gse250521/` | ~10 GB | No | Heavily used in spatial supplement |
| P1_GSE286332 | `project/data/external/gse286332/` | ~3 GB | No | Used in P3/P4 D-day analyses |
| P1_GSE151179 | `project/data/external/gse151179/` (verify) | ~3 GB | No | RAI-dediff axis dataset |
| P1_GSE248205 | `project/data/external/gse248205/` | ~5 GB | No | Spatial supplement |
| arcasHLA on K2 | `/data/thca/_repo_offload/arcasHLA/` | 8.4 GB | No | **Paper-2 zone only** — never join to Paper-1 outcomes |
| arcasHLA on GSE213647 | `/data/thca/_repo_offload/arcasHLA_GSE213647/` | 24 GB | No | **Paper-2 zone only** |

---

## B. Public open + small (metadata-first, then matrices on approval)

For each:
1. Pull metadata via GEOparse / GEOquery into `data_registry/<zone>/<dataset_id>/metadata/`.
2. Inspect phenotype labels.
3. If labels are clean → proceed with matrix download (per-dataset approval).
4. If labels are unclear → mark `NEEDS_MANUAL_REVIEW` and stop.

**Metadata-only commands (safe to run unattended):**

```bash
# Bulk pull SOFT family files for all GSE accessions in master catalog
python3 scripts/boundary_audit/fetch_geo_metadata.py \
    --catalog data_registry/manifests/master_dataset_catalog.csv \
    --out data_registry/_geo_metadata_cache/ \
    --metadata-only
```

(NOT YET WRITTEN — script to be authored once Seungho approves the metadata pull.)

| dataset_id | source | est. matrix size | Run? | Notes |
|---|---|---|---|---|
| P1_GSE33630 | GEO | ~150 MB | Approve | ATC microarray |
| P1_GSE29265 | GEO | ~80 MB | Approve | Verify subtype |
| P1_GSE60542 | GEO | ~100 MB | Approve | Verify subtype |
| P1_GSE76039 | GEO | ~5 GB | Approve | Cross-check vs Landa 2016 |
| P1_GSE3467 | GEO | ~30 MB | Approve | Small classic |
| P1_GSE53157 | GEO | ~80 MB | Approve | Verify subtype |
| P1_GSE27155 | GEO | ~120 MB | Approve | Multiple subtypes |
| P1_GSE3678 | GEO | ~30 MB | Approve | Paired |
| P1_GSE50901 | GEO | ~100 MB | Approve | Verify |
| P1_GSE65144 | GEO | ~40 MB | Approve | ATC |
| P1_GSE3950 | GEO | ~80 MB | Approve | Verify |
| P1_GSE6004 | GEO | ~40 MB | Approve | FTC/FA |
| P1_GSE54958 | GEO | ~3 GB | Approve | Methylation 450K |
| P1_GSE51090 | GEO | ~2 GB | Approve | Methylation |
| P1_GSE53051 | GEO | ~5 GB | Approve | Multi-cancer methylation — filter to thyroid |
| P1_GSE82208 | GEO | ~2 GB | Approve | FNA cytology |
| P1_GSE104260 | GEO | ~2 GB | Approve | RAI-refractory |
| P1_GSE158285 | GEO | ~3 GB | Approve | ATC scRNA |
| P1_GSE171472 | GEO | ~5 GB | Approve | PTC scRNA |
| P1_GSE193581 | GEO | ~5 GB | Approve | PTC scRNA |
| P1_GSE197211 | GEO | ~5 GB | Approve | Verify |
| P1_GSE214083 | GEO | ~5 GB | Approve | Verify |
| P1_GSE242530 | GEO | ~10 GB | Approve | Spatial |
| P1_GSE246297 | GEO | ~3 GB | Approve | Verify exists |
| P1_GSE40807 | GEO | ~80 MB | Approve | mPTC vs PTC |
| P1_GSE58545 | GEO | ~100 MB | Approve | Verify |
| P1_GSE154763 | GEO | ~30 GB raw / ~5 GB processed | Approve | Pan-cancer myeloid |
| P1_GSE150082 | GEO | ~3 GB | Approve | Verify |
| P1_GSE161656 | GEO | ~3 GB | Approve | Verify |
| P1_GSE76096 | GEO | ~2 GB | Approve | RAI-refractory |
| P1_GSE184362 | GEO | ~10 GB | Approve | PTC ecosystem scRNA |
| P1_GSE191288 | GEO | ~10 GB | Approve | PTC scRNA |
| P1_GSE213647 | GEO | NA (arcasHLA already pulled) | No | Verify primary expression matrix separately |
| BR_GSE138198 | GEO | ~150 MB | Approve | Bridge expression |
| BR_GSE163203 | GEO | ~10 GB | Approve | Bridge scRNA |
| BR_GSE230424 | GEO | ~50 GB | **Approve required** | Bridge spatial — large |
| BR_GSE29315 | GEO | ~80 MB | Approve | Verify phenotype |
| BR_GSE6339 | GEO | ~40 MB | Approve | AITD candidate |
| BR_GSE71956 | GEO | ~5 GB | Approve | GD T-cell methylation |
| BR_GSE136709 | GEO | ~3 GB | Approve | GD B-cell |
| BR_GSE112079 | GEO | ~3 GB | Approve | HT candidate |
| BR_GSE36841 | GEO | ~80 MB | Approve | HT-only candidate |
| BR_E_MEXP_2612 | ArrayExpress | ~80 MB | Approve | AITD tissue |

**Per-dataset GEO matrix command template (run after approval):**

```bash
DSET=GSE33630   # replace
DEST=data_registry/paper1_cancer_core/${DSET}
mkdir -p "$DEST"
wget -P "$DEST" "https://ftp.ncbi.nlm.nih.gov/geo/series/${DSET:0:${#DSET}-3}nnn/${DSET}/matrix/${DSET}_series_matrix.txt.gz"
wget -P "$DEST" "https://ftp.ncbi.nlm.nih.gov/geo/series/${DSET:0:${#DSET}-3}nnn/${DSET}/soft/${DSET}_family.soft.gz"
```

For SRA-FASTQ pulls (scRNA-seq raw):

```bash
PRJ=PRJNA685185
prefetch --max-size 100G "$PRJ" -O /data/thca/sra_cache/
fasterq-dump --split-files --threads 8 -O /data/thca/_pre_processing/${PRJ}/ /data/thca/sra_cache/${PRJ}/*.sra
```

**STAR / parallel pipelines must use `--tmpdir /data/thca/_tmp` per CLAUDE.md.**

---

## C. cBioPortal datasets (API only)

```bash
# Example: TCGA-THCA Pan-Cancer Atlas re-call
python3 -m cbioportal_api download \
    --study thca_tcga_pan_can_atlas_2018 \
    --out data_registry/paper1_cancer_core/thca_tcga_pan_can_atlas_2018/
```

| dataset_id | endpoint | Run? | Notes |
|---|---|---|---|
| P1_TCGA_PANCAN_2018 | cBioPortal | Approve | Re-harmonized |
| P1_TCGA_THCA_PUB | cBioPortal | No | Already-used (TERT mutations recovered) |
| P1_MSK_LANDA_2016 | cBioPortal | Approve | Confirm overlap with Landa 2016 supplementary |

---

## D. ICGC + GENIE + DepMap + CCLE (downloadable but large)

| dataset_id | size | Run? | Command |
|---|---|---|---|
| P1_ICGC_THCA_SA | ~10–50 GB | **Approve required** | `dcc download --project THCA-SA -o data_registry/paper1_cancer_core/icgc_thca_sa/` |
| P1_GENIE_THYROID | ~2 GB extracts | **Approve + Synapse login** | `synapse get -r syn7222066 --downloadLocation data_registry/paper1_cancer_core/genie_thyroid/` |
| P1_DEPMAP_THYROID | ~5 GB curated | Approve | https://depmap.org/portal/download/all/ |
| P1_CCLE_THYROID | ~10 GB | Approve | `wget` from CCLE downloads page |

---

## E. HLA reference + AITD literature (manual extraction)

| dataset_id | mode | Run? | Notes |
|---|---|---|---|
| P2_AFND_KOR | manual export from allelefrequencies.net | Manual | Save query CSVs |
| P2_AFND_EAS | manual export | Manual | Save query CSVs |
| P2_AFND_GLOBAL | manual export | Manual | Save query CSVs |
| P2_AFND_GRAVES_REF | manual export | Manual | Disease-cohort frequencies (cases) |
| P2_AFND_HASHIMOTO_REF | manual export | Manual | Disease-cohort frequencies (cases) |
| P2_KOR_5802 / P2_PARK_2010_KORHLA | PMID lookup → supplement extraction | Manual | Confirm PMID first; resolve duplicate |
| P2_SHIN_2019 | journals.plos.org supplement | Manual | DOI 10.1371/journal.pone.0223334 |
| P2_CHO_2011 | PubMed → PMC supplement | Manual | Locate exact PMID |
| P2_PARK_2005 | PubMed → manuscript table | Manual | GD-only |
| P2_JANG_2011 | PubMed → manuscript table | Manual | DRB1 SBT |
| P2_BAEK_2021 | PubMed → supplement | Manual | Class II NGS reference |
| P2_KOGES_KWAK_2014 | PubMed + KoGES portal | Manual | GWAS sum stats only |
| P2_TOMER_GD_HLA / P2_BEDNARCZUK_GD / P2_CHU_2018_GD / P2_HU_CHINESE_AITD / P2_HEWARD_UK_GD / P2_VAIDYA_UK_GD / P2_SIMMONDS_UK_AITD | PubMed/PMC | Manual | Extract primary-data tables only — ignore review-narrative |
| P2_IPD_IMGT_HLA | release tarball | Approve | https://github.com/ANHIG/IMGTHLA |
| P2_HLA_LIGAND_ATLAS | tarball download | Approve | hla-ligand-atlas.org |
| P2_NETMHC_REF | service registration | Manual | https://services.healthtech.dtu.dk/ |
| P2_FINNGEN_THYROID / P2_FINNGEN_HASHIMOTO / P2_FINNGEN_GRAVES | summary stats download | Approve | r10.finngen.fi (release 10) |
| P2_GWAS_CATALOG / P2_GWAS_HASHIMOTO / P2_GWAS_GRAVES / P2_GWAS_TPO_AB / P2_GWAS_TSH / P2_GWAS_HYPOTHY / P2_GWAS_HYPERTHY | EBI GWAS Catalog API | Approve | Filter HLA-region SNPs separately |
| P2_IEU_OPENGWAS | API queries | Approve | gwas.mrcieu.ac.uk |
| P2_PANUKBB | summary stats download | Approve | pan.ukbb.broadinstitute.org |

---

## F. Restricted / blocked

| dataset_id | mode | Run? | Notes |
|---|---|---|---|
| RR_DBGAP_TCGA_RAW_BAMS | dbGaP controlled | **No — even if approved, never run HLA imputation against these BAMs for Paper 1** | Hard rule per HLA_CANCER_SEPARATION_RULES.md §1.2 |
| RR_EGA_KOREAN_TBD | EGA controlled | **Blocked** | Catalog only |
| RR_KOREA_BIOBANK | KoGES/KBN application | **Blocked** | IRB + KCDC application required |
| P2_BBJ_THYROID | BBJ application | **Blocked** | Verify access route |
| P2_UKBB_THYROID | UKB application | **Blocked** | Application required for primary data |
| P2_KMDP_HLA | KMDP aggregate-only | **Blocked for individual data** | Aggregate frequency tables only |
| P2_KCDC_KORHLA | KCDC | **Blocked** | Aggregate only |

---

## G. Approval workflow

1. Seungho reads this plan.
2. For each dataset, sets a status: `approve_full` / `approve_metadata_only` / `defer` / `reject`.
3. The script `scripts/boundary_audit/run_downloads.py` (TBD) will read the approval file and execute only approved actions.
4. All downloaded data lands under `project/data/external/<dataset_id>/` (bind-mounted to `/data`) per CLAUDE.md disk-layout rules.

---

## H. Risk notes

- **Disk:** Bulk downloads (ICGC THCA-SA ~50 GB, BR_GSE230424 ~50 GB, P1_GSE154763 ~30 GB raw) MUST land on `/data`, not root. Confirm with `findmnt project/data` before pulling.
- **`/tmp`:** STAR / GNU parallel pipelines must use `--tmpdir /data/thca/_tmp` (per `v17_pod_D_k2_star_tmp_failure` lesson).
- **Boundary:** No download in this plan implies any HLA-cancer joining. arcasHLA outputs at `/data/thca/_repo_offload/arcasHLA*/` are sequestered to Paper 2 zone.
