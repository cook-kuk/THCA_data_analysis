# v14 strengthening — External cohort scan (Tier-1)

Author run: 2026-04-25
Working dir: `/home/seungho/personal/THCA_data_analysis/project`
Output TSV: `results/v14_strengthening/external_cohort_tacstd2_braf.tsv`
Output JSON: `results/v14_strengthening/external_cohort_tacstd2_braf_summary.json`
Script: `notebooks_or_scripts/v14_strengthening_external_cohort.py`

The aim of this pass was to test whether the v14 BRAF V600E ↔ TACSTD2 (TROP2)
high finding replicates in a third / fourth independent thyroid cohort beyond
TCGA-THCA + the GEO sets already used by v8 / v9 / v13 / v14
(GSE27155, GSE29265, GSE33630, GSE60542).

---

## 1. ICGC-THCA — status

The historic ICGC Data Portal (`dcc.icgc.org`) was retired in **June 2024**.
All requests to legacy endpoints (project listing, `/api/v1/projects/THCA-CN`,
`/api/v1/projects/THCA-SA`) now 200-redirect to the static retirement notice
HTML at `https://dcc.icgc.org/index.html`. There is therefore no live REST or
GraphQL API for ICGC-THCA. The replacement is **ICGC ARGO** + the legacy
**ICGC 25K** open S3 bucket.

The open data bucket is reachable without DACO at:
`https://object.genomeinformatics.org/icgc25k-open/release_28/`. I enumerated
all THCA-prefixed projects and their per-donor data types via S3
`?list-type=2&prefix=...`:

| Project | Donors found | Open data types present (per donor) | Independent of TCGA? |
|---|---|---|---|
| **THCA-CN** Thyroid Cancer — China | 50 | `donor`, `sample`, `specimen`, **`ssm_open`** | yes |
| **THCA-SA** Thyroid Cancer — Saudi Arabia | 243 | `donor`, `sample`, `specimen`, **`ssm_open`** | yes |
| **THCA-US** TCGA Head&Neck Thyroid Carcinoma | 506 | `donor`, `sample`, `specimen`, `ssm_open`, `cnsm`, **`exp_seq`** (91 donors), `meth_array`, `mirna_seq`, `pexp` | **no — same as TCGA-THCA already in v14** |

**Verdict for ICGC-THCA**:
- THCA-CN and THCA-SA are open-access but contain **only somatic mutation
  calls** (`ssm_open`). There is no `exp_seq` / `exp_array` / `pexp` per donor
  in release 28 for these two projects. Without expression we cannot test the
  TACSTD2 differential; only BRAF V600E prevalence is recoverable.
- THCA-US carries open-access RNA-seq for 91 donors, but it is derived from
  the same TCGA-THCA cohort that is already the discovery substrate for the
  v14 paper. Re-using it is not an independent replication.
- Raw RNA-seq for THCA-CN / THCA-SA is hosted on EGA `EGAC00001000010` and
  is **DACO-controlled**. Access requires institutional DACO approval,
  which is not feasible inside this 30 to 45 min budget.

Net: ICGC-THCA does not provide an independent open expression cohort for
this analysis.

---

## 2. CPTAC-THCA — status

The Proteomic Data Commons (PDC) GraphQL endpoint at
`https://pdc.cancer.gov/graphql` was queried with `allPrograms { projects {
studies } }` and the response (HTTP 200, ~38 KB) was scanned for any study
whose `submitter_id_name`, `disease_type`, or `primary_site` contains "THCA"
or "thyroid" (case-insensitive). The full CPTAC catalogue currently exposes
the following primary sites only:

> AML, BRCA, CCRCC, COAD, GBM, HNSCC, IPMN, LSCC, LUAD, non-ccRCC, OV, PDA,
> STAD, T-ALL, UCEC, plus methodology / cell-line / xenograft / serum sub-
> studies. **No thyroid / THCA study exists** in CPTAC2, CPTAC3, CPTAC-KF,
> CPTAC-Other, PTRC, APOLLO, ICPC, or the Georgetown Proteomics Research
> Program.

A web search for a 2024-2026 CPTAC thyroid release also returned nothing.

Net: **CPTAC-THCA does not exist** as of release 2026-04. No proteomic
TROP2 + BRAF table can be pulled. This is a hard "data not available"
result — not a permission barrier.

---

## 3. Alternative public GEO cohort — GSE58545

Because both putative external cohorts proved unusable, I selected an
additional, previously-unused public GEO dataset for the strengthening pass:

- **GSE58545** — Bauer et al. (2014), *"Gene expression signature associated
  with BRAFV600E mutation in human papillary thyroid carcinoma based on
  transgenic mouse model"*.
- Platform: GPL96 / Affymetrix HG-U133A (RMA-normalised log2 expression
  embedded in the GEO series_matrix file).
- n = **45**: 27 PTC tumours with explicit per-sample BRAF/RET/RAS status in
  `Sample_characteristics_ch1` (BRAF V600E n=18, RET fusion n=8, RAS mutant
  n=1) plus 18 normal thyroids.
- TACSTD2 (TROP2) is represented by probe **202286_s_at**.
- BRAF status method: extracted directly from the GEO sample-level
  characteristics field
  `"braf/ret/ras status; braf(+)- v600e mutation presence; ret(+)- ret
  rearrangement presence, ras(+) - ras mutation presence: <CALL>"`. No
  inference from the BRS / MAPK signature was used, so the test is not
  circular with v14's expression-based scoring.
- Novelty check: `grep -rE "GSE58545"` over `results/` and `reports/` returns
  no prior hits — this dataset has not been touched by v8 / v8p1 / v13 / v14.

### 3.1 Test, identical to v14

Following v14, BRAF-WT PTC was defined as the union of RET-fusion + RAS-mut
+ any genotype-other PTC (here that is the 8 RET + 1 RAS = 9 samples).
Z-scoring was performed across PTC tumours only (normals excluded from the z
distribution) so that the cohort-internal z directly mirrors v14's
tumour-level claim.

| Group | n | mean log2 TACSTD2 (RMA) | mean PTC-internal z |
|---|---|---|---|
| PTC, BRAF V600E | 18 | 10.978 | **+0.326** |
| PTC, BRAF-WT (RET + RAS) | 9 | 9.351 | **−0.652** |
| Normal thyroid | 18 | 5.055 | n/a (excluded from z) |

- Welch t = 1.93, df ≈ 8.8, two-sided **p = 0.086**
- Mann-Whitney U = 21, two-sided **p = 0.0020**
- BRAF V600E vs Normal thyroid: Welch t = 20.7, **p ≈ 1.6e-20**
- Direction predicted by v14 (BRAF V600E > BRAF-WT): **replicated** (effect
  size +0.98 z, ~3.2-fold linear difference at log2 scale)

### 3.2 Caveats

- Small WT arm (n=9) — Welch is borderline; Mann-Whitney is
  significant but is influenced by the BRAF V600E group's monotonic high
  expression.
- HG-U133A is a 22k-probe array. RMA is intra-platform comparable but the
  effect size in z units is not directly comparable to TCGA-THCA RNA-seq z
  (different dynamic range). Direction and rank-test significance are the
  load-bearing comparators.
- One PTC tumour (NIS135) is a clear BRAF V600E low-TACSTD2 outlier
  (log2 = 8.68, z = −1.06 within PTC), pulling the BRAF mean down. The
  Mann-Whitney is robust to it; the Welch is sensitive.

---

## 4. Verdict (Korean)

> **v14 BRAF↔TACSTD2 link이 외부 4번째 코호트에서 replicate되나?**
>
> **부분 (partial replication).**
>
> - ICGC-THCA-CN/SA: 발현 데이터 없음 (`ssm_open` 변이 콜만 공개) → 검증 불가.
> - ICGC-THCA-US: TCGA-THCA와 동일 코호트라 독립 검증 아님.
> - CPTAC-THCA: PDC에 갑상선 (THCA / thyroid) 연구 자체가 없음 → 검증 불가.
> - 대체 코호트 GSE58545 (n=27 PTC, HG-U133A, BRAF/RET/RAS metadata 명시):
>   방향성 일치 (BRAF V600E PTC가 BRAF-WT PTC 대비 TACSTD2 평균 +0.98 z 높음).
>   Welch t-test p=0.086 (경계), Mann-Whitney p=0.0020 (유의). BRAF vs 정상
>   조직 p≈1.6e-20.
>
> 따라서 v14의 종양 수준 BRAF V600E ↔ TROP2-high 주장은 GSE58545에서
> **방향성 + 비모수적으로 유의하게 재현**되며, 모수적 검정에서는 작은 WT
> 표본 때문에 borderline이다. ICGC와 CPTAC는 본 분석에서 검증 자원으로
> 활용 불가능함이 확인되었으므로, 추후 strengthening은 GSE213647 (632
> 샘플 RNA-seq, BRAF는 BRS 추정 필요) 또는 GSE76039 (n=37 PDTC/ATC,
> BRAF는 supplementary table 외부 매칭 필요)로 확장하는 것을 권장.

---

## Reproducibility

```bash
# 1. fetch series matrix (no auth needed)
curl -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58545/matrix/GSE58545_series_matrix.txt.gz" \
     -o /tmp/gse58545/GSE58545_series_matrix.txt.gz
# 2. run analysis
python3 notebooks_or_scripts/v14_strengthening_external_cohort.py
```

Outputs:
- `results/v14_strengthening/external_cohort_tacstd2_braf.tsv` — per-sample
  status, log2 TACSTD2, PTC-internal z, all-sample z, plus an embedded
  `# SUMMARY` block.
- `results/v14_strengthening/external_cohort_tacstd2_braf_summary.json` —
  same summary as machine-readable JSON.
