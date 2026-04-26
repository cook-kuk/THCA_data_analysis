#!/usr/bin/env python3
"""v8.1 — CCLE thyroid cell-line BRS52 validation.

Pulls expression for the 23 CCLE thyroid cell lines via cBioPortal REST API
(study ccle_broad_2019, profile ccle_broad_2019_rna_seq_mrna RPKM), applies
the BRS52 expression signature centroids built on TCGA-THCA, and computes
sign concordance with the v8.1 raw-count DESeq2 log2FC.

This is a fundamentally different substrate (cancer cell lines, in vitro,
RNA-seq RPKM) from the patient-tumour cohorts already validated, so a
positive concordance is meaningful additional support.

Outputs (results/v8p1_rigor/f_ccle_validation/):
  ccle_thyroid_metadata.tsv       — 23 lines × histology / oncotree / etc.
  ccle_thyroid_expression.tsv     — gene × line expression
  ccle_brs_labels.tsv             — per-line BRS score + label
  ccle_concordance.tsv            — sign concordance vs TCGA raw-DESeq2
  ccle_validation_analysis.md     — narrative

Log: logs/v8p1_ccle_brs.log
"""
from __future__ import annotations
import sys, time, json
from pathlib import Path
import numpy as np
import pandas as pd
import requests

ROOT = Path("/opt/thyroid-dash/project")
RES  = ROOT / "results" / "v8p1_rigor"
OUT  = RES / "f_ccle_validation"
OUT.mkdir(parents=True, exist_ok=True)
LOG  = ROOT / "logs" / "v8p1_ccle_brs.log"

API = "https://www.cbioportal.org/api"
STUDY = "ccle_broad_2019"
PROFILE = "ccle_broad_2019_rna_seq_mrna"   # RPKM

_log_f = open(LOG, "w")
def log(msg: str) -> None:
    stamp = time.strftime("[%Y-%m-%d %H:%M:%S]")
    print(f"{stamp} {msg}", flush=True)
    _log_f.write(f"{stamp} {msg}\n"); _log_f.flush()

log("v8.1 CCLE thyroid BRS52 validation")
log("=" * 60)

# --------------------------------------------------------------------
# 1. Pull thyroid sample metadata
# --------------------------------------------------------------------
r = requests.get(f"{API}/studies/{STUDY}/clinical-data",
                 params={"clinicalDataType":"SAMPLE","projection":"SUMMARY"},
                 timeout=60)
r.raise_for_status()
clin = r.json()
log(f"[clin] {len(clin)} clinical-data rows")

from collections import defaultdict
sample_meta = defaultdict(dict)
for row in clin:
    sample_meta[row["sampleId"]][row["clinicalAttributeId"]] = row["value"]

thyroid_lines = [s for s, a in sample_meta.items()
                 if "thyroid" in a.get("CANCER_TYPE","").lower()]
log(f"[thyroid] CCLE thyroid lines: {len(thyroid_lines)}")

meta_rows = []
for s in thyroid_lines:
    a = sample_meta[s]
    meta_rows.append(dict(
        sample_id=s,
        name=a.get("NAME","?"),
        cancer_type_detailed=a.get("CANCER_TYPE_DETAILED","?"),
        oncotree=a.get("ONCOTREE_CODE","?"),
        hist_subtype=a.get("HIST_SUBTYPE1","?"),
        depmap_id=a.get("DEPMAPID","?"),
        primary_site=a.get("PRIMARY_SITE","?"),
        sample_type=a.get("SAMPLE_TYPE","?"),
    ))
meta_df = pd.DataFrame(meta_rows)
meta_df.to_csv(OUT / "ccle_thyroid_metadata.tsv", sep="\t", index=False)
log(f"[meta] wrote {OUT / 'ccle_thyroid_metadata.tsv'}")
log(f"[meta] histology breakdown:\n{meta_df['cancer_type_detailed'].value_counts().to_string()}")

# --------------------------------------------------------------------
# 2. Load BRS52 gene list, fetch expression for those genes only
# --------------------------------------------------------------------
brs = pd.read_csv(ROOT/"results"/"v12_literature"/"crosscheck"/
                  "reference_gene_lists"/"Chakravarty_2011_BRS52.tsv",
                  sep="\t")["gene"].tolist()
log(f"[brs] loaded {len(brs)} BRS52 genes")

# Resolve gene symbols → Entrez IDs
r = requests.post(f"{API}/genes/fetch", params={"geneIdType":"HUGO_GENE_SYMBOL"},
                  json=brs, timeout=60)
r.raise_for_status()
gene_records = r.json()
log(f"[genes] resolved {len(gene_records)}/{len(brs)} BRS genes via cBioPortal")

# Also pull the v8.1 raw-DESeq2 sig genes — we need their expression in CCLE
# to compute the concordance test. To avoid pulling 6 605 individual genes,
# we'll just pull the BRS+all-sig union via a sample-list query.
de = pd.read_csv(RES/"a_deseq2_raw"/"deseq2_tcga_raw.tsv", sep="\t")
de = de.dropna(subset=["gene_name"]).drop_duplicates("gene_name")
de_sig_genes = de.loc[de["significant"], "gene_name"].tolist()
log(f"[de] raw-DESeq2 sig genes: {len(de_sig_genes)}")

target_genes = sorted(set(brs) | set(de_sig_genes))
log(f"[fetch] union of BRS + sig genes: {len(target_genes)}")

# Resolve all target genes to entrez IDs
r = requests.post(f"{API}/genes/fetch", params={"geneIdType":"HUGO_GENE_SYMBOL"},
                  json=target_genes, timeout=120)
r.raise_for_status()
gene_lookup = {g["hugoGeneSymbol"]: g for g in r.json()}
log(f"[fetch] resolved {len(gene_lookup)}/{len(target_genes)} target genes")

entrez_ids = [gene_lookup[g]["entrezGeneId"] for g in target_genes if g in gene_lookup]

# --------------------------------------------------------------------
# 3. Build a sample list and pull molecular data
# --------------------------------------------------------------------
# cBioPortal /molecular-data takes a sampleListId or sampleIds list.
# Use sampleIds explicitly.

# Endpoint: POST /molecular-profiles/{profileId}/molecular-data/fetch
chunk = 1500   # entrez IDs per request
all_data = []
log(f"[expr] pulling RPKM for {len(entrez_ids)} genes × {len(thyroid_lines)} lines, "
    f"chunked at {chunk}/req")
for i in range(0, len(entrez_ids), chunk):
    sub = entrez_ids[i:i+chunk]
    r = requests.post(
        f"{API}/molecular-profiles/{PROFILE}/molecular-data/fetch",
        json={"entrezGeneIds": sub, "sampleIds": thyroid_lines},
        timeout=180,
    )
    if r.status_code != 200:
        log(f"[expr] chunk {i//chunk} HTTP {r.status_code} body[:200]: {r.text[:200]}")
        continue
    all_data.extend(r.json())
    log(f"[expr]  chunk {i//chunk+1}: +{len(r.json())} rows  (total {len(all_data)})")

if not all_data:
    log("[expr] FAILED to fetch any data — aborting")
    sys.exit(1)

# Build a gene × sample matrix
df = pd.DataFrame(all_data)
log(f"[expr] dataframe cols: {list(df.columns)[:8]}")

# CCLE molecular-data rows have: entrezGeneId, sampleId, value (RPKM)
expr = df.pivot_table(index="entrezGeneId", columns="sampleId", values="value")
log(f"[expr] pivot: {expr.shape}")

# Map entrez → HGNC for downstream
entrez2hgnc = {g["entrezGeneId"]: g["hugoGeneSymbol"]
               for g in gene_lookup.values()}
expr.index = expr.index.map(entrez2hgnc)
expr = expr[~expr.index.duplicated(keep="first")]
expr.to_csv(OUT / "ccle_thyroid_expression.tsv", sep="\t")
log(f"[expr] wrote {OUT / 'ccle_thyroid_expression.tsv'}  shape {expr.shape}")

# RPKM → log2(RPKM + 1) for centroid-correlation comparability with TCGA log2-TPM
expr_log = np.log2(expr.fillna(0) + 1.0)

# --------------------------------------------------------------------
# 4. Build BRS centroids on TCGA-THCA (mirror BRS validation script)
# --------------------------------------------------------------------
TCGA_TPM = Path("/data/thca/v5_cross_cancer/raw/THCA/tcga_log2tpm.tsv.gz")
TCGA_LBL = Path("/data/thca/v5_cross_cancer/raw/THCA/tcga_samples_labels.tsv")

tpm = pd.read_csv(TCGA_TPM, sep="\t", index_col=0)
labels = pd.read_csv(TCGA_LBL, sep="\t").set_index("sample_id")["label"]

brs_in_tcga = [g for g in brs if g in tpm.index]
samples_in = [s for s in tpm.columns if s in labels.index]
tpm_brs = tpm.loc[brs_in_tcga, samples_in]
y = labels.loc[samples_in]
c_braf = tpm_brs.loc[:, (y=="BRAF").values].mean(axis=1)
c_ras  = tpm_brs.loc[:, (y=="RAS").values].mean(axis=1)
log(f"[centroid] BRAF n={(y=='BRAF').sum()}  RAS n={(y=='RAS').sum()}  "
    f"BRS-genes-in-TCGA={len(brs_in_tcga)}")

# --------------------------------------------------------------------
# 5. Apply BRS to CCLE thyroid lines
# --------------------------------------------------------------------
common = expr_log.index.intersection(c_braf.index)
log(f"[brs-ccle] BRS genes overlapping CCLE: {len(common)}")
expr_brs = expr_log.loc[common]
cb = c_braf.loc[common]
cr = c_ras.loc[common]

Xc = expr_brs.values - expr_brs.values.mean(axis=0, keepdims=True)
nx = np.sqrt((Xc**2).sum(axis=0))
nb = np.sqrt(((cb - cb.mean())**2).sum())
nr = np.sqrt(((cr - cr.mean())**2).sum())
corr_b = (Xc * (cb - cb.mean()).values[:, None]).sum(axis=0) / (nx*nb + 1e-12)
corr_r = (Xc * (cr - cr.mean()).values[:, None]).sum(axis=0) / (nx*nr + 1e-12)
score = corr_b - corr_r
brs_label = np.where(score > 0, "BRAF_like", "RAS_like")
brs_df = pd.DataFrame({
    "sample": expr_brs.columns,
    "corr_braf": corr_b,
    "corr_ras": corr_r,
    "brs_score": score,
    "brs_label": brs_label,
}).set_index("sample")
brs_df = brs_df.join(meta_df.set_index("sample_id")[
    ["name","cancer_type_detailed","oncotree","hist_subtype"]
], how="left")
brs_df.to_csv(OUT / "ccle_brs_labels.tsv", sep="\t")
log(f"[brs-ccle] BRS labels:\n{brs_df['brs_label'].value_counts().to_string()}")
log(f"[brs-ccle] By histology:\n"
    f"{pd.crosstab(brs_df['cancer_type_detailed'], brs_df['brs_label']).to_string()}")

# --------------------------------------------------------------------
# 6. Concordance test vs TCGA raw-DESeq2
# --------------------------------------------------------------------
de_sig = de.loc[de["significant"]].drop_duplicates("gene_name").set_index("gene_name")
log(f"[concord] raw-DESeq2 sig: {len(de_sig)}")

braf_lines = brs_df[brs_df["brs_label"]=="BRAF_like"].index
ras_lines  = brs_df[brs_df["brs_label"]=="RAS_like"].index
log(f"[concord] BRAF-like lines: {len(braf_lines)}, RAS-like: {len(ras_lines)}")

if len(braf_lines) >= 2 and len(ras_lines) >= 2:
    delta = (expr_log.loc[:, braf_lines].mean(axis=1) -
             expr_log.loc[:, ras_lines].mean(axis=1))
    overlap = de_sig.index.intersection(delta.index)
    log(f"[concord] sig × CCLE expr overlap: {len(overlap)}")
    if len(overlap):
        tcga_sgn = (de_sig.loc[overlap, "log2FoldChange"] > 0).values
        ccle_sgn = (delta.loc[overlap] > 0).values
        agree = int((tcga_sgn == ccle_sgn).sum())
        pct = 100*agree/len(overlap)
        log(f"[concord] BRAF-up sign concordance (TCGA vs CCLE thyroid): "
            f"{agree}/{len(overlap)} = {pct:.1f}%")
    else:
        agree = 0; pct = 0.0
else:
    delta = pd.Series(dtype=float)
    overlap = pd.Index([])
    agree = 0; pct = 0.0
    log("[concord] insufficient lines per BRS class for delta")

# Cross with the existing 3-cohort consensus to make a 4-cohort match-count
con3 = pd.read_csv(RES/"e_brs_validation"/"three_cohort_per_gene_consensus.tsv",
                   sep="\t").set_index("gene")
log(f"[concord] 3-cohort consensus rows: {len(con3)}")

if len(overlap):
    overlap_4 = con3.index.intersection(overlap)
    log(f"[concord] genes in all 4 cohorts: {len(overlap_4)}")
    if len(overlap_4):
        # Per-gene match count (3 GEO + CCLE) vs TCGA
        ccle_match = (((de_sig.loc[overlap_4, "log2FoldChange"] > 0).values ==
                       (delta.loc[overlap_4] > 0).values).astype(int))
        match4 = con3.loc[overlap_4, "match_count_3cohort"].values + ccle_match
        from collections import Counter
        cnt = Counter(match4.tolist())
        log(f"[concord] 4-cohort match-count distribution: "
            f"{dict(sorted(cnt.items()))}")
        log(f"[concord] genes matching TCGA in all 4 cohorts: "
            f"{(match4==4).sum()} ({100*(match4==4).sum()/len(match4):.1f}%)")
        log(f"[concord] genes matching in ≥ 3 of 4: "
            f"{(match4>=3).sum()} ({100*(match4>=3).sum()/len(match4):.1f}%)")
        out4 = pd.DataFrame({
            "gene": overlap_4,
            "match_count_4cohort": match4,
        })
        out4.to_csv(OUT / "four_cohort_per_gene_consensus.tsv", sep="\t", index=False)

# Final concordance summary
summary = pd.DataFrame({
    "metric": [
        "ccle_thyroid_lines_total",
        "ccle_brs_BRAF_like",
        "ccle_brs_RAS_like",
        "deseq2_sig_in_ccle_universe",
        "BRAF_up_concordance_TCGA_vs_CCLE",
        "ccle_brs_genes_used",
    ],
    "value": [
        len(thyroid_lines),
        int((brs_df["brs_label"]=="BRAF_like").sum()),
        int((brs_df["brs_label"]=="RAS_like").sum()),
        len(overlap),
        f"{pct:.2f}%" if len(overlap) else "n/a",
        len(common),
    ],
})
summary.to_csv(OUT / "ccle_concordance.tsv", sep="\t", index=False)
log(f"[write] {OUT / 'ccle_concordance.tsv'}")

log("v8.1 CCLE BRS validation complete")
_log_f.close()
