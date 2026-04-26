#!/usr/bin/env python3
"""v6 F12 — build a single AnnData from GSE193581 (Lu et al. JCI 2023).

Selects only patients with explicit BRAF/RAS mutation calls in
JCI169653 Supplemental Table S1. Each per-GSM file is a dense
gene x cell tab-separated matrix (rows=genes, first row=cell barcodes,
first col=gene symbol).

Output: /data/thca/scrna/raw/scrna_F12.h5ad with obs columns
  patient_id, gsm, sample, cancer_type, braf_status, ras_status,
  braf_vs_ras (BRAF / RAS / NA).
"""
from __future__ import annotations

import gzip
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix

RAW_DIR = Path("/data/thca/scrna/raw/gse193581")
OUT_H5AD = Path("/data/thca/scrna/raw/scrna_F12.h5ad")
OUT_META = Path("/data/thca/scrna/raw/metadata_F12.tsv")

# (gsm, sample_tag, patient_id, cancer_type, BRAF_status, RAS_status)
# from JCI169653 Supplemental Table S1.
SAMPLES = [
    # BRAF-mutant cohort
    ("GSM5814576", "PTC03", "PTC03", "PTC", "mut", "wt"),
    ("GSM5814579", "PTC05", "PTC05", "PTC", "mut", "wt"),
    ("GSM5814580", "PTC06", "PTC06", "PTC", "mut", "wt"),
    ("GSM5814583", "ATC08", "ATC08", "ATC", "mut", "wt"),
    ("GSM5814584", "ATC09", "ATC09", "ATC", "mut", "wt"),
    # RAS-mutant cohort
    ("GSM5814586", "ATC11", "ATC11", "ATC", "wt", "mut"),
    ("GSM5814587", "ATC12", "ATC12", "ATC", "wt", "mut"),
    ("GSM5814588", "ATC13", "ATC13", "ATC", "wt", "mut"),
    ("GSM5814591", "ATC17", "ATC17", "ATC", "wt", "mut"),
]


def read_one(gsm: str, tag: str) -> ad.AnnData:
    p = RAW_DIR / f"{gsm}_{tag}_UMI.txt.gz"
    print(f"[read] {gsm} {tag} from {p.name}", flush=True)
    # Robust read: header row = barcodes, first col = gene symbol.
    df = pd.read_csv(p, sep="\t", index_col=0, compression="gzip")
    # df: index=gene, columns=barcodes; transpose so rows=cells, cols=genes
    Xc = csr_matrix(df.values.T.astype(np.float32))
    var = pd.DataFrame({"gene_symbol": df.index.astype(str).values}, index=df.index.astype(str))
    obs_idx = [str(b) for b in df.columns]
    obs = pd.DataFrame(index=obs_idx)
    a = ad.AnnData(X=Xc, obs=obs, var=var)
    return a


def main():
    t0 = time.time()
    parts = []
    for gsm, tag, pat, cancer, braf, ras in SAMPLES:
        a = read_one(gsm, tag)
        a.obs["patient_id"] = pat
        a.obs["gsm"] = gsm
        a.obs["sample"] = tag
        a.obs["cancer_type"] = cancer
        a.obs["braf_status"] = braf
        a.obs["ras_status"] = ras
        a.obs["braf_vs_ras"] = "BRAF" if braf == "mut" else ("RAS" if ras == "mut" else "NA")
        parts.append(a)
        print(f"  {tag}: {a.n_obs} cells x {a.n_vars} genes (BRAF={braf}, RAS={ras})", flush=True)

    print(f"[concat] {len(parts)} samples", flush=True)
    full = ad.concat(parts, axis=0, join="inner", merge="first", index_unique="-")
    full.var_names_make_unique()
    full.obs_names_make_unique()

    full.obs["mutation_status"] = full.obs["braf_vs_ras"]

    print(f"[done] {full.n_obs} cells x {full.n_vars} genes", flush=True)
    print(f"  patients: {dict(full.obs['patient_id'].value_counts().sort_index())}", flush=True)
    print(f"  braf_vs_ras: {dict(full.obs['braf_vs_ras'].value_counts())}", flush=True)
    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    full.write_h5ad(OUT_H5AD, compression="gzip")
    print(f"[write] {OUT_H5AD}", flush=True)

    full.obs.to_csv(OUT_META, sep="\t")
    print(f"[write] {OUT_META}", flush=True)
    print(f"[time ] {time.time() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
