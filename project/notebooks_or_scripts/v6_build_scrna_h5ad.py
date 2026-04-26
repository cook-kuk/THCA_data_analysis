#!/usr/bin/env python3
"""v6 Wave 1 — build a single AnnData from the 7 GSE184362 PTC tumor 10X samples.

NO synthetic data. Reads filtered MTX triplets per sample, concatenates with
patient_id labels in obs.
"""
from __future__ import annotations

import os
import sys
import time
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.io import mmread
from scipy.sparse import csr_matrix

RAW_DIR = Path("/data/thca/scrna/raw/gse184362")
OUT_H5AD = Path("/data/thca/scrna/raw/scrna_raw.h5ad")
OUT_META = Path("/data/thca/scrna/raw/metadata.tsv")

SAMPLES = [
    ("GSM5585102", "PTC1_T", "Patient1"),
    ("GSM5585104", "PTC2_T", "Patient2"),
    ("GSM5585107", "PTC3_T", "Patient3"),
    ("GSM5585112", "PTC5_T", "Patient5"),
    ("GSM5585117", "PTC8_T", "Patient8"),
    ("GSM5585119", "PTC9_T", "Patient9"),
    ("GSM5585121", "PTC10_T", "Patient10"),
]


def read_one(gsm: str, tag: str) -> ad.AnnData:
    pfx = RAW_DIR / f"{gsm}_{tag}"
    mtx_p = Path(str(pfx) + "_matrix.mtx.gz")
    bc_p = Path(str(pfx) + "_barcodes.tsv.gz")
    ft_p = Path(str(pfx) + "_features.tsv.gz")
    print(f"[read] {gsm} {tag}", flush=True)
    M = mmread(str(mtx_p))
    M = csr_matrix(M).T  # mtx is (genes x cells); transpose -> (cells x genes)
    with gzip.open(bc_p, "rt") as fh:
        barcodes = [ln.strip() for ln in fh]
    with gzip.open(ft_p, "rt") as fh:
        feats = [ln.strip().split("\t") for ln in fh]
    # 10X features.tsv: gene_id  gene_symbol  feature_type
    if len(feats[0]) >= 2:
        gene_ids = [r[0] for r in feats]
        gene_syms = [r[1] for r in feats]
    else:
        gene_ids = [r[0] for r in feats]
        gene_syms = gene_ids
    var = pd.DataFrame({"gene_symbol": gene_syms}, index=gene_ids)
    obs_idx = [f"{tag}_{b}" for b in barcodes]
    obs = pd.DataFrame(index=obs_idx)
    a = ad.AnnData(X=M, obs=obs, var=var)
    return a


def main():
    t0 = time.time()
    parts = []
    for gsm, tag, patient in SAMPLES:
        a = read_one(gsm, tag)
        a.obs["patient_id"] = patient
        a.obs["gsm"] = gsm
        a.obs["sample"] = tag
        parts.append(a)
        print(f"  {tag}: {a.n_obs} cells x {a.n_vars} genes", flush=True)

    print(f"[concat] {len(parts)} samples", flush=True)
    # Use 'inner' join on var to keep gene set intersection (some samples have
    # 33538 features, others 36601 — different reference annotations).
    full = ad.concat(parts, axis=0, join="inner", merge="first", index_unique=None)
    full.var_names_make_unique()
    full.obs_names_make_unique()

    # mutation_status not in GSE184362 metadata — leave NA
    full.obs["mutation_status"] = "NA"

    print(f"[done] {full.n_obs} cells x {full.n_vars} genes", flush=True)
    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    full.write_h5ad(OUT_H5AD, compression="gzip")
    print(f"[write] {OUT_H5AD}", flush=True)

    # save obs metadata as TSV
    full.obs.to_csv(OUT_META, sep="\t")
    print(f"[write] {OUT_META}", flush=True)
    print(f"[time ] {time.time() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
