#!/usr/bin/env python3
"""v6 Wave 1 — classical CPU baseline pipeline.

QC -> normalize -> HVG -> PCA -> integration (Harmony, scVI optional) -> Leiden -> UMAP.
CPU-only environment; default integrator is Harmony (fast on CPU).
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

IN_H5AD = Path("/data/thca/scrna/raw/scrna_raw.h5ad")
OUT_H5AD = Path("/data/thca/scrna/processed/classical_baseline.h5ad")
LOG = Path("/opt/thyroid-dash/project/logs/v6_classical_baseline.log")


def main():
    t0 = time.time()
    sc.settings.verbosity = 1
    sc.settings.n_jobs = 4

    print(f"[load] {IN_H5AD}", flush=True)
    A = sc.read_h5ad(IN_H5AD)
    print(f"  raw: {A.n_obs} cells x {A.n_vars} genes", flush=True)

    # Use gene symbols for downstream MT detection
    if "gene_symbol" in A.var.columns:
        A.var["symbol"] = A.var["gene_symbol"].astype(str)
    else:
        A.var["symbol"] = A.var_names.astype(str)
    A.var["mt"] = A.var["symbol"].str.upper().str.startswith("MT-")

    sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    pre_qc = A.n_obs
    A = A[(A.obs["n_genes_by_counts"] > 200) &
          (A.obs["n_genes_by_counts"] < 6000) &
          (A.obs["pct_counts_mt"] < 25)].copy()
    print(f"[qc  ] {pre_qc} -> {A.n_obs} cells (200<n_genes<6000, pct_mt<25)", flush=True)

    sc.pp.filter_genes(A, min_cells=10)
    print(f"[gene-filter] {A.n_vars} genes (min_cells=10)", flush=True)

    # Stratified subsample by patient_id, cap at 20k
    target_total = 20000
    if A.n_obs > target_total:
        rng = np.random.default_rng(42)
        per_pat = A.obs.groupby("patient_id").size()
        n_pat = len(per_pat)
        per_quota = target_total // n_pat
        keep_idx = []
        for pid, n in per_pat.items():
            sel = np.where(A.obs["patient_id"].values == pid)[0]
            take = min(per_quota, len(sel))
            keep_idx.extend(rng.choice(sel, size=take, replace=False).tolist())
        A = A[sorted(keep_idx)].copy()
        print(f"[subsample] -> {A.n_obs} cells across {n_pat} patients", flush=True)

    A.layers["counts"] = A.X.copy()
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    print("[norm] normalize_total + log1p done", flush=True)

    sc.pp.highly_variable_genes(A, n_top_genes=3000, flavor="seurat", batch_key="patient_id")
    print(f"[hvg ] {int(A.var['highly_variable'].sum())} HVGs", flush=True)
    A_hvg = A[:, A.var["highly_variable"]].copy()

    sc.pp.scale(A_hvg, max_value=10)
    sc.tl.pca(A_hvg, n_comps=50, svd_solver="arpack", random_state=42)
    A.obsm["X_pca"] = A_hvg.obsm["X_pca"]
    print(f"[pca ] X_pca shape {A.obsm['X_pca'].shape}", flush=True)

    # Try Harmony first (fast on CPU). scVI on full matrix is slow on CPU.
    integrated = False
    try:
        import harmonypy as hm
        print("[harmony] running run_harmony …", flush=True)
        ho = hm.run_harmony(A.obsm["X_pca"], A.obs, "patient_id", max_iter_harmony=20)
        Z = np.asarray(ho.Z_corr)
        # harmonypy returns Z_corr with shape (n_cells, n_pcs); some versions
        # return (n_pcs, n_cells). Normalize so rows == n_cells.
        if Z.ndim == 2 and Z.shape[0] != A.n_obs and Z.shape[1] == A.n_obs:
            Z = Z.T
        A.obsm["X_int"] = Z
        integrated = True
        print(f"[harmony] X_int shape {A.obsm['X_int'].shape}", flush=True)
        A.uns["integration_method"] = "harmony"
    except Exception as e:
        print(f"[harmony-fail] {type(e).__name__}: {e}", flush=True)

    if not integrated:
        print("[fallback] using uncorrected PCA as latent", flush=True)
        A.obsm["X_int"] = A.obsm["X_pca"]
        A.uns["integration_method"] = "pca_only"

    sc.pp.neighbors(A, use_rep="X_int", n_neighbors=15, random_state=42)
    sc.tl.leiden(A, resolution=0.6, random_state=42, flavor="igraph", n_iterations=2, directed=False)
    n_clust = A.obs["leiden"].nunique()
    print(f"[leiden] {n_clust} clusters at resolution 0.6", flush=True)

    sc.tl.umap(A, random_state=42)
    print("[umap] done", flush=True)

    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    A.write_h5ad(OUT_H5AD, compression="gzip")
    print(f"[write] {OUT_H5AD}", flush=True)
    print(f"[time] total {time.time() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
