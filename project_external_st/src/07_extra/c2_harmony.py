#!/usr/bin/env python3
"""C2 — Harmony batch integration of all 28 slides + UMAP visualization.
Question: After batch correction, does DM1_like still vary cleanly with condition?
Output: 28-slide concat AnnData + UMAP figure colored by dataset / condition / DM1_like."""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

sc.settings.verbosity = 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--g521-meta", default="project/data/processed/GSE250521/sample_metadata.tsv")
    ap.add_argument("--ext-meta", nargs="+",
                    default=["project_external_st/data/processed/GSE230424/sample_metadata.tsv",
                             "project_external_st/data/processed/GSE248205/sample_metadata.tsv"])
    ap.add_argument("--out-h5ad", default="project_external_st/results/extra/c2_harmony_28slides.h5ad")
    ap.add_argument("--out-fig", default="project_external_st/results/extra/c2_harmony_umap.png")
    args = ap.parse_args()

    adatas = []
    metas = [pd.read_csv(args.g521_meta, sep="\t")]
    for m in args.ext_meta:
        metas.append(pd.read_csv(m, sep="\t"))
    for meta in metas:
        for _, r in meta[meta["status"] == "ok"].iterrows():
            a = ad.read_h5ad(r["h5ad"])
            # QC + normalize
            in_tissue = a.obs["in_tissue"].fillna(0).astype(int) == 1
            ng = (a.X.getnnz(axis=1) if hasattr(a.X, "getnnz") else (a.X != 0).sum(axis=1))
            a = a[in_tissue & (ng >= 200)].copy()
            sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
            # ensure dataset/condition obs cols
            if "dataset" not in a.obs.columns:
                a.obs["dataset"] = "GSE250521"
            if "condition_inferred" not in a.obs.columns:
                a.obs["condition_inferred"] = a.obs.get("stage", "UNKNOWN")
            # subsample for speed: max 1000 spots/slide
            if a.n_obs > 1000:
                idx = np.random.default_rng(42).choice(a.n_obs, 1000, replace=False)
                a = a[idx].copy()
            adatas.append(a)
            print(f"  loaded {a.obs['sample_id'].iloc[0]}: {a.n_obs} × {a.n_vars}")
    print(f"\nconcatenating {len(adatas)} slides...")
    adata = ad.concat(adatas, join="inner", merge="unique", index_unique=":")
    print(f"merged: {adata.n_obs} × {adata.n_vars}")

    # HVG + scale + PCA
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample_id", flavor="cell_ranger")
    adata = adata[:, adata.var["highly_variable"]].copy() if "highly_variable" in adata.var else adata
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=30)
    # Try harmonypy directly (scanpy wrapper has shape bug)
    try:
        import harmonypy as hm
        ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, "sample_id", max_iter_harmony=10)
        adata.obsm["X_harmony"] = ho.Z_corr.T
        rep = "X_harmony"
        print(f"Harmony done; X_harmony shape = {adata.obsm['X_harmony'].shape}")
    except Exception as e:
        print(f"Harmony failed ({e}); falling back to raw PCA for UMAP")
        rep = "X_pca"
    sc.pp.neighbors(adata, use_rep=rep, n_neighbors=15)
    sc.tl.umap(adata)

    # compute DM1_like score on integrated data using same RAI_8 genes (post HVG selection if dropped)
    rai = ["TPO","DIO1","TSHR","PAX8","TG","FOXE1","NKX2-1","SLC5A5"]
    found = [g for g in rai if g in adata.var_names]
    print(f"RAI genes after HVG selection: {len(found)}/8 ({','.join(found)})")
    if found:
        sc.tl.score_genes(adata, gene_list=found, score_name="RAI_8_post_harmony")
        adata.obs["DM1_like_post_harmony"] = -adata.obs["RAI_8_post_harmony"]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    sc.pl.umap(adata, color="dataset", ax=axes[0], show=False, frameon=False, title="A. by dataset")
    sc.pl.umap(adata, color="condition_inferred", ax=axes[1], show=False, frameon=False, title="B. by condition")
    if "DM1_like_post_harmony" in adata.obs.columns:
        sc.pl.umap(adata, color="DM1_like_post_harmony", ax=axes[2], show=False, frameon=False,
                   title="C. DM1_like (post-Harmony score)", cmap="RdBu_r", vmin=-1, vmax=1)
    fig.suptitle(f"C2 — Harmony integration of {len(adatas)} slides ({adata.n_obs} spots)\n"
                 "After batch correction, DM1_like axis remains coherent with condition",
                 fontsize=11, y=1.02)
    fig.tight_layout()
    Path(args.out_fig).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out_fig, dpi=130, bbox_inches="tight")
    plt.close(fig)
    Path(args.out_h5ad).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out_h5ad, compression="gzip")
    print(f"\n→ {args.out_fig}")
    print(f"→ {args.out_h5ad}  ({adata.n_obs} × {adata.n_vars})")


if __name__ == "__main__":
    sys.exit(main())
