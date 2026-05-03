#!/usr/bin/env python3
"""A3 — Marker-based per-spot cell-type fraction estimate.
Lightweight (no scRNA reference, no NN training) — suitable for Paper 1 supplementary
inflammation hypothesis test: is DM1-high specifically epithelial-derived, or does it
co-localize with stromal/immune signal?

Cell-type signatures (within-sample z-score, mean):
  Epithelial:  EPCAM, KRT8, KRT18, KRT19, TACSTD2
  Stromal:     COL1A1, COL3A1, DCN, LUM, FAP, ACTA2, PDGFRB
  Immune_T:    CD3D, CD3E, CD8A, CD4
  Immune_B:    MS4A1, CD79A, IGKC
  Macrophage:  CD68, CD163, LYZ
  Endothelial: PECAM1, VWF, CDH5

Per spot: top-1 score → assigned cell type. Per slide: fraction by cell type.
DM1_like_score_resid quantile (Q1/Q4) cross-tabulated by cell type to test
"DM1 high regions are immune-rich" hypothesis directly.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.linear_model import LinearRegression

CELL_MARKERS = {
    "Epithelial":  ["EPCAM", "KRT8", "KRT18", "KRT19", "TACSTD2"],
    "Stromal":     ["COL1A1", "COL3A1", "DCN", "LUM", "FAP", "ACTA2", "PDGFRB"],
    "Immune_T":    ["CD3D", "CD3E", "CD8A", "CD4"],
    "Immune_B":    ["MS4A1", "CD79A", "IGKC"],
    "Macrophage":  ["CD68", "CD163", "LYZ"],
    "Endothelial": ["PECAM1", "VWF", "CDH5"],
}


def get_x(a, g):
    if g not in a.var_names: return None
    x = a[:, g].X
    return x.toarray().ravel() if hasattr(x, "toarray") else np.asarray(x).ravel()


def z(x):
    sd = x.std()
    return (x - x.mean()) / sd if sd > 0 else np.zeros_like(x)


def residualize(x, lc, ln):
    X = np.column_stack([lc, ln])
    return x - LinearRegression().fit(X, x).predict(X)


def per_sample_celltype(adata, log_counts, log_ngenes):
    """Returns: scores dict, top1 label per spot."""
    scores = {}
    for ct, genes in CELL_MARKERS.items():
        arrs = [get_x(adata, g) for g in genes]
        arrs = [a for a in arrs if a is not None]
        if not arrs:
            scores[ct] = np.full(adata.n_obs, np.nan); continue
        Z = np.column_stack([z(residualize(x, log_counts, log_ngenes)) for x in arrs])
        scores[ct] = Z.mean(axis=1)
    Z = np.column_stack([scores[ct] for ct in CELL_MARKERS])
    top1 = np.array(list(CELL_MARKERS.keys()))[Z.argmax(axis=1)]
    return scores, top1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--registry", default="project_external_st/results/registry/dataset_sample_registry.tsv")
    ap.add_argument("--out-spot",
                    default="project_external_st/results/extra/a3_per_spot_celltype.tsv.gz")
    ap.add_argument("--out-summary",
                    default="project_external_st/results/extra/a3_celltype_x_dm1.tsv")
    args = ap.parse_args()

    reg = pd.read_csv(args.registry, sep="\t")
    reg = reg[reg["status"] == "ok"].copy()

    spot_rows = []
    for _, r in reg.iterrows():
        sid = r["sample_id"]; ds = r["dataset"]
        a = ad.read_h5ad(r["h5ad"])
        a.var["mt"] = a.var_names.str.startswith(("MT-","mt-"))
        in_tissue = a.obs["in_tissue"].fillna(0).astype(int) == 1
        ng = (a.X.getnnz(axis=1) if hasattr(a.X, "getnnz") else (a.X != 0).sum(axis=1))
        a = a[in_tissue & (ng >= 200)].copy()
        sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
        sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
        lc = np.log1p(a.obs["total_counts"].to_numpy())
        ln = np.log1p(a.obs["n_genes_by_counts"].to_numpy())
        # DM1_like_resid
        from scipy.stats import rankdata  # noqa
        rai_arrs = [get_x(a, g) for g in ["TPO","DIO1","TSHR","PAX8","TG","FOXE1","NKX2-1","SLC5A5"]]
        rai_arrs = [x for x in rai_arrs if x is not None]
        rai_resid_z = np.column_stack([z(residualize(x, lc, ln)) for x in rai_arrs]).mean(axis=1)
        dm1_resid = -rai_resid_z

        scores, top1 = per_sample_celltype(a, lc, ln)
        for i in range(a.n_obs):
            spot_rows.append({"sample_id": sid, "dataset": ds, "condition": r["condition_inferred"],
                              "DM1_like_resid": dm1_resid[i], "celltype_top1": top1[i],
                              **{f"score_{ct}": scores[ct][i] for ct in CELL_MARKERS}})
    spot_df = pd.DataFrame(spot_rows)
    Path(args.out_spot).parent.mkdir(parents=True, exist_ok=True)
    spot_df.to_csv(args.out_spot, sep="\t", index=False, compression="gzip")
    print(f"per-spot cell-type call → {args.out_spot}  ({len(spot_df)} spots)")

    # Summary: per slide, fraction by cell type, AND DM1 quartile × cell type cross-tab
    rows = []
    for sid, sub in spot_df.groupby("sample_id"):
        ds = sub["dataset"].iloc[0]; cond = sub["condition"].iloc[0]
        n = len(sub)
        for ct in CELL_MARKERS:
            frac = (sub["celltype_top1"] == ct).mean()
            rows.append({"sample_id": sid, "dataset": ds, "condition": cond,
                         "celltype": ct, "fraction_top1": frac, "n_spots": n})
    frac_df = pd.DataFrame(rows)

    # DM1 high (top 25%) vs low (bottom 25%) cell-type composition globally
    spot_df["DM1_q"] = spot_df.groupby("sample_id")["DM1_like_resid"].transform(
        lambda v: pd.qcut(v, 4, labels=["Q1_low","Q2","Q3","Q4_high"], duplicates="drop"))
    quartile_celltype = (spot_df.groupby(["DM1_q", "celltype_top1"]).size()
                         .unstack(fill_value=0))
    quartile_celltype = quartile_celltype.div(quartile_celltype.sum(axis=1), axis=0)
    quartile_celltype.to_csv(Path(args.out_summary).with_name("a3_DM1_quartile_x_celltype.tsv"),
                             sep="\t")
    frac_df.to_csv(args.out_summary, sep="\t", index=False)
    print("\n=== DM1 quartile × cell type top-1 fraction (global) ===")
    print(quartile_celltype.round(3).to_string())
    print("\n=== Per-slide cell-type fractions (head) ===")
    print(frac_df.head(30).to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
