#!/usr/bin/env python3
"""A1 — Drop-one-out per-gene RAI_8 robustness.
For each gene g in RAI_8 (8 genes), compute RAI_7 = mean z of (RAI_8 - g), then test
DM1_like_drop_g vs THYROID_NONOVERLAP correlation across 12 external samples.

Question: 8 genes 중 1개 빼도 cross-validation 살아있나?
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from scipy.stats import spearmanr, pearsonr

RAI_8 = ["TPO", "DIO1", "TSHR", "PAX8", "TG", "FOXE1", "NKX2-1", "SLC5A5"]
NONOVERLAP = ["SLC26A4", "IYD", "DUOX1", "DUOX2", "TFF3", "HHEX", "GLIS3", "DIO2"]


def get_gene(adata, name):
    for c in (name, name.replace("-", "_")):
        if c in adata.var_names:
            x = adata[:, c].X
            return x.toarray().ravel() if hasattr(x, "toarray") else np.asarray(x).ravel()
    return None


def z(x):
    sd = x.std()
    return (x - x.mean()) / sd if sd > 0 else np.zeros_like(x)


def residualize(x, log_counts, log_ngenes):
    X = np.column_stack([log_counts, log_ngenes])
    return x - LinearRegression().fit(X, x).predict(X)


def score_set(adata, genes, log_counts, log_ngenes):
    arrs = [get_gene(adata, g) for g in genes]
    arrs = [a for a in arrs if a is not None]
    if not arrs:
        return np.full(adata.n_obs, np.nan)
    Z = np.column_stack([z(residualize(a, log_counts, log_ngenes)) for a in arrs])
    return Z.mean(axis=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--registry", default="project_external_st/results/registry/dataset_sample_registry.tsv")
    ap.add_argument("--out", default="project_external_st/results/extra/a1_drop_one_out.tsv")
    args = ap.parse_args()

    reg = pd.read_csv(args.registry, sep="\t")
    reg = reg[reg["status"] == "ok"].copy()

    # for each sample compute drop-one-out RAI scores + NONOVERLAP, save sample-mean (epi top 25%)
    rows = []
    for _, r in reg.iterrows():
        sid = r["sample_id"]; ds = r["dataset"]
        h = Path(r["h5ad"])
        a = ad.read_h5ad(h)
        # QC
        a.var["mt"] = a.var_names.str.startswith(("MT-","mt-"))
        in_tissue = a.obs["in_tissue"].fillna(0).astype(int) == 1
        ng = (a.X.getnnz(axis=1) if hasattr(a.X, "getnnz") else (a.X != 0).sum(axis=1))
        a = a[in_tissue & (ng >= 200)].copy()
        # normalize
        import scanpy as sc
        sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
        sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
        log_counts = np.log1p(a.obs["total_counts"].to_numpy())
        log_ngenes = np.log1p(a.obs["n_genes_by_counts"].to_numpy())

        # Epithelial score for top-25% mask
        epi_genes = ["EPCAM", "KRT8", "KRT18", "KRT19", "TACSTD2"]
        epi_arrs = [get_gene(a, g) for g in epi_genes]; epi_arrs = [x for x in epi_arrs if x is not None]
        epi_score = np.column_stack([z(x) for x in epi_arrs]).mean(axis=1)
        epi_mask = epi_score >= np.quantile(epi_score, 0.75)

        # NONOVERLAP score (same for all drop-one)
        nonov = score_set(a, NONOVERLAP, log_counts, log_ngenes)
        nonov_epi25_mean = nonov[epi_mask].mean()

        # Drop-each: RAI_7 then DM1_drop = -RAI_7
        for g_drop in RAI_8:
            kept = [g for g in RAI_8 if g != g_drop]
            rai7 = score_set(a, kept, log_counts, log_ngenes)
            dm1_drop = -rai7
            rows.append({"sample_id": sid, "dataset": ds, "condition": r["condition_inferred"],
                         "gene_dropped": g_drop, "n_genes_used": len(kept),
                         "DM1_drop_epi25_mean": dm1_drop[epi_mask].mean(),
                         "NONOVERLAP_epi25_mean": nonov_epi25_mean,
                         "n_spots_post_qc": a.n_obs, "n_spots_epi25": int(epi_mask.sum())})

        # also full RAI_8 for reference
        rai8 = score_set(a, RAI_8, log_counts, log_ngenes)
        dm1_full = -rai8
        rows.append({"sample_id": sid, "dataset": ds, "condition": r["condition_inferred"],
                     "gene_dropped": "(none, full RAI_8)", "n_genes_used": 8,
                     "DM1_drop_epi25_mean": dm1_full[epi_mask].mean(),
                     "NONOVERLAP_epi25_mean": nonov_epi25_mean,
                     "n_spots_post_qc": a.n_obs, "n_spots_epi25": int(epi_mask.sum())})
    df = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)

    # summary: for each gene_dropped, sample-mean correlation with NONOVERLAP
    summary = []
    for g, sub in df.groupby("gene_dropped"):
        x = sub["DM1_drop_epi25_mean"].to_numpy()
        y = sub["NONOVERLAP_epi25_mean"].to_numpy()
        ok = ~(np.isnan(x) | np.isnan(y))
        if ok.sum() < 3: continue
        rho, p = spearmanr(x[ok], y[ok])
        r, pp = pearsonr(x[ok], y[ok])
        summary.append({"gene_dropped": g, "n_samples": int(ok.sum()),
                        "spearman_rho": rho, "spearman_p": p,
                        "pearson_r": r, "pearson_p": pp})
    s = pd.DataFrame(summary).sort_values("pearson_r")
    s.to_csv(Path(args.out).with_name("a1_drop_one_out_summary.tsv"), sep="\t", index=False)
    print(s.to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
