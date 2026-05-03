#!/usr/bin/env python3
"""External-ST unified scorer (both raw within-sample z + depth-residualized).

Per sample:
- Score sets A..I (RAI_8, DM1_like=-RAI_8, TDS_overlap, THYROID_NONOVERLAP, Epithelial,
  CAF_ECM, EMT, Hypoxia, Proliferation)
- Two versions: _raw (within-sample z of log1p) + _resid (z of residuals from gene ~ log_counts + log_ngenes)
- Output spot-level table → results/scores/{dataset}_{sample}_spot_scores.tsv.gz
"""
from __future__ import annotations
import argparse, logging, sys
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.linear_model import LinearRegression

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("score_ext")

GENESETS = {
    "RAI_8":      ["TPO","DIO1","TSHR","PAX8","TG","FOXE1","NKX2-1","SLC5A5"],
    "TDS_overlap":["TG","TPO","TSHR","SLC5A5","DIO1","PAX8"],
    "THYROID_NONOVERLAP": ["SLC26A4","IYD","DUOX1","DUOX2","TFF3","HHEX","GLIS3","DIO2"],
    "CAF_ECM":    ["COL1A1","COL1A2","COL3A1","COL6A1","DCN","LUM","FAP","ACTA2","PDGFRB","POSTN"],
    "EMT":        ["VIM","FN1","SNAI1","SNAI2","ZEB1","ZEB2","TWIST1","CDH2","ITGA5","ITGB1"],
    "Hypoxia":    ["VEGFA","CA9","SLC2A1","LDHA","PDK1","BNIP3","NDRG1","EGLN3","ADM"],
    "Proliferation": ["MKI67","TOP2A","PCNA","MCM2","MCM5","STMN1","CENPF"],
    "Epithelial": ["EPCAM","KRT8","KRT18","KRT19","TACSTD2"],
}
GENE_ALIASES = {"NKX2-1": ["NKX2-1","NKX2_1","TITF1"], "SLC5A5": ["SLC5A5","NIS"]}


def resolve(adata: ad.AnnData, name: str) -> str | None:
    cands = GENE_ALIASES.get(name, [name])
    for c in cands:
        if c in adata.var_names: return c
    return None


def get_x(adata: ad.AnnData, gene: str) -> np.ndarray:
    x = adata[:, gene].X
    return x.toarray().ravel() if hasattr(x, "toarray") else np.asarray(x).ravel()


def z(x: np.ndarray) -> np.ndarray:
    sd = x.std()
    return (x - x.mean()) / sd if sd > 0 else np.zeros_like(x)


def residualize(x: np.ndarray, log_counts: np.ndarray, log_ngenes: np.ndarray) -> np.ndarray:
    X = np.column_stack([log_counts, log_ngenes])
    m = LinearRegression().fit(X, x)
    return x - m.predict(X)


def score_one_sample(in_h5ad: Path, out_tsv: Path) -> dict:
    a = ad.read_h5ad(in_h5ad)
    sid = a.obs["sample_id"].iloc[0]
    log.info("[%s] %d × %d", sid, a.n_obs, a.n_vars)

    a.var["mt"] = a.var_names.str.startswith(("MT-","mt-"))
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    n0 = a.n_obs
    in_tissue = a.obs["in_tissue"].fillna(0).astype(int) == 1
    a = a[in_tissue & (a.obs["n_genes_by_counts"] >= 200)].copy()
    log.info("[%s] QC %d→%d", sid, n0, a.n_obs)

    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)

    log_counts = np.log1p(a.obs["total_counts"].to_numpy())
    log_ngenes = np.log1p(a.obs["n_genes_by_counts"].to_numpy())

    avail: dict[str, list[str]] = {}
    for set_name, genes in GENESETS.items():
        found = [g for g in (resolve(a, g) for g in genes) if g]
        avail[set_name] = found
        if not found:
            a.obs[f"{set_name}_score_raw"] = np.nan
            a.obs[f"{set_name}_score_resid"] = np.nan
            a.obs[f"{set_name}_n_genes"] = 0; continue
        gene_x = np.column_stack([get_x(a, g) for g in found])
        z_raw = np.column_stack([z(gene_x[:, i]) for i in range(gene_x.shape[1])])
        a.obs[f"{set_name}_score_raw"] = z_raw.mean(axis=1)
        z_resid = np.column_stack(
            [z(residualize(gene_x[:, i], log_counts, log_ngenes)) for i in range(gene_x.shape[1])]
        )
        a.obs[f"{set_name}_score_resid"] = z_resid.mean(axis=1)
        a.obs[f"{set_name}_n_genes"] = len(found)

    a.obs["DM1_like_score_raw"] = -a.obs["RAI_8_score_raw"]
    a.obs["DM1_like_score_resid"] = -a.obs["RAI_8_score_resid"]

    keep_cols = ["sample_id","dataset","source_accession","condition_raw","condition_inferred","disease_axis",
                 "total_counts","n_genes_by_counts","pct_counts_mt","array_row","array_col",
                 "pxl_row_in_fullres","pxl_col_in_fullres","in_tissue"]
    keep_cols += [c for c in a.obs.columns if c.endswith("_score_raw") or c.endswith("_score_resid") or c.endswith("_n_genes")]
    keep_cols += ["DM1_like_score_raw","DM1_like_score_resid"]
    keep_cols = [c for c in dict.fromkeys(keep_cols) if c in a.obs.columns]
    df = a.obs[keep_cols].copy()
    df["spot_id"] = a.obs.index
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False, compression="gzip")
    return {"sample_id": sid, "n_spots_post_qc": a.n_obs,
            **{f"{s}_n_genes": len(g) for s, g in avail.items()},
            **{f"{s}_genes": ",".join(g) for s, g in avail.items()}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta-files", nargs="+", required=True,
                    help="sample_metadata.tsv files from parse step")
    ap.add_argument("--scores-dir", default="project_external_st/results/scores")
    ap.add_argument("--qc-out", default="project_external_st/results/qc/gene_availability_by_dataset.tsv")
    ap.add_argument("--combined-out",
                    default="project_external_st/results/scores/all_external_spots_scored.tsv.gz")
    args = ap.parse_args()

    rows = []; combined = []
    for meta_f in args.meta_files:
        meta = pd.read_csv(meta_f, sep="\t")
        for _, r in meta[meta["status"] == "ok"].iterrows():
            sid = r["sample_id"]; ds = r["dataset"]
            in_h = Path(r["h5ad"])
            out_tsv = Path(args.scores_dir) / f"{ds}_{sid}_spot_scores.tsv.gz"
            try:
                info = score_one_sample(in_h, out_tsv)
                info["dataset"] = ds; info["condition"] = r["condition_inferred"]
                rows.append(info)
                df = pd.read_csv(out_tsv, sep="\t")
                combined.append(df)
            except Exception as e:
                log.exception("[%s] failed", sid)
                rows.append({"dataset": ds, "sample_id": sid, "status": "fail", "error": str(e)})
    qc = pd.DataFrame(rows)
    Path(args.qc_out).parent.mkdir(parents=True, exist_ok=True)
    qc.to_csv(args.qc_out, sep="\t", index=False)
    if combined:
        all_df = pd.concat(combined, ignore_index=True)
        Path(args.combined_out).parent.mkdir(parents=True, exist_ok=True)
        all_df.to_csv(args.combined_out, sep="\t", index=False, compression="gzip")
        log.info("combined: %d spots × %d cols → %s", len(all_df), all_df.shape[1], args.combined_out)
    print(qc.to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
