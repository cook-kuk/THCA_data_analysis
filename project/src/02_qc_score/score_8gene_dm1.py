#!/usr/bin/env python3
"""QC + within-sample z-score 기반 8-gene RAI / DM1-like / TDS / CAF/EMT/hypoxia/proliferation/epithelial scoring.

Input:  per-sample {sid}.raw.h5ad  (output of parse_gse250521.py)
Output: per-sample {sid}.scored.h5ad + combined obs table
"""
from __future__ import annotations
import argparse, logging, sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("score")

GENESETS = {
    "RAI_8":      ["TPO","DIO1","TSHR","PAX8","TG","FOXE1","NKX2-1","SLC5A5"],
    "TDS_like":   ["TG","TPO","TSHR","SLC5A5","DIO1","PAX8"],
    "CAF_ECM":    ["COL1A1","COL1A2","COL3A1","COL6A1","DCN","LUM","FAP","ACTA2","PDGFRB","POSTN"],
    "EMT":        ["VIM","FN1","SNAI1","SNAI2","ZEB1","ZEB2","TWIST1","CDH2","ITGA5","ITGB1"],
    "Hypoxia":    ["VEGFA","CA9","SLC2A1","LDHA","PDK1","BNIP3","NDRG1","EGLN3","ADM"],
    "Proliferation": ["MKI67","TOP2A","PCNA","MCM2","MCM5","STMN1","CENPF"],
    "Epithelial": ["EPCAM","KRT8","KRT18","KRT19","TACSTD2"],
}
GENE_ALIASES = {"NKX2-1": ["NKX2-1","NKX2_1","TITF1"], "SLC5A5": ["SLC5A5","NIS"]}


def resolve_gene(adata: ad.AnnData, name: str) -> str | None:
    cands = GENE_ALIASES.get(name, [name])
    var_names = set(adata.var_names)
    for c in cands:
        if c in var_names: return c
    return None


def within_sample_z(adata: ad.AnnData, gene: str) -> np.ndarray:
    x = adata[:, gene].X
    x = x.toarray().ravel() if hasattr(x, "toarray") else np.asarray(x).ravel()
    mu, sd = x.mean(), x.std()
    return (x - mu) / sd if sd > 0 else np.zeros_like(x)


def score_geneset(adata: ad.AnnData, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    found = [g for g in (resolve_gene(adata, g) for g in genes) if g]
    if not found:
        return np.full(adata.n_obs, np.nan), found
    Z = np.column_stack([within_sample_z(adata, g) for g in found])
    return Z.mean(axis=1), found


def qc_filter(adata: ad.AnnData, min_genes: int = 200, log_pfx: str = "") -> ad.AnnData:
    adata.var["mt"] = adata.var_names.str.startswith(("MT-","mt-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    n0 = adata.n_obs
    in_tissue = adata.obs["in_tissue"].fillna(0).astype(int) == 1
    pass_qc = (adata.obs["n_genes_by_counts"] >= min_genes) & in_tissue
    adata = adata[pass_qc].copy()
    log.info("%s QC: %d → %d (in_tissue + n_genes>=%d)", log_pfx, n0, adata.n_obs, min_genes)
    return adata


def normalize(adata: ad.AnnData) -> ad.AnnData:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def score_sample(in_h5ad: Path, out_h5ad: Path, qc_dir: Path) -> dict:
    adata = ad.read_h5ad(in_h5ad)
    sid = adata.obs["sample_id"].iloc[0]
    log.info("[%s] %d spots × %d genes", sid, adata.n_obs, adata.n_vars)

    adata = qc_filter(adata, log_pfx=f"[{sid}]")
    adata = normalize(adata)

    avail: dict[str, list[str]] = {}
    for set_name, genes in GENESETS.items():
        score, found = score_geneset(adata, genes)
        adata.obs[f"{set_name}_score"] = score
        adata.obs[f"{set_name}_n_genes"] = len(found)
        avail[set_name] = found
        log.info("[%s] %s: %d/%d genes (%s)", sid, set_name, len(found), len(genes), ",".join(found))

    if "RAI_8_score" in adata.obs:
        adata.obs["DM1_like_score"] = -adata.obs["RAI_8_score"]
        adata.obs["score_n_genes"] = adata.obs["RAI_8_n_genes"]

    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_h5ad, compression="gzip")

    keep_cols = ["sample_id","stage","stage_metadata_raw","gsm","total_counts","n_genes_by_counts",
                 "pct_counts_mt","array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres",
                 "in_tissue","RAI_8_score","DM1_like_score","TDS_like_score","CAF_ECM_score",
                 "EMT_score","Hypoxia_score","Proliferation_score","Epithelial_score","score_n_genes"]
    obs_out = adata.obs[[c for c in keep_cols if c in adata.obs.columns]].copy()
    obs_out["spot_id"] = adata.obs.index
    qc_dir.mkdir(parents=True, exist_ok=True)
    obs_out.to_csv(qc_dir / f"{sid}.spot_scores.tsv.gz", sep="\t", index=False, compression="gzip")
    return {"sample_id": sid, "n_spots_post_qc": adata.n_obs,
            **{f"{s}_n": len(g) for s, g in avail.items()},
            **{f"{s}_genes": ",".join(g) for s, g in avail.items()}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", default="project/data/processed/GSE250521/sample_metadata.tsv")
    ap.add_argument("--processed-dir", default="project/data/processed/GSE250521")
    ap.add_argument("--qc-dir", default="project/results/00_qc")
    ap.add_argument("--combined-out", default="project/results/01_spatial_score/all_spots_scored.tsv.gz")
    args = ap.parse_args()

    meta = pd.read_csv(args.meta, sep="\t")
    rows = []
    combined = []
    for _, r in meta[meta["status"] == "ok"].iterrows():
        sid = r["sample_id"]
        in_h = Path(r["h5ad"])
        out_h = Path(args.processed_dir) / sid / f"{sid}.scored.h5ad"
        try:
            info = score_sample(in_h, out_h, Path(args.qc_dir))
            rows.append(info)
            df = pd.read_csv(Path(args.qc_dir) / f"{sid}.spot_scores.tsv.gz", sep="\t")
            combined.append(df)
        except Exception as e:
            log.exception("[%s] failed", sid)
            rows.append({"sample_id": sid, "status": "fail", "error": str(e)})

    qc_summary = pd.DataFrame(rows)
    qc_summary.to_csv(Path(args.qc_dir) / "scoring_summary.tsv", sep="\t", index=False)
    if combined:
        all_spots = pd.concat(combined, ignore_index=True)
        Path(args.combined_out).parent.mkdir(parents=True, exist_ok=True)
        all_spots.to_csv(args.combined_out, sep="\t", index=False, compression="gzip")
        log.info("combined: %d spots × %d cols → %s", len(all_spots), all_spots.shape[1], args.combined_out)
    print(qc_summary.to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
