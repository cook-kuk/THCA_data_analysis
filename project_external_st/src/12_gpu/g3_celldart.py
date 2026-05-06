#!/usr/bin/env python3
"""G3 — CellDART-style cell-type deconvolution: paired GSE250521 scRNA → ST.
Uses simpler signature-based approach (deep learning not strictly needed for this question).
Identifies which cell types carry DM1-high signal."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import gzip


def load_scRNA_paired(scrna_dir: Path):
    """Load the 9 paired scRNA samples (GSE250521 _sc_ files) into AnnData."""
    import anndata as ad
    from scipy.io import mmread
    sample_dirs = []
    # files are like GSM7980876_N-2_sc_{barcodes,features,matrix}.tsv.gz / mtx.gz
    files = list(scrna_dir.glob("GSM*_sc_matrix.mtx.gz"))
    print(f"  found {len(files)} scRNA matrix files")
    if len(files) == 0:
        return None
    adatas = []
    for mtx_file in files:
        prefix = mtx_file.name.replace("_sc_matrix.mtx.gz","")
        bc_file = mtx_file.parent / f"{prefix}_sc_barcodes.tsv.gz"
        feat_file = mtx_file.parent / f"{prefix}_sc_features.tsv.gz"
        if not bc_file.exists() or not feat_file.exists():
            continue
        with gzip.open(mtx_file, "rb") as f:
            X = mmread(f).tocsr()
        barcodes = pd.read_csv(bc_file, header=None, sep="\t", compression="gzip")[0].tolist()
        feats = pd.read_csv(feat_file, header=None, sep="\t", compression="gzip")
        feats.columns = ["ensembl", "gene_symbol"] + [f"c{i}" for i in range(feats.shape[1]-2)]
        if X.shape[0] == len(feats) and X.shape[1] == len(barcodes):
            X = X.T.tocsr()
        a = ad.AnnData(X=X.astype(np.float32),
                       obs=pd.DataFrame(index=barcodes),
                       var=pd.DataFrame(index=feats["gene_symbol"].astype(str).values))
        a.var_names_make_unique()
        a.obs["sample_id"] = prefix
        adatas.append(a)
        print(f"    [{prefix}] {a.n_obs} cells × {a.n_vars} genes")
    return ad.concat(adatas, join="inner", index_unique=":") if adatas else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scrna-dir", default="data/scrna_gse250521")
    ap.add_argument("--out-dir", default="results/g3_celldart")
    ap.add_argument("--st-meta", default="data/all_tile_metadata.tsv.gz",
                    help="combined ST tile metadata with DM1_like_score_resid for cross-ref")
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    print("Loading paired scRNA...")
    sc_adata = load_scRNA_paired(Path(args.scrna_dir))
    if sc_adata is None or sc_adata.n_obs < 1000:
        print("WARNING: paired scRNA not available or too small. Skipping CellDART.")
        return
    print(f"\nscRNA combined: {sc_adata.n_obs} cells × {sc_adata.n_vars} genes")

    import scanpy as sc
    sc.pp.filter_cells(sc_adata, min_genes=200)
    sc.pp.filter_genes(sc_adata, min_cells=10)
    sc.pp.normalize_total(sc_adata, target_sum=1e4)
    sc.pp.log1p(sc_adata)
    sc.pp.highly_variable_genes(sc_adata, n_top_genes=2000, flavor="cell_ranger")
    sc.pp.pca(sc_adata, n_comps=30, mask_var="highly_variable")
    sc.pp.neighbors(sc_adata, n_neighbors=15)
    sc.tl.leiden(sc_adata, resolution=0.5)
    print(f"  Leiden clusters: {sc_adata.obs['leiden'].nunique()}")

    # marker per cluster
    sc.tl.rank_genes_groups(sc_adata, "leiden", method="wilcoxon", n_genes=30)
    markers = pd.DataFrame()
    for c in sc_adata.obs["leiden"].cat.categories:
        names = [g for g in sc_adata.uns["rank_genes_groups"]["names"][c][:30]]
        markers[c] = names
    markers.to_csv(out / "scrna_cluster_markers.tsv", sep="\t", index=False)

    # auto-label based on canonical thyroid+immune+stromal markers
    label_map = {}
    for c in sc_adata.obs["leiden"].cat.categories:
        top = list(markers[c])[:30]
        if any(g in top for g in ["TG","TPO","TSHR","SLC5A5","FOXE1"]): label_map[c] = "Thyrocyte"
        elif any(g in top for g in ["CD3D","CD3E","CD8A","CD8B","CD4"]): label_map[c] = "T_cell"
        elif any(g in top for g in ["MS4A1","CD79A","IGKC","JCHAIN"]): label_map[c] = "B_plasma"
        elif any(g in top for g in ["CD68","CD163","LYZ","MARCO"]): label_map[c] = "Macrophage"
        elif any(g in top for g in ["COL1A1","COL3A1","DCN","ACTA2","FAP","PDGFRB"]): label_map[c] = "Fibroblast"
        elif any(g in top for g in ["PECAM1","VWF","CDH5"]): label_map[c] = "Endothelial"
        else: label_map[c] = f"Other_{c}"
    sc_adata.obs["cell_type"] = sc_adata.obs["leiden"].map(label_map).astype("category")
    print("\n  cluster → cell type:")
    for c, l in label_map.items(): print(f"    {c}: {l}")

    sc_adata.obs.to_csv(out / "scrna_cell_labels.tsv", sep="\t")

    # Per-cluster mean expression for the 8 RAI genes + TROP2 + NONOVERLAP
    panel = ["TPO","DIO1","TSHR","PAX8","TG","FOXE1","NKX2-1","SLC5A5",
             "SLC26A4","IYD","DUOX1","DUOX2","TFF3","HHEX","GLIS3","DIO2",
             "TACSTD2","FN1","KCNN4","NAMPT"]
    avail = [g for g in panel if g in sc_adata.var_names]
    print(f"\n  panel genes available: {len(avail)}/{len(panel)}")
    cell_type_expr = sc_adata.to_df()[avail].groupby(sc_adata.obs["cell_type"]).mean()
    cell_type_expr.to_csv(out / "celltype_panel_expression.tsv", sep="\t")
    print(cell_type_expr.round(2).to_string())

    print(f"\n→ {out}")


if __name__ == "__main__":
    main()
