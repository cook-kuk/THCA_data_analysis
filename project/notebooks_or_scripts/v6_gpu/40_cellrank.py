#!/usr/bin/env python3
"""v6 Wave 2 — RNA velocity + CellRank 2 macrostates + fate driver genes.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="h5ad with spliced/unspliced layers")
    ap.add_argument("--output_dir", default="/opt/thyroid-dash/project/results/v6_scrna")
    ap.add_argument("--n_macrostates", type=int, default=3)
    args = ap.parse_args()
    out = Path(args.output_dir); (out/"trajectory").mkdir(parents=True, exist_ok=True)

    try:
        import scanpy as sc, scvelo as scv, cellrank as cr
        import numpy as np, pandas as pd, plotly.express as px
        from cellrank.kernels import VelocityKernel, ConnectivityKernel
        from cellrank.estimators import GPCCA
    except Exception as e:
        print(f"FATAL imports: {e}", file=sys.stderr); sys.exit(2)

    adata = sc.read_h5ad(args.input)
    if "spliced" not in adata.layers:
        print("FATAL: no spliced/unspliced layers — re-process raw fastq with velocyto first", file=sys.stderr); sys.exit(2)

    print("[cellrank] computing velocity moments + dynamics …")
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)

    print("[cellrank] kernels …")
    vk = VelocityKernel(adata).compute_transition_matrix()
    ck = ConnectivityKernel(adata).compute_transition_matrix()
    combined = 0.8 * vk + 0.2 * ck

    print("[cellrank] GPCCA macrostates …")
    g = GPCCA(combined)
    g.compute_schur(n_components=20)
    g.compute_macrostates(n_states=args.n_macrostates)
    g.compute_fate_probabilities()

    adata.write(out/"trajectory/cellrank_macrostates.h5ad")

    # Driver genes per fate
    ms_names = list(g.macrostates.cat.categories)
    driver_rows = []
    for ms in ms_names:
        try:
            d = g.compute_lineage_drivers(lineages=ms)
            for _, r in d.head(50).iterrows():
                driver_rows.append({"fate":ms, "gene":r.name, "corr":float(r.iloc[0]) if len(r) else None})
        except Exception as e:
            driver_rows.append({"fate":ms, "gene":None, "corr":None, "error":str(e)})
    pd.DataFrame(driver_rows).to_csv(out/"trajectory/fate_driver_genes.tsv", sep="\t", index=False)

    # Interactive UMAP (CellRank's plot is matplotlib; we make a simple plotly UMAP colored by macrostate)
    if "X_umap" in adata.obsm:
        umap = adata.obsm["X_umap"]
        df = pd.DataFrame({"UMAP1":umap[:,0],"UMAP2":umap[:,1],
                           "macrostate":adata.obs.get("macrostates", "unassigned")})
        fig = px.scatter(df, x="UMAP1", y="UMAP2", color="macrostate",
                         title="CellRank 2 macrostates",
                         template="plotly_dark", opacity=0.7)
        fig.write_html("/opt/thyroid-dash/project/reports/html/figs_interactive/v6/trajectory_umap.html",
                       include_plotlyjs="cdn")

    Path(f"{args.output_dir}/.stamp_40_cellrank").touch()
    print("[cellrank] DONE")

if __name__ == "__main__":
    main()
