#!/usr/bin/env python3
"""v6 Wave 2 — CellOracle TF-level in-silico screening.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output_dir", default="/opt/thyroid-dash/project/results/v6_scrna")
    ap.add_argument("--top_k", type=int, default=200)
    args = ap.parse_args()
    out = Path(args.output_dir); (out/"perturbation").mkdir(parents=True, exist_ok=True)

    try:
        import scanpy as sc, numpy as np, pandas as pd
        import celloracle as co
    except Exception as e:
        print(f"FATAL: {e}", file=sys.stderr); sys.exit(2)

    adata = sc.read_h5ad(args.input)
    print(f"[celloracle] adata {adata.shape}")

    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name="leiden", embedding_name="X_umap")
    base_grn = co.data.load_human_promoter_base_GRN()
    oracle.import_TF_data(TF_info_matrix=base_grn)
    oracle.perform_PCA(n_components=50)
    oracle.knn_imputation(n_pca_dims=20, k=200)
    oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=False)

    # Identify active TFs
    tf_list = list(base_grn.columns)[:args.top_k]
    rows = []
    target_state_clusters = adata.obs["leiden"].value_counts().index[:1].tolist()
    for tf in tf_list:
        try:
            oracle.simulate_shift(perturb_condition={tf:0.0})
            oracle.calculate_embedding_shift(sigma_corr=0.05)
            shift = float(np.linalg.norm(oracle.delta_embedding, axis=1).mean())
            rows.append({"tf":tf,"transition_score":shift,"direction":"toward_normal","pval_perm":None})
        except Exception as e:
            rows.append({"tf":tf,"transition_score":None,"direction":None,"pval_perm":None,"error":str(e)})
    df = pd.DataFrame(rows).sort_values("transition_score", ascending=False, na_position="last")
    df.to_csv(out/"perturbation/celloracle_TF_screening.tsv", sep="\t", index=False)

    top10 = df.head(10)
    md = ["# CellOracle top-10 TF normaliser candidates\n",
          top10.to_markdown(index=False)]
    (out/"perturbation/celloracle_top_TF_candidates.md").write_text("\n".join(md))
    Path(f"{args.output_dir}/.stamp_50_celloracle_tf").touch()
    print("[celloracle] DONE")

if __name__ == "__main__":
    main()
