#!/usr/bin/env python3
"""v6 Wave 2 — VEGA pathway-constrained latent + per-cell pathway DIAL.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output_dir", default="/opt/thyroid-dash/project/results/v6_scrna")
    ap.add_argument("--gmt", default="/data/thca/scrna/models/h.all.v2024.1.Hs.symbols.gmt")
    ap.add_argument("--epochs", type=int, default=300)
    args = ap.parse_args()
    out = Path(args.output_dir); (out/"interpretation").mkdir(parents=True, exist_ok=True)

    try:
        import torch, scanpy as sc, numpy as np, pandas as pd
        from vega import VEGA
    except Exception as e:
        print(f"FATAL: VEGA not importable. pip install vega-tools or git+https://github.com/LucasESBS/vega.git\n{e}",
              file=sys.stderr); sys.exit(2)

    if not Path(args.gmt).exists():
        try:
            import gseapy as gp
            lib = gp.get_library("MSigDB_Hallmark_2020")
            with open(args.gmt, "w") as fh:
                for term, genes in lib.items():
                    fh.write(f"{term}\thttp://msigdb\t" + "\t".join(genes) + "\n")
        except Exception as e:
            print(f"FATAL: cannot fetch hallmark gmt: {e}", file=sys.stderr); sys.exit(2)

    adata = sc.read_h5ad(args.input)
    model = VEGA(adata=adata, gmt_path=args.gmt)
    model.train_vega(n_epochs=args.epochs)
    Z = model.to_latent()  # samples × pathways
    np.save(out/"interpretation/vega_latent.npy", Z)
    pd.DataFrame(Z, index=adata.obs_names, columns=[f"PW_{i}" for i in range(Z.shape[1])]).to_csv(
        out/"interpretation/vega_pathway_activations.tsv", sep="\t")

    # Per-cell pathway DIAL
    try:
        sys.path.insert(0,"/opt/thyroid-dash/project/notebooks_or_scripts")
        from v5p1_common import compute_dial, get_classifier_factories
        from sklearn.model_selection import LeaveOneGroupOut
        Y_col = "subtype" if "subtype" in adata.obs else "mutation_status"
        B_col = "patient_id"
        Y = adata.obs[Y_col].astype(str).values
        B = adata.obs[B_col].astype(str).values
        logo = LeaveOneGroupOut()
        splits = list(logo.split(np.arange(len(Y)), Y, groups=B))
        rows = []
        for clf_name, fac in get_classifier_factories().items():
            res = compute_dial(Z.astype("float32"), Y, B, fac, splits)
            res.update({"layer":"vega_pathway_embedding","cancer":"THCA","classifier":clf_name})
            rows.append(res)
        pd.DataFrame(rows).to_csv(out/"interpretation/pathway_dial_per_cell.tsv", sep="\t", index=False)
    except Exception as e:
        print(f"[vega] WARN DIAL skipped: {e}")

    Path(f"{args.output_dir}/.stamp_60_vega").touch()
    print("[vega] DONE")

if __name__ == "__main__":
    main()
