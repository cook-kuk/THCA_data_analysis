"""
End-to-end PantheonOS pipeline runner: MOSCOT → mapped celltypes → Squidpy LR.

Self-contained runner that adapts the staged PantheonOS skills
(skill_sc_spatial_mapping.md, Gallery #6 LR) to our THCA cohort.

Wrapper around the 3 stub scripts but inlines the logic to avoid the abort guards
since data IS present (per user clarification 2026-05-08).

Inputs (defaults wired to GSE250521 PTC-1 + Lu 2023 GSE193581 atlas):
  --visium     project/data/processed/GSE250521/GSM7980864_PTC-1/GSM7980864_PTC-1.scored.h5ad
  --sc         project/results/v17_lu2023/GSE193581_hvg_adata.h5ad
  --celltype   author_celltype
  --out        project/results/pantheonos_demo/PTC-1/

Output:
  out/transport_matrix.npy       — MOSCOT OT solution
  out/visium_with_mapped.h5ad    — Visium + obs['mapped_celltype']
  out/lr_means.csv               — Squidpy LR means
  out/lr_pvalues.csv             — Squidpy LR p-values
  out/lr_top.csv                 — top significant LR pairs
  out/lr_dotplot.{png,pdf}       — LR dotplot
  out/run_meta.json
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--visium", default="project/data/processed/GSE250521/GSM7980864_PTC-1/GSM7980864_PTC-1.scored.h5ad")
    p.add_argument("--sc", default="project/results/v17_lu2023/GSE193581_hvg_adata.h5ad")
    p.add_argument("--celltype-col", default="author_celltype")
    p.add_argument("--histology-filter", default="PTC", help="Subset sc to one histology to match Visium")
    p.add_argument("--out", default="project/results/pantheonos_demo/PTC-1/")
    p.add_argument("--alpha", type=float, default=0.01,
                   help="moscot 0.5 requires (0,1]; small=near-linear (gene-expr dominant), 1=pure GW (structure-only)")
    p.add_argument("--tau-a", type=float, default=1.0)
    p.add_argument("--tau-b", type=float, default=0.8)
    p.add_argument("--n-perms", type=int, default=500)
    p.add_argument("--top-k", type=int, default=30)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--skip-lr", action="store_true")
    p.add_argument("--force", action="store_true", help="Re-run even if output exists")
    args = p.parse_args()

    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # Skip if already done
    if (out / "run_meta.json").exists() and not getattr(args, "force", False):
        print(f"[skip] {out}/run_meta.json already exists — use --force to redo")
        return

    print("="*60); print("Stage 1/3: load + normalize"); print("="*60)
    import anndata as ad
    sc_full = ad.read_h5ad(args.sc)
    sp = ad.read_h5ad(args.visium)
    print(f"sc total: {sc_full.shape}  histology: {sc_full.obs['histology'].value_counts().to_dict()}")
    print(f"sp:       {sp.shape}")

    # Subset sc to matching histology
    if args.histology_filter and args.histology_filter in sc_full.obs['histology'].unique():
        sc = sc_full[sc_full.obs['histology'] == args.histology_filter].copy()
        print(f"sc filtered to histology={args.histology_filter}: {sc.shape}")
    else:
        sc = sc_full

    # Common genes
    common = sc.var_names.intersection(sp.var_names)
    print(f"common genes: {len(common)}")
    sc = sc[:, common].copy()
    sp = sp[:, common].copy()

    # Normalize spatial: log1p (sc is already log1p+scaled HVG; spatial is raw)
    import scanpy as scn
    sp_x_raw = sp.X.copy()
    if sp.X.max() > 50:  # raw counts
        scn.pp.normalize_total(sp, target_sum=1e4)
        scn.pp.log1p(sp)
        print(f"sp normalized: log1p, max now {sp.X.max():.2f}")

    print("\n"+"="*60); print("Stage 2/3: MOSCOT optimal transport"); print("="*60)
    from moscot.problems.space import MappingProblem
    mp = MappingProblem(sc, sp)
    mp = mp.prepare(sc_attr={"attr": "X"})
    mp = mp.solve(alpha=args.alpha, tau_a=args.tau_a, tau_b=args.tau_b, epsilon=0.01)
    pi = np.array(mp.solutions[("src", "tgt")].transport_matrix)
    print(f"transport matrix: {pi.shape}, sum={pi.sum():.4f}")
    np.save(out / "transport_matrix.npy", pi)

    # Map celltype labels: for each spot, take argmax over weighted celltype distribution
    print("\nmapping celltypes from sc → spots...")
    labels = sc.obs[args.celltype_col].astype(str).values
    uniq, inv = np.unique(labels, return_inverse=True)
    M = np.eye(len(uniq))[inv]            # (n_sc, n_celltypes)
    # moscot 0.5: pi.shape = (n_sp_spots, n_sc_cells); pi @ M → (n_sp, n_celltypes)
    if pi.shape[1] == M.shape[0]:
        score_sums = pi.dot(M)
    else:
        score_sums = pi.T.dot(M)
    best_idx = score_sums.argmax(axis=1)
    sp.obs['mapped_celltype'] = pd.Categorical(uniq[best_idx])
    sp.obs['mapped_celltype_conf'] = score_sums.max(axis=1) / (score_sums.sum(axis=1) + 1e-12)
    print(f"mapped celltype distribution:\n{sp.obs['mapped_celltype'].value_counts()}")
    sp.write_h5ad(out / "visium_with_mapped.h5ad")

    # Spatial overlay plot — mapped celltype on tissue coords (reviewer validation)
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        coords = sp.obsm["spatial"]
        cats = sp.obs["mapped_celltype"].cat.categories
        palette = sns.color_palette("tab10", n_colors=len(cats))
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        # left: celltype
        for i, c in enumerate(cats):
            mask = sp.obs["mapped_celltype"] == c
            axes[0].scatter(coords[mask, 0], -coords[mask, 1], s=4, c=[palette[i]], label=f"{c} (n={mask.sum()})", alpha=0.8)
        axes[0].set_title(f"Mapped celltype — {Path(args.visium).stem.replace('.scored','')}")
        axes[0].legend(fontsize=8, loc="upper right")
        axes[0].set_aspect("equal"); axes[0].axis("off")
        # right: confidence
        sc_plot = axes[1].scatter(coords[:, 0], -coords[:, 1], s=4,
                                   c=sp.obs["mapped_celltype_conf"].values,
                                   cmap="viridis", alpha=0.8)
        axes[1].set_title("Mapping confidence (max score / total)")
        axes[1].set_aspect("equal"); axes[1].axis("off")
        fig.colorbar(sc_plot, ax=axes[1], shrink=0.8)
        fig.tight_layout()
        fig.savefig(out / "spatial_overlay.png", dpi=150, bbox_inches="tight")
        fig.savefig(out / "spatial_overlay.pdf", bbox_inches="tight")
        plt.close(fig)
        print(f"  spatial overlay → {out}/spatial_overlay.[png|pdf]")
    except Exception as e:
        print(f"  spatial plot failed (non-fatal): {e}")

    if args.skip_lr:
        print("[skip-lr] done at stage 2"); return

    print("\n"+"="*60); print("Stage 3/3: Squidpy ligand-receptor"); print("="*60)
    import squidpy as sq
    n_groups = sp.obs['mapped_celltype'].nunique()
    if n_groups < 2:
        print(f"WARN: only {n_groups} celltype after mapping — skipping LR")
        (out / "run_meta.json").write_text(json.dumps({
            "source": "PantheonOS port: MOSCOT skill + Gallery #6 LR",
            "visium": args.visium, "sc": args.sc,
            "histology_filter": args.histology_filter,
            "n_common_genes": int(len(common)),
            "n_sc_cells": int(sc.shape[0]),
            "n_sp_spots": int(sp.shape[0]),
            "celltype_col": args.celltype_col,
            "n_groups": int(n_groups),
            "n_significant_lr": 0,
            "lr_skipped": True,
            "elapsed_sec": round(time.time() - t0, 1),
        }, indent=2))
        return
    sq.gr.ligrec(
        sp,
        cluster_key="mapped_celltype",
        n_perms=args.n_perms,
        seed=args.seed,
        use_raw=False,
        copy=False,
    )
    res = sp.uns["mapped_celltype_ligrec"]
    means = res["means"]; pvals = res["pvalues"]
    means.to_csv(out / "lr_means.csv")
    pvals.to_csv(out / "lr_pvalues.csv")

    sig_mask = pvals.values < 0.05
    flat = []
    for i, src in enumerate(pvals.index):
        for j, tgt in enumerate(pvals.columns):
            if sig_mask[i, j]:
                flat.append({
                    "ligand_receptor": str(src),
                    "cluster_pair": str(tgt),
                    "mean": means.values[i, j],
                    "pval": pvals.values[i, j],
                })
    top = pd.DataFrame(flat).sort_values(["pval", "mean"], ascending=[True, False]).head(args.top_k)
    top.to_csv(out / "lr_top.csv", index=False)
    print(f"significant LR pairs: {len(flat)} | top-{args.top_k} → lr_top.csv")

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12, 9))
    try:
        sq.pl.ligrec(
            sp, cluster_key="mapped_celltype",
            means_range=(0.5, None),
            pvalue_threshold=0.05, ax=ax,
        )
        fig.savefig(out / "lr_dotplot.png", dpi=150, bbox_inches="tight")
        fig.savefig(out / "lr_dotplot.pdf", bbox_inches="tight")
        print(f"dotplot → {out}/lr_dotplot.[png|pdf]")
    except Exception as e:
        print(f"dotplot failed (non-fatal): {e}")
    plt.close(fig)

    (out / "run_meta.json").write_text(json.dumps({
        "source": "PantheonOS port: MOSCOT skill + Gallery #6 LR",
        "visium": args.visium,
        "sc": args.sc,
        "histology_filter": args.histology_filter,
        "n_common_genes": int(len(common)),
        "n_sc_cells": int(sc.shape[0]),
        "n_sp_spots": int(sp.shape[0]),
        "celltype_col": args.celltype_col,
        "n_groups": int(n_groups),
        "n_significant_lr": int(len(flat)),
        "elapsed_sec": round(time.time() - t0, 1),
    }, indent=2))
    print(f"\n[done] elapsed {time.time()-t0:.1f}s → {out}")


if __name__ == "__main__":
    main()
