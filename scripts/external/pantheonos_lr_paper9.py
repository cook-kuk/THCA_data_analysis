"""
PantheonOS Gallery #6 → Spatial ligand-receptor analysis on paper9 Visium.

Per Gallery description: "Ligand–receptor analysis of disease-associated tissue
microenvironments." We implement via Squidpy (scverse standard, no PantheonOS
runtime needed). The PantheonOS framing for this trajectory is: leader →
analysis_expert reads SKILL → biologist hypothesizes from LR table → reporter
writes PDF. We collapse to a single-shot script — agentic loop deferred.

Usage:
  python scripts/external/pantheonos_lr_paper9.py \\
      --visium project/results/paper9_perturbation/visium_sample01.h5ad \\
      --celltype-col niche_label \\
      --out project/results/paper9_lr/sample01/

Required dependencies:
  pip install squidpy omnipath

The cell-type column (--celltype-col) must already exist in adata.obs — typically
from the existing niche clustering (S_F31, S_F47). LR analysis interprets cell-cell
communication between those niches.

Marathon note: paper-blocking-status for paper9 is unconfirmed — confirm with author
before running. This script does NOT auto-execute.
"""
import argparse
import json
import sys
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--visium", required=True, help="Path to Visium .h5ad")
    p.add_argument("--celltype-col", required=True,
                   help="adata.obs column with discrete labels (niche cluster, cell type, etc.)")
    p.add_argument("--out", required=True, help="Output directory")
    p.add_argument("--n-perms", type=int, default=1000, help="Squidpy ligrec permutations")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--top-k", type=int, default=30, help="Top LR pairs to highlight in figures")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    visium = Path(args.visium)
    out = Path(args.out)
    if not visium.exists():
        sys.exit(f"ABORT: --visium {visium} does not exist")
    out.mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        print("[dry-run] inputs validated; LR analysis deferred")
        return

    try:
        import anndata as ad
        import scanpy as sc
        import squidpy as sq
        import pandas as pd
        import matplotlib.pyplot as plt
    except ImportError as e:
        sys.exit(f"ABORT: missing dep ({e}). Install: pip install squidpy omnipath scanpy")

    print(f"[load] {visium}")
    adata = ad.read_h5ad(visium)
    print(f"  shape: {adata.shape}")
    if args.celltype_col not in adata.obs.columns:
        sys.exit(
            f"ABORT: --celltype-col {args.celltype_col!r} not in adata.obs. "
            f"Available: {list(adata.obs.columns)[:20]}"
        )

    # Ensure clusters are categorical with deterministic order
    adata.obs[args.celltype_col] = adata.obs[args.celltype_col].astype("category")
    n_groups = adata.obs[args.celltype_col].nunique()
    print(f"[setup] {n_groups} groups in {args.celltype_col}")

    # Squidpy ligand-receptor (CellPhoneDB-style permutation test, Omnipath db)
    print(f"[ligrec] running {args.n_perms} permutations...")
    sq.gr.ligrec(
        adata,
        cluster_key=args.celltype_col,
        n_perms=args.n_perms,
        seed=args.seed,
        copy=False,
    )

    res_key = f"{args.celltype_col}_ligrec"
    if res_key not in adata.uns:
        sys.exit(f"ABORT: squidpy did not write {res_key} to adata.uns")
    res = adata.uns[res_key]

    # Means + p-values tables
    means = res["means"]
    pvals = res["pvalues"]
    means.to_csv(out / "lr_means.csv")
    pvals.to_csv(out / "lr_pvalues.csv")

    # Top-K significant pairs across all cluster pairs
    sig_mask = pvals.values < 0.05
    flat = []
    for i, src in enumerate(pvals.index):
        for j, tgt in enumerate(pvals.columns):
            if sig_mask[i, j]:
                flat.append({
                    "ligand": src[0] if isinstance(src, tuple) else src,
                    "receptor": src[1] if isinstance(src, tuple) else "",
                    "source_cluster": tgt[0] if isinstance(tgt, tuple) else "",
                    "target_cluster": tgt[1] if isinstance(tgt, tuple) else tgt,
                    "mean": means.values[i, j],
                    "pval": pvals.values[i, j],
                })
    top = pd.DataFrame(flat).sort_values(["pval", "mean"], ascending=[True, False]).head(args.top_k)
    top.to_csv(out / "lr_top_significant.csv", index=False)
    print(f"[out] {len(flat)} significant LR pairs; top-{args.top_k} → lr_top_significant.csv")

    # Heatmap dotplot
    fig, ax = plt.subplots(figsize=(10, 8))
    sq.pl.ligrec(
        adata,
        cluster_key=args.celltype_col,
        source_groups=None,
        target_groups=None,
        means_range=(0.5, None),
        pvalue_threshold=0.05,
        ax=ax,
    )
    fig.savefig(out / "lr_dotplot.png", dpi=150, bbox_inches="tight")
    fig.savefig(out / "lr_dotplot.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[out] dotplot → {out}/lr_dotplot.[png|pdf]")

    # Provenance
    (out / "run_meta.json").write_text(json.dumps({
        "source": "PantheonOS Gallery #6 spatial disease LR (squidpy port)",
        "visium": str(visium),
        "celltype_col": args.celltype_col,
        "n_perms": args.n_perms,
        "seed": args.seed,
        "n_groups": int(n_groups),
        "n_significant": int(len(flat)),
    }, indent=2))
    print(f"[done] {out}")


if __name__ == "__main__":
    main()
