"""
PantheonOS skill_sc_spatial_mapping → MOSCOT optimal-transport sc → spatial.

Stub-only: we do NOT currently have matched sc + spatial THCA pairs. This is
a reviewer-reserve scaffold. The recipe is verbatim from:
  project/external/pantheonos/single_cell_spatial_analysis/skill_sc_spatial_mapping.md

Usage (when matched data appears):
  pip install moscot
  python scripts/external/pantheonos_moscot_stub.py \\
      --sc project/data/sc/thca_sc.h5ad \\
      --sp project/data/spatial/thca_visium.h5ad \\
      --out project/results/paper9_moscot/

Marathon note: this script is a placeholder. Reviewer-reserve only.
"""
import argparse
import sys
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sc", required=True, help="Single-cell adata .h5ad")
    p.add_argument("--sp", required=True, help="Spatial adata .h5ad")
    p.add_argument("--out", required=True)
    p.add_argument("--alpha", type=float, default=0.0, help="0=linear OT, 1=Gromov-Wasserstein")
    p.add_argument("--tau-a", type=float, default=1.0)
    p.add_argument("--tau-b", type=float, default=0.8)
    p.add_argument("--celltype-col", default="celltype")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    sc_path = Path(args.sc)
    sp_path = Path(args.sp)
    out = Path(args.out)
    if not sc_path.exists():
        sys.exit(f"ABORT: --sc {sc_path} does not exist (THCA sc-RNA not yet acquired)")
    if not sp_path.exists():
        sys.exit(f"ABORT: --sp {sp_path} does not exist")

    if args.dry_run:
        print("[dry-run] paths exist; MOSCOT mapping wiring deferred to runtime")
        return

    try:
        import anndata as ad
        from moscot.problems.space import MappingProblem
        import numpy as np
    except ImportError as e:
        sys.exit(f"ABORT: missing dep ({e}). Install with: pip install moscot anndata")

    out.mkdir(parents=True, exist_ok=True)
    adata_sc = ad.read_h5ad(sc_path)
    adata_sp = ad.read_h5ad(sp_path)

    # Common-gene filter (cell cycle excluded per skill recipe)
    common = adata_sc.var_names.intersection(adata_sp.var_names)
    adata_sc = adata_sc[:, common].copy()
    adata_sp = adata_sp[:, common].copy()

    mp = MappingProblem(adata_sc, adata_sp)
    mp.solve(alpha=args.alpha, tau_a=args.tau_a, tau_b=args.tau_b)

    pi = np.array(mp.solutions[("src", "tgt")].transport_matrix)
    np.save(out / "transport_matrix.npy", pi)

    # Imputed gene expression on spatial coords
    imputed = pi.dot(adata_sc.X)
    adata_pred = ad.AnnData(X=imputed, obsm=adata_sp.obsm.copy(), obs=adata_sp.obs.copy())
    adata_pred.var_names = adata_sc.var_names
    adata_pred.write_h5ad(out / "spatial_imputed.h5ad")

    print(f"Wrote transport matrix + imputed h5ad to {out}")


if __name__ == "__main__":
    main()
