"""
PantheonOS evolution_topact → TopACT spatial cell-type classifier on our Visium.

Self-contained: imports the bare classifier from project/external/pantheonos/evolution_topact/
WITHOUT the GA wrapper or PantheonOS runtime.

Usage:
  python scripts/external/pantheonos_topact_visium.py \\
      --visium project/data/spatial/visium_thca_sample01.h5ad \\
      --sc-ref project/data/sc_ref/thca_atlas.h5ad \\
      --out project/results/paper9_topact/sample01/

Required user inputs (no defaults — script aborts if missing):
  --visium     adata.h5ad of one Visium sample
  --sc-ref     adata.h5ad single-cell reference with .obs['celltype']
  --out        output dir

Marathon note: paper-blocking-status of TopACT for paper9 is unconfirmed.
This script does NOT auto-run — explicit invocation required.
"""
import argparse
import sys
from pathlib import Path

# Ensure we can import the staged TopACT module
ROOT = Path(__file__).resolve().parents[2]
TOPACT_STAGED = ROOT / "project" / "external" / "pantheonos" / "evolution_topact"
sys.path.insert(0, str(TOPACT_STAGED))

# The staged files are flattened (topact_classifier.py instead of topact/classifier.py),
# so we'd need to either un-flatten or import individually. Simpler: copy to a real
# topact/ package on first run.
def ensure_topact_package() -> Path:
    pkg = TOPACT_STAGED / "topact_pkg"
    pkg.mkdir(exist_ok=True)
    (pkg / "__init__.py").write_text("")
    for stem in ["classifier", "constantlookuplist", "countdata",
                 "densetools", "filtering", "sparsetools", "spatial"]:
        src = TOPACT_STAGED / f"topact_{stem}.py"
        dst = pkg / f"{stem}.py"
        if src.exists() and not dst.exists():
            dst.write_text(src.read_text())
    return pkg


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--visium", required=True, help="Path to Visium .h5ad")
    p.add_argument("--sc-ref", required=True, help="Path to single-cell reference .h5ad")
    p.add_argument("--celltype-col", default="celltype")
    p.add_argument("--out", required=True, help="Output directory")
    p.add_argument("--dry-run", action="store_true", help="Setup only, no fitting")
    args = p.parse_args()

    visium = Path(args.visium)
    sc_ref = Path(args.sc_ref)
    out = Path(args.out)
    if not visium.exists():
        sys.exit(f"ABORT: --visium {visium} does not exist")
    if not sc_ref.exists():
        sys.exit(f"ABORT: --sc-ref {sc_ref} does not exist")
    out.mkdir(parents=True, exist_ok=True)

    pkg = ensure_topact_package()
    print(f"[setup] TopACT package at {pkg}", file=sys.stderr)
    sys.path.insert(0, str(pkg.parent))
    from topact_pkg.classifier import SCDataset  # noqa: F401
    from topact_pkg.spatial import CountGrid     # noqa: F401

    if args.dry_run:
        print("[dry-run] imports OK; exit before adata load")
        return

    import anndata as ad
    sp = ad.read_h5ad(visium)
    sc = ad.read_h5ad(sc_ref)
    print(f"[load] visium: {sp.shape}, sc-ref: {sc.shape}")

    if args.celltype_col not in sc.obs.columns:
        sys.exit(f"ABORT: --celltype-col {args.celltype_col!r} not in sc-ref.obs")

    # TopACT API skeleton — see project/external/pantheonos/evolution_topact/evaluator.py
    # for canonical fit/predict pattern. Implementation pending user confirmation
    # of which Visium sample(s) to run, scoring metric, and whether to compare to
    # existing niche clusters (S_F31/S_F47).
    raise NotImplementedError(
        "Full TopACT fit/predict not wired — needs user decision on:\n"
        "  1. Visium sample selection (paper9 has 28; pick representative or all?)\n"
        "  2. SC reference (we don't currently have THCA scRNA matched to Visium)\n"
        "  3. Comparison target (replace existing niche clusters? overlay?)\n"
        "Run with --dry-run to validate setup; remove this NotImplementedError once decided."
    )


if __name__ == "__main__":
    main()
