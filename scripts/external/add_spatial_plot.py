"""Add spatial overlay plot to existing pantheonos_demo runs (no MOSCOT redo)."""
import sys
from pathlib import Path
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns


def make_plot(sample_dir: Path) -> None:
    h5ad = sample_dir / "visium_with_mapped.h5ad"
    if not h5ad.exists():
        print(f"  skip {sample_dir.name}: no visium_with_mapped.h5ad")
        return
    sp = ad.read_h5ad(h5ad)
    coords = sp.obsm["spatial"]
    cats = sp.obs["mapped_celltype"].cat.categories
    palette = sns.color_palette("tab10", n_colors=len(cats))
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for i, c in enumerate(cats):
        mask = sp.obs["mapped_celltype"] == c
        axes[0].scatter(coords[mask, 0], -coords[mask, 1], s=4, c=[palette[i]],
                        label=f"{c} (n={mask.sum()})", alpha=0.8)
    axes[0].set_title(f"Mapped celltype — {sample_dir.name}")
    axes[0].legend(fontsize=8, loc="upper right")
    axes[0].set_aspect("equal"); axes[0].axis("off")
    sc_plot = axes[1].scatter(coords[:, 0], -coords[:, 1], s=4,
                               c=sp.obs["mapped_celltype_conf"].values,
                               cmap="viridis", alpha=0.8)
    axes[1].set_title("Mapping confidence")
    axes[1].set_aspect("equal"); axes[1].axis("off")
    fig.colorbar(sc_plot, ax=axes[1], shrink=0.8)
    fig.tight_layout()
    fig.savefig(sample_dir / "spatial_overlay.png", dpi=150, bbox_inches="tight")
    fig.savefig(sample_dir / "spatial_overlay.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  {sample_dir.name}: spatial_overlay.[png|pdf]")


if __name__ == "__main__":
    root = Path("project/results/pantheonos_demo")
    for d in sorted(root.iterdir()):
        if d.is_dir() and (d / "run_meta.json").exists():
            make_plot(d)
