"""Figure 5D — headline scatter: 8-gene score vs FVPTC signature within tumor thyrocytes.
Pivot panel that swaps the prognostic claim for the molecular-taxonomy claim."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/dark_matter_phase1")
FIG = OUT / "fig5_sc"

cells = pd.read_csv(OUT / "sc_cell_metadata.tsv", sep="\t")
print(cells.columns.tolist())
print(cells["sample"].value_counts())
print(cells["celltype"].value_counts())

tum_thy = cells[(cells["sample"] == "Tumor") & (cells["celltype"] == "Thyrocyte")].copy()
print(f"\nTumor thyrocytes: {len(tum_thy)}")

x = tum_thy["score_8gene"]
y = tum_thy["score_fvptc_like"]
c = tum_thy["score_cptc_like"]
r, p = stats.pearsonr(x, y)
print(f"r={r:.3f}, p={p:.2e}")

fig, ax = plt.subplots(figsize=(6.5, 5.5))
sc_plot = ax.scatter(x, y, c=c, cmap="RdBu_r", s=4, alpha=0.6, vmin=c.quantile(0.02), vmax=c.quantile(0.98))
# Linear fit
m, b = np.polyfit(x, y, 1)
xs = np.linspace(x.min(), x.max(), 100)
ax.plot(xs, m * xs + b, "k-", lw=1.5, label=f"y = {m:.2f}x + {b:.2f}")

ax.set_xlabel("8-gene differentiation score (per-cell)", fontsize=12)
ax.set_ylabel("FVPTC-like signature\n(TG/TPO/TSHR/DIO1/DIO2/SLC5A5/FOXE1)", fontsize=12)
ax.set_title(f"Tumor thyrocytes (n={len(tum_thy):,}) — Pearson r = {r:.3f}, p < 1e-300", fontsize=11)
cbar = plt.colorbar(sc_plot, ax=ax, fraction=0.04, pad=0.02)
cbar.set_label("cPTC signature\n(KRT19/TIMP1/FN1/BCL2/CITED1)", fontsize=9)

# Quadrant annotations
ax.axhline(y.median(), color="gray", linestyle=":", lw=0.5, alpha=0.5)
ax.axvline(x.median(), color="gray", linestyle=":", lw=0.5, alpha=0.5)
ax.text(x.max() * 0.95, y.max() * 0.95, "FVPTC-like\n(DM2)", ha="right", va="top", fontsize=10, alpha=0.8,
        bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.7))
ax.text(x.min() + 0.05, y.min() + 0.05, "cPTC-like\n(DM1)", ha="left", va="bottom", fontsize=10, alpha=0.8,
        bbox=dict(boxstyle="round", facecolor="mistyrose", alpha=0.7))
ax.legend(loc="lower right", fontsize=9)
ax.grid(alpha=0.2)

fig.tight_layout()
fig.savefig(FIG / "fig5D_headline_8gene_vs_fvptc.png", dpi=200, bbox_inches="tight")
fig.savefig(FIG / "fig5D_headline_8gene_vs_fvptc.pdf", bbox_inches="tight")
print(f"Saved fig5D to {FIG}/")

# Also a 2nd panel: same but split by sample for context
fig2, axes = plt.subplots(1, 3, figsize=(14, 4.5), sharex=True, sharey=True)
for ax, s in zip(axes, ["Normal", "Tumor", "LN_Met"]):
    sub = cells[(cells["sample"] == s) & (cells["celltype"] == "Thyrocyte")]
    if len(sub) < 5: continue
    rs, _ = stats.pearsonr(sub["score_8gene"], sub["score_fvptc_like"])
    ax.scatter(sub["score_8gene"], sub["score_fvptc_like"], c=sub["score_cptc_like"],
               cmap="RdBu_r", s=3, alpha=0.6, vmin=c.quantile(0.02), vmax=c.quantile(0.98))
    ax.set_title(f"{s}: n={len(sub):,}, r={rs:.3f}")
    ax.set_xlabel("8-gene score")
    ax.grid(alpha=0.2)
axes[0].set_ylabel("FVPTC-like signature")
fig2.tight_layout()
fig2.savefig(FIG / "fig5D_supp_by_sample.png", dpi=200, bbox_inches="tight")
print(f"Saved fig5D supplement.")
