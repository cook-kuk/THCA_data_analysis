"""Figure 5D v2 — external sc validation headline scatter, GSE193581 PTC+ATC pooled."""
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/dark_matter_phase2")
FIG = OUT / "fig_p2a"
FIG.mkdir(exist_ok=True)
a = ad.read_h5ad("/home/seungho/personal/THCA_data_analysis/project/results/v17_lu2023/GSE193581_hvg_adata.h5ad")

FVPTC_EXT = [g for g in ["TG", "TPO", "TSHR", "DIO2", "PAX8"] if g in a.var_names]
CPTC = [g for g in ["KRT19", "TIMP1", "FN1", "CITED1"] if g in a.var_names]
sc.tl.score_genes(a, gene_list=FVPTC_EXT, score_name="score_fvptc_ext")
sc.tl.score_genes(a, gene_list=CPTC, score_name="score_cptc")

mal = a[a.obs["author_celltype"] == "Malignant cell"].copy()

r_pool, _ = stats.pearsonr(mal.obs["DM_score"], mal.obs["score_fvptc_ext"])

fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), gridspec_kw={"width_ratios": [1, 1]})

# Left: PTC+ATC pooled scatter
ax = axes[0]
m, b = np.polyfit(mal.obs["DM_score"], mal.obs["score_fvptc_ext"], 1)
xs = np.linspace(mal.obs["DM_score"].min(), mal.obs["DM_score"].max(), 200)
sc_obj = ax.scatter(mal.obs["DM_score"], mal.obs["score_fvptc_ext"],
                    c=mal.obs["histology"].map({"PTC": "tab:blue", "ATC": "tab:red"}),
                    s=2, alpha=0.4)
ax.plot(xs, m*xs+b, "k-", lw=1.5, label=f"y={m:.2f}x+{b:.2f}")
ax.set_xlabel("8-gene DM score (per cell)", fontsize=11)
ax.set_ylabel("FVPTC signature\n(TG/TPO/TSHR/DIO2/PAX8)", fontsize=11)
ax.set_title(f"GSE193581 malignant cells (n={len(mal):,})\nPooled Pearson r = {r_pool:.3f}", fontsize=11)
# Custom legend: PTC vs ATC
from matplotlib.lines import Line2D
ax.legend(handles=[
    Line2D([0],[0], marker="o", color="w", markerfacecolor="tab:blue", markersize=8, label=f"PTC (n=8,590)"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor="tab:red", markersize=8, label=f"ATC (n=6,034)"),
    Line2D([0],[0], color="k", lw=1.5, label=f"y={m:.2f}x+{b:.2f}"),
], loc="upper left", fontsize=9)
ax.grid(alpha=0.2)

# Right: per-patient r forest plot
ax = axes[1]
per_pt = []
for s, g in mal.obs.groupby("sample", observed=True):
    if len(g) >= 30:
        r, p = stats.pearsonr(g["DM_score"], g["score_fvptc_ext"])
        per_pt.append({"sample": s, "n": len(g), "r": r,
                       "histology": g["histology"].iloc[0]})
pp = pd.DataFrame(per_pt).sort_values("r", ascending=True)
y_pos = np.arange(len(pp))
colors = pp["histology"].map({"PTC": "tab:blue", "ATC": "tab:red"})
ax.barh(y_pos, pp["r"], height=0.7, color=colors, alpha=0.7)
ax.axvline(0.7, color="green", linestyle="--", lw=1, alpha=0.5)
ax.axvline(0.5, color="orange", linestyle="--", lw=1, alpha=0.5)
ax.axvline(r_pool, color="black", lw=1.5, label=f"pooled = {r_pool:.3f}")
ax.set_yticks(y_pos)
ax.set_yticklabels([f"{s} ({h}, n={n})" for s, h, n in zip(pp["sample"], pp["histology"], pp["n"])], fontsize=9)
ax.set_xlabel("Per-patient Pearson r")
ax.set_xlim(-0.15, 1.0)
ax.set_title(f"Per-patient r distribution (n={len(pp)} samples)\n"
             f"PTC median r = {pp[pp['histology']=='PTC']['r'].median():.3f}", fontsize=11)
ax.legend(loc="lower right", fontsize=9)
ax.grid(alpha=0.2, axis="x")

fig.suptitle("Figure 5D v2 — External sc validation (GSE193581, Lu 2023 JCI)\n"
             "8-gene DM score recapitulates FVPTC differentiation axis across multi-patient cohort", fontsize=12, y=1.02)
fig.tight_layout()
fig.savefig(FIG / "fig5D_v2_external_validation.png", dpi=200, bbox_inches="tight")
fig.savefig(FIG / "fig5D_v2_external_validation.pdf", bbox_inches="tight")
print(f"Saved fig5D_v2_external_validation to {FIG}/")
print(f"Pooled r = {r_pool:.3f}")
