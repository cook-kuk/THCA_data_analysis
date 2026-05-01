"""GSE184362 multi-site (T/P/LN) trajectory expansion.
Extends Phase 1 finding (Normal→Tumor→LN_Met dedifferentiation) from 1 patient to 4-7 patients."""
from pathlib import Path
import re
import json
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = Path("/data/thca/external_sc/GSE184362/extracted")
OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/dark_matter_phase2")
FIG = OUT / "fig_multisite"
FIG.mkdir(exist_ok=True)
sc.settings.figdir = FIG

PANEL_8GENE = ["DIO1", "FOXE1", "NKX2-1", "PAX8", "SLC5A5", "TG", "TPO", "TSHR"]
FVPTC_GENES = ["TG", "TPO", "TSHR", "DIO1", "DIO2", "SLC5A5", "FOXE1"]
CPTC_GENES = ["KRT19", "TIMP1", "FN1", "BCL2", "CITED1"]
THY_MARKERS = ["TG", "TPO", "TSHR", "TFF3", "PAX8", "FOXE1"]
T_MARKERS = ["CD3D", "CD3E", "CD8A"]
MYE_MARKERS = ["LYZ", "CD68", "CD14"]
END_MARKERS = ["PECAM1", "VWF", "CDH5"]
FIB_MARKERS = ["COL1A1", "COL1A2", "DCN"]

# 1. Load all available samples
files = sorted(DATA.glob("*_matrix.mtx.gz"))
print(f"All matrix files: {len(files)}")

adatas = []
for fp in files:
    base = fp.name.replace("_matrix.mtx.gz", "")
    m = re.match(r"GSM\d+_(PTC\d+)_([A-Za-z]+)$", base)
    if not m:
        continue
    patient, tissue = m.group(1), m.group(2)
    a = sc.read_mtx(fp).T  # cells × genes
    barcodes = pd.read_csv(DATA / f"{base}_barcodes.tsv.gz", header=None, sep="\t")[0].values
    features = pd.read_csv(DATA / f"{base}_features.tsv.gz", header=None, sep="\t")
    a.obs_names = [f"{patient}_{tissue}_{b}" for b in barcodes]
    a.var_names = features[1].values
    a.var_names_make_unique()
    a.obs["patient"] = patient
    a.obs["tissue"] = tissue
    a.obs["sample"] = f"{patient}_{tissue}"
    adatas.append(a)
    print(f"  {patient} {tissue}: {a.shape}")
adata = adatas[0].concatenate(adatas[1:], batch_key="batch_id")
adata.obs["sample"] = adata.obs["patient"].astype(str) + "_" + adata.obs["tissue"].astype(str)
print(f"Concatenated: {adata.shape}")

# 2. QC + filter
mt_genes = adata.var_names.str.startswith("MT-")
adata.obs["pct_mt"] = (adata[:, mt_genes].X.sum(axis=1).A1 / adata.X.sum(axis=1).A1) * 100
adata.obs["n_genes"] = (adata.X > 0).sum(axis=1).A1
adata = adata[(adata.obs["n_genes"] >= 200) & (adata.obs["pct_mt"] < 25)].copy()
sc.pp.filter_genes(adata, min_cells=10)
print(f"Post-QC: {adata.shape}")

# 3. Normalize + HVG + cluster
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat", min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=15)
sc.tl.leiden(adata, resolution=0.6, flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(adata, min_dist=0.3)

# 4. Cell type
for ct, markers in [("Thyrocyte", THY_MARKERS), ("T_cell", T_MARKERS),
                    ("Myeloid", MYE_MARKERS), ("Endothelial", END_MARKERS),
                    ("Fibroblast", FIB_MARKERS)]:
    present = [g for g in markers if g in adata.raw.var_names]
    if present:
        sc.tl.score_genes(adata, gene_list=present, score_name=f"score_{ct}", use_raw=True)
score_cols = [c for c in adata.obs.columns if c.startswith("score_") and c.split("_", 1)[1] in
              {"Thyrocyte", "T_cell", "Myeloid", "Endothelial", "Fibroblast"}]
cluster_means = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
adata.obs["celltype"] = adata.obs["leiden"].map(cluster_means.idxmax(axis=1).str.replace("score_", ""))
print(f"\nCell types:\n{adata.obs['celltype'].value_counts()}")

# 5. Score signatures
panel_present = [g for g in PANEL_8GENE if g in adata.raw.var_names]
fvptc_present = [g for g in FVPTC_GENES if g in adata.raw.var_names]
cptc_present = [g for g in CPTC_GENES if g in adata.raw.var_names]
print(f"\n8-gene present: {panel_present}\nFVPTC: {fvptc_present}\ncPTC: {cptc_present}")
sc.tl.score_genes(adata, gene_list=panel_present, score_name="score_8gene", use_raw=True)
sc.tl.score_genes(adata, gene_list=fvptc_present, score_name="score_fvptc", use_raw=True)
sc.tl.score_genes(adata, gene_list=cptc_present, score_name="score_cptc", use_raw=True)

# 6. Filter to thyrocytes
thy = adata[adata.obs["celltype"] == "Thyrocyte"].copy()
print(f"\nThyrocytes: {thy.shape}")
print(f"Thyrocytes per sample:\n{thy.obs['sample'].value_counts()}")
print(f"\nThyrocytes per tissue:\n{thy.obs['tissue'].value_counts()}")

# 7. Trajectory: per-tissue 8-gene mean (per patient)
print("\n=== Per-patient × tissue 8-gene mean ===")
pivot = thy.obs.groupby(["patient", "tissue"], observed=True)["score_8gene"].mean().unstack()
print(pivot.round(3))
pivot.to_csv(OUT / "p2a2_multisite_8gene_means.tsv", sep="\t")

# 8. Per-tissue summary across patients
print("\n=== Per-tissue summary (across patients) ===")
tissue_summary = thy.obs.groupby("tissue", observed=True)[["score_8gene", "score_fvptc", "score_cptc"]].agg(
    ["mean", "std", "count"])
print(tissue_summary)

# 9. Statistical test: T vs LN dedifferentiation
t_thy = thy[thy.obs["tissue"].isin(["T"])].obs["score_8gene"]
ln_thy = thy[thy.obs["tissue"].isin(["LeftLN", "RightLN"])].obs["score_8gene"]
p_thy_score = thy[thy.obs["tissue"].isin(["P"])].obs["score_8gene"] if (thy.obs["tissue"]=="P").any() else None
print(f"\n=== Tissue-level 8-gene score ===")
print(f"T (tumor): mean={t_thy.mean():.3f}, n={len(t_thy)}")
if p_thy_score is not None:
    print(f"P (paratumor): mean={p_thy_score.mean():.3f}, n={len(p_thy_score)}")
print(f"LN (left+right): mean={ln_thy.mean():.3f}, n={len(ln_thy)}")
mw_T_LN = stats.mannwhitneyu(t_thy, ln_thy, alternative="greater")
print(f"T vs LN Mann-Whitney p={mw_T_LN.pvalue:.2e}")
if p_thy_score is not None:
    mw_P_T = stats.mannwhitneyu(p_thy_score, t_thy, alternative="greater")
    print(f"P vs T Mann-Whitney p={mw_P_T.pvalue:.2e}")

# 10. Pooled correlation 8-gene vs FVPTC across all thyrocytes
r_pool, _ = stats.pearsonr(thy.obs["score_8gene"], thy.obs["score_fvptc"])
print(f"\n=== POOLED 8-gene ↔ FVPTC r (all multi-site thyrocytes) = {r_pool:.3f} ===")

# 11. Per-patient correlation in tumor only (replication)
per_pt_T = []
for s, g in thy.obs[thy.obs["tissue"] == "T"].groupby("patient", observed=True):
    if len(g) >= 30:
        r, _ = stats.pearsonr(g["score_8gene"], g["score_fvptc"])
        per_pt_T.append({"patient": s, "n": len(g), "r": r})
per_pt_T = pd.DataFrame(per_pt_T).sort_values("r", ascending=False)
print(f"\nTumor-only per-patient r (replication of P2-A2):\n{per_pt_T.to_string(index=False)}")

# 12. Figures: trajectory + UMAP + per-tissue distributions
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
ax = axes[0]
# Per-tissue boxplot of 8-gene scores in thyrocytes
tissue_order = ["P", "T", "LeftLN", "RightLN"]
tissue_order = [t for t in tissue_order if t in thy.obs["tissue"].unique()]
data_for_box = [thy.obs.loc[thy.obs["tissue"] == t, "score_8gene"].values for t in tissue_order]
bp = ax.boxplot(data_for_box, labels=tissue_order, showfliers=False, patch_artist=True)
colors = ["lightgreen", "tab:orange", "lightcoral", "lightcoral"][:len(tissue_order)]
for patch, c in zip(bp["boxes"], colors):
    patch.set_facecolor(c); patch.set_alpha(0.7)
ax.set_ylabel("8-gene differentiation score")
ax.set_title(f"GSE184362 — Multi-site dedifferentiation\nThyrocytes per tissue, MW T>LN p={mw_T_LN.pvalue:.1e}")
ax.grid(alpha=0.2, axis="y")

# Per-patient connectivity
ax = axes[1]
for pt, sub in pivot.iterrows():
    valid_tissues = [t for t in tissue_order if t in sub.index and not pd.isna(sub.get(t))]
    if len(valid_tissues) >= 2:
        ax.plot([tissue_order.index(t) for t in valid_tissues], [sub[t] for t in valid_tissues],
                marker="o", label=pt, alpha=0.8, lw=2)
ax.set_xticks(range(len(tissue_order)))
ax.set_xticklabels(tissue_order)
ax.set_xlabel("Tissue type")
ax.set_ylabel("Mean 8-gene score (per patient)")
ax.set_title("Per-patient dedifferentiation paths")
ax.legend(loc="best", fontsize=8, ncol=2)
ax.grid(alpha=0.2)

fig.suptitle("Figure 5F (multi-site upgrade) — Tissue-level dedifferentiation trajectory in 4 PTC patients", y=1.02)
fig.tight_layout()
fig.savefig(FIG / "fig5F_multisite_trajectory.png", dpi=200, bbox_inches="tight")
fig.savefig(FIG / "fig5F_multisite_trajectory.pdf", bbox_inches="tight")

# UMAP plot of multi-site thyrocytes
sc.pl.umap(adata, color=["tissue", "celltype", "score_8gene", "score_fvptc"],
           save="_multisite_overview.png", show=False, ncols=2)
print(f"\nSaved figures to {FIG}/")

# 13. Save summary
result = {
    "n_total_cells_post_qc": int(adata.shape[0]),
    "n_thyrocytes": int(len(thy)),
    "n_patients_with_multisite": int(per_pt_T["patient"].nunique()),
    "tissue_distribution": thy.obs["tissue"].value_counts().to_dict(),
    "per_tissue_8gene_mean": thy.obs.groupby("tissue", observed=True)["score_8gene"].mean().round(3).to_dict(),
    "per_tissue_fvptc_mean": thy.obs.groupby("tissue", observed=True)["score_fvptc"].mean().round(3).to_dict(),
    "per_tissue_cptc_mean": thy.obs.groupby("tissue", observed=True)["score_cptc"].mean().round(3).to_dict(),
    "MW_T_vs_LN_p": float(mw_T_LN.pvalue),
    "pooled_r_8gene_fvptc_multisite": float(r_pool),
    "tumor_only_per_patient_r": per_pt_T.to_dict(orient="records"),
}
(OUT / "p2a2_multisite_summary.json").write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
