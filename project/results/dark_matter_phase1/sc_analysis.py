"""sc analysis on GSE241184 (Pu et al, 1 tumor + 1 normal + 1 LN met).

Goals:
  Step 4 — does the 8-gene signature show structure within sc thyrocytes?
           bimodality + intra- vs inter-tumor variance.
  Step 5 — pseudotime within thyrocytes; locate DM1/DM2-like cells.
  +     — DICER1/EIF1AX downstream proxies (basic miRNA/translation pathway scores).

Hard limit: only 3 patients, so we cannot do "DM1 patient vs DM2 patient" sc comparison.
We instead ask whether a SINGLE tumor's thyrocytes split into the cluster signal at
sc resolution — which is itself a meaningful claim if true.
"""
from pathlib import Path
import json
import warnings

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
DATA = Path("/data/thca/v17_gse241184")
OUT = ROOT / "results" / "dark_matter_phase1"
FIG = OUT / "fig5_sc"
FIG.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = FIG
sc.settings.set_figure_params(dpi=100, dpi_save=120, frameon=True)

PANEL_8GENE = ["DIO1", "FOXE1", "NKX2-1", "PAX8", "SLC5A5", "TG", "TPO", "TSHR"]
CELLTYPE_MARKERS = {
    "Thyrocyte": ["TG", "TPO", "TSHR", "TFF3", "PAX8", "FOXE1"],
    "T_cell": ["CD3D", "CD3E", "CD8A", "CD4"],
    "B_cell": ["MS4A1", "CD19", "CD79A"],
    "Myeloid": ["LYZ", "CD68", "CD14", "AIF1"],
    "Endothelial": ["PECAM1", "VWF", "CDH5"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN"],
}

# ============================================================================
# 1. Load all 3 samples
# ============================================================================
samples = {"Tumor": DATA / "Thyroid_tumor", "Normal": DATA / "Normal_thyroid", "LN_Met": DATA / "Lymph_node"}
adatas = []
for label, path in samples.items():
    a = sc.read_mtx(path / "matrix.mtx.gz").T  # cells × genes
    barcodes = pd.read_csv(path / "barcodes.tsv.gz", header=None, sep="\t")[0].values
    features = pd.read_csv(path / "features.tsv.gz", header=None, sep="\t")
    a.obs_names = [f"{label}_{b}" for b in barcodes]
    a.var_names = features[1].values  # symbols
    a.var["ensembl"] = features[0].values
    a.var_names_make_unique()
    a.obs["sample"] = label
    adatas.append(a)
    print(f"  {label}: {a.shape}")
adata = adatas[0].concatenate(adatas[1:], batch_key="batch", batch_categories=list(samples.keys()))
adata.obs["sample"] = adata.obs["batch"]
print(f"Concat: {adata.shape}")

# ============================================================================
# 2. QC + filter
# ============================================================================
sc.pp.calculate_qc_metrics(adata, qc_vars=[], percent_top=None, log1p=False, inplace=True)
mt_genes = adata.var_names.str.startswith("MT-")
adata.obs["pct_mt"] = (adata[:, mt_genes].X.sum(axis=1).A1 / adata.X.sum(axis=1).A1) * 100
print(f"\nQC: median n_genes_by_counts={adata.obs['n_genes_by_counts'].median():.0f}, "
      f"median pct_mt={adata.obs['pct_mt'].median():.1f}")

# Filter cells: at least 200 genes, less than 25% mt
adata = adata[(adata.obs["n_genes_by_counts"] >= 200) & (adata.obs["pct_mt"] < 25)].copy()
sc.pp.filter_genes(adata, min_cells=10)
print(f"After QC filter: {adata.shape}")

# ============================================================================
# 3. Normalize + HVG + PCA + neighbors + UMAP + leiden
# ============================================================================
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat", min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=15)
sc.tl.umap(adata, min_dist=0.3)
sc.tl.leiden(adata, resolution=0.6, flavor="igraph", n_iterations=2, directed=False)
print(f"Leiden clusters: {adata.obs['leiden'].nunique()}")

# ============================================================================
# 4. Cell-type annotation via marker scoring
# ============================================================================
for ct, markers in CELLTYPE_MARKERS.items():
    present = [g for g in markers if g in adata.raw.var_names]
    if present:
        sc.tl.score_genes(adata, gene_list=present, score_name=f"score_{ct}", use_raw=True)

# Assign cell type per cluster: argmax of mean score per cluster
score_cols = [f"score_{ct}" for ct in CELLTYPE_MARKERS]
cluster_means = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
cluster_assign = cluster_means.idxmax(axis=1).str.replace("score_", "")
adata.obs["celltype"] = adata.obs["leiden"].map(cluster_assign)
print(f"\nCell type counts:\n{adata.obs['celltype'].value_counts()}")
print(f"\nCell type × sample:\n{pd.crosstab(adata.obs['celltype'], adata.obs['sample'])}")

# ============================================================================
# 5. 8-gene panel signature score
# ============================================================================
panel_present = [g for g in PANEL_8GENE if g in adata.raw.var_names]
print(f"\n8-gene panel present in data: {panel_present}")
sc.tl.score_genes(adata, gene_list=panel_present, score_name="score_8gene", use_raw=True)

# ============================================================================
# 6. Intra-tumor heterogeneity in thyrocytes
# ============================================================================
thy = adata[adata.obs["celltype"] == "Thyrocyte"].copy()
print(f"\nThyrocytes by sample:\n{thy.obs['sample'].value_counts()}")

# 8-gene score distribution
s_summary = thy.obs.groupby("sample", observed=True)["score_8gene"].describe()
print(f"\n8-gene score in thyrocytes:\n{s_summary[['count', 'mean', 'std', 'min', '50%', 'max']]}")

# Variance decomposition: total var = within-sample + between-sample
overall_var = thy.obs["score_8gene"].var()
within_var = thy.obs.groupby("sample", observed=True)["score_8gene"].var().mean()
between_var = thy.obs.groupby("sample", observed=True)["score_8gene"].mean().var()
print(f"\nVariance decomposition (8-gene score in thyrocytes):")
print(f"  Total variance: {overall_var:.4f}")
print(f"  Mean within-sample variance: {within_var:.4f}")
print(f"  Between-sample variance: {between_var:.4f}")
print(f"  Within / Between ratio: {within_var/between_var:.2f}  (>1 = ITH dominates inter-tumor)")

# Hartigan-style dip / bimodality on tumor thyrocytes
from scipy import stats as sp_stats
tum_thy = thy[thy.obs["sample"] == "Tumor"].copy()
if len(tum_thy) >= 50:
    # Skewness and kurtosis as rough bimodality indicators; bimodality coefficient = (skew^2 + 1) / kurt
    sk = sp_stats.skew(tum_thy.obs["score_8gene"])
    ku = sp_stats.kurtosis(tum_thy.obs["score_8gene"]) + 3  # back to pearson kurtosis
    bc = (sk**2 + 1) / ku if ku > 0 else np.nan
    print(f"\nTumor thyrocytes (n={len(tum_thy)}): skew={sk:.3f}, kurtosis={ku:.3f}, bimodality_coef={bc:.3f}")
    print(f"  bimodality_coef > 0.555 suggests bimodal")

# Sub-cluster thyrocytes only
sc.pp.neighbors(thy, n_pcs=30, n_neighbors=15, use_rep="X_pca")
sc.tl.leiden(thy, resolution=0.5, flavor="igraph", n_iterations=2, directed=False, key_added="thy_subcluster")
sc.tl.umap(thy, min_dist=0.3)

# Per-subcluster 8-gene score, by sample
thy_sub_summary = thy.obs.groupby(["thy_subcluster", "sample"], observed=True)["score_8gene"].agg(["count", "mean"])
print(f"\nThyrocyte sub-cluster × sample:\n{thy_sub_summary}")

# ============================================================================
# 7. DICER1/EIF1AX downstream proxy
# ============================================================================
# DICER1 functional axis — miRNA biogenesis machinery + let-7 targets
DICER1_AXIS = ["DICER1", "DROSHA", "DGCR8", "AGO1", "AGO2", "TRBP2"]
# EIF1AX functional axis — translation initiation
EIF1AX_AXIS = ["EIF1AX", "EIF1", "EIF2S1", "EIF4E", "EIF4G1"]
# Follicular thyrocyte signature (vs papillary)
FVPTC_LIKE = ["TG", "TPO", "TSHR", "DIO1", "DIO2", "SLC5A5", "FOXE1"]  # high diff markers
CPTC_LIKE = ["KRT19", "TIMP1", "FN1", "BCL2", "CITED1"]  # cPTC markers

for name, gl in [("dicer1_axis", DICER1_AXIS), ("eif1ax_axis", EIF1AX_AXIS),
                 ("fvptc_like", FVPTC_LIKE), ("cptc_like", CPTC_LIKE)]:
    present = [g for g in gl if g in adata.raw.var_names]
    if present:
        sc.tl.score_genes(adata, gene_list=present, score_name=f"score_{name}", use_raw=True)
        print(f"  scored {name}: {len(present)}/{len(gl)} genes")

# Thyrocyte-only summary — copy new score columns from adata.obs without overwriting subcluster
for col in ["score_dicer1_axis", "score_eif1ax_axis", "score_fvptc_like", "score_cptc_like"]:
    if col in adata.obs.columns:
        thy.obs[col] = adata.obs.loc[thy.obs.index, col].values
thy_means = thy.obs.groupby("sample", observed=True)[
    ["score_8gene", "score_dicer1_axis", "score_eif1ax_axis", "score_fvptc_like", "score_cptc_like"]
].mean()
print(f"\nThyrocyte signature scores by sample:\n{thy_means}")

# Within tumor: relate 8gene score to FVPTC/cPTC scores
if len(tum_thy) > 0:
    for col in ["score_dicer1_axis", "score_eif1ax_axis", "score_fvptc_like", "score_cptc_like"]:
        if col in adata.obs.columns:
            tum_thy.obs[col] = adata.obs.loc[tum_thy.obs.index, col].values
    corr = tum_thy.obs[["score_8gene", "score_fvptc_like", "score_cptc_like",
                        "score_dicer1_axis", "score_eif1ax_axis"]].corr()
    print(f"\nIntra-tumor (thyrocytes) correlations:\n{corr.round(3)}")

# ============================================================================
# 8. Save figures
# ============================================================================
sc.pl.umap(adata, color=["sample", "celltype", "score_8gene"],
           save="_overview.png", show=False, ncols=3)
sc.pl.umap(adata, color=["score_fvptc_like", "score_cptc_like",
                          "score_dicer1_axis", "score_eif1ax_axis"],
           save="_signatures.png", show=False, ncols=2)
sc.pl.umap(thy, color=["sample", "thy_subcluster", "score_8gene",
                        "score_fvptc_like", "score_cptc_like"],
           save="_thyrocytes.png", show=False, ncols=3)

# Histograms of 8-gene score per sample (thyrocytes)
fig, ax = plt.subplots(1, 1, figsize=(6, 4))
for s in thy.obs["sample"].unique():
    vals = thy.obs.loc[thy.obs["sample"] == s, "score_8gene"]
    ax.hist(vals, bins=40, alpha=0.5, label=f"{s} (n={len(vals)})", density=True)
ax.set_xlabel("8-gene panel score (thyrocyte)")
ax.set_ylabel("density")
ax.legend()
ax.set_title("Intra-tumor heterogeneity of 8-gene signature\nGSE241184 thyrocytes")
fig.savefig(FIG / "hist_8gene_thyrocyte.png", dpi=120, bbox_inches="tight")

# Save key tables
adata.obs[["sample", "leiden", "celltype", "score_8gene", "score_fvptc_like", "score_cptc_like",
           "score_dicer1_axis", "score_eif1ax_axis"]].to_csv(OUT / "sc_cell_metadata.tsv", sep="\t")

results = {
    "n_cells_post_qc": int(adata.shape[0]),
    "n_thyrocytes": int(len(thy)),
    "thyrocytes_by_sample": thy.obs["sample"].value_counts().to_dict(),
    "8gene_score_thyrocyte_summary": s_summary.round(3).to_dict(),
    "variance_decomposition": {
        "total_var": float(overall_var),
        "mean_within_sample_var": float(within_var),
        "between_sample_var": float(between_var),
        "within_over_between": float(within_var / between_var) if between_var > 0 else None,
    },
    "tumor_thyrocyte_bimodality": {
        "n": int(len(tum_thy)),
        "skew": float(sk) if len(tum_thy) >= 50 else None,
        "kurtosis": float(ku) if len(tum_thy) >= 50 else None,
        "bimodality_coef": float(bc) if len(tum_thy) >= 50 else None,
        "interpretation": "BC > 0.555 = bimodal",
    },
    "thyrocyte_signature_means_by_sample": thy_means.round(3).to_dict(),
    "tumor_intra_correlations": corr.round(3).to_dict() if len(tum_thy) > 0 else None,
}
(OUT / "sc_step4_5_summary.json").write_text(json.dumps(results, indent=2, default=str))
print(f"\nSaved figures to {FIG}/, results to sc_step4_5_summary.json")
