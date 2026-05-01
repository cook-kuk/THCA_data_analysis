"""P2-A2: GSE184362 (11 PTC patients, Pan et al — 7 tumor samples available).
Test: 8-gene ↔ FVPTC correlation in adult multi-patient PTC tumor sc."""
from pathlib import Path
import json
import warnings
warnings.filterwarnings("ignore")
import re

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = Path("/data/thca/external_sc/GSE184362/extracted")
OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/dark_matter_phase2")
FIG = OUT / "fig_p2a2"
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

# 1. Load all 7 tumor samples
prefixes = sorted(set(re.match(r"GSM\d+_(PTC\d+)_T", f.name).group(1)
                       for f in DATA.glob("*_T_matrix.mtx.gz")))
print(f"Tumor samples: {prefixes}")

adatas = []
for p in prefixes:
    matches = list(DATA.glob(f"GSM*_{p}_T_matrix.mtx.gz"))
    if not matches:
        continue
    base = matches[0].name.replace("_matrix.mtx.gz", "")
    a = sc.read_mtx(DATA / f"{base}_matrix.mtx.gz").T  # cells × genes
    barcodes = pd.read_csv(DATA / f"{base}_barcodes.tsv.gz", header=None, sep="\t")[0].values
    features = pd.read_csv(DATA / f"{base}_features.tsv.gz", header=None, sep="\t")
    a.obs_names = [f"{p}_{b}" for b in barcodes]
    a.var_names = features[1].values
    a.var_names_make_unique()
    a.obs["sample"] = p
    adatas.append(a)
    print(f"  {p}: {a.shape}")
adata = adatas[0].concatenate(adatas[1:], batch_key="batch", batch_categories=prefixes)
adata.obs["sample"] = adata.obs["batch"]
print(f"Concatenated: {adata.shape}")

# 2. QC + filter
mt_genes = adata.var_names.str.startswith("MT-")
adata.obs["pct_mt"] = (adata[:, mt_genes].X.sum(axis=1).A1 / adata.X.sum(axis=1).A1) * 100
adata.obs["n_genes"] = (adata.X > 0).sum(axis=1).A1
adata = adata[(adata.obs["n_genes"] >= 200) & (adata.obs["pct_mt"] < 25)].copy()
sc.pp.filter_genes(adata, min_cells=10)
print(f"Post-QC: {adata.shape}")

# 3. Normalize + HVG + PCA + neighbors + leiden
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat", min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=15)
sc.tl.leiden(adata, resolution=0.6, flavor="igraph", n_iterations=2, directed=False)
print(f"Leiden clusters: {adata.obs['leiden'].nunique()}")

# 4. Cell-type annotation via marker scoring
for ct, markers in [("Thyrocyte", THY_MARKERS), ("T_cell", T_MARKERS),
                    ("Myeloid", MYE_MARKERS), ("Endothelial", END_MARKERS),
                    ("Fibroblast", FIB_MARKERS)]:
    present = [g for g in markers if g in adata.raw.var_names]
    if present:
        sc.tl.score_genes(adata, gene_list=present, score_name=f"score_{ct}", use_raw=True)
score_cols = [c for c in adata.obs.columns if c.startswith("score_") and c.split("_", 1)[1] in
              {"Thyrocyte", "T_cell", "Myeloid", "Endothelial", "Fibroblast"}]
cluster_means = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
cluster_assign = cluster_means.idxmax(axis=1).str.replace("score_", "")
adata.obs["celltype"] = adata.obs["leiden"].map(cluster_assign)
print(f"\nCell type counts:\n{adata.obs['celltype'].value_counts()}")

# 5. Score 8-gene + FVPTC + cPTC
panel_present = [g for g in PANEL_8GENE if g in adata.raw.var_names]
print(f"\n8-gene panel present: {panel_present}")
sc.tl.score_genes(adata, gene_list=panel_present, score_name="score_8gene", use_raw=True)
fvptc_present = [g for g in FVPTC_GENES if g in adata.raw.var_names]
print(f"FVPTC present: {fvptc_present}")
sc.tl.score_genes(adata, gene_list=fvptc_present, score_name="score_fvptc", use_raw=True)
cptc_present = [g for g in CPTC_GENES if g in adata.raw.var_names]
print(f"cPTC present: {cptc_present}")
sc.tl.score_genes(adata, gene_list=cptc_present, score_name="score_cptc", use_raw=True)

# 6. Filter to thyrocytes (tumor cells)
thy = adata[adata.obs["celltype"] == "Thyrocyte"].copy()
print(f"\nThyrocytes: {thy.shape}")
print(f"Thyrocytes per sample:\n{thy.obs['sample'].value_counts()}")

# 7. Per-patient correlation
per_pt = []
for s, g in thy.obs.groupby("sample", observed=True):
    if len(g) >= 30:
        r, p = stats.pearsonr(g["score_8gene"], g["score_fvptc"])
        per_pt.append({"sample": s, "n_cells": len(g), "r": r, "p": p})
pp = pd.DataFrame(per_pt).sort_values("r", ascending=False)
print(f"\n=== Per-patient r ===\n{pp.to_string(index=False)}")

r_pool, p_pool = stats.pearsonr(thy.obs["score_8gene"], thy.obs["score_fvptc"])
print(f"\n=== POOLED r = {r_pool:.3f}, p < 1e-300 ===")
print(f"Median patient r = {pp['r'].median():.3f}")
print(f"r > 0.7 patients: {(pp['r'] > 0.7).sum()} / {len(pp)}")
print(f"r > 0.5 patients: {(pp['r'] > 0.5).sum()} / {len(pp)}")

# 8. Bootstrap CI
rng = np.random.default_rng(42)
boot = []
x = thy.obs["score_8gene"].values; y = thy.obs["score_fvptc"].values
for _ in range(500):
    idx = rng.integers(0, len(x), len(x))
    boot.append(stats.pearsonr(x[idx], y[idx])[0])
ci_low, ci_high = np.percentile(boot, [2.5, 97.5])
print(f"Bootstrap 95% CI: [{ci_low:.3f}, {ci_high:.3f}]")

# 9. Decision rule
verdict = "FAIL"
median_r = pp["r"].median()
n_high = int((pp["r"] > 0.7).sum())
if r_pool > 0.7 and median_r > 0.7 and n_high >= 3:
    verdict = "PASS — Nat Comm reach"
elif r_pool > 0.5 and median_r > 0.5:
    verdict = "PARTIAL — Cell Rep Med / JCI Insight"
print(f"\n=== VERDICT: {verdict} ===")

# 10. Figure: per-patient forest + pooled scatter
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
ax = axes[0]
sub = pp.sort_values("r", ascending=True)
ax.barh(np.arange(len(sub)), sub["r"], color="steelblue", alpha=0.8)
ax.axvline(0.7, color="green", linestyle="--", lw=1, alpha=0.5, label="PASS (0.7)")
ax.axvline(0.5, color="orange", linestyle="--", lw=1, alpha=0.5, label="PARTIAL (0.5)")
ax.axvline(r_pool, color="red", lw=1.5, label=f"pooled = {r_pool:.3f}")
ax.set_yticks(np.arange(len(sub)))
ax.set_yticklabels([f"{s} (n={n})" for s, n in zip(sub["sample"], sub["n_cells"])])
ax.set_xlabel("Per-patient Pearson r")
ax.set_xlim(min(-0.1, sub['r'].min()-0.05), 1.0)
ax.set_title(f"GSE184362 PTC tumor thyrocytes ({len(sub)} samples)\nPooled r = {r_pool:.3f} [{ci_low:.3f}, {ci_high:.3f}]")
ax.legend(loc="lower right", fontsize=9)
ax.grid(alpha=0.2, axis="x")

ax = axes[1]
sc_obj = ax.scatter(thy.obs["score_8gene"], thy.obs["score_fvptc"],
                    c=thy.obs["score_cptc"], cmap="RdBu_r", s=2, alpha=0.5,
                    vmin=thy.obs["score_cptc"].quantile(0.02), vmax=thy.obs["score_cptc"].quantile(0.98))
m, b = np.polyfit(thy.obs["score_8gene"], thy.obs["score_fvptc"], 1)
xs = np.linspace(thy.obs["score_8gene"].min(), thy.obs["score_8gene"].max(), 100)
ax.plot(xs, m*xs + b, "k-", lw=1.5)
ax.set_xlabel("8-gene panel score")
ax.set_ylabel("FVPTC signature")
ax.set_title(f"GSE184362 thyrocytes (n={len(thy):,})\nr = {r_pool:.3f}")
plt.colorbar(sc_obj, ax=ax, label="cPTC score")

fig.suptitle(f"P2-A2 — External validation cohort 2 (GSE184362)\n{verdict}", fontsize=12, y=1.02)
fig.tight_layout()
fig.savefig(FIG / "p2a2_gse184362_summary.png", dpi=200, bbox_inches="tight")
fig.savefig(FIG / "p2a2_gse184362_summary.pdf", bbox_inches="tight")

# 11. Save
pp.to_csv(OUT / "p2a2_per_patient_r.tsv", sep="\t", index=False)
result = {
    "cohort": "GSE184362 (Pan/Lu 2021 7 PTC tumors)",
    "n_total_cells_post_qc": int(adata.shape[0]),
    "n_thyrocytes": int(len(thy)),
    "n_patients": int(len(pp)),
    "thyrocytes_per_sample": {s: int(n) for s, n in zip(pp["sample"], pp["n_cells"])},
    "pooled_r": float(r_pool),
    "bootstrap_95CI": [float(ci_low), float(ci_high)],
    "median_patient_r": float(median_r),
    "patients_r_gt_0.7": int(n_high),
    "patients_r_gt_0.5": int((pp["r"] > 0.5).sum()),
    "FVPTC_genes_used": fvptc_present,
    "phase1_GSE241184_r": 0.905,
    "p2a1_GSE193581_r_PTC": 0.687,
    "p2a1_GSE193581_r_PTC_ATC": 0.893,
    "verdict": verdict,
}
(OUT / "p2a2_summary.json").write_text(json.dumps(result, indent=2))
print(f"\nSaved p2a2_summary.json + figures.")
print(json.dumps(result, indent=2))
