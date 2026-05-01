"""Phase 7 ULTIMATE: TCGA bulk PCA/UMAP + per-gene LR explainability + stemness + drug targets + ensemble model."""
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict, cross_val_score
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
OUT = ROOT / "results/dark_matter_phase2"
FIG = OUT / "web/figures"
DATA_OUT = OUT / "web/data"

expr = pd.read_csv(ROOT / "data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv", sep="\t", index_col=0)
master = pd.read_csv(OUT / "../dark_matter_phase1/tcga_dm_master_with_pfi.tsv", sep="\t")
master_v17 = pd.read_csv(ROOT / "results/v17/tables/sample_master_v17_full.tsv", sep="\t", low_memory=False)
master_v17["tcga_short"] = master_v17["sample_id"].str[:12]
master_v17_tcga = master_v17[(master_v17["dataset"] == "TCGA-THCA") & (master_v17["normal_vs_tumor"] == "tumor")]
master = master.merge(master_v17_tcga[["tcga_short", "histology_subtype", "sample_id"]].drop_duplicates("tcga_short"),
                      on="tcga_short", how="left")
master["expr_col"] = master["sample_id"]
master = master[master["expr_col"].isin(expr.columns)]

PANEL = ["DIO1", "FOXE1", "NKX2-1", "PAX8", "SLC5A5", "TG", "TPO", "TSHR"]
panel_present = [g for g in PANEL if g in expr.index]

# ============================================================================
# 1. TCGA bulk PCA / UMAP (DM1 vs DM2 visual)
# ============================================================================
print("=== 1. TCGA bulk PCA on Dark Matter cohort ===")
dm = master[master["dm_status"] & master["v17_dark_cluster"].isin(["DM1", "DM2"])].copy()
print(f"DM with expr: {len(dm)}")

# Use top 2000 most variable genes for PCA
gene_var = expr.var(axis=1).sort_values(ascending=False)
top_var = gene_var.head(2000).index.tolist()
X_full = expr.loc[top_var, dm["expr_col"]].T.values
y = (dm["v17_dark_cluster"] == "DM2").astype(int).values

X_scaled = StandardScaler().fit_transform(X_full)
pca = PCA(n_components=10, random_state=42)
X_pca = pca.fit_transform(X_scaled)
print(f"PCA var explained PC1-PC5: {pca.explained_variance_ratio_[:5].round(3)}")
print(f"Cumulative: {np.cumsum(pca.explained_variance_ratio_[:5]).round(3)}")

# UMAP if available
try:
    import umap
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.3, random_state=42)
    X_umap = reducer.fit_transform(X_pca)
    has_umap = True
except ImportError:
    print("umap-learn not available, skipping UMAP")
    has_umap = False
    X_umap = None

# Plot
ncols = 3 if has_umap else 2
fig, axes = plt.subplots(1, ncols, figsize=(5*ncols, 5))
ax = axes[0]
colors_y = ["#4361ee" if v == 0 else "#f4a261" for v in y]
ax.scatter(X_pca[:, 0], X_pca[:, 1], c=colors_y, s=50, alpha=0.8, edgecolors="white", linewidths=0.5)
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([0],[0], marker="o", color="w", markerfacecolor="#4361ee", markersize=10, label=f"DM1 (n={int((y==0).sum())})"),
                    Line2D([0],[0], marker="o", color="w", markerfacecolor="#f4a261", markersize=10, label=f"DM2 (n={int((y==1).sum())})")])
ax.set_title("PCA (top 2000 variable genes, n=136)")
ax.grid(alpha=0.3)

# Variance scree
ax = axes[1]
ax.bar(range(1, 11), pca.explained_variance_ratio_*100, color="steelblue", alpha=0.85)
ax.set_xlabel("PC"); ax.set_ylabel("% variance explained")
ax.set_title("PCA scree plot")
ax.grid(alpha=0.3)

if has_umap:
    ax = axes[2]
    ax.scatter(X_umap[:, 0], X_umap[:, 1], c=colors_y, s=50, alpha=0.8, edgecolors="white", linewidths=0.5)
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    ax.set_title("UMAP from PCA(10)")
    ax.legend(handles=[Line2D([0],[0], marker="o", color="w", markerfacecolor="#4361ee", markersize=10, label="DM1"),
                       Line2D([0],[0], marker="o", color="w", markerfacecolor="#f4a261", markersize=10, label="DM2")])
    ax.grid(alpha=0.3)

fig.suptitle("Figure E22 — TCGA-THCA bulk dim-reduction within Dark Matter (DM1 vs DM2)", y=1.02)
fig.tight_layout()
fig.savefig(FIG / "figE22_pca_umap_dm.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved figE22")

# ============================================================================
# 2. Per-gene LR coefficients (explainability for FVPTC prediction)
# ============================================================================
print("\n=== 2. Per-gene LR explainability ===")
dm_h = dm[dm["histology_subtype"].isin(["cPTC", "FVPTC"])].copy()
X_panel = expr.loc[panel_present, dm_h["expr_col"]].T.values
y_h = (dm_h["histology_subtype"] == "FVPTC").astype(int).values

# Standardize
X_scaled_p = StandardScaler().fit_transform(X_panel)
clf = LogisticRegression(max_iter=2000, C=1.0)
clf.fit(X_scaled_p, y_h)
coefs = pd.DataFrame({"gene": panel_present, "coef_FVPTCup": clf.coef_[0],
                      "abs_coef": np.abs(clf.coef_[0])}).sort_values("abs_coef", ascending=False)
print(coefs.round(3))
coefs.to_csv(DATA_OUT / "lr_panel_coefficients.tsv", sep="\t", index=False)

# Plot
fig, ax = plt.subplots(figsize=(8, 4.5))
colors_c = ["tab:orange" if v > 0 else "tab:blue" for v in coefs["coef_FVPTCup"]]
ax.barh(coefs["gene"], coefs["coef_FVPTCup"], color=colors_c, alpha=0.85)
ax.axvline(0, color="black", lw=0.5)
ax.set_xlabel("LR coefficient (positive = pushes toward FVPTC)")
ax.set_title("Figure E23 — 8-gene logistic regression coefficients (FVPTC vs cPTC)\n(standardized features, all 8 genes contribute)")
for i, (_, r) in enumerate(coefs.iterrows()):
    ax.text(r["coef_FVPTCup"] + (0.02 if r["coef_FVPTCup"] >= 0 else -0.02), i,
            f"{r['coef_FVPTCup']:.2f}", va="center", fontsize=9,
            ha="left" if r["coef_FVPTCup"] >= 0 else "right")
ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
fig.savefig(FIG / "figE23_lr_coefficients.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved figE23")

# ============================================================================
# 3. Stemness score (mRNAsi-like) by cluster
# ============================================================================
print("\n=== 3. Stemness mRNAsi-like score by cluster ===")
# Embryonic stem cell signature (Malta 2018 Cell mRNAsi proxy genes)
STEMNESS = ["LIN28B", "LIN28A", "SOX2", "POU5F1", "NANOG", "ZFP42", "DPPA3", "DPPA5",
            "DNMT3B", "TET1", "BMI1", "MYC", "CDH1", "TERT", "CCND1"]
avail = [g for g in STEMNESS if g in expr.index]
print(f"Stemness available: {avail}")
if len(avail) >= 5:
    stem = expr.loc[avail].mean(axis=0)
    master["stemness"] = master["expr_col"].map(stem)
    s_dm1 = master.loc[master["v17_dark_cluster"] == "DM1", "stemness"].dropna()
    s_dm2 = master.loc[master["v17_dark_cluster"] == "DM2", "stemness"].dropna()
    s_braf = master.loc[master["fourway"] == "BRAF V600E", "stemness"].dropna() if "fourway" in master else None
    t, p = stats.ttest_ind(s_dm1, s_dm2, equal_var=False)
    d_stem = (s_dm2.mean() - s_dm1.mean()) / np.sqrt((s_dm1.std()**2 + s_dm2.std()**2)/2)
    print(f"DM1 stem mean={s_dm1.mean():.3f}, DM2={s_dm2.mean():.3f}, p={p:.4f}, d={d_stem:.2f}")

# ============================================================================
# 4. Drug target / repurposing — DM2-up DEG match
# ============================================================================
print("\n=== 4. Drug target match for DM2-up genes ===")
# Curated FDA-approved or clinical-trial drug-target list
DRUG_TARGETS = {
    "TPO": ["Methimazole (anti-thyroid, indirect inhibitor)"],
    "TG": ["Thyroglobulin antibody therapy (experimental)"],
    "TSHR": ["TSH inhibitors / receptor antagonists (CS-17, NCGC00161856 — preclinical)"],
    "DIO1": ["Iopanoic acid, propylthiouracil (deiodinase inhibitor)"],
    "SLC5A5": ["Iodide transport modulators (perchlorate)"],
    "BRAF": ["dabrafenib, vemurafenib"],
    "TERT": ["imetelstat (TERT inhibitor, experimental)"],
    "RET": ["selpercatinib, pralsetinib"],
    "NTRK1": ["larotrectinib"],
    "NTRK3": ["entrectinib"],
    "ALK": ["crizotinib, alectinib"],
    "EGFR": ["erlotinib, gefitinib"],
    "PIK3CA": ["alpelisib"],
    "MTOR": ["everolimus, sirolimus"],
    "MET": ["cabozantinib"],
    "CDK4": ["palbociclib, abemaciclib"],
    "CDK6": ["palbociclib, abemaciclib"],
    "BCL2": ["venetoclax"],
    "PARP1": ["olaparib, niraparib"],
    "ESR1": ["tamoxifen, fulvestrant"],
    "AR": ["enzalutamide"],
    "CDH1": ["(no direct, structural)"],
    "MYC": ["(no direct; BET inhibitors indirect)"],
    "HMGA2": ["(no direct; let-7 mimic experimental)"],
    "LIN28B": ["(no direct)"],
}

dge = pd.read_csv(DATA_OUT / "dge_dm1_vs_dm2.tsv", sep="\t")
top_up = dge.nlargest(30, "log2FC_DM2vDM1")
matches = []
for _, r in top_up.iterrows():
    g = r["gene"]
    if g in DRUG_TARGETS:
        matches.append({"gene": g, "log2FC": r["log2FC_DM2vDM1"], "p": r["p"], "drugs": "; ".join(DRUG_TARGETS[g])})
if matches:
    drugs_df = pd.DataFrame(matches)
    print(drugs_df.to_string(index=False))
    drugs_df.to_csv(DATA_OUT / "drug_targets_dm2_up.tsv", sep="\t", index=False)

# Top down genes (DM1-up = potential DM1-specific targets)
top_down = dge.nsmallest(30, "log2FC_DM2vDM1")
matches_dn = []
for _, r in top_down.iterrows():
    g = r["gene"]
    if g in DRUG_TARGETS:
        matches_dn.append({"gene": g, "log2FC": r["log2FC_DM2vDM1"], "p": r["p"], "drugs": "; ".join(DRUG_TARGETS[g])})

# ============================================================================
# 5. Ensemble model (8-gene + age + stage_ord)
# ============================================================================
print("\n=== 5. Ensemble model with age + stage ===")
def stage_ord(s):
    if not isinstance(s, str): return np.nan
    s = s.strip().upper()
    if "IV" in s: return 4
    if "III" in s: return 3
    if "II" in s and "III" not in s and "IV" not in s: return 2
    if "I" in s: return 1
    return np.nan

dm_h["age_yrs"] = pd.to_numeric(dm_h["age"], errors="coerce")
if "ajcc_pathologic_tumor_stage" in dm_h.columns:
    dm_h["stage_ord"] = dm_h["ajcc_pathologic_tumor_stage"].apply(stage_ord)
else:
    dm_h["stage_ord"] = dm_h["stage"].apply(stage_ord)
mask = dm_h[["age_yrs", "stage_ord"]].notna().all(axis=1)
print(f"With age+stage: {int(mask.sum())} / {len(dm_h)}")

# Compare AUC: 8-gene vs 8-gene+age vs 8-gene+age+stage
results = {}
for label, feat in [("8-gene only", panel_present),
                     ("8-gene + age", panel_present + ["age_yrs"]),
                     ("8-gene + age + stage", panel_present + ["age_yrs", "stage_ord"])]:
    sub = dm_h[mask].copy()
    if "age_yrs" in feat:
        X_arr = np.column_stack([expr.loc[panel_present, sub["expr_col"]].T.values, sub["age_yrs"].values.reshape(-1, 1)])
        if "stage_ord" in feat:
            X_arr = np.column_stack([X_arr, sub["stage_ord"].values.reshape(-1, 1)])
    else:
        X_arr = expr.loc[panel_present, sub["expr_col"]].T.values
    y_e = (sub["histology_subtype"] == "FVPTC").astype(int).values
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    clf = LogisticRegression(max_iter=2000, C=1.0)
    auc_scores = cross_val_score(clf, StandardScaler().fit_transform(X_arr), y_e, cv=skf, scoring="roc_auc")
    auc_mean = auc_scores.mean()
    auc_std = auc_scores.std()
    results[label] = {"AUC_mean": float(auc_mean), "AUC_std": float(auc_std), "n_features": X_arr.shape[1], "n_samples": len(y_e)}
    print(f"  {label}: AUC = {auc_mean:.3f} ± {auc_std:.3f} (n={len(y_e)} patients)")

# ============================================================================
# Final figure — ensemble + drugs combined
# ============================================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 4.5))
ax = axes[0]
labels = list(results.keys())
means = [results[l]["AUC_mean"] for l in labels]
stds = [results[l]["AUC_std"] for l in labels]
colors_e = ["#4361ee", "#9d4edd", "#e63946"]
ax.bar(labels, means, yerr=stds, color=colors_e, alpha=0.85, capsize=8)
for i, (m, s) in enumerate(zip(means, stds)):
    ax.text(i, m + s + 0.01, f"{m:.3f} ± {s:.3f}", ha="center", fontsize=10, fontweight="bold")
ax.set_ylabel("ROC AUC (5-fold CV)")
ax.set_ylim(0.5, 1.0)
ax.set_title("Figure E24A — Ensemble model: 8-gene + clinical covariates")
ax.grid(alpha=0.3, axis="y")
plt.setp(ax.get_xticklabels(), rotation=15, ha="right")

ax = axes[1]
if matches:
    drug_top = pd.DataFrame(matches).head(10)
    ax.barh(drug_top["gene"], drug_top["log2FC"], color="tab:orange", alpha=0.85)
    for i, (_, r) in enumerate(drug_top.iterrows()):
        ax.text(r["log2FC"] + 0.05, i, "  " + r["drugs"][:35] + ("..." if len(r["drugs"]) > 35 else ""),
                va="center", fontsize=8)
    ax.set_xlabel("log2FC (DM2-up)")
    ax.set_title("Figure E24B — DM2-up druggable targets")
    ax.grid(alpha=0.3, axis="x")
fig.tight_layout()
fig.savefig(FIG / "figE24_ensemble_drugs.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved figE24")

# Save summary
final = {
    "PCA_explained_var_PC1_5": pca.explained_variance_ratio_[:5].round(4).tolist(),
    "PCA_cumulative_PC1_5": np.cumsum(pca.explained_variance_ratio_[:5]).round(4).tolist(),
    "lr_coefficients": coefs.round(3).to_dict(orient="records"),
    "stemness_DM1_vs_DM2": {"DM1_mean": float(s_dm1.mean()), "DM2_mean": float(s_dm2.mean()),
                            "delta": float(s_dm2.mean() - s_dm1.mean()), "p": float(p), "Cohen_d": float(d_stem)},
    "drug_targets_DM2_up": matches,
    "drug_targets_DM1_up": matches_dn,
    "ensemble_AUC": results,
}
(DATA_OUT / "p8_ultimate_summary.json").write_text(json.dumps(final, indent=2, default=str))
print(f"\n=== ALL P8 DONE ===")
