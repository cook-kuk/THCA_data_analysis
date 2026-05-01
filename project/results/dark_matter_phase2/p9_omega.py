"""Phase 8 OMEGA: ROC overlay + gene-gene corr + dot plot + Cox interaction + confusion matrix + Sankey."""
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix, brier_score_loss
from sklearn.preprocessing import StandardScaler
from lifelines import CoxPHFitter

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

dm = master[master["dm_status"] & master["v17_dark_cluster"].isin(["DM1", "DM2"])].copy()
dm_h = dm[dm["histology_subtype"].isin(["cPTC", "FVPTC"])].copy()
y = (dm_h["histology_subtype"] == "FVPTC").astype(int).values

# ============================================================================
# 1. Multi-model ROC overlay (8-gene LR vs TPO alone vs DIO1 alone vs random)
# ============================================================================
print("=== 1. Multi-model ROC overlay ===")
X_full = expr.loc[panel_present, dm_h["expr_col"]].T.values
X_scaled = StandardScaler().fit_transform(X_full)
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
clf = LogisticRegression(max_iter=2000, C=1.0)
y_proba_8gene = cross_val_predict(clf, X_scaled, y, cv=skf, method="predict_proba")[:, 1]
fpr_8, tpr_8, _ = roc_curve(y, y_proba_8gene)
auc_8 = roc_auc_score(y, y_proba_8gene)

# Single genes
single_results = {}
fig, ax = plt.subplots(figsize=(8, 7))
ax.plot(fpr_8, tpr_8, lw=3, color="#4361ee", label=f"8-gene panel CV AUC = {auc_8:.3f}")
colors_g = plt.cm.Set2(np.linspace(0, 1, len(panel_present)))
for i, g in enumerate(panel_present):
    gx = expr.loc[g, dm_h["expr_col"]].values
    auc_g = roc_auc_score(y, gx)
    fpr_g, tpr_g, _ = roc_curve(y, gx)
    if auc_g < 0.5:
        # Reverse direction
        auc_g = 1 - auc_g
        fpr_g, tpr_g, _ = roc_curve(y, -gx)
    single_results[g] = float(auc_g)
    if g in ["TPO", "DIO1", "TG"]:
        ax.plot(fpr_g, tpr_g, lw=1.5, color=colors_g[i], alpha=0.8, label=f"{g} alone AUC = {auc_g:.3f}")

# Random null (3 random gene sets size 8)
rng = np.random.default_rng(42)
null_aucs = []
for i in range(3):
    rgenes = rng.choice(expr.index, size=len(panel_present), replace=False)
    Xr = StandardScaler().fit_transform(expr.loc[rgenes, dm_h["expr_col"]].T.values)
    yp = cross_val_predict(clf, Xr, y, cv=skf, method="predict_proba")[:, 1]
    auc_r = roc_auc_score(y, yp)
    fpr_r, tpr_r, _ = roc_curve(y, yp)
    null_aucs.append(float(auc_r))
    if i == 0:
        ax.plot(fpr_r, tpr_r, lw=1.2, color="gray", alpha=0.6, ls="--", label=f"Random 8 genes AUC = {auc_r:.3f}")

ax.plot([0, 1], [0, 1], "--", color="gray", lw=1, alpha=0.5)
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
ax.set_title(f"Figure E25 — Multi-model ROC (cPTC vs FVPTC within DM, n={len(y)})\nPanel vs single-gene baselines vs random null")
ax.legend(loc="lower right", fontsize=10)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(FIG / "figE25_multimodel_ROC.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print(f"  8-gene AUC = {auc_8:.3f}")
print(f"  Single gene AUCs: {single_results}")
print(f"  Random null AUCs: {null_aucs}")
print("Saved figE25")

# ============================================================================
# 2. 8-gene gene-gene correlation matrix
# ============================================================================
print("\n=== 2. 8-gene gene-gene correlation matrix ===")
panel_expr = expr.loc[panel_present, dm_h["expr_col"]]
gene_corr = panel_expr.T.corr()
print(gene_corr.round(2))
gene_corr.to_csv(DATA_OUT / "panel_gene_gene_corr.tsv", sep="\t")

fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(gene_corr.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
ax.set_xticks(range(len(gene_corr))); ax.set_yticks(range(len(gene_corr)))
ax.set_xticklabels(gene_corr.columns, rotation=30, ha="right")
ax.set_yticklabels(gene_corr.index)
for i in range(len(gene_corr)):
    for j in range(len(gene_corr)):
        v = gene_corr.values[i, j]
        ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                color="white" if abs(v) > 0.6 else "black", fontsize=10)
plt.colorbar(im, ax=ax, label="Pearson r")
ax.set_title("Figure E26 — 8-gene panel gene-gene correlation\n(within DM cohort cPTC + FVPTC samples, n=125)")
fig.tight_layout()
fig.savefig(FIG / "figE26_panel_gene_corr.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved figE26")

# ============================================================================
# 3. Dot plot: 8-gene × 4 driver groups (bubble)
# ============================================================================
print("\n=== 3. Dot plot 8-gene × 4-way ===")
def fourway(row):
    if row.get("has_braf_v600e"): return "BRAF V600E"
    if row.get("has_ras_mut"): return "RAS hotspot"
    if row.get("v17_dark_cluster") == "DM1": return "DM1"
    if row.get("v17_dark_cluster") == "DM2": return "DM2"
    return "Other"
master["fourway"] = master.apply(fourway, axis=1)
groups = ["BRAF V600E", "RAS hotspot", "DM1", "DM2"]

# For each gene × group, compute mean expression + % cells expressing (>= 25th pctile)
dot_data = []
all_expr = expr.loc[panel_present]
threshold = all_expr.quantile(0.25, axis=1)
for g in panel_present:
    for grp in groups:
        cols = master.loc[master["fourway"] == grp, "expr_col"].tolist()
        if not cols: continue
        vals = expr.loc[g, cols]
        dot_data.append({"gene": g, "group": grp, "mean": float(vals.mean()),
                         "pct_above_q25": float((vals >= threshold[g]).mean() * 100),
                         "n": len(cols)})
dot_df = pd.DataFrame(dot_data)

fig, ax = plt.subplots(figsize=(9, 6))
gene_order = panel_present
group_order = groups
mean_norm = dot_df.copy()
# z-score per gene
mean_norm["mean_z"] = mean_norm.groupby("gene")["mean"].transform(lambda x: (x - x.mean()) / x.std())
for _, r in mean_norm.iterrows():
    yp = gene_order.index(r["gene"])
    xp = group_order.index(r["group"])
    color_val = r["mean_z"]
    size = r["pct_above_q25"] * 4
    ax.scatter(xp, yp, s=size, c=color_val, cmap="RdBu_r", vmin=-1.5, vmax=1.5,
               edgecolors="black", lw=0.5)
ax.set_xticks(range(len(group_order))); ax.set_xticklabels(group_order, rotation=15, ha="right")
ax.set_yticks(range(len(gene_order))); ax.set_yticklabels(gene_order)
ax.set_title(f"Figure E27 — 8-gene expression dot plot × 4 driver groups\nsize = % samples above panel-25th-pctile, color = row-z mean expression")
ax.grid(alpha=0.3)
# Color legend
sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(vmin=-1.5, vmax=1.5))
plt.colorbar(sm, ax=ax, label="row-z mean expr")
fig.tight_layout()
fig.savefig(FIG / "figE27_dotplot_4way.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved figE27")

# ============================================================================
# 4. Cox cluster × age interaction
# ============================================================================
print("\n=== 4. Cox cluster × age interaction ===")
master["age_yrs"] = pd.to_numeric(master["age"], errors="coerce")
dm_cox = master[master["dm_status"] & master["v17_dark_cluster"].isin(["DM1", "DM2"])].copy()
dm_cox["cluster_dm2"] = (dm_cox["v17_dark_cluster"] == "DM2").astype(int)
dm_cox["cluster_age_int"] = dm_cox["cluster_dm2"] * dm_cox["age_yrs"]
sub = dm_cox[["PFI", "PFI.time", "cluster_dm2", "age_yrs", "cluster_age_int"]].dropna()
print(f"n={len(sub)}, events={int(sub['PFI'].sum())}")
if int(sub["PFI"].sum()) >= 3:
    cph = CoxPHFitter(penalizer=0.05)
    try:
        cph.fit(sub.rename(columns={"PFI": "event", "PFI.time": "time"}),
                duration_col="time", event_col="event")
        s = cph.summary[["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].round(3)
        print(s.to_string())
        cox_int = s.to_dict()
    except Exception as e:
        print(f"Cox failed: {e}")
        cox_int = {"error": str(e)}

# ============================================================================
# 5. Confusion matrix + calibration at threshold 0.5
# ============================================================================
print("\n=== 5. Confusion matrix + calibration ===")
y_pred = (y_proba_8gene >= 0.5).astype(int)
cm = confusion_matrix(y, y_pred)
print(f"Confusion matrix:\n{cm}")
brier = brier_score_loss(y, y_proba_8gene)
print(f"Brier score: {brier:.3f}")

# Calibration plot — bin predictions, plot true rate
nbins = 8
bin_edges = np.linspace(0, 1, nbins + 1)
bins = np.digitize(y_proba_8gene, bin_edges) - 1
bins = np.clip(bins, 0, nbins - 1)
calib_data = []
for b in range(nbins):
    mask = bins == b
    if mask.sum() > 0:
        calib_data.append({
            "bin_center": (bin_edges[b] + bin_edges[b+1]) / 2,
            "mean_pred": float(y_proba_8gene[mask].mean()),
            "true_rate": float(y[mask].mean()),
            "n": int(mask.sum()),
        })
calib_df = pd.DataFrame(calib_data)
print(calib_df.round(3))

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
# Confusion matrix
ax = axes[0]
im = ax.imshow(cm, cmap="Blues", aspect="equal")
ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
ax.set_xticklabels(["Pred cPTC", "Pred FVPTC"])
ax.set_yticklabels(["True cPTC", "True FVPTC"])
for i in range(2):
    for j in range(2):
        ax.text(j, i, str(cm[i, j]), ha="center", va="center",
                color="white" if cm[i, j] > cm.max()*0.5 else "black", fontsize=18, fontweight="bold")
sens = cm[1, 1] / cm[1].sum() if cm[1].sum() > 0 else 0
spec = cm[0, 0] / cm[0].sum() if cm[0].sum() > 0 else 0
ppv = cm[1, 1] / cm[:, 1].sum() if cm[:, 1].sum() > 0 else 0
npv = cm[0, 0] / cm[:, 0].sum() if cm[:, 0].sum() > 0 else 0
ax.set_title(f"Figure E28A — Confusion matrix at threshold 0.5\nSens={sens:.2f}, Spec={spec:.2f}, PPV={ppv:.2f}, NPV={npv:.2f}")

# Calibration
ax = axes[1]
ax.plot([0, 1], [0, 1], "k--", lw=1, label="Perfect calibration")
ax.scatter(calib_df["mean_pred"], calib_df["true_rate"], s=calib_df["n"]*8, alpha=0.7, color="tab:blue", edgecolors="black")
for _, r in calib_df.iterrows():
    ax.annotate(f"n={r['n']}", (r["mean_pred"], r["true_rate"]), xytext=(5, 5), textcoords="offset points", fontsize=8)
ax.set_xlabel("Mean predicted probability")
ax.set_ylabel("Observed FVPTC rate")
ax.set_title(f"Figure E28B — Calibration plot (Brier = {brier:.3f})")
ax.set_xlim(-0.05, 1.05); ax.set_ylim(-0.05, 1.05)
ax.legend(loc="upper left")
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(FIG / "figE28_confusion_calibration.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved figE28")

# ============================================================================
# 6. Pathway × cluster mean heatmap (raw values, not z-scored)
# ============================================================================
print("\n=== 6. Pathway raw mean heatmap ===")
pm = pd.read_csv(DATA_OUT / "fourway_pathway_means.tsv", sep="\t", index_col=0)
fig, ax = plt.subplots(figsize=(7, 8))
im = ax.imshow(pm.values, cmap="viridis", aspect="auto")
ax.set_xticks(range(len(pm.columns))); ax.set_yticks(range(len(pm)))
ax.set_xticklabels(pm.columns, rotation=20, ha="right")
ax.set_yticklabels(pm.index)
for i in range(len(pm)):
    for j in range(len(pm.columns)):
        v = pm.values[i, j]
        ax.text(j, i, f"{v:.1f}", ha="center", va="center",
                color="white" if v < pm.values.mean() else "black", fontsize=9)
plt.colorbar(im, ax=ax, label="mean log2(TPM+1)")
ax.set_title("Figure E29 — Pathway raw mean heatmap × 4 driver groups\n(raw, not z-scored — absolute value comparison)")
fig.tight_layout()
fig.savefig(FIG / "figE29_pathway_raw_heatmap.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved figE29")

# Save final
final = {
    "ROC_8gene_AUC": float(auc_8),
    "single_gene_AUCs": single_results,
    "random_null_AUCs": null_aucs,
    "panel_gene_corr_min": float(gene_corr.values[~np.eye(len(gene_corr), dtype=bool)].min()),
    "panel_gene_corr_max": float(gene_corr.values[~np.eye(len(gene_corr), dtype=bool)].max()),
    "cox_cluster_age_interaction": cox_int if "cox_int" in dir() else None,
    "confusion_at_0.5": {"cm": cm.tolist(), "sens": float(sens), "spec": float(spec), "ppv": float(ppv), "npv": float(npv)},
    "brier_score": float(brier),
    "calibration_bins": calib_df.to_dict(orient="records"),
}
(DATA_OUT / "p9_omega_summary.json").write_text(json.dumps(final, indent=2, default=str))
print(f"\n=== ALL P9 OMEGA DONE ===")
