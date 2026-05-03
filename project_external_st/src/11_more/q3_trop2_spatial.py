#!/usr/bin/env python3
"""Q3 — TROP2 (TACSTD2) spatial localization in ST data.
Closes therapeutic loop: DM1-high regions = TROP2-high regions = sacituzumab targets."""
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

OUT = Path("project_external_st/results/extra")

# Load all 28 slides + extract TACSTD2 per spot + correlate with DM1_like_score_resid
rows = []
# external 12
ext_reg = pd.read_csv("project_external_st/results/registry/dataset_sample_registry.tsv", sep="\t")
for _, r in ext_reg[ext_reg["status"]=="ok"].iterrows():
    h = Path(r["h5ad"])
    if not h.exists(): continue
    a = ad.read_h5ad(h)
    in_t = a.obs["in_tissue"].fillna(0).astype(int) == 1
    a = a[in_t].copy()
    if "TACSTD2" not in a.var_names:
        continue
    import scanpy as sc
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    trop2 = a[:, "TACSTD2"].X.toarray().ravel() if hasattr(a[:, "TACSTD2"].X, "toarray") else np.asarray(a[:, "TACSTD2"].X).ravel()
    spot = pd.read_csv(Path("project_external_st/results/scores") / f"{r['dataset']}_{r['sample_id']}_spot_scores.tsv.gz", sep="\t")
    n = min(len(trop2), len(spot))
    if n < 50: continue
    dm1 = spot["DM1_like_score_resid"].iloc[:n].to_numpy()
    epi = spot["Epithelial_score_resid"].iloc[:n].to_numpy()
    trop2 = trop2[:n]
    ok = ~(np.isnan(dm1) | np.isnan(trop2))
    if ok.sum() < 50: continue
    rho_dm1, p_dm1 = spearmanr(trop2[ok], dm1[ok])
    rho_epi, p_epi = spearmanr(trop2[ok], epi[ok])
    rows.append({"sample_id": r["sample_id"], "dataset": r["dataset"], "condition": r["condition_inferred"],
                 "n_spots": int(ok.sum()),
                 "TROP2_mean": float(trop2[ok].mean()),
                 "rho_TROP2_DM1": rho_dm1, "p_TROP2_DM1": p_dm1,
                 "rho_TROP2_Epi": rho_epi, "p_TROP2_Epi": p_epi})

# GSE250521
g521_reg = pd.read_csv("project/data/processed/GSE250521/sample_metadata.tsv", sep="\t")
for _, r in g521_reg[g521_reg["status"]=="ok"].iterrows():
    h = Path(r["h5ad"])
    if not h.exists(): continue
    a = ad.read_h5ad(h)
    in_t = a.obs["in_tissue"].fillna(0).astype(int) == 1
    a = a[in_t].copy()
    if "TACSTD2" not in a.var_names:
        continue
    import scanpy as sc
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    trop2 = a[:, "TACSTD2"].X.toarray().ravel() if hasattr(a[:, "TACSTD2"].X, "toarray") else np.asarray(a[:, "TACSTD2"].X).ravel()
    # GSE250521 raw spot scores
    spot = pd.read_csv("project/results/01_spatial_score/all_spots_scored.tsv.gz", sep="\t")
    spot = spot[spot["sample_id"] == r["sample_id"]]
    n = min(len(trop2), len(spot))
    if n < 50: continue
    dm1 = spot["DM1_like_score"].iloc[:n].to_numpy()
    epi = spot["Epithelial_score"].iloc[:n].to_numpy()
    trop2 = trop2[:n]
    ok = ~(np.isnan(dm1) | np.isnan(trop2))
    if ok.sum() < 50: continue
    rho_dm1, p_dm1 = spearmanr(trop2[ok], dm1[ok])
    rho_epi, p_epi = spearmanr(trop2[ok], epi[ok])
    rows.append({"sample_id": r["sample_id"], "dataset": "GSE250521", "condition": r["stage_inferred"],
                 "n_spots": int(ok.sum()),
                 "TROP2_mean": float(trop2[ok].mean()),
                 "rho_TROP2_DM1": rho_dm1, "p_TROP2_DM1": p_dm1,
                 "rho_TROP2_Epi": rho_epi, "p_TROP2_Epi": p_epi})

df = pd.DataFrame(rows)
df.to_csv(OUT / "q3_trop2_spatial.tsv", sep="\t", index=False)
print(f"  {len(df)} slides analyzed")
print("\n=== TROP2 ↔ DM1_like spatial correlation per slide ===")
print(df.to_string(index=False))
print(f"\n  mean ρ(TROP2, DM1) = {df['rho_TROP2_DM1'].mean():.3f}")
print(f"  slides with positive (TROP2-high = DM1-high colocalize) ρ: {(df['rho_TROP2_DM1'] > 0).sum()} / {len(df)}")

# Figure
fig, axes = plt.subplots(1, 2, figsize=(15, 6))
import seaborn as sns
ax = axes[0]
order = ["PT","PTC","LPTC","ATC","CONTROL","HT","GD","PTC_HT"]
df["condition"] = df["condition"].fillna("unknown")
df_plot = df[df["condition"].isin(order)].copy()
sns.boxplot(data=df_plot, x="condition", y="rho_TROP2_DM1", order=order, ax=ax,
            palette="RdBu_r")
sns.stripplot(data=df_plot, x="condition", y="rho_TROP2_DM1", order=order, ax=ax,
              color="black", size=8, alpha=0.7)
ax.axhline(0, color="grey", lw=0.5, ls=":")
ax.set_title(f"Q3.A — TROP2 ↔ DM1 spatial correlation per slide ({len(df_plot)} slides)\n"
             f"Positive ρ = TROP2-high spots colocalize with DM1-high spots",
             fontsize=11)
ax.set_xlabel(""); ax.set_ylabel("Spearman ρ (TROP2 vs DM1) per slide")

ax = axes[1]
ax.scatter(df["TROP2_mean"], df["rho_TROP2_DM1"], s=80,
           c=[{"PT":"#3C6B4F","PTC":"#B8893C","LPTC":"#A04451","ATC":"#962E2E",
                "CONTROL":"#4C72B0","HT":"#DD8452","GD":"#55A467","PTC_HT":"#C44E52"}.get(c,"grey")
              for c in df["condition"]],
           edgecolor="black")
ax.axhline(0, color="grey", lw=0.5, ls=":")
ax.set_xlabel("Mean TROP2 (log-norm) per slide")
ax.set_ylabel("ρ(TROP2, DM1) per slide")
ax.set_title("Q3.B — TROP2 abundance vs colocalization", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "q3_trop2_spatial.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"→ {OUT / 'q3_trop2_spatial.png'}")
