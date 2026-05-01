"""P2-D Full driver map (TCGA quick win) — 7-class hierarchy × cluster × phenotype × Cox."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
OUT = ROOT / "results/dark_matter_phase2"
OUT.mkdir(exist_ok=True)

# Reuse Phase 1 master table
df = pd.read_csv(ROOT / "results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv", sep="\t")
print(f"Loaded master: {len(df)} rows, columns include: {list(df.columns)[:15]}")

# Load fusion table from earlier inventory
fusion = pd.read_csv(ROOT / "metadata/v3_fusion_anchor_tcga.tsv", sep="\t")
print(f"Fusion table: {fusion.shape}, cols: {list(fusion.columns)[:8]}")
fusion["tcga_short"] = fusion["sample_id"].str[:12] if "sample_id" in fusion.columns else fusion["tcga12"]

# Hierarchical 7-class classification (mutually exclusive, priority order)
def classify(row):
    if row.get("has_braf_v600e"):
        return "Class1_BRAF_V600E"
    if row.get("has_ras_mut"):
        return "Class2_RAS_hotspot"
    # Need fusion data merged
    fc = str(row.get("fusion_classes", ""))
    if "fusion" in fc.lower() and "other" not in fc.lower():
        return "Class3_Fusion"
    da = str(row.get("driver_anchor_v17", ""))
    if da == "DICER1_EIF1AX_PPM1D":
        return "Class4_DICER1_EIF1AX"
    if row.get("tert_pos"):
        return "Class5_TERT_only"
    return "Class6_True_driver_neg"

# Merge fusion classes from sample_master
master_v17 = pd.read_csv(ROOT / "results/v17/tables/sample_master_v17_full.tsv", sep="\t", low_memory=False)
master_v17["tcga_short"] = master_v17["sample_id"].str[:12]
master_v17_tcga = master_v17[(master_v17["dataset"] == "TCGA-THCA") & (master_v17["normal_vs_tumor"] == "tumor")]
df = df.merge(master_v17_tcga[["tcga_short", "fusion_classes"]].drop_duplicates("tcga_short"),
              on="tcga_short", how="left", suffixes=("", "_m"))

df["driver_class"] = df.apply(classify, axis=1)
print(f"\n=== Driver class distribution (TCGA n={len(df)}) ===")
class_counts = df["driver_class"].value_counts()
print(class_counts)

# Class × cluster
ct = pd.crosstab(df["driver_class"], df["v17_dark_cluster"], dropna=False)
print(f"\n=== Driver class × 8-gene cluster ===\n{ct}")

# Per-class clinical phenotype
print("\n=== Per-class phenotype ===")
phen_summary = []
for cls, sub in df.groupby("driver_class", observed=True):
    s = {
        "class": cls, "n": len(sub),
        "DM1_n": int((sub["v17_dark_cluster"] == "DM1").sum()),
        "DM2_n": int((sub["v17_dark_cluster"] == "DM2").sum()),
        "age_mean": float(pd.to_numeric(sub["age"], errors="coerce").mean()),
        "PFI_events": int(sub["PFI"].sum()) if "PFI" in sub.columns else None,
        "PFI_n": int(sub["PFI"].notna().sum()) if "PFI" in sub.columns else None,
        "PFI_event_rate": float(sub["PFI"].sum() / sub["PFI"].notna().sum()) if "PFI" in sub.columns and sub["PFI"].notna().sum() > 0 else None,
    }
    phen_summary.append(s)
phen_df = pd.DataFrame(phen_summary)
print(phen_df.to_string(index=False))

# PFI Cox with class as covariate (BRAF V600E reference)
cph_df = df[["PFI", "PFI.time", "driver_class"]].dropna()
cph_df["class_idx"] = cph_df["driver_class"].astype("category").cat.codes
ref_class = "Class1_BRAF_V600E"
cph_df["is_BRAF"] = (cph_df["driver_class"] == ref_class).astype(int)
print(f"\n=== Cox HR (PFI, by class vs BRAF V600E) ===")
for cls in cph_df["driver_class"].unique():
    if cls == ref_class: continue
    sub = cph_df[cph_df["driver_class"].isin([ref_class, cls])].copy()
    sub["is_other"] = (sub["driver_class"] == cls).astype(int)
    if int(sub["PFI"].sum()) < 3 or sub["is_other"].nunique() < 2:
        continue
    try:
        cph = CoxPHFitter(penalizer=0.01)
        cph.fit(sub.rename(columns={"PFI": "event", "PFI.time": "time"})[["time", "event", "is_other"]],
                duration_col="time", event_col="event")
        s_row = cph.summary.iloc[0]
        print(f"  {cls} vs BRAF: HR={s_row['exp(coef)']:.2f} [{s_row['exp(coef) lower 95%']:.2f}-{s_row['exp(coef) upper 95%']:.2f}], p={s_row['p']:.3f}, n={len(sub)}, events={int(sub['PFI'].sum())}")
    except Exception as e:
        print(f"  {cls}: {e}")

# Save
phen_df.to_csv(OUT / "p2d_driver_class_phenotype.tsv", sep="\t", index=False)
ct.to_csv(OUT / "p2d_driver_class_x_cluster.tsv", sep="\t")
df[["tcga_short", "driver_class", "v17_dark_cluster", "PFI", "PFI.time"]].to_csv(
    OUT / "p2d_per_sample_classification.tsv", sep="\t", index=False)

# Figure: stacked bar — class × cluster
fig, ax = plt.subplots(figsize=(10, 5))
class_order = sorted(df["driver_class"].dropna().unique())
dm1_pct = []; dm2_pct = []; un_pct = []
for c in class_order:
    sub = df[df["driver_class"] == c]
    n = len(sub)
    dm1_pct.append((sub["v17_dark_cluster"] == "DM1").sum() / n * 100 if n else 0)
    dm2_pct.append((sub["v17_dark_cluster"] == "DM2").sum() / n * 100 if n else 0)
    un_pct.append(sub["v17_dark_cluster"].isna().sum() / n * 100 if n else 0)
x = np.arange(len(class_order))
ax.bar(x, dm1_pct, label="DM1 (cPTC-like)", color="tab:blue", alpha=0.8)
ax.bar(x, dm2_pct, bottom=dm1_pct, label="DM2 (FVPTC-like)", color="tab:orange", alpha=0.8)
ax.bar(x, un_pct, bottom=np.array(dm1_pct)+np.array(dm2_pct), label="No cluster", color="lightgray", alpha=0.5)
ax.set_xticks(x)
ax.set_xticklabels([c.replace("Class", "") for c in class_order], rotation=20, ha="right", fontsize=9)
ax.set_ylabel("% of class")
ax.set_title(f"Figure S1 — Driver class × 8-gene cluster (TCGA-THCA, n={len(df)})\n"
             f"Class N: {dict(class_counts)}")
ax.legend(loc="upper right", fontsize=9)
ax.grid(alpha=0.2, axis="y")
fig.tight_layout()
fig.savefig(OUT / "fig_S1_driver_class_x_cluster.png", dpi=200, bbox_inches="tight")
fig.savefig(OUT / "fig_S1_driver_class_x_cluster.pdf", bbox_inches="tight")
print(f"\nSaved figure S1 + tables.")

summary = {
    "n_total": int(len(df)),
    "class_counts": {k: int(v) for k, v in class_counts.to_dict().items()},
    "DM1_DM2_by_class": ct.to_dict(),
    "phenotype": phen_df.to_dict(orient="records"),
}
(OUT / "p2d_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print("Saved p2d_summary.json")
