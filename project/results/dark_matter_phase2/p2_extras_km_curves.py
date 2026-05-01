"""Additional figures: Kaplan-Meier curves + DM1/DM2 stratification visuals + DICER1+ outcome."""
from pathlib import Path
import pandas as pd
import numpy as np
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
OUT = ROOT / "results/dark_matter_phase2"
FIG = OUT / "web/figures"

df = pd.read_csv(ROOT / "results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv", sep="\t")

# Add driver class from P2-D
master_v17 = pd.read_csv(ROOT / "results/v17/tables/sample_master_v17_full.tsv", sep="\t", low_memory=False)
master_v17["tcga_short"] = master_v17["sample_id"].str[:12]
master_v17_tcga = master_v17[(master_v17["dataset"] == "TCGA-THCA") & (master_v17["normal_vs_tumor"] == "tumor")]
df = df.merge(master_v17_tcga[["tcga_short", "fusion_classes"]].drop_duplicates("tcga_short"),
              on="tcga_short", how="left", suffixes=("", "_m"))

def classify(row):
    if row.get("has_braf_v600e"): return "BRAF V600E"
    if row.get("has_ras_mut"): return "RAS hotspot"
    da = str(row.get("driver_anchor_v17", ""))
    if da == "DICER1_EIF1AX_PPM1D": return "DICER1/EIF1AX"
    if row.get("tert_pos"): return "TERT-only"
    return "True driver-neg"
df["driver_class"] = df.apply(classify, axis=1)

# === Figure A: KM by driver class ===
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
ax = axes[0]
class_colors = {"BRAF V600E": "#e63946", "RAS hotspot": "#f4a261",
                "DICER1/EIF1AX": "#9d4edd", "TERT-only": "#1a1a2e",
                "True driver-neg": "#2a9d8f"}
for cls in ["BRAF V600E", "RAS hotspot", "True driver-neg", "DICER1/EIF1AX", "TERT-only"]:
    sub = df[(df["driver_class"] == cls) & df["PFI"].notna()]
    if len(sub) < 3: continue
    kmf = KaplanMeierFitter()
    kmf.fit(sub["PFI.time"], sub["PFI"], label=f"{cls} (n={len(sub)}, ev={int(sub['PFI'].sum())})")
    kmf.plot_survival_function(ax=ax, ci_show=False, color=class_colors[cls], lw=2)

# Logrank
classes_for_logrank = df[df["PFI"].notna() & df["driver_class"].isin(["BRAF V600E", "RAS hotspot", "True driver-neg", "DICER1/EIF1AX", "TERT-only"])].copy()
lr = multivariate_logrank_test(classes_for_logrank["PFI.time"], classes_for_logrank["driver_class"], classes_for_logrank["PFI"])
ax.set_title(f"Figure 8A — PFI by driver class (TCGA-THCA n=482)\nMultivariate logrank p = {lr.p_value:.4f}")
ax.set_xlabel("PFI time (days)")
ax.set_ylabel("Progression-Free Survival probability")
ax.set_ylim(0.5, 1.02)
ax.legend(loc="lower left", fontsize=8)
ax.grid(alpha=0.3)

# === Figure B: KM DM1 vs DM2 within Dark Matter ===
ax = axes[1]
dm = df[df["dm_status"] & df["v17_dark_cluster"].isin(["DM1", "DM2"]) & df["PFI"].notna()]
for cls, color in [("DM1", "#4361ee"), ("DM2", "#f4a261")]:
    sub = dm[dm["v17_dark_cluster"] == cls]
    kmf = KaplanMeierFitter()
    kmf.fit(sub["PFI.time"], sub["PFI"], label=f"{cls} (n={len(sub)}, ev={int(sub['PFI'].sum())})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color=color, lw=2)

dm1 = dm[dm["v17_dark_cluster"] == "DM1"]
dm2 = dm[dm["v17_dark_cluster"] == "DM2"]
lr2 = logrank_test(dm1["PFI.time"], dm2["PFI.time"], dm1["PFI"], dm2["PFI"])
ax.set_title(f"Figure 8B — PFI within Dark Matter: DM1 vs DM2\nlogrank p = {lr2.p_value:.3f} (TCGA underpowered)")
ax.set_xlabel("PFI time (days)")
ax.set_ylabel("Progression-Free Survival probability")
ax.set_ylim(0.7, 1.02)
ax.legend(loc="lower left", fontsize=9)
ax.grid(alpha=0.3)

fig.suptitle("Kaplan-Meier survival curves — TCGA-THCA PFI (Liu 2018 CDR)", fontsize=12, y=1.02)
fig.tight_layout()
fig.savefig(FIG / "fig8_KM_curves.png", dpi=200, bbox_inches="tight")
fig.savefig(FIG / "fig8_KM_curves.pdf", bbox_inches="tight")
print(f"Saved fig8_KM_curves.png/.pdf to {FIG}")
print(f"Multivariate logrank p (5 classes) = {lr.p_value:.4f}")
print(f"DM1 vs DM2 logrank p = {lr2.p_value:.3f}")

# === Figure C: PFI rate by class — bar chart with n labels ===
class_stats = []
for cls, sub in df.groupby("driver_class", observed=True):
    cs = sub.dropna(subset=["PFI"])
    n = len(cs); evt = int(cs["PFI"].sum())
    class_stats.append({"class": cls, "n": n, "events": evt, "rate": evt/n*100 if n else 0})
class_stats = pd.DataFrame(class_stats).sort_values("rate")
print("\nClass PFI rates:")
print(class_stats.to_string(index=False))
class_stats.to_csv(OUT / "p2_extras_class_pfi_rates.tsv", sep="\t", index=False)
