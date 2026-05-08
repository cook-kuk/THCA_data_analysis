"""
Quick gene-panel combination scan — alternative compact RAI-lineage readouts.

Tests panels of various sizes/compositions vs the discovery DM1/DM2 reference partition
in TCGA-THCA. AUC + ARI for each combination.

Outputs:
  - panel_combos.tsv  (full table)
  - fig_panel_combos.png  (forest plot of AUC + ARI)
"""
from __future__ import annotations
import gzip
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import roc_auc_score, adjusted_rand_score
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent

# ---------- gene definitions ----------
RAI_8        = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]
NONOVERLAP_8 = ["SLC26A4","IYD","DUOX1","DUOX2","TFF3","HHEX","GLIS3","DIO2"]
TF_4         = ["FOXE1","NKX2-1","PAX8","HHEX"]
RAI_5_eff    = ["SLC5A5","TPO","TG","TSHR","DIO1"]   # effectors only
TDS_16       = sorted(set(RAI_8 + NONOVERLAP_8))      # full TDS-core
TIERA_extra  = ["DIO2","SLC5A8","THRA","THRB"]        # additional TDS members
MECHANISM    = ["STAT3","FOSL1","JUNB","DNMT1","DNMT3B"]

PANELS = {
    # baselines
    "RAI_8 (canonical)":              RAI_8,
    "TDS-16 (full)":                  TDS_16,
    "NONOVERLAP_8":                   NONOVERLAP_8,
    "TF_4_backbone":                  TF_4,
    # variants
    "RAI_5_effectors_only":           RAI_5_eff,
    "TF_3_within_RAI8":               ["FOXE1","NKX2-1","PAX8"],
    "RAI_8 + HHEX (9)":               RAI_8 + ["HHEX"],
    "RAI_8 + TF_4 (12 unique)":       sorted(set(RAI_8 + TF_4)),
    "RAI_8 + NONOVERLAP_8 (16)":      RAI_8 + NONOVERLAP_8,
    "TF_4 + RAI_5 (9)":               TF_4 + RAI_5_eff,
    "TF_3 + iodide (NIS+TPO+TG+SLC26A4+IYD)": ["FOXE1","NKX2-1","PAX8","SLC5A5","TPO","TG","SLC26A4","IYD"],
    "Iodide-handling 6":              ["SLC5A5","SLC26A4","TPO","DUOX1","DUOX2","IYD"],
    "Lineage TF + Effector minimal 6":["FOXE1","NKX2-1","PAX8","TG","TPO","DIO1"],
    "Lineage TF + Effector minimal 4":["FOXE1","NKX2-1","TG","TPO"],
    "TIERA TDS-extended 12":          sorted(set(RAI_8 + ["DIO2","SLC5A8","THRA","THRB"])),
    "Mechanism arm 5 (STAT3+AP1+DNMT)":MECHANISM,
}

# ---------- load TCGA-THCA expression ----------
NEEDED = set()
for g_list in PANELS.values():
    NEEDED |= set(g_list)
NEEDED |= {"TITF1","NKX2_1"}

print(f"[load] streaming pancan for {len(NEEDED)} genes...")
rows = {}
with gzip.open("project/data/raw/TCGA_pancan/pancan_geneExp.gz","rt") as f:
    header = next(f).strip().split("\t")
    samples = header[1:]
    for line in f:
        parts = line.rstrip("\n").split("\t")
        sym = parts[0].strip()
        if sym in NEEDED:
            rows[sym] = [float(x) if x not in ("","NA","NaN") else np.nan for x in parts[1:]]
expr = pd.DataFrame(rows, index=samples).T
if "NKX2-1" not in expr.index:
    if "NKX2_1" in expr.index: expr = expr.rename(index={"NKX2_1":"NKX2-1"})
    elif "TITF1" in expr.index: expr = expr.rename(index={"TITF1":"NKX2-1"})
expr = expr[~expr.index.duplicated(keep="first")]

# ---------- THCA only + DM call from master ----------
master = pd.read_csv("project/results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv", sep="\t")
master["short"] = master["sample_id"].str[:15]
expr.columns = [c[:15] for c in expr.columns]
common = sorted(set(expr.columns) & set(master["short"]))
expr_thca = expr[common]
dm_lookup = master.set_index("short")["dm"]
dm_thca = dm_lookup.loc[common]
keep = dm_thca.isin(["DM1","DM2"])
print(f"[THCA] expr {expr_thca.shape}; DM1/DM2 labeled n={keep.sum()}")

# ---------- score each panel ----------
results = []
ref = (dm_thca[keep]=="DM1").astype(int).values  # 1 = DM1, 0 = DM2
for name, gene_list in PANELS.items():
    found = [g for g in gene_list if g in expr_thca.index]
    if len(found) < 2:
        results.append({"panel": name, "n_genes": len(found), "n_found": len(found),
                        "AUC": np.nan, "ARI": np.nan, "missing": ",".join(set(gene_list)-set(found))})
        continue
    sub = expr_thca.loc[found, keep[keep].index].fillna(method="ffill", axis=1)
    sub_z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    score = sub_z.mean(axis=0)  # mean z per sample → high = differentiated; we want low for DM1
    # AUC: predict DM1 (1=DM1). Score lower = DM1, so we use -score.
    auc = roc_auc_score(ref, -score.values)
    # ARI: KMeans k=2 on z-vector × samples
    km = KMeans(n_clusters=2, n_init=10, random_state=0).fit(sub_z.T.values)
    labels = km.labels_
    # align labels: which cluster has lower mean score → assign as DM1=1
    if score.values[labels==0].mean() < score.values[labels==1].mean():
        cl = (labels==0).astype(int)
    else:
        cl = (labels==1).astype(int)
    ari = adjusted_rand_score(ref, cl)
    results.append({"panel": name, "n_genes": len(gene_list), "n_found": len(found),
                    "AUC": round(float(auc),3), "ARI": round(float(ari),3),
                    "missing": ",".join(set(gene_list)-set(found))})

res_df = pd.DataFrame(results).sort_values("AUC", ascending=False)
res_df.to_csv(OUT/"panel_combos.tsv", sep="\t", index=False)
print("\n[results]")
print(res_df.to_string(index=False))

# ---------- forest plot ----------
fig, axes = plt.subplots(1, 2, figsize=(14, 7))
df_plot = res_df.dropna(subset=["AUC","ARI"]).iloc[::-1]  # reverse for plot
y = np.arange(len(df_plot))
axes[0].barh(y, df_plot["AUC"], color="#2c5e9c", edgecolor="black")
axes[0].set_yticks(y); axes[0].set_yticklabels(df_plot["panel"], fontsize=9)
axes[0].set_xlabel("AUC (DM1 vs DM2 classification, supervised)")
axes[0].axvline(0.5, color="grey", linestyle=":")
axes[0].axvline(0.962, color="red", linestyle="--", lw=1.0, label="RAI_8 anchor 0.962")
axes[0].legend(fontsize=8); axes[0].set_xlim(0.5, 1.0)
axes[0].set_title("AUC per panel (deployment metric)")
for i, v in enumerate(df_plot["AUC"]):
    axes[0].text(v+0.005, i, f"{v}", va="center", fontsize=8)

axes[1].barh(y, df_plot["ARI"], color="#c0392b", edgecolor="black")
axes[1].set_yticks(y); axes[1].set_yticklabels([])
axes[1].set_xlabel("ARI (unsupervised KMeans vs DM1/DM2 reference)")
axes[1].axvline(0.49, color="red", linestyle="--", lw=1.0, label="RAI_8 anchor 0.49")
axes[1].axvline(0.90, color="grey", linestyle=":", lw=1.0, label="TIERA67 anchor 0.90")
axes[1].legend(fontsize=8); axes[1].set_xlim(0, 1.0)
axes[1].set_title("ARI per panel (discovery metric)")
for i, v in enumerate(df_plot["ARI"]):
    axes[1].text(v+0.01, i, f"{v}", va="center", fontsize=8)

plt.suptitle("Gene-panel combination scan — AUC (deployment) vs ARI (discovery)\nReference: TCGA-THCA DM1 vs DM2 master partition",
             fontsize=12)
plt.tight_layout()
plt.savefig(OUT/"fig_panel_combos.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"[plot] {OUT}/fig_panel_combos.png")
