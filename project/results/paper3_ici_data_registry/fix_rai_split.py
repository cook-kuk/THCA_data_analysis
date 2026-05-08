#!/usr/bin/env python3
"""Fix the RAI before/after split bug. Re-compute the GSE151179 RAI module shift."""
import pandas as pd, numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/paper3_ici_track_b_lite")
FIG = OUT / "figures_png"
REG = Path("/home/seungho/personal/THCA_data_analysis/project/reports/paper3_ici")

mods = pd.read_csv(REG/"paper3_ici_module_gene_list.tsv", sep="\t")
MODULES = sorted(set(mods.module))

# Reload GSE151179 score (with metadata) — file has cohort + rai_collection cols
z = pd.read_csv(OUT/"scores_per_cohort/GSE151179_module_scores.tsv", sep="\t", index_col=0)
print(f"loaded n={len(z)}")
print("rai_collection raw values:", z.rai_collection.value_counts().to_dict())

# Parse the value AFTER ':' to get true Before/After token
def rai_status(s):
    s = (s or "").lower()
    # Format: "collection before/after rai: Before" or "collection before/after rai: After"
    parts = s.split(":")
    if len(parts) < 2: return "unknown"
    last = parts[-1].strip()
    if last == "before": return "Before"
    if last == "after":  return "After"
    return "unknown"

z["rai"] = z.rai_collection.apply(rai_status)
print("RAI parsed:", z.rai.value_counts().to_dict())

pre = z[z.rai=="Before"]
post = z[z.rai=="After"]
print(f"\npre n={len(pre)}  post n={len(post)}")

rows=[]
for mod in MODULES:
    a = pre[mod].dropna().values
    b = post[mod].dropna().values
    if len(a)<3 or len(b)<3: continue
    d = (b.mean()-a.mean()) / np.sqrt((a.std(ddof=1)**2 + b.std(ddof=1)**2)/2 + 1e-12)
    t,p = stats.ttest_ind(b, a, equal_var=False)
    rows.append({
        "module":mod,"n_before":len(a),"n_after":len(b),
        "mean_before":float(a.mean()),"mean_after":float(b.mean()),
        "cohens_d_after_minus_before":float(d),"t":float(t),"p":float(p),
        "sign":"+" if d>0 else "-"
    })

rai_df = pd.DataFrame(rows)
rai_df.to_csv(OUT/"gse151179_rai_after_vs_before.tsv", sep="\t", index=False)
print("\nGSE151179 RAI shift (post-RAI vs pre-RAI):")
print(rai_df.round(3).to_string())

# Also stratify by tissue type — primary tumor vs lymph-node-met
def tissue_class(s):
    s = (s or "").lower()
    if "primary" in s: return "primary_tumor"
    if "non-neoplastic" in s or "normal" in s: return "non_neoplastic"
    if "lymph node" in s or "metastasis" in s: return "ln_met"
    return "unknown"
z["tissue_class"] = z.tissue_type.apply(tissue_class)
print("\ntissue class:", z.tissue_class.value_counts().to_dict())

# Primary vs LN_met (post-RAI 9 + 4 + 4 = 17 LN; primary 17)
prim = z[z.tissue_class=="primary_tumor"]
ln = z[z.tissue_class=="ln_met"]
nrm = z[z.tissue_class=="non_neoplastic"]
print(f"primary {len(prim)}, ln_met {len(ln)}, non_neoplastic {len(nrm)}")

ti_rows=[]
for mod in MODULES:
    a = prim[mod].dropna().values; b = ln[mod].dropna().values
    if len(a)<3 or len(b)<3: continue
    d = (b.mean()-a.mean())/np.sqrt((a.std(ddof=1)**2+b.std(ddof=1)**2)/2+1e-12)
    t,p = stats.ttest_ind(b,a,equal_var=False)
    ti_rows.append({"module":mod,"d_LN_minus_primary":float(d),"p":float(p),"sign":"+" if d>0 else "-"})
ti_df = pd.DataFrame(ti_rows)
ti_df.to_csv(OUT/"gse151179_lnmet_vs_primary.tsv", sep="\t", index=False)
print("\nLN met (post-RAI predominantly) vs primary tumor:")
print(ti_df.round(3).to_string())

# Update fig7
plt.figure(figsize=(8,5))
mods_l = list(MODULES)
rai_d = rai_df.set_index("module").reindex(mods_l)
plt.barh(mods_l, rai_d.cohens_d_after_minus_before,
         color=["red" if x>0 else "steelblue" for x in rai_d.cohens_d_after_minus_before])
plt.axvline(0, color="k", lw=0.8)
plt.xlabel("Cohen's d (post-RAI − pre-RAI)")
plt.title(f"GSE151179 — RAI-refractory progression: post-RAI vs pre-RAI module shift\n"
          f"(n_before={int(rai_df.n_before.iloc[0])}, n_after={int(rai_df.n_after.iloc[0])})")
plt.tight_layout()
plt.savefig(FIG/"fig7_rai_after_vs_before.png", dpi=150)
plt.close()
print("\nfig7 updated")

# Add fig9: LN_met vs primary
plt.figure(figsize=(8,5))
ti_d = ti_df.set_index("module").reindex(mods_l)
plt.barh(mods_l, ti_d.d_LN_minus_primary,
         color=["red" if x>0 else "steelblue" for x in ti_d.d_LN_minus_primary])
plt.axvline(0, color="k", lw=0.8)
plt.xlabel("Cohen's d (LN-met − primary tumor)")
plt.title(f"GSE151179 — Lymph-node metastases vs primary tumor module shift\n"
          f"(LN-met n={len(ln)}, primary n={len(prim)})")
plt.tight_layout()
plt.savefig(FIG/"fig9_lnmet_vs_primary.png", dpi=150)
plt.close()
print("fig9 saved")
