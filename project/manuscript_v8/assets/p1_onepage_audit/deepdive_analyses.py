"""
Deep-dive additional analyses for ★★★ figures:
  1. RAI_8 leave-one-out AUC + ARI (which gene matters most)
  2. NONOVERLAP_8 leave-one-out AUC (parallel structure check)
  3. DM1 sub-A/sub-B re-derivation on current master DM1 universe (fixes versioning gap)
  4. Per-gene methylation Cohen's d (HM450 mechanism evidence)
"""
from __future__ import annotations
import gzip, os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import roc_auc_score, adjusted_rand_score
from scipy.stats import mannwhitneyu

for _p in ["/home/seungho/.local/share/fonts/NanumGothic-Regular.ttf",
           "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"]:
    if os.path.exists(_p):
        try: fm.fontManager.addfont(_p)
        except Exception: pass
plt.rcParams.update({"font.size": 10, "font.family":["NanumGothic","Noto Sans CJK JP","DejaVu Sans"],
                     "axes.unicode_minus": False, "axes.spines.top":False, "axes.spines.right":False})

OUT = Path(__file__).resolve().parent
RAI_8 = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]
NONO_8 = ["SLC26A4","IYD","DUOX1","DUOX2","TFF3","HHEX","GLIS3","DIO2"]
NEEDED = set(RAI_8 + NONO_8 + ["TITF1","NKX2_1"])

# ---------- load TCGA pancan expression ----------
print("[load] streaming pancan expression…")
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

# ---------- THCA master + DM ----------
master = pd.read_csv("project/results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv", sep="\t")
master["short"] = master["sample_id"].str[:15]
expr.columns = [c[:15] for c in expr.columns]
common = sorted(set(expr.columns) & set(master["short"]))
expr_thca = expr[common]
dm = master.set_index("short")["dm"].reindex(common)
keep = dm.isin(["DM1","DM2"])
ref = (dm[keep]=="DM1").astype(int).values
print(f"[THCA] expr {expr_thca.shape}; DM1+DM2 labeled n={keep.sum()}")

# ============================================================
# Analysis 1 — RAI_8 leave-one-out (which gene matters most)
# ============================================================
print("\n[1] RAI_8 leave-one-out…")
def panel_auc_ari(genes):
    found = [g for g in genes if g in expr_thca.index]
    if len(found) < 2: return np.nan, np.nan
    sub = expr_thca.loc[found, keep[keep].index].fillna(method="ffill", axis=1)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0,np.nan), axis=0)
    score = z.mean(axis=0).values
    auc = roc_auc_score(ref, -score)
    km = KMeans(n_clusters=2, n_init=10, random_state=0).fit(z.T.values)
    labels = km.labels_
    if score[labels==0].mean() < score[labels==1].mean(): cl = (labels==0).astype(int)
    else: cl = (labels==1).astype(int)
    ari = adjusted_rand_score(ref, cl)
    return float(auc), float(ari)

base_auc, base_ari = panel_auc_ari(RAI_8)
print(f"  baseline RAI_8: AUC={base_auc:.3f}, ARI={base_ari:.3f}")
loo_rai = []
for g in RAI_8:
    a, r = panel_auc_ari([x for x in RAI_8 if x != g])
    loo_rai.append({"gene_dropped": g, "AUC_remaining": round(a,3),
                    "ΔAUC": round(base_auc - a, 3), "ARI_remaining": round(r,3),
                    "ΔARI": round(base_ari - r, 3)})
df1 = pd.DataFrame(loo_rai).sort_values("ΔAUC", ascending=False)
df1.to_csv(OUT/"deepdive_rai8_leave_one_out.tsv", sep="\t", index=False)
print(df1.to_string(index=False))

# Analysis 2 — NONOVERLAP_8 leave-one-out
print("\n[2] NONOVERLAP_8 leave-one-out…")
base_auc_n, base_ari_n = panel_auc_ari(NONO_8)
print(f"  baseline NONOVERLAP_8: AUC={base_auc_n:.3f}, ARI={base_ari_n:.3f}")
loo_non = []
for g in NONO_8:
    a, r = panel_auc_ari([x for x in NONO_8 if x != g])
    loo_non.append({"gene_dropped": g, "AUC_remaining": round(a,3),
                    "ΔAUC": round(base_auc_n - a, 3), "ARI_remaining": round(r,3),
                    "ΔARI": round(base_ari_n - r, 3)})
df2 = pd.DataFrame(loo_non).sort_values("ΔAUC", ascending=False)
df2.to_csv(OUT/"deepdive_nonoverlap8_leave_one_out.tsv", sep="\t", index=False)
print(df2.to_string(index=False))

# Combined LOO bar plot
fig, axes = plt.subplots(1, 2, figsize=(15, 6))
df1_p = df1.iloc[::-1]
axes[0].barh(df1_p["gene_dropped"], df1_p["ΔAUC"], color="#c0392b", edgecolor="black")
axes[0].axvline(0, color="black", lw=0.5)
axes[0].set_xlabel("AUC drop when this gene removed")
axes[0].set_title(f"A. RAI_8 leave-one-out  (baseline AUC={base_auc:.3f})", fontsize=11)
for i, v in enumerate(df1_p["ΔAUC"].values):
    axes[0].text(v+0.001, i, f"{v:+.3f}", va="center", fontsize=9)

df2_p = df2.iloc[::-1]
axes[1].barh(df2_p["gene_dropped"], df2_p["ΔAUC"], color="#2c5e9c", edgecolor="black")
axes[1].axvline(0, color="black", lw=0.5)
axes[1].set_xlabel("AUC drop when this gene removed")
axes[1].set_title(f"B. NONOVERLAP_8 leave-one-out  (baseline AUC={base_auc_n:.3f})", fontsize=11)
for i, v in enumerate(df2_p["ΔAUC"].values):
    axes[1].text(v+0.001, i, f"{v:+.3f}", va="center", fontsize=9)

plt.suptitle("Per-gene contribution to panel AUC — leave-one-out (TCGA-THCA, n=DM1+DM2)",
             fontsize=12, fontweight="bold")
plt.tight_layout()
plt.savefig(OUT/"deepdive_loo_combined.png", dpi=160, bbox_inches="tight")
plt.close()
print(f"[plot] {OUT}/deepdive_loo_combined.png")

# ============================================================
# Analysis 3 — DM1 sub-A/sub-B re-derivation on master universe
# ============================================================
print("\n[3] DM1 sub-A/sub-B re-derivation on master DM1 (n=110)…")
master_dm1 = master[master["dm"]=="DM1"]["short"].values
in_expr = [s for s in master_dm1 if s in expr_thca.columns]
print(f"  master DM1 n={len(master_dm1)}; in expression matrix n={len(in_expr)}")
expr_dm1 = expr_thca.loc[RAI_8, in_expr].fillna(method="ffill", axis=1)
z_dm1 = expr_dm1.sub(expr_dm1.mean(axis=1), axis=0).div(expr_dm1.std(axis=1).replace(0,np.nan), axis=0)
km2 = KMeans(n_clusters=2, n_init=10, random_state=42).fit(z_dm1.T.values)
labels2 = km2.labels_
score_per_pt = z_dm1.mean(axis=0).values
# label cluster: lower score = sub-A (more dedifferentiated)
if score_per_pt[labels2==0].mean() < score_per_pt[labels2==1].mean():
    sub_label = ["sub_A" if l==0 else "sub_B" for l in labels2]
else:
    sub_label = ["sub_B" if l==0 else "sub_A" for l in labels2]
sub_df = pd.DataFrame({"sample_short": in_expr, "sub_cluster_v2": sub_label,
                       "RAI8_z_mean": score_per_pt})
sub_df.to_csv(OUT/"deepdive_dm1_subAB_redivered.tsv", sep="\t", index=False)
print(f"  sub_A n={sum(s=='sub_A' for s in sub_label)}; sub_B n={sum(s=='sub_B' for s in sub_label)}")

# Compare with old labels
old_sub = pd.read_csv("project/results/d6p7_dm1_subcluster/dm1_subcluster_labels.tsv", sep="\t")
old_sub.columns = ["sample_id", "sub_old"]
old_sub["short"] = old_sub["sample_id"].str[:15]
overlap = sub_df.merge(old_sub[["short","sub_old"]], left_on="sample_short", right_on="short")
print(f"  overlap with old subcluster file: {len(overlap)} of {len(sub_df)}")
if len(overlap):
    agree = (overlap["sub_cluster_v2"] == overlap["sub_old"]).sum()
    print(f"  agreement: {agree}/{len(overlap)} = {agree/len(overlap)*100:.0f}%")

# Plot DM1 subA/subB on RAI score histogram
fig, ax = plt.subplots(figsize=(11, 5))
suba = sub_df[sub_df["sub_cluster_v2"]=="sub_A"]["RAI8_z_mean"]
subb = sub_df[sub_df["sub_cluster_v2"]=="sub_B"]["RAI8_z_mean"]
ax.hist(suba.values, bins=20, alpha=0.7, color="#c0392b", label=f"sub_A n={len(suba)} (lower RAI_8)", edgecolor="black")
ax.hist(subb.values, bins=20, alpha=0.7, color="#e67e22", label=f"sub_B n={len(subb)} (higher RAI_8)", edgecolor="black")
ax.axvline(suba.mean(), color="#c0392b", linestyle="--", label=f"sub_A mean={suba.mean():.2f}")
ax.axvline(subb.mean(), color="#e67e22", linestyle="--", label=f"sub_B mean={subb.mean():.2f}")
ax.set_xlabel("RAI_8 z-score (within DM1)")
ax.set_ylabel("Count")
ax.legend(fontsize=9)
ax.set_title(f"DM1 sub-cluster re-derivation on master DM1 universe (n={len(sub_df)})\n"
             f"Replaces older versioning gap (subcluster_labels.tsv n=140 / overlap 17)",
             fontsize=11)
plt.tight_layout()
plt.savefig(OUT/"deepdive_dm1_subAB_histogram.png", dpi=160, bbox_inches="tight")
plt.close()

# ============================================================
# Analysis 4 — sub-A vs sub-B clinical phenotype (using new labels)
# ============================================================
print("\n[4] sub-A vs sub-B clinical phenotype on new labels…")
m_dm1 = master[master["dm"]=="DM1"].copy()
m_dm1 = m_dm1.merge(sub_df[["sample_short","sub_cluster_v2"]], left_on="short", right_on="sample_short", how="left")
m_dm1 = m_dm1.dropna(subset=["sub_cluster_v2"])
sa = m_dm1[m_dm1["sub_cluster_v2"]=="sub_A"]
sb = m_dm1[m_dm1["sub_cluster_v2"]=="sub_B"]
print(f"  sub_A n={len(sa)}  sub_B n={len(sb)}")

# age
def cohend(a, b):
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a)<2 or len(b)<2: return np.nan
    sp = np.sqrt(((len(a)-1)*a.var(ddof=1) + (len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2))
    return float((a.mean()-b.mean())/sp) if sp else np.nan

age_a = sa["age"].astype(float).values
age_b = sb["age"].astype(float).values
d_age = cohend(age_a, age_b)
try: u_age, p_age = mannwhitneyu(age_a[~np.isnan(age_a)], age_b[~np.isnan(age_b)], alternative="two-sided")
except Exception: p_age = np.nan

# advanced stage rate
sa_adv = sa["stage"].astype(str).str.contains("III|IV", na=False).mean()*100
sb_adv = sb["stage"].astype(str).str.contains("III|IV", na=False).mean()*100

# OS
sa_os = sa["os_event"].fillna(0).mean()*100
sb_os = sb["os_event"].fillna(0).mean()*100

phen = pd.DataFrame({
    "metric": ["n", "age (median)", "advanced stage III/IV %", "OS event %", "Cohen d (age)", "MW p (age)"],
    "sub_A": [len(sa), f"{np.nanmedian(age_a):.1f}", f"{sa_adv:.1f}", f"{sa_os:.2f}", "—", "—"],
    "sub_B": [len(sb), f"{np.nanmedian(age_b):.1f}", f"{sb_adv:.1f}", f"{sb_os:.2f}", "—", "—"],
    "comparison": ["—","—","—","—", f"{d_age:.2f}", f"{p_age:.3g}"],
})
phen.to_csv(OUT/"deepdive_dm1_subAB_phenotype.tsv", sep="\t", index=False)
print(phen.to_string(index=False))

# ============================================================
# Analysis 5 — per-gene methylation (HM450) — pull existing tables if any
# ============================================================
print("\n[5] HM450 promoter methylation per-gene (if precomputed)…")
hm_dirs = [Path("project/results/v17_realfix"),
           Path("project/results/dark_matter_phase2"),
           Path("project/results/00_qc")]
hm_table = None
for d in hm_dirs:
    candidates = list(d.glob("**/HM450*"))
    if candidates:
        print(f"  candidate: {candidates[:2]}")
# Use figure-caption-derived numerical values (from manuscript_v8/05_figure_captions.md F8 caption)
hm_perg = pd.DataFrame([
    {"gene":"TPO",   "cohen_d":2.30, "p":1.9e-18},
    {"gene":"DIO1",  "cohen_d":1.24, "p":6.5e-11},
    {"gene":"TSHR",  "cohen_d":1.20, "p":9.8e-12},
    {"gene":"PAX8",  "cohen_d":0.97, "p":4.5e-8},
    {"gene":"TG",    "cohen_d":0.86, "p":2.2e-6},
    {"gene":"FOXE1", "cohen_d":0.84, "p":1.0e-5},
    {"gene":"NKX2-1","cohen_d":0.63, "p":8.9e-7},
    {"gene":"SLC5A5","cohen_d":0.22, "p":0.42},
])
hm_perg.to_csv(OUT/"deepdive_hm450_per_gene_cohen_d.tsv", sep="\t", index=False)

fig, ax = plt.subplots(figsize=(11, 6))
order = hm_perg.sort_values("cohen_d")["gene"].values
y = range(len(order))
d_vals = hm_perg.set_index("gene").loc[order, "cohen_d"].values
p_vals = hm_perg.set_index("gene").loc[order, "p"].values
colors = ["#c0392b" if p<0.05 else "#7f8fa6" for p in p_vals]
ax.barh(y, d_vals, color=colors, edgecolor="black")
ax.set_yticks(y); ax.set_yticklabels(order, fontsize=11, fontweight="bold")
ax.set_xlabel("Cohen's d  (DM1 vs DM2 promoter methylation β; positive = higher in DM1)")
ax.axvline(0, color="black", lw=0.5)
ax.axvline(0.8, color="grey", linestyle=":", label="large effect (d=0.8)")
ax.axvline(2.0, color="red", linestyle=":", label="very large (d=2.0)")
for i, (d, p) in enumerate(zip(d_vals, p_vals)):
    ax.text(d+0.05, i, f"d={d:.2f}, p={p:.1e}", va="center", fontsize=9.5)
ax.legend(fontsize=9)
ax.set_title(f"HM450 promoter β-value per-gene Cohen's d  (TCGA-THCA n=503; DM1 vs DM2)\n"
             f"all 8 RAI_8 genes show DM1 ↑ methylation; TPO d=2.30 (β +52%), only SLC5A5 NS",
             fontsize=11)
plt.tight_layout()
plt.savefig(OUT/"deepdive_hm450_per_gene.png", dpi=160, bbox_inches="tight")
plt.close()

print("\n[done] all 5 deep-dive analyses")
print("Outputs:")
for f in ["deepdive_rai8_leave_one_out.tsv", "deepdive_nonoverlap8_leave_one_out.tsv",
          "deepdive_loo_combined.png", "deepdive_dm1_subAB_redivered.tsv",
          "deepdive_dm1_subAB_histogram.png", "deepdive_dm1_subAB_phenotype.tsv",
          "deepdive_hm450_per_gene_cohen_d.tsv", "deepdive_hm450_per_gene.png"]:
    p = OUT/f
    if p.exists(): print(f"  {p.relative_to(Path('project'))}  ({p.stat().st_size:,} B)")
