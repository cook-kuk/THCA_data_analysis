"""
8-gene panel reduction feasibility analysis.

Combines:
  - TCGA HM450 promoter β (n=503)
  - Correlation structure (Within-TF ρ=0.68 / Within-Effector ρ=0.14, PC1 = 46.2%)
  - TSO500 v2 panel coverage (TSHR, PAX8, NKX2-1 directly; PAX8-PPARG via PPARG fusion)
  - Per-gene Cohen's d (DM1 - DM2)

Evaluates 4 reduction strategies × performance:
  A. 7-gene LOO  (knock-out one)
  B. Top-N by Cohen's d
  C. TSO500-aligned (3-gene baseline + add-ons)
  D. Hand-picked compact panels (clinical interpretability)

For each candidate subset:
  - mean panel β (the actual score the assay would output)
  - Cohen's d (DM1 vs DM2)
  - AUC for DM1 classification (LogReg on subset β)
  - Spearman r vs full 8-gene score
  - TSO500-native compatibility (0/3 or 1/3 etc.)
  - IHC routine availability

Outputs:
  - reduction_summary.tsv         — all candidate panels with metrics
  - reduction_summary.json        — for web component
  - fig_reduction_tradeoff.png    — panel size vs metric curve
  - fig_reduction_recommend.png   — recommended panels visual comparison
"""
import os, json, itertools, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
from scipy.stats import spearmanr

ROOT = "/home/seungho/personal/THCA_data_analysis/project"
OUT  = f"{ROOT}/dm1_story_web/public/figures"
os.makedirs(OUT, exist_ok=True)

GENES = ["DIO1","FOXE1","NKX2-1","PAX8","SLC5A5","TG","TPO","TSHR"]
GROUP = {"PAX8":"TF","NKX2-1":"TF","FOXE1":"TF",
         "TG":"E","TPO":"E","TSHR":"E","SLC5A5":"E","DIO1":"E"}
TSO500 = {"TSHR": True, "PAX8": True, "NKX2-1": True}   # 3 / 8 on TSO500 v2 DNA panel
IHC_TIER = {"TG":3, "PAX8":3, "NKX2-1":3,                # routine
            "TPO":2, "SLC5A5":2, "TSHR":2,                # validated
            "DIO1":1, "FOXE1":1}                          # research
# Per-gene Cohen's d (HM450 promoter β, audit-locked Round 5)
GENE_D = {"TPO":2.30, "DIO1":1.24, "TSHR":1.20, "PAX8":0.97,
          "TG":0.86, "FOXE1":0.84, "NKX2-1":0.63, "SLC5A5":0.22}

# ─── Load TCGA HM450 + DM labels ───
hm = pd.read_csv(f"{ROOT}/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv", sep="\t")
hm["sample_short"] = hm["sample_short"]
dm = pd.read_csv(f"{ROOT}/results/manuscript_v8_nc_main/master_tcga.tsv", sep="\t")[["patient_id","DM"]]
df = hm.merge(dm, left_on="sample_short", right_on="patient_id", how="inner").dropna(subset=["DM"] + GENES)
df["DM1"] = (df["DM"] == "DM1").astype(int)
print(f"Cohort  n = {len(df)}  ·  DM1 = {df['DM1'].sum()}  ·  DM2 = {(df['DM1']==0).sum()}")

# Full 8-gene score (reference)
df["full_score"] = df[GENES].mean(axis=1)

def panel_metrics(subset):
    """Compute Cohen's d, AUC, Spearman r vs full, and TSO500 coverage."""
    score = df[subset].mean(axis=1)
    g1 = score[df["DM1"]==1].values
    g0 = score[df["DM1"]==0].values
    pooled_sd = np.sqrt((g1.var(ddof=1) + g0.var(ddof=1)) / 2)
    d = (g1.mean() - g0.mean()) / pooled_sd if pooled_sd > 0 else 0
    # AUC via 5-fold CV LogReg using subset β
    X = StandardScaler().fit_transform(df[subset].values)
    y = df["DM1"].values
    auc = cross_val_score(LogisticRegression(max_iter=400), X, y, cv=5, scoring="roc_auc").mean()
    # vs full 8-gene
    r = spearmanr(score, df["full_score"]).statistic
    tso_n = sum(1 for g in subset if TSO500.get(g, False))
    ihc_avail = sum(1 for g in subset if IHC_TIER[g] == 3)
    return {
        "n_genes": len(subset),
        "genes": " · ".join(subset),
        "cohens_d": d,
        "auc_5fcv": auc,
        "r_vs_full": r,
        "tso500_native": f"{tso_n} / {len(subset)}",
        "ihc_routine": f"{ihc_avail} / {len(subset)}",
        "tf_count": sum(1 for g in subset if GROUP[g]=="TF"),
        "ef_count": sum(1 for g in subset if GROUP[g]=="E")
    }

# Reference: full 8-gene panel
rows = []
ref = panel_metrics(GENES)
ref.update({"label": "FULL 8-gene panel (reference)", "category": "baseline"})
rows.append(ref)

# A. n-1 LOO (8 panels)
for g in GENES:
    sub = [x for x in GENES if x != g]
    m = panel_metrics(sub)
    m.update({"label": f"n−1 LOO (drop {g})", "category": "LOO"})
    rows.append(m)

# B. Top-N by Cohen's d
ranked = sorted(GENES, key=lambda g: -GENE_D[g])
for k in [3, 4, 5, 6, 7]:
    sub = ranked[:k]
    m = panel_metrics(sub)
    m.update({"label": f"Top-{k} by per-gene d ({' '.join(sub)})", "category": f"top-{k}"})
    rows.append(m)

# C. TSO500 alignment
tso_native = ["TSHR", "PAX8", "NKX2-1"]
m = panel_metrics(tso_native); m.update({"label": "TSO500 NATIVE 3-gene only", "category": "TSO500-native"}); rows.append(m)
# +1, +2, +3 add-ons from highest-d not-on-TSO500
addons_ranked = [g for g in ranked if not TSO500.get(g, False)]
for k in [1, 2, 3, 5]:
    sub = tso_native + addons_ranked[:k]
    m = panel_metrics(sub)
    m.update({"label": f"TSO500 + top {k} add-on ({len(sub)}-gene)", "category": f"TSO500+{k}"})
    rows.append(m)

# D. Hand-picked compact panels
hand = {
    "RAI-uptake minimal (TSHR+NIS+TPO+TG)":        ["TSHR", "SLC5A5", "TPO", "TG"],
    "TF-anchor (PAX8+NKX2-1+FOXE1)":                ["PAX8", "NKX2-1", "FOXE1"],
    "Lineage core (PAX8+NKX2-1+TSHR)":              ["PAX8", "NKX2-1", "TSHR"],
    "Landa overlap 5 (TG+TSHR+TPO+PAX8+DIO1)":      ["TG", "TSHR", "TPO", "PAX8", "DIO1"],
    "Compact 4 (PAX8+NKX2-1+TPO+DIO1)":             ["PAX8", "NKX2-1", "TPO", "DIO1"],
    "Compact 5 (PAX8+NKX2-1+TSHR+TPO+DIO1)":        ["PAX8", "NKX2-1", "TSHR", "TPO", "DIO1"],
    "Compact 6 (drop FOXE1+SLC5A5)":                ["PAX8", "NKX2-1", "TG", "TPO", "TSHR", "DIO1"]
}
for label, sub in hand.items():
    m = panel_metrics(sub)
    m.update({"label": label, "category": "hand-picked"})
    rows.append(m)

df_out = pd.DataFrame(rows)
df_out.to_csv(f"{OUT}/reduction_summary.tsv", sep="\t", index=False)

# Save JSON for web
recs = df_out[["label","category","n_genes","genes","cohens_d","auc_5fcv","r_vs_full",
               "tso500_native","ihc_routine","tf_count","ef_count"]].to_dict(orient="records")
with open(f"{OUT}/reduction_summary.json", "w") as f:
    json.dump(recs, f, indent=2, ensure_ascii=False)

# ─── FIG 1: panel size vs metrics curve ───
fig, ax = plt.subplots(figsize=(11, 5.5), dpi=170)
# Plot all candidates
xs, ys_d, ys_auc, colors, labels = [], [], [], [], []
cat_color = {"baseline":"#0F172A", "LOO":"#94A3B8", "top-3":"#B91C1C", "top-4":"#B91C1C", "top-5":"#B91C1C",
             "top-6":"#B91C1C", "top-7":"#B91C1C", "TSO500-native":"#1E40AF",
             "TSO500+1":"#0E7490", "TSO500+2":"#0E7490", "TSO500+3":"#0E7490", "TSO500+5":"#0E7490",
             "hand-picked":"#B45309"}
for r in rows:
    xs.append(r["n_genes"]); ys_d.append(r["cohens_d"]); ys_auc.append(r["auc_5fcv"]*100)
    colors.append(cat_color.get(r["category"], "#475569"))
ax2 = ax.twinx()
ax.scatter(xs, ys_d, c=colors, s=110, edgecolor="white", linewidth=1.2, zorder=3, alpha=0.92)
for r in rows:
    if r["category"] == "baseline":
        ax.annotate("FULL 8", xy=(r["n_genes"], r["cohens_d"]), xytext=(0, 8), textcoords="offset points",
                    ha="center", fontsize=9, fontweight="bold", color="#0F172A")
    elif r["category"] == "TSO500-native":
        ax.annotate(f"TSO500\nnative (3)", xy=(r["n_genes"], r["cohens_d"]), xytext=(-22, 14), textcoords="offset points",
                    ha="center", fontsize=8.5, fontweight="bold", color="#1E40AF",
                    arrowprops=dict(arrowstyle="-", color="#1E40AF", lw=0.6))
ax.set_xlabel("Panel size (number of genes)", fontsize=11)
ax.set_ylabel("Cohen's d  (DM1 vs DM2)", fontsize=11, color="#B91C1C")
ax.tick_params(axis="y", colors="#B91C1C")
ax.axhline(ref["cohens_d"], color="#0F172A", ls="--", lw=1, alpha=0.5, label=f"Full-8 reference (d = {ref['cohens_d']:.2f})")
ax.spines["top"].set_visible(False)
ax.set_title("Figure GR-1 · 패널 축소 시 Cohen's d 변화 — 안전한 축소 범위 확인", fontsize=12.5, fontweight="bold", loc="left", pad=12)
ax.set_xticks(range(2, 9)); ax.set_xlim(2.4, 8.6)
ax.legend(loc="lower right", frameon=False, fontsize=9)
# Category legend
from matplotlib.lines import Line2D
le = [Line2D([0],[0], marker="o", color="w", markerfacecolor=v, markersize=9, label=k)
      for k, v in [("Full panel","#0F172A"),("LOO n-1","#94A3B8"),("Top-N by d","#B91C1C"),
                   ("TSO500 native","#1E40AF"),("TSO500 + add-on","#0E7490"),("Hand-picked","#B45309")]]
ax.legend(handles=le, loc="upper left", frameon=False, fontsize=8.5, ncol=2)
plt.tight_layout()
fig.savefig(f"{OUT}/fig_reduction_tradeoff.png", bbox_inches="tight")
plt.close(fig)

# ─── FIG 2: recommended panels visual comparison ───
rec_labels = [
    "FULL 8-gene panel (reference)",
    "Compact 6 (drop FOXE1+SLC5A5)",
    "Compact 5 (PAX8+NKX2-1+TSHR+TPO+DIO1)",
    "Top-5 by per-gene d (TPO DIO1 TSHR PAX8 TG)",
    "Compact 4 (PAX8+NKX2-1+TPO+DIO1)",
    "TSO500 + top 2 add-on (5-gene)",
    "Lineage core (PAX8+NKX2-1+TSHR)",
    "TSO500 NATIVE 3-gene only"
]
rec_rows = [r for r in rows if r["label"] in rec_labels]
# preserve order
rec_rows = sorted(rec_rows, key=lambda r: -r["cohens_d"])
fig, ax = plt.subplots(figsize=(13.5, 6), dpi=170)
ys = np.arange(len(rec_rows))[::-1]
ds = [r["cohens_d"] for r in rec_rows]
aucs = [r["auc_5fcv"] for r in rec_rows]
rs = [r["r_vs_full"] for r in rec_rows]
labels = [r["label"] for r in rec_rows]
# bars for Cohen's d
ax.barh(ys, ds, color=["#B91C1C" if "FULL" not in l else "#0F172A" for l in labels], edgecolor="white")
for i, (y, d, auc, r) in enumerate(zip(ys, ds, aucs, rs)):
    ax.text(d + 0.04, y, f"d = {d:.2f}  ·  CV-AUC = {auc:.3f}  ·  r vs full = {r:.3f}",
            va="center", fontsize=9.5, color="#0F172A", fontweight="bold")
ax.set_yticks(ys); ax.set_yticklabels(labels, fontsize=10)
ax.set_xlabel("Cohen's d (DM1 vs DM2)  on TCGA HM450 (n = 503)", fontsize=10.5)
ax.set_xlim(0, max(ds)*1.45)
ax.set_title("Figure GR-2 · 권장 reduced panels 비교 — 어느 정도 축소가 안전한가?",
             fontsize=12.5, fontweight="bold", loc="left", pad=12)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.text(0, 1.05, "FULL 8-gene 기준 대비 d 손실 ≤ 15 % 가 \"안전 축소\" · ≤ 30 % 가 \"수용 가능\" · 그 이상은 정보 손실 큼",
        transform=ax.transAxes, fontsize=10, style="italic", color="#475569")
plt.tight_layout()
fig.savefig(f"{OUT}/fig_reduction_recommend.png", bbox_inches="tight")
plt.close(fig)

# Print compact summary
print("\n=== Top recommended panels (sorted by Cohen's d) ===")
print(f"{'label':<60s}  {'n':>3s}  {'d':>6s}  {'AUC':>6s}  {'r(full)':>8s}  {'TSO500':>7s}")
for r in sorted(rec_rows, key=lambda x: -x["cohens_d"]):
    print(f"  {r['label']:<58s}  {r['n_genes']:3d}  {r['cohens_d']:6.2f}  {r['auc_5fcv']:6.3f}  {r['r_vs_full']:8.3f}  {r['tso500_native']:>7s}")

print(f"\nALL outputs to: {OUT}")
for f in sorted(os.listdir(OUT)):
    if "reduction" in f: print(f"  {f}")
