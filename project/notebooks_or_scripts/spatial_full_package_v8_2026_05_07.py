#!/usr/bin/env python3
"""Spatial v8 — 8 more figures (S_F57 ~ S_F64) + comprehensive tables."""
from __future__ import annotations
from pathlib import Path
import numpy as np, pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import spearmanr, mannwhitneyu, kendalltau
import warnings; warnings.filterwarnings("ignore")

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT/"project/results/spatial_full_2026_05_06"
ASSETS = ROOT/"project/papers_hub_2026_05_04/assets/spatial_full"
G521H = ROOT/"project/data/processed/GSE250521"
EXTH = ROOT/"project_external_st/data/processed"

PAL = {"PT":"#3C6B4F","PTC":"#34547A","LPTC":"#B8893C","ATC":"#7B1F2A",
       "PTC_HT":"#962E2E","CONTROL":"#52525a","HT":"#34547A","GD":"#B8893C"}
ORDER = ["PT","PTC","LPTC","ATC"]
trop_cmap = LinearSegmentedColormap.from_list("trop2",["#fffaf2","#ffd58a","#d97c2e","#7B1F2A"],N=256)
dm_cmap   = LinearSegmentedColormap.from_list("dm",["#34547A","#fffaf2","#7B1F2A"],N=256)

def score_h5ad(p, gs_list):
    a = ad.read_h5ad(p)
    rvar = a.raw.var_names.astype(str) if a.raw is not None else a.var_names.astype(str)
    X = a.raw.X if a.raw is not None else a.X
    idx = [list(rvar.values).index(g) for g in gs_list if g in rvar.values]
    if not idx: return None, None
    sub = X[:, idx]
    if hasattr(sub,"toarray"): sub = sub.toarray()
    tot = np.array(a.obs.get("total_counts", np.ones(a.n_obs))).flatten()
    tot = np.where(tot==0,1,tot)
    return np.log1p(sub/tot[:,None]*1e4), (a.obs["array_row"].values, a.obs["array_col"].values) if "array_row" in a.obs.columns else (None, None)

print("[load]")
g521 = pd.read_csv(ROOT/"project/results/01_spatial_score/all_spots_scored.tsv.gz", sep="\t").rename(columns={"stage":"condition"})
trop_all = pd.read_csv(OUT/"spatial_per_spot_TROP2.tsv.gz", sep="\t")
df_E = pd.read_csv(OUT/"spatial_E_TROP2_spot_distribution.tsv", sep="\t")
df_C = pd.read_csv(OUT/"spatial_C_HT_TLS_spatial.tsv", sep="\t")
df_moran = pd.read_csv(OUT/"spatial_F40_per_sample_per_axis_morans.tsv", sep="\t")

# ===== S_F57: Per-sample TROP2 violin per stage =====
print("[F57] Per-sample TROP2 violin")
fig, ax = plt.subplots(figsize=(15, 5.5))
positions = []; data_v = []; colors_v = []
g521_samples_by_stage = {}
for s in ORDER:
    g521_samples_by_stage[s] = sorted([h for h in trop_all[(trop_all.dataset=="GSE250521")&(trop_all.condition==s)]["sample_id"].unique()])
x = 0
xticks = []; xlabels = []
for stage in ORDER:
    for sample in g521_samples_by_stage[stage]:
        sub = trop_all[trop_all.sample_id == sample]
        nonzero = sub[sub.TROP2 > 0.001]["TROP2"].values
        if len(nonzero) < 5: continue
        positions.append(x); data_v.append(nonzero); colors_v.append(PAL[stage])
        xticks.append(x); xlabels.append(f"{sample.split('_')[-1]}\n[{stage}]"); x += 1
    x += 0.7
parts = ax.violinplot(data_v, positions=positions, widths=0.8, showmeans=True, showmedians=True)
for body, c in zip(parts['bodies'], colors_v):
    body.set_facecolor(c); body.set_alpha(0.7); body.set_edgecolor("black")
parts['cmeans'].set_color("black"); parts['cmedians'].set_color("white")
ax.set_xticks(xticks); ax.set_xticklabels(xlabels, fontsize=8, rotation=70, ha="right")
ax.set_ylabel("TROP2 (log-norm, non-zero spots only)")
ax.set_title("S_F57  Per-sample TROP2 expression distribution — 16 GSE250521 slides (non-zero spots only)\n"
             "PTC+LPTC samples: distribution shifted to higher mean + wider variance (niche zones)",
             fontsize=12, fontweight="bold")
ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(ASSETS/"S_F57_TROP2_per_sample_violin.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F58: GSE230424 PTC+HT TROP2 + DM1 RGB overlay =====
print("[F58] GSE230424 RGB overlay")
g230_paths = sorted((EXTH/"GSE230424").glob("*/GSM*.raw.h5ad"))
ext_dfs = [pd.read_csv(f, sep="\t") for f in (ROOT/"project_external_st/results/scores").glob("GSE230424_*.tsv.gz")]
ext_g230 = pd.concat(ext_dfs, ignore_index=True).rename(columns={"DM1_like_score_raw":"DM1_like_score","Epithelial_score_raw":"Epithelial_score"})
fig = plt.figure(figsize=(15, 4.4))
g = gs.GridSpec(1, 4, wspace=0.18, top=0.85, bottom=0.04, left=0.04, right=0.97)
def to01(x):
    x = np.array(x, dtype=float)
    x = x - np.nanmin(x); x = x/np.nanmax(x) if np.nanmax(x)>0 else x
    return np.nan_to_num(x, 0)
for c, p in enumerate(g230_paths):
    sample = p.parent.name
    ax = fig.add_subplot(g[0, c])
    sub_t = trop_all[trop_all.sample_id == sample]
    sub_g = ext_g230[ext_g230.sample_id == sample]
    if sub_t.empty or sub_g.empty: ax.axis("off"); continue
    sub_t["row_col"] = sub_t["sample_id"]+"_"+sub_t["array_row"].astype(str)+"_"+sub_t["array_col"].astype(str)
    sub_g["row_col"] = sub_g["sample_id"]+"_"+sub_g["array_row"].astype(str)+"_"+sub_g["array_col"].astype(str)
    mg = sub_t.merge(sub_g[["row_col","DM1_like_score","Epithelial_score"]], on="row_col")
    if mg.empty: ax.axis("off"); continue
    R = to01(mg["TROP2"].values)
    G = to01(mg["Epithelial_score"].values)
    B = to01(np.maximum(mg["DM1_like_score"].values, 0))
    rgb = np.stack([R, G, B], axis=1)
    ax.scatter(mg["array_col"].values, -mg["array_row"].values, c=rgb, s=4, alpha=0.92)
    ax.set_title(f"PTC+HT · {sample.split('_')[-1]}\nR=TROP2 G=Epi B=DM1+", fontsize=10, color="#962E2E", fontweight="bold")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
fig.suptitle("S_F58  GSE230424 PTC+HT — 3-channel RGB overlay (S_F27 parallel for cancer cohort)\n"
             "TROP2 weak in HT-overlap PTC; Epi+DM1 dominant; HT-specific TLS niche은 immune marker로 (S_F8)",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F58_GSE230424_RGB_overlay.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F59: Per-stage TROP2 enrichment (CDF + KDE) =====
print("[F59] Per-stage TROP2 CDF")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
ax = axes[0]
for stage, color in zip(ORDER, ["#3C6B4F","#34547A","#B8893C","#7B1F2A"]):
    sub = trop_all[(trop_all.dataset=="GSE250521") & (trop_all.condition==stage)]
    vals = sub["TROP2"].values
    sorted_vals = np.sort(vals)
    cdf = np.arange(1, len(sorted_vals)+1) / len(sorted_vals)
    ax.plot(sorted_vals, cdf, color=color, lw=2, label=f"{stage} (n_spots={len(vals)})")
ax.set_xlabel("TROP2 (log-norm)"); ax.set_ylabel("Cumulative probability")
ax.set_title("A. TROP2 expression CDF per stage\n(curve가 우측으로 shift = TROP2-high spot 비율 증가)", fontsize=11, fontweight="bold")
ax.legend(); ax.grid(alpha=0.3)
ax = axes[1]
# KDE-like density for non-zero only
for stage, color in zip(ORDER, ["#3C6B4F","#34547A","#B8893C","#7B1F2A"]):
    sub = trop_all[(trop_all.dataset=="GSE250521") & (trop_all.condition==stage) & (trop_all.TROP2 > 0.001)]
    vals = sub["TROP2"].values
    if len(vals) > 50:
        ax.hist(vals, bins=40, density=True, alpha=0.4, color=color, label=f"{stage}", edgecolor="black", linewidth=0.3)
ax.set_xlabel("TROP2 (log-norm, non-zero)"); ax.set_ylabel("Density")
ax.set_title("B. TROP2 expression density per stage (non-zero spots)\n(분포 shape difference visible)", fontsize=11, fontweight="bold")
ax.legend(); ax.grid(axis="y", alpha=0.3)
fig.suptitle("S_F59  Per-stage TROP2 cumulative + density — distribution shape differences",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.93])
fig.savefig(ASSETS/"S_F59_TROP2_CDF_density.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F60: TROP2 niche cluster ID color visualization (single PTC sample) =====
print("[F60] Niche cluster ID visualization")
def visium_hex_clusters_with_ids(rows, cols, mask):
    pos_to_idx = {}
    for i, (r, c) in enumerate(zip(rows, cols)):
        if mask[i]: pos_to_idx[(int(r), int(c))] = i
    parent = {p: p for p in pos_to_idx.values()}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    for (r, c), i in pos_to_idx.items():
        for dr, dc in [(-1,-1),(-1,1),(0,-2),(0,2),(1,-1),(1,1)]:
            j = pos_to_idx.get((r+dr, c+dc))
            if j is not None: union(i, j)
    cluster_ids = np.full(len(rows), -1, dtype=int)
    root_to_id = {}
    next_id = 0
    for (r, c), i in pos_to_idx.items():
        root = find(i)
        if root not in root_to_id:
            root_to_id[root] = next_id; next_id += 1
        cluster_ids[i] = root_to_id[root]
    return cluster_ids

# Pick 4 PTC+LPTC samples with highest Moran I
top_samples = df_E.sort_values("morans_I_TROP2", ascending=False).head(4)["sample_id"].tolist()
fig = plt.figure(figsize=(15, 4.4))
g = gs.GridSpec(1, 4, wspace=0.15, top=0.85, bottom=0.04, left=0.04, right=0.97)
import matplotlib.colors as mcolors
for c, sample in enumerate(top_samples):
    ax = fig.add_subplot(g[0, c])
    sub = trop_all[trop_all.sample_id == sample]
    if "array_row" not in sub.columns: continue
    rows_arr = sub["array_row"].astype(int).values
    cols_arr = sub["array_col"].astype(int).values
    threshold = sub["TROP2"].quantile(0.85)
    mask = (sub["TROP2"] > threshold).values
    cluster_ids = visium_hex_clusters_with_ids(rows_arr, cols_arr, mask)
    # background spots gray
    bg = ~mask
    ax.scatter(cols_arr[bg], -rows_arr[bg], c="#dddddd", s=2, alpha=0.5)
    # cluster spots colored by cluster ID (use tab20)
    if cluster_ids[mask].size:
        n_clusters = cluster_ids[mask].max() + 1
        cmap = plt.get_cmap("tab20", max(20, n_clusters))
        ax.scatter(cols_arr[mask], -rows_arr[mask], c=cluster_ids[mask], cmap=cmap, s=8, alpha=0.95)
        cond = sub["condition"].iloc[0]
        ax.set_title(f"{cond} · {sample.split('_')[-1]}\n{n_clusters} clusters identified", fontsize=10, color=PAL.get(cond,"#666"), fontweight="bold")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
fig.suptitle("S_F60  TROP2-high niche cluster IDs visualized (top 4 niche-organized samples)\n"
             "Each color = 1 contiguous niche cluster (Visium hex 6-neighbor adjacency); 회색 = background",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F60_niche_cluster_IDs.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F61: Closure battery target × model heatmap =====
print("[F61] Closure battery target × model")
clos_metric = pd.read_csv(ROOT/"project/results/03_pathology_poc/closure_battery_metrics.tsv", sep="\t")
sec_a = clos_metric[clos_metric["section"].astype(str).str.startswith("A")].copy()
sec_a["target"] = sec_a["experiment"].astype(str).str.replace("raw_","")
fig, ax = plt.subplots(figsize=(11, 4))
targets = sorted(sec_a["target"].unique())
models = sorted(sec_a["model"].unique())
mat = np.full((len(targets), len(models)), np.nan)
for i, t in enumerate(targets):
    for j, m in enumerate(models):
        sub = sec_a[(sec_a.target == t) & (sec_a.model == m)]
        if not sub.empty: mat[i, j] = sub["pooled_spearman_r"].iloc[0]
im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-0.3, vmax=0.3)
for i in range(len(targets)):
    for j in range(len(models)):
        if not np.isnan(mat[i,j]):
            color = "white" if abs(mat[i,j]) > 0.18 else "black"
            ax.text(j, i, f"{mat[i,j]:.3f}", ha="center", va="center", color=color, fontsize=11, fontweight="bold")
ax.set_xticks(range(len(models))); ax.set_xticklabels(models, fontsize=11)
ax.set_yticks(range(len(targets))); ax.set_yticklabels(targets, fontsize=10.5, family="monospace")
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Pooled Spearman ρ")
ax.set_title("S_F61  Closure battery section A — target × model pooled ρ heatmap\n"
             "All cells fall within ±0.1 — H&E ResNet50 cannot predict any DM1/RAI/TDS axis at useful threshold (0.3)",
             fontsize=11, fontweight="bold")
fig.tight_layout()
fig.savefig(ASSETS/"S_F61_closure_target_model_heatmap.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F62: GSE230424 spot-level DM1 vs TROP2 vs immune-set (joint scatter) =====
print("[F62] GSE230424 joint axis scatter")
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
sample = "GSM7221917_P3"  # P3 has strongest TLS niche
sub_g = ext_g230[ext_g230.sample_id == sample]
sub_g["row_col"] = sub_g["sample_id"]+"_"+sub_g["array_row"].astype(str)+"_"+sub_g["array_col"].astype(str)
sub_t = trop_all[trop_all.sample_id == sample]
sub_t["row_col"] = sub_t["sample_id"]+"_"+sub_t["array_row"].astype(str)+"_"+sub_t["array_col"].astype(str)
mg_p3 = sub_g.merge(sub_t[["row_col","TROP2"]], on="row_col")
# also get HLA-II per spot
hla_path = next((p for p in g230_paths if p.parent.name == sample), None)
if hla_path:
    HLA_II = ["HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1"]
    scores, _ = score_h5ad(hla_path, HLA_II)
    if scores is not None:
        a = ad.read_h5ad(hla_path)
        df_hla = pd.DataFrame({"spot_id": a.obs_names.astype(str), "HLA_II": scores.mean(axis=1)})
        # match by spot_id directly
        if "spot_id" in mg_p3.columns:
            mg_p3 = mg_p3.merge(df_hla, on="spot_id", how="left")
# 4 panels: TROP2 vs DM1, TROP2 vs HLA-II, DM1 vs HLA-II, DM1 vs Epi
plots = [
    ("DM1_like_score","TROP2","DM1_like × TROP2"),
    ("HLA_II","TROP2","HLA-II × TROP2"),
    ("HLA_II","DM1_like_score","HLA-II × DM1_like"),
    ("Epithelial_score","DM1_like_score","Epithelial × DM1_like"),
]
for ax, (xcol, ycol, label) in zip(axes.flat, plots):
    if xcol not in mg_p3.columns or ycol not in mg_p3.columns: ax.axis("off"); continue
    sub_d = mg_p3.dropna(subset=[xcol, ycol])
    if sub_d.empty: ax.axis("off"); continue
    ax.hexbin(sub_d[xcol], sub_d[ycol], gridsize=40, cmap="OrRd", mincnt=2)
    rho, p = spearmanr(sub_d[xcol], sub_d[ycol], nan_policy="omit")
    ax.set_xlabel(xcol); ax.set_ylabel(ycol)
    ax.set_title(f"{label}\nrho = {rho:.3f}, p = {p:.1e}", fontsize=10.5)
    ax.axhline(0, ls=":", lw=0.4, color="gray"); ax.axvline(0, ls=":", lw=0.4, color="gray")
fig.suptitle(f"S_F62  GSE230424 {sample.split('_')[-1]} (P3, strongest TLS niche) — joint axis spot scatter\n"
             "DM1 × HLA-II 강한 양의 상관 = HT-overlap PTC의 immune-active dark matter 영역 collocation",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig(ASSETS/"S_F62_GSE230424_joint_axes_P3.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F63: 28-sample comprehensive summary table figure =====
print("[F63] 28-sample comprehensive summary table")
df_summary = df_E.copy()
df_summary = df_summary.sort_values(["dataset","condition","sample_id"]).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(14, 11))
ax.axis("off")
# build table
cols_show = ["sample_id","dataset","condition","n_spots","TROP2_mean","TROP2_p90","TROP2_high_spot_frac","morans_I_TROP2"]
df_short = df_summary[cols_show].copy()
df_short["sample_id"] = df_short["sample_id"].apply(lambda x: x.split("_")[-1])
df_short["dataset"] = df_short["dataset"].apply(lambda x: x.replace("GSE","").replace("GSE",""))
df_short["TROP2_mean"] = df_short["TROP2_mean"].apply(lambda x: f"{x:.2f}")
df_short["TROP2_p90"] = df_short["TROP2_p90"].apply(lambda x: f"{x:.1f}")
df_short["TROP2_high_spot_frac"] = df_short["TROP2_high_spot_frac"].apply(lambda x: f"{x*100:.1f}%")
df_short["morans_I_TROP2"] = df_short["morans_I_TROP2"].apply(lambda x: f"{x:+.3f}")
df_short.columns = ["sample","dataset","condition","n_spots","mean","p90","high%","Moran I"]

# color per row by condition
row_colors = [PAL.get(c, "#fff") for c in df_summary["condition"]]
row_colors_alpha = []
for c in df_summary["condition"]:
    h = PAL.get(c, "#999")
    # add alpha
    row_colors_alpha.append((int(h[1:3],16)/255, int(h[3:5],16)/255, int(h[5:7],16)/255, 0.18))

table = ax.table(cellText=df_short.values, colLabels=df_short.columns, loc="center", cellLoc="center")
table.auto_set_font_size(False); table.set_fontsize(9)
table.scale(1, 1.4)
# color header
for j in range(len(df_short.columns)):
    table[(0, j)].set_facecolor("#17212f")
    table[(0, j)].set_text_props(color="white", weight="bold")
# color rows
for i in range(len(df_short)):
    for j in range(len(df_short.columns)):
        table[(i+1, j)].set_facecolor(row_colors_alpha[i])
        # bold Moran I if niche
        if j == len(df_short.columns)-1:
            moran_v = df_summary.iloc[i]["morans_I_TROP2"]
            if not np.isnan(moran_v) and moran_v > 0.25:
                table[(i+1, j)].set_text_props(weight="bold", color="#7B1F2A")

fig.suptitle("S_F63  Comprehensive 28-sample TROP2 summary table\n"
             "Bold red Moran I = niche-organized (>0.25); rows colored by condition",
             fontsize=12, fontweight="bold", y=0.97)
fig.savefig(ASSETS/"S_F63_28_sample_summary_table.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F64: Final summary heatmap — sample × marker Moran I =====
print("[F64] Final 28-sample × multi-axis Moran I heatmap")
# Combine: for each sample, get all available Moran I per axis
combined_morans = []
for _, r in df_E.iterrows():
    rec = {"sample": r["sample_id"], "dataset": r["dataset"], "condition": r["condition"],
           "TROP2": r["morans_I_TROP2"]}
    # GSE250521: pull from df_moran
    if r["dataset"] == "GSE250521":
        m_row = df_moran[df_moran["sample"] == r["sample_id"]]
        if not m_row.empty:
            for axis in ["DM1_like_score","Epithelial_score","Proliferation_score","CAF_ECM_score","EMT_score","Hypoxia_score","RAI_8_score"]:
                if axis in m_row.columns:
                    rec[axis.replace("_score","")] = m_row.iloc[0][axis]
    # GSE230424: HT TLS markers
    if r["dataset"] == "GSE230424":
        for marker in ["HLA_II","B_cell","TLS","IGHV_AICDA"]:
            cs = df_C[(df_C["sample_id"]==r["sample_id"]) & (df_C["gene_set"]==marker)]
            if not cs.empty: rec[marker] = cs["morans_I"].iloc[0]
    combined_morans.append(rec)
df_combined = pd.DataFrame(combined_morans)
df_combined.to_csv(OUT/"spatial_F64_28sample_multiaxis_morans.tsv", sep="\t", index=False)

fig, ax = plt.subplots(figsize=(14, 12))
df_combined_sorted = df_combined.sort_values(["dataset","condition","sample"]).reset_index(drop=True)
axes_show = ["TROP2","DM1_like","Epithelial","Proliferation","CAF_ECM","EMT","Hypoxia","HLA_II","B_cell","TLS","IGHV_AICDA"]
mat = np.full((len(df_combined_sorted), len(axes_show)), np.nan)
for i in range(len(df_combined_sorted)):
    for j, a in enumerate(axes_show):
        v = df_combined_sorted.iloc[i].get(a, np.nan)
        if not pd.isna(v): mat[i, j] = float(v)
im = ax.imshow(mat, aspect="auto", cmap="RdYlGn_r", vmin=-0.1, vmax=0.85)
for i in range(len(df_combined_sorted)):
    for j in range(len(axes_show)):
        v = mat[i,j]
        if not np.isnan(v):
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    color="white" if v>0.5 else "black", fontsize=8.5, fontweight="bold")
ax.set_xticks(range(len(axes_show))); ax.set_xticklabels(axes_show, fontsize=10, family="monospace", rotation=30, ha="right")
ax.set_yticks(range(len(df_combined_sorted)))
ax.set_yticklabels([f"{r['condition']}·{r['sample'].split('_')[-1]}" for _,r in df_combined_sorted.iterrows()], fontsize=9, family="monospace")
plt.colorbar(im, ax=ax, fraction=0.04, pad=0.04, label="Moran's I")
ax.set_title("S_F64  Master Moran's I matrix — 28 samples × 11 axes\n"
             "(빈칸 = 해당 cohort에서 axis 미측정; 색 = niche organization 강도)",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(ASSETS/"S_F64_master_moran_matrix.png", dpi=160, bbox_inches="tight"); plt.close(fig)

print("\nDONE — 8 new figures (S_F57 ~ S_F64) saved")
