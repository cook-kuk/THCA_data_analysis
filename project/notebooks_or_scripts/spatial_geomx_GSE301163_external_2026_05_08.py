#!/usr/bin/env python3
"""GSE301163 GeoMx DSP cross-platform external validation for Paper 2 spatial.
n=78 ROIs (Patient 1 + 2) · PanCK+ (32) / VIM+ (46) compartment · histology trajectory
Normal → Normal_adjacent → Non-neoplastic nodular → microPTC → PDTC + Tumor Capsule
+ DICER1 / DGCR8 / non-mutated genotype layer.

Output: 6 figures S_F101 ~ S_F106 + 1 summary TSV.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import mannwhitneyu, spearmanr
import warnings; warnings.filterwarnings("ignore")

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
DATA = ROOT/"project/data/external_geomx_GSE301163"
META = ROOT/"project/results/geo_search_2026_05_08/GSE301163_dsp_metadata.tsv"
OUT = ROOT/"project/results/spatial_full_2026_05_06"
ASSETS = ROOT/"project/papers_hub_2026_05_04/assets/spatial_full"

# ---------- Load matrix + metadata ----------
print("[load] GSE301163 matrix + metadata")
mat = pd.read_csv(DATA/"GSE301163_Normalized_ST_matrix.txt.gz", sep="\t", index_col=0)
print(f"  matrix: {mat.shape}  (genes x ROIs)")
meta = pd.read_csv(META, sep="\t")
# matrix col names: Sample_DSP-...; meta dsp col: Sample_DSP-...
sample_to_dsp = dict(zip(meta["dsp"], meta.to_dict("records")))
mat_cols = list(mat.columns)
matched = [c for c in mat_cols if c in sample_to_dsp]
print(f"  matrix ROIs matched to metadata: {len(matched)} / {len(mat_cols)}")
mat = mat[matched]

# Build sample-level metadata aligned to matrix columns
meta_aligned = pd.DataFrame([sample_to_dsp[c] for c in matched])
meta_aligned["dsp"] = matched
print(f"  meta_aligned shape: {meta_aligned.shape}")
meta_aligned.to_csv(OUT/"GSE301163_meta_aligned.tsv", sep="\t", index=False)

# ---------- Gene sets (matched to our spatial-full pipeline) ----------
GENE_SETS = {
    "TROP2":            ["TACSTD2"],
    "DM1_axis":         ["TG","TPO","TSHR","SLC5A5","DUOX1","DUOX2","IYD","DIO1","DIO2"],
    "RAI_8":            ["TG","TPO","SLC5A5","FOXE1","TSHR","NKX2-1","PAX8","DIO1"],
    "HLA_II":           ["HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1"],
    "B_cell":           ["CD19","MS4A1","CD79A","CD79B","SDC1","JCHAIN"],
    "TLS":              ["CXCL13","CCL19","CCL21","CXCR5","CCR7","SELL"],
    "IGHV_AICDA":       ["AICDA","IGHM","IGHG1","IGHA1","IGKC","IGLC2"],
    "Thyroid_TF":       ["PAX8","NKX2-1","FOXE1","HHEX","TITF2"],
    "Cell_cycle":       ["MKI67","CCNB1","CDK1","MCM2","TOP2A","BIRC5","UBE2C","CCNE1"],
    "Epithelial":       ["KRT8","KRT18","KRT19","EPCAM","CDH1"],
    "Stromal":          ["VIM","FAP","ACTA2","COL1A1","COL3A1","PDGFRB"],
    "EMT":              ["CDH1","CDH2","VIM","ZEB1","ZEB2","SNAI1","SNAI2","TWIST1","FN1"],
    "Hypoxia":          ["HIF1A","EPAS1","VEGFA","CA9","BNIP3","SLC2A1","NDRG1","LDHA"],
    "Eight_gene_DM1":   ["TIMP1","FN1","COL1A1","COL3A1","SPARC","POSTN","SERPINH1","COL5A1"],
}

# ---------- Per-ROI score ----------
print("[score] per-ROI gene-set means (log1p of normalized count)")
scores = {}
for set_name, genes in GENE_SETS.items():
    avail = [g for g in genes if g in mat.index]
    if not avail:
        scores[set_name] = pd.Series(np.nan, index=mat.columns); continue
    sub = mat.loc[avail]
    log_norm = np.log1p(sub).mean(axis=0)
    scores[set_name] = log_norm
    print(f"  {set_name:18s} {len(avail)}/{len(genes)} genes available")

scores_df = pd.DataFrame(scores).reset_index().rename(columns={"index":"dsp"})
df = meta_aligned.merge(scores_df, on="dsp")
df.to_csv(OUT/"GSE301163_per_ROI_scores.tsv", sep="\t", index=False)
print(f"[write] per-ROI scores: {OUT/'GSE301163_per_ROI_scores.tsv'}  shape={df.shape}")

# ---------- Comparison stats ----------
HIST_ORDER = ["Normal thyroid","Normal adjacent","Non-neoplastic nodular area",
              "microPTC","PDTC","Tumor Capsule"]
PAL = {"Normal thyroid":"#3C6B4F","Normal adjacent":"#74A88A",
       "Non-neoplastic nodular area":"#888","microPTC":"#34547A",
       "PDTC":"#7B1F2A","Tumor Capsule":"#962E2E"}
GENO_PAL = {"non-DICER1-mutated":"#888","DICER1-mutated":"#7B1F2A","DGCR8-E518K-mutated":"#34547A"}

def mw_p(a, b):
    a = np.asarray(a); b = np.asarray(b)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2: return np.nan, np.nan
    try:
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        return u, p
    except Exception:
        return np.nan, np.nan

# ---------- S_F101 — ROI heatmap (78 × axes) ----------
print("[F101] ROI heatmap")
HEATMAP_AXES = ["TROP2","DM1_axis","RAI_8","Thyroid_TF","Epithelial","Stromal",
                "Cell_cycle","EMT","Hypoxia","HLA_II","B_cell","TLS","IGHV_AICDA","Eight_gene_DM1"]
df_sorted = df.copy()
df_sorted["hist_rank"] = df_sorted["histology"].map({h:i for i,h in enumerate(HIST_ORDER)}).fillna(99)
df_sorted = df_sorted.sort_values(["hist_rank","label","patient","tissue"])
M = df_sorted[HEATMAP_AXES].values
# z-score per axis
M_z = (M - np.nanmean(M,0)) / (np.nanstd(M,0) + 1e-9)
fig = plt.figure(figsize=(15, 14))
gss = gs.GridSpec(1, 2, width_ratios=[14, 1], wspace=0.02)
ax = fig.add_subplot(gss[0])
cmap_div = LinearSegmentedColormap.from_list("div", ["#34547A","#fffaf2","#7B1F2A"], N=256)
im = ax.imshow(M_z, cmap=cmap_div, vmin=-2.4, vmax=2.4, aspect="auto")
ax.set_xticks(range(len(HEATMAP_AXES)))
ax.set_xticklabels(HEATMAP_AXES, rotation=45, ha="right", fontsize=9.5)
ax.set_yticks(range(len(df_sorted)))
labels_left = [f"{r['patient'][-1]} · {r['histology'][:25]:25s} · {r['label']}" for _, r in df_sorted.iterrows()]
ax.set_yticklabels(labels_left, fontsize=6.0)
for tick, (_, r) in zip(ax.get_yticklabels(), df_sorted.iterrows()):
    tick.set_color(PAL.get(r["histology"], "#0F1A2E"))
plt.colorbar(im, cax=fig.add_subplot(gss[1]), label="z-score (per axis)")
ax.set_title("S_F101  GSE301163 GeoMx DSP — 78 ROIs × 14 axes z-scored heatmap (cross-platform external validation)\n"
             "n=78: Patient1+2 · Normal/NormAdj/Nodular/microPTC/PDTC/Capsule × PanCK+ / VIM+ · DICER1 / DGCR8 / non-mutated",
             fontsize=11.5, fontweight="bold")
fig.savefig(ASSETS/"S_F101_GSE301163_ROI_heatmap.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------- S_F102 — TROP2 by compartment × histology ----------
print("[F102] TROP2 PanCK+ vs VIM+ by histology")
fig = plt.figure(figsize=(13, 5.6))
ax = fig.add_subplot(111)
positions = []; data = []; colors = []; lbls = []
xpos = 0
for h in HIST_ORDER:
    for cmp in ["PanCK+","VIM+"]:
        sub = df[(df.histology==h) & (df.label==cmp)]
        if len(sub) < 1: continue
        v = sub["TROP2"].dropna().values
        positions.append(xpos); data.append(v); lbls.append(f"{h[:18]}·{cmp}")
        colors.append(PAL[h] if cmp=="PanCK+" else "#cccccc")
        xpos += 1
    xpos += 0.6  # gap between histology blocks
bp = ax.boxplot(data, positions=positions, widths=0.55, patch_artist=True, showfliers=False,
                medianprops=dict(color="#0F1A2E", lw=1.5))
for b, c in zip(bp["boxes"], colors):
    b.set_facecolor(c); b.set_edgecolor("#0F1A2E"); b.set_alpha(0.85)
for i, (vals, c) in enumerate(zip(data, colors)):
    j = (np.random.RandomState(i).rand(len(vals))-0.5)*0.3
    ax.scatter(np.full(len(vals), positions[i]) + j, vals, color=c,
               edgecolor="black", lw=0.4, s=42, zorder=3)
ax.set_xticks(positions); ax.set_xticklabels(lbls, rotation=70, ha="right", fontsize=8.5)
ax.set_ylabel("TROP2 (TACSTD2 log1p normalized)", fontsize=11)
ax.set_title("S_F102  GSE301163 — TROP2 PanCK+ epithelial vs VIM+ stromal by histology  ·  cross-platform validation\n"
             "Hypothesis (signature-level only): TROP2 enriched in PanCK+ epithelial of tumor histologies (microPTC / PDTC)",
             fontsize=11.5, fontweight="bold")
ax.grid(True, ls=":", alpha=0.4, axis="y")

# Add MW p for tumor compartments
mptc_p = df[(df.histology=="microPTC") & (df.label=="PanCK+")]["TROP2"].values
mptc_v = df[(df.histology=="microPTC") & (df.label=="VIM+")]["TROP2"].values
pdtc_p = df[(df.histology=="PDTC") & (df.label=="PanCK+")]["TROP2"].values
pdtc_v = df[(df.histology=="PDTC") & (df.label=="VIM+")]["TROP2"].values
norm_p = df[(df.histology=="Normal thyroid") & (df.label=="PanCK+")]["TROP2"].values
_, p_mptc_panck_vs_vim = mw_p(mptc_p, mptc_v)
_, p_pdtc_panck_vs_vim = mw_p(pdtc_p, pdtc_v)
_, p_mptc_panck_vs_norm = mw_p(mptc_p, norm_p)
ax.text(0.005, 0.98, (f"microPTC PanCK+ vs VIM+   MW p = {p_mptc_panck_vs_vim:.3g}\n"
                       f"PDTC PanCK+ vs VIM+         MW p = {p_pdtc_panck_vs_vim:.3g}\n"
                       f"microPTC PanCK+ vs Normal PanCK+  MW p = {p_mptc_panck_vs_norm:.3g}"),
        transform=ax.transAxes, va="top", ha="left", fontsize=9, family="monospace",
        bbox=dict(boxstyle="round,pad=0.4", fc="#fff5e0", ec="#B8893C", lw=1))
fig.tight_layout()
fig.savefig(ASSETS/"S_F102_GSE301163_TROP2_compartment.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------- S_F103 — DM1 axis by histology ----------
print("[F103] DM1 axis by histology")
fig = plt.figure(figsize=(11, 5.2))
ax = fig.add_subplot(111)
data = []; positions = []; clrs = []; lbls = []
for i,h in enumerate(HIST_ORDER):
    sub = df[df.histology==h]
    if len(sub) < 2: continue
    data.append(sub["DM1_axis"].dropna().values)
    positions.append(i); clrs.append(PAL[h]); lbls.append(f"{h[:22]}\n(n={len(sub)})")
bp = ax.boxplot(data, positions=positions, widths=0.55, patch_artist=True, showfliers=False)
for b, c in zip(bp["boxes"], clrs):
    b.set_facecolor(c); b.set_edgecolor("#0F1A2E"); b.set_alpha(0.85)
for i,(vals,c) in enumerate(zip(data, clrs)):
    j = (np.random.RandomState(i).rand(len(vals))-0.5)*0.3
    ax.scatter(np.full(len(vals), positions[i]) + j, vals, color=c, edgecolor="black", lw=0.4, s=46, zorder=3)
ax.set_xticks(positions); ax.set_xticklabels(lbls, fontsize=9.5)
for tick, c in zip(ax.get_xticklabels(), clrs):
    tick.set_color(c); tick.set_fontweight("bold")
ax.set_ylabel("DM1 axis (TG·TPO·TSHR·SLC5A5·DUOX1·DUOX2·IYD·DIO1/2 log1p mean)", fontsize=10)
ax.set_title("S_F103  GSE301163 — DM1 axis (thyroid differentiation set) by histology · cross-platform consistency check\n"
             "Visium 28-sample finding (post-fix): PT/PTC/LPTC ~0.4-0.6, ATC collapse to ~0.12 (Moran I). GeoMx = ROI-mean test.",
             fontsize=11, fontweight="bold")
ax.grid(True, ls=":", alpha=0.4, axis="y")
fig.tight_layout()
fig.savefig(ASSETS/"S_F103_GSE301163_DM1_by_histology.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------- S_F104 — Eight-gene DM1 mini-index by histology + genotype ----------
print("[F104] Eight-gene DM1 by histology / genotype")
fig = plt.figure(figsize=(13, 5.4))
g104 = gs.GridSpec(1, 2, wspace=0.28)
ax = fig.add_subplot(g104[0,0])
data=[]; positions=[]; clrs=[]; lbls=[]
for i,h in enumerate(HIST_ORDER):
    sub = df[df.histology==h]
    if len(sub)<2: continue
    data.append(sub["Eight_gene_DM1"].dropna().values); positions.append(i)
    clrs.append(PAL[h]); lbls.append(f"{h[:18]}\n(n={len(sub)})")
bp = ax.boxplot(data, positions=positions, widths=0.55, patch_artist=True, showfliers=False)
for b,c in zip(bp["boxes"],clrs): b.set_facecolor(c); b.set_edgecolor("#0F1A2E"); b.set_alpha(0.85)
for i,(vals,c) in enumerate(zip(data, clrs)):
    j = (np.random.RandomState(i).rand(len(vals))-0.5)*0.3
    ax.scatter(np.full(len(vals), positions[i]) + j, vals, color=c, edgecolor="black", lw=0.4, s=46)
ax.set_xticks(positions); ax.set_xticklabels(lbls, fontsize=9)
ax.set_ylabel("8-gene DM1 (TIMP1·FN1·COL1A1·COL3A1·SPARC·POSTN·SERPINH1·COL5A1) log1p mean", fontsize=9)
ax.set_title("A. by histology", fontsize=10.5, fontweight="bold")
ax.grid(True, ls=":", alpha=0.4, axis="y")

ax2 = fig.add_subplot(g104[0,1])
data=[]; positions=[]; clrs=[]; lbls=[]
for i,(g,c) in enumerate(GENO_PAL.items()):
    sub = df[df.genotype==g]
    if len(sub)<2: continue
    data.append(sub["Eight_gene_DM1"].dropna().values); positions.append(i)
    clrs.append(c); lbls.append(f"{g[:22]}\n(n={len(sub)})")
bp = ax2.boxplot(data, positions=positions, widths=0.55, patch_artist=True, showfliers=False)
for b,c in zip(bp["boxes"],clrs): b.set_facecolor(c); b.set_edgecolor("#0F1A2E"); b.set_alpha(0.85)
for i,(vals,c) in enumerate(zip(data, clrs)):
    j = (np.random.RandomState(i).rand(len(vals))-0.5)*0.3
    ax2.scatter(np.full(len(vals), positions[i]) + j, vals, color=c, edgecolor="black", lw=0.4, s=46)
ax2.set_xticks(positions); ax2.set_xticklabels(lbls, fontsize=9)
ax2.set_ylabel("8-gene DM1 mini-index", fontsize=9)
ax2.set_title("B. by genotype", fontsize=10.5, fontweight="bold")
ax2.grid(True, ls=":", alpha=0.4, axis="y")

fig.suptitle("S_F104  GSE301163 — 8-gene DM1 mini-index (Paper 1 main signature) cross-platform check",
             fontsize=11.5, fontweight="bold", y=1.02)
fig.savefig(ASSETS/"S_F104_GSE301163_8gene_DM1.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------- S_F105 — Thyroid TF / RAI_8 / Differentiation by histology ----------
print("[F105] Differentiation collapse cross-platform check")
fig = plt.figure(figsize=(15, 5.4))
g105 = gs.GridSpec(1, 3, wspace=0.30)
DIFF_PANELS = [("Thyroid_TF","PAX8·NKX2-1·FOXE1·HHEX"),
               ("RAI_8","TG·TPO·SLC5A5·FOXE1·TSHR·NKX2-1·PAX8·DIO1"),
               ("DM1_axis","TG·TPO·TSHR·SLC5A5·DUOX1/2·IYD·DIO1/2")]
for k,(axis,desc) in enumerate(DIFF_PANELS):
    ax = fig.add_subplot(g105[0,k])
    data=[]; positions=[]; clrs=[]; lbls=[]
    for i,h in enumerate(HIST_ORDER):
        sub = df[df.histology==h]
        if len(sub)<2: continue
        data.append(sub[axis].dropna().values); positions.append(i)
        clrs.append(PAL[h]); lbls.append(f"{h[:14]}")
    bp = ax.boxplot(data, positions=positions, widths=0.55, patch_artist=True, showfliers=False)
    for b,c in zip(bp["boxes"],clrs): b.set_facecolor(c); b.set_edgecolor("#0F1A2E"); b.set_alpha(0.85)
    for i,(vals,c) in enumerate(zip(data, clrs)):
        j = (np.random.RandomState(i+k).rand(len(vals))-0.5)*0.3
        ax.scatter(np.full(len(vals), positions[i]) + j, vals, color=c, edgecolor="black", lw=0.4, s=42)
    ax.set_xticks(positions); ax.set_xticklabels(lbls, rotation=45, ha="right", fontsize=8.5)
    ax.set_title(f"{axis}\n{desc[:55]}", fontsize=10, fontweight="bold")
    ax.set_ylabel("log1p mean (normalized count)", fontsize=9)
    ax.grid(True, ls=":", alpha=0.4, axis="y")
fig.suptitle("S_F105  GSE301163 — Differentiation collapse cross-platform check (3 axes × 6 histology) · PDTC = collapse",
             fontsize=11.5, fontweight="bold", y=1.02)
fig.savefig(ASSETS/"S_F105_GSE301163_differentiation.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------- S_F106 — Cross-platform consistency summary ----------
print("[F106] Cross-platform consistency summary")
# Visium GSE250521 Moran I means (post-fix)
visium_mor = pd.read_csv(OUT/"spatial_morans_v7plus_2026_05_07.tsv", sep="\t")
visium_g521 = visium_mor[visium_mor["dataset"]=="GSE250521"]
visium_means = visium_g521.groupby(["condition","axis"])["morans_I"].mean().unstack()

# GeoMx GSE301163 ROI means by histology
geomx_means = df.groupby("histology")[
    ["TROP2","DM1_axis","RAI_8","HLA_II","TLS","Cell_cycle","EMT","Hypoxia","Eight_gene_DM1"]
].mean()

fig = plt.figure(figsize=(16, 9))
g106 = gs.GridSpec(2, 2, hspace=0.45, wspace=0.30, top=0.93, bottom=0.06, left=0.07, right=0.97)

# Panel A: Visium Moran I means PT→ATC trajectory
axA = fig.add_subplot(g106[0,0])
v_axes = ["TROP2","DM1_like_score","RAI_8_score","HLA_II","TLS","Cell_cycle"]
v_axes_present = [a for a in v_axes if a in visium_means.columns]
for a in v_axes_present:
    pt_atc = ["PT","PTC","LPTC","ATC"]
    ys = [visium_means.loc[c,a] if c in visium_means.index else np.nan for c in pt_atc]
    axA.plot(pt_atc, ys, "o-", lw=1.6, ms=6, label=a)
axA.set_xlabel("Visium GSE250521 stage", fontsize=10); axA.set_ylabel("mean Moran's I (post-fix)", fontsize=10)
axA.set_title("A. Visium GSE250521 — niche organization trajectory (Moran's I)", fontsize=10.5, fontweight="bold")
axA.legend(fontsize=8, loc="upper left", ncol=2, frameon=False)
axA.grid(True, ls=":", alpha=0.4); axA.axhline(0, color="#888", lw=0.6, ls="--")

# Panel B: GeoMx ROI means histology trajectory
axB = fig.add_subplot(g106[0,1])
g_axes = ["TROP2","DM1_axis","RAI_8","HLA_II","TLS","Cell_cycle"]
hist_traj = [h for h in HIST_ORDER if h in geomx_means.index]
for a in g_axes:
    if a not in geomx_means.columns: continue
    ys = [geomx_means.loc[h,a] if h in geomx_means.index else np.nan for h in hist_traj]
    axB.plot(range(len(hist_traj)), ys, "o-", lw=1.6, ms=6, label=a)
axB.set_xticks(range(len(hist_traj))); axB.set_xticklabels(hist_traj, rotation=30, ha="right", fontsize=8.5)
axB.set_ylabel("ROI-mean log1p (normalized count)", fontsize=10)
axB.set_title("B. GeoMx GSE301163 — ROI-mean trajectory by histology (cross-platform)", fontsize=10.5, fontweight="bold")
axB.legend(fontsize=8, loc="upper right", ncol=2, frameon=False)
axB.grid(True, ls=":", alpha=0.4)

# Panel C: TROP2 specifically — Visium Moran I (per condition) vs GeoMx ROI mean (per histology)
axC = fig.add_subplot(g106[1,0])
v_trop = visium_g521[visium_g521.axis=="TROP2"].groupby("condition")["morans_I"].mean()
v_order = ["PT","PTC","LPTC","ATC"]
xs = np.arange(len(v_order))
ys = [v_trop.get(c, np.nan) for c in v_order]
axC.bar(xs, ys, color=["#3C6B4F","#34547A","#B8893C","#7B1F2A"], edgecolor="black", alpha=0.85)
axC.set_xticks(xs); axC.set_xticklabels(v_order, fontweight="bold")
axC.set_ylabel("Moran's I (Visium hex 6-nbr)", fontsize=10)
axC.set_title("C. Visium TROP2 niche Moran's I — confirmed PTC+LPTC tumor-specific", fontsize=10.5, fontweight="bold")
axC.axhline(0.13, color="#888", lw=0.6, ls="--", label="threshold")
axC.grid(True, ls=":", alpha=0.4, axis="y")

axD = fig.add_subplot(g106[1,1])
trop_panck = df[df.label=="PanCK+"].groupby("histology")["TROP2"].mean()
trop_vim = df[df.label=="VIM+"].groupby("histology")["TROP2"].mean()
hist_show = [h for h in HIST_ORDER if h in trop_panck.index or h in trop_vim.index]
xs = np.arange(len(hist_show))
panck_y = [trop_panck.get(h, np.nan) for h in hist_show]
vim_y = [trop_vim.get(h, np.nan) for h in hist_show]
w = 0.36
axD.bar(xs - w/2, panck_y, width=w, color="#34547A", edgecolor="black", alpha=0.85, label="PanCK+ epithelial")
axD.bar(xs + w/2, vim_y, width=w, color="#B8893C", edgecolor="black", alpha=0.85, label="VIM+ stromal")
axD.set_xticks(xs); axD.set_xticklabels([h[:18] for h in hist_show], rotation=30, ha="right", fontsize=8.5)
axD.set_ylabel("TROP2 ROI-mean log1p", fontsize=10)
axD.set_title("D. GeoMx TROP2 ROI-mean — PanCK+ vs VIM+ across histology", fontsize=10.5, fontweight="bold")
axD.legend(fontsize=9, frameon=False)
axD.grid(True, ls=":", alpha=0.4, axis="y")

fig.suptitle("S_F106  Cross-platform consistency — Visium spot-level Moran's I (n=28) ↔ GeoMx DSP ROI-level mean (n=78)\n"
             "Hypothesis (signature-level only): TROP2 enriched in PanCK+ epithelial of microPTC/PDTC consistent with Visium tumor-specific niche claim",
             fontsize=12.0, fontweight="bold", y=0.995)
fig.savefig(ASSETS/"S_F106_cross_platform_summary.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------- Summary stats TSV ----------
print("[stats] cross-platform comparison stats")
stats_recs = []
# 1. TROP2 PanCK+ vs VIM+ in tumor histologies
for h in ["microPTC","PDTC","Normal thyroid","Non-neoplastic nodular area"]:
    p = df[(df.histology==h)&(df.label=="PanCK+")]["TROP2"].values
    v = df[(df.histology==h)&(df.label=="VIM+")]["TROP2"].values
    _, pp = mw_p(p, v)
    stats_recs.append({"comparison": f"TROP2 PanCK+ vs VIM+ ({h})",
                       "n_panck": len(p), "n_vim": len(v),
                       "panck_mean": float(np.nanmean(p)) if len(p) else np.nan,
                       "vim_mean": float(np.nanmean(v)) if len(v) else np.nan,
                       "MW_p": pp})
# 2. TROP2 microPTC PanCK+ vs Normal PanCK+
for ax_name in ["TROP2","DM1_axis","RAI_8","Eight_gene_DM1","HLA_II","Cell_cycle"]:
    p = df[(df.histology=="microPTC")&(df.label=="PanCK+")][ax_name].values
    n = df[(df.histology=="Normal thyroid")&(df.label=="PanCK+")][ax_name].values
    _, pp = mw_p(p, n)
    stats_recs.append({"comparison": f"{ax_name} microPTC PanCK+ vs Normal PanCK+",
                       "n_panck": len(p), "n_vim": len(n),
                       "panck_mean": float(np.nanmean(p)) if len(p) else np.nan,
                       "vim_mean": float(np.nanmean(n)) if len(n) else np.nan,
                       "MW_p": pp})
# 3. PDTC vs microPTC differentiation
for ax_name in ["DM1_axis","RAI_8","Thyroid_TF"]:
    a = df[(df.histology=="PDTC")&(df.label=="PanCK+")][ax_name].values
    b = df[(df.histology=="microPTC")&(df.label=="PanCK+")][ax_name].values
    _, pp = mw_p(a, b)
    stats_recs.append({"comparison": f"{ax_name} PDTC PanCK+ vs microPTC PanCK+",
                       "n_panck": len(a), "n_vim": len(b),
                       "panck_mean": float(np.nanmean(a)) if len(a) else np.nan,
                       "vim_mean": float(np.nanmean(b)) if len(b) else np.nan,
                       "MW_p": pp})
stats_df = pd.DataFrame(stats_recs)
stats_df.to_csv(OUT/"GSE301163_crossplatform_stats.tsv", sep="\t", index=False)
print(stats_df.round(3).to_string(index=False))
print(f"[write] {OUT/'GSE301163_crossplatform_stats.tsv'}")

print("\n=== DONE — 6 figures (S_F101 ~ S_F106) generated ===")
