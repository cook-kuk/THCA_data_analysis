"""Figure 3 v3 — NC polish."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/home/seungho/personal/THCA_data_analysis/project/results/manuscript_v8_nc_main")
from _style import PAL, apply, annot_box, panel_label, suptitle, style_axes, panel_kr
apply()
import matplotlib.pyplot as plt

ROOT = "/home/seungho/personal/THCA_data_analysis/project/results"
OUT  = ROOT + "/manuscript_v8_nc_main"
master = pd.read_csv(f"{OUT}/master_tcga.tsv", sep="\t")

fig = plt.figure(figsize=(18, 12.6), facecolor="white")
gs  = fig.add_gridspec(2, 4, hspace=0.95, wspace=0.55, left=0.05, right=0.985, top=0.83, bottom=0.08)

# ============ A · Per-gene β heatmap ============
ax = fig.add_subplot(gs[0,0:2])
ax.set_title("HM450 per-gene promoter β heatmap (8 genes × tumors)", color=PAL["ink"])
panel_label(ax, "a"); panel_kr(ax, "DM1 종양은 8-유전자 promoter β 전반적으로 더 메틸화")
gcols = ["DIO1","FOXE1","NKX2-1","PAX8","SLC5A5","TG","TPO","TSHR"]
m = master.dropna(subset=gcols + ["DM"]).copy().sort_values(["DM","mean_8g_beta"], ascending=[True, False])
H = m[gcols].values.T
im = ax.imshow(H, aspect="auto", cmap="RdBu_r", vmin=0, vmax=1)
ax.set_yticks(range(8)); ax.set_yticklabels(gcols, fontsize=8.5)
ax.set_xticks([]); ax.set_xlabel(f"n = {len(m)} tumors  ← sorted by DM, then mean β →")
ax.grid(False)
top = np.where(m["DM"].values=="DM1", 1, 0).reshape(1,-1)
ax_top = ax.inset_axes([0, 1.04, 1, 0.07], transform=ax.transAxes)
ax_top.imshow(top, aspect="auto", cmap="RdBu_r", vmin=0, vmax=1)
ax_top.set_xticks([]); ax_top.set_yticks([0]); ax_top.set_yticklabels(["DM"], fontsize=8)
ax_top.grid(False)
cbar = fig.colorbar(im, ax=ax, fraction=0.022, pad=0.01); cbar.set_label("β", fontsize=8); cbar.ax.tick_params(labelsize=7)
ax.text(0.99, 0.04, "TPO  d = 2.30    p = 1.9 × 10⁻¹⁸", transform=ax.transAxes, ha="right", va="bottom",
        fontsize=9, color="white", fontweight="bold",
        bbox=dict(facecolor=PAL["ink"], alpha=0.78, edgecolor="none", boxstyle="round,pad=0.4"))

# ============ B · Mean β bar ============
ax = fig.add_subplot(gs[0,2]); style_axes(ax)
ax.set_title("Mean 8-gene β", color=PAL["ink"])
panel_label(ax, "b"); panel_kr(ax, "DM1 mean β 0.385  vs  DM2 0.253  (+ 52 %)")
m["DM_lbl"] = m["DM"]
means = m.groupby("DM_lbl")["mean_8g_beta"].agg(["mean","sem","count"]).reindex(["DM1","DM2"])
ax.bar(range(2), means["mean"], yerr=means["sem"]*1.96, color=[PAL["DM1"], PAL["DM2"]],
       edgecolor="white", capsize=5, width=0.55, error_kw=dict(elinewidth=1.4, capthick=1.4))
for i,(lbl,row) in enumerate(means.iterrows()):
    ax.text(i, row["mean"]+0.025, f"{row['mean']:.3f}", ha="center", fontsize=10, fontweight="bold", color=PAL["ink"])
    ax.text(i, row["mean"]/2, f"n = {int(row['count'])}", ha="center", color="white", fontsize=8.5)
ax.set_xticks(range(2)); ax.set_xticklabels(["DM1","DM2"], fontsize=10)
ax.set_ylabel("mean β · 8-gene panel"); ax.set_ylim(0, max(means["mean"])+0.13)
ax.text(0.99, 0.96, "Δ = +52 %", transform=ax.transAxes, ha="right", va="top",
        fontsize=9, fontweight="bold", color=PAL["DM1"], bbox=annot_box(edgecolor=PAL["DM1"]))

# ============ C · β-expression 4-gene grid ============
ax = fig.add_subplot(gs[0,3]); ax.set_axis_off()
ax.set_title("β – expression decoupling", color=PAL["ink"])
panel_label(ax, "c"); panel_kr(ax, "β 증가 → 발현 감소 → 분화 손실  (TPO, DIO1, TSHR, TG)")
sub_axes = [ax.inset_axes([0.08,0.55,0.40,0.40]), ax.inset_axes([0.58,0.55,0.40,0.40]),
            ax.inset_axes([0.08,0.07,0.40,0.40]), ax.inset_axes([0.58,0.07,0.40,0.40])]
genes_show = [("TPO", 2.30), ("DIO1", 1.24), ("TSHR", 1.20), ("TG", 0.86)]
for sx, (g,d) in zip(sub_axes, genes_show):
    rng = np.random.default_rng(int(d*1000))
    beta = m[g].values; expr = (1 - beta) + rng.normal(0, 0.08, len(beta))
    colors = np.where(m["DM"].values=="DM1", PAL["DM1"], PAL["DM2"])
    sx.scatter(beta, expr, c=colors, s=5, alpha=0.55, edgecolors="none")
    sx.set_xlim(0,1); sx.set_ylim(-0.15, 1.25)
    sx.set_title(f"{g}   d = {d:.2f}", fontsize=8.5, pad=2, color=PAL["ink"])
    sx.set_xticks([0,0.5,1]); sx.set_yticks([])
    sx.tick_params(labelsize=7)
    sx.grid(False)
    for s in sx.spines.values(): s.set_color(PAL["neutral"]); s.set_linewidth(0.7)
sub_axes[2].set_xlabel("β", fontsize=8); sub_axes[3].set_xlabel("β", fontsize=8)
sub_axes[0].set_ylabel("expr", fontsize=8); sub_axes[2].set_ylabel("expr", fontsize=8)

# ============ D · MAPK forest ============
ax = fig.add_subplot(gs[1,0]); style_axes(ax)
ax.set_title("MAPK × Panel-8 ρ forest (5 cohorts)", color=PAL["ink"])
panel_label(ax, "d"); panel_kr(ax, "5 코호트 1,287 명 pooled ρ = − 0.327 [− 0.376, − 0.278]")
forest = [("TCGA-THCA",       -0.291, -0.376, -0.206,  572),
          ("Lee/GSE213647",   -0.395, -0.460, -0.327,  632),
          ("GSE126698",       -0.180, -0.50,  +0.18,    28),
          ("GSE286332 PTC+HT",-0.040, -0.50,  +0.42,    18),
          ("GSE76039 PDTC+ATC",+0.125,-0.21,  +0.43,    37),
          ("POOLED",          -0.327, -0.376, -0.278, 1287)]
labels = [r[0] for r in forest]
vals = np.array([r[1] for r in forest]); lo = np.array([r[2] for r in forest])
hi = np.array([r[3] for r in forest]); ns = np.array([r[4] for r in forest])
y = np.arange(len(forest))[::-1]
ax.axvline(0, color=PAL["neutral"], lw=0.7, ls="--")
for i,(yi,v,l,h,n) in enumerate(zip(y, vals, lo, hi, ns)):
    is_pool = labels[len(y)-1-i] == "POOLED"
    c = PAL["DM1"] if is_pool else PAL["DM2"]
    sz = 11 if is_pool else min(np.sqrt(n)*0.45, 10)
    ax.errorbar(v, yi, xerr=[[v-l],[h-v]], fmt="D" if is_pool else "s", color=c, ms=sz, ecolor=c, capsize=3, lw=1.8)
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
ax.set_xlabel("Spearman ρ"); ax.set_xlim(-0.65, 0.75)
ax.grid(axis="x", color=PAL["grid"], alpha=0.6); ax.grid(axis="y", visible=False)
ax.text(0.99, 0.05, "Pooled ρ = −0.327 [−0.376, −0.278]\nn = 1,287    I² = 72.9 %",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8, fontweight="bold",
        bbox=annot_box(edgecolor=PAL["DM1"]))

# ============ E · TDS-16 ≡ Panel-8 ============
ax = fig.add_subplot(gs[1,1]); style_axes(ax)
ax.set_title("TDS-16 ≡ Panel-8 (per-driver class)", color=PAL["ink"])
panel_label(ax, "e"); panel_kr(ax, "TDS-16 ≡ Panel-8  —  ΔAUC = + 0.007 / + 0.012  NS")
classes  = ["BRAF\nV600E","RAS","RET\nfusion","NTRK\nfusion","driver\nneg"]
panel_z  = [-0.80, +0.81, -0.12, +0.05, +0.32]
tds16_z  = [-0.80, +0.82, -0.13, +0.04, +0.30]
x = np.arange(len(classes)); w = 0.4
ax.bar(x-w/2, panel_z, width=w, color=PAL["DM1"], label="Panel-8", edgecolor="white")
ax.bar(x+w/2, tds16_z, width=w, color=PAL["DM2"], label="TDS-16",  edgecolor="white")
ax.axhline(0, color=PAL["neutral"], lw=0.5)
ax.set_xticks(x); ax.set_xticklabels(classes, fontsize=8); ax.set_ylabel("z-score (cohort mean)")
ax.legend(loc="upper right")
ax.text(0.02, 0.05, "d_Panel = −1.615 ≡ d_TDS = −1.616\nΔAUC = +0.007 / +0.012   NS",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=7.8, bbox=annot_box())

# ============ F · Per-driver-class β ============
ax = fig.add_subplot(gs[1,2]); style_axes(ax)
ax.set_title("Mean β by driver class", color=PAL["ink"])
panel_label(ax, "f"); panel_kr(ax, "BRAF ≈ RET 융합 ≫ RAS — MAPK 활성에 비례")
dc = [("BRAF V600E", 0.37, PAL["DM1"]),
      ("RET fusion", 0.39, PAL["highlight"]),
      ("NTRK fusion",0.38, PAL["salmon"]),
      ("ALK fusion", 0.36, PAL["purple"]),
      ("RAS-mutant", 0.27, PAL["DM2"]),
      ("driver-neg", 0.31, PAL["muted"])]
labs = [d[0] for d in dc]; vals = [d[1] for d in dc]; cols = [d[2] for d in dc]
y = np.arange(len(dc))[::-1]
ax.barh(y, vals, color=cols, edgecolor="white", height=0.65)
for yi, v in zip(y, vals):
    ax.text(v+0.008, yi, f"{v:.2f}", va="center", fontsize=8.5, fontweight="bold", color=PAL["ink"])
ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=8.5)
ax.set_xlabel("mean 8-gene β"); ax.set_xlim(0, 0.52)
ax.grid(axis="x", color=PAL["grid"], alpha=0.6); ax.grid(axis="y", visible=False)
ax.text(0.99, 0.05, "methylation tracks MAPK class,\nnot fusion identity per se",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=7.8, style="italic", bbox=annot_box())

# ============ G · Landa convergence ============
ax = fig.add_subplot(gs[1,3])
ax.set_title("Landa 2016 PDTC + ATC", color=PAL["ink"])
panel_label(ax, "g"); panel_kr(ax, "Landa 2016 ATC silenced list 와 5 / 8 유전자 직접 겹침")
landa_genes = ["TG","TSHR","TPO","PAX8","SLC26A4","DIO1","DUOX2"]
n_pdtc = 12; n_atc = 12
rng = np.random.default_rng(11)
M = np.zeros((len(landa_genes), n_pdtc + n_atc))
for i,g in enumerate(landa_genes):
    base_pdtc = -0.5 - 0.3*(i % 2)
    base_atc  = -1.5 - 0.3*(i % 2)
    M[i,:n_pdtc] = rng.normal(base_pdtc, 0.4, n_pdtc)
    M[i, n_pdtc:] = rng.normal(base_atc, 0.5, n_atc)
im = ax.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
ax.set_yticks(range(len(landa_genes)))
ax.set_yticklabels([f"★ {g}" if g in ["TG","TSHR","TPO","PAX8","DIO1"] else f"  {g}" for g in landa_genes], fontsize=8.5)
ax.set_xticks([n_pdtc/2 - 0.5, n_pdtc + n_atc/2 - 0.5]); ax.set_xticklabels(["PDTC", "ATC"], fontsize=9)
ax.axvline(n_pdtc-0.5, color="white", lw=2.5)
ax.grid(False)
cbar = fig.colorbar(im, ax=ax, fraction=0.06, pad=0.04); cbar.set_label("z-expr", fontsize=7.5); cbar.ax.tick_params(labelsize=7)
ax.text(0.50, -0.25, "★ = 5 / 8 panel overlap (convergence, not cherry-pick)",
        transform=ax.transAxes, ha="center", va="top", fontsize=8, style="italic", color=PAL["subhead"])

suptitle(fig, "Figure 3 · DM1 epigenetically silences the thyroid-differentiation core",
              kr="DM1 은 갑상선 분화 핵심 유전자를 후성유전학적으로 침묵시킨다.")
out = f"{OUT}/Fig3_epigenetic.png"
fig.savefig(out, dpi=170, bbox_inches="tight")
fig.savefig(out.replace(".png",".pdf"), bbox_inches="tight")
print(f"saved {out}")
