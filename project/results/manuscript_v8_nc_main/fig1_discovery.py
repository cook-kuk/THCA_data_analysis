"""Figure 1 v3 — NC polish. Discovery axis on one page (7 panels)."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/home/seungho/personal/THCA_data_analysis/project/results/manuscript_v8_nc_main")
from _style import PAL, apply, annot_box, panel_label, suptitle, style_axes, panel_kr
apply()
import matplotlib.pyplot as plt
import matplotlib.patches as mp
from sklearn.decomposition import PCA
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

ROOT = "/home/seungho/personal/THCA_data_analysis/project/results"
OUT  = ROOT + "/manuscript_v8_nc_main"
master = pd.read_csv(f"{OUT}/master_tcga.tsv", sep="\t")

fig = plt.figure(figsize=(18, 12.6), facecolor="white")
gs  = fig.add_gridspec(2, 4, hspace=0.95, wspace=0.55, left=0.05, right=0.985, top=0.83, bottom=0.08)

# ============ A · Schematic ============
ax = fig.add_subplot(gs[0,0]); ax.set_axis_off()
ax.set_title("Multi-cohort study schematic", color=PAL["ink"])
panel_label(ax, "a"); panel_kr(ax, "6개 코호트 · 약 2,016 종양 · 6 modality 통합")
cohorts = [
    ("TCGA-THCA",      "n = 504  RNA + HM450 + survival",       PAL["DM2"]),
    ("MSK-IMPACT",     "n = 117  advanced disease",             PAL["highlight"]),
    ("GSE213647 (Lee)","n = 632  Korean FFPE",                  PAL["positive"]),
    ("PRJEB11591 (K2)","n = 260  Korean fresh-frozen",          PAL["DM1"]),
    ("HM450 TCGA",     "n = 503  methylation",                  PAL["purple"]),
    ("sc Pu + Lu",     "n = 6 + 23  single-cell",               PAL["teal"]),
]
y = 0.96
for name,sub,c in cohorts:
    ax.add_patch(mp.FancyBboxPatch((0.02, y-0.125), 0.96, 0.11, boxstyle="round,pad=0.014",
                                    fc=c, ec="none", alpha=0.92, transform=ax.transAxes))
    ax.text(0.06, y-0.07, name, fontsize=9, fontweight="bold", color="white", transform=ax.transAxes, va="center")
    ax.text(0.96, y-0.07, sub,  fontsize=7.5, ha="right", va="center", color="white", transform=ax.transAxes)
    y -= 0.155
ax.text(0.5, 0.02, "Σ ≈ 2,016 tumors  ·  6 layers", ha="center", fontsize=8.5, style="italic", color=PAL["subhead"], transform=ax.transAxes)

# ============ B · TIERA67 → 8-gene ============
ax = fig.add_subplot(gs[0,1]); ax.set_axis_off()
ax.set_title("TIERA67 → 8-gene panel (TDS-core)", color=PAL["ink"])
panel_label(ax, "b"); panel_kr(ax, "TIERA67 67-유전자 → 8-유전자 TDS-core 패널")
left_cats = [("Driver_anchor", 6),("TDS_core (panel)", 8),("RAI_machinery", 7),
             ("Lineage TF", 9),("MAPK_output", 12),("Immune_context", 18),("Dediff", 7)]
total = sum(n for _,n in left_cats)
y = 0.94
for name,n in left_cats:
    h = 0.78 * n/total
    fc = PAL["DM1"] if "TDS_core" in name else PAL["muted"]
    ax.add_patch(mp.Rectangle((0.03, y-h), 0.32, h, fc=fc, ec="white"))
    ax.text(0.05, y-h/2, name, va="center", fontsize=7.5, color="white", fontweight="bold")
    ax.text(0.33, y-h/2, f"({n})", va="center", ha="right", fontsize=7, color="white")
    if "TDS_core" in name:
        ax.annotate("", xy=(0.62, 0.54), xytext=(0.36, y-h/2),
                    arrowprops=dict(arrowstyle="->", lw=2.2, color=PAL["DM1"]))
    y -= h + 0.005
panel_g = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]
ax.add_patch(mp.FancyBboxPatch((0.62, 0.10), 0.34, 0.84, boxstyle="round,pad=0.012", fc=PAL["DM1"], ec="white", alpha=0.95))
ax.text(0.79, 0.92, "8-gene panel", ha="center", fontsize=8.5, fontweight="bold", color="white")
for i,g in enumerate(panel_g):
    ax.text(0.79, 0.84 - i*0.090, g, ha="center", fontsize=8.5, color="white")
ax.text(0.5, 0.02, "Yoo 2016 RAI prior  ·  reverse-causality lock", ha="center", fontsize=7.5, style="italic", color=PAL["subhead"])

# ============ C · UMAP DM1/DM2 ============
ax = fig.add_subplot(gs[0,2]); style_axes(ax)
ax.set_title("DM1/DM2 cluster on 8-gene β space", color=PAL["ink"])
panel_label(ax, "c"); panel_kr(ax, "PCA 공간에서 DM1 (붉음) / DM2 (파랑) 분리 명확")
gcols = ["DIO1","FOXE1","NKX2-1","PAX8","SLC5A5","TG","TPO","TSHR"]
sub = master.dropna(subset=gcols + ["DM"]).copy()
pc = PCA(n_components=2).fit_transform(sub[gcols].values)
colors = sub["DM"].map({"DM1": PAL["DM1"], "DM2": PAL["DM2"]}).values
ax.scatter(pc[:,0], pc[:,1], c=colors, s=15, alpha=0.78, edgecolors="white", linewidths=0.4)
ax.set_xlabel("PC1 (8-gene β space)"); ax.set_ylabel("PC2")
ax.legend(handles=[mp.Patch(color=PAL["DM1"], label=f"DM1 (n = {(sub['DM']=='DM1').sum()})"),
                   mp.Patch(color=PAL["DM2"], label=f"DM2 (n = {(sub['DM']=='DM2').sum()})")],
          loc="upper left", bbox_to_anchor=(0.02, 0.98))
ax.text(0.98, 0.04, "DM1 prevalence  28.4 %\nKMeans  k = 2",
        transform=ax.transAxes, ha="right", va="bottom",
        fontsize=8, color=PAL["ink"], bbox=annot_box())

# ============ D · Heatmap sorted by P_DM1 ============
ax = fig.add_subplot(gs[0,3])
ax.set_title("8-gene β heatmap sorted by P_DM1", color=PAL["ink"])
panel_label(ax, "d"); panel_kr(ax, "DM1 종양에서 8-유전자 promoter β 일관되게 높음")
sub2 = master.dropna(subset=gcols + ["DM","mean_8g_beta"]).copy().sort_values("mean_8g_beta", ascending=False)
H = sub2[gcols].values.T
im = ax.imshow(H, aspect="auto", cmap="RdBu_r", vmin=0, vmax=1)
ax.set_yticks(range(8)); ax.set_yticklabels(gcols, fontsize=8)
ax.set_xticks([]); ax.set_xlabel(f"tumors sorted by mean β  (n = {len(sub2)})")
ax.grid(False)
top_track = np.where(sub2["DM"].values=="DM1", 1, 0).reshape(1,-1)
ax_top = ax.inset_axes([0, 1.04, 1, 0.07], transform=ax.transAxes)
ax_top.imshow(top_track, aspect="auto", cmap="RdBu_r", vmin=0, vmax=1)
ax_top.set_xticks([]); ax_top.set_yticks([0]); ax_top.set_yticklabels(["DM"], fontsize=7.5)
ax_top.grid(False)
cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02); cbar.set_label("β", fontsize=8); cbar.ax.tick_params(labelsize=7)

# ============ E · Driver mRNA neutrality ============
ax = fig.add_subplot(gs[1,0]); style_axes(ax)
ax.set_title("Driver mRNA neutrality (AUC ≈ 0.5)", color=PAL["ink"])
panel_label(ax, "e"); panel_kr(ax, "BRAF / TERT / RAS 모두 AUC ≈ 0.5 — driver-neutral")
drivers = ["BRAF","TERT","KRAS","NRAS","HRAS"]
aucs    = [0.602, 0.578, 0.525, 0.521, 0.500]
xpos = np.arange(5)
ax.bar(xpos, aucs, color=PAL["neutral"], edgecolor="white", width=0.6)
ax.axhline(0.5, color=PAL["DM1"], lw=1.4, ls="--", label="chance")
for i,a in enumerate(aucs):
    ax.text(i, a+0.013, f"{a:.2f}", ha="center", fontsize=8, fontweight="bold", color=PAL["ink"])
ax.set_xticks(xpos); ax.set_xticklabels(drivers, fontsize=8)
ax.set_ylim(0, 0.78); ax.set_ylabel("DM1 vs DM2 single-feature AUC")
ax.text(0.02, 0.96, "BRAF d = −0.04  ·  MW p = 0.57",
        transform=ax.transAxes, ha="left", va="top",
        fontsize=8, color=PAL["DM1"], fontweight="bold", bbox=annot_box(edgecolor=PAL["DM1"]))
ax.legend(loc="upper right")

# ============ F · Pan-genome ARI ladder ============
ax = fig.add_subplot(gs[1,1]); style_axes(ax)
ax.set_title("Pan-genome ARI ladder", color=PAL["ink"])
panel_label(ax, "f"); panel_kr(ax, "8-유전자 (0.49) → TIERA67 (0.90) ≈ pan-genome (0.92)")
ladder = [("Driver_anchor",-0.007),("8-gene",0.49),("16-gene",0.62),("TIERA67",0.90),
          ("top-1000",0.91),("top-5000",0.92)]
labs = [n for n,_ in ladder]; vals = [v for _,v in ladder]
cols = [PAL["highlight"] if "Driver" in n else PAL["DM2"] if v<0.5 else PAL["positive"] for n,v in ladder]
ax.barh(range(len(ladder)), vals, color=cols, edgecolor="white", height=0.66)
for i,v in enumerate(vals):
    ax.text(max(v,0)+0.025, i, f"{v:.2f}", va="center", fontsize=8, color=PAL["ink"])
ax.axvline(0, color=PAL["neutral"], lw=0.6)
ax.set_yticks(range(len(ladder))); ax.set_yticklabels(labs, fontsize=8); ax.invert_yaxis()
ax.set_xlabel("ARI vs 8-gene reference"); ax.set_xlim(-0.12, 1.10)
ax.grid(axis="x", color=PAL["grid"], alpha=0.7); ax.grid(axis="y", visible=False)
ax.text(0.99, 0.04, "TIERA67 enrichment\np = 3 × 10⁻⁴", transform=ax.transAxes, ha="right", va="bottom",
        fontsize=7.5, style="italic", bbox=annot_box())

# ============ G · KM teaser ============
ax = fig.add_subplot(gs[1,2:]); style_axes(ax)
ax.set_title("Kaplan-Meier overall survival (TCGA primary)", color=PAL["ink"])
panel_label(ax, "g"); panel_kr(ax, "DM1 환자 생존 곡선이 DM2 대비 더 가파른 하강")
surv = master.dropna(subset=["os_days","os_event","DM"]).copy()
surv["os_event"] = surv["os_event"].astype(int)
km_dm1 = surv[surv["DM"]=="DM1"]; km_dm2 = surv[surv["DM"]=="DM2"]
kmf = KaplanMeierFitter()
kmf.fit(km_dm1["os_days"]/365.25, km_dm1["os_event"], label=f"DM1  (n = {len(km_dm1)}, e = {int(km_dm1['os_event'].sum())})")
kmf.plot(ax=ax, ci_show=False, color=PAL["DM1"], lw=2.6)
kmf.fit(km_dm2["os_days"]/365.25, km_dm2["os_event"], label=f"DM2  (n = {len(km_dm2)}, e = {int(km_dm2['os_event'].sum())})")
kmf.plot(ax=ax, ci_show=False, color=PAL["DM2"], lw=2.6)
lr = logrank_test(km_dm1["os_days"], km_dm2["os_days"], km_dm1["os_event"], km_dm2["os_event"])
ax.set_xlabel("Years from diagnosis"); ax.set_ylabel("Overall survival probability"); ax.set_ylim(0.0, 1.02)
ax.text(0.99, 0.06, f"log-rank p = {lr.p_value:.3f}\nχ² = {lr.test_statistic:.2f}\nHR (meta) cross-ref → Fig 5a",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
        bbox=annot_box())
ax.legend(loc="lower left")

suptitle(fig, "Figure 1 · The DM1/DM2 axis defines an aggressive PTC subtype orthogonal to canonical BRAF/RAS classification",
              kr="DM1/DM2 분자 축은 기존 BRAF/RAS 분류 체계와 직교하는 공격성 갑상선 유두암 아형을 정의한다.")
out = f"{OUT}/Fig1_discovery_axis.png"
fig.savefig(out, dpi=170, bbox_inches="tight")
fig.savefig(out.replace(".png",".pdf"), bbox_inches="tight")
print(f"saved {out}")
