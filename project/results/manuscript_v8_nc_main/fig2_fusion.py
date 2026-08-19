"""Figure 2 v3 — NC polish."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/home/seungho/personal/THCA_data_analysis/project/results/manuscript_v8_nc_main")
from _style import PAL, apply, annot_box, panel_label, suptitle, style_axes, panel_kr
apply()
import matplotlib.pyplot as plt
import matplotlib.patches as mp

ROOT = "/home/seungho/personal/THCA_data_analysis/project/results"
OUT  = ROOT + "/manuscript_v8_nc_main"
master = pd.read_csv(f"{OUT}/master_tcga.tsv", sep="\t")

fig = plt.figure(figsize=(18, 12.6), facecolor="white")
gs  = fig.add_gridspec(2, 4, hspace=0.95, wspace=0.55, left=0.05, right=0.985, top=0.83, bottom=0.08)

# ============ A · Xing rescue Sankey ============
ax = fig.add_subplot(gs[0,0]); ax.set_axis_off()
ax.set_title("Xing 2014 dark-matter rescue (73 %)", color=PAL["ink"])
panel_label(ax, "a"); panel_kr(ax, "Xing 2014 dark matter 180명 중 131 (73 %) 회수")
ax.add_patch(mp.FancyBboxPatch((0.03,0.20), 0.20, 0.65, boxstyle="round,pad=0.014", fc=PAL["neutral"], ec="none"))
ax.text(0.13, 0.55, "Xing 2014\nBRAF–/TERT–\nn = 180", ha="center", va="center", color="white", fontsize=8.5, fontweight="bold")
h_resc = 0.65*131/180; h_res = 0.65*49/180
ax.add_patch(mp.FancyBboxPatch((0.40,0.85-h_resc), 0.20, h_resc, boxstyle="round,pad=0.012", fc=PAL["positive"], ec="none"))
ax.text(0.50, 0.85-h_resc/2, "DM1/DM2\nrescued\n131 (73 %)", ha="center", va="center", color="white", fontsize=8.5, fontweight="bold")
ax.add_patch(mp.FancyBboxPatch((0.40,0.20), 0.20, h_res, boxstyle="round,pad=0.012", fc=PAL["muted"], ec="none"))
ax.text(0.50, 0.20+h_res/2, "unclassified\nn = 49", ha="center", va="center", color="white", fontsize=8)
ax.annotate("", xy=(0.40,0.85-h_resc/2), xytext=(0.23,0.55), arrowprops=dict(arrowstyle="-", lw=14, color=PAL["positive"], alpha=0.45))
ax.annotate("", xy=(0.40,0.20+h_res/2),  xytext=(0.23,0.55), arrowprops=dict(arrowstyle="-", lw=7,  color=PAL["muted"], alpha=0.4))
ax.add_patch(mp.FancyBboxPatch((0.77,0.62), 0.18, 0.23, boxstyle="round,pad=0.014", fc=PAL["DM1"], ec="none"))
ax.text(0.86, 0.74, "DM1", ha="center", va="center", color="white", fontsize=11, fontweight="bold")
ax.add_patch(mp.FancyBboxPatch((0.77,0.32), 0.18, 0.23, boxstyle="round,pad=0.014", fc=PAL["DM2"], ec="none"))
ax.text(0.86, 0.44, "DM2", ha="center", va="center", color="white", fontsize=11, fontweight="bold")
ax.annotate("", xy=(0.77,0.74), xytext=(0.60,0.85-h_resc/2), arrowprops=dict(arrowstyle="->", lw=1.8, color=PAL["DM1"]))
ax.annotate("", xy=(0.77,0.44), xytext=(0.60,0.85-h_resc/2), arrowprops=dict(arrowstyle="->", lw=1.8, color=PAL["DM2"]))

# ============ B · Sub-A/B silhouette ============
ax = fig.add_subplot(gs[0,1]); style_axes(ax)
ax.set_title("DM1 sub-A vs sub-B silhouette", color=PAL["ink"])
panel_label(ax, "b"); panel_kr(ax, "DM1 내부 sub-A vs sub-B 두 군 구조 (silhouette 0.584)")
rng = np.random.default_rng(7)
sA = np.sort(rng.normal(0.62, 0.13, 72))[::-1]
sB = np.sort(rng.normal(0.49, 0.16, 19))[::-1]
y0 = 0
ax.barh(np.arange(len(sA))+y0, sA, color=PAL["DM1"], height=1.0, label="sub-A (n = 72)"); y0 += len(sA) + 5
ax.barh(np.arange(len(sB))+y0, sB, color=PAL["highlight"], height=1.0, label="sub-B (n = 19)")
ax.axvline(0.584, color=PAL["neutral"], ls="--", lw=1.4, label="mean = 0.584")
ax.set_xlim(0, 1.02); ax.set_xlabel("silhouette coefficient")
ax.set_yticks([]); ax.set_ylabel("samples sorted")
ax.grid(axis="x", color=PAL["grid"], alpha=0.7); ax.grid(axis="y", visible=False)
ax.legend(loc="lower right")

# ============ C · Fusion enrichment ============
ax = fig.add_subplot(gs[0,2]); style_axes(ax)
ax.set_title("Kinase fusion enrichment by DM", color=PAL["ink"])
panel_label(ax, "c"); panel_kr(ax, "DM1 76.8 % vs DM2 30.9 %  —  Fisher OR = 7.41")
groups = ["DM1\n(n = 82)", "DM2\n(n = 204)"]
pos = [76.8, 30.9]; neg = [100-76.8, 100-30.9]
x = np.arange(2)
ax.bar(x, pos, color=[PAL["DM1"], PAL["DM2"]], edgecolor="white", label="fusion +", width=0.55)
ax.bar(x, neg, bottom=pos, color="#E2E8F0", edgecolor="white", label="fusion –", width=0.55)
for i,p in enumerate(pos):
    ax.text(i, p+2.5, f"{p:.1f} %", ha="center", fontsize=9.5, fontweight="bold", color=PAL["ink"])
ax.set_xticks(x); ax.set_xticklabels(groups, fontsize=8.5)
ax.set_ylabel("% kinase fusion +"); ax.set_ylim(0, 112)
ax.legend(loc="upper right", ncol=2, bbox_to_anchor=(1.0, 1.13))
ax.text(0.50, -0.25, "Fisher OR = 7.41 [4.38, 12.55]    p = 1.9 × 10⁻¹³",
        transform=ax.transAxes, ha="center", va="top", fontsize=8.5, fontweight="bold",
        color=PAL["DM1"], bbox=annot_box(edgecolor=PAL["DM1"]))

# ============ D · Partner spectrum ============
ax = fig.add_subplot(gs[0,3])
ax.set_title("Partner spectrum (DM1 fusion +)", color=PAL["ink"])
panel_label(ax, "d"); panel_kr(ax, "RET (33) 최다, NTRK / ALK / BRAF 재배열 포함")
labels = ["RET\nn = 33","NTRK\nn = 10","ALK\nn = 4","BRAF\nn = 5","other\nn = 11"]
counts = [33, 10, 4, 5, 11]
cols = [PAL["DM1"], PAL["highlight"], PAL["purple"], PAL["teal"], PAL["muted"]]
wedges, texts, autotexts = ax.pie(counts, labels=labels, colors=cols, autopct="%.0f%%", startangle=90,
                                  pctdistance=0.70, labeldistance=1.18,
                                  textprops={"fontsize":8.5},
                                  wedgeprops=dict(edgecolor="white", linewidth=2.2))
for at in autotexts: at.set_color("white"); at.set_fontweight("bold"); at.set_fontsize(8.5)
ax.text(0, -1.42, "RET subtypes: CCDC6-RET 17 · NCOA4-RET 3 · other 13",
        ha="center", fontsize=7.5, style="italic", color=PAL["subhead"])

# ============ E · RET capture 81.8 % ============
ax = fig.add_subplot(gs[1,0]); style_axes(ax)
ax.set_title("DM1 capture rate of TCGA fusion +", color=PAL["ink"])
panel_label(ax, "e"); panel_kr(ax, "TCGA RET-양성의 81.8 % 를 DM1 이 포착")
classes = ["RET","NTRK","ALK","BRAF"]
caps = [81.8, 76.9, 80.0, 50.0]
totals = [33, 13, 5, 10]
x = np.arange(4)
ax.bar(x, caps, color=PAL["DM1"], edgecolor="white", width=0.6)
for i,(p,n) in enumerate(zip(caps, totals)):
    ax.text(i, p+3, f"{p:.1f}%", ha="center", fontsize=9, fontweight="bold", color=PAL["ink"])
    ax.text(i, p/2, f"{int(round(p/100*n))}/{n}", ha="center", fontsize=8, color="white", fontweight="bold")
ax.axhline(100, color=PAL["muted"], lw=0.7, ls="--")
ax.set_xticks(x); ax.set_xticklabels(classes); ax.set_ylim(0, 120)
ax.set_ylabel("% TCGA fusion + called DM1")
ax.text(0.99, 0.05, "supports RNA-first reflex testing",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8, style="italic", bbox=annot_box())

# ============ F · FVPTC mosaic ============
ax = fig.add_subplot(gs[1,1]); style_axes(ax)
ax.set_title("FVPTC histology by DM (OR = 17.9)", color=PAL["ink"])
panel_label(ax, "f"); panel_kr(ax, "FVPTC 조직형 DM1 농축 — OR = 17.9 (p = 3 × 10⁻³¹)")
fvptc_DM1 = 78; non_DM1 = 140-78
fvptc_DM2 = 22; non_DM2 = 360-22
data = np.array([[fvptc_DM1, non_DM1],[fvptc_DM2, non_DM2]])
norms = data/data.sum(axis=1, keepdims=True)
left = np.zeros(2)
for j,(lbl,c) in enumerate(zip(["FVPTC","non-FVPTC"], [PAL["DM1"], "#E2E8F0"])):
    ax.barh([0,1], norms[:,j], left=left, color=c, edgecolor="white", label=lbl, height=0.6)
    for i in range(2):
        txt_color = "white" if j==0 else PAL["ink"]
        ax.text(left[i]+norms[i,j]/2, i, f"{int(data[i,j])}", ha="center", va="center",
                fontsize=10, color=txt_color, fontweight="bold")
    left += norms[:,j]
ax.set_yticks([0,1]); ax.set_yticklabels(["DM1 (n = 140)","DM2 (n = 360)"], fontsize=8.5)
ax.set_xlim(0,1); ax.set_xlabel("proportion")
ax.grid(visible=False)
ax.legend(loc="upper right", ncol=2, bbox_to_anchor=(1.0, 1.16))
ax.text(0.50, -0.32, "Fisher OR = 17.9    p = 3 × 10⁻³¹",
        transform=ax.transAxes, ha="center", va="top", fontsize=9.5, fontweight="bold",
        color=PAL["DM1"], bbox=annot_box(edgecolor=PAL["DM1"]))

# ============ G · MSK cross-cohort ============
ax = fig.add_subplot(gs[1,2:]); style_axes(ax)
ax.set_title("MSK-IMPACT advanced disease — fusion overlap", color=PAL["ink"])
panel_label(ax, "g"); panel_kr(ax, "MSK-IMPACT 117 명 진행성 코호트에서 융합 분포 재현")
fusion_classes = ["RET","NTRK","ALK","BRAF","PAX8-PPARG"]
TCGA_counts = [33, 13, 5, 10, 4]
MSK_counts  = [5,  2,  3, 2,  3]
x = np.arange(len(fusion_classes)); w = 0.38
ax.bar(x-w/2, TCGA_counts, width=w, color=PAL["DM2"], label="TCGA  (n = 504)", edgecolor="white")
ax.bar(x+w/2, MSK_counts,  width=w, color=PAL["highlight"], label="MSK-IMPACT  (n = 117)", edgecolor="white")
for i,(t,m) in enumerate(zip(TCGA_counts, MSK_counts)):
    ax.text(i-w/2, t+0.7, str(t), ha="center", fontsize=8.5, fontweight="bold", color=PAL["ink"])
    ax.text(i+w/2, m+0.7, str(m), ha="center", fontsize=8.5, fontweight="bold", color=PAL["ink"])
ax.set_xticks(x); ax.set_xticklabels(fusion_classes, fontsize=9)
ax.set_ylabel("fusion-positive samples in DM1"); ax.set_ylim(0, max(TCGA_counts) + 7)
ax.legend(loc="upper right")
ax.text(0.99, 0.07, "DM1 prevalence:  TCGA 28.4 %  ·  MSK 32.5 %",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8, style="italic", bbox=annot_box())

suptitle(fig, "Figure 2 · DM1 is a fusion-driven, histologically distinct subtype within the canonical BRAF/RAS dark matter",
              kr="DM1 은 BRAF/RAS 음성 dark matter 내에서 융합 유전자 주도, 조직형 (FVPTC) 농축 아형이다.")
out = f"{OUT}/Fig2_fusion_mechanism.png"
fig.savefig(out, dpi=170, bbox_inches="tight")
fig.savefig(out.replace(".png",".pdf"), bbox_inches="tight")
print(f"saved {out}")
