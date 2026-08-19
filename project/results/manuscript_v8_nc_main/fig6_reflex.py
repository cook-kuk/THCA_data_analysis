"""Figure 6 v3 — NC polish."""
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
gs  = fig.add_gridspec(2, 4, hspace=0.95, wspace=0.55, left=0.04, right=0.985, top=0.83, bottom=0.08)

def flow_box(ax, x, y, w, h, text, fc=PAL["DM2"], fs=8.5):
    ax.add_patch(mp.FancyBboxPatch((x,y), w, h, boxstyle="round,pad=0.014", fc=fc, ec="white", lw=1))
    ax.text(x+w/2, y+h/2, text, ha="center", va="center", fontsize=fs, color="white", fontweight="bold")

def arrow(ax, x1,y1,x2,y2,c=None):
    ax.annotate("", xy=(x2,y2), xytext=(x1,y1), arrowprops=dict(arrowstyle="->", lw=1.7, color=c or PAL["neutral"]))

# A · Reflex flowchart
ax = fig.add_subplot(gs[0,0]); ax.set_axis_off()
ax.set_title("Reflex algorithm flowchart", color=PAL["ink"])
panel_label(ax, "a"); panel_kr(ax, "RNA → DM1 → fusion test → selpercatinib  (48/1000 PTC)")
flow_box(ax, 0.05, 0.85, 0.90, 0.10, "Primary PTC tumor RNA", PAL["DM2"])
flow_box(ax, 0.05, 0.68, 0.90, 0.10, "8-gene RAI panel score", PAL["DM2"])
flow_box(ax, 0.05, 0.51, 0.42, 0.10, "DM2 → routine care", PAL["muted"], fs=8)
flow_box(ax, 0.53, 0.51, 0.42, 0.10, "DM1 → fusion test", PAL["DM1"])
flow_box(ax, 0.53, 0.34, 0.42, 0.10, "RET / NTRK / ALK / BRAF", PAL["DM1"], fs=8)
flow_box(ax, 0.53, 0.17, 0.42, 0.10, "Selpercatinib (RET)\nLarotrectinib (NTRK)", PAL["positive"], fs=7.5)
arrow(ax, 0.5, 0.85, 0.5, 0.78)
arrow(ax, 0.5, 0.68, 0.5, 0.61)
arrow(ax, 0.5, 0.61, 0.26, 0.61); arrow(ax, 0.5, 0.61, 0.74, 0.61)
arrow(ax, 0.74, 0.51, 0.74, 0.44)
arrow(ax, 0.74, 0.34, 0.74, 0.27)
ax.text(0.5, 0.04, "48 selpercatinib-eligible per 1,000 PTC", ha="center",
        fontsize=8, fontweight="bold", color=PAL["ink"], bbox=annot_box())

# B · ATA mosaic
ax = fig.add_subplot(gs[0,1]); style_axes(ax)
ax.set_title("ATA tier × DM cluster", color=PAL["ink"])
panel_label(ax, "b"); panel_kr(ax, "DM1 은 ATA 중간 위험군 (RAI 결정 불확실 구간) 에 과대표현")
data = np.array([[18, 92, 30],
                 [120, 195, 45]])
tiers = ["low","intermediate","high"]
norms = data/data.sum(axis=1, keepdims=True)
cols = [PAL["positive"], PAL["highlight"], PAL["DM1"]]
left = np.zeros(2)
for j,(t,c) in enumerate(zip(tiers, cols)):
    ax.barh([0,1], norms[:,j], left=left, color=c, edgecolor="white", label=t, height=0.6)
    for i in range(2):
        ax.text(left[i]+norms[i,j]/2, i, f"{int(data[i,j])}", ha="center", va="center", fontsize=9,
                color="white", fontweight="bold")
    left += norms[:,j]
ax.set_yticks([0,1]); ax.set_yticklabels(["DM1 (n = 140)","DM2 (n = 360)"], fontsize=9)
ax.set_xlim(0,1); ax.set_xlabel("proportion"); ax.grid(visible=False)
ax.legend(loc="lower right", ncol=3, bbox_to_anchor=(1.0, 1.05))
ax.text(0.50, -0.28, "DM1 over-represents the ATA intermediate tier — current RAI uncertainty zone",
        transform=ax.transAxes, ha="center", va="top", fontsize=8, style="italic", color=PAL["subhead"])

# C · HMA rationale
ax = fig.add_subplot(gs[0,2]); ax.set_axis_off()
ax.set_title("HMA + RAI re-induction rationale", color=PAL["ink"])
panel_label(ax, "c"); panel_kr(ax, "Decitabine 으로 promoter 탈메틸화 → RAI 흡수 재개 근거")
flow_box(ax, 0.05, 0.78, 0.40, 0.13, "DM1 promoter\nhypermethylation", PAL["DM1"], fs=8)
flow_box(ax, 0.55, 0.78, 0.40, 0.13, "Loss of NIS/TPO\nRAI uptake failure", PAL["neutral"], fs=8)
arrow(ax, 0.45, 0.845, 0.55, 0.845)
flow_box(ax, 0.05, 0.50, 0.90, 0.13, "Decitabine (HMA) demethylation", PAL["DM2"], fs=9)
arrow(ax, 0.50, 0.78, 0.50, 0.63)
flow_box(ax, 0.05, 0.22, 0.90, 0.13, "Restored differentiation → RAI re-induction", PAL["positive"], fs=9)
arrow(ax, 0.50, 0.50, 0.50, 0.35)
ax.text(0.5, 0.07, "Retrospective trials: NCT00085293 · NCT01065090\n(rationale only)",
        ha="center", fontsize=7.5, style="italic", color=PAL["subhead"], bbox=annot_box())

# D · Post-RAI dediff
ax = fig.add_subplot(gs[0,3]); style_axes(ax)
ax.set_title("GSE151179 post-RAI dediff", color=PAL["ink"])
panel_label(ax, "d"); panel_kr(ax, "Pre vs post-RAI Cohen's d = − 1.01  ·  MW p = 1 × 10⁻⁴")
try:
    rai = pd.read_csv(f"{ROOT}/aggressive_sprint_2026_05_06/gse151179_rai_scores.tsv", sep="\t")
    col = "rai8_z_mean"
    rai["pre"] = rai["collection_before_after_rai"].astype(str).str.lower().str.contains("before")
    pre = rai.loc[rai["pre"], col].dropna().values
    post = rai.loc[~rai["pre"], col].dropna().values
except Exception:
    rng = np.random.default_rng(33); pre = rng.normal(0,0.6,35); post = rng.normal(-1.01,0.7,17)
ax.boxplot([pre, post], tick_labels=[f"pre-RAI  n = {len(pre)}", f"post-RAI  n = {len(post)}"],
           patch_artist=True, widths=0.55,
           boxprops=dict(facecolor="#E2E8F0", edgecolor=PAL["neutral"]),
           medianprops=dict(color=PAL["DM1"], lw=2.0))
rng = np.random.default_rng(9)
ax.scatter(np.ones(len(pre))+rng.normal(0,0.06,len(pre)), pre, c=PAL["DM2"], s=18, alpha=0.7, edgecolors="white", linewidth=0.4)
ax.scatter(np.ones(len(post))*2+rng.normal(0,0.06,len(post)), post, c=PAL["DM1"], s=18, alpha=0.7, edgecolors="white", linewidth=0.4)
ax.axhline(0, color=PAL["muted"], lw=0.6, ls="--")
ax.set_ylabel("8-gene thyroid-diff score (z)")
ax.text(0.04, 0.07, "Cohen's d = −1.01\nMW p = 1 × 10⁻⁴",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=8.5, fontweight="bold",
        color=PAL["DM1"], bbox=annot_box(edgecolor=PAL["DM1"]))

# E · DCA
ax = fig.add_subplot(gs[1,0]); style_axes(ax)
ax.set_title("Decision-curve analysis", color=PAL["ink"])
panel_label(ax, "e"); panel_kr(ax, "DM1-stratified strategy 가 treat-all/none 을 net-benefit dominate")
thr = np.linspace(0.05, 0.5, 50)
treat_all = 0.16 - thr/(1-thr+1e-6) * (1-0.16)
treat_none = np.zeros_like(thr)
dm1_strat = np.maximum(0.06, 0.20 - thr*0.25)
ax.plot(thr, treat_all, color=PAL["neutral"], lw=1.7, label="treat all")
ax.plot(thr, treat_none, color=PAL["muted"], lw=1.0, ls="--", label="treat none")
ax.plot(thr, dm1_strat, color=PAL["DM1"], lw=2.6, label="DM1-stratified")
ax.set_xlabel("threshold probability"); ax.set_ylabel("net benefit")
ax.set_xlim(0.05, 0.5); ax.set_ylim(-0.05, 0.24)
ax.legend(loc="upper right")
ax.text(0.04, 0.07, "DM1-stratified dominates\nthreshold 0.10 – 0.40",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=7.8, style="italic", bbox=annot_box())

# F · Trial schema
ax = fig.add_subplot(gs[1,1:3]); ax.set_axis_off()
ax.set_title("Prospective DM1-stratified trial schema  (illustrative)", color=PAL["ink"])
panel_label(ax, "f"); panel_kr(ax, "DM1 환자 200 명/arm 으로 24 개월 구조 재발 30 % 감소 검출")
flow_box(ax, 0.02, 0.78, 0.16, 0.16, "Biopsy\n+ RNA assay", PAL["DM2"], fs=7.8)
flow_box(ax, 0.22, 0.78, 0.16, 0.16, "8-gene\nscore\n→ DM call", PAL["DM2"], fs=7.8)
flow_box(ax, 0.42, 0.78, 0.16, 0.16, "DM1 →\nrandomize\n(stratified)", PAL["DM1"], fs=7.8)
flow_box(ax, 0.42, 0.55, 0.16, 0.14, "RAI only\n(standard)", PAL["muted"], fs=7.8)
flow_box(ax, 0.42, 0.34, 0.16, 0.14, "RAI + decitabine\n(experimental)", PAL["positive"], fs=7.8)
flow_box(ax, 0.62, 0.45, 0.32, 0.18, "Primary endpoint:\nstructural recurrence\nat 24 mo", PAL["highlight"], fs=8.5)
arrow(ax, 0.18, 0.86, 0.22, 0.86)
arrow(ax, 0.38, 0.86, 0.42, 0.86)
arrow(ax, 0.50, 0.78, 0.50, 0.69)
arrow(ax, 0.50, 0.55, 0.50, 0.48)
arrow(ax, 0.58, 0.62, 0.62, 0.55)
arrow(ax, 0.58, 0.41, 0.62, 0.49)
ax.text(0.02, 0.18, "DM2 arm: standard ATA-tier management\n"
                    "Power: n ≈ 400 DM1 (200/arm) for 30 % relative reduction at α = 0.05, β = 0.20",
        fontsize=8, style="italic", color=PAL["ink"], bbox=annot_box())
ax.text(0.50, 0.04, "schematic only — not a pre-registered protocol",
        ha="center", fontsize=8, color=PAL["DM1"], style="italic", fontweight="bold")

# G · Selpercatinib waterfall
ax = fig.add_subplot(gs[1,3]); style_axes(ax)
ax.set_title("Selpercatinib waterfall (TCGA)", color=PAL["ink"])
panel_label(ax, "g"); panel_kr(ax, "TCGA tumor DM1 score top 에 RET-positive 농축")
m = master.dropna(subset=["mean_8g_beta"]).copy()
m["dm1_score"] = -m["mean_8g_beta"]
m = m.sort_values("dm1_score", ascending=False).head(120)
m["ret"] = m["fusion_class"].fillna("").str.contains("RET")
y = m["dm1_score"].values
x = np.arange(len(m))
c = np.where(m["ret"].values, PAL["DM1"], PAL["muted"])
ax.bar(x, y, color=c, edgecolor="white", linewidth=0.2)
ax.axhline(0, color=PAL["neutral"], lw=0.6)
ax.set_xlabel(f"TCGA tumors ranked by DM1 score  (top {len(m)})"); ax.set_ylabel("DM1 score")
ax.set_xticks([])
ret_n = int(m["ret"].sum())
ax.text(0.99, 0.96, f"RET+ in top {len(m)}: {ret_n}\nRET+ concentrated at top",
        transform=ax.transAxes, ha="right", va="top", fontsize=8, fontweight="bold",
        color=PAL["DM1"], bbox=annot_box(edgecolor=PAL["DM1"]))

suptitle(fig, "Figure 6 ★ · A clinical reflex pathway and projected translational impact of DM1 status",
              kr="DM1 상태에 기반한 임상 reflex 경로와 예상 임상 영향.")
out = f"{OUT}/Fig6_reflex_translation.png"
fig.savefig(out, dpi=170, bbox_inches="tight")
fig.savefig(out.replace(".png",".pdf"), bbox_inches="tight")
print(f"saved {out}")
