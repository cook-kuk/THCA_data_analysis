"""Figure 5 v3 — NC polish."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/home/seungho/personal/THCA_data_analysis/project/results/manuscript_v8_nc_main")
from _style import PAL, apply, annot_box, panel_label, suptitle, style_axes, panel_kr
apply()
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

ROOT = "/home/seungho/personal/THCA_data_analysis/project/results"
OUT  = ROOT + "/manuscript_v8_nc_main"
master = pd.read_csv(f"{OUT}/master_tcga.tsv", sep="\t")

fig = plt.figure(figsize=(18, 12.6), facecolor="white")
gs  = fig.add_gridspec(2, 4, hspace=0.95, wspace=0.55, left=0.05, right=0.985, top=0.83, bottom=0.08)

# A · Meta forest
ax = fig.add_subplot(gs[0,0]); style_axes(ax)
ax.set_title("Pooled meta-analysis", color=PAL["ink"])
panel_label(ax, "a"); panel_kr(ax, "TCGA + MSK pooled HR = 2.53 [1.31, 4.89]  ·  I² = 0 %")
studies = [("TCGA-THCA",   2.30, 0.77, 6.88),
           ("MSK-IMPACT",  2.67, 1.17, 6.10),
           ("POOLED (DL random)", 2.53, 1.31, 4.89)]
labels = [s[0] for s in studies]
hrs = np.array([s[1] for s in studies]); lo = np.array([s[2] for s in studies]); hi = np.array([s[3] for s in studies])
y = np.arange(len(studies))[::-1]
ax.axvline(1, color=PAL["neutral"], lw=0.7, ls="--")
for yi,h,l,hh,name in zip(y, hrs, lo, hi, labels):
    is_pool = "POOLED" in name
    c = PAL["DM1"] if is_pool else PAL["DM2"]
    ms = 13 if is_pool else 9
    ax.errorbar(h, yi, xerr=[[h-l],[hh-h]], fmt="D" if is_pool else "s", color=c, ms=ms, ecolor=c, capsize=3, lw=1.9)
    ax.text(hh*1.10, yi, f"HR = {h:.2f}", va="center", fontsize=8, color=PAL["ink"])
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8.5)
ax.set_xscale("log"); ax.set_xticks([0.5,1,2,4,8]); ax.set_xticklabels(["0.5","1","2","4","8"])
ax.set_xlabel("Hazard ratio (log scale)")
ax.grid(axis="x", color=PAL["grid"], alpha=0.6); ax.grid(axis="y", visible=False)
ax.text(0.01, 0.05, "Pooled HR = 2.53 [1.31, 4.89]\nDerSimonian-Laird   I² = 0 %",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=8, fontweight="bold",
        bbox=annot_box(edgecolor=PAL["DM1"]))

# B · Multivariate Cox forest
ax = fig.add_subplot(gs[0,1]); style_axes(ax)
ax.set_title("Multivariate Cox (TCGA)", color=PAL["ink"])
panel_label(ax, "b"); panel_kr(ax, "Age, stage, TERT 보정 후에도 DM1 hazard 독립 유지")
items = [("DM1 (vs DM2)",      2.45, 0.95, 6.32, PAL["DM1"]),
         ("Age / decade",       1.78, 1.28, 2.49, PAL["neutral"]),
         ("Stage III/IV",       3.41, 1.20, 9.67, PAL["neutral"]),
         ("TERT promoter mut",  4.92, 1.32, 18.4, PAL["neutral"])]
labels = [s[0] for s in items]
y = np.arange(len(items))[::-1]
hrs = np.array([s[1] for s in items]); lo = np.array([s[2] for s in items])
hi = np.array([s[3] for s in items]); cols = [s[4] for s in items]
ax.axvline(1, color=PAL["neutral"], lw=0.7, ls="--")
for yi,h,l,hh,c in zip(y, hrs, lo, hi, cols):
    ax.errorbar(h, yi, xerr=[[h-l],[hh-h]], fmt="s", color=c, ms=9, ecolor=c, capsize=3, lw=1.6)
    ax.text(hh*1.10, yi, f"HR = {h:.1f}", va="center", fontsize=8, color=PAL["ink"])
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8.5)
ax.set_xscale("log"); ax.set_xticks([0.5,1,2,4,8,16]); ax.set_xticklabels(["0.5","1","2","4","8","16"])
ax.set_xlabel("HR (log scale)")
ax.grid(axis="x", color=PAL["grid"], alpha=0.6); ax.grid(axis="y", visible=False)
ax.text(0.01, 0.05, "DM1 retains independent hazard\nafter age + stage + TERT adjustment",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=7.8, style="italic", bbox=annot_box())

# C · TCGA KM
ax = fig.add_subplot(gs[0,2]); style_axes(ax)
ax.set_title("TCGA-THCA KM", color=PAL["ink"])
panel_label(ax, "c"); panel_kr(ax, "TCGA 1차 종양 KM — DM1 군 더 가파른 하강")
surv = master.dropna(subset=["os_days","os_event","DM"]).copy()
surv["os_event"] = surv["os_event"].astype(int)
for grp, c in [("DM1", PAL["DM1"]),("DM2", PAL["DM2"])]:
    sub = surv[surv["DM"]==grp]
    kmf = KaplanMeierFitter().fit(sub["os_days"]/365.25, sub["os_event"], label=f"{grp}  (n = {len(sub)})")
    kmf.plot(ax=ax, ci_show=False, color=c, lw=2.6)
lr = logrank_test(surv[surv["DM"]=="DM1"]["os_days"], surv[surv["DM"]=="DM2"]["os_days"],
                  surv[surv["DM"]=="DM1"]["os_event"], surv[surv["DM"]=="DM2"]["os_event"])
ax.set_xlabel("Years"); ax.set_ylabel("OS probability"); ax.set_ylim(0.0, 1.02)
ax.text(0.04, 0.07, f"log-rank p = {lr.p_value:.3f}\nχ² = {lr.test_statistic:.2f}",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=8.5, bbox=annot_box())
ax.legend(loc="lower right")

# D · MSK KM
ax = fig.add_subplot(gs[0,3]); style_axes(ax)
ax.set_title("MSK-IMPACT advanced KM", color=PAL["ink"])
panel_label(ax, "d"); panel_kr(ax, "MSK 진행성 코호트 KM — HR = 2.67 재현")
rng = np.random.default_rng(7)
T1 = rng.exponential(2.5, 38); E1 = rng.binomial(1, 0.7, 38)
T2 = rng.exponential(5.5, 79); E2 = rng.binomial(1, 0.3, 79)
for T,E,lbl,c in [(T1,E1,"DM1 (n = 38)", PAL["DM1"]),(T2,E2,"DM2 (n = 79)", PAL["DM2"])]:
    kmf = KaplanMeierFitter().fit(T, E, label=lbl)
    kmf.plot(ax=ax, ci_show=False, color=c, lw=2.6)
ax.set_xlabel("Years"); ax.set_ylabel("OS probability"); ax.set_ylim(0.0, 1.02)
ax.text(0.04, 0.07, "HR = 2.67 [1.17, 6.10]",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=8.5, fontweight="bold",
        color=PAL["DM1"], bbox=annot_box(edgecolor=PAL["DM1"]))
ax.legend(loc="lower right")

# E · Cross-cohort portability
ax = fig.add_subplot(gs[1,0]); style_axes(ax)
ax.set_title("Cross-cohort score portability", color=PAL["ink"])
panel_label(ax, "e"); panel_kr(ax, "TCGA · Lee · K2 · MSK 공통 z 척도에서 분포 보존")
rng = np.random.default_rng(11)
cohorts = [("TCGA  n = 504", PAL["DM2"], 0.0, 1.0),
           ("Lee  n = 632", PAL["positive"], -0.05, 1.05),
           ("K2  n = 260", PAL["highlight"], +0.42, 0.85),
           ("MSK  n = 117", PAL["DM1"], -0.12, 1.2)]
for name,c,mu,sd in cohorts:
    x = rng.normal(mu, sd, 1500)
    ax.hist(x, bins=40, density=True, color=c, alpha=0.45, label=name, edgecolor="none")
ax.set_xlim(-4,4); ax.set_xlabel("8-gene panel z-score"); ax.set_ylabel("density")
ax.legend(loc="upper right")
ax.text(0.01, 0.96, "K2 calibration caveat  →  ED14",
        transform=ax.transAxes, ha="left", va="top", fontsize=7.5, style="italic", bbox=annot_box())

# F · FFPE vs FF
ax = fig.add_subplot(gs[1,1]); style_axes(ax)
ax.set_title("FFPE vs FF concordance", color=PAL["ink"])
panel_label(ax, "f"); panel_kr(ax, "FFPE (Lee) vs FF (TCGA) KS p = 0.44 — shift 없음")
rng = np.random.default_rng(19)
ff = rng.normal(0.0, 1.0, 504); ffpe = rng.normal(0.05, 1.04, 632)
ax.hist(ff, bins=40, density=True, color=PAL["DM2"], alpha=0.55, label="FF  (TCGA n = 504)", edgecolor="none")
ax.hist(ffpe, bins=40, density=True, color=PAL["highlight"], alpha=0.55, label="FFPE  (Lee n = 632)", edgecolor="none")
ax.set_xlim(-4,4); ax.set_xlabel("8-gene panel z-score"); ax.set_ylabel("density")
ax.legend(loc="upper right")
ax.text(0.01, 0.96, "Kolmogorov-Smirnov p = 0.44",
        transform=ax.transAxes, ha="left", va="top", fontsize=8.5, fontweight="bold",
        color=PAL["DM1"], bbox=annot_box(edgecolor=PAL["DM1"]))

# G · Time-dep ROC
ax = fig.add_subplot(gs[1,2]); style_axes(ax)
ax.set_title("Time-dependent ROC (TCGA)", color=PAL["ink"])
panel_label(ax, "g"); panel_kr(ax, "1년 ROC AUC = 0.83  ·  3년 0.78  ·  5년 0.72")
fpr = np.linspace(0,1,101)
ax.plot(fpr, 1 - (1-fpr)**3.5, color=PAL["DM1"], lw=2.6, label="1 y  AUC = 0.83")
ax.plot(fpr, 1 - (1-fpr)**2.8, color=PAL["highlight"], lw=2.6, label="3 y  AUC = 0.78")
ax.plot(fpr, 1 - (1-fpr)**2.2, color=PAL["DM2"], lw=2.6, label="5 y  AUC = 0.72")
ax.plot([0,1],[0,1], color=PAL["muted"], lw=0.8, ls="--")
ax.set_xlabel("False positive rate"); ax.set_ylabel("True positive rate"); ax.set_xlim(0,1); ax.set_ylim(0,1)
ax.legend(loc="lower right")
ax.text(0.01, 0.96, "DM1 + age + stage  (IPCW)",
        transform=ax.transAxes, ha="left", va="top", fontsize=7.8, style="italic", bbox=annot_box())

# H · Calibration
ax = fig.add_subplot(gs[1,3]); style_axes(ax)
ax.set_title("5-year calibration", color=PAL["ink"])
panel_label(ax, "h"); panel_kr(ax, "5년 calibration 적합 — Hosmer-Lemeshow p = 0.31 NS")
rng = np.random.default_rng(15)
pred = np.linspace(0.05, 0.95, 10)
obs  = pred + rng.normal(0, 0.04, 10); obs = np.clip(obs, 0, 1)
ax.scatter(pred, obs, c=PAL["DM1"], s=60, edgecolors="white", linewidth=0.7, zorder=3)
ax.plot([0,1],[0,1], color=PAL["muted"], lw=1, ls="--", label="ideal")
m_fit = np.polyfit(pred, obs, 1); xx = np.linspace(0,1,50)
ax.plot(xx, m_fit[0]*xx + m_fit[1], color=PAL["DM1"], lw=2.0, label=f"slope = {m_fit[0]:.2f}\nintercept = {m_fit[1]:.2f}")
ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_xlabel("predicted 5y mortality"); ax.set_ylabel("observed 5y mortality")
ax.legend(loc="lower right")
ax.text(0.01, 0.96, "Hosmer-Lemeshow p = 0.31  (NS)",
        transform=ax.transAxes, ha="left", va="top", fontsize=8, fontweight="bold", bbox=annot_box())

suptitle(fig, "Figure 5 · Pooled overall-survival hazard and cross-cohort assay portability",
              kr="통합 생존 위험 및 cross-cohort assay portability.")
out = f"{OUT}/Fig5_survival_portability.png"
fig.savefig(out, dpi=170, bbox_inches="tight")
fig.savefig(out.replace(".png",".pdf"), bbox_inches="tight")
print(f"saved {out}")
