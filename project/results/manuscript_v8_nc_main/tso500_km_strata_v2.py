"""
Regenerate Fig SR-1 with cleaner 3x2 layout · much larger per-panel · no overlapping lines.
"""
import os, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

for f in ["/usr/share/fonts/truetype/nanum/NanumGothic.ttf"]:
    if os.path.exists(f): font_manager.fontManager.addfont(f); break
plt.rcParams["font.family"] = ["NanumGothic","DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

ROOT = "/home/seungho/personal/THCA_data_analysis/project"
OUT  = f"{ROOT}/dm1_story_web/public/figures"

hm = pd.read_csv(f"{ROOT}/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv", sep="\t")
cl = pd.read_csv(f"{ROOT}/results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv", sep="\t")
GENES_F8 = ["DIO1","FOXE1","NKX2-1","PAX8","SLC5A5","TG","TPO","TSHR"]
df = hm.merge(cl, left_on="sample_short", right_on="tcga_short", how="inner").dropna(subset=GENES_F8)
df["t3"] = df[["TSHR","PAX8","NKX2-1"]].mean(axis=1)
df["f8"] = df[GENES_F8].mean(axis=1)
df["ihc3"] = df[["TG","PAX8","NKX2-1"]].mean(axis=1)
df["age"] = pd.to_numeric(df.get("age"), errors="coerce")
df["stage_hi"] = df["stage"].astype(str).str.contains("III|IV", regex=True, na=False)
df["braf_pos"] = df["has_braf_v600e"].astype(str).str.lower().isin(["1","true","yes"])
df["tert_pos"] = df["tert_pos"].astype(str).str.lower().isin(["1","true","yes"])
df["age_55"]   = df["age"] > 55

STRATA = {
    "Full cohort":               np.ones(len(df), bool),
    "BRAF V600E+":               df["braf_pos"].values,
    "Stage III/IV":              df["stage_hi"].values,
    "TERT+  (small n)":          df["tert_pos"].values,
    "Age > 55":                  df["age_55"].fillna(False).values,
    "BRAF+ & (StgIII/IV or TERT+)": (df["braf_pos"] & (df["stage_hi"] | df["tert_pos"])).values
}
strata_show = list(STRATA.keys())

def km_stats(sub, col):
    m = sub[[col,"PFI.time","PFI"]].dropna()
    m = m[m["PFI.time"] > 0]
    if len(m) < 12 or m["PFI"].sum() < 4: return None
    q1, q2 = m[col].quantile([1/3, 2/3])
    hi = m[m[col] >= q2]; lo = m[m[col] <= q1]
    if len(hi) < 5 or len(lo) < 5: return None
    lr = logrank_test(hi["PFI.time"], lo["PFI.time"], hi["PFI"], lo["PFI"])
    z = (m[col]-m[col].mean())/m[col].std()
    try:
        cph = CoxPHFitter().fit(pd.DataFrame({"t":m["PFI.time"], "e":m["PFI"], "z":z}),
                                duration_col="t", event_col="e")
        hr    = float(np.exp(cph.params_["z"]))
        cox_p = float(cph.summary.loc["z","p"])
        ci_lo = float(np.exp(cph.confidence_intervals_.loc["z"].iloc[0]))
        ci_hi = float(np.exp(cph.confidence_intervals_.loc["z"].iloc[1]))
    except Exception:
        hr = cox_p = ci_lo = ci_hi = np.nan
    return {"hi":hi, "lo":lo, "lr_p":lr.p_value, "hr":hr, "ci_lo":ci_lo, "ci_hi":ci_hi, "cox_p":cox_p, "n":len(m), "ev":int(m["PFI"].sum())}

PANELS = [("TSO500-3", "t3", "#B91C1C"), ("Full-8", "f8", "#0F172A")]

# ─── Grid layout: 6 strata x 2 panels = 12 subplots (6 rows × 2 cols)
fig = plt.figure(figsize=(22, 32), dpi=140)
gs = fig.add_gridspec(6, 2, hspace=0.55, wspace=0.22, top=0.965, bottom=0.03, left=0.06, right=0.98)

for i, lab in enumerate(strata_show):
    sub = df[STRATA[lab]]
    for j, (pname, col, color) in enumerate(PANELS):
        ax = fig.add_subplot(gs[i, j])
        res = km_stats(sub, col)
        if res is None:
            ax.text(0.5, 0.5, "n / events 부족", transform=ax.transAxes, ha="center", va="center", fontsize=14)
            ax.set_axis_off()
            continue
        hi, lo = res["hi"], res["lo"]
        kmf_hi = KaplanMeierFitter().fit(hi["PFI.time"]/365.25, hi["PFI"], label=f"DM1-like top-33 %  n={len(hi)}  ev={int(hi['PFI'].sum())}")
        kmf_lo = KaplanMeierFitter().fit(lo["PFI.time"]/365.25, lo["PFI"], label=f"non-DM1 bot-33 %  n={len(lo)}  ev={int(lo['PFI'].sum())}")
        kmf_hi.plot_survival_function(ax=ax, color=color, ci_alpha=0.14, lw=3.5)
        kmf_lo.plot_survival_function(ax=ax, color="#94A3B8", ci_alpha=0.14, lw=3.5)
        sig = "★★★" if res["lr_p"]<0.001 else "★★" if res["lr_p"]<0.01 else "★" if res["lr_p"]<0.05 else "n.s."
        col_frame = color if res["lr_p"]<0.05 else "#CBD5E1"
        # Stats box
        stat_text = (f"Log-rank p = {res['lr_p']:.3e}   {sig}\n"
                     f"Cox HR = {res['hr']:.2f}  [{res['ci_lo']:.2f}-{res['ci_hi']:.2f}]\n"
                     f"Cox p = {res['cox_p']:.2e}")
        ax.text(0.05, 0.12, stat_text, transform=ax.transAxes, fontsize=13,
                family="monospace", va="bottom",
                bbox=dict(facecolor="white", edgecolor=col_frame, linewidth=2.5, boxstyle="round,pad=0.7"))
        ax.set_xlabel("Years  (progression-free)", fontsize=13)
        ax.set_ylabel("PFI probability" if j == 0 else "", fontsize=13)
        ax.set_ylim(0.55, 1.02); ax.set_xlim(0, 15)
        ax.grid(True, ls=":", alpha=0.35)
        # Title
        head = f"{lab}  ·  {pname}"
        ax.set_title(head, fontsize=16, fontweight="bold", color=color, loc="left", pad=8)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.legend(loc="lower left", fontsize=11, frameon=True, edgecolor="#CBD5E1")

fig.suptitle("Fig SR-1 · PFI Kaplan-Meier — TSO500-3 (TSHR·PAX8·NKX2-1) vs Full-8, 6 개 high-risk stratum",
             fontsize=20, fontweight="bold", y=0.988)
fig.savefig(f"{OUT}/fig_tso500_rai_failure_km_strata.png", bbox_inches="tight")
plt.close(fig)
print(f"Wrote {OUT}/fig_tso500_rai_failure_km_strata.png")
