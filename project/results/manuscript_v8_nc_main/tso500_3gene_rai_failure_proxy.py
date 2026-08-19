"""
TSO500-3 (TSHR · PAX8 · NKX2-1) — RAI failure proxy KM test.

임상적으로 RAI failure 가 우려되는 high-risk stratum 에서만 검정 :
  1. BRAF V600E positive  (RAI failure 의 가장 흔한 분자 origin)
  2. Stage III/IV  (RAI 가 더 적극 적용되는 군)
  3. Age > 55     (예후 risk factor)
  4. TERT positive (가장 aggressive)
  5. Joint: BRAF+ AND (stage III/IV OR TERT+)  =  "RAI failure risk pool"

각 stratum 안에서 PFI 가 가장 events 많고 임상 의미 가까움.
3-gene mean score → tertile (top-33% vs bottom-33%) split → log-rank + KM curves.
"""
import os, json, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

for f in ["/usr/share/fonts/truetype/nanum/NanumGothic.ttf"]:
    if os.path.exists(f):
        font_manager.fontManager.addfont(f); break
plt.rcParams["font.family"] = ["NanumGothic","DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

ROOT = "/home/seungho/personal/THCA_data_analysis/project"
OUT  = f"{ROOT}/dm1_story_web/public/figures"

# Load β + clinical
hm = pd.read_csv(f"{ROOT}/results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv", sep="\t")
cl = pd.read_csv(f"{ROOT}/results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv", sep="\t")
GENES_T3 = ["TSHR","PAX8","NKX2-1"]
GENES_F8 = ["DIO1","FOXE1","NKX2-1","PAX8","SLC5A5","TG","TPO","TSHR"]
df = hm.merge(cl, left_on="sample_short", right_on="tcga_short", how="inner").dropna(subset=GENES_F8)
# Score: high β = DM1-like = worse
df["t3"] = df[GENES_T3].mean(axis=1)
df["f8"] = df[GENES_F8].mean(axis=1)
print(f"Merged n = {len(df)}")

df["age"] = pd.to_numeric(df.get("age"), errors="coerce")
df["stage_hi"]  = df["stage"].astype(str).str.contains("III|IV", regex=True, na=False)
df["braf_pos"]  = df["has_braf_v600e"].astype(str).str.lower().isin(["1","true","yes"])
df["tert_pos"]  = df["tert_pos"].astype(str).str.lower().isin(["1","true","yes"])
df["age_55"]    = df["age"] > 55
df["rai_risk_pool"] = df["braf_pos"] & (df["stage_hi"] | df["tert_pos"])

STRATA = {
    "Full cohort":          np.ones(len(df), bool),
    "BRAF V600E+":          df["braf_pos"].values,
    "Stage III/IV":         df["stage_hi"].values,
    "TERT-positive":        df["tert_pos"].values,
    "Age > 55":             df["age_55"].fillna(False).values,
    "BRAF+ & (StgIII/IV or TERT+)  ★ RAI-risk pool": df["rai_risk_pool"].values
}

def kmcox_test(sub, score_col="t3", tcol="PFI.time", ecol="PFI", split="tertile"):
    m = sub[[score_col, tcol, ecol]].dropna()
    m = m[m[tcol] > 0]
    n = len(m); ev = int(m[ecol].sum())
    if n < 12 or ev < 4: return dict(n=n, ev=ev)
    # Tertile split
    if split == "tertile":
        q1, q2 = m[score_col].quantile([1/3, 2/3])
        hi = m[m[score_col] >= q2]; lo = m[m[score_col] <= q1]
    elif split == "median":
        med = m[score_col].median()
        hi = m[m[score_col] > med]; lo = m[m[score_col] <= med]
    if len(hi) < 5 or len(lo) < 5 or hi[ecol].sum() < 2 or lo[ecol].sum() < 2:
        return dict(n=n, ev=ev, n_hi=len(hi), n_lo=len(lo), ev_hi=int(hi[ecol].sum()), ev_lo=int(lo[ecol].sum()))
    lr = logrank_test(hi[tcol], lo[tcol], hi[ecol], lo[ecol])
    # Cox HR on z-score
    mz = m.copy()
    mz["z"] = (mz[score_col]-mz[score_col].mean())/mz[score_col].std()
    try:
        cph = CoxPHFitter().fit(mz[[tcol,ecol,"z"]], duration_col=tcol, event_col=ecol)
        hr = float(np.exp(cph.params_["z"]))
        cox_p = float(cph.summary.loc["z","p"])
        ci_lo = float(np.exp(cph.confidence_intervals_.loc["z"].iloc[0]))
        ci_hi = float(np.exp(cph.confidence_intervals_.loc["z"].iloc[1]))
    except Exception:
        hr=cox_p=ci_lo=ci_hi=np.nan
    return dict(n=n, ev=ev, n_hi=len(hi), n_lo=len(lo),
                ev_hi=int(hi[ecol].sum()), ev_lo=int(lo[ecol].sum()),
                logrank_p=float(lr.p_value), hr=hr, hr_lo=ci_lo, hr_hi=ci_hi, cox_p=cox_p,
                hi_subset=hi, lo_subset=lo)

# Run for all strata × T3 vs F8 × PFI/OS
rows = []
for label, mask in STRATA.items():
    sub = df[mask]
    for ep, (tcol, ecol) in [("PFI", ("PFI.time","PFI")), ("OS", ("OS.time","OS"))]:
        for panel, col in [("TSO500-3", "t3"), ("Full-8", "f8")]:
            res = kmcox_test(sub, col, tcol, ecol, split="tertile")
            rows.append({"stratum": label, "endpoint": ep, "panel": panel, **{k:v for k,v in res.items() if k not in ["hi_subset","lo_subset"]}})

df_out = pd.DataFrame(rows)
df_out.to_csv(f"{OUT}/tso500_rai_failure_proxy.tsv", sep="\t", index=False)
with open(f"{OUT}/tso500_rai_failure_proxy.json","w") as f:
    json.dump(rows, f, indent=2, default=str)

print(f"\n{'Stratum':<55s}  {'EP':<4s}  {'Panel':<10s}  {'n':>4s}/{'ev':<3s}  {'hi/lo':>9s}  {'HR':>6s}[{'95CI':<11s}] {'CoxP':>9s}  {'LR-p':>9s}  sig")
for r in rows:
    if r.get("logrank_p") is None: continue
    if "n" not in r: continue
    sig = "★★★" if r["logrank_p"]<0.001 else "★★" if r["logrank_p"]<0.01 else "★" if r["logrank_p"]<0.05 else ""
    print(f"  {r['stratum']:<53s}  {r['endpoint']:<3s}  {r['panel']:<9s}  {r['n']:4d}/{r['ev']:<3d}  {r['n_hi']:3d}/{r['n_lo']:<3d}  {r['hr']:6.2f}[{r['hr_lo']:.2f}-{r['hr_hi']:.2f}] {r['cox_p']:9.3e}  {r['logrank_p']:9.3e} {sig}")

# ─── KM grid: rows = strata, cols = panels (PFI only since OS too few events) ───
strata_show = ["Full cohort","BRAF V600E+","Stage III/IV","TERT-positive","Age > 55","BRAF+ & (StgIII/IV or TERT+)  ★ RAI-risk pool"]
fig, axes = plt.subplots(len(strata_show), 2, figsize=(13.5, 3.0*len(strata_show)), dpi=150)
for i, lab in enumerate(strata_show):
    mask = STRATA[lab]
    sub = df[mask]
    for j, (panel, col, color) in enumerate([("TSO500-3 (TSHR · PAX8 · NKX2-1)", "t3", "#B91C1C"),
                                              ("Full-8 (reference)",              "f8", "#0F172A")]):
        ax = axes[i,j]
        res = kmcox_test(sub, col, "PFI.time", "PFI", split="tertile")
        if res.get("logrank_p") is None or "hi_subset" not in res:
            ax.text(0.5, 0.5, f"n = {res.get('n','?')} · events = {res.get('ev','?')}\n불충분", transform=ax.transAxes, ha="center", va="center", fontsize=11)
            ax.set_axis_off()
            continue
        hi = res["hi_subset"]; lo = res["lo_subset"]
        kmf_hi = KaplanMeierFitter().fit(hi["PFI.time"]/365.25, hi["PFI"], label=f"DM1-like (top T3, n={len(hi)}, ev={int(hi['PFI'].sum())})")
        kmf_lo = KaplanMeierFitter().fit(lo["PFI.time"]/365.25, lo["PFI"], label=f"non-DM1 (bot T3, n={len(lo)}, ev={int(lo['PFI'].sum())})")
        kmf_hi.plot_survival_function(ax=ax, color=color, ci_alpha=0.13, lw=2.4)
        kmf_lo.plot_survival_function(ax=ax, color="#94A3B8", ci_alpha=0.13, lw=2.4)
        sig = "★★★" if res["logrank_p"]<0.001 else "★★" if res["logrank_p"]<0.01 else "★" if res["logrank_p"]<0.05 else "ns"
        ax.text(0.04, 0.10,
                f"Log-rank p = {res['logrank_p']:.3e}  {sig}\n"
                f"Cox HR     = {res['hr']:.2f} [{res['hr_lo']:.2f}-{res['hr_hi']:.2f}]\n"
                f"Cox p      = {res['cox_p']:.2e}",
                transform=ax.transAxes, fontsize=9, family="monospace",
                bbox=dict(facecolor="white", edgecolor=color if res["logrank_p"]<0.05 else "#CBD5E1", linewidth=1.2, boxstyle="round,pad=0.4"))
        ax.set_xlabel("Years (PFI)"); ax.set_ylabel("PFI probability" if j==0 else "")
        ax.set_ylim(0.55, 1.02); ax.set_xlim(0, 15)
        ax.set_title(f"{lab}  ·  {panel}", fontsize=10.5, fontweight="bold", color=color, loc="left")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.legend(loc="lower left", fontsize=8, frameon=False)
fig.suptitle("TSO500-3 vs Full-8 — PFI KM (tertile split) across high-risk strata",
             fontsize=14, fontweight="bold", y=1.005)
plt.tight_layout()
fig.savefig(f"{OUT}/fig_tso500_rai_failure_km_strata.png", bbox_inches="tight")
plt.close(fig)

print(f"\nWrote KM grid:  {OUT}/fig_tso500_rai_failure_km_strata.png")
