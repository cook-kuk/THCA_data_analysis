#!/usr/bin/env python3
"""Sprint O — Multivariate Cox + dichotomized survival on TCGA-THCA.
Tests if DM1_like effect survives after adjusting for BRAF, RAS, age, stage."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

OUT = Path("project_external_st/results/extra")

df = pd.read_csv(OUT / "s_tcga_thca_scored.tsv", sep="\t")
print(f"loaded {len(df)} TCGA-THCA samples")

# encode covariates
stage_map = {"Stage I":1,"Stage II":2,"Stage III":3,"Stage IVA":4,"Stage IVB":5,"Stage IVC":6}
df["stage_ord"] = df["ajcc_pathologic_tumor_stage"].map(stage_map)
df["age_z"] = (df["age"] - df["age"].mean()) / df["age"].std()
df["braf"] = df["has_braf_v600e"].astype(float)
df["ras"]  = df["has_ras_mut"].astype(float)
df["dm1_z"] = (df["DM1_like"] - df["DM1_like"].mean()) / df["DM1_like"].std()

# === 1. Univariate (already known but for table) ===
def cox_univariate(df, covar, time_col, event_col):
    sub = df[[time_col, event_col, covar]].dropna()
    if len(sub) < 30: return None
    cph = CoxPHFitter().fit(sub.rename(columns={time_col:"time", event_col:"event"}),
                             duration_col="time", event_col="event")
    return {"covariate": covar, "n": len(sub),
            "HR": cph.hazard_ratios_[covar],
            "HR_ci_lo": np.exp(cph.confidence_intervals_.loc[covar, "95% lower-bound"]),
            "HR_ci_hi": np.exp(cph.confidence_intervals_.loc[covar, "95% upper-bound"]),
            "p": cph.summary.loc[covar, "p"]}

uni_rows = []
for outcome in [("PFI.time","PFI"),("DFI.time","DFI"),("OS.time","OS"),("DSS.time","DSS")]:
    for cov in ["dm1_z","braf","ras","age_z","stage_ord"]:
        r = cox_univariate(df, cov, outcome[0], outcome[1])
        if r:
            r["outcome"] = outcome[1]; r["model"] = "univariate"; uni_rows.append(r)

# === 2. Multivariate ===
multi_rows = []
for outcome in [("PFI.time","PFI"),("DFI.time","DFI"),("OS.time","OS"),("DSS.time","DSS")]:
    sub = df[[outcome[0], outcome[1], "dm1_z","braf","ras","age_z","stage_ord"]].dropna()
    if len(sub) < 30: continue
    sub = sub.rename(columns={outcome[0]:"time", outcome[1]:"event"})
    try:
        cph = CoxPHFitter(penalizer=0.01).fit(sub, duration_col="time", event_col="event")
        for cov in ["dm1_z","braf","ras","age_z","stage_ord"]:
            multi_rows.append({"outcome": outcome[1], "model": "multivariate",
                               "covariate": cov, "n": len(sub),
                               "HR": cph.hazard_ratios_[cov],
                               "HR_ci_lo": np.exp(cph.confidence_intervals_.loc[cov, "95% lower-bound"]),
                               "HR_ci_hi": np.exp(cph.confidence_intervals_.loc[cov, "95% upper-bound"]),
                               "p": cph.summary.loc[cov, "p"]})
    except Exception as e:
        print(f"  [skip] {outcome[1]} multivariate: {e}")

cox_df = pd.DataFrame(uni_rows + multi_rows)
cox_df.to_csv(OUT / "o_multivariate_cox.tsv", sep="\t", index=False)
print("\n=== Multivariate Cox (DM1_like covariate p across outcomes) ===")
m = cox_df[(cox_df["model"]=="multivariate") & (cox_df["covariate"]=="dm1_z")]
print(m[["outcome","n","HR","HR_ci_lo","HR_ci_hi","p"]].to_string(index=False))

# === 3. Dichotomized: top 30% vs bottom 30% PFI ===
ok = df.dropna(subset=["PFI","PFI.time","DM1_like"]).copy()
q30, q70 = np.percentile(ok["DM1_like"], [30, 70])
ok["DM1_dichotomy"] = pd.cut(ok["DM1_like"], bins=[-np.inf, q30, q70, np.inf],
                              labels=["DM1_low (bot30)","mid","DM1_high (top30)"])
lo = ok[ok["DM1_dichotomy"]=="DM1_low (bot30)"]
hi = ok[ok["DM1_dichotomy"]=="DM1_high (top30)"]
lr = logrank_test(lo["PFI.time"], hi["PFI.time"], lo["PFI"], hi["PFI"])
print(f"\n=== Dichotomized PFI (DM1 top30 vs bot30) ===")
print(f"  log-rank p = {lr.p_value:.2e}")

# Cox HR for dichotomized
ok["dm1_high"] = (ok["DM1_dichotomy"] == "DM1_high (top30)").astype(int)
ok2 = ok[ok["DM1_dichotomy"].isin(["DM1_low (bot30)","DM1_high (top30)"])]
ok2 = ok2[["PFI.time","PFI","dm1_high"]].rename(columns={"PFI.time":"time","PFI":"event"})
cph_d = CoxPHFitter().fit(ok2, duration_col="time", event_col="event")
hr_d = cph_d.hazard_ratios_["dm1_high"]
ci_lo = np.exp(cph_d.confidence_intervals_.loc["dm1_high","95% lower-bound"])
ci_hi = np.exp(cph_d.confidence_intervals_.loc["dm1_high","95% upper-bound"])
p_d = cph_d.summary.loc["dm1_high","p"]
print(f"  Dichotomized Cox HR (top30 vs bot30) = {hr_d:.2f} [{ci_lo:.2f}, {ci_hi:.2f}], p = {p_d:.2e}")

# === 4. Figure: forest plot of multivariate Cox ===
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# panel A: forest plot
ax = axes[0]
m = cox_df[(cox_df["model"]=="multivariate") & (cox_df["outcome"]=="PFI")].sort_values("HR")
y = np.arange(len(m))
ax.errorbar(m["HR"], y, xerr=[m["HR"]-m["HR_ci_lo"], m["HR_ci_hi"]-m["HR"]],
            fmt="o", color="#962E2E", capsize=4, markersize=10)
ax.axvline(1, color="grey", lw=0.5, ls=":")
ax.set_yticks(y); ax.set_yticklabels(m["covariate"])
ax.set_xlabel("Hazard ratio (PFI)")
ax.set_xscale("log")
for i, (_, r) in enumerate(m.iterrows()):
    ax.text(r["HR_ci_hi"]*1.1, i,
            f"HR {r['HR']:.2f}, p={r['p']:.2e}", va="center", fontsize=9)
ax.set_title(f"O.A — Multivariate Cox PFI (n={int(m['n'].iloc[0])})\n"
             f"DM1_like (per 1 SD) HR = {m[m['covariate']=='dm1_z']['HR'].iloc[0]:.2f}, "
             f"p = {m[m['covariate']=='dm1_z']['p'].iloc[0]:.2e}", fontsize=11)

# panel B: KM dichotomized
ax = axes[1]
for label, sub_df, color in [("DM1_low (bot 30%)", lo, "#3C6B4F"),
                              ("DM1_high (top 30%)", hi, "#962E2E")]:
    kmf = KaplanMeierFitter().fit(sub_df["PFI.time"], sub_df["PFI"], label=f"{label} (n={len(sub_df)})")
    kmf.plot_survival_function(ax=ax, color=color, ci_show=True)
ax.set_title(f"O.B — Dichotomized DM1 PFI (top 30% vs bot 30%)\n"
             f"HR = {hr_d:.2f} [{ci_lo:.2f}, {ci_hi:.2f}], p = {p_d:.2e}, log-rank p = {lr.p_value:.2e}",
             fontsize=11)
ax.set_xlabel("PFI time (days)"); ax.set_ylabel("Progression-free probability")
ax.legend(loc="lower left", frameon=False)
fig.tight_layout()
fig.savefig(OUT / "o_multivariate_cox.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"→ {OUT / 'o_multivariate_cox.png'}")
