#!/usr/bin/env python3
"""N1 — Dark-matter subset analysis (Paper 1 target population).
Filter TCGA-THCA to BRAF-negative AND RAS-negative (= dark-matter), then
recompute stage trend + survival. Effect should be STRONGER in this subset
since Paper 1's primary thesis is about dark-matter biology specifically."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

OUT = Path("project_external_st/results/extra")
df = pd.read_csv(OUT / "s_tcga_thca_scored.tsv", sep="\t")
print(f"All TCGA-THCA: {len(df)}")

# Define dark matter: BRAF-negative AND RAS-negative
df["is_dark_matter"] = (~df["has_braf_v600e"].astype(bool)) & (~df["has_ras_mut"].astype(bool))
dm = df[df["is_dark_matter"]].copy()
print(f"Dark matter (BRAF-neg AND RAS-neg): n = {len(dm)}")
print(f"BRAF+: n = {df['has_braf_v600e'].sum()}, RAS+: n = {df['has_ras_mut'].sum()}")

# Stage trend within dark matter
stage_map = {"Stage I":1,"Stage II":2,"Stage III":3,"Stage IVA":4,"Stage IVB":5,"Stage IVC":6}
dm["stage_ord"] = dm["ajcc_pathologic_tumor_stage"].map(stage_map)
ok = dm.dropna(subset=["stage_ord","DM1_like"])
rho, p = spearmanr(ok["stage_ord"], ok["DM1_like"])
print(f"\n=== Dark matter stage trend (n={len(ok)}) ===")
print(f"  DM1_like vs stage_ord ρ = {rho:.3f}, p = {p:.2e}")
print("\n  per-stage mean:")
print(ok.groupby("ajcc_pathologic_tumor_stage")["DM1_like"].agg(["mean","std","count"]).to_string())

# Compare to all-cohort effect
all_ok = df.dropna(subset=["ajcc_pathologic_tumor_stage","DM1_like"])
all_ok["stage_ord"] = all_ok["ajcc_pathologic_tumor_stage"].map(stage_map)
rho_all, p_all = spearmanr(all_ok["stage_ord"], all_ok["DM1_like"])
print(f"\n  Compared to ALL TCGA-THCA: ρ_all = {rho_all:.3f}, p_all = {p_all:.2e}")
print(f"  Effect size in dark matter / all = {rho/rho_all:.2f}x")

# Survival in dark matter
dm_surv = dm.dropna(subset=["PFI","PFI.time","DM1_like"]).copy()
print(f"\n=== Dark matter PFI (n={len(dm_surv)}) ===")
cph = CoxPHFitter().fit(
    dm_surv[["PFI.time","PFI","DM1_like"]].rename(columns={"PFI.time":"time","PFI":"event"}),
    duration_col="time", event_col="event")
hr = cph.hazard_ratios_["DM1_like"]
ci_lo = np.exp(cph.confidence_intervals_.loc["DM1_like","95% lower-bound"])
ci_hi = np.exp(cph.confidence_intervals_.loc["DM1_like","95% upper-bound"])
p_cox = cph.summary.loc["DM1_like","p"]
print(f"  Continuous Cox HR = {hr:.2f} [{ci_lo:.2f}, {ci_hi:.2f}], p = {p_cox:.2e}")

# Dichotomized HR within dark matter
vals = dm_surv["DM1_like"].astype(float).to_numpy()
q30, q70 = np.percentile(vals, [30, 70])
lo = dm_surv[dm_surv["DM1_like"] <= q30]
hi = dm_surv[dm_surv["DM1_like"] >= q70]
lr = logrank_test(lo["PFI.time"], hi["PFI.time"], lo["PFI"], hi["PFI"])
print(f"  Log-rank top30 vs bot30: p = {lr.p_value:.2e}")
ok2 = pd.concat([lo.assign(grp=0), hi.assign(grp=1)])
cph_d = CoxPHFitter().fit(
    ok2[["PFI.time","PFI","grp"]].rename(columns={"PFI.time":"time","PFI":"event"}),
    duration_col="time", event_col="event")
hr_d = cph_d.hazard_ratios_["grp"]
ci_lo_d = np.exp(cph_d.confidence_intervals_.loc["grp","95% lower-bound"])
ci_hi_d = np.exp(cph_d.confidence_intervals_.loc["grp","95% upper-bound"])
p_d = cph_d.summary.loc["grp","p"]
print(f"  Dichotomized Cox HR (top30 vs bot30) = {hr_d:.2f} [{ci_lo_d:.2f}, {ci_hi_d:.2f}], p = {p_d:.2e}")

# Save summary
summary = pd.DataFrame([
    {"cohort":"All TCGA-THCA", "n":len(all_ok), "stage_rho":rho_all, "stage_p":p_all},
    {"cohort":"Dark matter (BRAF-neg AND RAS-neg)", "n":len(ok), "stage_rho":rho, "stage_p":p,
     "PFI_HR_continuous":hr, "PFI_HR_continuous_p":p_cox,
     "PFI_HR_dichotomy":hr_d, "PFI_HR_dichotomy_p":p_d}
])
summary.to_csv(OUT / "n1_dark_matter_subset.tsv", sep="\t", index=False)

# Figure
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
import seaborn as sns
ax = axes[0]
order = ["Stage I","Stage II","Stage III","Stage IVA","Stage IVC"]
sns.boxplot(data=ok[ok["ajcc_pathologic_tumor_stage"].isin(order)],
            x="ajcc_pathologic_tumor_stage", y="DM1_like", order=order, ax=ax,
            palette=["#3C6B4F","#5C8B6E","#B8893C","#A04451","#962E2E"])
ax.set_title(f"N1.A — Dark-matter stage trend (n={len(ok)})\n"
             f"Spearman ρ = {rho:.3f}, p = {p:.2e}\n"
             f"All-cohort ρ = {rho_all:.3f} (compare)",
             fontsize=11)
ax.set_xlabel(""); ax.set_ylabel("DM1_like")
ax.axhline(0, color="grey", lw=.4, ls=":")

ax = axes[1]
for label, sub_df, color in [("DM1_low (bot 30%)", lo, "#3C6B4F"),
                              ("DM1_high (top 30%)", hi, "#962E2E")]:
    kmf = KaplanMeierFitter().fit(sub_df["PFI.time"], sub_df["PFI"], label=f"{label} (n={len(sub_df)})")
    kmf.plot_survival_function(ax=ax, color=color, ci_show=True)
ax.set_title(f"N1.B — Dark-matter PFI dichotomized (n={len(dm_surv)})\n"
             f"HR = {hr_d:.2f} [{ci_lo_d:.2f}, {ci_hi_d:.2f}], p = {p_d:.2e}\n"
             f"log-rank p = {lr.p_value:.2e}",
             fontsize=11)
ax.set_xlabel("PFI time (days)"); ax.set_ylabel("PFI probability")
ax.legend(loc="lower left", frameon=False, fontsize=9)

ax = axes[2]
# DM1 distribution by mutation group
df["mut_group"] = "BRAF+/RAS+"
df.loc[df["has_braf_v600e"].astype(bool) & ~df["has_ras_mut"].astype(bool), "mut_group"] = "BRAF+ only"
df.loc[~df["has_braf_v600e"].astype(bool) & df["has_ras_mut"].astype(bool), "mut_group"] = "RAS+ only"
df.loc[~df["has_braf_v600e"].astype(bool) & ~df["has_ras_mut"].astype(bool), "mut_group"] = "Dark matter"
sns.boxplot(data=df, x="mut_group", y="DM1_like",
            order=["BRAF+ only","RAS+ only","BRAF+/RAS+","Dark matter"], ax=ax,
            palette=["#962E2E","#B8893C","#7B1F2A","#3C6B4F"])
ax.set_title("N1.C — DM1_like by mutation status\n"
             "Dark matter = paper 1 target subset", fontsize=11)
ax.set_xlabel(""); ax.set_ylabel("DM1_like")
for i, g in enumerate(["BRAF+ only","RAS+ only","BRAF+/RAS+","Dark matter"]):
    n = (df["mut_group"]==g).sum()
    ax.text(i, ax.get_ylim()[0]+0.05, f"n={n}", ha="center", fontsize=9)

fig.tight_layout()
fig.savefig(OUT / "n1_dark_matter.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"→ {OUT / 'n1_dark_matter.png'}")
