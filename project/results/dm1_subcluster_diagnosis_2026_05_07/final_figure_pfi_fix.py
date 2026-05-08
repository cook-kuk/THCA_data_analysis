"""Final 3-panel manuscript-ready figure: OS-before / PFI-after / stage-explain."""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

ROOT = Path(__file__).resolve().parent

# load merged data with PFI etc.
master = pd.read_csv("project/results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv", sep="\t")
master["patient_id"] = master["sample_id"].str[:12]
pancan = pd.read_csv("project/data/raw/TCGA_pancan/survival.tsv", sep="\t")
thca = pancan[pancan["cancer type abbreviation"] == "THCA"].copy()
thca["patient_id"] = thca["_PATIENT"]

m = master.merge(thca[["patient_id","OS","OS.time","PFI","PFI.time","ajcc_pathologic_tumor_stage"]],
                 on="patient_id", how="left")

sub = pd.read_csv("project/results/d6p7_dm1_subcluster/dm1_subcluster_labels.tsv", sep="\t")
sub.columns = ["sample_id","sub_cluster"]
m = m.merge(sub, on="sample_id", how="left")

m["dm4"] = m["dm"].astype(str)
isdm1 = m["dm"]=="DM1"
m.loc[isdm1 & (m["sub_cluster"]=="sub_A"), "dm4"] = "DM1_subA"
m.loc[isdm1 & (m["sub_cluster"]=="sub_B"), "dm4"] = "DM1_subB"
m.loc[isdm1 & m["sub_cluster"].isna(), "dm4"] = "DM1_other"

def stage_lvl(s):
    s = str(s)
    if "IV" in s: return "IV"
    if "III" in s: return "III"
    if " II" in s or s.endswith("II"): return "II"
    if " I" in s or s.endswith("I"): return "I"
    return None
m["stage_lvl"] = m["ajcc_pathologic_tumor_stage"].apply(stage_lvl)

GROUP4 = ["DM1_subA", "DM1_subB", "DM1_other", "DM2", "not_DM"]
COLORS4 = {"DM1_subA":"#c0392b","DM1_subB":"#e67e22","DM1_other":"#a93226",
           "DM2":"#2c5e9c","not_DM":"#7f8fa6"}

fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

# ---------- Panel A — OS 4-way (BEFORE — sparse) ----------
ax = axes[0]
sd = m.dropna(subset=["OS","OS.time"]).copy(); sd = sd[sd["OS.time"]>0]; sd["OS"]=sd["OS"].astype(int)
for g in GROUP4:
    s = sd[sd["dm4"]==g]
    if len(s)<2: continue
    kmf = KaplanMeierFitter()
    kmf.fit(s["OS.time"]/365.25, event_observed=s["OS"],
            label=f"{g} (n={len(s)}, e={int(s['OS'].sum())})")
    kmf.plot_survival_function(ax=ax, ci_show=False, color=COLORS4[g], lw=2)
sd_in = sd[sd["dm4"].isin(GROUP4)]
lr = multivariate_logrank_test(sd_in["OS.time"], sd_in["dm4"], sd_in["OS"])
ax.set_title(f"A. OS (BEFORE — uninformative)\nlog-rank χ²={lr.test_statistic:.2f}, p={lr.p_value:.3f}    total events={int(sd['OS'].sum())} / n={len(sd)}",
             fontsize=10)
ax.set_xlabel("Years"); ax.set_ylabel("Overall survival probability")
ax.set_ylim(0.7, 1.02); ax.grid(alpha=0.3); ax.legend(fontsize=7, loc="lower left")

# ---------- Panel B — PFI 4-way (FIX) ----------
ax = axes[1]
sd = m.dropna(subset=["PFI","PFI.time"]).copy(); sd = sd[sd["PFI.time"]>0]; sd["PFI"]=sd["PFI"].astype(int)
for g in GROUP4:
    s = sd[sd["dm4"]==g]
    if len(s)<2: continue
    kmf = KaplanMeierFitter()
    kmf.fit(s["PFI.time"]/365.25, event_observed=s["PFI"],
            label=f"{g} (n={len(s)}, e={int(s['PFI'].sum())})")
    kmf.plot_survival_function(ax=ax, ci_show=False, color=COLORS4[g], lw=2)
sd_in = sd[sd["dm4"].isin(GROUP4)]
lr = multivariate_logrank_test(sd_in["PFI.time"], sd_in["dm4"], sd_in["PFI"])
ax.set_title(f"B. PFI (FIX — canonical THCA endpoint)\nlog-rank χ²={lr.test_statistic:.2f}, p={lr.p_value:.3f}    total events={int(sd['PFI'].sum())} / n={len(sd)}",
             fontsize=10)
ax.set_xlabel("Years"); ax.set_ylabel("Progression-free probability")
ax.set_ylim(0.5, 1.02); ax.grid(alpha=0.3); ax.legend(fontsize=7, loc="lower left")

# ---------- Panel C — Stage distribution + age annotation ----------
ax = axes[2]
stage_colors = {"I":"#94d4a4","II":"#f9c863","III":"#e67e22","IV":"#c0392b"}
counts = pd.DataFrame({g: m[m["dm4"]==g]["stage_lvl"].value_counts() for g in GROUP4}).fillna(0)
counts = counts.reindex(["I","II","III","IV"]).fillna(0)
pct = counts / counts.sum() * 100
bottom = np.zeros(len(GROUP4))
for st in ["I","II","III","IV"]:
    vals = pct.loc[st].values
    ax.bar(GROUP4, vals, bottom=bottom, color=stage_colors[st], label=f"Stage {st}", edgecolor="white")
    bottom += vals
# Annotate median age above each bar
ages = m.groupby("dm4")["age"].median()
for i, g in enumerate(GROUP4):
    ax.text(i, 102, f"age {ages[g]:.0f}", ha="center", fontsize=9, color="#333", fontweight="bold")
ax.set_ylabel("% of group"); ax.set_ylim(0, 110)
ax.set_title("C. Stage % per group + median age\n(stage axis ⊥ DM axis; DM2 older → age confounds OS)", fontsize=10)
ax.legend(loc="lower right", fontsize=8)
ax.tick_params(axis="x", labelrotation=20)

plt.suptitle(
    "DM1 / DM2 / not_DM survival inversion — diagnosis (A) + fix (B) + explanation (C)\n"
    "OS too sparse (events <4%) → PFI separates groups; AJCC stage measures presentation extent ≠ DM dedifferentiation axis",
    fontsize=12)
plt.tight_layout()
plt.savefig(ROOT / "FINAL_dm_survival_pfi_fix.png", dpi=160, bbox_inches="tight")
plt.close()
print(f"[plot] {ROOT / 'FINAL_dm_survival_pfi_fix.png'}")
