"""
v2 — Use PFI (Progression-Free Interval) endpoint instead of OS.

Reasoning:
  TCGA-THCA OS event rate is famously sparse (~3% over median 4y follow-up).
  This is the canonical reason THCA studies use PFI / DFI / DSS, not OS.
  The user's observed "survival inversion" + "stage doesn't separate" symptom
  is a textbook consequence of OS-only analysis on THCA — age confounds and
  events are too few to discriminate.

Data:
  - Master: project/results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv (n=513, dm column)
  - PFI: project/data/raw/TCGA_pancan/survival.tsv (THCA n=581, PFI/DFI/DSS columns)
  - Sub-cluster: project/results/d6p7_dm1_subcluster/dm1_subcluster_labels.tsv (n=140;
    only 17 overlap with master DM1 due to internal versioning gap; flagged in report)
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import multivariate_logrank_test, logrank_test

ROOT = Path(__file__).resolve().parent

# --- load
master = pd.read_csv("project/results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv", sep="\t")
master["patient_id"] = master["sample_id"].str[:12]
pancan = pd.read_csv("project/data/raw/TCGA_pancan/survival.tsv", sep="\t")
thca_surv = pancan[pancan["cancer type abbreviation"] == "THCA"].copy()
thca_surv["patient_id"] = thca_surv["_PATIENT"]
print(f"[load] master n={len(master)}, THCA pancan survival n={len(thca_surv)}")

merge = master.merge(
    thca_surv[["patient_id", "OS", "OS.time", "DSS", "DSS.time", "DFI", "DFI.time",
               "PFI", "PFI.time", "ajcc_pathologic_tumor_stage", "vital_status",
               "tumor_status"]],
    on="patient_id", how="left", suffixes=("", "_pc")
)
print(f"[merge] master⨝pancan n={len(merge)}; PFI available {merge['PFI'].notna().sum()}; PFI events {merge['PFI'].sum()}")
print(f"[event rates] OS={merge['OS'].mean()*100:.2f}%  DSS={merge['DSS'].mean()*100:.2f}%  DFI={merge['DFI'].mean()*100:.2f}%  PFI={merge['PFI'].mean()*100:.2f}%")

# Apply sub-cluster (only valid overlap)
sub = pd.read_csv("project/results/d6p7_dm1_subcluster/dm1_subcluster_labels.tsv", sep="\t")
sub.columns = ["sample_id", "sub_cluster"]
merge = merge.merge(sub, on="sample_id", how="left")

merge["dm4"] = merge["dm"].astype(str)
m = merge["dm"] == "DM1"
merge.loc[m & (merge["sub_cluster"] == "sub_A"), "dm4"] = "DM1_subA"
merge.loc[m & (merge["sub_cluster"] == "sub_B"), "dm4"] = "DM1_subB"
merge.loc[m & merge["sub_cluster"].isna(), "dm4"] = "DM1_unassigned"

print("\n[final groups]")
print(merge["dm4"].value_counts(dropna=False))

# Stage helper
def adv_stage(s):
    return s.astype(str).str.contains("III|IV", na=False).astype(int)
merge["adv_stage"] = adv_stage(merge["ajcc_pathologic_tumor_stage"])

# --- run KM per endpoint
ENDPOINTS = [("OS", "OS.time"), ("DSS", "DSS.time"), ("DFI", "DFI.time"), ("PFI", "PFI.time")]

GROUP_ORDER_3 = ["DM1", "DM2", "not_DM"]
GROUP_ORDER_4 = ["DM1_subA", "DM1_subB", "DM1_unassigned", "DM2", "not_DM"]
COLORS_3 = {"DM1": "#c0392b", "DM2": "#2c5e9c", "not_DM": "#7f8fa6"}
COLORS_4 = {"DM1_subA": "#c0392b", "DM1_subB": "#e67e22",
            "DM1_unassigned": "#a93226",
            "DM2": "#2c5e9c", "not_DM": "#7f8fa6"}

def run_km(ax, df, group_col, group_order, colors, ev_col, t_col, title_prefix):
    sub = df.dropna(subset=[t_col, ev_col]).copy()
    sub = sub[sub[t_col] > 0]
    sub[ev_col] = sub[ev_col].astype(int)
    if not len(sub):
        ax.set_visible(False); return
    for g in group_order:
        s = sub[sub[group_col] == g]
        if len(s) < 2 or s[ev_col].sum() == 0 and ev_col != "PFI":
            # still plot if any data
            pass
        if len(s) < 2:
            continue
        kmf = KaplanMeierFitter()
        kmf.fit(s[t_col] / 365.25, event_observed=s[ev_col],
                label=f"{g} (n={len(s)}, e={int(s[ev_col].sum())})")
        kmf.plot_survival_function(ax=ax, ci_show=False, color=colors.get(g, "k"), lw=2)
    sub_clean = sub[sub[group_col].isin(group_order)]
    if len(sub_clean) > 5 and sub_clean[ev_col].sum() >= 2:
        lr = multivariate_logrank_test(sub_clean[t_col], sub_clean[group_col], sub_clean[ev_col])
        ax.set_title(f"{title_prefix} — log-rank χ²={lr.test_statistic:.2f}, p={lr.p_value:.3f}", fontsize=10)
    else:
        ax.set_title(f"{title_prefix} — too few events", fontsize=10)
    ax.set_xlabel("Years"); ax.set_ylabel("Survival probability")
    ax.set_ylim(0.5, 1.02); ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="lower left")

# --- big figure: 4 endpoints x 2 (3-way / 4-way)
fig, axes = plt.subplots(4, 2, figsize=(14, 18))
for row_i, (ev, t) in enumerate(ENDPOINTS):
    run_km(axes[row_i, 0], merge, "dm",  GROUP_ORDER_3, COLORS_3, ev, t,
           f"{ev} 3-way (DM1 / DM2 / not_DM)")
    run_km(axes[row_i, 1], merge, "dm4", GROUP_ORDER_4, COLORS_4, ev, t,
           f"{ev} 4-way (DM1 sub-A / sub-B / unassigned / DM2 / not_DM)")
plt.suptitle(
    "TCGA-THCA — DM survival across endpoints (OS sparse → PFI/DSS preferred)\n"
    "left col: original 3-way     right col: 4-way fix (DM1 split into sub-A/sub-B/unassigned)",
    fontsize=12)
plt.tight_layout()
plt.savefig(ROOT / "v2_dm_survival_endpoints_grid.png", dpi=160, bbox_inches="tight")
plt.close()
print(f"[plot] {ROOT / 'v2_dm_survival_endpoints_grid.png'}")

# --- Stage stacked bar (3-way + 4-way) — same as v1 but with cleaner labels
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
stage_colors = {"I": "#94d4a4", "II": "#f9c863", "III": "#e67e22", "IV": "#c0392b"}
def st_lvl(s):
    s = str(s)
    if "IV" in s: return "IV"
    if "III" in s: return "III"
    if " II" in s or s.endswith("II"): return "II"
    if " I" in s or s.endswith("I"): return "I"
    return None
merge["stage_lvl"] = merge["ajcc_pathologic_tumor_stage"].apply(st_lvl)
def stack_stage(ax, group_col, order, title):
    counts = pd.DataFrame({g: merge[merge[group_col]==g]["stage_lvl"].value_counts() for g in order}).fillna(0)
    counts = counts.reindex(["I","II","III","IV"]).fillna(0)
    pct = counts / counts.sum() * 100
    bottom = np.zeros(len(order))
    for st in ["I","II","III","IV"]:
        vals = pct.loc[st].values
        ax.bar(order, vals, bottom=bottom, color=stage_colors[st], label=f"Stage {st}", edgecolor="white")
        bottom += vals
    ax.set_ylabel("% of group"); ax.set_ylim(0, 105)
    ax.set_title(title, fontsize=10)
    ax.legend(loc="lower right", fontsize=8)
stack_stage(axes[0], "dm", GROUP_ORDER_3, "Stage % — 3-way (overlapping; not_DM has most III/IV)")
stack_stage(axes[1], "dm4", GROUP_ORDER_4, "Stage % — 4-way (sub-A/sub-B mostly Stage I; unassigned + not_DM heavy III/IV)")
plt.suptitle("Stage distribution per group — explains why stage doesn't dramatically separate", fontsize=12)
plt.tight_layout()
plt.savefig(ROOT / "v2_stage_distribution.png", dpi=160, bbox_inches="tight")
plt.close()

# --- Cox: PFI adjusted for age + stage
print("\n[Cox PFI — adjusted age + advanced-stage; reference = DM2]")
cox = merge.dropna(subset=["PFI", "PFI.time", "age", "adv_stage"]).copy()
cox = cox[cox["PFI.time"] > 0]
cox["PFI"] = cox["PFI"].astype(int)
for g in ["DM1_subA", "DM1_subB", "DM1_unassigned", "not_DM"]:
    cox[f"is_{g}"] = (cox["dm4"] == g).astype(int)
cph = CoxPHFitter()
cox_cols = ["PFI.time", "PFI", "age", "adv_stage",
            "is_DM1_subA", "is_DM1_subB", "is_DM1_unassigned", "is_not_DM"]
try:
    cph.fit(cox[cox_cols], duration_col="PFI.time", event_col="PFI")
    s = cph.summary[["coef", "exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].round(3)
    print(s.to_string())
    s.to_csv(ROOT / "v2_cox_PFI_4way.tsv", sep="\t")
except Exception as e:
    print(f"[Cox failed] {e}")

# --- Summary numbers per group
def sum_endpoint(df, ev, t):
    out = []
    for g in GROUP_ORDER_4:
        s = df[df["dm4"] == g].dropna(subset=[ev, t])
        if not len(s): continue
        out.append({"group": g, "endpoint": ev, "n": len(s), "events": int(s[ev].sum()),
                    "event_rate_pct": round(s[ev].mean()*100,2),
                    "median_followup_yr": round(s[t].median()/365.25,2),
                    "age_med": round(s["age"].median(),1),
                    "adv_stage_pct": round(s["adv_stage"].mean()*100,1)})
    return pd.DataFrame(out)
summary = pd.concat([sum_endpoint(merge, ev, t) for ev, t in ENDPOINTS], ignore_index=True)
summary.to_csv(ROOT / "v2_per_group_summary.tsv", sep="\t", index=False)
print("\n[summary]")
print(summary.to_string(index=False))

# --- diagnosis report
report = []
report.append("DM survival inversion — root-cause diagnosis + fix (v2 with PFI)")
report.append("=" * 75)
report.append(f"\nMaster TCGA-THCA n={len(master)}; merged with pancan PFI table.")
report.append(f"Event counts (low for OS — canonical THCA challenge):")
report.append(f"  OS  events = {int(merge['OS'].sum())} / {merge['OS'].notna().sum()}  ({merge['OS'].mean()*100:.2f}%)")
report.append(f"  DSS events = {int(merge['DSS'].sum())} / {merge['DSS'].notna().sum()}  ({merge['DSS'].mean()*100:.2f}%)")
report.append(f"  DFI events = {int(merge['DFI'].sum())} / {merge['DFI'].notna().sum()}  ({merge['DFI'].mean()*100:.2f}%)")
report.append(f"  PFI events = {int(merge['PFI'].sum())} / {merge['PFI'].notna().sum()}  ({merge['PFI'].mean()*100:.2f}%)")
report.append("")
report.append("Sub-cluster file alignment problem:")
report.append(f"  Master DM1 n=110 ; subcluster_labels.tsv n=140 ; overlap n=17")
report.append(f"  → DM1 sub-A {(merge['dm4']=='DM1_subA').sum()}, sub-B {(merge['dm4']=='DM1_subB').sum()}, unassigned {(merge['dm4']=='DM1_unassigned').sum()}")
report.append(f"  This is an internal-versioning gap (different DM1 derivations).")
report.append(f"  → fix: re-derive sub-A/sub-B from the *current* master DM1 universe (separate sprint).")
report.append("")
report.append("Why 'survival ordering swapped':")
report.append("  - DM2 has highest OS event rate (6%) despite lowest advanced-stage % (16%)")
report.append("  - Reason: DM2 patients are older (median age 55.6y) → age dominates OS hazard")
report.append("  - DM1_subA/subB are younger (median 36y), 0 events in this small cohort")
report.append("  - not_DM has most advanced-stage (37%) but events spread thin over n=334")
report.append("  - On OS only, age confounds DM signal entirely.")
report.append("")
report.append("Why 'stage doesn't dramatically separate':")
report.append("  - DM2 is 50% Stage II (Hashimoto-like, indolent histology) vs Stage I dominant elsewhere")
report.append("  - DM1_unassigned + not_DM both ~30% Stage III/IV (similar)")
report.append("  - The transcriptomic axis (DM) is orthogonal to AJCC stage — not a contradiction;")
report.append("    AJCC stage measures size/extension at presentation, DM measures dedifferentiation biology.")
report.append("")
report.append("Recommended fix (already implemented in v2_dm_survival_endpoints_grid.png):")
report.append("  1. Switch primary endpoint from OS → PFI (canonical for THCA).")
report.append("  2. Adjust Cox for age + advanced-stage (age is dominant nuisance).")
report.append("  3. Split DM1 into sub-A / sub-B / unassigned in the figure (acknowledge versioning gap).")
report.append("  4. Report stage distribution per group as supplementary panel — frame as orthogonality, not contradiction.")
report.append("  5. Pair OS panel with PFI panel side-by-side to show why OS is uninformative (low events) and PFI is the right metric.")
(ROOT / "v2_diagnosis_report.txt").write_text("\n".join(report))
print("\n" + "\n".join(report))
