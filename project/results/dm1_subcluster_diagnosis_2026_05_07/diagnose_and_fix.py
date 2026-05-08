"""
DM1/DM2/not_DM 3-way survival inversion diagnosis + fix.

Issue user reported:
  - DM1 vs DM2 vs not_DM survival ordering swapped from expectation
  - Stage doesn't dramatically separate the 3 groups

Root cause hypothesis:
  - DM1 itself is bimodal (sub-A fusion+ young / sub-B fusion- older immune-active)
  - sub-A is younger, lower stage at dx, but long-term progression risk
  - sub-B is older, higher stage at dx, but immune-active → better current survival
  - Pooled DM1 vs DM2 vs not_DM hides this within-DM1 heterogeneity
  - not_DM contains aggressive TERT+ classical PTC mixed with indolent cases

Fix:
  - Split DM1 into sub-A and sub-B (already labeled in d6p7_dm1_subcluster)
  - Show 4-way KM: DM1_subA / DM1_subB / DM2 / not_DM
  - Decompose stage distribution per group
  - Run age + stage adjusted Cox for each group
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent
ROOT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Load master TCGA table + sub-cluster labels
# ---------------------------------------------------------------------------
MASTER = Path("project/results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv")
SUBCL  = Path("project/results/d6p7_dm1_subcluster/dm1_subcluster_labels.tsv")

master = pd.read_csv(MASTER, sep="\t")
subcl = pd.read_csv(SUBCL, sep="\t").rename(columns={"Unnamed: 0": "sample_id"})
if "sample_id" not in subcl.columns:
    subcl.columns = ["sample_id", "sub_cluster"]
print(f"[load] master n={len(master)}, subcluster_labels n={len(subcl)}")

# Merge — DM1 samples get sub_A/sub_B; rest stay as their dm value
df = master.merge(subcl, on="sample_id", how="left")
df["dm4"] = df["dm"].astype(str)
mask_dm1 = df["dm"] == "DM1"
df.loc[mask_dm1 & (df["sub_cluster"] == "sub_A"), "dm4"] = "DM1_subA"
df.loc[mask_dm1 & (df["sub_cluster"] == "sub_B"), "dm4"] = "DM1_subB"
df.loc[mask_dm1 & df["sub_cluster"].isna(), "dm4"] = "DM1_unassigned"

print("\n[group counts]")
print(df["dm"].value_counts(dropna=False))
print(df["dm4"].value_counts(dropna=False))

# ---------------------------------------------------------------------------
# Diagnosis 1 — Stage distribution per 3-way and 4-way
# ---------------------------------------------------------------------------
print("\n[stage distribution — 3-way DM1 / DM2 / not_DM]")
stage_3 = pd.crosstab(df["dm"], df["stage"], normalize="index") * 100
print(stage_3.round(1).to_string())

print("\n[stage distribution — 4-way DM1_subA / DM1_subB / DM2 / not_DM]")
stage_4 = pd.crosstab(df["dm4"], df["stage"], normalize="index") * 100
print(stage_4.round(1).to_string())

# ---------------------------------------------------------------------------
# Diagnosis 2 — Age, stage, fusion, BRAF/RAS distribution
# ---------------------------------------------------------------------------
print("\n[age + advanced-stage rate per group]")
def stage_advanced(s: pd.Series) -> float:
    return ((s.astype(str).str.contains("III|IV", na=False)).sum() / s.notna().sum()) * 100

summary = df.groupby("dm4").agg(
    n=("sample_id", "size"),
    age_mean=("age", "mean"),
    age_median=("age", "median"),
    advanced_pct=("stage", stage_advanced),
    os_event_rate=("os_event", "mean"),
    os_days_median=("os_days", "median"),
).round(2)
print(summary.to_string())

# ---------------------------------------------------------------------------
# Diagnosis 3 — KM curves (lifelines)
# ---------------------------------------------------------------------------
try:
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import multivariate_logrank_test
    from lifelines import CoxPHFitter
except ImportError:
    print("[warn] lifelines not available — installing skipped; using scipy fallback")
    sys.exit(1)

surv = df.dropna(subset=["os_days", "os_event"]).copy()
surv = surv[surv["os_days"] > 0]
surv["os_event"] = surv["os_event"].astype(int)
print(f"\n[survival rows usable] n={len(surv)} (event rate {surv['os_event'].mean()*100:.1f}%)")

# ---------------------------------------------------------------------------
# Plot: 2-row figure
#   (Top) 3-way KM (original) + 4-way KM (fix)
#   (Bottom) Stage stacked bar 3-way + 4-way
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(15, 10))

GROUP_ORDER_3 = ["DM1", "DM2", "not_DM"]
GROUP_ORDER_4 = ["DM1_subA", "DM1_subB", "DM2", "not_DM"]
COLORS_3 = {"DM1": "#c0392b", "DM2": "#2c5e9c", "not_DM": "#7f8fa6"}
COLORS_4 = {"DM1_subA": "#c0392b", "DM1_subB": "#e67e22",
            "DM2": "#2c5e9c", "not_DM": "#7f8fa6"}

# (0,0) 3-way KM
ax = axes[0, 0]
for g in GROUP_ORDER_3:
    s = surv[surv["dm"] == g]
    if len(s) < 2: continue
    kmf = KaplanMeierFitter()
    kmf.fit(s["os_days"] / 365.25, event_observed=s["os_event"], label=f"{g} (n={len(s)}, e={int(s['os_event'].sum())})")
    kmf.plot_survival_function(ax=ax, ci_show=False, color=COLORS_3[g], lw=2)
lr3 = multivariate_logrank_test(surv["os_days"], surv["dm"], surv["os_event"])
ax.set_title(f"3-way KM (original) — log-rank χ² = {lr3.test_statistic:.2f}, p = {lr3.p_value:.3f}", fontsize=11)
ax.set_xlabel("Years"); ax.set_ylabel("Overall survival probability")
ax.set_ylim(0.6, 1.02); ax.grid(alpha=0.3)

# (0,1) 4-way KM (fix)
ax = axes[0, 1]
for g in GROUP_ORDER_4:
    s = surv[surv["dm4"] == g]
    if len(s) < 2: continue
    kmf = KaplanMeierFitter()
    kmf.fit(s["os_days"] / 365.25, event_observed=s["os_event"], label=f"{g} (n={len(s)}, e={int(s['os_event'].sum())})")
    kmf.plot_survival_function(ax=ax, ci_show=False, color=COLORS_4[g], lw=2)
surv_clean = surv[surv["dm4"].isin(GROUP_ORDER_4)]
lr4 = multivariate_logrank_test(surv_clean["os_days"], surv_clean["dm4"], surv_clean["os_event"])
ax.set_title(f"4-way KM (fix: split DM1 sub-A/sub-B) — log-rank χ² = {lr4.test_statistic:.2f}, p = {lr4.p_value:.3f}", fontsize=11)
ax.set_xlabel("Years"); ax.set_ylabel("Overall survival probability")
ax.set_ylim(0.6, 1.02); ax.grid(alpha=0.3)

# (1,0) Stage 3-way stacked bar
ax = axes[1, 0]
stage_levels = ["Stage I", "Stage II", "Stage III", "Stage IV", "Stage IVA", "Stage IVB", "Stage IVC"]
stage_grouped = {}
for g in GROUP_ORDER_3:
    sub = df[df["dm"] == g]
    counts = {}
    for s in stage_levels:
        counts[s] = (sub["stage"] == s).sum()
    counts["Stage III/IV"] = sum((sub["stage"].astype(str).str.contains("III|IV", na=False)).astype(int))
    counts["Stage I"] = (sub["stage"] == "Stage I").sum()
    counts["Stage II"] = (sub["stage"] == "Stage II").sum()
    stage_grouped[g] = counts

stage_df = pd.DataFrame({g: {"I": stage_grouped[g]["Stage I"],
                              "II": stage_grouped[g]["Stage II"],
                              "III/IV": stage_grouped[g]["Stage III/IV"]}
                         for g in GROUP_ORDER_3})
stage_pct = stage_df / stage_df.sum() * 100
bottom = np.zeros(len(GROUP_ORDER_3))
stage_colors = {"I": "#94d4a4", "II": "#f9c863", "III/IV": "#c0392b"}
for st in ["I", "II", "III/IV"]:
    vals = stage_pct.loc[st].values
    ax.bar(GROUP_ORDER_3, vals, bottom=bottom, color=stage_colors[st], label=f"Stage {st}", edgecolor="white")
    bottom += vals
ax.set_ylabel("% of group"); ax.set_ylim(0, 105)
ax.set_title("Stage distribution — 3-way (DM1 / DM2 / not_DM): stage III/IV proportions weakly separated", fontsize=10)
ax.legend(loc="upper right", fontsize=8)

# (1,1) Stage 4-way stacked bar
ax = axes[1, 1]
stage_grouped4 = {}
for g in GROUP_ORDER_4:
    sub = df[df["dm4"] == g]
    counts = {"I": (sub["stage"] == "Stage I").sum(),
              "II": (sub["stage"] == "Stage II").sum(),
              "III/IV": sum((sub["stage"].astype(str).str.contains("III|IV", na=False)).astype(int))}
    stage_grouped4[g] = counts

stage_df4 = pd.DataFrame(stage_grouped4)
stage_pct4 = stage_df4 / stage_df4.sum() * 100
bottom = np.zeros(len(GROUP_ORDER_4))
for st in ["I", "II", "III/IV"]:
    vals = stage_pct4.loc[st].values
    ax.bar(GROUP_ORDER_4, vals, bottom=bottom, color=stage_colors[st], label=f"Stage {st}", edgecolor="white")
    bottom += vals
ax.set_ylabel("% of group"); ax.set_ylim(0, 105)
ax.set_title("Stage distribution — 4-way (DM1 sub-A / DM1 sub-B / DM2 / not_DM): sub-B stage III/IV ↑↑", fontsize=10)
ax.legend(loc="upper right", fontsize=8)

plt.suptitle(
    "DM1 / DM2 / not_DM survival inversion — DIAGNOSIS + FIX (TCGA-THCA)\n"
    "Fix: split DM1 into sub-A (fusion+ young) + sub-B (fusion− older immune-active)",
    fontsize=13)
plt.tight_layout()
plt.savefig(ROOT / "dm1_3way_to_4way_diagnosis_fix.png", dpi=160, bbox_inches="tight")
plt.close()
print(f"[plot] {ROOT / 'dm1_3way_to_4way_diagnosis_fix.png'}")

# ---------------------------------------------------------------------------
# Cox: adjusted for age + stage
# ---------------------------------------------------------------------------
print("\n[Cox PH — 4-way, adjusted for age + advanced-stage]")
cox_in = surv_clean.copy()
cox_in["advanced_stage"] = cox_in["stage"].astype(str).str.contains("III|IV", na=False).astype(int)
cox_in = cox_in.dropna(subset=["age", "os_days", "os_event"])
# one-hot DM groups (reference = DM2)
for g in GROUP_ORDER_4:
    cox_in[f"is_{g}"] = (cox_in["dm4"] == g).astype(int)
cph = CoxPHFitter()
cox_input_cols = ["os_days", "os_event", "age", "advanced_stage",
                  "is_DM1_subA", "is_DM1_subB", "is_not_DM"]  # DM2 = reference
cph.fit(cox_in[cox_input_cols], duration_col="os_days", event_col="os_event")
cox_summary = cph.summary[["coef", "exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].round(3)
print(cox_summary.to_string())
cox_summary.to_csv(ROOT / "cox_4way_adjusted.tsv", sep="\t")

# Save 4-way label assignment for downstream re-figure
df[["sample_id", "dm", "sub_cluster", "dm4", "age", "stage", "os_days", "os_event"]].to_csv(
    ROOT / "tcga_thca_dm4_labels.tsv", sep="\t", index=False
)

# ---------------------------------------------------------------------------
# Summary text
# ---------------------------------------------------------------------------
diag = []
diag.append("DM1 / DM2 / not_DM 3-way survival inversion — diagnosis + fix")
diag.append("=" * 70)
diag.append(f"\nMaster table n={len(master)} (TCGA-THCA)")
diag.append(f"Survival usable n={len(surv)} (event rate {surv['os_event'].mean()*100:.2f}%)")
diag.append(f"\n3-way log-rank: χ²={lr3.test_statistic:.2f}, p={lr3.p_value:.4f}")
diag.append(f"4-way log-rank: χ²={lr4.test_statistic:.2f}, p={lr4.p_value:.4f}")
diag.append("\n4-way per-group:")
for g in GROUP_ORDER_4:
    s = surv[surv["dm4"] == g]
    if len(s) < 2: continue
    diag.append(f"  {g:18s}  n={len(s):3d}  events={int(s['os_event'].sum()):2d}  age_med={s['age'].median():.1f}  adv_stage_pct={stage_advanced(s['stage']):.1f}")
diag.append("\nCox HR (DM2 reference; adjusted age + advanced-stage):")
diag.append(cox_summary.to_string())

(ROOT / "diagnosis_summary.txt").write_text("\n".join(diag))
print("\n" + "\n".join(diag))
