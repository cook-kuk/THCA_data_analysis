"""Speed batch: pediatric subset + DICER1+ specific PFI + age-stratified r."""
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
OUT = ROOT / "results/dark_matter_phase2"
FIG = OUT / "web/figures"
DATA_OUT = OUT / "web/data"

df = pd.read_csv(OUT / "../dark_matter_phase1/tcga_dm_master_with_pfi.tsv", sep="\t")
df["age_num"] = pd.to_numeric(df["age"], errors="coerce")
print(f"Total: {len(df)} TCGA samples, age range {df['age_num'].min():.0f}-{df['age_num'].max():.0f}")

# ============================================================================
# 1. Pediatric / young-adult (<35 yo) subset
# ============================================================================
print("\n=== 1. Young-adult <35 subset ===")
young = df[df["age_num"] < 35].copy()
old = df[df["age_num"] >= 55].copy()
print(f"<35: n={len(young)}, ≥55: n={len(old)}")
print(f"<35 driver class:\n{young['mutation_group'].value_counts(dropna=False)}")
print(f"<35 dm_status %: {young['dm_status'].mean()*100:.1f}%")
print(f"≥55 dm_status %: {old['dm_status'].mean()*100:.1f}%")

# Within <35, DM cluster distribution
y_dm = young[young["dm_status"] & young["v17_dark_cluster"].isin(["DM1", "DM2"])].copy()
print(f"<35 + DM with cluster: n={len(y_dm)}, DM1={int((y_dm['v17_dark_cluster']=='DM1').sum())}, DM2={int((y_dm['v17_dark_cluster']=='DM2').sum())}")

# PFI events young vs old
print(f"<35 PFI events: {int(young['PFI'].sum())}/{int(young['PFI'].notna().sum())}")
print(f"≥55 PFI events: {int(old['PFI'].sum())}/{int(old['PFI'].notna().sum())}")

# Age-stratified DM% across stratums
strata = [(0, 35, "<35"), (35, 55, "35-54"), (55, 100, "≥55")]
strat_results = []
for lo, hi, lbl in strata:
    sub = df[(df["age_num"] >= lo) & (df["age_num"] < hi)]
    n = len(sub)
    if n == 0: continue
    n_dm = int(sub["dm_status"].sum())
    n_dm1 = int((sub["v17_dark_cluster"] == "DM1").sum())
    n_dm2 = int((sub["v17_dark_cluster"] == "DM2").sum())
    n_braf = int(sub["has_braf_v600e"].sum())
    n_ras = int(sub["has_ras_mut"].sum())
    pfi_n = int(sub["PFI"].notna().sum())
    pfi_e = int(sub["PFI"].sum())
    strat_results.append({
        "age_band": lbl, "n": n, "DM_n": n_dm, "DM_pct": round(n_dm/n*100, 1),
        "BRAF_pct": round(n_braf/n*100, 1), "RAS_pct": round(n_ras/n*100, 1),
        "DM1_n": n_dm1, "DM2_n": n_dm2,
        "PFI_n": pfi_n, "PFI_events": pfi_e,
        "PFI_rate_pct": round(pfi_e/pfi_n*100, 2) if pfi_n else None,
    })
strat_df = pd.DataFrame(strat_results)
print(f"\n=== Age-stratified driver landscape ===\n{strat_df.to_string(index=False)}")
strat_df.to_csv(DATA_OUT / "age_stratified_landscape.tsv", sep="\t", index=False)

# ============================================================================
# 2. DICER1/EIF1AX-specific PFI (whole cohort + within DM)
# ============================================================================
print("\n=== 2. DICER1/EIF1AX+ specific PFI ===")
df["alt_pos"] = df["driver_anchor_v17"].isin(["DICER1_EIF1AX_PPM1D"])
print(f"DICER1/EIF1AX/PPM1D+ in TCGA: {int(df['alt_pos'].sum())}")
print(f"DICER1+ PFI events: {int(df.loc[df['alt_pos'], 'PFI'].sum())}/{int(df.loc[df['alt_pos'], 'PFI'].notna().sum())}")

# Compare DICER1+ vs DICER1- within DM
dm_full = df[df["dm_status"]].copy()
ev_alt_pos = int(dm_full.loc[dm_full["alt_pos"], "PFI"].sum())
n_alt_pos = int(dm_full["alt_pos"].sum())
ev_alt_neg = int(dm_full.loc[~dm_full["alt_pos"], "PFI"].sum())
n_alt_neg = int((~dm_full["alt_pos"]).sum())
print(f"DICER1+ in DM: {ev_alt_pos}/{n_alt_pos} ({ev_alt_pos/n_alt_pos*100 if n_alt_pos else 0:.1f}%)")
print(f"DICER1- in DM: {ev_alt_neg}/{n_alt_neg} ({ev_alt_neg/n_alt_neg*100 if n_alt_neg else 0:.1f}%)")

# Logrank
sub_alt = dm_full[dm_full["PFI"].notna()].copy()
if int(sub_alt["alt_pos"].sum()) >= 2 and int((~sub_alt["alt_pos"]).sum()) >= 5:
    ap = sub_alt[sub_alt["alt_pos"]]; an = sub_alt[~sub_alt["alt_pos"]]
    lr = logrank_test(ap["PFI.time"], an["PFI.time"], ap["PFI"], an["PFI"])
    print(f"DICER1+ vs - within DM PFI logrank p = {lr.p_value:.4f}")

# ============================================================================
# 3. Age-stratified figure
# ============================================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
ax = axes[0]
ax2 = ax.twinx()
x = np.arange(len(strat_df))
ax.bar(x - 0.2, strat_df["BRAF_pct"], 0.18, label="BRAF V600E %", color="#e63946", alpha=0.85)
ax.bar(x, strat_df["RAS_pct"], 0.18, label="RAS %", color="#f4a261", alpha=0.85)
ax.bar(x + 0.2, strat_df["DM_pct"], 0.18, label="Dark Matter %", color="#2a9d8f", alpha=0.85)
ax.set_xticks(x); ax.set_xticklabels(strat_df["age_band"])
ax.set_xlabel("Age band (years at diagnosis)")
ax.set_ylabel("% of patients")
ax.set_title(f"Figure E15A — Driver landscape by age (TCGA-THCA n={len(df)})")
ax.legend(loc="upper left"); ax.grid(alpha=0.3, axis="y")
ax2.plot(x, strat_df["DM1_n"] / (strat_df["DM1_n"] + strat_df["DM2_n"]).replace(0, np.nan) * 100,
         "o--", color="#4361ee", lw=2, markersize=8, label="DM1 % within DM")
ax2.set_ylabel("DM1 % within DM (line)", color="#4361ee")
ax2.tick_params(axis="y", labelcolor="#4361ee")
ax2.set_ylim(0, 100)
ax2.legend(loc="upper right")

# DICER1 split bar
ax = axes[1]
ax.bar(["DICER1+ (n=" + str(n_alt_pos) + ")", "DICER1- (n=" + str(n_alt_neg) + ")"],
       [ev_alt_pos / n_alt_pos * 100 if n_alt_pos else 0,
        ev_alt_neg / n_alt_neg * 100 if n_alt_neg else 0],
       color=["#9d4edd", "#dee2e6"], alpha=0.85)
ax.set_ylabel("PFI event rate (%)")
ax.set_title(f"Figure E15B — DICER1+ vs DICER1- PFI rate within Dark Matter")
ax.grid(alpha=0.3, axis="y")
fig.tight_layout()
fig.savefig(FIG / "figE15_age_dicer1_subgroup.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved figE15")

# Save summary
result = {
    "young_lt35": {"n": len(young), "DM_pct": float(young["dm_status"].mean()*100), "PFI_events": int(young["PFI"].sum())},
    "old_ge55": {"n": len(old), "DM_pct": float(old["dm_status"].mean()*100), "PFI_events": int(old["PFI"].sum())},
    "age_stratified": strat_df.to_dict(orient="records"),
    "alt_driver_pfi": {"DICER1_pos_n": n_alt_pos, "DICER1_pos_events": ev_alt_pos,
                       "DICER1_neg_n": n_alt_neg, "DICER1_neg_events": ev_alt_neg},
}
(DATA_OUT / "p5_speed_summary.json").write_text(json.dumps(result, indent=2, default=str))
print(json.dumps(result, indent=2))
