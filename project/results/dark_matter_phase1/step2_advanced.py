"""Step 2 advanced rescues: multivariate Cox + Xing-rescue PFI events + KM curves."""
from pathlib import Path
import json

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy import stats

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
OUT = ROOT / "results" / "dark_matter_phase1"

df = pd.read_csv(OUT / "tcga_dm_master_with_pfi.tsv", sep="\t")
print(f"Loaded {len(df)} rows. Stage column: {df['ajcc_pathologic_tumor_stage'].notna().sum() if 'ajcc_pathologic_tumor_stage' in df else 'absent'}")

# Encode stage
def stage_ord(s):
    if not isinstance(s, str): return np.nan
    s = s.strip().upper()
    if "IV" in s: return 4
    if "III" in s: return 3
    if "II" in s and "III" not in s and "IV" not in s: return 2
    if "I" in s: return 1
    return np.nan
df["stage_ord"] = df["ajcc_pathologic_tumor_stage"].apply(stage_ord)
df["age_yrs"] = pd.to_numeric(df["age"], errors="coerce")

# ============================================================================
# 1. Multivariate Cox: cluster + age + stage on PFI within DM
# ============================================================================
print("\n=== MULTIVARIATE COX (DM-only, PFI ~ cluster + age + stage) ===")
dm = df[df["dm_status"] & df["v17_dark_cluster"].isin(["DM1", "DM2"])].copy()
dm["cluster_dm2"] = (dm["v17_dark_cluster"] == "DM2").astype(int)
sub = dm[["PFI", "PFI.time", "cluster_dm2", "age_yrs", "stage_ord"]].dropna()
print(f"n={len(sub)}, events={int(sub['PFI'].sum())}")
if int(sub["PFI"].sum()) >= 3:
    cph = CoxPHFitter(penalizer=0.01)
    cph.fit(sub.rename(columns={"PFI": "event", "PFI.time": "time"}),
            duration_col="time", event_col="event")
    print(cph.summary[["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].round(3).to_string())

# ============================================================================
# 2. Xing-rescue PFI events: among Xing low-risk patients, compare DM2 vs DM1/none
# ============================================================================
print("\n=== XING-RESCUE PFI: low-risk Xing × DM-cluster ===")
df["xing_low_risk"] = (~df["has_braf_v600e"].fillna(False).astype(bool)) | (~df["tert_pos"].fillna(False).astype(bool))
# Actually Xing low-risk = NOT both BRAF+TERT+
df["xing_low_risk"] = ~(df["has_braf_v600e"].fillna(False).astype(bool) & df["tert_pos"].fillna(False).astype(bool))

# Pick Xing-low-risk patients with cluster + PFI
low = df[df["xing_low_risk"] & df["v17_dark_cluster"].isin(["DM1", "DM2"]) & df["PFI"].notna()].copy()
low["cluster_dm2"] = (low["v17_dark_cluster"] == "DM2").astype(int)
n_low = len(low); ev_low = int(low["PFI"].sum())
n_dm2_in_low = int(low["cluster_dm2"].sum())
ev_dm1 = int(low.loc[low["cluster_dm2"] == 0, "PFI"].sum())
ev_dm2 = int(low.loc[low["cluster_dm2"] == 1, "PFI"].sum())
print(f"Xing low-risk + clustered: n={n_low}, total PFI events={ev_low}")
print(f"  DM1: n={n_low - n_dm2_in_low}, events={ev_dm1} ({ev_dm1/(n_low - n_dm2_in_low)*100:.1f}%)")
print(f"  DM2: n={n_dm2_in_low}, events={ev_dm2} ({ev_dm2/n_dm2_in_low*100:.1f}%)" if n_dm2_in_low else "")

# Logrank
if ev_low >= 3:
    lr = logrank_test(
        low.loc[low["cluster_dm2"] == 0, "PFI.time"], low.loc[low["cluster_dm2"] == 1, "PFI.time"],
        event_observed_A=low.loc[low["cluster_dm2"] == 0, "PFI"], event_observed_B=low.loc[low["cluster_dm2"] == 1, "PFI"]
    )
    print(f"  Logrank p = {lr.p_value:.4f}")

# Cox HR within Xing-low-risk
if ev_low >= 3 and low["cluster_dm2"].nunique() == 2:
    cph2 = CoxPHFitter()
    cph2.fit(low[["PFI.time", "PFI", "cluster_dm2"]].rename(columns={"PFI": "event", "PFI.time": "time"}),
             duration_col="time", event_col="event")
    s2 = cph2.summary.iloc[0]
    print(f"  Cox HR(DM2/DM1 within Xing low-risk): {s2['exp(coef)']:.2f} "
          f"[{s2['exp(coef) lower 95%']:.2f}-{s2['exp(coef) upper 95%']:.2f}], p={s2['p']:.4f}")

# ============================================================================
# 3. DICER1/EIF1AX vs DM2 — does the DICER1+ subset within DM2 drive any signal?
# ============================================================================
print("\n=== DICER1/EIF1AX-positive within DM cohort: PFI events ===")
dm_full = df[df["dm_status"]].copy()
alt = ["DICER1_EIF1AX_PPM1D", "RET", "NTRK_fusion", "RET_fusion", "ALK_fusion", "TP53"]
dm_full["alt_driver_pos"] = dm_full["driver_anchor_v17"].isin(alt)
adv = dm_full[dm_full["PFI"].notna()].copy()
n_alt = int(adv["alt_driver_pos"].sum())
ev_alt_pos = int(adv.loc[adv["alt_driver_pos"], "PFI"].sum())
ev_alt_neg = int(adv.loc[~adv["alt_driver_pos"], "PFI"].sum())
n_alt_neg = int((~adv["alt_driver_pos"]).sum())
print(f"  Alt-driver+ (n={n_alt}): events={ev_alt_pos} ({ev_alt_pos/n_alt*100 if n_alt else 0:.1f}%)")
print(f"  Alt-driver− (n={n_alt_neg}): events={ev_alt_neg} ({ev_alt_neg/n_alt_neg*100 if n_alt_neg else 0:.1f}%)")
if n_alt and n_alt_neg:
    table = [[ev_alt_pos, n_alt - ev_alt_pos], [ev_alt_neg, n_alt_neg - ev_alt_neg]]
    odds, p = stats.fisher_exact(table)
    print(f"  Fisher OR={odds:.2f}, p={p:.4f}")

# ============================================================================
# 4. Save everything to JSON
# ============================================================================
result = {
    "dm_pfi_multivariate_done": int(sub["PFI"].sum()) >= 3 if "sub" in dir() else False,
    "xing_lowrisk_DM_split": {
        "n_total": n_low, "events_total": ev_low,
        "DM1_n": n_low - n_dm2_in_low, "DM1_events": ev_dm1,
        "DM2_n": n_dm2_in_low, "DM2_events": ev_dm2,
    },
    "alt_driver_pfi": {
        "alt_pos_n": n_alt, "alt_pos_events": ev_alt_pos,
        "alt_neg_n": n_alt_neg, "alt_neg_events": ev_alt_neg,
    },
}
(OUT / "step2_advanced.json").write_text(json.dumps(result, indent=2, default=str))
print(f"\nSaved step2_advanced.json")
