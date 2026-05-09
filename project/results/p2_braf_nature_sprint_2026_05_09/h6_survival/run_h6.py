"""
H6 — Clinical outcome translation of DM1/DM2 within BRAF-cPTC stratum.

Endpoints: OS, PFI, DFI, DSS (TCGA pancan canonical), recurrence (new tumor event),
and RAI-related fields (i_131_total_administered_dose, additional_radiation_therapy).

Strata tested (in priority order):
  S1: BRAF_like / cPTC                 [headline]
  S2: BRAF_like (any histology)
  S3: full TCGA-THCA
  S4: BRAF_like / cPTC × TERT⁺/⁻       [4-way / interaction]

Reference DM2 (small RAI-refractory-leaning n) → contrasts focus on DM1 vs DM2 and DM1 vs not_DM.
"""
from __future__ import annotations
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy import stats as sstats

warnings.filterwarnings("ignore")

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p2_braf_nature_sprint_2026_05_09/h6_survival")
OUT.mkdir(parents=True, exist_ok=True)

MASTER = "/data/thca/repo_results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv"
PANCAN = "/home/seungho/personal/THCA_data_analysis/project/data/raw/TCGA_pancan/survival.tsv"

# ----- load -----
master = pd.read_csv(MASTER, sep="\t")
pan = pd.read_csv(PANCAN, sep="\t")
pan = pan[pan["cancer type abbreviation"] == "THCA"].copy()
pan["sample_short"] = pan["sample"].str[:15]  # TCGA-XX-XXXX-01

# canonical sample short
master["sample_short"] = master["sample_id"].str[:15]
mg = master.merge(
    pan[["sample_short", "OS", "OS.time", "DSS", "DSS.time",
         "DFI", "DFI.time", "PFI", "PFI.time", "ajcc_pathologic_tumor_stage"]],
    on="sample_short", how="left", suffixes=("", "_pan"),
)
print(f"[merge] master n={len(master)} ⨝ pancan THCA n={len(pan)} → joined n={len(mg)}")
print(f"  PFI events: {int(mg['PFI'].sum())} / {mg['PFI'].notna().sum()} | "
      f"DFI events: {int(mg['DFI'].sum())} / {mg['DFI'].notna().sum()} | "
      f"DSS events: {int(mg['DSS'].sum())} / {mg['DSS'].notna().sum()} | "
      f"OS events (pancan): {int(mg['OS'].sum())} / {mg['OS'].notna().sum()}")

# stage harmonization → adv = III/IV
def adv_stage(s):
    if pd.isna(s): return np.nan
    s = str(s)
    if "IV" in s or "III" in s: return 1
    return 0
mg["adv_stage"] = mg["ajcc_stage_group"].map(adv_stage)
mg["age_yr"] = pd.to_numeric(mg["age_at_diagnosis"], errors="coerce")
mg["sex_m"] = (mg["gender"].astype(str).str.upper() == "MALE").astype(int)
mg["tert_mut"] = (mg["tert_promoter_integrated"] == "mutated").astype(int)
mg["recurrence"] = mg["new_tumor_event_after_initial_treatment"].map({"YES": 1, "NO": 0})

# binary DM contrasts
mg["is_DM1"] = (mg["dm"] == "DM1").astype(int)
mg["is_DM2"] = (mg["dm"] == "DM2").astype(int)
mg["dm_DM1_vs_DM2"] = np.where(mg["dm"]=="DM1", 1, np.where(mg["dm"]=="DM2", 0, np.nan))
mg["dm_DM1_vs_other"] = np.where(mg["dm"]=="DM1", 1, np.where(mg["dm"].isin(["DM2","not_DM"]), 0, np.nan))

# strata masks
S1 = (mg["molecular_subtype"]=="BRAF_like") & (mg["histology_subtype"]=="cPTC")
S2 = (mg["molecular_subtype"]=="BRAF_like")
S3 = pd.Series(True, index=mg.index)

print(f"\n[strata] S1 BRAF_like/cPTC n={int(S1.sum())} | S2 BRAF_like n={int(S2.sum())} | S3 full n={int(S3.sum())}")
for tag, mk in [("S1", S1), ("S2", S2), ("S3", S3)]:
    sub = mg[mk]
    print(f"  {tag}: dm = {sub['dm'].value_counts().to_dict()} | "
          f"PFI events = {int(sub['PFI'].sum() or 0)} | "
          f"recurrence = {(sub['recurrence']==1).sum()}")

# ============================================================
# Cox table
# ============================================================
ROWS = []

def cox_fit(df, duration, event, covars, label, stratum, endpoint, contrast):
    """Fit a Cox model and emit one row per primary covariate."""
    use = df.dropna(subset=[duration, event] + covars).copy()
    use = use[use[duration] > 0]
    n = len(use); ev = int(use[event].sum())
    if n < 8 or ev < 3:
        ROWS.append({
            "stratum": stratum, "endpoint": endpoint, "model": label,
            "covariate": contrast, "n": n, "events": ev,
            "HR": np.nan, "CI_low": np.nan, "CI_high": np.nan, "p": np.nan,
            "note": "insufficient n/events",
        })
        return None
    cph = CoxPHFitter(penalizer=0.0)
    try:
        cph.fit(use[[duration, event] + covars],
                duration_col=duration, event_col=event, robust=True)
    except Exception as e:
        # try with small penalizer
        try:
            cph = CoxPHFitter(penalizer=0.05)
            cph.fit(use[[duration, event] + covars],
                    duration_col=duration, event_col=event, robust=True)
        except Exception as e2:
            ROWS.append({
                "stratum": stratum, "endpoint": endpoint, "model": label,
                "covariate": contrast, "n": n, "events": ev,
                "HR": np.nan, "CI_low": np.nan, "CI_high": np.nan, "p": np.nan,
                "note": f"fit_failed: {type(e2).__name__}",
            })
            return None
    s = cph.summary
    for cv in covars:
        if cv not in s.index: continue
        ROWS.append({
            "stratum": stratum, "endpoint": endpoint, "model": label,
            "covariate": cv, "n": n, "events": ev,
            "HR": float(s.loc[cv, "exp(coef)"]),
            "CI_low": float(s.loc[cv, "exp(coef) lower 95%"]),
            "CI_high": float(s.loc[cv, "exp(coef) upper 95%"]),
            "p": float(s.loc[cv, "p"]),
            "note": "",
        })
    return cph

# Endpoints to run
ENDPOINTS = [
    ("OS",  "OS.time",  "OS"),
    ("PFI", "PFI.time", "PFI"),
    ("DFI", "DFI.time", "DFI"),
    ("DSS", "DSS.time", "DSS"),
]

STRATA = [("S1_BRAF_cPTC", S1), ("S2_BRAF_any", S2), ("S3_full", S3)]

# DM1 vs DM2 (within DM-positive only)
for sname, mk in STRATA:
    sub = mg[mk & mg["dm"].isin(["DM1","DM2"])].copy()
    for ev_col, t_col, label in ENDPOINTS:
        cox_fit(sub, t_col, ev_col,
                ["is_DM1", "age_yr", "adv_stage", "sex_m"],
                label="DM1_vs_DM2_adj_age_stage_sex", stratum=sname,
                endpoint=label, contrast="is_DM1")

# DM1 vs (DM2 + not_DM)  — wider baseline
for sname, mk in STRATA:
    sub = mg[mk].copy()
    sub = sub.dropna(subset=["dm_DM1_vs_other"])
    sub["is_DM1_wide"] = sub["dm_DM1_vs_other"].astype(int)
    for ev_col, t_col, label in ENDPOINTS:
        cox_fit(sub, t_col, ev_col,
                ["is_DM1_wide", "age_yr", "adv_stage", "sex_m"],
                label="DM1_vs_others_adj_age_stage_sex", stratum=sname,
                endpoint=label, contrast="is_DM1_wide")

# DM2 vs not_DM (positive control: rai_score_v17 says DM2 is high-RAI silencing)
for sname, mk in STRATA:
    sub = mg[mk & mg["dm"].isin(["DM2","not_DM"])].copy()
    sub["is_DM2_only"] = (sub["dm"]=="DM2").astype(int)
    for ev_col, t_col, label in ENDPOINTS:
        cox_fit(sub, t_col, ev_col,
                ["is_DM2_only", "age_yr", "adv_stage", "sex_m"],
                label="DM2_vs_notDM_adj_age_stage_sex", stratum=sname,
                endpoint=label, contrast="is_DM2_only")

# TERT independence test in S1 (BRAF-cPTC): is_DM1 + tert_mut both
for ev_col, t_col, label in ENDPOINTS:
    sub = mg[S1 & mg["dm"].isin(["DM1","DM2","not_DM"])].copy()
    sub["is_DM1_wide"] = (sub["dm"]=="DM1").astype(int)
    cox_fit(sub, t_col, ev_col,
            ["is_DM1_wide", "tert_mut", "age_yr", "adv_stage", "sex_m"],
            label="DM1_AND_TERT_jointly_adj_S1", stratum="S1_BRAF_cPTC",
            endpoint=label, contrast="is_DM1_wide")
    # TERT alone
    cox_fit(sub, t_col, ev_col,
            ["tert_mut", "age_yr", "adv_stage", "sex_m"],
            label="TERT_alone_S1_BRAF_cPTC", stratum="S1_BRAF_cPTC",
            endpoint=label, contrast="tert_mut")

# TERT alone full stratum (replicate the +memory HR=7.57 anchor on OS)
for ev_col, t_col, label in ENDPOINTS:
    cox_fit(mg.copy(), t_col, ev_col,
            ["tert_mut", "age_yr", "adv_stage", "sex_m"],
            label="TERT_alone_full_TCGA", stratum="S3_full",
            endpoint=label, contrast="tert_mut")

# Recurrence Cox on PFI.time as exposure proxy is already PFI; do a logistic-style chi-square +
# Recurrence: use new_tumor_event flag against DM groups.
def fisher_summary(tab, label, stratum, contrast):
    """tab: 2x2 [[a,b],[c,d]]"""
    try:
        odds, p = sstats.fisher_exact(tab)
    except Exception:
        odds, p = (np.nan, np.nan)
    a,b = tab[0]; c,d = tab[1]
    # Wald CI on log odds
    if a>0 and b>0 and c>0 and d>0:
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
        loghr = np.log(odds)
        lo, hi = np.exp(loghr - 1.96*se), np.exp(loghr + 1.96*se)
    else:
        lo = hi = np.nan
    n = a+b+c+d; ev = a+c
    ROWS.append({
        "stratum": stratum, "endpoint": label, "model": "fisher_exact",
        "covariate": contrast, "n": int(n), "events": int(ev),
        "HR": float(odds), "CI_low": float(lo), "CI_high": float(hi),
        "p": float(p), "note": f"a/b/c/d={a}/{b}/{c}/{d}",
    })

# Recurrence
for sname, mk in STRATA:
    sub = mg[mk].dropna(subset=["recurrence"])
    # DM1 vs DM2
    s12 = sub[sub["dm"].isin(["DM1","DM2"])]
    if len(s12)>0:
        a = ((s12["dm"]=="DM1") & (s12["recurrence"]==1)).sum()
        b = ((s12["dm"]=="DM1") & (s12["recurrence"]==0)).sum()
        c = ((s12["dm"]=="DM2") & (s12["recurrence"]==1)).sum()
        d = ((s12["dm"]=="DM2") & (s12["recurrence"]==0)).sum()
        fisher_summary([[a,b],[c,d]], "recurrence", sname, "DM1_vs_DM2_OR")
    # DM1 vs (DM2+not_DM)
    s1w = sub[sub["dm"].isin(["DM1","DM2","not_DM"])]
    a = ((s1w["dm"]=="DM1") & (s1w["recurrence"]==1)).sum()
    b = ((s1w["dm"]=="DM1") & (s1w["recurrence"]==0)).sum()
    c = ((s1w["dm"]!="DM1") & (s1w["recurrence"]==1)).sum()
    d = ((s1w["dm"]!="DM1") & (s1w["recurrence"]==0)).sum()
    fisher_summary([[a,b],[c,d]], "recurrence", sname, "DM1_vs_others_OR")

# 4-way DM × TERT interaction in S1
sub_int = mg[S1 & mg["dm"].isin(["DM1","DM2","not_DM"])].copy()
sub_int["dm_x_tert"] = sub_int["is_DM1"] * sub_int["tert_mut"]
sub_int["is_DM1_int"] = sub_int["is_DM1"]
for ev_col, t_col, label in ENDPOINTS:
    cox_fit(sub_int, t_col, ev_col,
            ["is_DM1_int", "tert_mut", "dm_x_tert", "age_yr", "adv_stage", "sex_m"],
            label="DM1_x_TERT_interaction_S1", stratum="S1_BRAF_cPTC",
            endpoint=label, contrast="dm_x_tert")

# write Cox table
ctab = pd.DataFrame(ROWS)
ctab.to_csv(OUT/"h6_cox_table.tsv", sep="\t", index=False)
print(f"\n[ok] wrote {OUT/'h6_cox_table.tsv'} | rows={len(ctab)}")

# ============================================================
# KM curve — headline: BRAF-cPTC × DM (DM1 vs DM2 vs not_DM) on PFI
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))

# Panel A: PFI in BRAF-cPTC stratified by DM
ax = axes[0]
sub = mg[S1].dropna(subset=["PFI", "PFI.time", "dm"]).copy()
sub = sub[sub["PFI.time"]>0]
colors = {"DM1": "#d62728", "DM2": "#1f77b4", "not_DM": "#7f7f7f"}
for grp in ["DM1", "DM2", "not_DM"]:
    g = sub[sub["dm"]==grp]
    if len(g)<3: continue
    kmf = KaplanMeierFitter()
    kmf.fit(g["PFI.time"]/365.25, g["PFI"],
            label=f"{grp} (n={len(g)}, ev={int(g['PFI'].sum())})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color=colors[grp])
# logrank DM1 vs DM2
g1 = sub[sub["dm"]=="DM1"]; g2 = sub[sub["dm"]=="DM2"]
if len(g1)>2 and len(g2)>2:
    lr = logrank_test(g1["PFI.time"], g2["PFI.time"], g1["PFI"], g2["PFI"])
    ax.text(0.02, 0.04, f"DM1 vs DM2 log-rank p={lr.p_value:.3g}",
            transform=ax.transAxes, fontsize=9)
ax.set_title("BRAF-cPTC stratum — PFI by DM status (TCGA-THCA)")
ax.set_xlabel("Years"); ax.set_ylabel("Progression-free probability")
ax.set_ylim(0.5, 1.02); ax.grid(alpha=0.3)

# Panel B: PFI in BRAF-cPTC × TERT × DM (4-way)
ax = axes[1]
sub2 = mg[S1].dropna(subset=["PFI", "PFI.time", "dm"]).copy()
sub2 = sub2[sub2["PFI.time"]>0]
sub2["grp4"] = (
    sub2["dm"].map({"DM1":"DM1","DM2":"DM2","not_DM":"notDM"}).astype(str)
    + "_" + sub2["tert_mut"].map({1:"TERT+", 0:"TERT-"})
)
combo_colors = {
    "DM1_TERT+": "#8b0000", "DM1_TERT-": "#d62728",
    "DM2_TERT+": "#08306b", "DM2_TERT-": "#1f77b4",
    "notDM_TERT+": "#525252", "notDM_TERT-": "#bdbdbd",
}
for grp in ["DM1_TERT+", "DM1_TERT-", "DM2_TERT-", "notDM_TERT+", "notDM_TERT-"]:
    g = sub2[sub2["grp4"]==grp]
    if len(g)<3: continue
    kmf = KaplanMeierFitter()
    kmf.fit(g["PFI.time"]/365.25, g["PFI"],
            label=f"{grp} (n={len(g)}, ev={int(g['PFI'].sum())})")
    kmf.plot_survival_function(ax=ax, ci_show=False, color=combo_colors.get(grp,"black"))
ax.set_title("BRAF-cPTC — PFI × DM × TERT (4-way)")
ax.set_xlabel("Years"); ax.set_ylabel("Progression-free probability")
ax.set_ylim(0.4, 1.02); ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(OUT/"h6_km_braf_cptc_dm.png", dpi=180)
plt.close()
print(f"[ok] wrote {OUT/'h6_km_braf_cptc_dm.png'}")

# Save raw merged data for record
mg[["sample_id","sample_short","dm","molecular_subtype","histology_subtype",
    "tert_promoter_integrated","tert_mut",
    "OS","OS.time","PFI","PFI.time","DFI","DFI.time","DSS","DSS.time",
    "recurrence","age_yr","adv_stage","sex_m","rai_score_v17","tds16_score_v17"]].to_csv(
    OUT/"h6_merged_clinical.tsv", sep="\t", index=False)

# Headline numbers for report
def get_row(stratum, endpoint, model, cov):
    m = ctab[(ctab["stratum"]==stratum)&(ctab["endpoint"]==endpoint)
            &(ctab["model"]==model)&(ctab["covariate"]==cov)]
    if len(m)==0: return None
    return m.iloc[0]

print("\n[HEADLINE]")
for ep in ["OS","PFI","DFI","DSS"]:
    for stratum in ["S1_BRAF_cPTC","S2_BRAF_any","S3_full"]:
        r = get_row(stratum, ep, "DM1_vs_DM2_adj_age_stage_sex", "is_DM1")
        if r is None: continue
        print(f"  [{stratum:18s}] {ep:4s} DM1_vs_DM2: HR={r['HR']:.2f} "
              f"[{r['CI_low']:.2f},{r['CI_high']:.2f}] p={r['p']:.3g} "
              f"n={r['n']} ev={r['events']} {r['note']}")
print()
for ep in ["OS","PFI","DFI","DSS"]:
    r = get_row("S1_BRAF_cPTC", ep, "DM1_AND_TERT_jointly_adj_S1", "is_DM1_wide")
    rT = get_row("S1_BRAF_cPTC", ep, "DM1_AND_TERT_jointly_adj_S1", "tert_mut")
    if r is None: continue
    print(f"  [S1 joint     ] {ep:4s} DM1_wide: HR={r['HR']:.2f} [{r['CI_low']:.2f},{r['CI_high']:.2f}] p={r['p']:.3g}  "
          f"|  TERT: HR={rT['HR']:.2f} [{rT['CI_low']:.2f},{rT['CI_high']:.2f}] p={rT['p']:.3g}  "
          f"n={r['n']} ev={r['events']}")

# write json with numeric headlines
hl = {}
for ep in ["OS","PFI","DFI","DSS"]:
    for stratum in ["S1_BRAF_cPTC","S2_BRAF_any","S3_full"]:
        r = get_row(stratum, ep, "DM1_vs_DM2_adj_age_stage_sex", "is_DM1")
        if r is not None:
            hl[f"{stratum}_{ep}_DM1vsDM2"] = {
                "HR": r["HR"], "CI_low": r["CI_low"], "CI_high": r["CI_high"],
                "p": r["p"], "n": int(r["n"]), "events": int(r["events"]),
                "note": r["note"],
            }
    r = get_row("S1_BRAF_cPTC", ep, "DM1_AND_TERT_jointly_adj_S1", "is_DM1_wide")
    rT = get_row("S1_BRAF_cPTC", ep, "DM1_AND_TERT_jointly_adj_S1", "tert_mut")
    if r is not None:
        hl[f"S1_{ep}_joint_DM1"] = {
            "HR": r["HR"], "CI_low": r["CI_low"], "CI_high": r["CI_high"],
            "p": r["p"], "n": int(r["n"]), "events": int(r["events"]), "note": r["note"]}
    if rT is not None:
        hl[f"S1_{ep}_joint_TERT"] = {
            "HR": rT["HR"], "CI_low": rT["CI_low"], "CI_high": rT["CI_high"],
            "p": rT["p"], "n": int(rT["n"]), "events": int(rT["events"]), "note": rT["note"]}

with open(OUT/"h6_headline.json","w") as f:
    json.dump(hl, f, indent=2, default=lambda x: None if pd.isna(x) else float(x))

print(f"\n[ok] wrote {OUT/'h6_headline.json'}")
print("DONE")
