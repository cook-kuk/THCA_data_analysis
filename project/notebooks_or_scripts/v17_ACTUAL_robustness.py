#!/usr/bin/env python3
"""v17 ACTUAL ship — 3 critical robustness tasks combined.

A1-A Stage-adjusted Cox + Firth-like (lifelines penalizer)
A1-B Bootstrap CI for 4-group logrank (1000 iter)
A1-C Leave-one-out sensitivity (36 LOO + 6-event LOO)

Outputs to results/v17_actual/.
"""
from __future__ import annotations
import json
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
warnings.filterwarnings("ignore")

ROOT = Path("/opt/thyroid-dash/project")
SRC  = ROOT / "results" / "v17_tert_recovery" / "v2" / "sample_master_v17_tert_v2.tsv"
OUT  = ROOT / "results" / "v17_actual"
OUT.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(42)

# --------------------------------------------------------------------
# Load + prep
# --------------------------------------------------------------------
df = pd.read_csv(SRC, sep="\t")
df["tert_pos"] = (df["tert_promoter_integrated"].astype(str).str.lower() == "mutated").astype(int)
df["has_surv"] = df["os_event"].notna() & df["os_days"].notna()
df_surv = df[df["has_surv"]].copy()
df_surv["event"] = df_surv["os_event"].astype(int)
df_surv["time"] = df_surv["os_days"].astype(float)
df_surv["age"] = pd.to_numeric(df_surv["age_clinical"].fillna(df_surv.get("age_at_diagnosis", np.nan)), errors="coerce")

# Stage encoded: I=1, II=2, III=3, IV=4 (Stage 0a + NaN dropped to NaN)
def stage_to_int(s):
    s = str(s)
    if "IV" in s: return 4
    if "III" in s: return 3
    if "II" in s: return 2
    if "I" in s and "0" not in s: return 1
    return np.nan
df_surv["stage_int"] = df_surv["ajcc_stage_group"].map(stage_to_int)
df_surv["stage_advanced"] = (df_surv["stage_int"] >= 3).astype(float)
df_surv["sex_m"] = (df_surv["sex_clinical"].astype(str).str.upper() == "MALE").astype(int)

# Paper 4-group: A_BRAF_only / B_RAS_only / C_TERT_pos / D_triple_neg
def paper_group(row):
    if row["tert_pos"]: return "C_TERT_pos"
    q = row["quad_group"]
    if q == "A_braf_only": return "A_BRAF_only"
    if q == "B_ras_only":  return "B_RAS_only"
    if q == "D_triple_negative": return "D_triple_neg"
    return None
df_surv["paper_group"] = df_surv.apply(paper_group, axis=1)
df_surv = df_surv.dropna(subset=["paper_group"])

print(f"=== Cohort summary ===")
print(f"  N (with survival): {len(df_surv)}")
print(f"  Total events:      {df_surv['event'].sum()}")
print(f"  TERT+ N:           {df_surv['tert_pos'].sum()}")
print(f"  TERT+ events:      {((df_surv['tert_pos']==1) & (df_surv['event']==1)).sum()}")
print(f"  TERT+ Stage III/IV:{((df_surv['tert_pos']==1) & (df_surv['stage_advanced']==1)).sum()} / {df_surv['tert_pos'].sum()}")
print(f"  Median follow-up:  {df_surv['time'].median():.0f} days")
print()
print(df_surv.groupby("paper_group").agg(N=("event","size"), events=("event","sum"),
       med_age=("age","median"), pct_advanced=("stage_advanced","mean")).round(2))

# ====================================================================
# A1-A: Cox models
# ====================================================================
def fit_cox(data: pd.DataFrame, covars: list[str], penalizer: float = 0.01,
            label: str = "") -> dict:
    """Fit Cox with given covariates; lifelines penalizer ≈ Firth-like
    correction for small-event count. Return summary of `tert_pos` row (or first)."""
    cph = CoxPHFitter(penalizer=penalizer)
    cols = ["time", "event"] + covars
    sub = data[cols].dropna()
    try:
        cph.fit(sub, duration_col="time", event_col="event")
        s = cph.summary
        out = {"label": label, "n": len(sub), "events": int(sub["event"].sum())}
        # Get TERT row if exists, else first var
        target = "tert_pos" if "tert_pos" in covars else covars[0]
        if target in s.index:
            out["coef"] = float(s.loc[target, "coef"])
            out["HR"] = float(s.loc[target, "exp(coef)"])
            out["HR_lo"] = float(s.loc[target, "exp(coef) lower 95%"])
            out["HR_hi"] = float(s.loc[target, "exp(coef) upper 95%"])
            out["p"] = float(s.loc[target, "p"])
            out["se"] = float(s.loc[target, "se(coef)"])
        out["covars"] = covars
        out["concordance"] = float(cph.concordance_index_)
        return out
    except Exception as e:
        return {"label": label, "n": len(sub), "error": f"{type(e).__name__}: {e}"}

print("\n=== A1-A: Cox models on TERT (penalizer=0.01) ===")
cox_results = []
cox_results.append(fit_cox(df_surv, ["tert_pos"],                                 label="TERT only"))
cox_results.append(fit_cox(df_surv, ["tert_pos", "age"],                          label="TERT + age"))
cox_results.append(fit_cox(df_surv, ["tert_pos", "stage_int"],                    label="TERT + stage"))
cox_results.append(fit_cox(df_surv, ["tert_pos", "stage_int", "age", "sex_m"],    label="TERT + stage + age + sex"))

cox_df = pd.DataFrame(cox_results)
cox_df.to_csv(OUT / "A1A_cox_unadjusted.tsv", sep="\t", index=False)
print(cox_df[["label","n","events","HR","HR_lo","HR_hi","p","concordance"]].round(4).to_string(index=False))

# 4-group analysis
print("\n=== A1-A: 4-group dummy-encoded Cox ===")
groups_dummy = pd.get_dummies(df_surv["paper_group"], prefix="g").astype(int)
dfg = pd.concat([df_surv[["time","event","stage_int","age","sex_m"]], groups_dummy], axis=1)
# Drop reference (BRAF only) to avoid singular
ref_col = "g_A_BRAF_only"
dfg2 = dfg.drop(columns=[ref_col])
covars_4g = [c for c in dfg2.columns if c.startswith("g_")] + ["stage_int","age","sex_m"]
cph = CoxPHFitter(penalizer=0.01)
sub = dfg2[["time","event"] + covars_4g].dropna()
cph.fit(sub, duration_col="time", event_col="event")
g4_summary = cph.summary[["coef","exp(coef)","exp(coef) lower 95%","exp(coef) upper 95%","p"]]
g4_summary.columns = ["coef","HR","HR_lo","HR_hi","p"]
g4_summary.index.name = "covariate"
g4_summary.reset_index().to_csv(OUT / "A1A_4group_cox_adjusted.tsv", sep="\t", index=False, float_format="%.4g")
print(g4_summary.round(4))

# Save full multivariate
cox_full = cox_results[3]  # TERT + stage + age + sex
with open(OUT / "A1A_cox_full_multivariate.tsv", "w") as f:
    f.write("metric\tvalue\n")
    for k, v in cox_full.items():
        f.write(f"{k}\t{v}\n")

# Decision logic
print("\n=== A1-A: Decision logic ===")
adj = cox_results[3]  # TERT + stage + age + sex
hr = adj.get("HR", float("nan"))
hr_lo = adj.get("HR_lo", float("nan"))
hr_hi = adj.get("HR_hi", float("nan"))
p = adj.get("p", float("nan"))
if hr_lo > 1.0 and hr > 2.0 and p < 0.05:
    scenario = "A"
    decision = ("TERT HR > 2 with 95% CI excluding 1, even after stage + age + sex adjustment. "
                "TERT adds independent prognostic value over stage. Paper narrative survives strongly.")
elif 1.5 <= hr < 3 and (hr_lo < 1.0 or 0.05 <= p < 0.2):
    scenario = "B"
    decision = ("TERT HR ~ 1.5-3 but 95% CI includes 1 or p borderline (0.05-0.2). "
                "TERT correlates with stage. Reframe as molecular handle for Stage III/IV.")
else:
    scenario = "C"
    decision = ("TERT HR ~ 1 or 95% CI very wide / p >> 0.05 after stage adjustment. "
                "TERT signal mostly captured by stage. Recommend TERT 4-group → supplementary.")

print(f"  Scenario: {scenario}")
print(f"  Adjusted HR: {hr:.2f} (95% CI {hr_lo:.2f}–{hr_hi:.2f}), p = {p:.4g}")
print(f"  Decision: {decision}")

# Save A1-A summary
A1A = {
    "tert_n": int(df_surv["tert_pos"].sum()),
    "tert_events": int(((df_surv["tert_pos"]==1) & (df_surv["event"]==1)).sum()),
    "tert_stage_advanced_pct": float(((df_surv["tert_pos"]==1) & (df_surv["stage_advanced"]==1)).sum() / df_surv["tert_pos"].sum()),
    "univariate_HR": cox_results[0].get("HR"),
    "univariate_p": cox_results[0].get("p"),
    "multivariate_HR": hr,
    "multivariate_HR_95CI": [hr_lo, hr_hi],
    "multivariate_p": p,
    "scenario": scenario,
    "decision": decision,
    "concordance_full": cox_results[3].get("concordance"),
}
(OUT / "A1A_summary.json").write_text(json.dumps(A1A, indent=2))
(OUT / "A1A_decision.md").write_text(f"""# A1-A Stage-Adjusted Cox — Decision

**Scenario: {scenario}**

| Model | HR (TERT+) | 95% CI | p-value | Concordance |
|---|---:|---|---:|---:|
| TERT only (univariate) | {cox_results[0].get('HR', float('nan')):.2f} | {cox_results[0].get('HR_lo', float('nan')):.2f}–{cox_results[0].get('HR_hi', float('nan')):.2f} | {cox_results[0].get('p', float('nan')):.4g} | {cox_results[0].get('concordance', float('nan')):.3f} |
| TERT + age | {cox_results[1].get('HR', float('nan')):.2f} | {cox_results[1].get('HR_lo', float('nan')):.2f}–{cox_results[1].get('HR_hi', float('nan')):.2f} | {cox_results[1].get('p', float('nan')):.4g} | {cox_results[1].get('concordance', float('nan')):.3f} |
| TERT + stage | {cox_results[2].get('HR', float('nan')):.2f} | {cox_results[2].get('HR_lo', float('nan')):.2f}–{cox_results[2].get('HR_hi', float('nan')):.2f} | {cox_results[2].get('p', float('nan')):.4g} | {cox_results[2].get('concordance', float('nan')):.3f} |
| **TERT + stage + age + sex (multivariate)** | **{hr:.2f}** | **{hr_lo:.2f}–{hr_hi:.2f}** | **{p:.4g}** | {cox_results[3].get('concordance', float('nan')):.3f} |

## Decision

{decision}

## TERT cohort baseline

- N (TERT+): 36 / 504 (7.1 %)
- Events: 6 (16.7 % of TERT+, vs 10/468 = 2.1 % WT)
- Stage III/IV in TERT+: {A1A['tert_stage_advanced_pct']*100:.0f} %

## Implication for Manuscript v3 (R8 paragraph)

See `submission/npj/manuscript_v3.md` — narrative is reframed per scenario {scenario}.
""")

print(f"  → wrote {OUT/'A1A_summary.json'} and A1A_decision.md")
print(f"  → scenario: {scenario}")

# ====================================================================
# A1-B: Bootstrap CI for 4-group logrank
# ====================================================================
print("\n=== A1-B: Bootstrap 4-group logrank (1000 iter) ===")
N_BOOT = 1000
groups = df_surv["paper_group"].values
times = df_surv["time"].values
events = df_surv["event"].values

def bootstrap_iter(seed: int) -> dict:
    rng = np.random.default_rng(seed)
    # Stratified resample per group
    idx = []
    for g in np.unique(groups):
        gi = np.where(groups == g)[0]
        idx.extend(rng.choice(gi, size=len(gi), replace=True))
    idx = np.array(idx)
    try:
        res = multivariate_logrank_test(times[idx], groups[idx], events[idx])
        p = float(res.p_value)
    except Exception:
        p = np.nan
    # Effect-size proxy: TERT+ event rate / triple-neg event rate
    g_tert = (groups[idx] == "C_TERT_pos")
    g_tn   = (groups[idx] == "D_triple_neg")
    er_tert = events[idx][g_tert].mean() if g_tert.any() else np.nan
    er_tn   = events[idx][g_tn].mean() if g_tn.any() else np.nan
    rr = er_tert / er_tn if er_tn > 0 else np.nan
    return {"seed": seed, "p": p, "rate_tert": er_tert, "rate_tn": er_tn, "rate_ratio": rr}

boot = pd.DataFrame([bootstrap_iter(i) for i in range(N_BOOT)])
boot.to_csv(OUT / "A1B_bootstrap_logrank.tsv", sep="\t", index=False, float_format="%.4g")

p_med = boot["p"].median()
p_lo  = boot["p"].quantile(0.025)
p_hi  = boot["p"].quantile(0.975)
pct_below_05 = (boot["p"] < 0.05).mean() * 100
pct_below_1e4 = (boot["p"] < 1e-4).mean() * 100

rr_med = boot["rate_ratio"].median()
rr_lo = boot["rate_ratio"].quantile(0.025)
rr_hi = boot["rate_ratio"].quantile(0.975)

A1B = {
    "n_boot": N_BOOT,
    "p_median": float(p_med),
    "p_95ci": [float(p_lo), float(p_hi)],
    "pct_p_lt_0p05": float(pct_below_05),
    "pct_p_lt_1e4": float(pct_below_1e4),
    "rate_ratio_median": float(rr_med),
    "rate_ratio_95ci": [float(rr_lo), float(rr_hi)],
    "robust": bool(p_hi < 0.01 and rr_lo > 2),
}
(OUT / "A1B_summary.json").write_text(json.dumps(A1B, indent=2))
print(f"  median p = {p_med:.2g}, 95% CI [{p_lo:.2g}, {p_hi:.2g}]")
print(f"  % iter p<0.05: {pct_below_05:.1f}, % iter p<1e-4: {pct_below_1e4:.1f}")
print(f"  TERT+/triple-neg event-rate ratio: median {rr_med:.2f}, 95% CI [{rr_lo:.2f}, {rr_hi:.2f}]")
print(f"  Robust (CI hi p < 0.01 AND CI lo rate-ratio > 2): {A1B['robust']}")

# Effect-size CI tsv
pd.DataFrame({"metric": ["rate_ratio_median","rate_ratio_lo95","rate_ratio_hi95",
                          "p_median","p_lo95","p_hi95","pct_p_lt_05"],
              "value":  [rr_med, rr_lo, rr_hi, p_med, p_lo, p_hi, pct_below_05]
             }).to_csv(OUT / "A1B_effect_size_ci.tsv", sep="\t", index=False,
                       float_format="%.4g")

# ====================================================================
# A1-C: Leave-one-out sensitivity
# ====================================================================
print("\n=== A1-C: Leave-one-out sensitivity ===")
loo_results = []
tert_pos_idx = df_surv.index[df_surv["tert_pos"] == 1].tolist()
event_in_tert_idx = df_surv.index[(df_surv["tert_pos"] == 1) & (df_surv["event"] == 1)].tolist()

# 36 LOO: drop one TERT+ patient at a time, redo logrank
for i, pid in enumerate(tert_pos_idx):
    sub = df_surv.drop(index=pid)
    res = multivariate_logrank_test(sub["time"], sub["paper_group"], sub["event"])
    is_event_loo = pid in event_in_tert_idx
    loo_results.append({
        "iteration": i,
        "patient_idx": int(pid),
        "is_event": int(is_event_loo),
        "logrank_p": float(res.p_value),
        "logrank_chi2": float(res.test_statistic),
    })
loo_df = pd.DataFrame(loo_results)
loo_df.to_csv(OUT / "A1C_loo_logrank.tsv", sep="\t", index=False, float_format="%.4g")

p_full = float(multivariate_logrank_test(df_surv["time"], df_surv["paper_group"], df_surv["event"]).p_value)
p_max = loo_df["logrank_p"].max()
p_max_event_only = loo_df[loo_df["is_event"]==1]["logrank_p"].max()
n_below_05 = (loo_df["logrank_p"] < 0.05).sum()
n_below_001 = (loo_df["logrank_p"] < 0.001).sum()

# 6 event LOO subset
event_loo = loo_df[loo_df["is_event"]==1].copy().reset_index(drop=True)
event_loo.to_csv(OUT / "A1C_event_loo.tsv", sep="\t", index=False, float_format="%.4g")

# Influence score: |log10 p_LOO - log10 p_full|
loo_df["influence"] = (np.log10(loo_df["logrank_p"]) - np.log10(p_full)).abs()
loo_df.sort_values("influence", ascending=False).to_csv(
    OUT / "A1C_influence_score.tsv", sep="\t", index=False, float_format="%.4g")

A1C = {
    "p_full":     p_full,
    "n_loo":      len(loo_df),
    "p_max_loo":  float(p_max),
    "p_max_event_only_loo": float(p_max_event_only),
    "n_loo_p_lt_0p05":  int(n_below_05),
    "n_loo_p_lt_0p001": int(n_below_001),
    "robust_36":  bool(p_max < 0.05),
    "robust_event":  bool(p_max_event_only < 0.01),
    "fragile_warn": bool(p_max_event_only > 0.05),
}
(OUT / "A1C_summary.json").write_text(json.dumps(A1C, indent=2))
print(f"  Full-cohort logrank p: {p_full:.2g}")
print(f"  Worst LOO p (any 1 patient dropped): {p_max:.4g}")
print(f"  Worst LOO p (event patients only):   {p_max_event_only:.4g}")
print(f"  LOO iter with p<0.05 / p<0.001: {n_below_05}/{n_below_001} of {len(loo_df)}")
print(f"  Robust (all 36 LOO p<0.05): {A1C['robust_36']}")
print(f"  Robust (all event LOO p<0.01): {A1C['robust_event']}")
print(f"  Fragile warn (any event LOO p>0.05): {A1C['fragile_warn']}")

# ====================================================================
# Final scenario decision combining all 3
# ====================================================================
final = {
    "A1A": A1A, "A1B": A1B, "A1C": A1C,
    "scenario": A1A["scenario"],
    "submit_verdict": ("GO_strong" if (A1A["scenario"] == "A" and A1B["robust"] and A1C["robust_event"])
                       else "GO_with_reframe" if A1A["scenario"] in ("A", "B")
                       else "REVISE_to_supplementary"),
}
(OUT / "FINAL_decision.json").write_text(json.dumps(final, indent=2, default=str))

print("\n" + "="*60)
print(f"FINAL VERDICT: {final['submit_verdict']}  (scenario {A1A['scenario']})")
print("="*60)
