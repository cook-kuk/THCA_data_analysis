#!/usr/bin/env python3
"""D3-P5 — GSE286332 18-sample P(DM1) gradient driver analysis.

Tests whether the P(DM1) 0.05–0.30 spectrum within the 18-sample GSE286332 cohort
is driven by HLA-II infiltration, 8-gene RAI dedifferentiation, or generic immune,
and whether the PTC+HT → P_DM1 ↓ relationship is mediated by HLA-II (Baron-Kenny).
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

PROJ = Path("/opt/thyroid-dash/project")
RES = PROJ / "results/d3p5_pdm1_gradient"
RES.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. Load
# ============================================================
print("=== 1. Load GSE286332 P(DM1), 8-gene, HLA-I/II, immune ===")
dm = pd.read_csv(PROJ / "results/p3_gse286332/dm12_predictions.tsv", sep="\t", index_col=0)
panel = pd.read_csv(PROJ / "results/p3_gse286332/8gene_panel_per_sample.tsv", sep="\t", index_col=0)
hla = pd.read_csv(PROJ / "results/p3_gse286332/hla_module_scores.tsv", sep="\t", index_col=0)
scores = pd.read_csv(PROJ / "results/p5_8gene_vs_hla_autocorr/scores_GSE286332.tsv", sep="\t", index_col=0)

df = pd.DataFrame({
    "P_DM1": dm["P_DM1"],
    "group": dm["group"],
    "g8_RAI": panel["RAI_score_8gene"],
    "HLA_I": hla["HLA_I"],
    "HLA_II": hla["HLA_II"],
    "immune": scores["immune"],
}, index=dm.index)
df["is_PTC_HT"] = (df["group"] == "PTC_HT").astype(int)
print(df.to_string())
df.to_csv(RES / "merged_18sample.tsv", sep="\t")

# ============================================================
# 2. P_DM1 distribution (per group)
# ============================================================
print("\n=== 2. P_DM1 distribution per group ===")
for g in ["PTC", "PTC_HT"]:
    s = df.loc[df["group"] == g, "P_DM1"]
    print(f"  {g}: n={len(s)}, range=[{s.min():.3f}, {s.max():.3f}], "
          f"median={s.median():.3f}, mean={s.mean():.3f}, sd={s.std(ddof=1):.3f}")

# Mann-Whitney
ng = df.loc[df["group"] == "PTC", "P_DM1"].values
th = df.loc[df["group"] == "PTC_HT", "P_DM1"].values
mw_u, mw_p = stats.mannwhitneyu(ng, th, alternative="two-sided")
print(f"  MW PTC vs PTC+HT: U={mw_u:.1f}, p={mw_p:.4g}")

# ============================================================
# 3. P_DM1 vs continuous covariates (Spearman)
# ============================================================
print("\n=== 3. Spearman ρ between P_DM1 and covariates ===")
covs = ["g8_RAI", "HLA_I", "HLA_II", "immune"]
spear_rows = []
for c in covs:
    rho, p = stats.spearmanr(df["P_DM1"], df[c])
    print(f"  P_DM1 ~ {c}: ρ={rho:+.3f}, p={p:.4g}")
    spear_rows.append(dict(covariate=c, spearman_rho=round(rho, 3), p=p))
pd.DataFrame(spear_rows).to_csv(RES / "spearman_pdm1_vs_covariates.tsv", sep="\t", index=False)

# ============================================================
# 4. Linear regression decomposition (P_DM1 ~ HLA-II + g8 + immune)
# ============================================================
print("\n=== 4. Linear regression decomposition ===")
# Standardize each predictor for comparability of coefficients
z = df[covs].apply(lambda c: (c - c.mean()) / (c.std(ddof=1) + 1e-9))
X = sm.add_constant(z[["HLA_II", "g8_RAI", "immune"]])
y = df["P_DM1"]
fit = sm.OLS(y, X).fit()
print(fit.summary())
decomp = {
    "model": "P_DM1 ~ HLA_II + g8_RAI + immune (z-standardized)",
    "R2": round(fit.rsquared, 3),
    "R2_adj": round(fit.rsquared_adj, 3),
    "n": int(fit.nobs),
    "coefficients": {k: dict(beta=round(v, 4), se=round(fit.bse[k], 4),
                              t=round(fit.tvalues[k], 3), p=float(fit.pvalues[k]))
                     for k, v in fit.params.items()},
}
print(f"\nR² = {decomp['R2']:.3f}, R² adj = {decomp['R2_adj']:.3f}")
print(f"Dominant predictor: {max(['HLA_II', 'g8_RAI', 'immune'], key=lambda k: abs(fit.params[k]))}")

# Single-predictor R² for each
single_r2 = {}
for c in ["HLA_II", "g8_RAI", "immune"]:
    Xs = sm.add_constant(z[[c]])
    f = sm.OLS(y, Xs).fit()
    single_r2[c] = round(f.rsquared, 3)
    print(f"  Single-predictor R²: P_DM1 ~ {c}: {f.rsquared:.3f}")
decomp["single_predictor_R2"] = single_r2

# Hierarchical: HLA_II first, then add g8, then immune
print("\n  Hierarchical R² addition:")
prev_r2 = 0.0
hier = []
for col in ["HLA_II", "g8_RAI", "immune"]:
    used = [c for c in ["HLA_II", "g8_RAI", "immune"] if c == col or c in [x[0] for x in hier]]
    Xh = sm.add_constant(z[used])
    fh = sm.OLS(y, Xh).fit()
    delta = fh.rsquared - prev_r2
    print(f"    +{col}: R²={fh.rsquared:.3f} (Δ={delta:.3f})")
    hier.append((col, round(fh.rsquared, 3), round(delta, 3)))
    prev_r2 = fh.rsquared
decomp["hierarchical_R2"] = hier

(RES / "decomposition.json").write_text(json.dumps(decomp, indent=2, default=str))

# ============================================================
# 5. Mediation: PTC+HT → HLA_II → P_DM1 (Baron-Kenny + bootstrap)
# ============================================================
print("\n=== 5. Mediation analysis (Baron-Kenny + bootstrap) ===")

def baron_kenny(treat, mediator, outcome):
    """Returns total, direct, indirect, %mediated."""
    Xc = sm.add_constant(treat)
    fit_c = sm.OLS(outcome, Xc).fit()  # total effect c
    c_total = fit_c.params[1]

    Xa = sm.add_constant(treat)
    fit_a = sm.OLS(mediator, Xa).fit()  # treat → mediator
    a = fit_a.params[1]

    Xb = sm.add_constant(np.column_stack([treat, mediator]))
    fit_b = sm.OLS(outcome, Xb).fit()  # treat + mediator → outcome
    b = fit_b.params[2]
    c_direct = fit_b.params[1]

    indirect = a * b
    pct_mediated = 100 * indirect / c_total if abs(c_total) > 1e-9 else np.nan
    return c_total, c_direct, indirect, pct_mediated, a, b

def boot_mediation(treat, mediator, outcome, n_boot=5000, seed=42):
    rng = np.random.default_rng(seed)
    n = len(treat)
    indirects = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        try:
            _, _, ind, _, _, _ = baron_kenny(treat[idx], mediator[idx], outcome[idx])
            indirects.append(ind)
        except Exception:
            continue
    indirects = np.array(indirects)
    return dict(mean=np.mean(indirects),
                ci_lo=np.percentile(indirects, 2.5),
                ci_hi=np.percentile(indirects, 97.5),
                p_emp=2 * min(np.mean(indirects > 0), np.mean(indirects < 0)))

treat = df["is_PTC_HT"].values.astype(float)
y_outcome = df["P_DM1"].values
mediation_results = {}

for med_name in ["HLA_II", "g8_RAI", "immune", "HLA_I"]:
    med = df[med_name].values
    c_total, c_direct, indirect, pct, a, b = baron_kenny(treat, med, y_outcome)
    bb = boot_mediation(treat, med, y_outcome, n_boot=5000, seed=42)
    print(f"\n  Mediator: {med_name}")
    print(f"    a (treat → mediator) = {a:+.4f}")
    print(f"    b (mediator → outcome | treat) = {b:+.4f}")
    print(f"    c_total = {c_total:+.4f}")
    print(f"    c_direct = {c_direct:+.4f}")
    print(f"    indirect (a*b) = {indirect:+.4f}")
    print(f"    %mediated = {pct:+.1f}%")
    print(f"    bootstrap indirect mean={bb['mean']:+.4f}, 95% CI [{bb['ci_lo']:+.4f}, {bb['ci_hi']:+.4f}], p_emp={bb['p_emp']:.4g}")
    mediation_results[med_name] = dict(
        a=round(float(a), 4), b=round(float(b), 4),
        c_total=round(float(c_total), 4), c_direct=round(float(c_direct), 4),
        indirect=round(float(indirect), 4), pct_mediated=round(float(pct), 1),
        boot_indirect_mean=round(float(bb["mean"]), 4),
        boot_ci_lo=round(float(bb["ci_lo"]), 4),
        boot_ci_hi=round(float(bb["ci_hi"]), 4),
        boot_p_emp=float(bb["p_emp"]),
    )

(RES / "mediation_results.json").write_text(json.dumps(mediation_results, indent=2))

# ============================================================
# 6. Continuous PTC+HT severity proxy (HLA-II quartile within full 18)
# ============================================================
print("\n=== 6. PTC+HT severity gradient (within-PTC test) ===")
# Within 9 PTC samples, is there pre-clinical Hashimoto signal?
ptc_only = df[df["group"] == "PTC"].copy()
print(f"  PTC subset n={len(ptc_only)}")
print(f"  PTC HLA-II range: [{ptc_only['HLA_II'].min():.3f}, {ptc_only['HLA_II'].max():.3f}]")
print(f"  PTC P_DM1 range: [{ptc_only['P_DM1'].min():.3f}, {ptc_only['P_DM1'].max():.3f}]")

if ptc_only["HLA_II"].std() > 0.05:
    rho_in, p_in = stats.spearmanr(ptc_only["HLA_II"], ptc_only["P_DM1"])
    print(f"  Within-PTC Spearman ρ(HLA_II, P_DM1) = {rho_in:+.3f}, p={p_in:.4g}")
    rho_in2, p_in2 = stats.spearmanr(ptc_only["g8_RAI"], ptc_only["P_DM1"])
    print(f"  Within-PTC Spearman ρ(g8_RAI, P_DM1) = {rho_in2:+.3f}, p={p_in2:.4g}")
    within_ptc = dict(spearman_HLA2_P_DM1=round(float(rho_in), 3), p=float(p_in),
                       spearman_g8_P_DM1=round(float(rho_in2), 3), p2=float(p_in2))
else:
    within_ptc = dict(note="HLA_II variation within PTC too small")

(RES / "within_ptc_severity.json").write_text(json.dumps(within_ptc, indent=2))

# ============================================================
# 7. Decision summary
# ============================================================
print("\n=== 7. DECISION ===")
hla2_r2 = single_r2["HLA_II"]
g8_r2 = single_r2["g8_RAI"]
hla2_pct = mediation_results["HLA_II"]["pct_mediated"]
g8_pct = mediation_results["g8_RAI"]["pct_mediated"]
hla2_boot_sig = (mediation_results["HLA_II"]["boot_ci_lo"] * mediation_results["HLA_II"]["boot_ci_hi"]) > 0

if hla2_r2 > 0.5 and abs(hla2_pct) > 50 and hla2_boot_sig:
    decision = "STRONG"
    msg = (f"HLA-II is dominant predictor of P(DM1) gradient (single-predictor R²={hla2_r2:.2f}, "
           f"%mediated={hla2_pct:.0f}%, bootstrap CI excludes 0). "
           f"PTC+HT → HLA-II infiltration → P_DM1 ↓ mechanistic claim supported.")
elif hla2_r2 > 0.3 or g8_r2 > 0.3:
    decision = "MODERATE"
    msg = (f"Multi-axis convergence (HLA-II R²={hla2_r2:.2f}, 8-gene R²={g8_r2:.2f}). "
           f"P_DM1 gradient driven by both immune and dedifferentiation axes simultaneously.")
else:
    decision = "WEAK"
    msg = "P_DM1 gradient weakly correlated with covariates at n=18 — binary framework retained."

print(f"  Decision: {decision}")
print(f"  {msg}")

summary = {
    "n_samples": int(len(df)),
    "PTC_n": int((df["group"] == "PTC").sum()),
    "PTC_HT_n": int((df["group"] == "PTC_HT").sum()),
    "P_DM1_PTC_mean": round(float(ng.mean()), 3),
    "P_DM1_PTC_HT_mean": round(float(th.mean()), 3),
    "MW_p": float(mw_p),
    "spearman_pdm1_vs_covariates": spear_rows,
    "decomposition": decomp,
    "mediation": mediation_results,
    "within_ptc_severity": within_ptc,
    "decision": decision,
    "message": msg,
}
(RES / "D3P5_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print(f"\n✓ Outputs to {RES}")
