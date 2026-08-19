#!/usr/bin/env python3
"""AUDIT 09 — robustness of the RAI-anchored recurrence result.

AUDIT 08 produced a split verdict that must not be left as-is:

    Cox, panel score continuous  : HR 0.71 [0.40, 1.27], P = 0.246   -> null
    log-rank, median split       : P = 0.043                          -> nominal

Reporting the second without interrogating it is exactly the selective reading
Prof. Kang warned against. This script asks whether the median split survives
scrutiny, and adds the adjusted and dose analyses that the cohort now permits.

Analyses
  A  maximally selected log-rank across all cut-points, with a permutation-corrected
     p-value (the median is one of many cut-points; choosing it is a free parameter)
  B  dose-response across score tertiles and quartiles (is there monotonicity?)
  C  multivariable Cox: score + stage + age + composition
  D  proportional-hazards check (Schoenfeld-type correlation of scaled residuals with time)
  E  radioiodine DOSE analyses: does the score track cumulative mCi, and does dose confound?
  F  sensitivity: drop patients who also received external-beam radiation
  G  minimum detectable effect for the treatment x score interaction

Outputs
  results/tables/audit09_*_2026_08_06.tsv
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
RADIO = Path("/data/rai_atlas/raw/TCGA_THCA_biotab"
             "/nationwidechildrens.org_clinical_radiation_thca.txt")
SURV = Path("/data/thca/repo_data/raw/TCGA_pancan/survival.tsv")
OUT = ROOT / "results" / "tables"
SEED = 20260806
NPERM = 5000

spec = importlib.util.spec_from_file_location(
    "base", ROOT / "scripts" / "tcga_rai_best_response_2026_08_06.py")
_base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(_base)

a8spec = importlib.util.spec_from_file_location(
    "a8", ROOT / "scripts" / "audit" / "08_time_from_rai_to_recurrence_2026_08_06.py")
_a8 = importlib.util.module_from_spec(a8spec)
a8spec.loader.exec_module(_a8)

logrank = _a8.logrank


# ----------------------------------------------------------------- Cox, multivariable
def cox_fit(t, e, X, max_iter=100, tol=1e-9):
    """Breslow-tie Cox by Newton-Raphson. X is n x p (no intercept)."""
    t = np.asarray(t, float)
    e = np.asarray(e, float)
    X = np.atleast_2d(np.asarray(X, float))
    if X.shape[0] != len(t):
        X = X.T
    order = np.argsort(t)
    t, e, X = t[order], e[order], X[order]
    n, p = X.shape
    b = np.zeros(p)
    ev = np.where(e == 1)[0]
    for _ in range(max_iter):
        eta = X @ b
        eta = np.clip(eta, -50, 50)
        r = np.exp(eta)
        g = np.zeros(p)
        H = np.zeros((p, p))
        for i in ev:
            at = t >= t[i]
            ra = r[at]
            Xa = X[at]
            s0 = ra.sum()
            s1 = (ra[:, None] * Xa).sum(axis=0)
            s2 = (ra[:, None, None] * Xa[:, :, None] * Xa[:, None, :]).sum(axis=0)
            m = s1 / s0
            g += X[i] - m
            H -= s2 / s0 - np.outer(m, m)
        try:
            step = np.linalg.solve(H, g)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(H, g, rcond=None)[0]
        b_new = b - step
        if np.max(np.abs(b_new - b)) < tol:
            b = b_new
            break
        b = b_new
    try:
        cov = np.linalg.inv(-H)
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(-H)
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(se > 0, b / se, np.nan)
    pv = 2 * (1 - stats.norm.cdf(np.abs(z)))
    return b, se, pv


def martingale_like_resid(t, e, x, beta):
    """Scaled score residuals at each event time, for a crude PH check."""
    order = np.argsort(t)
    t, e, x = np.asarray(t)[order], np.asarray(e)[order], np.asarray(x)[order]
    r = np.exp(beta * x)
    out = []
    for i in np.where(e == 1)[0]:
        at = t >= t[i]
        s0 = r[at].sum()
        s1 = (r[at] * x[at]).sum()
        out.append((t[i], x[i] - s1 / s0))
    return np.array(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    # ------------------------------------------------ rebuild the RAI-anchored frame
    df = _base.build()
    rad = pd.read_csv(RADIO, sep="\t", skiprows=[1, 2], low_memory=False)
    units = rad["radiation_adjuvant_units"].astype(str).str.lower()
    rai = rad[units.eq("mci")].copy()
    rai["start"] = pd.to_numeric(rai["radiation_therapy_started_days_to"], errors="coerce")
    first = (rai.groupby("bcr_patient_barcode")["start"].min().reset_index()
             .rename(columns={"bcr_patient_barcode": "patient", "start": "rai_start_day"}))
    df = df.merge(first, on="patient", how="left")

    sv = pd.read_csv(SURV, sep="\t", low_memory=False)
    sv["patient"] = sv["_PATIENT"]
    sv = sv[["patient", "new_tumor_event_dx_days_to", "last_contact_days_to",
             "age_at_initial_pathologic_diagnosis"]].drop_duplicates("patient")
    for c in sv.columns[1:]:
        sv[c] = pd.to_numeric(sv[c], errors="coerce")
    df = df.merge(sv, on="patient", how="left")

    d = df[df["rai_start_day"].notna()].copy()
    has = d["new_tumor_event_dx_days_to"].notna()
    pre = has & (d["new_tumor_event_dx_days_to"] < d["rai_start_day"])
    d["time"] = np.where(has & ~pre,
                         d["new_tumor_event_dx_days_to"] - d["rai_start_day"],
                         d["last_contact_days_to"] - d["rai_start_day"])
    d["event"] = np.where(has & ~pre, 1, 0)
    d["RAI_8"] = pd.to_numeric(d["RAI_8"], errors="coerce")
    d = d[d["time"].notna() & (d["time"] > 0) & d["RAI_8"].notna()].copy()
    d["stage_high"] = d["AJCC_PATHOLOGIC_TUMOR_STAGE"].astype(str).str.contains(
        "III|IV", regex=True).astype(int)
    d["age"] = pd.to_numeric(d["age_at_initial_pathologic_diagnosis"], errors="coerce")
    d["comp"] = pd.to_numeric(d.get("leuko_frac"), errors="coerce")
    d["dose"] = pd.to_numeric(d["rai_dose_mci"], errors="coerce")

    t, e, x = d["time"].values, d["event"].values, d["RAI_8"].values
    n, nev = len(d), int(e.sum())
    print(f"\nanalysis set n={n}, events={nev}\n")

    rows = []

    # ------------------------------------------------------------------ A  cut-point
    # Vectorised log-rank over all candidate cut-points at once. Patients are sorted
    # by time; at each event the at-risk counts are reverse cumulative sums, so the
    # whole cut-point sweep is a single matrix operation and the permutation loop is
    # feasible. Breslow handling of ties (few ties here).
    lo_q, hi_q = np.quantile(x, [0.10, 0.90])
    cand = np.unique(x[(x >= lo_q) & (x <= hi_q)])

    _o = np.argsort(t, kind="mergesort")
    t_s, e_s = t[_o], e[_o]
    ev_pos = np.where(e_s == 1)[0]
    n_risk = (len(t_s) - np.arange(len(t_s)))[ev_pos].astype(float)

    def max_chi_vec(xv):
        """Max log-rank chi2 over all candidate cut-points. Returns (chi2, cut)."""
        G = (xv[_o][:, None] > cand[None, :])            # n x K, sorted by time
        n1 = np.cumsum(G[::-1], axis=0)[::-1]            # at-risk in group 1
        keep = (G.sum(axis=0) >= 15) & ((~G).sum(axis=0) >= 15)
        if not keep.any():
            return -1.0, np.nan
        n1e = n1[ev_pos][:, keep].astype(float)
        p = n1e / n_risk[:, None]
        O1 = G[ev_pos][:, keep].sum(axis=0).astype(float)
        E1 = p.sum(axis=0)
        V = (p * (1 - p)).sum(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            chi = np.where(V > 0, (O1 - E1) ** 2 / V, np.nan)
        if not np.isfinite(chi).any():
            return -1.0, np.nan
        j = int(np.nanargmax(chi))
        return float(chi[j]), float(cand[keep][j])

    obs_chi, obs_cut = max_chi_vec(x)
    max_chi = max_chi_vec
    med = float(np.median(x))
    med_chi, med_p = logrank(t[x <= med], e[x <= med], t[x > med], e[x > med])

    perm_max = np.empty(NPERM)
    for i in range(NPERM):
        perm_max[i] = max_chi(rng.permutation(x))[0]
    p_corr = float((perm_max >= obs_chi).mean())

    rows.append(dict(analysis="A · median split (as reported in AUDIT 08)",
                     n=n, events=nev, statistic=f"chi2={med_chi:.3f}",
                     p_value=med_p, p_corrected=np.nan,
                     verdict="nominal, uncorrected",
                     note=f"cut at median {med:.4f}"))
    rows.append(dict(analysis="A · MAXIMALLY SELECTED log-rank over all cut-points",
                     n=n, events=nev, statistic=f"chi2_max={obs_chi:.3f}",
                     p_value=float(1 - stats.chi2.cdf(obs_chi, 1)),
                     p_corrected=p_corr,
                     verdict="NULL after correction" if p_corr > 0.05 else "survives correction",
                     note=(f"best cut {obs_cut:.4f}; {NPERM} permutations of the score; "
                           f"naive p ignores that the cut-point was chosen")))
    print(f"A  median chi2={med_chi:.3f} p={med_p:.4f} | "
          f"max chi2={obs_chi:.3f} at cut {obs_cut:.4f} | permutation-corrected p={p_corr:.4f}")

    # ------------------------------------------------------------------ B  monotonic?
    for k, lab in [(3, "tertile"), (4, "quartile")]:
        qs = pd.qcut(d["RAI_8"], k, labels=False, duplicates="drop")
        parts = []
        for g in sorted(pd.unique(qs.dropna())):
            m = qs == g
            parts.append(f"g{int(g)+1}: n={int(m.sum())} ev={int(e[m.values].sum())} "
                         f"({100*e[m.values].sum()/max(m.sum(),1):.1f}%)")
        # trend test: Cox on the ordinal group index
        b, se, pv = cox_fit(t, e, qs.values.astype(float).reshape(-1, 1))
        rows.append(dict(analysis=f"B · {lab} dose-response (Cox trend on ordinal group)",
                         n=n, events=nev,
                         statistic=f"HR/step={np.exp(b[0]):.3f}",
                         p_value=float(pv[0]), p_corrected=np.nan,
                         verdict="monotone trend" if pv[0] < 0.05 else "no monotone trend",
                         note=" | ".join(parts)))
        print(f"B  {lab}: " + " | ".join(parts) + f"  -> trend p={pv[0]:.4f}")

    # ------------------------------------------------------------------ C  adjusted
    for terms, lab in [(["RAI_8"], "score only"),
                       (["RAI_8", "stage_high"], "score + stage"),
                       (["RAI_8", "age"], "score + age"),
                       (["RAI_8", "comp"], "score + composition"),
                       (["RAI_8", "stage_high", "age", "comp"], "score + stage + age + comp")]:
        sub = d[["time", "event"] + terms].dropna()
        if len(sub) < 30 or sub["event"].sum() < 5:
            continue
        b, se, pv = cox_fit(sub["time"].values, sub["event"].values,
                            sub[terms].values)
        rows.append(dict(analysis=f"C · Cox {lab}", n=len(sub),
                         events=int(sub["event"].sum()),
                         statistic=(f"HR={np.exp(b[0]):.3f} "
                                    f"[{np.exp(b[0]-1.96*se[0]):.3f},"
                                    f"{np.exp(b[0]+1.96*se[0]):.3f}]"),
                         p_value=float(pv[0]), p_corrected=np.nan,
                         verdict="significant" if pv[0] < 0.05 else "null",
                         note=f"EPV={sub['event'].sum()/len(terms):.1f}"))
        print(f"C  {lab:32} HR={np.exp(b[0]):.3f} P={pv[0]:.4f} (EPV={sub['event'].sum()/len(terms):.1f})")

    # ------------------------------------------------------------------ D  PH check
    b1, se1, _ = cox_fit(t, e, x.reshape(-1, 1))
    res = martingale_like_resid(t, e, x, b1[0])
    if len(res) > 5:
        rho, pph = stats.spearmanr(res[:, 0], res[:, 1])
        rows.append(dict(analysis="D · proportional-hazards check",
                         n=n, events=nev,
                         statistic=f"rho(resid, time)={rho:.3f}",
                         p_value=float(pph), p_corrected=np.nan,
                         verdict="PH violated" if pph < 0.05 else "no PH violation detected",
                         note="Spearman correlation of score residuals at event times with time; "
                              "a crude Schoenfeld-type check, not a formal test"))
        print(f"D  PH check rho={rho:.3f} P={pph:.4f}")

    # ------------------------------------------------------------------ E  dose
    dd = d[d["dose"].notna()]
    rho_d, p_d = stats.spearmanr(dd["RAI_8"], dd["dose"])
    rows.append(dict(analysis="E · panel score vs cumulative radioiodine dose",
                     n=len(dd), events=int(dd["event"].sum()),
                     statistic=f"Spearman rho={rho_d:.3f}", p_value=float(p_d),
                     p_corrected=np.nan,
                     verdict="score tracks dose" if p_d < 0.05 else "score independent of dose",
                     note=f"median dose {dd['dose'].median():.1f} mCi "
                          f"(IQR {dd['dose'].quantile(.25):.1f}-{dd['dose'].quantile(.75):.1f})"))
    print(f"E  score vs dose rho={rho_d:.3f} P={p_d:.4f}")

    sub = d[["time", "event", "RAI_8", "dose"]].dropna()
    b, se, pv = cox_fit(sub["time"].values, sub["event"].values,
                        sub[["RAI_8", "dose"]].values)
    rows.append(dict(analysis="E · Cox score + cumulative dose", n=len(sub),
                     events=int(sub["event"].sum()),
                     statistic=(f"HR_score={np.exp(b[0]):.3f} "
                                f"[{np.exp(b[0]-1.96*se[0]):.3f},{np.exp(b[0]+1.96*se[0]):.3f}]; "
                                f"HR_dose/100mCi={np.exp(b[1]*100):.3f}"),
                     p_value=float(pv[0]), p_corrected=np.nan,
                     verdict="significant" if pv[0] < 0.05 else "null",
                     note=f"dose P={pv[1]:.4f}"))
    print(f"E  Cox score+dose: HR_score={np.exp(b[0]):.3f} P={pv[0]:.4f}; dose P={pv[1]:.4f}")

    # ------------------------------------------------------------------ F  drop EBRT
    if "had_ebrt" in d.columns:
        s = d[~d["had_ebrt"].astype(bool)]
        b, se, pv = cox_fit(s["time"].values, s["event"].values,
                            s["RAI_8"].values.reshape(-1, 1))
        rows.append(dict(analysis="F · sensitivity, exclude external-beam co-treated",
                         n=len(s), events=int(s["event"].sum()),
                         statistic=(f"HR={np.exp(b[0]):.3f} "
                                    f"[{np.exp(b[0]-1.96*se[0]):.3f},{np.exp(b[0]+1.96*se[0]):.3f}]"),
                         p_value=float(pv[0]), p_corrected=np.nan,
                         verdict="significant" if pv[0] < 0.05 else "null",
                         note=f"dropped {int(d['had_ebrt'].astype(bool).sum())} patients "
                              f"who also received Gy/cGy radiation"))
        print(f"F  drop EBRT (n={len(s)}): HR={np.exp(b[0]):.3f} P={pv[0]:.4f}")

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "audit09_rai_anchored_robustness_2026_08_06.tsv",
               sep="\t", index=False)
    print(f"\nwrote {OUT / 'audit09_rai_anchored_robustness_2026_08_06.tsv'}")
    return out


if __name__ == "__main__":
    main()
